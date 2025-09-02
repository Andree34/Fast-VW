#include "Slow_simplifier_PSLG.hpp"

using std::string;
using std::vector;
using std::ifstream;
using std::istringstream;
using std::getline;
using std::list;
using std::map;
using std::pair;
using std::make_pair;
using std::to_string;
using std::cout;
using std::endl;
using std::runtime_error; // TODO temporarily ???

Slow_simplifier_PSLG::Slow_simplifier_PSLG(const std::string& input_folder_name, bool gen_test, bool auto_simplify) :
    name(
        input_folder_name)
{
    start_time = std::chrono::high_resolution_clock::now();
    end_time = start_time;

    // parse input file into chains (flexible parsing; blank lines separate chains; 2 numbers per line are points)
    const auto chains_points = parse_input_file_as_chains("../data/" + input_folder_name + "/data.in");

    // build internal node list, merge duplicate coordinates and create per-node structures
    build_internal_structures_from_chains(chains_points);

    // initial vertex count is number of node occurrences
    init_vertex_count = static_cast<int>(PI.size());

    // run simplification if auto
    if (auto_simplify)
    {
        simplify();
        end_time = std::chrono::high_resolution_clock::now();
    }

    if (gen_test)
        generate_test_output();
}

long long Slow_simplifier_PSLG::get_PITC() const
{
    return point_in_triangle_checks;
}


// .in file contains lines with 2 integers with an empty line between chains
std::vector<std::vector<Slow_simplifier_PSLG::Point>>
Slow_simplifier_PSLG::parse_input_file_as_chains(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
    {
        throw std::runtime_error("Failed to open input file: " + path);
    }

    vector<vector<Point>> chains_points;
    vector<Point> current_chain;
    string line;

    auto flush_current_chain = [&]()
    {
        if (!current_chain.empty())
        {
            chains_points.push_back(current_chain);
            current_chain.clear();
        }
    };

    while (std::getline(in, line))
    {
        // trim-as-blank check
        bool all_ws = true;
        for (char c : line)
        {
            if (!std::isspace((unsigned char)c))
            {
                all_ws = false;
                break;
            }
        }

        if (all_ws)
        {
            // chain separator
            flush_current_chain();
            continue;
        }

        std::istringstream iss(line);
        std::vector<K::FT> nums;
        K::FT v;
        while (iss >> v) nums.push_back(v);

        if (nums.size() >= 2)
        {
            for (size_t i = 0; i + 1 < nums.size(); i += 2)
            {
                current_chain.emplace_back(nums[i], nums[i + 1]);
            }
        }
        // silently ignore odd trailing number in a line, malformed input
    }
    flush_current_chain();

    return chains_points;
}

// build internal structures (merge duplicates, create node occurrences)
void Slow_simplifier_PSLG::build_internal_structures_from_chains(const std::vector<std::vector<Point>>& chains_points)
{
    // prepare member chains (chain sizes)
    chains.clear();
    chains.resize(chains_points.size());

    points.clear();
    PI.clear();
    node_to_global_vid.clear();
    global_coord_to_vid.clear();
    global_vid_counts.clear();
    chain_pos.clear();
    chain_closed.clear();

    int node_index = 0;
    for (size_t cid = 0; cid < chains_points.size(); ++cid)
    {
        const auto& chain_pts = chains_points[cid];
        chains[cid].reserve(chain_pts.size());
        bool closed = false;
        if (chain_pts.size() >= 2)
        {
            // closed if first equals last
            if (chain_pts.front() == chain_pts.back())
                closed = true;
        }
        chain_closed.push_back(closed ? 1 : 0);

        // if chain is closed by repeating first point at end, skip final duplicate occurrence.
        size_t limit = chain_pts.size();
        if (closed && limit > 0) limit = limit - 1;

        for (size_t pi_idx = 0; pi_idx < limit; ++pi_idx)
        {
            auto p = chain_pts[pi_idx];

            // merge duplicates by coordinate
            pair<K::FT, K::FT> coord = {p.x(), p.y()};
            auto it = global_coord_to_vid.find(coord);
            int gid;
            if (it == global_coord_to_vid.end())
            {
                gid = static_cast<int>(global_coord_to_vid.size());
                global_coord_to_vid[coord] = gid;
                global_vid_counts.push_back(0);
                // ensure chain-count arrays grow in parallel if they exist
                if (global_vid_chain_count.size() < global_coord_to_vid.size())
                {
                    global_vid_chain_count.push_back(0);
                    global_vid_chain_last_seen.push_back(-1);
                }
            }
            else
            {
                gid = it->second;
            }

            // if the chain-count arrays exist, update them (count distinct chains)
            if (!global_vid_chain_last_seen.empty())
            {
                if (global_vid_chain_last_seen[gid] != static_cast<int>(cid))
                {
                    global_vid_chain_count[gid] += 1;
                    global_vid_chain_last_seen[gid] = static_cast<int>(cid);
                }
            }

            // append node entry to points list
            NodeEntry ne;
            ne.p = p;
            ne.meta.first = node_index;
            ne.meta.second = 0; // block cursor initialised
            ne.chain_id = static_cast<int>(cid);
            points.push_back(ne);
            auto itlist = points.end();
            --itlist;
            PI.push_back(itlist);

            // bookkeeping
            node_to_global_vid.push_back(gid);
            global_vid_counts[gid] += 1;
            chains[cid].push_back(node_index);
            chain_pos.push_back(static_cast<int>(chains[cid].size()) - 1);

            ++node_index;
        }
    }

    // now we have:
    // - points (list of NodeEntry occurrences)
    // - PI mapping node_index -> iterator into points list
    // - chains vector: per-chain node indices in order
    // - chain_pos mapping node_index -> position in its chain
    // - node_to_global_vid mapping node_index -> global coord id
    // - global_vid_counts mapping global coord id -> number occurrences across all chains

    // init some members used by simplify
    init_vertex_count = static_cast<int>(PI.size());
}


// return (prev_node_index, next_node_index) in the same chain as node index
std::pair<int, int> Slow_simplifier_PSLG::get_neighbours(const int ind) const
{
    // TODO remove
    {
        if (ind < 0 || ind >= (int)PI.size())
        {
            throw runtime_error(
                "get_neighbours: invalid node index: " + std::to_string(ind) + " (PI.size=" + std::to_string(PI.size())
                +
                ")");
        }
        if (chain_pos.size() != PI.size())
        {
            // defensive: chain_pos should align with PI
            throw runtime_error(
                "get_neighbours: chain_pos.size() mismatch (chain_pos.size=" + std::to_string(chain_pos.size()) +
                ", PI.size=" + std::to_string(PI.size()) + ")");
        }
    }

    // we avoid directly dereferencing PI[ind] here beyond reading chain_id; if the iterator is invalid this will still throw
    const int cid = PI[ind]->chain_id;
    // TODO remove
    {
        if (cid < 0 || cid >= (int)chains.size())
        {
            throw runtime_error(
                "get_neighbours: invalid chain id for node " + std::to_string(ind) + ": " + std::to_string(cid));
        }
    }

    const int pos = chain_pos[ind];
    const auto& chain_nodes = chains[cid];

    int prev = -1;
    int next = -1;

    if (chain_nodes.empty()) return {-1, -1};

    // previous
    if (pos > 0)
        prev = chain_nodes[pos - 1];
    else
    {
        if (chain_closed[cid])
            prev = chain_nodes.back();
        else
            prev = -1;
    }

    // next
    if (pos + 1 < (int)chain_nodes.size())
        next = chain_nodes[pos + 1];
    else
    {
        if (chain_closed[cid])
            next = chain_nodes.front();
        else
            next = -1;
    }

    return {prev, next};
}

bool Slow_simplifier_PSLG::is_in_triangle(Point p, Point tr1, Point tr2, Point tr3) const
{
    // handle degenerate case: collinear points (never blocked)
    if (CGAL::collinear(tr1, tr2, tr3))
        return false;

    CGAL::Orientation ori1 = CGAL::orientation(tr1, tr2, p);
    CGAL::Orientation ori2 = CGAL::orientation(tr2, tr3, p);
    CGAL::Orientation ori3 = CGAL::orientation(tr3, tr1, p);

    return !(ori1 == CGAL::RIGHT_TURN || ori2 == CGAL::RIGHT_TURN || ori3 == CGAL::RIGHT_TURN);
}

K::FT Slow_simplifier_PSLG::get_area(int ind)
{
    // TODO remove
    {
        if (ind < 0 || ind >= (int)PI.size())
        {
            throw runtime_error("Index for area invalid" + std::to_string(ind));
        }
    }

    auto [nb1, nb2] = get_neighbours(ind);
    if (nb1 == -1 || nb2 == -1)
        return K::FT(0); // endpoints / undefined triangle => area 0 (not removable)
    Point p1 = get_pi(ind)->p;
    Point p2 = nb1 >= 0 ? get_pi(nb1)->p : Point(0, 0);
    Point p3 = nb2 >= 0 ? get_pi(nb2)->p : Point(0, 0);

    return abs(K::Triangle_2(p1, p2, p3).area());
}

void Slow_simplifier_PSLG::handle_point(Point_iterator pi,
                                        std::map<std::pair<K::FT, int>, Point>& ordered_triangles,
                                        const std::vector<char>& removed,
                                        std::map<std::pair<K::FT, int>, Point>::iterator& mi)
{
    // TODO remove
    {
        if (pi == points.end())
        {
            // Should never happen if caller only uses valid node indices
            return;
        }
    }

    auto& node = *pi;
    int ind = node.meta.first;
    // TODO remove
    {
        if (ind < 0 || ind >= (int)PI.size()) return;
    }
    int& block = node.meta.second;

    auto [nb1, nb2] = get_neighbours(ind);

    // if missing neighbors, it's an endpoint on an open chain => not removable
    if (nb1 == -1 || nb2 == -1)
        return;

    // swap points 2 and 3 to get ccw triangle if needed
    if (CGAL::orientation(pi->p, get_pi(nb1)->p, get_pi(nb2)->p) != CGAL::LEFT_TURN)
        std::swap(nb1, nb2);

    for (; block < init_vertex_count; block++)
    {
        if (removed[block])
            continue;

        if (block == ind || block == nb1 || block == nb2)
            continue;

        Point oth;
        // TODO remove
        {
            // if the iterator for 'block' was invalid/dangling, this could throw: wrap in try-catch to convert to informative message
            try
            {
                oth = get_pi(block)->p;
            }
            catch (const std::exception& e)
            {
                throw runtime_error(
                    std::string("handle_point: get_pi(block) threw for block=") + std::to_string(block) + "; what(): " +
                    e.
                    what());
            }
        }
        ++point_in_triangle_checks;
        if (is_in_triangle(oth, pi->p, get_pi(nb1)->p, get_pi(nb2)->p))
            break;
    }

    if (block == init_vertex_count)
    {
        // insert into ordered_triangles and record iterator in mi
        auto key = std::make_pair(get_area(ind), ind);
        auto itpair = ordered_triangles.insert(std::make_pair(key, pi->p));
        mi = itpair.first;
    }
}

void Slow_simplifier_PSLG::handle_neighbour(Point_iterator pi,
                                            std::map<std::pair<K::FT, int>, Point>& ordered_triangles,
                                            std::map<std::pair<K::FT, int>, Point>::iterator& mi)
{
    // remove from map if present
    if (mi != ordered_triangles.end())
        ordered_triangles.erase(mi);
    mi = ordered_triangles.end();

    // TODO: REMOVE
    {
        if (pi == points.end())
        {
            // odd: caller passed end() or iterator invalid; do not dereference
            return;
        }
    }
    pi->meta.second = 0;
}

// helper to check whether a node occurrence is a candidate for removal:
// must have exactly 1 distinct chain (i.e. not a junction)
bool Slow_simplifier_PSLG::is_node_candidate_removable(int node_index) const
{
    if (node_index < 0 || node_index >= (int)node_to_global_vid.size())
    {
        // TODO: remove this potentially, should never be possible
        throw runtime_error("Invalid candidate index");
        return false;
    }
    int gid = node_to_global_vid[node_index];

    if (gid < 0 || gid >= (int)global_vid_chain_count.size()) return false;
    if (global_vid_chain_count[gid] != 1) return false; // junction (multiple distinct chains) or isolated weirdness


    // must have both neighbours
    auto [nb1, nb2] = get_neighbours(node_index);
    if (nb1 == -1 || nb2 == -1) return false; // open endpoint
    // otherwise candidate
    return true;
}

void Slow_simplifier_PSLG::simplify(int remaining_vertices)
{
    int vertices_left = init_vertex_count - (int)result.size();
    std::cout << "SIMPLIFYING FROM " << vertices_left << " TO " << remaining_vertices << std::endl;

    // sorts triangles by area, and then by vertex index (in order to make the algo predictable)
    // this also ensures unique keys
    std::map<std::pair<K::FT, int>, Point> ordered_triangles;

    // keeps track of vertices that were removed (char is used since vector<bool> is bad practice)
    // seed the removed mask from 'result' so repeated simplify() calls are safe
    std::vector<char> removed(init_vertex_count);
    if (!result.empty())
    {
        for (int idx : result)
        {
            if (idx >= 0 && idx < init_vertex_count)
                removed[idx] = 1;
        }
    }

    // maps each vertex index to the corresponding iterator in the map
    std::vector<decltype(ordered_triangles.begin())> index_to_MI(init_vertex_count);
    for (auto& it : index_to_MI)
        it = ordered_triangles.end();

    // initialize candidate nodes (not junctions or endpoints)
    {
        for (int i = 0; i < init_vertex_count; ++i)
        {
            if (removed[i]) continue; // skip previously removed nodes

            bool candidate = false;
            try
            {
                candidate = is_node_candidate_removable(i);
            }
            catch (const std::exception& e)
            {
                throw runtime_error(
                    std::string("simplify: is_node_candidate_removable threw for i=") + std::to_string(i) + "; what(): "
                    + e.what());
            }

            if (!candidate)
            {
                cout << "  node " << i << " not candidate" << endl;
                continue;
            }

            auto pi = get_pi(i);
            handle_point(pi, ordered_triangles, removed, index_to_MI[i]);
        }
    }

    cout << "DEBUG: ordered_triangles size after init = " << ordered_triangles.size() << endl;

    int bound = init_vertex_count - remaining_vertices;
    for (int it_count = 0; it_count < bound; it_count++)
    {
        if (ordered_triangles.empty())
        {
            // no more removable candidates (likely because remaining graph contains only endpoints/junctions or degenerate triangles).
            // we stop early (cannot remove more without breaking topology).
            break;
        }

        // get vertex handle of next vertex that is removed (smallest area non-blocked)
        auto best_it = ordered_triangles.begin();
        Point best_point = best_it->second;
        int index = best_it->first.second;
        ordered_triangles.erase(best_it);
        // mark its map-iterator slot as end (it is no longer in the map)
        if (index >= 0 && index < (int)index_to_MI.size())
        {
            index_to_MI[index] = ordered_triangles.end();
        }

        // add index of vertex to result
        result.push_back(index);

        // handle neighbor vertices: remove them from ordered set (they will be re-evaluated)
        auto [index_nb1, index_nb2] = get_neighbours(index);
        if (index_nb1 != -1)
            handle_neighbour(get_pi(index_nb1), ordered_triangles, index_to_MI[index_nb1]);
        if (index_nb2 != -1)
            handle_neighbour(get_pi(index_nb2), ordered_triangles, index_to_MI[index_nb2]);

        // remove vertex occurrence from chain and points list
        // erase from chain vector and update chain_pos for subsequent nodes
        int cid = get_pi(index)->chain_id;
        int pos = chain_pos[index];
        auto& vec = chains[cid];
        vec.erase(vec.begin() + pos);

        // update chain_pos for nodes after pos
        for (int j = pos; j < (int)vec.size(); ++j)
        {
            chain_pos[vec[j]] = j;
        }

        // erase from points list
        points.erase(get_pi(index));
        // mark removed and keep PI[index] as dangling (we won't use it again)
        removed[index] = 1;

        // keep global_vid_counts consistent so candidacy tests remain correct
        {
            int gid = node_to_global_vid[index];
            if (gid >= 0 && gid < (int)global_vid_counts.size())
                global_vid_counts[gid] -= 1;
        }

        // unblock points: re-evaluate all nodes that were candidates but currently not in ordered_triangles
        for (int j = 0; j < init_vertex_count; j++)
        {
            if (removed[j]) continue;
            if (!is_node_candidate_removable(j)) continue; // not a candidate (endpoints / junctions)
            if (index_to_MI[j] != ordered_triangles.end()) continue; // already in ordered_triangles

            handle_point(get_pi(j), ordered_triangles, removed, index_to_MI[j]);
        }
    }
}


void Slow_simplifier_PSLG::generate_test_output()
{
    std::ofstream fout("../data/" + name + "/data.out");
    for (auto index : result)
        fout << index << std::endl;
}

void Slow_simplifier_PSLG::chain_to_ipe(bool original)
{
    // build ipe chains from current chains (skip empty chains)
    std::vector<IPE::Chain> out_chains;
    out_chains.reserve(chains.size());

    for (size_t cid = 0; cid < chains.size(); ++cid)
    {
        const auto &chain_nodes = chains[cid];
        if (chain_nodes.empty()) continue;

        IPE::Chain chain;
        chain.reserve(chain_nodes.size());

        for (int node_idx : chain_nodes)
        {
            const Point &pp = get_pi(node_idx)->p;
            chain.emplace_back(CGAL::to_double(pp.x()), CGAL::to_double(pp.y()));
        }

        out_chains.push_back(std::move(chain));
    }

    if (out_chains.empty()) return;

    // write all chains into one ipe file
    IPE::chains_to_IPE(name, out_chains, original);
}


// generate ipe output for multiple target sizes (descending), simplifying before each write
void Slow_simplifier_PSLG::create_ipe_chains(std::vector<int> save_sizes)
{
    std::sort(save_sizes.begin(), save_sizes.end(), std::greater<int>());

    for (int vertex_count : save_sizes)
    {
        // compute current remaining nodes (sum of chain lengths)
        int current_remaining = 0;
        for (const auto &ch : chains) current_remaining += static_cast<int>(ch.size());

        if (vertex_count > current_remaining)
            continue;

        // simplify down to the requested total remaining count (across all chains)
        if (vertex_count < current_remaining)
            simplify(vertex_count);

        // write all chains (in same file) to ipe
        chain_to_ipe();
    }
}

int Slow_simplifier_PSLG::get_initial_vertex_count() const
{
    return init_vertex_count;
}
