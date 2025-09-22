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
    start_time = std::chrono::steady_clock::now();

    // parse input file into chains (flexible parsing; blank lines separate chains; 2 numbers per line are points)
    const auto chains_points = parse_input_file_as_chains("../data/" + input_folder_name + "/data.in");

    // build internal node list, merge duplicate coordinates and create per-node structures
    build_internal_structures_from_chains(chains_points);

    // run simplification if auto
    if (auto_simplify)
    {
        simplify();
        end_time = std::chrono::steady_clock::now();
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

    global_removed.clear();
    points.clear();
    PI.clear();
    node_to_gid.clear();
    global_coord_to_vid.clear();
    chain_pos.clear();
    chain_closed.clear();
    gid_to_neighbors.clear();

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
            const auto& p = chain_pts[pi_idx];

            // merge duplicates by coordinate
            pair<K::FT, K::FT> coord = {p.x(), p.y()};
            auto it = global_coord_to_vid.find(coord);
            int gid;
            if (it == global_coord_to_vid.end())
            {
                gid = static_cast<int>(global_coord_to_vid.size());
                global_coord_to_vid[coord] = gid;
            }
            else
            {
                gid = it->second;
            }

            // append node entry to points list
            NodeEntry ne;
            ne.p = p;
            ne.chain_id = static_cast<int>(cid);
            points.push_back(ne);
            auto itlist = points.end();
            --itlist;
            PI.push_back(itlist);

            // bookkeeping
            node_to_gid.push_back(gid);
            chains[cid].push_back(node_index);
            chain_pos.push_back(static_cast<int>(chains[cid].size()) - 1);

            ++node_index;
        }
    }

    // build global neighbor sets for junction detection (initial state)
    for (size_t cid = 0; cid < chains.size(); ++cid)
    {
        const auto& chain_nodes = chains[cid];
        const size_t chain_len = chain_nodes.size();
        if (chain_len == 0) continue;

        for (size_t i = 0; i < chain_len; ++i)
        {
            int node_i = chain_nodes[i];
            int node_prev = (i == 0) ? (chain_closed[cid] ? chain_nodes.back() : -1) : chain_nodes[i - 1];
            int node_next = (i == chain_len - 1) ? (chain_closed[cid] ? chain_nodes.front() : -1) : chain_nodes[i + 1];

            int gid_i = node_to_gid[node_i];
            if (node_prev != -1)
            {
                int gid_prev = node_to_gid[node_prev];
                gid_to_neighbors[gid_i].insert(gid_prev);
            }
            if (node_next != -1)
            {
                int gid_next = node_to_gid[node_next];
                gid_to_neighbors[gid_i].insert(gid_next);
            }
        }
    }

    // now we have:
    // - points (list of NodeEntry occurrences)
    // - PI mapping node_index -> iterator into points list
    // - chains vector: per-chain node indices in order
    // - chain_pos mapping node_index -> position in its chain
    // - node_to_global_vid mapping node_index -> global coord id
    // - global_vid_to_original_neighbors mapping global coord id -> set of neighboring global coord ids (not updated during simplification)

    init_global_vertex_count = static_cast<int>(global_coord_to_vid.size());

    // setup removed vector for global vertices
    global_removed.resize(init_global_vertex_count);

    // build reverse mapping global_vid -> occurrences and canonical point per global vid
    gid_to_nodes.clear();
    gid_to_nodes.resize(init_global_vertex_count);
    gid_to_point.clear();
    gid_to_point.resize(init_global_vertex_count);

    for (int node = 0; node < (int)node_to_gid.size(); ++node)
    {
        int gid = node_to_gid[node];
        gid_to_nodes[gid].push_back(node);
        // set canonical point if not yet set (we can overwrite repeatedly, it's the same coordinate)
        gid_to_point[gid] = PI[node]->p;
    }
}


// return (prev_node_index, next_node_index) in the same chain as node index
std::pair<int, int> Slow_simplifier_PSLG::get_neighbours(const int ind) const
{
    // defensive checks
    if (ind < 0 || ind >= (int)PI.size())
    {
        throw runtime_error(
            "get_neighbours: invalid node index: " + std::to_string(ind) + " (PI.size=" + std::to_string(PI.size())
            +
            ")");
    }
    if (chain_pos.size() != PI.size())
    {
        throw runtime_error(
            "get_neighbours: chain_pos.size() mismatch (chain_pos.size=" + std::to_string(chain_pos.size()) +
            ", PI.size=" + std::to_string(PI.size()) + ")");
    }

    // if this occurrence was erased, PI[ind] will be points.end()
    if (PI[ind] == points.end())
        return {-1, -1};

    const int cid = PI[ind]->chain_id;
    if (cid < 0 || cid >= (int)chains.size())
    {
        throw runtime_error(
            "get_neighbours: invalid chain id for node " + std::to_string(ind) + ": " + std::to_string(cid));
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

bool Slow_simplifier_PSLG::is_in_triangle(const Point& p, const Point& tr1, const Point& tr2, const Point& tr3)
{
    // handle degenerate case: collinear points (never blocked)
    if (CGAL::collinear(tr1, tr2, tr3))
        return false;

    const CGAL::Orientation ori1 = CGAL::orientation(tr1, tr2, p);
    const CGAL::Orientation ori2 = CGAL::orientation(tr2, tr3, p);
    const CGAL::Orientation ori3 = CGAL::orientation(tr3, tr1, p);

    return !(ori1 == CGAL::RIGHT_TURN || ori2 == CGAL::RIGHT_TURN || ori3 == CGAL::RIGHT_TURN);
}

K::FT Slow_simplifier_PSLG::get_area(int ind)
{
    auto [nb1, nb2] = get_neighbours(ind);
    if (nb1 == -1 || nb2 == -1)
        return K::FT(0); // endpoints / undefined triangle => area 0 (not removable)
    Point p1 = get_pi(ind)->p;
    Point p2 = nb1 >= 0 ? get_pi(nb1)->p : Point(0, 0);
    Point p3 = nb2 >= 0 ? get_pi(nb2)->p : Point(0, 0);

    return abs(K::Triangle_2(p1, p2, p3).area());
}

void Slow_simplifier_PSLG::handle_point(
    int gid,
    std::map<std::pair<K::FT, int>, Point>& ordered_triangles,
    std::map<std::pair<K::FT, int>, Point>::iterator& mi)
{
    if (gid < 0 || gid >= (int)gid_to_point.size())
        throw runtime_error("handle_point: invalid gid");

    if (global_removed[gid])
        return;

    // quick candidate check
    if (!is_node_candidate_removable(gid))
        return;

    // get the two current neighbour gids (must exist because is_node_candidate_removable passed)
    const auto& ngset = gid_to_neighbors.at(gid);
    auto sit = ngset.begin();
    int gid_nb1 = *sit;
    ++sit;
    int gid_nb2 = *sit;

    // canonical coordinates
    Point p = gid_to_point[gid];
    Point p1 = gid_to_point[gid_nb1];
    Point p2 = gid_to_point[gid_nb2];

    // ensure CCW (triangle order must be CCW)
    if (CGAL::orientation(p, p1, p2) != CGAL::LEFT_TURN)
    {
        std::swap(p1, p2);
        std::swap(gid_nb1, gid_nb2);
    }

    // scan for blocking vertices
    const int global_count = static_cast<int>(gid_to_point.size());
    bool blocked = false;

    for (auto block = 0; block < global_count; ++block)
    {
        if (block == gid) continue;
        if (global_removed[block]) continue;
        if (block == gid_nb1 || block == gid_nb2) continue;

        Point oth = gid_to_point[block];
        ++point_in_triangle_checks;

        if (is_in_triangle(oth, p, p1, p2))
        {
            blocked = true;
            break;
        }
    }

    if (!blocked)
    {
        K::FT area = abs(K::Triangle_2(p, p1, p2).area());
        auto key = std::make_pair(area, gid);
        auto itpair = ordered_triangles.insert({key, p});
        mi = itpair.first;
    }
}


// remove global id from ordered_triangles
void Slow_simplifier_PSLG::handle_neighbour_global(int gid, std::map<std::pair<K::FT, int>, Point>& ordered_triangles,
                                                   std::vector<std::map<std::pair<K::FT, int>, Point>::iterator>&
                                                   index_to_MI)
{
    if (gid < 0 || gid >= (int)index_to_MI.size()) return;
    auto& mi = index_to_MI[gid];
    if (mi != ordered_triangles.end())
        ordered_triangles.erase(mi);
    mi = ordered_triangles.end();
}

// helper to check whether a node occurrence is a candidate for removal:
// must have exactly 1 distinct chain (i.e. not a junction)
[[nodiscard]] bool Slow_simplifier_PSLG::is_node_candidate_removable(int gid) const
{
    if (gid >= (int)gid_to_point.size())
        throw runtime_error("Invalid candidate gid");

    // junction/open end detection: candidates must have exactly 2 neighbors globally
    auto it = gid_to_neighbors.find(gid);
    if (it == gid_to_neighbors.end() || it->second.size() != 2)
        return false;

    // disallow removing third-last node from chain when it is closed
    for (int node_index : gid_to_nodes[gid])
    {
        const int cid = PI[node_index]->chain_id;
        if (chain_closed[cid] && static_cast<int>(chains[cid].size()) <= 3)
            return false;
    }

    // otherwise candidate
    return true;
}


unsigned long long Slow_simplifier_PSLG::get_vertices_left() const
{
    return init_global_vertex_count - result.size();
}

void Slow_simplifier_PSLG::simplify(const int remaining_vertices)
{
    // compute how many vertices currently left as sum of chain sizes
    auto vertices_left = get_vertices_left();
    std::cout << "SIMPLIFYING FROM " << vertices_left << " TO " << remaining_vertices << std::endl;
    if (remaining_vertices >= vertices_left)
        return; // nothing to do

    // sorts triangles by area, and then by vertex index (in order to make the algo predictable)
    // this also ensures unique keys
    std::map<std::pair<K::FT, int>, Point> ordered_triangles;

    // keeps track of vertices that were removed (char is used since vector<bool> is bad practice)
    // seed the removed mask from 'result' so repeated simplify() calls are safe
    if (!result.empty())
    {
        for (int v : result)
        {
            if (v >= 0 && v < (int)global_removed.size())
            {
                // v is a gid
                global_removed[v] = 1;
            }
            else
            {
                throw runtime_error("Out of range node index in result: " + std::to_string(v));
            }
        }
    }

    // maps each global vertex id to the corresponding iterator in the map
    const int global_count = static_cast<int>(gid_to_point.size());
    std::vector<decltype(ordered_triangles.begin())> index_to_MI(global_count);
    for (auto& it : index_to_MI)
        it = ordered_triangles.end();

    // initialize member gid_to_neighbors from current chains/occurrences
    gid_to_neighbors.clear();
    for (int gid = 0; gid < global_count; ++gid)
    {
        if (global_removed[gid]) continue;
        for (int node_idx : gid_to_nodes[gid])
        {
            if (node_idx < 0 || node_idx >= (int)PI.size()) continue;
            int cid = PI[node_idx]->chain_id;
            int pos = chain_pos[node_idx];
            if (!(pos >= 0 && cid >= 0 && cid < (int)chains.size())) continue;
            const auto& vec = chains[cid];
            if (!(pos < (int)vec.size() && vec[pos] == node_idx)) continue; // removed occurrence
            auto [nb1, nb2] = get_neighbours(node_idx);
            if (nb1 != -1) gid_to_neighbors[gid].insert(node_to_gid[nb1]);
            if (nb2 != -1) gid_to_neighbors[gid].insert(node_to_gid[nb2]);
        }
    }

    // initialize candidate global vertices (not junctions or endpoints) using member gid_to_neighbors
    for (int gid = 0; gid < global_count; ++gid)
    {
        if (global_removed[gid]) continue; // already removed
        auto it = gid_to_neighbors.find(gid);
        if (it == gid_to_neighbors.end() || it->second.size() != 2)
        {
            continue;
        }

        auto mi = index_to_MI[gid];
        handle_point(gid, ordered_triangles, mi);
        index_to_MI[gid] = mi;
    }

    auto current_remaining = vertices_left;
    auto bound = current_remaining - remaining_vertices;
    for (int it_count = 0; it_count < bound; it_count++)
    {
        if (ordered_triangles.empty())
        {
            // no more removable candidates (likely because remaining graph contains only endpoints/junctions or degenerate triangles).
            // we stop early (cannot remove more without breaking topology).
            break;
        }

        // get global vertex id of next vertex that is removed (smallest area non-blocked)
        auto best_it = ordered_triangles.begin();

        // Check PSLG conditions
        if (!is_node_candidate_removable(best_it->first.second))
        {
            // no longer a candidate remove from map and continue
            int gid = best_it->first.second;
            ordered_triangles.erase(best_it);
            index_to_MI[gid] = ordered_triangles.end();
            it_count--;
            continue;
        }

        Point best_point = best_it->second;
        int gid = best_it->first.second;
        ordered_triangles.erase(best_it);
        // mark its map-iterator slot as end (it is no longer in the map)
        if (gid >= 0 && gid < (int)index_to_MI.size())
        {
            index_to_MI[gid] = ordered_triangles.end();
        }

        // add global id of vertex to result
        result.push_back(gid);

        // snapshot the current neighbours (before we mutate data)
        std::set<int> old_neighbours;
        auto itg = gid_to_neighbors.find(gid);
        if (itg != gid_to_neighbors.end())
            old_neighbours = itg->second;

        // remove the gid from ordered_triangles entries of its neighbours (they must be re-evaluated)
        for (int neigh_gid : old_neighbours)
        {
            if (neigh_gid < 0 || neigh_gid >= (int)index_to_MI.size()) continue;
            handle_neighbour_global(neigh_gid, ordered_triangles, index_to_MI);
        }

        // remove all occurrences of this global vertex from chains and points list
        // TODO check if necessary
        auto occurrences = gid_to_nodes[gid]; // copy
        std::sort(occurrences.begin(), occurrences.end(), [&](int a, int b)
        {
            int pa = (a >= 0 && a < (int)chain_pos.size()) ? chain_pos[a] : -1;
            int pb = (b >= 0 && b < (int)chain_pos.size()) ? chain_pos[b] : -1;
            if (pa != pb) return pa > pb;
            return a > b;
        });

        for (int node_idx : occurrences)
        {
            if (node_idx < 0 || node_idx >= (int)PI.size()) continue;
            if (PI[node_idx] == points.end()) continue; // already erased

            int cid = PI[node_idx]->chain_id;
            int pos = chain_pos[node_idx];
            auto& vec = chains[cid];

            if (pos >= 0 && pos < (int)vec.size() && vec[pos] == node_idx)
            {
                vec.erase(vec.begin() + pos);
                // update chain_pos for nodes after pos
                for (int j = pos; j < (int)vec.size(); ++j)
                {
                    chain_pos[vec[j]] = j;
                }
            }

            // erase from points list and mark occurrence as erased
            points.erase(PI[node_idx]);
            PI[node_idx] = points.end();
            chain_pos[node_idx] = -1;
        }


        // mark removed globally
        global_removed[gid] = 1;

        // clear out the removed gid's neighbour set and occurrences (tidy)
        gid_to_neighbors[gid].clear();
        gid_to_nodes[gid].clear();

        // Now recompute neighbour sets for the gids that might have changed (old_neighbours),
        // because removing gid changes adjacency for those vertices.
        for (int neigh_gid : old_neighbours)
        {
            if (neigh_gid < 0 || neigh_gid >= (int)gid_to_point.size()) continue;
            if (global_removed[neigh_gid])
            {
                gid_to_neighbors[neigh_gid].clear();
                continue;
            }
            std::set<int> recomputed;
            for (int node_idx : gid_to_nodes[neigh_gid])
            {
                if (node_idx < 0 || node_idx >= (int)PI.size()) continue;
                int cid = PI[node_idx]->chain_id;
                int pos = chain_pos[node_idx];
                if (!(pos >= 0 && cid >= 0 && cid < (int)chains.size() && chains[cid][pos] == node_idx))
                    continue; // occurrence removed or stale
                auto [nb1_node, nb2_node] = get_neighbours(node_idx);
                if (nb1_node != -1) recomputed.insert(node_to_gid[nb1_node]);
                if (nb2_node != -1) recomputed.insert(node_to_gid[nb2_node]);
            }
            gid_to_neighbors[neigh_gid] = std::move(recomputed);
        }

        // reset its iterator slot for the removed gid
        if (gid >= 0 && gid < (int)index_to_MI.size())
            index_to_MI[gid] = ordered_triangles.end();

        // unblock points: re-evaluate affected neighbours (they were invalidated earlier)
        for (int gid_check : old_neighbours)
        {
            if (gid_check < 0 || gid_check >= (int)gid_to_point.size()) continue;
            if (global_removed[gid_check]) continue; // already removed
            auto itn = gid_to_neighbors.find(gid_check);
            if (itn == gid_to_neighbors.end() || itn->second.size() != 2) continue; // not a candidate
            if (index_to_MI[gid_check] != ordered_triangles.end()) continue; // already in ordered_triangles
            auto mi = index_to_MI[gid_check];
            handle_point(gid_check, ordered_triangles, mi);
            index_to_MI[gid_check] = mi;
        }
    }
}


void Slow_simplifier_PSLG::generate_test_output()
{
    // TODO: might be bugged
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
        const auto& chain_nodes = chains[cid];
        if (chain_nodes.empty()) continue;

        IPE::Chain chain;
        chain.reserve(chain_nodes.size());

        for (int node_idx : chain_nodes)
        {
            const Point& pp = get_pi(node_idx)->p;
            chain.emplace_back(CGAL::to_double(pp.x()), CGAL::to_double(pp.y()));
        }

        out_chains.push_back(std::move(chain));
    }

    if (out_chains.empty())
    {
        std::cerr << "No IPE chains found" << std::endl;
        return;
    }

    // write all chains into one ipe file
    IPE::chains_to_IPE(name, out_chains, chain_closed, get_vertices_left(), original);
}


// generate ipe output for multiple target sizes (descending), simplifying before each write
void Slow_simplifier_PSLG::create_ipe_chains(std::vector<int> save_sizes)
{
    std::sort(save_sizes.begin(), save_sizes.end(), std::greater<int>());

    chain_to_ipe(true);
    std::cout << "WROTE ORIGINAL" << std::endl;

    auto last_save_size = -1;
    for (int vertex_count : save_sizes)
    {
        int current_remaining = get_vertices_left();
        if (current_remaining == last_save_size)
        {
            std::cout << "Stopped simplifcation as no further simplification possible" << std::endl;
            break; // stop if we cannot simplify further
        }

        if (vertex_count > current_remaining)
            continue;

        // simplify down to the requested total remaining count (across all chains)
        if (vertex_count < current_remaining)
            simplify(vertex_count);

        // write all chains (in same file) to ipe
        chain_to_ipe();

        last_save_size = current_remaining;
    }
}

int Slow_simplifier_PSLG::get_initial_vertex_count() const
{
    return init_global_vertex_count;
}
