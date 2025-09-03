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

    global_removed.clear();
    points.clear();
    PI.clear();
    node_to_gid.clear();
    global_coord_to_vid.clear();
    chain_pos.clear();
    chain_closed.clear();
    gid_to_original_neighbors.clear();

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

    // build global neighbor sets for junction detection
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
                gid_to_original_neighbors[gid_i].insert(gid_prev);
            }
            if (node_next != -1)
            {
                int gid_next = node_to_gid[node_next];
                gid_to_original_neighbors[gid_i].insert(gid_next);
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

    // init some members used by simplify
    init_vertex_count = static_cast<int>(PI.size());

    // setup removed vector for global vertices
    global_removed.resize(global_coord_to_vid.size());

    // build reverse mapping global_vid -> occurrences and canonical point per global vid
    gid_to_nodes.clear();
    gid_to_nodes.resize(global_coord_to_vid.size());
    gid_to_point.clear();
    gid_to_point.resize(global_coord_to_vid.size());

    for (int node = 0; node < (int)node_to_gid.size(); ++node)
    {
        int gid = node_to_gid[node];
        gid_to_nodes[gid].push_back(node);
        // set canonical point if not yet set (we can overwrite repeatedly, it's the same coordinate)
        gid_to_point[gid] = PI[node]->p;
    }

    // block cursors per global id
    global_block_cursor.clear();
    global_block_cursor.resize(global_coord_to_vid.size(), 0);
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


    // print all oris
    std::cout << "ORIS: " << ori1 << " " << ori2 << " " << ori3 << std::endl;

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

void Slow_simplifier_PSLG::handle_point(int gid,
                                        std::map<std::pair<K::FT, int>, Point>& ordered_triangles,
                                        std::map<std::pair<K::FT, int>, Point>::iterator& mi)
{
    // TODO remove
    {
        if (gid < 0 || gid >= (int)gid_to_point.size())
        {
            throw runtime_error(
                "handle_point: invalid gid: " + std::to_string(gid) + " (gid_to_point.size=" +
                std::to_string(gid_to_point.size()) + ")");
        }
    }

    // TODO: remove
    {
        if (global_removed[gid])
        {
            throw runtime_error("handle_point: gid " + std::to_string(gid) + " already removed");
        }
    }

    // check candidate properties: not a junction and none of its occurrences are endpoints
    if (!is_node_candidate_removable(gid))
    {
        // TODO: remove cout, maybe remove this part altogether since I think we only call this when it is valid
        std::cout << "NOT REMOVAL IN HANDLE_POINT: " << gid << std::endl;
        return;
    }


    // For each occurrence compute triangle and check blocking. We accept the candidate if at least one occurrence
    // yields a valid (non-endpoint) triangle and is not blocked by any other global vertex.
    K::FT best_area = K::FT(0);
    bool any_ok = false;

    const int global_count = static_cast<int>(gid_to_point.size());

    for (int node_index : gid_to_nodes[gid])
    //    TODO: rework
    {
        // ensure occurrence still present (PI[node_index] might be dangling if removed previously)
        // attempt to detect removal by checking chain_pos validity: if chain_id invalid or pos out-of-range skip
        if (node_index < 0 || node_index >= (int)PI.size()) continue;
        try
        {
            auto & [p, chain_id] = *PI[node_index]; // may throw if iterator invalid

            // get occurrence-level neighbors (these are node indices)
            auto [nb1, nb2] = get_neighbours(node_index);

            // global neighbor ids (for comparing to 'block')
            int gid_nb1 = -1, gid_nb2 = -1;
            // neighbor Points (for triangle tests)
            Point pnb1, pnb2;

            // try to use occurrence neighbors first (both must exist)
            if (nb1 != -1 && nb2 != -1 && nb1 >= 0 &&
                nb1 < static_cast<int>(PI.size()) && nb2 >= 0 && nb2 < static_cast<int>(PI.size()))
            {
                gid_nb1 = node_to_gid[nb1];
                gid_nb2 = node_to_gid[nb2];
                pnb1 = get_pi(nb1)->p;
                pnb2 = get_pi(nb2)->p;
            }
            else
            {
                // fallback to global neighbors (from the set). Must be exactly two to be usable.
                auto it = gid_to_original_neighbors.find(gid);
                if (it == gid_to_original_neighbors.end() || it->second.size() != 2)
                    continue; // cannot reconstruct -> skip this occurrence

                auto sit = it->second.begin();
                gid_nb1 = *sit;
                ++sit;
                gid_nb2 = *sit;

                // use canonical points for those global neighbors
                pnb1 = gid_to_point[gid_nb1];
                pnb2 = gid_to_point[gid_nb2];
            }

            // ensure CCW: if not, swap both the points and the corresponding global ids
            if (CGAL::orientation(p, pnb1, pnb2) != CGAL::LEFT_TURN)
            {
                std::swap(pnb1, pnb2);
                std::swap(gid_nb1, gid_nb2);
            }

            // scan other global vertices (cursor stored in global_block_cursor[gid])
            int &block = global_block_cursor[gid];
            bool blocked = false;
            for (; block < global_count; ++block)
            {
                if (block == gid) continue;
                if (global_removed[block]) continue;

                // skip the immediate global neighbors (compare against gids)
                if (block == gid_nb1 || block == gid_nb2) continue;

                Point oth = gid_to_point[block];
                ++point_in_triangle_checks;

                if (is_in_triangle(oth, p, pnb1, pnb2))
                {
                    blocked = true;
                    ++block; // next time start after this blocking global
                    std::cout << "BLOCKED: " << gid << " by " << block - 1
                              << " with triangle (" << p << ", " << pnb1 << ", " << pnb2 << ")"
                              << " oth=" << oth << std::endl;
                    break;
                }
            }

            if (!blocked)
            {
                // compute area for this occurrence
                K::FT area = abs(K::Triangle_2(p, pnb1, pnb2).area());
                if (!any_ok || area < best_area)
                    best_area = area;
                any_ok = true;
            }
        }
        catch (const std::exception &e)
        {
            // occurrence iterator invalid/dangling -> skip this occurrence
            std ::cout << "handle_point: skipping occurrence " << node_index << " of gid " << gid
                      << " due to exception: " << e.what() << std::endl;
            continue;
        }
    }

    if (any_ok)
    {
        auto key = std::make_pair(best_area, gid);
        auto itpair = ordered_triangles.insert(std::make_pair(key, gid_to_point[gid]));
        mi = itpair.first;
    }
}



// remove global id from ordered_triangles and reset its block cursor
void Slow_simplifier_PSLG::handle_neighbour_global(int gid, std::map<std::pair<K::FT, int>, Point>& ordered_triangles, std::vector<std::map<std::pair<K::FT, int>, Point>::iterator>& index_to_MI)
{
    if (gid < 0 || gid >= (int)index_to_MI.size()) return;
    auto &mi = index_to_MI[gid];
    if (mi != ordered_triangles.end())
        ordered_triangles.erase(mi);
    mi = ordered_triangles.end();
    if (gid >= 0 && gid < (int)global_block_cursor.size())
        global_block_cursor[gid] = 0;
}

// helper to check whether a node occurrence is a candidate for removal:
// must have exactly 1 distinct chain (i.e. not a junction)
bool Slow_simplifier_PSLG::is_node_candidate_removable(int gid) const
{
    if (gid < 0 || gid >= (int)gid_to_point.size())
        throw runtime_error("Invalid candidate gid");

    // junction/open end detection: candidates must have exactly 2 neighbors globally
    auto it = gid_to_original_neighbors.find(gid);
    if (it == gid_to_original_neighbors.end() || it->second.size() != 2)
        return false;

    // ensure none of its occurrences are endpoints (we require every occurrence has both neighbours)
    // This is only necessary if we carea bout preserving original chain endpoints
    // for (int node_index : gid_to_nodes[gid])
    // {
    //     if (node_index < 0 || node_index >= (int)PI.size()) return false; // removed/dangling treated as not removable
    //     auto [nb1, nb2] = get_neighbours(node_index);
    //     if (nb1 == -1 || nb2 == -1) return false; // endpoint occurrence -> not removable globally
    // }

    // otherwise candidate
    return true;
}

void Slow_simplifier_PSLG::simplify(int remaining_vertices)
{
    // compute how many vertices currently left as sum of chain sizes
    int vertices_left = 0;
    for (const auto &ch : chains) vertices_left += static_cast<int>(ch.size());
    // Note: previously printed init-based count; preserve similar message but use current remaining
    std::cout << "SIMPLIFYING FROM " << vertices_left << " TO " << remaining_vertices << std::endl;

    // sorts triangles by area, and then by vertex index (in order to make the algo predictable)
    // this also ensures unique keys
    std::map<std::pair<K::FT, int>, Point> ordered_triangles;

    // keeps track of vertices that were removed (char is used since vector<bool> is bad practice)
    // seed the removed mask from 'result' so repeated simplify() calls are safe
    if (!result.empty())
    {
        for (int idx : result)
        {
            if (idx >= 0 && idx < init_vertex_count)
            {
                int gid = node_to_gid[idx];
                global_removed[gid] = 1;
            }
        }
    }

    // maps each global vertex id to the corresponding iterator in the map
    const int global_count = static_cast<int>(gid_to_point.size());
    std::vector<decltype(ordered_triangles.begin())> index_to_MI(global_count);
    for (auto& it : index_to_MI)
        it = ordered_triangles.end();

    // initialize candidate global vertices (not junctions or endpoints)
    for (int gid = 0; gid < global_count; ++gid)
    {
        if (global_removed[gid]) continue; // already removed

        bool candidate = false;
        try
        {
            candidate = is_node_candidate_removable(gid);
        }
        catch (const std::exception& e)
        {
            throw runtime_error(
                std::string("simplify: is_node_candidate_removable threw for gid=") + std::to_string(gid) + "; what(): "
                + e.what());
        }

        if (!candidate)
        {
            cout << "  global node " << gid << " not candidate" << endl;
            continue;
        }

        auto mi = index_to_MI[gid];
        handle_point(gid, ordered_triangles, mi);
        index_to_MI[gid] = mi;
    }


    cout << "DEBUG: ordered_triangles size after init = " << ordered_triangles.size() << endl;

    int current_remaining = 0;
    for (const auto &ch : chains) current_remaining += static_cast<int>(ch.size());
    int bound = current_remaining - remaining_vertices;
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

        // handle neighbor vertices: remove them from ordered set (they will be re-evaluated)
        // note that this does not need to account for >2 neighbors since junctions are never candidates
        auto it_neighbors = gid_to_original_neighbors.find(gid);
        if (it_neighbors != gid_to_original_neighbors.end())
        {
            for (int neigh_gid : it_neighbors->second)
            {
                if (neigh_gid < 0 || neigh_gid >= (int)index_to_MI.size()) continue;
                // use helper to remove neighbour entry from ordered_triangles and reset cursor
                handle_neighbour_global(neigh_gid, ordered_triangles, index_to_MI);
            }
        }

        // remove all occurrences of this global vertex from chains and points list
        for (int node_idx : gid_to_nodes[gid])
        {
            if (node_idx < 0 || node_idx >= (int)PI.size()) continue;
            // get chain id and pos
            int cid = PI[node_idx]->chain_id;
            int pos = chain_pos[node_idx];
            auto& vec = chains[cid];

            // try-catch in case pos is stale, but normally pos should be valid
            if (pos >= 0 && pos < (int)vec.size() && vec[pos] == node_idx)
            {
                vec.erase(vec.begin() + pos);
                // update chain_pos for nodes after pos
                for (int j = pos; j < (int)vec.size(); ++j)
                {
                    chain_pos[vec[j]] = j;
                }
            }
            // erase from points list
            points.erase(PI[node_idx]);
            // leave PI[node_idx] dangling (we won't use it again)
        }

        // mark removed globally
        global_removed[gid] = 1;
        // reset its iterator slot
        if (gid >= 0 && gid < (int)index_to_MI.size())
            index_to_MI[gid] = ordered_triangles.end();

        // unblock points: re-evaluate all global nodes that were candidates but currently not in ordered_triangles
        for (int j = 0; j < global_count; j++)
        {
            if (global_removed[j]) continue; // already removed
            if (!is_node_candidate_removable(j)) continue; // not a candidate (endpoints / junctions)
            if (index_to_MI[j] != ordered_triangles.end()) continue; // already in ordered_triangles

            auto mi = index_to_MI[j];
            handle_point(j, ordered_triangles, mi);
            index_to_MI[j] = mi;
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
    IPE::chains_to_IPE(name, out_chains, chain_closed, original);
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
    return init_vertex_count; // TODO: note this likely can just go for the more global version of the code
}
