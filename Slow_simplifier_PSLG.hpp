#pragma once

#include "Core.hpp"
#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <list>
#include "Utils.hpp"

class Slow_simplifier_PSLG
{
public:
	using Point = CGAL::Point_2<K>;
	using Timestamp = std::chrono::steady_clock::time_point;

	// name of the folder where the data was taken from
	std::string name;

	// vector that contains the sequence of vertex IDs,
	// in order of vertex removal given by the VW algorithm
	std::vector<int> result;

	/// <summary>
	/// PRE: given folder contains a correctly formatted data.in file
	/// POST: points and ct are initialised according to the given input
	///
	/// Instructions on how *.in and *.out files shall be formatted are provided in the data folder
	/// </summary>
	/// <param name="input_folder_name">the name of the input folder that can be found in the data folder</param>
	/// <param name="gen_test">if true, result will be stored as expected output for given test</param>
	/// <param name="auto_simplify">if true, shape will be simplified to 3 vertices after constructor finishes executing</param>
	explicit Slow_simplifier_PSLG(const std::string& input_folder_name, bool gen_test = false, bool auto_simplify = true);

	// returns the number point in triangle checks
	[[nodiscard]] long long get_PITC() const;

	// returns the time it took to run the algorithm on the given input
	// only works if shape was auto_simplified
	template<typename Time_unit = std::chrono::milliseconds>
	Metric::Runtime get_runtime()
	{
		return { std::chrono::duration_cast<Time_unit>(end_time - start_time).count(), time_type<Time_unit>() };
	}

	template<typename Time_unit = std::chrono::milliseconds>
	Metric get_metrics()
	{
		return Metric(init_global_vertex_count, name, get_runtime<Time_unit>(), get_PITC(), 0, 0);
	}

	// prints all registered metrics
	template<typename Time_unit = std::chrono::milliseconds>
	void print_all_metrics()
	{
		std::cout << get_metrics<Time_unit>();
	}

	void chain_to_ipe(bool original = false);

	void create_ipe_chains(std::vector<int> polygon_sizes);

	[[nodiscard]] int get_initial_vertex_count() const;

private:
	long long point_in_triangle_checks = 0;
	int init_global_vertex_count;

	Timestamp start_time;
	Timestamp end_time;

	// node list entry representing a single occurrence of a (possibly shared) vertex in a chain
	struct NodeEntry {
		Point p;                          // geometric point
		int chain_id = -1;                // which chain this node occurrence belongs to
	};

	// points given in the input
	std::list<NodeEntry> points;

	// array of point iterators (in this->points), indexed by node-index (occurrence index)
	using Point_iterator = std::list<NodeEntry>::iterator;
	std::vector<Point_iterator> PI;

	// chain data: each chain is a list of node-indices (node indices are indices into the 'points' occurrences)
	std::vector<std::vector<int>> chains;

	// map node-index -> position in its chain (so we can get prev/next quickly)
	std::vector<int> chain_pos;

	// chain closed flags
	std::vector<char> chain_closed;

	// occurrence bookkeeping for global vertex merging: maps node-index -> global coord id
	std::vector<int> node_to_gid;
	// reverse mapping: global vid -> occurrence node indices
	std::vector<std::vector<int>> gid_to_nodes;

	// location to the global id
	std::map<std::pair<K::FT, K::FT>, int> global_coord_to_vid;
	// canonical point per global vid (coordinate)
	std::vector<Point> gid_to_point;

	// global neighbors of a vertex, mapped using global vertex ids. This mapping is updated
	// during simplification and always reflects the current immediate neighbours (for junction detection).
	std::map<int, std::set<int>> gid_to_neighbors;

	// Maps global vertex id -> whether it was removed
	std::vector<char> global_removed;

	// get point iterator in the vertex list from a vertex id
	[[nodiscard]] Point_iterator& get_pi(const int id) { return PI[id]; }
	[[nodiscard]] const Point_iterator& get_pi(const int id) const { return PI[id]; }

	// get neighbours of the point (on the current chain)
	// returns pair(prev_node_index, next_node_index). If neighbour doesn't exist (open endpoint), returns -1 for that neighbour.
	[[nodiscard]] std::pair<int, int> get_neighbours(int ind) const;

	/// <summary>
	/// PRE: (tr1, tr2, tr3) is given in CCW order
	/// </summary>
	/// <param name="p">point that is checked</param>
	/// <param name="tr1">first triangle point</param>
	/// <param name="tr2">second triangle point</param>
	/// <param name="tr3">third triangle point</param>
	/// <returns> true if p is in (tr1, tr2, tr3), or on the boundary, false otherwise</returns>
	[[nodiscard]] bool is_in_triangle(Point p, Point tr1, Point tr2, Point tr3) const;

	// returns the area of the triangle corresponding to vh in the polygon/chain
	[[nodiscard]] K::FT get_area(int ind);

	/// <summary>
	/// PRE: 'gid' is a valid global vertex id (exists in global_vid_to_point)
	/// POST: If 'gid' corresponds to a removable global vertex and is not blocked by another point,
	///       then an entry (area, global_vid) -> point is inserted into 'ordered_triangles'.
	///       Otherwise, it is left untouched.
	/// </summary>
	void handle_point(int gid,
	                  std::map<std::pair<K::FT, int>, Point>& ordered_triangles,
	                  std::map<std::pair<K::FT, int>, Point>::iterator& mi);

	/// <summary>
	/// Remove the global vertex 'gid' from ordered_triangles (if present) and reset its blocking cursor.
	/// </summary>
	void handle_neighbour_global(int gid,
	                             std::map<std::pair<K::FT, int>, Point>& ordered_triangles,
	                             std::vector<std::map<std::pair<K::FT, int>, Point>::iterator>& index_to_MI);

	/// <summary>
	/// PRE: "remaining_vertices" <= vertices.size()
	///
	/// POST: Appends the VW index sequence computed by simplifying the shape down to the number of remaining verticies specified.
	///		  Chains are also simplified to the specified number of vertices where appropriate.
	///
	/// Note: It is considerably faster to call this function only once compared to multiple times.
	/// </summary>
	/// <param name="remaining_vertices">The number of remaining vertices in the resulting polygon(s)</param>
	void simplify(int remaining_vertices = 3);

	void generate_test_output();

	// helpers for PSLG parsing & bookkeeping
	std::vector<std::vector<Point>> parse_input_file_as_chains(const std::string& path);
	void build_internal_structures_from_chains(const std::vector<std::vector<Point>>& chains_points);

	/// <summary>
	/// PRE: gid is a valid global vertex id.
	/// POST: Returns true if the global vertex is a candidate for removal,
	///       i.e. not a junction (appears with exactly two distinct global neighbors) and none of its occurrences are endpoints.
	/// </summary>
	[[nodiscard]] bool is_node_candidate_removable(int gid) const;
	unsigned long long Slow_simplifier_PSLG::get_vertices_left() const;
};
