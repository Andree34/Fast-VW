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

class Slow_simplifier
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
	explicit Slow_simplifier(const std::string& input_folder_name, bool gen_test = false, bool auto_simplify = true);

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
		return Metric(init_vertex_count, name, get_runtime<Time_unit>(), get_PITC(), 0, 0);
	}

	// prints all registered metrics
	template<typename Time_unit = std::chrono::milliseconds>
	void print_all_metrics()
	{
		std::cout << get_metrics<Time_unit>();
	}

private:
	long long point_in_triangle_checks = 0;
	int init_vertex_count;

	Timestamp start_time;
	Timestamp end_time;

	// node list entry representing a single occurrence of a (possibly shared) vertex in a chain
	struct NodeEntry {
		Point p;                          // geometric point
		std::pair<int, int> meta;         // meta.first = node index, meta.second = 'block' cursor used by handle_point
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

	// occurrence bookkeeping for global vertex merging: maps global coordinate id -> list of node indices
	// (global ids are only used to detect junction/shared vertices as we cannot remove those)
	std::vector<int> node_to_global_vid;
	std::map<std::pair<K::FT, K::FT>, int> global_coord_to_vid;
	std::vector<int> global_vid_counts;

	// for junction detection
	std::vector<int> global_vid_chain_count;
	// builder helper to avoid counting a chain multiple times per gid
	std::vector<int> global_vid_chain_last_seen;

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
	/// PRE: 'pi' points to a valid NodeEntry (not removed).
	///      'removed' marks which node indices are already removed.
	///      'ordered_triangles' may or may not contain an entry for 'pi'.
	/// POST: If 'pi' corresponds to a removable vertex and is not blocked by another point,
	///       then an entry (area, node_index) -> point is inserted into 'ordered_triangles'.
	///       Otherwise, 'pi' is left untouched.
	/// </summary>
	void handle_point(Point_iterator pi,
	                  std::map<std::pair<K::FT, int>, Point>& ordered_triangles,
	                  const std::vector<char>& removed,
	                  std::map<std::pair<K::FT, int>, Point>::iterator& mi);

	/// <summary>
	/// PRE: 'pi' points to a valid NodeEntry (not removed).
	///      'mi' is the iterator to the entry in 'ordered_triangles' corresponding to 'pi',
	///      or 'ordered_triangles.end()' if not currently present.
	/// POST: Removes 'pi' from 'ordered_triangles' if present.
	///       Resets 'pi->meta.second' (block cursor) so the point will be re-evaluated later.
	///       Sets 'mi = ordered_triangles.end()'.
	/// </summary>
	void handle_neighbour(Point_iterator pi, std::map<std::pair<K::FT, int>, Point>& ordered_triangles, std::map<std::pair<K::FT, int>, Point>::iterator& mi);
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
	/// PRE: node_index is a valid index into PI.
	/// POST: Returns true if the node occurrence is a candidate for removal,
	///       i.e. not a junction (appears in exactly one distinct chain) and has both neighbours.
	/// </summary>
	[[nodiscard]] bool is_node_candidate_removable(int node_index) const;
};
