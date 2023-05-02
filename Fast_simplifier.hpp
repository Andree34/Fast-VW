#pragma once

#include "Core.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <iomanip>

template<typename T = CDT>
class Fast_simplifier
{
public:
	using Edge_iterator = typename T::Edge_iterator;
	using Vertex_handle = typename T::Vertex_handle;
	using Point = typename T::Point;
	using Edge = typename T::Edge;
	using Vertex_circulator = typename T::Vertex_circulator;
	using Vertex_iterator = typename std::list<std::pair<Vertex_handle, int>>::iterator;
	using Map_iterator = typename std::map<std::pair<K::FT, int>, Vertex_handle>::iterator;
	using Neighbours = typename std::pair<Vertex_handle, Vertex_handle>;
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
	/// <param name="auto_simplify">if true, shape will be simplied to 3 vertices after constructor finishes executing</param>
	Fast_simplifier(std::string input_folder_name, bool auto_simplify = true);

	// prints the resulting sequence of vertices removed
	void print_result() const;

	// prints the current polygon being stored
	void print_current_polygon() const;

	// returns the number point in triangle checks
	long long get_PITC() const;

	// returns the number of points added by CGAL detected during the algorithm run
	long long get_UPD() const;

	// returns the average degree of the vertices that are being checked for blocks
	double get_avg_degree() const;

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
		return Metric(init_vertex_count, name, get_runtime<Time_unit>(), get_PITC(), get_UPD(), get_avg_degree());
	}

	// prints all registered metrics
	template<typename Time_unit = std::chrono::milliseconds>
	void print_all_metrics()
	{
		std::cout << get_metrics<Time_unit>();
	}

	void polygon_to_ipe(bool original = false);

	void create_ipe_polygons(std::vector<int> polygon_sizes);

	/// <summary>
	/// PRE: "remaining_vertices" <= vertices.size()
	/// 
	/// POST: Appends the VW index sequence computed by simplifying the shape down to the number of remaining verticies specified.
	///		  Polygon is also simplified to the specified number of vertices.
	/// 
	/// Note: It is considerably faster to call this function only once compared to multiple times.
	/// </summary>
	/// <param name="remaining_vertices">The number of remaining vertices in the resulting polygon</param>
	void simplify(int remaining_vertices = 3);

private:
	long long point_in_triangle_checks = 0;
	long long unknown_point_detections = 0;
	int init_vertex_count;

	int total_degree = 0;
	int degree_registrations = 0;

	Timestamp start_time;
	Timestamp end_time;

	// points given in the input
	std::list<std::pair<Vertex_handle, int>> vertices;

	// array of vertex iterators (in this->vertices)
	std::vector<Vertex_iterator> VI;

	// vertex handle to vertex ID
	std::unordered_map<Vertex_handle, int> VH_to_id;

	// constrained triangulation of the polygon defined by points, 
	// where points contains the vertices of the polygon given in CCW order
	T ct;

	// get vertex iterator in the vertex list from a vertex handle
	Vertex_iterator get_vi(Vertex_handle vh);

	// get vertex iterator in the vertex list from a vertex id
	Vertex_iterator get_vi(int id) const;

	// get neighbours of the vertex (on the current polygon)
	Neighbours get_neighbours(Vertex_handle vh);

	/// <summary>
	/// PRE: (tr1, tr2, tr3) is given in CCW order
	/// </summary>
	/// <param name="p">point that is checked</param>
	/// <param name="tr1">first triangle vertex</param>
	/// <param name="tr2">second triangle vertex</param>
	/// <param name="tr3">third triangle vertex</param>
	/// <returns> true if p is in (tr1, tr2, tr3), or on the boundary, false otherwise</returns>
	bool is_in_triangle(Point p, Point tr1, Point tr2, Point tr3) const;

	/// <param name="vh">Vertex handle that is checked for blocks</param>
	/// <returns>index of first blocking vertex found, or -1 if there is none</returns>
	int get_block(Vertex_handle vh);
	
	// returns the area of the triangle corresponding to vh in the polygon
	K::FT get_area(Vertex_handle vh);

	// adds the given vertex from the map depending on whether it's blocked or not
	void handle_vertex(Vertex_handle vh, std::map<std::pair<K::FT, int>, Vertex_handle>& ordered_triangles,
		std::vector<std::vector<int>>& blocked, Map_iterator& it);

	// take out neighbours from the map, invalidate their index_to_SI entries, 
	// and then call handle_vertex on them
	void handle_neighbour(Vertex_handle vh, std::map<std::pair<K::FT, int>, Vertex_handle>& ordered_triangles,
		std::vector<std::vector<int>>& blocked, Map_iterator& it);
};