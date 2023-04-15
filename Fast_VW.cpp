#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <chrono>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::No_constraint_intersection_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CGAL::Constrained_triangulation_2<K, CGAL::Default, Itag> CT;

// the type T can (probably) be any type of constrained triangulation that is provided by CGAL
template<typename T = CDT>
struct VW_computator
{
	using Edge_iterator = typename T::Edge_iterator;
	using Vertex_handle = typename T::Vertex_handle;
	using Point = typename T::Point;
	using Edge = typename T::Edge;
	using Vertex_circulator = typename T::Vertex_circulator;
	using Vertex_iterator = typename std::list<std::pair<Vertex_handle, int>>::iterator;
	using Map_iterator = typename std::map<std::pair<K::FT, int>, Vertex_handle>::iterator;
	using Neighbours = typename std::pair<Vertex_handle, Vertex_handle>;

	// name of the folder where the data was taken from
	std::string name;

	int point_in_triangle_checks = 0;

	// points given in the input
	std::list<std::pair<Vertex_handle, int>> vertices;

	// array of vertex iterators (in this->vertices)
	std::vector<Vertex_iterator> VI;

	// vertex handle to vertex ID
	std::unordered_map<Vertex_handle, int> VH_to_id;
	
	// constrained triangulation of the polygon defined by points, 
	// where points contains the vertices of the polygon given in CCW order
	T ct;

	// vector that contains the sequence of vertex IDs, 
	// in order of vertex removal given by the VW algorithm
	std::vector<int> result;

	// get vertex iterator in the vertex list from a vertex handle
	Vertex_iterator get_vi(Vertex_handle vh)
	{
		return VI[VH_to_id[vh]];
	}

	// get vertex iterator in the vertex list from a vertex id
	Vertex_iterator get_vi(int id)
	{
		return VI[id];
	}

	// get neighbours of the vertex (on the current polygon)
	Neighbours get_neighbours(Vertex_handle vh)
	{
		Neighbours nb;
		Vertex_iterator curr = get_vi(vh);

		if (curr != vertices.begin())
		{
			auto prev = curr;
			--prev;
			nb.first = prev->first;
		}
		else nb.first = vertices.back().first;

		auto next = curr;
		++next;
		if (next != vertices.end())
			nb.second = next->first;
		else nb.second = vertices.front().first;

		return nb;
	}

	/// <summary>
	/// PRE: (tr1, tr2, tr3) is given in CCW order
	/// </summary>
	/// <param name="p">point that is checked</param>
	/// <param name="tr1">first triangle vertex</param>
	/// <param name="tr2">second triangle vertex</param>
	/// <param name="tr3">third triangle vertex</param>
	/// <returns> true if p is in (tr1, tr2, tr3), or on the boundary, false otherwise</returns>
	bool is_in_triangle(Point p, Point tr1, Point tr2, Point tr3)
	{
		CGAL::Orientation ori1 = CGAL::orientation(tr1, tr2, p);
		CGAL::Orientation ori2 = CGAL::orientation(tr2, tr3, p);
		CGAL::Orientation ori3 = CGAL::orientation(tr3, tr1, p);

		return !(ori1 == CGAL::RIGHT_TURN || ori2 == CGAL::RIGHT_TURN || ori3 == CGAL::RIGHT_TURN);
	}

	/// <param name="vh">Vertex handle that is checked for blocks</param>
	/// <returns>index of first blocking vertex found, or -1 if there is none</returns>
	int get_block(Vertex_handle vh)
	{
		// get the point's neighbour triangle
		auto [nb1, nb2] = get_neighbours(vh);

		// Swap points 2 and 3 to get ccw triangle
		if (CGAL::orientation(vh->point(), nb1->point(), nb2->point()) != CGAL::LEFT_TURN)
			std::swap(nb1, nb2);

		Point p1 = nb1->point();
		Point p2 = nb2->point();

		Vertex_circulator vc = ct.incident_vertices(vh), done(vc);
		do
		{
			// if the current vertex is a neighbour, skip it
			if (vc == nb1 || vc == nb2)
			{
				++vc;
				continue;
			}

			Point p_vc = vc->point();

			// check if point is in triangle here using vertex orientations
			// TODO: Also check if barycentric coordinate test is faster
			++point_in_triangle_checks;
			if(is_in_triangle(vc->point(), vh->point(), nb1->point(), nb2->point()))
				return VH_to_id[vc->handle()];
			++vc;

		} while (vc != done);

		return -1;
	}

	K::FT get_area(Vertex_handle vh)
	{
		Neighbours nb = get_neighbours(vh);
		K::Point_2 p1 = vh->point();
		K::Point_2 p2 = nb.first->point();
		K::Point_2 p3 = nb.second->point();

		return abs(K::Triangle_2(p1, p2, p3).area());
	}

	// adds the given vertex from the map depending on whether it's blocked or not
	void handle_vertex(Vertex_handle vh, std::map<std::pair<K::FT, int>, Vertex_handle>& ordered_triangles,
		std::vector<std::vector<int>>& blocked, Map_iterator& it)
	{
		int index = get_vi(vh)->second;
		int block = get_block(vh);
		if (block >= 0)
		{
			if(it == ordered_triangles.end())
				blocked[block].push_back(index);
			return;
		}

		if (it != ordered_triangles.end())
			return;

		it = ordered_triangles.insert(std::make_pair(std::make_pair(get_area(vh), index ), vh)).first;
	}

	// take out neighbours from the map, invalidate their index_to_SI entries, 
	// and then call handle_vertex on them
	void handle_neighbour(Vertex_handle vh, std::map<std::pair<K::FT, int>, Vertex_handle>& ordered_triangles,
		std::vector<std::vector<int>>& blocked, Map_iterator& it)
	{
		int index = VH_to_id[vh];

		// take out the vertex from the map
		if (it != ordered_triangles.end())
		{
			ordered_triangles.erase(it);
			it = ordered_triangles.end();
		}

		handle_vertex(vh, ordered_triangles, blocked, it);
	}

	void compute_result()
	{
		// sorts triangles by area, and then by vertex index (in order to make the algo predictable)
		// this also ensures unique keys, no multimap is not needed
		std::map<std::pair<K::FT, int>, Vertex_handle> ordered_triangles;

		// blocked[i] contains the indices of all vertices that are blocked
		std::vector<std::vector<int>> blocked(vertices.size());

		// keeps track of vertices that were removed (char is used since vector<bool> is bad practice)
		std::vector<char> removed(vertices.size());

		// maps each vertex index to the corresponding iterator in the map
		std::vector<Map_iterator> index_to_SI(vertices.size());
		for (auto& it : index_to_SI)
			it = ordered_triangles.end();

		// initialise ordered_triangles from initial triangulation
		for (int i = 0; i < vertices.size(); i++)
			handle_vertex(get_vi(i)->first, ordered_triangles, blocked, index_to_SI[i]);

		int remaining_vertices = 3;
		int bound = (int)vertices.size() - remaining_vertices;
		for (int i = 0; i < bound; i++)
		{
			while (get_block(ordered_triangles.begin()->second) >= 0)
			{
				assert(ordered_triangles.size() >= 1);

				Vertex_handle vh = ordered_triangles.begin()->second;
				ordered_triangles.erase(ordered_triangles.begin());

				int index = VH_to_id[vh];
				index_to_SI[index] = ordered_triangles.end();
				handle_vertex(vh, ordered_triangles, blocked, index_to_SI[index]);
			}

			assert(ordered_triangles.size() >= 1);

			// get vertex handle of next vertex that is removed (smallest area non-blocked)
			Vertex_handle best_vertex = ordered_triangles.begin()->second;
			ordered_triangles.erase(ordered_triangles.begin());

			// add index of vertex to result
			int index = VH_to_id[best_vertex];
			result.push_back(index);

			// get neighbours of vertex that will be removed
			auto [nb1, nb2] = get_neighbours(best_vertex);

			// make all the incident constraints regular triangulation edges
			ct.remove_incident_constraints(best_vertex);

			// remove vertex from polygon
			vertices.erase(get_vi(best_vertex));

			// track removal
			removed[index] = 1;

			// Remove the vertex from the CT.
			ct.remove(best_vertex);

			// add back the constraint to close the shape
			ct.insert_constraint(nb1, nb2);

			// unblock vertices that were blocked by best_vertex
			for (int to_unblock : blocked[index])
				if(!removed[to_unblock])
					handle_vertex(get_vi(to_unblock)->first, ordered_triangles, blocked, index_to_SI[to_unblock]);

			// handle neighbour vertices
			int index_nb1 = VH_to_id[nb1];
			int index_nb2 = VH_to_id[nb2];
			handle_neighbour(nb1, ordered_triangles, blocked, index_to_SI[index_nb1]);
			handle_neighbour(nb2, ordered_triangles, blocked, index_to_SI[index_nb2]);
		}
	}

	/// <summary>
	/// PRE: given folder contains a correctly formatted data.in file
	/// POST: points and ct are initialised according to the given input
	/// 
	/// Instructions on how *.in and *.out files shall be formatted are provided in the data folder
	/// </summary>
	/// <param name="input_folder_name">PRE: the name of the input folder that can be found in the data folder</param>
	VW_computator(std::string input_folder_name) :name(input_folder_name)
	{
		std::ifstream in("../data/" + input_folder_name + "/data.in");
		std::istream_iterator<Point> begin(in);
		std::istream_iterator<Point> end;

		std::vector<Point> points;
		int index = 0;
		for (auto p = begin; p != end; ++p, ++index)
		{
			// init vertices
			vertices.push_back({ ct.insert(*p), index });
			points.push_back(*p);
			VH_to_id[vertices.back().first] = index;

			// init vi
			Vertex_iterator vi = vertices.end();
			--vi;
			VI.push_back(vi);
		}

		// init CT
		ct.insert_constraint(points.begin(), points.end(), true);

		compute_result();
	}

	void print_result() const
	{
		for (auto index : result)
			std::cout << index << std::endl;
	}
};

namespace Test
{
	// add test here
	std::vector<std::string> test_names = {
		"BoundaryBlock",
		"Arrow"
	};

	bool run_test(std::string test_name, bool verbose = 1)
	{
		std::cout << "Running test: " << test_name << std::endl;

		std::vector<int> result = VW_computator(test_name).result;
		std::ifstream in("../data/" + test_name + "/data.out");
		std::istream_iterator<int> begin(in);
		std::istream_iterator<int> end;
		std::vector<int> expected_indices(begin, end);

		// Test is index count produced by algo matches expected index count
		assert((result.size() == expected_indices.size()) && "Resulting index count does not match test case");

		for (int i = 0; i < result.size(); i++)
			if (result[i] != expected_indices[i])
			{
				if (verbose)
				{
					std::cout << "Failed on test case: " << test_name << std::endl;
					std::cout << "Mismatch on sequence index: " << i << " (0-based)" << std::endl;
					std::cout << "Expected: " << expected_indices[i] << std::endl;
					std::cout << "Received: " << result[i] << std::endl;
				}
				else std::cout << "Failed" << std::endl;
				return false;
			}

		std::cout << "Passed" << std::endl;
		return true;
	}

	void run_tests()
	{
		std::cout << "Running tests..." << std::endl;
		for (auto& name : test_names)
		{
			std::cout << std::endl;
			if (!run_test(name))
				break;
		}
	}
}

namespace IPE
{
	void process_file(std::string name)
	{
		using Row = std::pair<std::pair<std::string, std::string>, std::string>;
		std::ifstream fin("../data/" + name + "/data.in");
		
		std::vector<std::pair<std::string, std::string>> rows;
		std::string str1, str2, str3;
		while (fin >> str1 >> str2 >> str3)
			rows.push_back({ str1, str2 });

		std::ofstream fout("../data/" + name + "/data.in");
		for (auto [str1, str2] : rows)
			fout << str1 << " " << str2 << std::endl;
	}
}

int main()
{
	// Sample way of generating a big test case
	/*int n = 1000000;
	std::ofstream fout("../data/HugeZigZag/data.in");
	for (int i = 0; i < n; ++i)
		fout << i << " " << (i & 1) << std::endl;
	fout << poly_size << " " << 2 * poly_size << std::endl;*/

	// Sample way of instantiating VW_computator
	auto start_time_full_algo = std::chrono::high_resolution_clock::now();

	VW_computator vw("India");

	auto end_time_full_algo = std::chrono::high_resolution_clock::now();
	auto duration_full_algo = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_full_algo - start_time_full_algo);

	std::cout << "Milliseconds passed entire algorithm (test case " + vw.name + "): " << duration_full_algo.count() << " milliseconds." << std::endl;
	std::cout << "Point in triangle checks (test case " + vw.name + "): " << vw.point_in_triangle_checks << std::endl;
	//vw.print_result();

	//Test::run_tests();
	//IPE::process_file("India");
}
