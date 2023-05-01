#include "Fast_simplifier.hpp"
#include "Utils.hpp"

template<typename T>
Fast_simplifier<T>::Fast_simplifier(std::string input_folder_name, bool auto_simplify) :name(input_folder_name)
{
	start_time = std::chrono::high_resolution_clock::now();
	end_time = start_time;

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
	init_vertex_count = (int)vertices.size();

	if (auto_simplify)
	{
		simplify();
		end_time = std::chrono::high_resolution_clock::now();
	}
}

template<typename T>
inline void Fast_simplifier<T>::print_result() const
{
	for (auto index : result)
		std::cout << index << std::endl;
}

template<typename T>
void Fast_simplifier<T>::print_current_polygon() const
{
	std::cout << std::setprecision(4) << std::fixed;
	for (auto& [vh, index] : vertices)
		std::cout << vh->point().x() << " " << vh->point().y() << std::endl;
}

template<typename T>
inline long long Fast_simplifier<T>::get_PITC() const
{
	return point_in_triangle_checks;
}

template<typename T>
inline long long Fast_simplifier<T>::get_UPD() const
{
	return unknown_point_detections;
}

template<typename T>
double Fast_simplifier<T>::get_avg_degree() const
{
	return (double)total_degree / (double)degree_registrations;
}

template<typename T>
void Fast_simplifier<T>::polygon_to_ipe(bool original)
{
	IPE::Polygon polygon;
	for (auto [vh, index] : vertices)
	{
		auto x = vh->point().x();
		auto y = vh->point().y();
		polygon.push_back({ CGAL::to_double(x), CGAL::to_double(y) });
	}

	IPE::polygon_to_IPE(name, polygon, original);
}

template<typename T>
void Fast_simplifier<T>::create_ipe_polygons(std::vector<int> polygon_sizes)
{
	sort(polygon_sizes.begin(), polygon_sizes.end(), std::greater<>());

	for (auto vertex_count : polygon_sizes)
	{
		if (vertex_count > vertices.size())
			continue;

		simplify(vertex_count);
		polygon_to_ipe();
	}
}

template<typename T>
void Fast_simplifier<T>::simplify(int remaining_vertices)
{
	assert(remaining_vertices <= (int)vertices.size() && "Shape is being simplified to more vertices than currently possible");

	// sorts triangles by area, and then by vertex index (in order to make the algo predictable)
	// this also ensures unique keys, no multimap is not needed
	std::map<std::pair<K::FT, int>, Vertex_handle> ordered_triangles;

	// blocked[i] contains the indices of all vertices that are blocked
	std::vector<std::vector<int>> blocked(init_vertex_count);

	// keeps track of vertices that were removed (char is used since vector<bool> is bad practice)
	std::vector<char> removed(init_vertex_count);

	// maps each vertex index to the corresponding iterator in the map
	std::vector<Map_iterator> index_to_SI(init_vertex_count);
	for (auto& it : index_to_SI)
		it = ordered_triangles.end();

	// initialise ordered_triangles from initial triangulation
	for (auto& [vh, i] : vertices)
		handle_vertex(get_vi(i)->first, ordered_triangles, blocked, index_to_SI[i]);

	int bound = (int)vertices.size() - remaining_vertices;
	for (int i = 0; i < bound; i++)
	{
		/*std::cerr << "Iteration " << i << std::endl;
		for (auto [comp, vh] : ordered_triangles)
		{
			std::cerr << comp.second << std::endl;
		}*/

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
			if (!removed[to_unblock])
				handle_vertex(get_vi(to_unblock)->first, ordered_triangles, blocked, index_to_SI[to_unblock]);

		// handle neighbour vertices
		int index_nb1 = VH_to_id[nb1];
		int index_nb2 = VH_to_id[nb2];
		handle_neighbour(nb1, ordered_triangles, blocked, index_to_SI[index_nb1]);
		handle_neighbour(nb2, ordered_triangles, blocked, index_to_SI[index_nb2]);
	}
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

template<typename T>
inline typename Fast_simplifier<T>::Vertex_iterator Fast_simplifier<T>::get_vi(Vertex_handle vh)
{
	return VI[VH_to_id[vh]];
}

template<typename T>
inline typename Fast_simplifier<T>::Vertex_iterator Fast_simplifier<T>::get_vi(int id) const
{
	return VI[id];
}

template<typename T>
inline typename Fast_simplifier<T>::Neighbours Fast_simplifier<T>::get_neighbours(Vertex_handle vh)
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

template<typename T>
inline bool Fast_simplifier<T>::is_in_triangle(Point p, Point tr1, Point tr2, Point tr3) const
{
	CGAL::Orientation ori1 = CGAL::orientation(tr1, tr2, p);
	CGAL::Orientation ori2 = CGAL::orientation(tr2, tr3, p);
	CGAL::Orientation ori3 = CGAL::orientation(tr3, tr1, p);

	return !(ori1 == CGAL::RIGHT_TURN || ori2 == CGAL::RIGHT_TURN || ori3 == CGAL::RIGHT_TURN);
}

template<typename T>
inline int Fast_simplifier<T>::get_block(Vertex_handle vh)
{
	// log degree of vertex
	++degree_registrations;
	total_degree += ct.degree(vh);

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
		// skip points added to the triangulation by CGAL
		if (VH_to_id.find(vc->handle()) == VH_to_id.end())
		{
			++unknown_point_detections;
			//std::cerr << vc->point().x() << " " << vc->point().y() << std::endl;
			++vc;
			continue;
		}

		// if the current vertex is a neighbour, skip it
		if (vc->handle() == nb1 || vc->handle() == nb2)
		{
			++vc;
			continue;
		}

		Point p_vc = vc->point();

		// check if point is in triangle here using vertex orientations
		// TODO: Also check if barycentric coordinate test is faster
		++point_in_triangle_checks;
		if (is_in_triangle(vc->point(), vh->point(), nb1->point(), nb2->point()))
			return VH_to_id[vc->handle()];
		++vc;

	} while (vc != done);

	return -1;
}

template<typename T>
inline K::FT Fast_simplifier<T>::get_area(Vertex_handle vh)
{
	Neighbours nb = get_neighbours(vh);
	K::Point_2 p1 = vh->point();
	K::Point_2 p2 = nb.first->point();
	K::Point_2 p3 = nb.second->point();

	return abs(K::Triangle_2(p1, p2, p3).area());
}

template<typename T>
inline void Fast_simplifier<T>::handle_vertex(Vertex_handle vh, std::map<std::pair<K::FT, int>, Vertex_handle>& ordered_triangles,
	std::vector<std::vector<int>>& blocked, Map_iterator& it)
{
	int index = get_vi(vh)->second;
	int block = get_block(vh);
	if (block >= 0)
	{
		if (it == ordered_triangles.end())
			blocked[block].push_back(index);
		return;
	}

	if (it != ordered_triangles.end())
		return;

	it = ordered_triangles.insert(std::make_pair(std::make_pair(get_area(vh), index), vh)).first;
}

template<typename T>
inline void Fast_simplifier<T>::handle_neighbour(Vertex_handle vh, std::map<std::pair<K::FT, int>, Vertex_handle>& ordered_triangles,
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

template class Fast_simplifier<CT>;
template class Fast_simplifier<CDT>;