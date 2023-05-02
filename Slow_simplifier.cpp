#include "Slow_simplifier.hpp"

Slow_simplifier::Slow_simplifier(std::string input_folder_name, bool gen_test) :name(input_folder_name)
{
	start_time = std::chrono::high_resolution_clock::now();
	end_time = start_time;

	std::ifstream in("../data/" + input_folder_name + "/data.in");
	std::istream_iterator<Point> begin(in);
	std::istream_iterator<Point> end;

	int index = 0;
	for (auto p = begin; p != end; ++p, ++index)
	{
		// init vertices
		points.push_back({ *p, {index, 0} });

		// init vi
		Point_iterator pi = points.end();
		--pi;
		PI.push_back(pi);
	}

	// init CT
	init_vertex_count = (int)points.size();

	simplify();
	end_time = std::chrono::high_resolution_clock::now();

	if (gen_test)
		generate_test_output();
}

long long Slow_simplifier::get_PITC() const
{
	return point_in_triangle_checks;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

Slow_simplifier::Point_iterator Slow_simplifier::get_pi(int id) const
{
	return PI[id];
}

std::pair<int, int> Slow_simplifier::get_neighbours(int ind)
{
	std::pair<int, int> nb;
	Point_iterator curr = get_pi(ind);

	if (curr != points.begin())
	{
		auto prev = curr;
		--prev;
		nb.first = prev->second.first;
	}
	else nb.first = points.back().second.first;

	auto next = curr;
	++next;
	if (next != points.end())
		nb.second = next->second.first;
	else nb.second = points.front().second.first;

	return nb;
}

bool Slow_simplifier::is_in_triangle(Point p, Point tr1, Point tr2, Point tr3) const
{
	// handle degenerate case: collinear points (never blocked since the shape is a simple polygon)
	if(CGAL::collinear(tr1, tr2, tr3))
		return false;

	CGAL::Orientation ori1 = CGAL::orientation(tr1, tr2, p);
	CGAL::Orientation ori2 = CGAL::orientation(tr2, tr3, p);
	CGAL::Orientation ori3 = CGAL::orientation(tr3, tr1, p);

	return !(ori1 == CGAL::RIGHT_TURN || ori2 == CGAL::RIGHT_TURN || ori3 == CGAL::RIGHT_TURN);
}

K::FT Slow_simplifier::get_area(int ind)
{
	auto [nb1, nb2] = get_neighbours(ind);
	Point p1 = get_pi(ind)->first;
	Point p2 = get_pi(nb1)->first;
	Point p3 = get_pi(nb2)->first;

	return abs(K::Triangle_2(p1, p2, p3).area());
}

void Slow_simplifier::handle_point(Point_iterator& pi, 
	std::map<std::pair<K::FT, int>, Point>& ordered_triangles, std::vector<char> removed, Map_iterator& mi)
{
	auto& [p, vals] = *pi;
	auto& [ind, block] = vals;
	auto [nb1, nb2] = get_neighbours(ind);

	// Swap points 2 and 3 to get ccw triangle
	if (CGAL::orientation(pi->first, get_pi(nb1)->first, get_pi(nb2)->first) != CGAL::LEFT_TURN)
		std::swap(nb1, nb2);

	for (; block < init_vertex_count; block++)
	{
		if (removed[block])
			continue;

		if (block == ind || block == nb1 || block == nb2)
			continue;

		Point oth = get_pi(block)->first;
		++point_in_triangle_checks;
		if (is_in_triangle(oth, p, get_pi(nb1)->first, get_pi(nb2)->first))
			break;
	}

	if (block == init_vertex_count)
		mi = ordered_triangles.insert(std::make_pair(std::make_pair(get_area(ind), ind), p)).first;
}

void Slow_simplifier::handle_neighbour(Point_iterator& pi,
	std::map<std::pair<K::FT, int>, Point>& ordered_triangles, Map_iterator& mi)
{
	if(mi != ordered_triangles.end())
		ordered_triangles.erase(mi);
	mi = ordered_triangles.end();
	pi->second.second = 0;
}

void Slow_simplifier::simplify(int remaining_vertices)
{
	// sorts triangles by area, and then by vertex index (in order to make the algo predictable)
	// this also ensures unique keys, no multimap is not needed
	std::map<std::pair<K::FT, int>, Point> ordered_triangles;

	// keeps track of vertices that were removed (char is used since vector<bool> is bad practice)
	std::vector<char> removed(init_vertex_count);

	// maps each vertex index to the corresponding iterator in the map
	std::vector<Map_iterator> index_to_MI(init_vertex_count);
	for (auto& it : index_to_MI)
		it = ordered_triangles.end();

	// initialise ordered_triangles from initial triangulation
	for (size_t i = 0; i < init_vertex_count; i++)
		handle_point(get_pi(i), ordered_triangles, removed, index_to_MI[i]);

	int bound = init_vertex_count - remaining_vertices;
	for (int i = 0; i < bound; i++)
	{
		assert(ordered_triangles.size() >= 1);

		// get vertex handle of next vertex that is removed (smallest area non-blocked)
		Point best_point = ordered_triangles.begin()->second;
		int index = ordered_triangles.begin()->first.second;
		ordered_triangles.erase(ordered_triangles.begin());

		// add index of vertex to result
		result.push_back(index);

		// handle neighbour vertices
		auto [index_nb1, index_nb2] = get_neighbours(index);
		handle_neighbour(get_pi(index_nb1), ordered_triangles, index_to_MI[index_nb1]);
		handle_neighbour(get_pi(index_nb2), ordered_triangles, index_to_MI[index_nb2]);

		// remove vertex from polygon
		points.erase(get_pi(index));

		// track removal
		removed[index] = 1;

		// unblock points
		for (int i = 0; i < init_vertex_count; i++)
		{
			if (removed[i])
				continue;

			if (index_to_MI[i] != ordered_triangles.end())
				continue;

			handle_point(get_pi(i), ordered_triangles, removed, index_to_MI[i]);
		}
	}
}

void Slow_simplifier::generate_test_output()
{
	std::ofstream fout("../data/" + name + "/data.out");
	for (auto index : result)
		fout << index << std::endl;
}