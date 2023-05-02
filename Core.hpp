#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <string>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::No_constraint_intersection_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CGAL::Constrained_triangulation_2<K, CGAL::Default, Itag> CT;

// stores all relevant metrics of a simplification run
struct Metric
{
	using Runtime = std::pair<long long, std::string>;

	long long init_vertex_count;

	std::string test_name;
	Runtime runtime;

	long long pitc;
	long long upd;		// might not be useful, this was more important during testing
	double avg_degree;

	Metric(long long init_vertex_count, std::string test_name, Runtime runtime, long long pitc, long long upd, double avg_degree) :
		init_vertex_count(init_vertex_count), test_name(test_name), runtime(runtime), pitc(pitc), upd(upd), avg_degree(avg_degree) {}

	friend std::ostream& operator<<(std::ostream& os, const Metric& metric);
};