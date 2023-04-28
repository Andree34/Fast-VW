#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include "Fast_simplifier.hpp"
#include "Slow_simplifier.hpp"
#include "Test.hpp"
#include "Utils.hpp"

int main()
{
	/////////////////////////////////////////////////

	// Sample way of generating a big test case
	/*int n = 1000000;
	std::ofstream fout("../data/HugeZigZag/data.in");
	for (int i = 0; i < n; ++i)
		fout << i << " " << (i & 1) << std::endl;
	fout << poly_size << " " << 2 * poly_size << std::endl;*/

	/////////////////////////////////////////////////

	// Sample way of instantiating VW_computator
	auto start_time_full_algo = std::chrono::high_resolution_clock::now();

	Fast_simplifier<CDT> vw("BoundaryBlock");
	std::cout << "Unknown point detections: " << vw.get_UPD() << std::endl;
	std::cout << "Number of vertices removes: " << vw.result.size() << std::endl;

	//std::cout << "First 100 vertices removed: " << std::endl;
	//for (size_t i = 0; i < 100; i++)
	//	std::cout << vw.result[i] << std::endl;

	auto end_time_full_algo = std::chrono::high_resolution_clock::now();
	auto duration_full_algo = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_full_algo - start_time_full_algo);

	std::cout << "Milliseconds passed entire algorithm (test case " + vw.name + "): " << duration_full_algo.count() << " milliseconds." << std::endl;
	std::cout << "Point in triangle checks (test case " + vw.name + "): " << vw.get_PITC() << std::endl;
	vw.print_current_polygon();

	/////////////////////////////////////////////////

	//Test<>::run_tests();
	//std::cout << "Next batch:" << std::endl;
	//Test<Fast_simplifier<CT>>::run_tests();
	//IPE::process_file("Arrow");
}
