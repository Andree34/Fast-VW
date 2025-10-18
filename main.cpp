#include "Core.hpp"
#include "Fast_simplifier.hpp"
#include "Fast_simplifier_PSLG.hpp"
#include "Slow_simplifier_PSLG.hpp"
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

	// Sample way of using a simplifier

	// Fast_simplifier<CDT> vw("StressTests/536", false);
	// vw.polygon_to_ipe(true);
	// vw.print_all_metrics();
	// vw.create_ipe_polygons({ 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000});
	// vw.print_current_polygon();


	// Fast_simplifier<CDT> vw3("BoundaryBlock", false);
	// vw3.polygon_to_ipe(true);
	// vw3.create_ipe_polygons({ 1,2,3,4,5,6,7,8,9});


	// Slow_simplifier_PSLG vw("BoundaryBlock", false, false);
	// Slow_simplifier_PSLG vw("IPE/extracted-test", false, false);
	Fast_simplifier_PSLG vw("IPE/gemeenten-2022_92952vtcs", false, false);
	// Slow_simplifier_PSLG vw("IPE/gemeenten-2022_19282vtcs", false, false);
	// Slow_simplifier_PSLG vw("IPE/gemeenten-2022_1929vtcs", false, false);
	// Slow_simplifier_PSLG vw("IPE/merge_test", false, false);
	// Slow_simplifier_PSLG vw("IPE/testing4", false, false);
	auto init_vertex_count = vw.get_initial_vertex_count();
	// do equispaced points of vertices for testing
	std::vector<int> sizes = { };
	const int steps = 20;
	for (int i = 1; i <= steps; ++i)
	{
		const int new_size = init_vertex_count - i * init_vertex_count / steps;
		if (sizes.empty() || sizes.back() != new_size)
			sizes.emplace_back( new_size );
	}
	// Filter out first size being same as input (can happen for small cases)
	if (!sizes.empty() && sizes.front() == init_vertex_count)
		sizes.erase(sizes.begin());
	vw.create_ipe_chains(sizes);



	// vw.print_all_metrics();

	//for (int i = STRESS_TESTS; i >=  500; i-= 100)
	//{
	//	Slow_simplifier vw("StressTests/" + std::to_string(i), true);
	//	vw.print_all_metrics();
	//}
	

	/////////////////////////////////////////////////

	// Sample way of using util functions
	//run_stress_tests({ 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000});
	//generate_metrics_csv<Fast_simplifier<>>(true, 4000);

	/////////////////////////////////////////////////

	// Use exact constructions when testing
	
	// Test<Slow_simplifier>::run_tests();
	// Test<Slow_simplifier_PSLG>::run_tests();
	//std::cout << "Next batch:" << std::endl;
	// Test<Fast_simplifier<CDT>>::run_tests();
	// std::cout << "Next batch:" << std::endl;
	// Test<Fast_simplifier<CT>>::run_tests();
	//IPE::process_file("StressTests/test");
}
