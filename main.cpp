#include "Core.hpp"
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

	// Sample way of using a simplifier

	//Fast_simplifier<CDT> vw("StressTests/536", false);
	//vw.polygon_to_ipe(true);
	//vw.print_all_metrics();
	//vw.create_ipe_polygons({ 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000});
	//vw.print_current_polygon();

	Slow_simplifier vw("StressTests/100");
	vw.print_all_metrics();

	//for (int i = STRESS_TESTS; i >= 500; i-= 100)
	//{
	//	Slow_simplifier vw("StressTests/" + std::to_string(i), true);
	//	vw.print_all_metrics();
	//}
	

	/////////////////////////////////////////////////

	// Sample way of using util functions
	//run_stress_tests({ 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000});
	//generate_metrics_csv<Slow_simplifier>(true, 100);

	/////////////////////////////////////////////////

	// Use exact constructions when testing
	
	//Test<Slow_simplifier>::run_tests();
	//std::cout << "Next batch:" << std::endl;
	//Test<Fast_simplifier<CDT>>::run_tests();
	//std::cout << "Next batch:" << std::endl;
	//Test<Fast_simplifier<CT>>::run_tests();
	//IPE::process_file("StressTests/test");
}
