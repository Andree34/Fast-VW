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

	/*Fast_simplifier<CDT> vw("StressTests/5");
	vw.print_all_metrics();*/

	////std::cout << "First 100 vertices removed: " << std::endl;
	////for (size_t i = 0; i < 100; i++)
	////	std::cout << vw.result[i] << std::endl;

	//vw.print_current_polygon();

	run_stress_tests();

	/////////////////////////////////////////////////

	/*Test<>::run_tests();
	std::cout << "Next batch:" << std::endl;
	Test<Fast_simplifier<CT>>::run_tests();*/
	//IPE::process_file("Arrow");
}
