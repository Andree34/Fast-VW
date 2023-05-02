#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "Fast_simplifier.hpp"
#include "Slow_simplifier.hpp"

template<typename T = Fast_simplifier<>>
class Test
{
public:
	// add test here
	static inline std::vector<std::string> test_names = {
		"Arrow",
		"TallArrow",
		"SimpleCollinear",
		"RationalNumbers",
		"BoundaryBlock",
		"BoundaryBlock2",
		"RandomTests/0",
		"RandomTests/1",
		"RandomTests/2",
		"RandomTests/3",
		"RandomTests/4",
		"RandomTests/5",
		"RandomTests/6",
		"RandomTests/7",
		"RandomTests/8",
		"RandomTests/9",
		"RandomTests/10"
	};

	// run all tests
	static void run_tests();

private:

	// run specific test
	static bool run_test(std::string test_name, bool verbose = true);
};