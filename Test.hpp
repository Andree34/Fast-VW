#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "Fast_simplifier.hpp"

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
		"BoundaryBlock2"
	};

	// run all tests
	static void run_tests();

private:

	// run specific test
	static bool run_test(std::string test_name, bool verbose = true);
};