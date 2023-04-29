#pragma once

#include <string>
#include <fstream>
#include <vector>

class IPE
{
public:
	static void process_file(std::string name);
};

// runs all stress tests
void run_stress_tests(int count = 33);

// returns the time type as a string
template<typename T> std::string time_type();