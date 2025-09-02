#pragma once

#include "Core.hpp"
#include <string>
#include <fstream>
#include <vector>
#include <chrono>
#define STRESS_TESTS 5236

class IPE
{
public:
	using Polygon = std::vector<std::pair<double, double>>;
	using Chain = std::vector<std::pair<double, double>>; // not necessarily closed
	static void process_file(std::string name);
	static void normalize_polygon(Polygon& polygon, double bound = 500, double margin = 10);
	static void IPE::normalize_chains(std::vector<Chain>& chains, double bound = 500, double margin = 10);
	static void polygon_to_IPE(std::string name, Polygon polygon, bool original);
	static void chains_to_IPE(std::string name, const std::vector<Chain>& chains, bool original = false);
};

// runs all stress tests
template<typename T = CDT>
void run_stress_tests(std::vector<int> to_store, int start = 0);

template<typename T = CDT, typename Time_unit = std::chrono::milliseconds>
void generate_metrics_csv(bool console_output = true, int start = 0);

// returns the time type as a string
template<typename T> 
std::string time_type();