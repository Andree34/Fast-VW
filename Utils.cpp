#include "Utils.hpp"
#include "Fast_simplifier.hpp"
#include <string>
#include <chrono>

void IPE::process_file(std::string name)
{
	using Row = std::pair<std::pair<std::string, std::string>, std::string>;
	std::ifstream fin("../data/" + name + "/data.in");

	std::vector<std::pair<std::string, std::string>> rows;
	std::string str1, str2, str3;
	while (fin >> str1 >> str2 >> str3)
		rows.push_back({ str1, str2 });

	std::ofstream fout("../data/" + name + "/data.in");
	for (auto [str1, str2] : rows)
		fout << str1 << " " << str2 << std::endl;
}

// PRE: count <= 33
void run_stress_tests(int count)
{
	for (size_t i = 0; i < count; i++)
	{
		std::string name = std::to_string(i);
		Fast_simplifier simplifier("StressTests/" + name);
		simplifier.print_all_metrics();
	}
}

template<typename T> std::string time_type() { return "unknown"; }
template<> std::string time_type<std::chrono::nanoseconds >() { return "nanoseconds"; }
template<> std::string time_type<std::chrono::microseconds>() { return "microseconds"; }
template<> std::string time_type<std::chrono::milliseconds>() { return "milliseconds"; }
template<> std::string time_type<std::chrono::seconds     >() { return "seconds"; }
template<> std::string time_type<std::chrono::minutes     >() { return "minutes"; }
template<> std::string time_type<std::chrono::hours       >() { return "hours"; }