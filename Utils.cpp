#include "Utils.hpp"

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