#include "Utils.hpp"
#include "Fast_simplifier.hpp"
#include <chrono>
#include <iomanip>
#include <algorithm>

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

void IPE::normalize_polygon(Polygon& polygon, double bound, double margin)
{
	using Point = std::pair<double, double>;

	// make all coordinates positive
	double minx = min_element(polygon.begin(), polygon.end(), [](Point p1, Point p2)
		{return p1.first < p2.first; })->first;
	double miny = min_element(polygon.begin(), polygon.end(), [](Point p1, Point p2)
		{return p1.second < p2.second; })->second;
	for (auto& [x, y] : polygon)
	{
		x -= minx;
		y -= miny;
	}

	// scale all coordinates to be in [-bound, bound] on both axes
	double maxx = max_element(polygon.begin(), polygon.end(), [](Point p1, Point p2)
		{return abs(p1.first) < abs(p2.first); })->first;
	double maxy = max_element(polygon.begin(), polygon.end(), [](Point p1, Point p2) 
		{return abs(p1.second) < abs(p2.second); })->second;
	double div_ratio = std::max(maxx, maxy) / bound;

	for (auto& [x, y] : polygon)
	{
		x /= div_ratio;
		y /= div_ratio;

		x += margin;
		y += margin;
	}
}

void IPE::polygon_to_IPE(std::string name, Polygon polygon, bool original)
{
	assert(polygon.size() >= 3 && "Simplified polygon has less than 3 vertices.");

	// bounds the polygon to [-bound, bound] on both axes
	IPE::normalize_polygon(polygon, 5000);

	std::string file_name = (original ? "original_" : "simplified_") + std::to_string(polygon.size());
	std::ofstream fout("../data/" + name + "/" + file_name + ".ipe");
	fout << "<?xml version='1.0' encoding='utf-8'?>" << std::endl;
	fout << "<ipe version=\"70212\" creator=\"miniipe\"><ipestyle name=\"miniipe\" /><page><layer name=\"my_layer\" /><path stroke=\"black\" fill=\"#000\" layer=\"my layer\">";
	fout << std::setprecision(4) << std::fixed;
	fout << polygon[0].first << " " << polygon[0].second << " m ";
	for (size_t i = 1; i < polygon.size(); i++)
	{
		auto [x, y] = polygon[i];
		fout << x << " " << y << " l  ";
	}
	fout << "h </path></page></ipe>";
	fout.close();
}

// PRE: count <= 5236
template<typename T>
void run_stress_tests(std::vector<int> to_store, int count)
{
	for (size_t i = 0; i < count; i++)
	{
		std::string name = std::to_string(i);
		Fast_simplifier<T> simplifier("StressTests/" + name, false);
		simplifier.polygon_to_ipe(true);
		simplifier.create_ipe_polygons(to_store);
		std::cout << "Polygon " << i << " done" << std::endl;
	}
}

template<typename T, typename Time_unit>
void generate_metrics_csv(bool console_output, int count)
{
	std::ofstream fout("../data/metrics.csv");
	fout << "TestName,VertexCount,Runtime(" + time_type<Time_unit>() + "),PITC,UPD,AvgDegree" << std::endl;
	for (size_t i = 0; i < count; i++)
	{
		std::string name = std::to_string(i);
		Fast_simplifier<T> simplifier("StressTests/" + name);

		if (console_output)
		{
			simplifier.print_all_metrics<Time_unit>();
			std::cout << std::endl;
		}
		
		Metric metric = simplifier.get_metrics<Time_unit>();
		fout << metric.test_name << ","
			 << metric.init_vertex_count << ","
			 << metric.runtime.first << ","
			 << metric.pitc << ","
			 << metric.upd << ","
			 << metric.avg_degree << ","
			 << std::endl;
	}

	fout.close();
}

// yay code bloat
template void run_stress_tests<CDT>(std::vector<int> to_store, int count);
template void run_stress_tests<CT>(std::vector<int> to_store, int count);

template void generate_metrics_csv<CDT, std::chrono::nanoseconds >(bool console_output, int count);
template void generate_metrics_csv<CDT, std::chrono::microseconds>(bool console_output, int count);
template void generate_metrics_csv<CDT, std::chrono::milliseconds>(bool console_output, int count);
template void generate_metrics_csv<CDT, std::chrono::seconds     >(bool console_output, int count);
template void generate_metrics_csv<CDT, std::chrono::minutes     >(bool console_output, int count);
template void generate_metrics_csv<CDT, std::chrono::hours       >(bool console_output, int count);
template void generate_metrics_csv<CT, std::chrono::nanoseconds  >(bool console_output, int count);
template void generate_metrics_csv<CT, std::chrono::microseconds >(bool console_output, int count);
template void generate_metrics_csv<CT, std::chrono::milliseconds >(bool console_output, int count);
template void generate_metrics_csv<CT, std::chrono::seconds      >(bool console_output, int count);
template void generate_metrics_csv<CT, std::chrono::minutes      >(bool console_output, int count);
template void generate_metrics_csv<CT, std::chrono::hours        >(bool console_output, int count);

template<typename T> std::string time_type() { return "unknown"; }
template<> std::string time_type<std::chrono::nanoseconds >() { return "nanoseconds"; }
template<> std::string time_type<std::chrono::microseconds>() { return "microseconds"; }
template<> std::string time_type<std::chrono::milliseconds>() { return "milliseconds"; }
template<> std::string time_type<std::chrono::seconds     >() { return "seconds"; }
template<> std::string time_type<std::chrono::minutes     >() { return "minutes"; }
template<> std::string time_type<std::chrono::hours       >() { return "hours"; }