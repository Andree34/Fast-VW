#include "Utils.hpp"
#include "Fast_simplifier.hpp"
#include "Slow_simplifier_PSLG.hpp"
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

void IPE::normalize_chains(std::vector<Chain>& chains, double bound, double margin)
{
	// find global bounding box across all chains
	double minx = std::numeric_limits<double>::infinity();
	double miny = std::numeric_limits<double>::infinity();
	double maxx = -std::numeric_limits<double>::infinity();
	double maxy = -std::numeric_limits<double>::infinity();
	bool any = false;

	for (const auto &c : chains)
	{
		for (const auto &pt : c)
		{
			any = true;
			minx = std::min(minx, pt.first);
			miny = std::min(miny, pt.second);
			maxx = std::max(maxx, pt.first);
			maxy = std::max(maxy, pt.second);
		}
	}

	if (!any) return; // nothing to normalize

	// shift so min becomes zero (same idea as normalize_polygon)
	for (auto &c : chains)
	{
		for (auto &pt : c)
		{
			pt.first -= minx;
			pt.second -= miny;
		}
	}

	// compute maximum absolute values after shift
	double max_abs_x = 0.0;
	double max_abs_y = 0.0;
	for (const auto &c : chains)
	{
		for (const auto &pt : c)
		{
			max_abs_x = std::max(max_abs_x, std::abs(pt.first));
			max_abs_y = std::max(max_abs_y, std::abs(pt.second));
		}
	}

	double div_ratio = 1.0;
	double maxdim = std::max(max_abs_x, max_abs_y);
	if (maxdim > 0.0)
		div_ratio = maxdim / bound;
	else
		div_ratio = 1.0;

	// scale and add margin
	for (auto &c : chains)
	{
		for (auto &pt : c)
		{
			pt.first /= div_ratio;
			pt.second /= div_ratio;
			pt.first += margin;
			pt.second += margin;
		}
	}
}

void IPE::polygon_to_IPE(std::string name, Polygon polygon, bool original)
{
	assert(polygon.size() >= 3 && "Simplified polygon has less than 3 vertices.");

	// bounds the polygon to [-bound, bound] on both axes
	IPE::normalize_polygon(polygon, 500);

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

void IPE::chains_to_IPE(std::string name, const std::vector<Chain>& chains_in, std::vector<char> chain_closed, int global_vertices_left, bool original)
{
	// collect non-empty chains and copy them so we can normalise
	std::vector<Chain> chains;
	chains.reserve(chains_in.size());;

	for (const auto &c : chains_in)
	{
		if (c.empty()) continue;
		chains.push_back(c);
	}

	if (chains.empty()) return; // nothing to write

	normalize_chains(chains);

	// create output filename
	std::string file_name = (original ? "original_chains_" : "simplified_chains_") + std::to_string( global_vertices_left );
	std::ofstream fout("../data/" + name + "/" + file_name + ".ipe");

	fout << "<?xml version='1.0' encoding='utf-8'?>" << std::endl;
	fout << "<ipe version=\"70212\" creator=\"miniipe\"><ipestyle name=\"miniipe\" /><page><layer name=\"chains\" />";
	fout << std::setprecision(4) << std::fixed;

	// emit each chain as its own path element; do not fill (polylines)
	for (const auto &tc : chains)
	{
		fout << "<path stroke=\"black\" layer=\"chain\">";
		// move to first point
		fout << tc[0].first << " " << tc[0].second << " m ";
		// subsequent points as lines
		for (size_t i = 1; i < tc.size(); ++i)
		{
			fout << tc[i].first << " " << tc[i].second << " l ";
		}
		if (chain_closed[&tc - &chains[0]] && tc.size() > 1) // IPE does not like single points being closed shapes
			fout << " h "; // close if requested
		fout << "</path>";
	}

	fout << "</page></ipe>";
	fout.close();
}

// PRE: count <= 5236
template<typename T>
void run_stress_tests(std::vector<int> to_store, int start)
{
	for (size_t i = start; i < STRESS_TESTS; i++)
	{
		std::string name = std::to_string(i);
		Fast_simplifier<T> simplifier("StressTests/" + name, false);
		simplifier.polygon_to_ipe(true);
		simplifier.create_ipe_polygons(to_store);
		std::cout << "Polygon " << i << " done" << std::endl;
	}
}

template<typename T, typename Time_unit>
void generate_metrics_csv(bool console_output, int start)
{
	std::ofstream fout("../data/metrics.csv");
	fout << "TestName,VertexCount,Runtime(" + time_type<Time_unit>() + "),PITC,PITC/VertexCount,UPD,AvgDegree," << std::endl;
	for (size_t i = start; i < STRESS_TESTS; i++)
	{
		std::string name = std::to_string(i);
		T simplifier("StressTests/" + name);

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
			 << metric.pitc_to_vertex_count_ratio << ","
			 << metric.upd << ","
			 << metric.avg_degree << ","
			 << std::endl;
	}

	fout.close();
}

// yay code bloat
template void run_stress_tests<CDT>(std::vector<int> to_store, int count);
template void run_stress_tests<CT>(std::vector<int> to_store, int count);

template void generate_metrics_csv<Fast_simplifier<CDT>, std::chrono::nanoseconds >(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CDT>, std::chrono::microseconds>(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CDT>, std::chrono::milliseconds>(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CDT>, std::chrono::seconds     >(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CDT>, std::chrono::minutes     >(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CDT>, std::chrono::hours       >(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CT>, std::chrono::nanoseconds  >(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CT>, std::chrono::microseconds >(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CT>, std::chrono::milliseconds >(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CT>, std::chrono::seconds      >(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CT>, std::chrono::minutes      >(bool console_output, int count);
template void generate_metrics_csv<Fast_simplifier<CT>, std::chrono::hours        >(bool console_output, int count);
template void generate_metrics_csv<Slow_simplifier_PSLG, std::chrono::nanoseconds	  >(bool console_output, int count);
template void generate_metrics_csv<Slow_simplifier_PSLG, std::chrono::microseconds     >(bool console_output, int count);
template void generate_metrics_csv<Slow_simplifier_PSLG, std::chrono::milliseconds     >(bool console_output, int count);
template void generate_metrics_csv<Slow_simplifier_PSLG, std::chrono::seconds          >(bool console_output, int count);
template void generate_metrics_csv<Slow_simplifier_PSLG, std::chrono::minutes          >(bool console_output, int count);
template void generate_metrics_csv<Slow_simplifier_PSLG, std::chrono::hours            >(bool console_output, int count);

template<typename T> std::string time_type() { return "unknown"; }
template<> std::string time_type<std::chrono::nanoseconds >() { return "nanoseconds"; }
template<> std::string time_type<std::chrono::microseconds>() { return "microseconds"; }
template<> std::string time_type<std::chrono::milliseconds>() { return "milliseconds"; }
template<> std::string time_type<std::chrono::seconds     >() { return "seconds"; }
template<> std::string time_type<std::chrono::minutes     >() { return "minutes"; }
template<> std::string time_type<std::chrono::hours       >() { return "hours"; }