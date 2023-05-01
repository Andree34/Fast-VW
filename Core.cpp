#include "Core.hpp"

std::ostream& operator<<(std::ostream& os, const Metric& metric)
{
	os << "Number of vertices initial polygon (test case " + metric.test_name + "): " << metric.init_vertex_count << std::endl;
	os << "Runtime of algorithm (test case " + metric.test_name + "): " << metric.runtime.first << " " + metric.runtime.second << std::endl;
	os << "Point in triangle checks (test case " + metric.test_name + "): " << metric.pitc << std::endl;
	os << "Unknown point detections (test case " + metric.test_name + "): " << metric.upd << std::endl;
	os << "Average vertex degree (test case " + metric.test_name + "): " << metric.avg_degree << std::endl;
	return os;
}