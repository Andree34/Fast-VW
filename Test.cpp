#include "Test.hpp"

template<typename T>
void Test<T>::run_tests()
{
	std::cout << "Running tests..." << std::endl;
	for (auto& name : test_names)
	{
		std::cout << "+-+-+-+-+-+-+-+-+-+-+-+" << std::endl;
		bool result = run_test(name);
		std::cout << "+-+-+-+-+-+-+-+-+-+-+-+" << std::endl;
		std::cout << std::endl;

		if (!result)
			return;
	}

	std::cout << std::endl << "All tests passed!" << std::endl;
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

template<typename T>
bool Test<T>::run_test(std::string test_name, bool verbose)
{
	std::cout << "Running test: " << test_name << std::endl;

	std::vector<int> result = T(test_name).result;
	std::ifstream in("../data/" + test_name + "/data.out");
	std::istream_iterator<int> begin(in);
	std::istream_iterator<int> end;
	std::vector<int> expected_indices(begin, end);

	// Test if index count produced by algo matches expected index count
	if (result.size() != expected_indices.size())
	{
		if (verbose)
		{
			std::cout << "Failed on test case: " << test_name << std::endl;
			std::cout << "Resulting index count does not match test case" << std::endl;
		}
		else std::cout << "Failed" << std::endl;
		return false;
	}

	for (int i = 0; i < result.size(); i++)
		if (result[i] != expected_indices[i])
		{
			if (verbose)
			{
				std::cout << "Failed on test case: " << test_name << std::endl;
				std::cout << "Mismatch on sequence index: " << i << " (0-based)" << std::endl;
				std::cout << "Expected: " << expected_indices[i] << std::endl;
				std::cout << "Received: " << result[i] << std::endl;
			}
			else std::cout << "Failed" << std::endl;
			return false;
		}

	std::cout << "Passed" << std::endl;
	return true;
}

template class Test<Fast_simplifier<CDT>>;
template class Test<Fast_simplifier<CT>>;
