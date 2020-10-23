#include <iostream>
#include <iomanip>
#include <vector>


/**
  *	Generates all combinations of the symmetryIndices list, stored in a list of lists. O(N*(2^N)), where N is the 
  * number of symmetric axes.
  */
std::vector<std::vector<int>> getReflectionIndices(const std::vector<int>& symmetryIndices)
{
	int N = symmetryIndices.size();
	int numSubLists = 1 << N;
	std::vector<std::vector<int>> reflectionLists;

	for(int i = 0; i < numSubLists; ++i)
	{
		std::vector<int> reflectionList = {};

		for(int j = 0; j < N; ++j)
		{
			if((i >> j) & 1)
			{
				reflectionList.push_back(symmetryIndices[j]);
			}
		}

		reflectionLists.push_back(reflectionList);
	}

	return reflectionLists;
}

int main(int argc, char* argv[])
{
	std::vector<int> symmetryIndices = {0, 2, 3};
	std::vector<std::vector<int>> reflectionLists = getReflectionIndices(symmetryIndices);

	for(int i = 0; i < reflectionLists.size(); ++i)
	{
		for(auto const& c : reflectionLists[i])
		{
			std::cout << c << ' ';
		}

		std::cout << std::endl;
	}
}