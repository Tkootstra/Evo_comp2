// graph_partitioning.cpp : Defines the entry point for the application.
//

#include "graph_partitioning.h"
#include <iostream>
#include <list>
#include <vector>
#include <fstream>
#include <string>


using namespace std;



int main()
{
	
	cout << "Hello CMake." << endl;
	std::ifstream dataFile("Graph500.txt");
	std::string line;
	while (std::getline(dataFile, line))
	{
		cout << line << endl;
	}



	// hier data inladen -> Graph500.txt
	return 0;
}


class Node
{
public:
	int indexLocation;
	std::list<int> ConnectionLocations;
	Node(int indexLoc, std::list<int> ConLocations)
	{
		indexLocation = indexLoc;
		ConnectionLocations = ConLocations;
	}
};

class Solution 
{
	
public:
	std::vector<int> valueVector;
	int values[];
	

	Solution(int values[])
	{	
		values = values;
		valueVector.insert(valueVector.begin(), values, values + (sizeof(values) / sizeof(*values)));
	}
};








