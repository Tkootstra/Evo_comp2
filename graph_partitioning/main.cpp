// kwk3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <iterator>
#include <random>
#include <map>
// #include "Node.h"

using namespace std;

// Node class
class Node
{
public:
    int indexLocation;
    std::list<int> ConnectionLocations;
    int numberOfConnections;
    Node(int indexLoc, int nCon, std::list<int> ConLocations)
    {   
        indexLocation = indexLoc;
        ConnectionLocations = ConLocations;
        numberOfConnections = nCon;
    }
};

class Solution 
{
	
public:
	std::vector<int> valueVector;
    int length;
	Solution(std::vector<int> values)
	{	
		valueVector = values;
        length = values.size();	
	}  
};


std::list<Solution> makeMultipleRandomSolution(int bitLength, int amount)
{
    std::list<Solution> allSolutions;
    for (size_t i = 0; i < amount; i++)
    {
        std::vector<int> solution(bitLength);
         for (size_t j = 0; j < bitLength; j++)
         {
             int bit = rand() % 2;
             solution[j] = bit;
         }
        allSolutions.push_back(Solution(solution));
    }
    return allSolutions;
}






int main()
{
   
    std::ifstream dataFile("Graph500.txt");
    std::string line;
    std::list<Node> nodeList;
    

    while (std::getline(dataFile, line))
    {
        
        // parse all info from .txt file line by line
        std::istringstream iss(line);
        std::vector<std::string> singleSplittedLine((std::istream_iterator<std::string>(iss)),
            std::istream_iterator<std::string>());

        int vertexNumber = std::stoi(singleSplittedLine[0]);
        int nConnections = std::stoi(singleSplittedLine[2]);
        std::list<int> otherLocs;

        // parse all other indice locations
        for (size_t i = 3; i < singleSplittedLine.size(); i++)
        {   
            int loc = std::stoi(singleSplittedLine[i]);
            otherLocs.push_back(loc);    
        }
        // make new nodes for each input line
        Node newNode = Node(vertexNumber, nConnections, otherLocs);
        nodeList.push_back(newNode);
        
    }
    std::cout << "Parsed data from txt file. Total amount of nodes: \n";
    std::cout << nodeList.size() << "\n";

 

}
