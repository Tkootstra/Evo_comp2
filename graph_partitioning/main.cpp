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
#include <tuple>

// c++ -std=c++17 -stdlib=libc++ main.cpp -o prog 
// ./prog


using namespace std;

// Node class

class Node
{
public:
    int indexLocation;
    std::list<int> ConnectionLocations;
    int numberOfConnections;
    int belongsToWhichPartition;
    Node(int indexLoc, int nCon, std::list<int> ConLocations)
    {   
        indexLocation = indexLoc;
        ConnectionLocations = ConLocations;
        numberOfConnections = nCon;
    }

    void setPartition(int partitionNumber)
    {
        belongsToWhichPartition = partitionNumber;
    }

    void flipPartition()
    {
        if (belongsToWhichPartition == 1)
        {
            belongsToWhichPartition = 0;
        }
        else
        {
            belongsToWhichPartition = 1;
        }
        
        
    }
};



class Graph
{
    public:
    int baseGainPartition0;
    int baseGainPartition1;

    std::vector<Node> Nodes;
    std::vector<int> gains;
    Graph(std::vector<Node> nodeList, vector<int> solution)
    {
        for (size_t i = 0; i < nodeList.size(); i++)
        {
            int idx = solution[i];
            nodeList[i].setPartition(idx);  
        }

        Nodes = nodeList;
        
        
    }
    int countCon(int partition)
    // for each node in the graph, if the node belongs to partition 0, 
    // check  all connections from this node and check if they belong to partition 1.
    // if so, 
    {
        
        int totalValue = 0;
        for (Node checkNode: Nodes)
        {
            if (checkNode.belongsToWhichPartition == partition)
            {
                std::list<int> connections = checkNode.ConnectionLocations;
                for (int index:connections)
                {
                    if (Nodes[index].belongsToWhichPartition != partition)
                    {
                        totalValue ++; 
                        // std::cout << totalValue << endl;
                    }
                }
            } 
        }
        return totalValue;
    }

    void countConnections()
    {
        baseGainPartition0 = countCon(0);
        baseGainPartition1 = countCon(1);
    }

    
    
    
};






std::tuple<std::map<int, std::list<int> >, std::map<int, std::list<int> >, int, int>computeGain(Graph graph)
{
    std::vector<int> gains;
    std::cout << "begin" << endl;
    graph.countConnections();
    std::cout << "begin2" << endl;
    int originalBaseGain0 = graph.baseGainPartition0;
    int originalBaseGain1 = graph.baseGainPartition1;
    int maxCon = 0;
    
    // std::list<std::map<int, std::list<int> > > results;
    std::map<int, std::list<int> > bucketA;
    std::map<int, std::list<int> > bucketB;
    int bucketAmaxPointer = -999999;
    int bucketBmaxPointer = -999999;
    

    
    // elke bit index flippen
    // voor elke graph die hieruit komt doe je count connections en compare je die met de base gain. 
    //dit moet je opslaan per node in een list
    for (size_t i = 0; i < graph.Nodes.size(); i++)
    {
        int diff;
        Node current = graph.Nodes[i];
        int whichPart = current.belongsToWhichPartition;
        graph.Nodes[i].flipPartition();
        graph.countConnections();
        int gain0 = graph.baseGainPartition0;
        int gain1 = graph.baseGainPartition1;
        if (whichPart == 0)
        {
            int diff = originalBaseGain0 - gain0;

        }
        else
        {
            int diff = originalBaseGain1 - gain1;
        }
        
         // welke kant moet dit op?
        std::cout << "difference" << endl;
        std::cout << diff << endl;
      

        if (whichPart == 0)
        {
            if (diff > bucketAmaxPointer)
            {
                bucketAmaxPointer = diff;
            }
            bucketA[diff].push_back(current.indexLocation);
        }
        else
        {
            if (diff > bucketBmaxPointer)
            {
                bucketBmaxPointer = diff;
            }
            bucketB[diff].push_back(current.indexLocation);
        } 
    }
    //TODO: bug: difference is altijd 1 tussen base graph en alle 500 opties 
    std::tuple<std::map<int, std::list<int> >, std::map<int, std::list<int> >, int, int> results \
    (bucketA, bucketB, bucketAmaxPointer, bucketBmaxPointer);
    return results;
    
}




std::list<vector<int> > makeMultipleRandomSolution(int bitLength, int amount)
{
    //TODO: moeten gelijk aantal 0 en 1 zijn
    std::list<vector<int> > allSolutions;
    for (size_t i = 0; i < amount; i++)
    {
        std::vector<int> solution(bitLength);
         for (size_t j = 0; j < bitLength; j++)
         {
             int bit = rand() % 2;
             solution[j] = bit;
         }
        allSolutions.push_back(solution);
    }
    return allSolutions;
}







std::vector<Node> parseGraph()
{
    std::ifstream dataFile("Graph500.txt");
    std::string line;
    std::vector<Node> nodeList;
    int maxcon = 0;

    while (std::getline(dataFile, line))
    {
        
        // parse all info from .txt file line by line
        std::istringstream iss(line);
        std::vector<std::string> singleSplittedLine((std::istream_iterator<std::string>(iss)),
            std::istream_iterator<std::string>());

        int vertexNumber = std::stoi(singleSplittedLine[0]);
        int nConnections = std::stoi(singleSplittedLine[2]);
        if (nConnections > maxcon)
        {
            maxcon = nConnections;
        }
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
    // std::cout << "Parsed data from txt file. Total amount of nodes: \n";
    // std::cout << nodeList.size() << "\n";
    // std::cout << "max number of connections is" << endl;
    // std::cout << maxcon << endl;
    return nodeList;
}

int main()
{   // parse nodes from txt file
    std::vector<Node> NodeList = parseGraph();
    std::list<vector<int> > solution = makeMultipleRandomSolution(500, 1);
    // std::cout << "kek" << endl;
    vector<int> sol = solution.front();
    // std::cout << "kek" << endl;
    Graph graaf = Graph(NodeList, sol);
    // std::cout << "kek" << endl;
    std::map<int, std::list<int> > bucketA;
    std::map<int, std::list<int> > bucketB;
    int bucketAmaxPointer;
    int bucketBmaxPointer;

    std::tuple<std::map<int, std::list<int> >, std::map<int, \
    std::list<int> >, int, int>results = computeGain(graaf);
    std::tie (bucketA, bucketB, bucketAmaxPointer, bucketBmaxPointer) = results;
    std::cout << bucketAmaxPointer << endl;
}

