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
    int cutStatePartition0;  // The amount of connections between partition 0 and 1
    int cutStatePartition1;

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

    int getPartitionSize(int partition)
    {
        int partitionSize = 0;

        for (size_t i = 0; i < Nodes.size(); i++)
        {
            if (Nodes[i].belongsToWhichPartition == partition) 
            {
                partitionSize++;
            }
        }

        return partitionSize;
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
        cutStatePartition0 = countCon(0);
        cutStatePartition1 = countCon(1);
    }
 
};

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


std::tuple<std::map<int, std::list<int> >, std::map<int, std::list<int> >, int, int>computeGain(Graph graph)
{
    std::vector<int> gains;

    graph.countConnections(); 

    int originalCutState0 = graph.cutStatePartition0;
    int originalCutState1 = graph.cutStatePartition1;
    // int maxCon = 0;
    
    // std::list<std::map<int, std::list<int> > > results;
    std::map<int, std::list<int> > bucket0;
    std::map<int, std::list<int> > bucket1;
    int bucket0maxPointer = -999;
    int bucket1maxPointer = -999;
    
    // elke bit index flippen
    // voor elke graph die hieruit komt doe je count connections en compare je die met de base gain. 
    //dit moet je opslaan per node in een list
    for (size_t i = 0; i < graph.Nodes.size(); i++)
    {
        
        Graph tempGraph = graph;
        Node current = tempGraph.Nodes[i];
        int whichPart = current.belongsToWhichPartition;
        
        tempGraph.Nodes[i].flipPartition();
        tempGraph.countConnections();
        
        int gain0 = tempGraph.cutStatePartition0;
        int gain1 = tempGraph.cutStatePartition1;
        int diff;
        
        if (whichPart == 0)
        {
            diff = originalCutState0 - gain0;
            
            if (diff > bucket0maxPointer)
            {
                bucket0maxPointer = diff;
            }
            bucket0[diff].push_back(current.indexLocation);

        }
        else
        {
            diff = originalCutState1 - gain1;

            if (diff > bucket1maxPointer)
            {
                bucket1maxPointer = diff;
            }
            bucket1[diff].push_back(current.indexLocation);
        }
        
        std::cout << "difference: ";
        std::cout << diff << endl;      
    }

    //TODO: bug: difference is altijd 1 tussen base graph en alle 500 opties 
    std::tuple<std::map<int, std::list<int> >, std::map<int, std::list<int> >, int, int> results \
    (bucket0, bucket1, bucket0maxPointer, bucket1maxPointer);
    return results;
    
}

std::vector<int> makeRandomSolution(int stringLength)
{
    // Generate random indices
    std::vector<int> indices(stringLength / 2);
    for (size_t j = 0; j < indices.size(); j++)
    {
        int bit = rand() % indices.size();
        indices[j] = bit;
    }

    // Fill in the solution vector
    std::vector<int> solution(stringLength, 0);
    for (auto index: indices)
    {
        solution[index] = 1;
    }
    return solution;
}

std::vector<int> makeRandomSolutionOld(int stringLength)
{
    std::vector<int> solution(stringLength);
        for (size_t j = 0; j < stringLength; j++)
        {
            int bit = rand() % 2;
            solution[j] = bit;
        }
    return solution;
}


std::list<vector<int> > makeMultipleRandomSolution(int stringLength, int amount)
{
    //TODO: moeten gelijk aantal 0 en 1 zijn
    std::list<vector<int> > allSolutions;
    for (size_t i = 0; i < amount; i++)
    {
        allSolutions.push_back(makeRandomSolution(stringLength));
    }
    return allSolutions;
}

void singleFidMath(Graph g)
{
    // Compute gain
    std::tuple<std::map<int, std::list<int> >, std::map<int, \
    std::list<int> >, int, int> results;
    results = computeGain(g);
    
    // Make buckets and pointers
    std::map<int, std::list<int> > bucket0;
    std::map<int, std::list<int> > bucket1;
    int bucket0maxPointer;
    int bucket1maxPointer;

    // Unpack tuple
    std::tie (bucket0, bucket1, bucket0maxPointer, bucket1maxPointer) = results;
    
    std::cout << "Max pointer 0: ";
    std::cout << bucket0maxPointer << endl;

    std::cout << "Max pointer 1: ";
    std::cout << bucket1maxPointer << endl;

    // Do until both buckets are empty:
        // If bucket0.size() > bucket1.size():              (And the other way around)
            // Move node with highest gain from 0 to 1
            // Remove node from bucket (the node becomes LOCKED)

            // Recompute the gain for all FREE cells

    
}

int main()
{   // parse nodes from txt file
    std::vector<Node> NodeList = parseGraph();
    std::vector<int> solution = makeRandomSolution(500);

    Graph graaf = Graph(NodeList, solution);
    singleFidMath(graaf);
    

}

