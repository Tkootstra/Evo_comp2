// Graph.hpp

// Custom classes
#include "Node.hpp"

// Base libraries
#include <vector>
#include <list>

using namespace std;

#ifndef Graph_hpp
#define Graph_hpp

class Graph
{
    public:
    // Functions
    void initializeGraph(std::vector<Node> nodeList, vector<int> solution);
    int getPartitionSize(int partition);
    void countConnections(int partition);

    // Variables
    int cutStatePartition0;  // The amount of connections between partition 0 and 1

    std::vector<Node> Nodes;
    
    Graph()
    {   
    }

};

#endif

void Graph::initializeGraph(std::vector<Node> nodeList, vector<int> solution)
{
    for (size_t i = 0; i < nodeList.size(); i++)
    {
        int idx = solution[i];
        nodeList[i].setPartition(idx);  
    }

    Nodes = nodeList;   
}

int Graph::getPartitionSize(int partition)
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

void Graph::countConnections(int partition)
// for each node in the graph, if the node belongs to partition 0, 
// check  all connections from this node and check if they belong to partition 1.
{
    int totalValue = 0;
    std::list<int> connections;

    for (Node const &checkNode: Nodes)
    {
        if (checkNode.belongsToWhichPartition == partition)
        {
            connections = checkNode.ConnectionLocations;
            for (int const &index: connections)
            {
                if (Nodes[index].belongsToWhichPartition != partition)
                {
                    totalValue ++; 
                }
            }
        } 
    }
    // return totalValue;
    cutStatePartition0 = totalValue;
}
 
