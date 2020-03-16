// Graph.hpp

// Custom classes
#include "test-Node.hpp"
#include "test-Edge.hpp"

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
    int countConnections(int partition) const;
    int countSingleCellConnections(const int nodeIndex) const;

    // Variables
    std::vector<Node> Nodes;
    
    Graph()
    {   
    }

};

#endif

void Graph::initializeGraph(std::vector<Node> nodeList, const vector<int> solution)
{
    for (size_t i = 0; i < nodeList.size(); i++)
    {
        int idx = solution[i];
        nodeList[i].setPartition(idx);  
    }

    Nodes = nodeList;   
}

int Graph::getPartitionSize(const int partition)
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

int Graph::countConnections(const int partition) const
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

    return totalValue;
}

int Graph::countSingleCellConnections(const int nodeIndex) const
{
    int totalValue = 0;
    Node checkNode = Nodes[nodeIndex - 1];
    int partition = checkNode.belongsToWhichPartition;
    std::list<int> connections = checkNode.ConnectionLocations;

    for (int const &index: connections)
    {
        if (Nodes[index - 1].belongsToWhichPartition != partition)
        {
            totalValue++;
        }
    }

    return totalValue;

}
 
