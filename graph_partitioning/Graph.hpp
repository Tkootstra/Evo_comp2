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
    void initializeGraph(std::vector<Node> nodeList, const vector<int> solution);
    int countConnections() const;
    int countSingleCellConnections(const int nodeIndex) const;

    // Variables
    std::vector<Node> Nodes;
    
    Graph() {}
};

#endif

void Graph::initializeGraph(std::vector<Node> nodeList, const vector<int> solution)
{
    for (size_t i = 0; i < nodeList.size(); i++)
    {
        // Set partition to value of that location in solution
        nodeList[i].belongsToWhichPartition = solution[i];  
    }

    Nodes = nodeList;   
}

int Graph::countConnections() const
{ // for each node in the graph, if the node belongs to partition 0, 
  // check  all connections from this node and check if they belong to partition 1.
    int totalValue = 0;

    for (Node const &checkNode: Nodes)
    {
        if (checkNode.belongsToWhichPartition == 0)
        {
            for (int const &index: checkNode.ConnectionLocations)
            {
                if (Nodes[index].belongsToWhichPartition == 1)
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

    // Loop through node's connections and count which one are in a a different partition
    for (int const &index: checkNode.ConnectionLocations)
    {
        if (Nodes[index - 1].belongsToWhichPartition != checkNode.belongsToWhichPartition)
        {
            totalValue++;
        }
    }

    return totalValue;

}
 
