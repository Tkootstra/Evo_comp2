// Node.hpp

#include <iostream>
#include <list>

#ifndef Node_hpp
#define Node_hpp

using namespace std;

// Node class
class Node
{
    public:
    // Functions
    void initializeNode(int indexLoc, int nCon, std::list<int> ConLocations);
    void setPartition(int partitionNumber);
    void flipPartition();

    // Variables
    std::list<int> ConnectionLocations;
    int indexLocation, numberOfConnections, belongsToWhichPartition;

    // For comparing Nodes
    bool operator == (const Node& s) const { return indexLocation == s.indexLocation; }
    bool operator != (const Node& s) const { return !operator==(s); }
    
    Node()
    {   
    }
};

#endif 

void Node::initializeNode(int indexLoc, int nCon, std::list<int> ConLocations)
{
    indexLocation = indexLoc;
    ConnectionLocations = ConLocations;
    numberOfConnections = nCon;
}

void Node::setPartition(int partitionNumber)
{
    belongsToWhichPartition = partitionNumber;
}

void Node::flipPartition()
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
