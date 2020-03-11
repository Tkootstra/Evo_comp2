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
    void setPartition(int partitionNumber);
    void flipPartition();

    // Variables
    int indexLocation;
    std::list<int> ConnectionLocations;
    int numberOfConnections;
    int belongsToWhichPartition;

    // For comparing Nodes
    bool operator == (const Node& s) const { return indexLocation == s.indexLocation; }
    bool operator != (const Node& s) const { return !operator==(s); }
    
    Node(int indexLoc, int nCon, std::list<int> ConLocations)
    {   
        indexLocation = indexLoc;
        ConnectionLocations = ConLocations;
        numberOfConnections = nCon;
    }
};

#endif 

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
