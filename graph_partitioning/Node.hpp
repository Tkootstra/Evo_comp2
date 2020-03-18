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
    void initializeNode(int indexLoc, int nCon, const std::list<int> ConLocations);
    void flipPartition();

    // Variables
    std::list<int> ConnectionLocations;
    int gain, indexLocation, numberOfConnections, belongsToWhichPartition;
    bool isFixed;

    // For comparing Nodes
    bool operator == (const Node& s) const { return indexLocation == s.indexLocation; }
    bool operator != (const Node& s) const { return !operator==(s); }
    
    Node() {}
};

#endif 

void Node::initializeNode(int indexLoc, int nCon, const std::list<int> ConLocations)
{
    indexLocation = indexLoc;
    ConnectionLocations = ConLocations;
    numberOfConnections = nCon;
    isFixed = false;
    gain = 0;
}

void Node::flipPartition()
{
    if (belongsToWhichPartition == 1) belongsToWhichPartition = 0;
    else belongsToWhichPartition = 1;
}
