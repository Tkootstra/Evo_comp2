// Edge.hpp

#include <list>

#ifndef Edge_hpp
#define Edge_hpp

using namespace std;

// Node class
class Edge
{
    public:
    // Functions
    void initializeEdge(const int ind);

    // Variables
    std::list<std::pair<int, int> > ConnectionLocations;
    int index;
    bool cutState;

    // For comparing Nodes
    bool operator == (const Edge& s) const { return index == s.index; }
    bool operator != (const Edge& s) const { return !operator==(s); }
    
    Edge()
    {   
    }
};

#endif 

void Edge::initializeEdge(const int ind)
{
    index = ind;
}

