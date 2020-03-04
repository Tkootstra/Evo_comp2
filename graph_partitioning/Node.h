#pragma once

#include <iostream>
#include <list>




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
