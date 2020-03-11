
// Custom classes
// #include "Node.hpp"   <- Don't include because they're included in Bucket.hpp
// #incluse "Graph.hpp"  <-|
#include "Bucket.hpp"

// Base libraries
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <iterator>
#include <random>
#include <map>
#include <algorithm>

// c++ -std=c++17 -stdlib=libc++ main.cpp -o prog 
// ./prog

using namespace std;

// Node class
// class Node
// {
//     public:
//     int indexLocation;
//     std::list<int> ConnectionLocations;
//     int numberOfConnections;
//     int belongsToWhichPartition;

//     bool operator == (const Node& s) const { return indexLocation == s.indexLocation; }
//     bool operator != (const Node& s) const { return !operator==(s); }
    
//     Node(int indexLoc, int nCon, std::list<int> ConLocations)
//     {   
//         indexLocation = indexLoc;
//         ConnectionLocations = ConLocations;
//         numberOfConnections = nCon;
//     }

//     void setPartition(int partitionNumber)
//     {
//         belongsToWhichPartition = partitionNumber;
//     }

//     void flipPartition()
//     {
//         if (belongsToWhichPartition == 1)
//         {
//             belongsToWhichPartition = 0;
//         }
//         else
//         {
//             belongsToWhichPartition = 1;
//         }  
        
//     }
// };

// class Graph
// {
//     public:
//     int cutStatePartition0;  // The amount of connections between partition 0 and 1

//     std::vector<Node> Nodes;
//     std::vector<int> gains;
    
//     Graph(std::vector<Node> nodeList, vector<int> solution)
//     {
//         for (size_t i = 0; i < nodeList.size(); i++)
//         {
//             int idx = solution[i];
//             nodeList[i].setPartition(idx);  
//         }

//         Nodes = nodeList;      
//     }

//     int getPartitionSize(int partition)
//     {
//         int partitionSize = 0;

//         for (size_t i = 0; i < Nodes.size(); i++)
//         {
//             if (Nodes[i].belongsToWhichPartition == partition) 
//             {
//                 partitionSize++;
//             }
//         }

//         return partitionSize;
//     }

//     void countConnections(int partition)
//     // for each node in the graph, if the node belongs to partition 0, 
//     // check  all connections from this node and check if they belong to partition 1.
//     // if so, 
//     {
//         int totalValue = 0;
//         std::list<int> connections;

//         for (Node checkNode: Nodes)
//         {
//             if (checkNode.belongsToWhichPartition == partition)
//             {
//                 connections = checkNode.ConnectionLocations;
//                 for (int index:connections)
//                 {
//                     if (Nodes[index].belongsToWhichPartition != partition)
//                     {
//                         totalValue ++; 
//                     }
//                 }
//             } 
//         }
//         // return totalValue;
//         cutStatePartition0 = totalValue;
//     }

//     // void countConnections()
//     // {
//     //     cutStatePartition0 = countCon(0);
//     // }
 
// };

// class Bucket
// {
//     public:
//     std::map<int, std::list<Node> > bucket0;
//     std::map<int, std::list<Node> > bucket1;
//     int bucket0maxPointer;
//     int bucket1maxPointer;
//     int bucket0Size;
//     int bucket1Size;
    
//     std::list<int> fixedNodes;
//     int currentSolution;

//     Bucket(std::map<int, std::list<Node> > b0, std::map<int, std::list<Node> > b1, std::list<int> fN, Graph g)
//     {
//         bucket0 = b0;
//         bucket1 = b1;
//         bucket0maxPointer = -999;
//         bucket1maxPointer = -999;
//         bucket0Size = g.getPartitionSize(0);
//         bucket1Size = 500 - bucket0Size;

//         fixedNodes = fN;
//     }

//     void addToBucket(int partition, int key, Node item)
//     {
//         if (partition == 0)
//         {
//             bucket0[key].push_back(item);

//             if (key > bucket0maxPointer)
//             {
//                 bucket0maxPointer = key;
//             }
//         }
//         else
//         {
//             bucket1[key].push_back(item);

//             if (key > bucket1maxPointer)
//             {
//                 bucket1maxPointer = key;
//             }
//         }
//     }

//     void updateBucket(int partition, int key, Node item)
//     {
//         // cout << "Updating bucket... ";
//         // 'Yank' item from one list to another
//         if (partition == 0)
//         {
//             for (auto &k: bucket0)
//             {
//                 // k.second gets bucket0[k] list
//                 std::list<Node> nodeList = k.second;

//                 for (auto const &i: nodeList)
//                 {
//                     if (i == item)
//                     {
//                         // If nodes match, remove from this bucket and add to other bucket
//                         nodeList.remove(i);
//                         addToBucket(partition, key, item);
//                         // std::cout << "now in bucket " << key << endl;
//                         return;
//                     }
//                 }
//             }
//         }
//         else
//         {
//             for (auto &k: bucket1)
//             {
//                 std::list<Node> nodeList = k.second;

//                 for (auto const &i: nodeList)
//                 {
//                     if (i == item)
//                     {
//                         nodeList.remove(i);
//                         addToBucket(partition, key, item);
//                         // std::cout << "now in bucket " << key << endl;
//                         return;
//                     }
//                 }
//             }  
//         }
//     }

//     Node popFromBucketKey(int partition, int key)
//     {
//         if (partition == 0)
//         {
//             if (bucket0[key].size() > 0)
//             {
//                 // Remove from bucket array
//                 Node item = bucket0[key].front();
//                 bucket0[key].pop_front();
//                 bucket0Size--;

//                 // Push to fixedNodes list
//                 fixedNodes.push_back(item.indexLocation);
//                 return item;
//             }
//             else
//             {
//                 // If bucket is empty, decrease maxPointer and run again
//                 bucket0maxPointer--;
//                 return popFromBucketKey(partition, key - 1);
//             } 
//         }
//         else 
//         {
//             if (bucket1[key].size() > 0)
//             {
//                 Node item = bucket1[key].front();
//                 bucket1[key].pop_front();
//                 bucket1Size--;

//                 fixedNodes.push_back(item.indexLocation);
//                 return item;
//             }
//             else
//             {
//                 bucket1maxPointer--;
//                 return popFromBucketKey(partition, key - 1);
//             } 
//         }
//     }
// };

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

std::vector<int> makeRandomSolution(int stringLength)
{
    bool equal = false;
    std::vector<int> solution(stringLength);
    while (!equal)
    {
        // std::cout << "not equal" << std::endl;
        for (size_t j = 0; j < stringLength; j++)
        {
            int bit = rand() % 2;
            solution[j] = bit;
        }
        int sum = 0;
        for (auto& n:solution) sum += n;
        if (sum == (stringLength/2)) {equal = true;}
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

Bucket computeGain(Graph graph, Bucket currentBucket)
{ //TODO: deze methode moet bij de bucket class horen <- Why?
    std::list<int> fixedNodes;
    int whichPart;
    bool notInFixed;
    int newCutState;

    int gain;

    graph.countConnections(0); 
    int originalCutState = graph.cutStatePartition0;
    currentBucket.currentSolution = originalCutState;
    
    for (size_t i = 0; i < graph.Nodes.size(); i++)
    {
        Node current = graph.Nodes[i];
        fixedNodes = currentBucket.fixedNodes;

        // Check if node not in fixedNodes, otherwise skip – no need to update gain
        notInFixed = std::find(fixedNodes.begin(), fixedNodes.end(), current.indexLocation) == fixedNodes.end();

        if (notInFixed)
        {
            Graph tempGraph = graph;
            Node current = tempGraph.Nodes[i];
            whichPart = current.belongsToWhichPartition;
            
            tempGraph.Nodes[current.indexLocation].flipPartition();
            tempGraph.countConnections(0);
            
            newCutState = tempGraph.cutStatePartition0;
            
            gain = originalCutState - newCutState;
            
            currentBucket.addToBucket(whichPart, gain, current);
        }    
    }

    return currentBucket;    
}

Bucket updateGain(Graph graph, Bucket currentBucket, std::list<int> connections)
{   // Do the same as computeGain but only for certain nodes and run 'updateBucket'
    // instead of 'addToBucket'
    std::list<int> fixedNodes;
    int whichPart;
    bool notInFixed;
    int newCutState;

    int gain;

    graph.countConnections(0); 
    int originalCutState = graph.cutStatePartition0;
    currentBucket.currentSolution = originalCutState;

    // std::cout << "Updating gain... ";
    // Only loop through necessary connections
    for (int i: connections)
    {
        Node current = graph.Nodes[i];
        fixedNodes = currentBucket.fixedNodes;

        notInFixed = std::find(fixedNodes.begin(), fixedNodes.end(), current.indexLocation) == fixedNodes.end();

        if (notInFixed)
        {
            Graph tempGraph = graph;
            Node current = tempGraph.Nodes[i];
            whichPart = current.belongsToWhichPartition;
            
            tempGraph.Nodes[current.indexLocation].flipPartition();
            tempGraph.countConnections(0);
            
            newCutState = tempGraph.cutStatePartition0;
            
            gain = originalCutState - newCutState;
            
            currentBucket.updateBucket(whichPart, gain, current);
        } 
    }

    // std::cout << "Updated gain!" << endl;
    return currentBucket;
}

void singleFidMath(Graph g)
{   
    // std::list<Node> freeNodes = getFreeNodes(g);
    std::list<int> fixedNodes;

    std::map<int, std::list<Node> > bucket0;
    std::map<int, std::list<Node> > bucket1;

    Bucket startBucket = Bucket(bucket0, bucket1, fixedNodes, g);

    Bucket results = computeGain(g, startBucket);
    
    // Make buckets and pointers
    int b0max = results.bucket0maxPointer;
    int b1max = results.bucket1maxPointer;
    int b0size = results.bucket0Size;
    int b1size = results.bucket1Size;

    int score = results.currentSolution;

    int nodeToChangeIndex0;
    int nodeToChangeIndex1;

    std::cout << "b0size: " << b0size << endl;
    std::cout << "Current solution: " << score << endl;

    // Do until both buckets are empty:
        // Move node with highest gain from 0 to 1
            // Remove node from bucket (the node becomes LOCKED)
        // Do same from 1 to 0, otherwise solution isn't valid

        // Recompute the gain for all neighboring nodes of the nodes that have been moved

    while (b0size > 0 && b1size > 0) 
    {
        // TODO: Segmentation fault occurs somewhere after score print 
        // but before first print in loop??
        // Seems to always occur around the size=230-240 mark

        std::cout << "Moving... ";
        // Move node from 0 to 1
        Node nodeToChange0 = results.popFromBucketKey(0, b0max);
        nodeToChangeIndex0 = nodeToChange0.indexLocation;
        g.Nodes[nodeToChangeIndex0].flipPartition();
        // std::cout << "1... ";

        // Move node from 1 to 0
        Node nodeToChange1 = results.popFromBucketKey(1, b1max);
        nodeToChangeIndex1 = nodeToChange1.indexLocation;
        g.Nodes[nodeToChangeIndex1].flipPartition();
        // std::cout << "2... ";

        // We now have a new valid partition; update gains for neighbors
        results = updateGain(g, results, nodeToChange0.ConnectionLocations);
        results = updateGain(g, results, nodeToChange1.ConnectionLocations);

        score = results.currentSolution;
        // std::cout << "Score!" << endl;
        
        // Update metadata
        b0max = results.bucket0maxPointer;
        b1max = results.bucket1maxPointer;
        b0size = results.bucket0Size;
        b1size = results.bucket1Size;

        std::cout << "b0: " << b0size << " | b1: " << b1size << endl;
        std::cout << "Current solution: " << score << endl;
        std::cout << "Pointer0: " << b0max << " | Pointer1: " << b1max << endl << endl;
    }

}

int main()
{   // parse nodes from txt file
    std::vector<Node> NodeList = parseGraph();
    std::vector<int> solution = makeRandomSolution(500);

    Graph graaf = Graph(NodeList, solution);
    singleFidMath(graaf);
} 

 