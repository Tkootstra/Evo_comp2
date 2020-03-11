
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

    return nodeList;
}

std::vector<int> makeRandomSolution(int stringLength)
{
    bool equal = false;
    std::vector<int> solution(stringLength);
    while (!equal)
    {
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

std::list<vector<int> > makeMultipleRandomSolutions(int stringLength, int amount)
{
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

Bucket updateGain(Graph graph, Bucket currentBucket, std::list<int> nodeConnections)
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

    // Only loop through necessary connections
    for (auto const &i: nodeConnections)
    {
        Node current = graph.Nodes[i - 1];
        fixedNodes = currentBucket.fixedNodes;

        notInFixed = std::find(fixedNodes.begin(), fixedNodes.end(), current.indexLocation) == fixedNodes.end();

        if (notInFixed)
        {
            Graph tempGraph = graph;
            Node current = tempGraph.Nodes[i - 1];
            whichPart = current.belongsToWhichPartition;
            
            tempGraph.Nodes[current.indexLocation].flipPartition();
            tempGraph.countConnections(0);
            
            newCutState = tempGraph.cutStatePartition0;
            
            gain = originalCutState - newCutState;
            
            currentBucket.updateBucket(whichPart, gain, current);
        } 
    }
    return currentBucket;
}

int singleFidMath(Graph g)
{   
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

    // std::cout << "b0size: " << b0size << endl;
    // std::cout << "Current solution: " << score << endl;

    // Do until both buckets are empty:
        // Move node with highest gain from 0 to 1
            // Remove node from bucket (the node becomes LOCKED)
        // Do same from 1 to 0, otherwise solution isn't valid

        // Recompute the gain for all neighboring nodes of the nodes that have been moved

    while (b0size > 0 && b1size > 0) 
    {
        // Move node from 0 to 1
        Node nodeToChange0 = results.popFromBucketKey(0, b0max);
        nodeToChangeIndex0 = nodeToChange0.indexLocation;
        g.Nodes[nodeToChangeIndex0].flipPartition();

        // Move node from 1 to 0
        Node nodeToChange1 = results.popFromBucketKey(1, b1max);
        nodeToChangeIndex1 = nodeToChange1.indexLocation;
        g.Nodes[nodeToChangeIndex1].flipPartition();

        // We now have a new valid partition; update gains for neighbors
        results = updateGain(g, results, nodeToChange0.ConnectionLocations);
        results = updateGain(g, results, nodeToChange1.ConnectionLocations);

        score = results.currentSolution;
        
        // Update metadata
        b0max = results.bucket0maxPointer;
        b1max = results.bucket1maxPointer;
        b0size = results.bucket0Size;
        b1size = results.bucket1Size;

        // std::cout << "b0: " << b0size << " | b1: " << b1size << endl;
        // std::cout << "Current solution: " << score << endl;
        // std::cout << "Pointer0: " << b0max << " | Pointer1: " << b1max << endl << endl;
    }

    return score;

}

std::vector<int> multipleFidMath(std::vector<Node> nodeList, int iterations)
{
    std::vector<int> results(iterations);
    std::vector<int> solution;
    int r;

    for (size_t i = 0; i < iterations; i++)
    {
        solution = makeRandomSolution(500);
        Graph graaf = Graph(nodeList, solution);

        r = singleFidMath(graaf);
        results[i] = r;

        std::cout << "Iteration " << i + 1 << ", Score " << r << endl;
    }

    return results;
}

int main()
{   // parse nodes from txt file
    std::vector<Node> NodeList = parseGraph();
    
    // Run multiple solutions
    std::vector<int> results = multipleFidMath(NodeList, 10);

} 

 