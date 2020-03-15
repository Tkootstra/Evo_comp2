
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
#include <chrono>
#include <time.h>

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
        Node newNode;
        newNode.initializeNode(vertexNumber, nConnections, otherLocs);
        nodeList.push_back(newNode);
    }

    return nodeList;
}

std::vector<int> makeRandomSolution(int stringLength)
{
    bool equal = false;
    std::vector<int> solution(stringLength);

    srand(time(NULL));

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

std::vector<int> perturbSolution(std::vector<int> solution, float ratio)
{
    std::vector<int> tempSolution = solution;
    int amountToPerturb = (int)(ratio * solution.size());
    
    for (size_t i = 0; i < amountToPerturb; i++)
    {
        int mutationLoc1 = rand() % solution.size();
        int toMutate1 = solution[mutationLoc1];

        int mutationLoc2 = rand() % solution.size();
        int toMutate2 = solution[mutationLoc2];

        // If bit-values are the same, try again (must be balanced)
        while (toMutate2 == toMutate1)
        {
            mutationLoc2 = rand() % solution.size();
            toMutate2 = solution[mutationLoc2];
        }

        if (toMutate1 == 0)
        {
            tempSolution[mutationLoc1] = 1;
            tempSolution[mutationLoc2] = 0;
        }
        else
        {
            tempSolution[mutationLoc1] = 0;
            tempSolution[mutationLoc2] = 1;
        }
    }
    
    return tempSolution;
}

Bucket computeGain(Graph graph, Bucket currentBucket)
{ //TODO: deze methode moet bij de bucket class horen <- Why?
    std::list<int> fixedNodes;
    int whichPart;
    bool notInFixed;
    int newCutState;

    int gain;
    Graph tempGraph;
    Node currentNode, current;

    // Count the current score
    graph.countConnections(0); 
    int originalCutState = graph.cutStatePartition0;
    currentBucket.currentSolution = originalCutState;
    
    for (size_t i = 0; i < graph.Nodes.size(); i++)
    {
        currentNode = graph.Nodes[i];
        fixedNodes = currentBucket.fixedNodes;

        // Check if node not in fixedNodes, otherwise skip – no need to update gain
        notInFixed = std::find(fixedNodes.begin(), fixedNodes.end(), currentNode.indexLocation) == fixedNodes.end();

        if (notInFixed)
        {
            tempGraph = graph;
            current = tempGraph.Nodes[i];
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
    Graph tempGraph;
    Node currentNode, current;

    // Count current score
    graph.countConnections(0); 
    int originalCutState = graph.cutStatePartition0;
    currentBucket.currentSolution = originalCutState;

    // Only loop through necessary connections
    for (auto const &i: nodeConnections)
    {
        currentNode = graph.Nodes[i - 1];
        fixedNodes = currentBucket.fixedNodes;

        notInFixed = std::find(fixedNodes.begin(), fixedNodes.end(), currentNode.indexLocation) == fixedNodes.end();

        if (notInFixed)
        {
            tempGraph = graph;
            current = tempGraph.Nodes[i - 1];
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

    // Do until both buckets are empty:
        // Move node with highest gain from 0 to 1
            // Remove node from bucket (the node becomes LOCKED)
        // Do same from 1 to 0, otherwise solution isn't valid

        // Recompute the gain for all neighboring nodes of the nodes that have been moved

    int updateFunctions = 0;

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

        updateFunctions = updateFunctions + nodeToChange0.numberOfConnections;
        updateFunctions = updateFunctions + nodeToChange1.numberOfConnections;
        // std::cout << "0: " << nodeToChange0.numberOfConnections << "| 1: " << nodeToChange1.numberOfConnections << endl;

        score = results.currentSolution;
        
        // Update metadata
        b0max = results.bucket0maxPointer;
        b1max = results.bucket1maxPointer;
        b0size = results.bucket0Size;
        b1size = results.bucket1Size;
    }

    return score;

}

std::vector<vector<double> > multiStartLocalSearch(std::vector<Node> nodeList, int iterations)
{
    std::vector<vector<double> > combinedResults(iterations);
    std::vector<double> cR(2);
    std::vector<int> solution;
    
    int r;
    std::chrono::steady_clock::time_point begin, end;
    std::chrono::duration<double> dur;
    Graph graaf;

    for (size_t i = 0; i < iterations; i++)
    {
        begin = std::chrono::steady_clock::now();
        solution = makeRandomSolution(500);
        graaf.initializeGraph(nodeList, solution);

        // Run one local search
        r = singleFidMath(graaf);

        // Calculate elapsed time
        end = std::chrono::steady_clock::now();
        dur = (end - begin);

        cR[0] = r;
        cR[1] = dur.count();
        combinedResults[i] = cR;
        
        std::cout << "Iteration " << i + 1 << ": Score " << r << " | Time: " << cR[1] << "s." << endl;
    }

    return combinedResults;
}

std::vector<vector<double> > iterativeLocalSearch(std::vector<Node> nodeList, int iterations, float perturbationRatio)
{
    std::vector<vector<double> > combinedResults(iterations);
    std::vector<double> cR(2);
    std::vector<int> solution = makeRandomSolution(500);
    std::vector<int> tempSolution = solution;
    
    int prevResult = 0;
    int thisResult;
    std::chrono::steady_clock::time_point begin, end;
    std::chrono::duration<double> dur;
    Graph graaf;

    for (size_t i = 0; i < iterations; i++)
    {
        begin = std::chrono::steady_clock::now();

        tempSolution = perturbSolution(solution, perturbationRatio);
        graaf.initializeGraph(nodeList, solution);

        // Run one local search
        thisResult = singleFidMath(graaf);

        if (thisResult > prevResult)
        {
            solution = tempSolution;
            prevResult = thisResult;
        }

        // Calculate elapsed time
        end = std::chrono::steady_clock::now();
        dur = (end - begin);

        cR[0] = thisResult;
        cR[1] = dur.count();
        combinedResults[i] = cR;
        
        std::cout << "Iteration " << i + 1 << ": Score " << thisResult << " | Time: " << cR[1] << "s." << endl;
    }

    return combinedResults;
}

void writeToFile(std::vector<vector<double> > results, std::string fileName)
{
    std::ofstream output_file(fileName);
    output_file << "Score Time" << endl;

    for (auto const &r: results)
    {
        for (auto const &v: r)
        {
            output_file << v << " ";
        }
        
        output_file << endl;
    }

    output_file.close();
}

int main()
{   // https://www.codeproject.com/Articles/1271904/Programming-Concurrency-in-Cplusplus-Part-1
    
    // parse nodes from txt file
    std::vector<Node> nodeList = parseGraph();
    
    // Run MLS
    std::vector<vector<double> > results = multiStartLocalSearch(nodeList, 10000);
    writeToFile(results, "MLS.txt");

    // Run ILS
    // std::vector<vector<double> > results = iterativeLocalSearch(nodeList, 100, 0.05);
    // writeToFile(results, "ILS.txt");

} 

 