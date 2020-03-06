// kwk3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

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
class Node
{
    public:
    int indexLocation;
    std::list<int> ConnectionLocations;
    int numberOfConnections;
    int belongsToWhichPartition;
    Node(int indexLoc, int nCon, std::list<int> ConLocations)
    {   
        indexLocation = indexLoc;
        ConnectionLocations = ConLocations;
        numberOfConnections = nCon;
    }

    void setPartition(int partitionNumber)
    {
        belongsToWhichPartition = partitionNumber;
    }

    void flipPartition()
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
};

class Graph
{
    public:
    int cutStatePartition0;  // The amount of connections between partition 0 and 1

    std::vector<Node> Nodes;
    std::vector<int> gains;
    
    Graph(std::vector<Node> nodeList, vector<int> solution)
    {
        for (size_t i = 0; i < nodeList.size(); i++)
        {
            int idx = solution[i];
            nodeList[i].setPartition(idx);  
        }

        Nodes = nodeList;      
    }

    int getPartitionSize(int partition)
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

    void countConnections(int partition)
    // for each node in the graph, if the node belongs to partition 0, 
    // check  all connections from this node and check if they belong to partition 1.
    // if so, 
    {
        int totalValue = 0;
        std::list<int> connections;

        for (Node checkNode: Nodes)
        {
            if (checkNode.belongsToWhichPartition == partition)
            {
                connections = checkNode.ConnectionLocations;
                for (int index:connections)
                {
                    if (Nodes[index].belongsToWhichPartition != partition)
                    {
                        totalValue ++; 
                    }
                }
            } 
        }
        // return totalValue;
        cutStatePartition0 = totalValue;
    }

    // void countConnections()
    // {
    //     cutStatePartition0 = countCon(0);
    // }
 
};

class Bucket
{
    public:
    std::map<int, std::list<Node> > bucket0;
    std::map<int, std::list<Node> > bucket1;
    int bucket0maxPointer;
    int bucket1maxPointer;
    int bucket0Size;
    int bucket1Size;
    std::list<int> fixedNodes;
    int currentSolution;

    Bucket(std::map<int, std::list<Node> > b0, std::map<int, std::list<Node> > b1, std::list<int> fN, Graph g)
    {
        bucket0 = b0;
        bucket1 = b1;
        bucket0maxPointer = -999;
        bucket1maxPointer = -999;
        bucket0Size = g.getPartitionSize(0);
        bucket1Size = 500 - bucket0Size;

        fixedNodes = fN;
        
    }

    void addToBucket(int partition, int key, Node item)
    {
        if (partition == 0)
        {
            bucket0[key].push_back(item);

            if (key > bucket0maxPointer)
            {
                bucket0maxPointer = key;
            }

        }
        else
        {
            bucket1[key].push_back(item);

            if (key > bucket1maxPointer)
            {
                bucket1maxPointer = key;
            }

        }
    }

    Node popFromBucketKey(int partition, int key)
    {
        // TODO: Buckets now get longer with duplicate elements. When a node's gain is updated, it should be just moved to another key
        if (partition == 0)
        {
            if (bucket0[key].size() > 0)
            {
                Node item = bucket0[key].front();
                bucket0[key].pop_front();
                bucket0Size--;

                fixedNodes.push_back(item.indexLocation);
                return item;
            }
            else
            {
                // std::cout << "decreasing from " << key << " to " << key - 1 << endl;
                bucket0maxPointer--;
                return popFromBucketKey(partition, key - 1);
            }
            
        }
        else 
        {
            if (bucket1[key].size() > 0)
            {
                Node item = bucket1[key].front();
                bucket1[key].pop_front();
                bucket1Size--;

                fixedNodes.push_back(item.indexLocation);
                return item;
            }
            else
            {
                // std::cout << "decreasing from " << key << " to " << key - 1 << endl;
                bucket1maxPointer--;
                return popFromBucketKey(partition, key - 1);
            }
            
        }

    }

};

Bucket computeGain(Graph graph, Bucket currentBucket)
{ //TODO: deze methode moet bij de bucket class horen
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

            // std::cout << "difference: ";
            // std::cout << gain << endl;  

        }    
    }

    return currentBucket;    
}

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

    int nodeToChangeIndex;

    // std::cout << "b0size: ";
    // std::cout << b0size << endl;

    // std::cout << "Current solution: ";
    // std::cout << score << endl;

    // Do until both buckets are empty:
        // If bucket0.size() > bucket1.size():              (And the other way around)
            // Move node with highest gain from 0 to 1
            // Remove node from bucket (the node becomes LOCKED)

        // Recompute the gain for all FREE cells

    while (b0size > 0 && b1size > 0) 
    {
        if (b0size >= b1size) 
        {
            Node nodeToChange = results.popFromBucketKey(0, b0max);
            nodeToChangeIndex = nodeToChange.indexLocation;
        }
        else
        {
            Node nodeToChange = results.popFromBucketKey(1, b1max);
            nodeToChangeIndex = nodeToChange.indexLocation;
        }

        // Flip bit
        // Graph tempGraph = g;
        // tempGraph.Nodes[nodeToChangeIndex].flipPartition();
        // g = tempGraph;

        g.Nodes[nodeToChangeIndex].flipPartition();
        results = computeGain(g, results);
        score = results.currentSolution;
        
        // Update metadata
        b0max = results.bucket0maxPointer;
        b1max = results.bucket1maxPointer;
        b0size = results.bucket0Size;
        b1size = results.bucket1Size;

        std::cout << "b0size: " << b0size << endl;
        std::cout << "Current solution: " << score << endl;
    }

}

int main()
{   // parse nodes from txt file
    std::vector<Node> NodeList = parseGraph();
    std::vector<int> solution = makeRandomSolution(500);

    Graph graaf = Graph(NodeList, solution);
    singleFidMath(graaf);
    

} 

 