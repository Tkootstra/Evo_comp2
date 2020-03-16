// Bucket.hpp

// Custom classes
#include "test-Graph.hpp"

// Base libraries
#include <list>
#include <vector>
#include <array>
#include <fstream>

using namespace std;

#ifndef Bucket_hpp
#define Bucket_hpp

class Bucket
{
    public:
    // Functions
    void addToBucket(int partition, int key, Node item);
    void updateBucket(int partition, int key, Node item);
    Node popFromBucketKey(int partition, int key);

    // Variables
    std::vector<std::list<Node> > bucket0, bucket1;
    std::list<int> fixedNodes;
    
    int bucket0maxPointer, bucket1maxPointer;
    int bucket0Size, bucket1Size;
    int currentSolution;

    Bucket(std::vector<std::list<Node> > b0, std::vector<std::list<Node> > b1, std::list<int> fN, Graph g)
    {
        bucket0 = b0;
        bucket0.resize(32);
        bucket1 = b1;
        bucket1.resize(32);
        bucket0maxPointer = -999;
        bucket1maxPointer = -999;
        bucket0Size = g.getPartitionSize(0);
        bucket1Size = 500 - bucket0Size;

        fixedNodes = fN;
    }
};

#endif

void Bucket::addToBucket(const int partition, const int key, const Node item)
{
    // std::cout << "add " << key << endl;
    if (partition == 0)
    {
        bucket0[key + 16].push_back(item);

        if (key > bucket0maxPointer)
        {
            bucket0maxPointer = key;
        }
    }
    else
    {
        bucket1[key + 16].push_back(item);

        if (key > bucket1maxPointer)
        {
            bucket1maxPointer = key;
        }
    }
}

void Bucket::updateBucket(const int partition, const int key, const Node item)
{
    // std::cout << "upd " << key << endl;
    // 'Yank' item from one list to another
    if (partition == 0)
    {
        std::list<Node> nodeList = bucket0[key + 16];
        nodeList.remove(item);
        addToBucket(partition, key, item);
    }
    else
    {
        std::list<Node> nodeList = bucket0[key + 16];
        nodeList.remove(item);
        addToBucket(partition, key, item); 
    }
}

Node Bucket::popFromBucketKey(const int partition, const int key)
{
    // std::cout << "pop " << key << endl;
    if (partition == 0)
    {
        if (bucket0[key + 16].size() > 0)
        {
            // Remove from bucket array
            Node item = bucket0[key + 16].front();
            bucket0[key + 16].remove(item);
            bucket0Size--;

            // Push to fixedNodes list
            fixedNodes.push_back(item.indexLocation);
            return item;
        }
        else
        {
            // If bucket is empty, decrease maxPointer and run again
            bucket0maxPointer--;
            return popFromBucketKey(partition, key - 1);
        } 
    }
    else 
    {
        if (bucket1[key + 16].size() > 0)
        {
            Node item = bucket1[key + 16].front();
            bucket1[key + 16].remove(item);
            bucket1Size--;

            fixedNodes.push_back(item.indexLocation);
            return item;
        }
        else
        {
            bucket1maxPointer--;
            return popFromBucketKey(partition, key - 1);
        } 
    }
}