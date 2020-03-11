// Bucket.hpp

// Custom classes
#include "Graph.hpp"

// Base libraries
#include <list>
#include <vector>
#include <map>

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
};

#endif

void Bucket::addToBucket(int partition, int key, Node item)
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

void Bucket::updateBucket(int partition, int key, Node item)
{
    // cout << "Updating bucket... ";
    // 'Yank' item from one list to another
    if (partition == 0)
    {
        for (auto &k: bucket0)
        {
            // k.second gets bucket0[k] list
            std::list<Node> nodeList = k.second;

            for (auto const &i: nodeList)
            {
                if (i == item)
                {
                    // If nodes match, remove from this bucket and add to other bucket
                    nodeList.remove(i);
                    addToBucket(partition, key, item);
                    // std::cout << "now in bucket " << key << endl;
                    return;
                }
            }
        }
    }
    else
    {
        for (auto &k: bucket1)
        {
            std::list<Node> nodeList = k.second;

            for (auto const &i: nodeList)
            {
                if (i == item)
                {
                    nodeList.remove(i);
                    addToBucket(partition, key, item);
                    // std::cout << "now in bucket " << key << endl;
                    return;
                }
            }
        }  
    }
}

Node Bucket::popFromBucketKey(int partition, int key)
{
    if (partition == 0)
    {
        if (bucket0[key].size() > 0)
        {
            // Remove from bucket array
            Node item = bucket0[key].front();
            bucket0[key].pop_front();
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
            bucket1maxPointer--;
            return popFromBucketKey(partition, key - 1);
        } 
    }
}