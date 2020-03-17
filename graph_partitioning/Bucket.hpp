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
    Node popFromBucketKey(int partition);
    int gainSum() const;

    // Variables
    std::map<int, std::list<Node> > bucket0, bucket1;

    int bucket0maxPointer, bucket1maxPointer;
    int bucket0Size, bucket1Size;
    int currentSolution;

    Bucket(std::map<int, std::list<Node> > b0, std::map<int, std::list<Node> > b1, const Graph g)
    {
        bucket0 = b0;
        bucket1 = b1;
        bucket0maxPointer = -999;
        bucket1maxPointer = -999;
        bucket0Size = 250;
        bucket1Size = 250;
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
    // 'Yank' item from one list to another
    if (partition == 0)
    {
        // Loop through dict keys
        for (auto &k: bucket0)
        {
            // Loop through list items (k.second retrieves list)
            for (auto const &i: k.second)
            {
                if (i == item)
                {
                    // If nodes match, remove from this bucket and add to other bucket
                    k.second.remove(i);
                    addToBucket(partition, key, item);
                    return;
                }
            }
        }
    }
    else
    {
        for (auto &k: bucket1)
        {
            for (auto const &i: k.second)
            {
                if (i == item)
                {
                    k.second.remove(i);
                    addToBucket(partition, key, item);
                    return;
                }
            }
        }  
    }
}

Node Bucket::popFromBucketKey(int partition)
{
    // std::cout << bucket0maxPointer << " " << bucket1maxPointer << endl;
    if (partition == 0)
    {
        if (bucket0[bucket0maxPointer].size() > 0)
        {
            // Remove from bucket array
            Node item = bucket0[bucket0maxPointer].front();
            bucket0[bucket0maxPointer].remove(item);
            bucket0Size--;

            // Push to fixedNodes list
            item.isFixed = true;
            return item;
        }
        else
        {
            // If bucket is empty, decrease maxPointer and run again
            bucket0maxPointer--;
            return popFromBucketKey(partition);
        } 
    }
    else 
    {
        if (bucket1[bucket1maxPointer].size() > 0)
        {
            Node item = bucket1[bucket1maxPointer].front();
            bucket1[bucket1maxPointer].remove(item);
            bucket1Size--;

            item.isFixed = true;
            return item;
        }
        else
        {
            bucket1maxPointer--;
            return popFromBucketKey(partition);
        } 
    }
}

int Bucket::gainSum() const
{
    int gain = 0;

    for (auto const &k: bucket0)
    {
        gain = gain + (k.first * k.second.size());
    }

    for (auto const &k: bucket1)
    {
        gain = gain + (k.first * k.second.size());
    }

    return gain;
}
