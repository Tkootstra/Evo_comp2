
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

// c++ -std=c++17 -stdlib=libc++ main.cpp -o prog 
// ./prog

using namespace std;

std::vector<Node> parseGraph()
{// Parses the Graph500.txt file and returns a vector of Nodes
    std::ifstream dataFile("Graph500.txt");
    std::string line;
    std::vector<Node> nodeList;
    // std::vector<std::pair<int, int> > netList;
    // int neti = 0;

    while (std::getline(dataFile, line))
    {
        
        // parse all info from .txt file line by line
        std::istringstream iss(line);
        std::vector<std::string> singleSplittedLine((std::istream_iterator<std::string>(iss)),
            std::istream_iterator<std::string>());

        int vertexNumber = std::stoi(singleSplittedLine[0]);
        int nConnections = std::stoi(singleSplittedLine[2]);
        std::list<int> otherLocs; 

        // parse all other indice locations
        for (size_t i = 3; i < singleSplittedLine.size(); i++)
        {   
            int loc = std::stoi(singleSplittedLine[i]);
            otherLocs.push_back(loc - 1); 
            // netList.push_back(std::make_pair(vertexNumber, loc));   
        }
        // make new nodes for each input line
        Node newNode;
        newNode.initializeNode(vertexNumber - 1, nConnections, otherLocs);
        nodeList.push_back(newNode);
    }

    return nodeList;
}

std::vector<int> makeRandomSolution(int stringLength)
{// Creates random solution with equal 0s and 1s
    bool equal = false;
    std::vector<int> solution(stringLength);

    // srand(time(NULL));

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
{ // Calls makeRandomSolution 'amount' times
    srand(time(NULL));
    std::list<vector<int> > allSolutions;
    for (size_t i = 0; i < amount; i++)
    {
        
        allSolutions.push_back(makeRandomSolution(stringLength));
    }
    return allSolutions;
}

std::vector<int> perturbSolution(std::vector<int> solution, const float ratio)
{ // Perturb/mutate solution by a certain ratio
    std::vector<int> tempSolution = solution;
    int amountToPerturb = (int)(ratio * solution.size());
    
    for (size_t i = 0; i < amountToPerturb; i++)
    {
        int mutationLoc1 = rand() % solution.size();
        int toMutate1 = tempSolution[mutationLoc1];

        int mutationLoc2 = rand() % solution.size();
        int toMutate2 = tempSolution[mutationLoc2];

        // If bit-values are the same, try again (must be balanced)
        while (toMutate2 == toMutate1)
        {
            mutationLoc2 = rand() % solution.size();
            toMutate2 = tempSolution[mutationLoc2];
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

std::vector<int> getBestSolution(const Graph &g)
{ // Retrieve the vector representation of the partitioning of graph g
    std::vector<int> solution(g.Nodes.size());

    for (size_t i = 0; i < g.Nodes.size(); i++)
    {
        solution[i] = g.Nodes[i].belongsToWhichPartition;
    }

    return solution;
}

std::vector<int> uniformCrossOver(std::vector<int> parent1, std::vector<int> parent2)
{
    // TODO: aanpassen zoals in de opdracht
    std::vector<int> child(parent1.size());
    for (size_t i = 0; i < parent1.size(); i++)
    {
        std::vector<int> options(2);
        options[0] = parent1[i];
        options[1] = parent2[i];
        int choiceIndex = rand() % 2;
        child[i] = options[choiceIndex];
    }
    return child;
}

Bucket initGain(Graph &graph, Bucket &currentBucket)
{ // Compute gain for all free nodes
    Node current;

    int gain;
    
    for (size_t i = 0; i < graph.Nodes.size(); i++)
    {
        // Check if node not in fixedNodes, otherwise skip
        if (!graph.Nodes[i].isFixed)
        {
            // Create a temporary graph to manipulate
            current = graph.Nodes[i];
            
            // Flip partition and retrieve how it affects the gain
            graph.Nodes[current.indexLocation].flipPartition();

            // Calculate gain and adjust bucket accordingly
            gain = graph.countSingleCellConnections(current.indexLocation, 0);
            graph.Nodes[current.indexLocation].flipPartition();

            currentBucket.addToBucket(current.belongsToWhichPartition, gain, current);
            graph.Nodes[i].gain = gain;
        }    
    }

    return currentBucket;    
}

Bucket updateGain(Graph &graph, Bucket &currentBucket, const std::list<int> neighbors)
{ // Do the same as computeGain but only for certain nodes and run 'updateBucket'
    // instead of 'addToBucket'

    Node current;
    
    int gain;

    // Only loop through necessary connections
    for (auto const &i: neighbors)
    {
        // Check if current node i is not in fixedNodes
        if (!graph.Nodes[i].isFixed)
        {
            // Create a temporary graph to manipulate
            current = graph.Nodes[i];
            
            // Flip partition and retrieve how it affects the gains
            graph.Nodes[current.indexLocation].flipPartition();
            gain = graph.countSingleCellConnections(current.indexLocation, current.gain);
            graph.Nodes[current.indexLocation].flipPartition();

            // Update bucket to reflect gain
            currentBucket.updateBucket(current.belongsToWhichPartition, gain, current);
            graph.Nodes[i].gain = gain;
        } 
    }
    return currentBucket;
}

std::pair<int, std::vector<int> > singleFMrun(Graph g)
{ // Run a single FM pass (250 moves back and forth). Returns the best solution. 
    // Some pre-work for initializing the Bucket
    std::map<int, std::list<Node> > bucket0;
    std::map<int, std::list<Node> > bucket1;
    
    // Create bucket and get size of partition 0 (1 is equally large)
    Bucket results = Bucket(bucket0, bucket1);
    int b0size = results.bucket0Size;
    
    // Compute the initial score and gains
    results = initGain(g, results);
    int score = results.gainSum();
    
    // Keep track of the best score and accompanying graph
    int bestGainSum = score;
    Graph bestScoreGraph = g;

    // Do until both buckets are empty:
        // Move node with highest gain from 0 to 1
            // Remove node from bucket (the node becomes LOCKED)
        // Do same from 1 to 0, otherwise solution isn't valid

        // Recompute the gain for all neighboring nodes of the nodes that have been moved

    while (b0size > 0) 
    {
        // Move node from 0 to 1
        Node nodeToChange0 = results.popFromBucketKey(0);
        g.Nodes[nodeToChange0.indexLocation].flipPartition();

        // Move node from 1 to 0
        Node nodeToChange1 = results.popFromBucketKey(1);
        g.Nodes[nodeToChange1.indexLocation].flipPartition();

        // We now have a new valid partition; update gains for neighbors
        results = updateGain(g, results, nodeToChange0.ConnectionLocations);
        results = updateGain(g, results, nodeToChange1.ConnectionLocations);

        // Count the score again
        score = results.gainSum();
        // std::cout << score << endl;
        
        // Update metadata
        b0size = results.bucket0Size;

        // Keep track of best score
        if (score > bestGainSum)
        {
            bestGainSum = score;
            bestScoreGraph = g;
        }
    }

    // Get the solution that matches the graph with the best score
    int bestScore = bestScoreGraph.countConnections();
    std::vector<int> bestScoreSolution = getBestSolution(bestScoreGraph);
    std::pair<int, std::vector<int> > returnPair = std::make_pair(bestScore, bestScoreSolution);

    // int sum = 0;
    // for (auto& n:bestScoreSolution) sum += n;
    // std::cout << sum << endl;
    
    return returnPair;
}

std::vector<vector<double> > multiStartLocalSearch(const std::vector<Node> nodeList, int iterations)
{ // Run a simple MLS. Just makes a random solution every time.

    std::cout << "Running MLS..." << endl;

    std::vector<double> cR(2);
    std::vector<vector<double> > combinedResults(iterations);
    
    std::pair<int, std::vector<int> > r;
    std::chrono::steady_clock::time_point begin;
    std::chrono::duration<double> dur;
    Graph graaf;

    for (size_t i = 0; i < iterations; i++)
    {
        // Start timer and init graph
        begin = std::chrono::steady_clock::now();
        graaf.initializeGraph(nodeList, makeRandomSolution(500));

        // Run one local search
        r = singleFMrun(graaf);

        // Calculate elapsed time
        dur = (std::chrono::steady_clock::now() - begin);

        // Save results
        cR[0] = r.first;
        cR[1] = dur.count();
        combinedResults[i] = cR;
        
        std::cout << "Iteration " << i + 1 << ": Score " << r.first << " | Time: " << cR[1] << "s." << endl;
    }

    return combinedResults;
}

std::vector<vector<double> > iterativeLocalSearch(const std::vector<Node> nodeList, int iterations, float perturbationRatio)
{ // Run ILS. Mutates a solution and carries on if better, reverts to previous if not.

    std::cout << "Running ILS..." << endl;

    std::vector<double> cR(2);
    std::vector<vector<double> > combinedResults(iterations);

    // Keep track of results
    std::pair<int, std::vector<int> > firstResultPair;
    std::pair<int, std::vector<int> > secondResultPair;
    int firstResult;
    int secondResult;
    
    std::chrono::steady_clock::time_point begin;
    std::chrono::duration<double> dur;

    // Create a random solution
    std::vector<int> solution = makeRandomSolution(500);
    std::vector<int> tempSolution(500);

    /* RUN ONE LOCAL SEARCH */
    begin = std::chrono::steady_clock::now();

    Graph graaf;
    graaf.initializeGraph(nodeList, solution);
    firstResultPair = singleFMrun(graaf);
    firstResult = firstResultPair.first;

    // Calculate elapsed time
    dur = (std::chrono::steady_clock::now() - begin);

    // Save results
    cR[0] = firstResult;
    cR[1] = dur.count();
    combinedResults[0] = cR;
    
    std::cout << "Iteration " << 0 + 1 << ": Score " << firstResult << " | Time: " << cR[1] << "s." << endl;

    // Perturb inital solution
    tempSolution = perturbSolution(solution, perturbationRatio);

    for (size_t i = 1; i < iterations; i++)
    {
        // int sum = 0;
        // for (auto& n:tempSolution) sum += n;
        // std::cout << sum << endl;

        begin = std::chrono::steady_clock::now();

        // Run second local search
        Graph g;
        g.initializeGraph(nodeList, tempSolution);

        secondResultPair = singleFMrun(g);
        secondResult = secondResultPair.first;

        // If this result is better (lower score), the new solution becomes 
        // our optimal result from this pass
        if (secondResult < firstResult)
        {
            firstResult = secondResult;
            solution = secondResultPair.second;
            tempSolution = solution;
        }
        else
        {
            tempSolution = perturbSolution(solution, perturbationRatio);
        }
        
        // Calculate elapsed time
        dur = (std::chrono::steady_clock::now() - begin);

        // Save results
        cR[0] = firstResult;
        cR[1] = dur.count();
        combinedResults[i] = cR;
        
        std::cout << "Iteration " << i + 1 << ": Score " << firstResult << " | Time: " << cR[1] << "s." << endl;
    }

    return combinedResults;
}

std::pair<vector<int>, vector<int> > sampleSolutions(std::list<vector<int> > solutions)
{
    
    int first = rand() % solutions.size();
    int second = rand() % solutions.size();
    while (first == second)
    {
        int second = rand() % solutions.size();
    }


    auto firstSolutionPointer = solutions.begin();
    std::advance(firstSolutionPointer, first);
    std::vector<int> firstSol = *firstSolutionPointer;


    auto secondSolutionPointer = solutions.begin();
    std::advance(secondSolutionPointer, second);
    std::vector<int> secondSol = *secondSolutionPointer;
   
    std::pair<vector<int>, vector<int>> returnValue = std::make_pair(firstSol, secondSol);
    

    return returnValue;
}

void geneticLocalSearch(const std::vector<Node> nodeList, int iterations, int populationSize)
{
    //TODO: implement GLS:
    // 1. population size 50
    // 2. no generations: each iteration, select 2 random parents, use uniform crossover to generate one child.
    // 3. FM local search on the child
    // 4. Let this child compete with the worst solution (dus alle solutions checken)
    // 5. 

 
    std::cout << "Running Genetic Local Search..." << endl;

    // make n random solutions

    std::list<vector<int> > startingSolutions = makeMultipleRandomSolutions(500, populationSize);
    
    std::vector<double> cR(2);
    std::vector<vector<double> > combinedResults(iterations);

    std::chrono::steady_clock::time_point begin;
    std::chrono::duration<double> dur;
    
    int iter = 0;
    while (iter <= iterations)
    {
        
        std::cout << "doing new gen" << endl;
        // search for worst solution
        std::vector<int> worstSol;
        int worstScore = 0;
        std::list<vector<int>> allSolutions;

        for (auto& solution: startingSolutions)
        {   
            std::cout << "initializing graph" << endl;
            Graph g;
            g.initializeGraph(nodeList, solution);
            std::cout << "doing single FM" << endl;
            std::pair<int, vector<int> > resultPair = singleFMrun(g);
            int score = resultPair.first;
            std::vector<int> currentSol = resultPair.second;
            if (score > worstScore) 
            {
                worstScore = score;
                worstSol = currentSol;
            }
            allSolutions.push_back(currentSol);
            iter++;
            std::cout << iter << ' ' << score << endl;
        }
        
        // sample 2 random parents
        std::cout << "Doing crossover" << endl;
        std::pair<vector<int>, vector<int>> randomSample = sampleSolutions(allSolutions);
        std::vector<int> first  = randomSample.first;
        std::vector<int> second = randomSample.second;
        std::cout << "doing uniform xover" << endl;
        std::vector<int> generatedChild = uniformCrossOver(first, second);

        int sum = 0;
        for (auto& n:generatedChild) sum += n;
        std::cout << "sum of child: "  << sum << endl;


        // run single FM with child
        std::cout << "doing single FM" << endl;
        Graph child;
        child.initializeGraph(nodeList, generatedChild);
        std::pair<int, vector<int> > childRes = singleFMrun(child);
        std::cout << "comparing sols" << endl;
        iter++;
        int childScore = childRes.first;
        std::vector<int> childSolution = childRes.second;
        
        // compare chil solution with worst one. if child better than worst, delete worst from list and insert child
        if (childScore <= worstScore) 
        {
            allSolutions.remove(worstSol);
            allSolutions.push_back(childSolution);
            
        }

        startingSolutions = allSolutions;

        std::cout << "childscore:" << childScore << endl;

    }

}


void writeToFile(const std::vector<vector<double> > results, const std::string fileName)
{ // Write results to .txt
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
    
    int runs = 1000;

    // parse nodes from txt file
    std::vector<Node> nodeList = parseGraph();
    
    // // Run MLS
    // std::vector<vector<double> > resultsMLS = multiStartLocalSearch(nodeList, runs);
    // writeToFile(resultsMLS, "MLS.txt");

    // // Run ILS
    // std::vector<vector<double> > resultsILS = iterativeLocalSearch(nodeList, runs, 0.1);
    // writeToFile(resultsILS, "ILS.txt");

    geneticLocalSearch(nodeList,runs, 50);

} 
