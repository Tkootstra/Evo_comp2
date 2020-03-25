
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
    int edgeCounter = 0;
    while (std::getline(dataFile, line))
    {
        
        // parse all info from .txt file line by line
        std::istringstream iss(line);
        std::vector<std::string> singleSplittedLine((std::istream_iterator<std::string>(iss)),
            std::istream_iterator<std::string>());

        int vertexNumber = std::stoi(singleSplittedLine[0]);
        int nConnections = std::stoi(singleSplittedLine[2]);
        std::list<int> otherLocs; 
        edgeCounter += nConnections;

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
void printVector(std::vector<int> sol)
{
    for (size_t i = 0; i <10; i++)
    {
        std::cout << sol[i]  << " "; 
    }
    std::cout << endl;
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
    // printVector(solution);
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

int hammingDistance(std::vector<int> parent1, std::vector<int> parent2)
{
    int distance = 0;
    for (size_t i = 0; i < parent1.size(); i++)
    {
        if (parent1[i] != parent2[i]) distance++;
    }
    return distance;
}

std::vector<int> invertSolution(std::vector<int> solution)
{
    for (size_t i = 0; i < solution.size(); i++)
    {
        if (solution[i] == 0) solution[i] = 1;
        else solution[i] = 0;
    }
    return solution;
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
    
    // std::cout << "doing hamming" << endl;
    // std:: cout << "pa1 size: " << parent1.size() << " pa2 size: " << parent2.size() << endl;
    // check the hamming distance between the parents. If it is larger than l/2, invert parent1.
    int distance = hammingDistance(parent1, parent2);
    if (distance > (parent2.size() /2)) parent1 = invertSolution(parent1);

    // std::cout << "computing balance" << endl;

    std::vector<int> child(parent1.size());
    int zeroCounter = 0;
    int oneCounter = 0;
    int similarBits = 0;
    // compute balance of missing ones
    for (size_t i = 0; i < child.size(); i++)
    {
        if (parent2[i] == parent1[i])
        {
            similarBits++;
            if (parent2[i] == 0) zeroCounter++;
            if (parent2[i] == 1) oneCounter++;
        }

    }
    // std::cout << "amount of zeroes: " << zeroCounter << endl;
    // std::cout << "amount of ones: " << oneCounter << endl;
    // std::cout << "amount of disagreebles: " << parent1.size() - similarBits << endl;

    int amountOfOnesToSample = (parent1.size() /2) - oneCounter;
    int amountOfZeroesToSample = (parent1.size() /2) - zeroCounter;
    int notEqualBits = parent2.size() - similarBits /2;
    int balance = zeroCounter - oneCounter ;
    std::vector<int> zeroes(amountOfZeroesToSample, 0);
    std::vector<int> ones(amountOfOnesToSample, 1);
    zeroes.insert(zeroes.end(), ones.begin(), ones.end());
    std::vector<int> samplePool = zeroes;
    std::random_shuffle(samplePool.begin(), samplePool.end());

    // std::cout << "making child" << endl;

    int sampleIdx = 0;
    for (size_t i = 0; i < parent1.size(); i++)
    {   // constraints: 1. same bit value if indice values are equal. 
        //              2. if bit not equal, fill with random 1 or 0, but these must be equal in number for the TOTAL solution (samplen zonder teruglegging)
        int bit;
        if (parent2[i] == parent1[i]) bit = parent1[i];
        else 
        {   
            bit = samplePool[sampleIdx];
            sampleIdx ++;
        }
        child[i] = bit;
    }
    // std::cout << "sampleidx" << sampleIdx << endl;
    // std::cout << "pool size: "<< samplePool.size() << endl;
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
        
        // dit 5-10 keer draaien, of totdat de returned solution niet meer beter is dan de vorige
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

// std::vector<vector<double> > MLS_test(const std::vector<Node> nodeList, int iterations, float perturbationRatio)
// { // Run ILS. Mutates a solution and carries on if better, reverts to previous if not.

//     std::cout << "Running ILS..." << endl;

//     std::vector<double> cR(2);
//     std::vector<vector<double> > combinedResults(iterations);

//     // Keep track of results
//     std::pair<int, std::vector<int> > firstResultPair;
//     std::pair<int, std::vector<int> > secondResultPair;
//     int firstResult;
//     int secondResult;
    
//     std::chrono::steady_clock::time_point begin;
//     std::chrono::duration<double> dur;

//     /* RUN ONE LOCAL SEARCH */
//     begin = std::chrono::steady_clock::now();

//     // Create a random solution
//     std::vector<int> solution = makeRandomSolution(500);
//     std::vector<int> tempSolution(500);

//     Graph graaf;
//     graaf.initializeGraph(nodeList, solution);
//     firstResultPair = singleFMrun(graaf);
//     firstResult = firstResultPair.first;

//     // Calculate elapsed time
//     dur = (std::chrono::steady_clock::now() - begin);

//     // Save results
//     cR[0] = firstResult;
//     cR[1] = dur.count();
//     combinedResults[0] = cR;
    
//     std::cout << "Iteration " << 0 + 1 << ": Score " << firstResult << " | Time: " << cR[1] << "s." << endl;

//     // Perturb inital solution
//     tempSolution = firstResultPair.second;

//     for (size_t i = 1; i < iterations; i++)
//     {
//         // int sum = 0;
//         // for (auto& n:tempSolution) sum += n;
//         // std::cout << sum << endl;

//         begin = std::chrono::steady_clock::now();

//         // Run second local search
//         Graph g;
//         g.initializeGraph(nodeList, tempSolution);

//         secondResultPair = singleFMrun(g);
//         secondResult = secondResultPair.first;

//         // If this result is better (lower score), the new solution becomes 
//         // our optimal result from this pass
//         // int counter = 0;
//         // std::cout << secondResult - firstResult << endl;
//         if (secondResult < firstResult)
//         {
//             std::cout << "hoi" << endl;
//             firstResult = secondResult;
//             solution = secondResultPair.second;
//             tempSolution = solution;
//             // counter ++;
//         }
//         // std::cout << "number of FM passes: " << counter << endl;
//         else    
//         {
//             // std::cout << "making new random sol" << endl;
//             g.initializeGraph(nodeList, makeRandomSolution(500));

//             firstResultPair = singleFMrun(g);
//             firstResult = firstResultPair.first;

//             tempSolution = firstResultPair.second;
//         }
        
//         // Calculate elapsed time
//         dur = (std::chrono::steady_clock::now() - begin);

//         // Save results
//         cR[0] = firstResult;
//         cR[1] = dur.count();
//         combinedResults[i] = cR;
        
//         std::cout << "Iteration " << i + 1 << ": Score " << firstResult << " | Time: " << cR[1] << "s." << endl;
//     }

//     return combinedResults;
// }


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

int randomNumberGenerator(int max)
{
    //TODO: hier zit een grote bug in, werkt totaal niet zoals het moet.
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<int> dist(0, max);
    // std::cout << "sampling random number" << endl;
    int randomNumber = dist(rng);
    return randomNumber;

    // srand((unsigned) time(0));
    // int randomNumber;
    // randomNumber = (rand() % max) +1;
    // return randomNumber;

    
}

std::vector<int> sampleFromList(std::list<vector<int>> sols, int idx)
{
    int index = 0;
    std::vector<int> match;
    for (auto sol: sols)
    {
        if (index == idx) match = sol;
        index++;

    }
    return match;
}
std::pair<vector<int>, vector<int> > sampleSolutions(std::list<vector<int> > solutions)
{
    // std::cout << "solution size: " << solutions.size() << endl;

    int first = randomNumberGenerator(solutions.size()-1);
    int second = randomNumberGenerator(solutions.size()-1);
    // std::cout << "succesful sampling :" << first << "  " << second << endl;
    while (first == second)
    {   
        // std::cout << "resampling" << endl;
        second = randomNumberGenerator(solutions.size());
        
    }

    // std::cout << "pointer shizzle" << endl;
    // auto firstSolutionPointer = solutions.begin();
    // std::advance(firstSolutionPointer, first);
    // std::vector<int> firstSol = *firstSolutionPointer;


    // auto secondSolutionPointer = solutions.begin();
    // std::advance(secondSolutionPointer, second);
    // std::vector<int> secondSol = *secondSolutionPointer;
    // std::cout << "pointer shizzle passed" << endl;

    std::vector<int> firstSol = sampleFromList(solutions, first);
    std::vector<int> secondSol = sampleFromList(solutions, second);

    // std::cout << "first sample size: " << firstSol.size() << " second sample size: " << secondSol.size() << endl;
    std::pair<vector<int>, vector<int>> returnValue = std::make_pair(firstSol, secondSol);
    return returnValue;
}



std::vector<vector<double> > geneticLocalSearch(const std::vector<Node> nodeList, int iterations, int populationSize)
{
    //TODO: implement GLS:
    // 1. population size 50
    // 2. no generations: each iteration, select 2 random parents, use uniform crossover to generate one child.
    // 3. FM local search on the child
    // 4. Let this child compete with the worst solution (dus alle solutions checken)
    // 5. 


    // TODO: bug => volgens mij gaat het met de indexen niet goed omdat je 50x normaal FM doet en daarna met crossover nog een keer los;
    // Segmentation fault (core dumped)
 
    std::cout << "Running Genetic Local Search..." << endl;

    // make n random solutions

    std::list<vector<int> > startingSolutions = makeMultipleRandomSolutions(500, populationSize);
    
    std::vector<double> cR(2);
    std::vector<vector<double> > combinedResults(iterations);

    std::chrono::steady_clock::time_point begin;
    std::chrono::duration<double> dur;
   
    
    bool cont = true;
    int iter = 0;
    while (cont)
    {
        
        begin = std::chrono::steady_clock::now();
        // std::cout << "doing new gen" << endl;
        // search for worst solution
        std::vector<int> worstSol;
        int worstScore = 0;
        std::list<vector<int>> allSolutions;

        for (auto& solution: startingSolutions)
        {   
            if (iter >= iterations)
            {   cont = false;
                break;
            }
            // std::cout << "initializing graph" << endl;
            Graph g;
            g.initializeGraph(nodeList, solution);
            

            // std::cout << "doing single FM" << endl;
            std::pair<int, vector<int> > resultPair = singleFMrun(g);
            int score = resultPair.first;
            std::vector<int> currentSol = resultPair.second;
            if (score > worstScore) 
            {
                worstScore = score;
                worstSol = currentSol;
            }
            allSolutions.push_back(currentSol);
           
            dur = (std::chrono::steady_clock::now() - begin);
            cR[0] = resultPair.first;
            cR[1] = dur.count();
            // std::cout << "size of combinedresults:" << combinedResults.size() << endl;
            combinedResults[iter] = cR;
            // std::cout << iter << ' ' << score << endl;
            iter++;
        }
        
        if (iter >= iterations)
        {
             cont = false;
             break;
        }
        
        // sample 2 random parents
        // std::cout << "sampling solutions" << endl;
        // std::cout << "allsolutions size: " << allSolutions.size() << endl;
        std::pair<vector<int>, vector<int>> randomSample = sampleSolutions(allSolutions);
        std::vector<int> first  = randomSample.first;
        std::vector<int> second = randomSample.second;
        // std::cout << "doing uniform xover" << endl;
        std::vector<int> generatedChild = uniformCrossOver(first, second);

        int sum = 0;
        for (auto& n:generatedChild) sum += n;
        // std::cout << "sum of child: "  << sum << endl;


        // run single FM with child
        
        Graph child;
        child.initializeGraph(nodeList, generatedChild);
        // std::cout << "doing single FM" << endl;
        std::pair<int, vector<int> > childRes = singleFMrun(child);
        // std::cout << "comparing sols" << endl;
        
        int childScore = childRes.first;
        std::vector<int> childSolution = childRes.second;
        
        // compare child solution with worst one. if child better than worst, delete worst from list and insert child
        if (childScore <= worstScore) 
        {
            allSolutions.remove(worstSol);
            allSolutions.push_back(childSolution);
            
        }
        // std::cout << "vector size:" << combinedResults.size() <<" index: " << iter << endl;
        startingSolutions = allSolutions;
        dur = (std::chrono::steady_clock::now() - begin);
        cR[0] = childRes.first;
        cR[1] = dur.count();
        combinedResults[iter] = cR;
        std::cout << "Iteration " << iter << ": Score " << childRes.first << " | Time: " << cR[1] << "s." << endl;
        iter++;
    }

    return combinedResults;

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
    srand(time(NULL));
    int runs = 1000;

    // parse nodes from txt file
    std::vector<Node> nodeList = parseGraph();
    
    // Run MLS
    // std::vector<vector<double> > resultsMLS = multiStartLocalSearch(nodeList, runs);
    // writeToFile(resultsMLS, "MLS.txt");

    // std::vector<vector<double> > kek = MLS_test(nodeList, runs, 0.1);
    // // Run ILS
    std::vector<vector<double> > resultsILS = iterativeLocalSearch(nodeList, runs, 0.1);
    // writeToFile(resultsILS, "ILS.txt");

    std::vector<vector<double> > resultsGLS = geneticLocalSearch(nodeList,runs, 50);
    // writeToFile(resultsGLS, "GLS.txt");

} 
