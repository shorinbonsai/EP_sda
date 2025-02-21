#include "Ep.h"
#include <algorithm>
#include <random>
// #include <unordered_set>

const int BOUTS_PER_INDIVIDUAL = 10;



vector<Individual> calculateRelativeFiteness(vector<Individual> &population, const vector<int> &target, unsigned seed = 0)
{
    vector<Individual> individuals;
    individuals.reserve(population.size());
    for (Individual &indiv : population)
    {
        indiv.hammingFitness = hammingFitness(indiv.sda, target);
        indiv.boutWins = 0;
        individuals.push_back(indiv);
    }

    // Setup random number generator with provided seed or time-based seed
    if (seed == 0)
    {
        seed = static_cast<unsigned>(time(nullptr));
    }
    mt19937 gen(seed);
    uniform_int_distribution<> distrib(0, population.size() - 1);

    for (Individual &individual : individuals)
    {
        for (int i = 0; i < BOUTS_PER_INDIVIDUAL; i++)
        {
            int opponentIdx = distrib(gen);
            // Make sure opponent is not the same individual
            while (opponentIdx == &individual - &individuals[0])
            {
                opponentIdx = distrib(gen);
            }

            if (individual.hammingFitness < individuals[opponentIdx].hammingFitness)
            {
                individual.boutWins++;
            }
        }
    }

    // // Sort by fitness
    // sort(individuals.begin(), individuals.end(), [](const Individual &a, const Individual &b)
    //      { return a.hammingFitness < b.hammingFitness; });

    // // Calculate relative fitness
    // double totalFitness = 0.0;
    // for (const Individual &individual : individuals)
    // {
    //     totalFitness += 1.0 / (1.0 + individual.hammingFitness);
    // }

    // for (Individual &individual : individuals)
    // {
    //     individual.relativeFitness = (1.0 / (1.0 + individual.hammingFitness)) / totalFitness;
    // }

    return individuals;
}

Ep::Ep(int SDANumStates, int SDAOutputLen, vector<int> &sequence, int numGens, ostream &MyFile, int numChars, int popSize, int boutSize)
{
    this->numChars = numChars;
    this->popSize = popSize;
    this->boutSize = boutSize;
    // SDA *currentPop, *newPop;
    Individual *currentPop;
    currentPop = new Individual[popSize];
    
    initFits.reserve(popSize);
    // init population
    for (int i = 0; i < popSize; ++i)
    {
        SDA newSDA = SDA(SDANumStates, numChars, responseLength, SDAOutputLen);
        Individual newInd;
        newInd.hammingFitness =  hammingFitness(newSDA, sequence);
        newInd.boutWins = 0;
        currentPop[i] = newInd;
        newPop[i] = currentPop[i];
        initFits.push_back(currentPop[i].hammingFitness);
        
    }
    // Evolve(SDANumStates, SDAOutputLen, numGens, MyFile);
}

/**
 * Destructor.
 */
Ep::~Ep() = default;

double fitness(SDA &sda, vector<int> &sequence)
{
    return 0.0;
}

/**
 * Calculates fitness based on Hamming distance between SDA output and target sequence.
 * Hamming distance counts positions where characters differ.
 * Returns raw Hamming distance (lower is better).
 *
 * @param sda The SDA to evaluate
 * @param target The target DNA sequence (using 0=A, 1=C, 2=G, 3=T)
 * @return Raw Hamming distance (number of differences)
 */
double hammingFitness(SDA &sda, const vector<int> &sequence)
{
    // Get SDA output
    vector<int> output = sda.rtnOutput(false, cout);

    if (output[0] == -1)
    {
        return INT_MAX; // SDA failed to generate output, return worst possible fitness
    }

    double differences = 0.0;
    int sequenceLength = min(output.size(), sequence.size());

    // Count positions that don't match
    for (int i = 0; i < sequenceLength; i++)
    {
        if (output[i] != sequence[i])
        {
            differences++;
        }
    }

    // Add penalties for length differences (each extra/missing character counts as a mismatch)
    differences += abs((int)(output.size() - sequence.size()));

    return differences; // Lower is better
}

int Ep::Evolve(vector<SDA> &population, const vector<int> &target, int numGens, ostream &MyFile, unsigned seed = 0)
{
    
    newPop = new Individual[popSize * 2];

    MyFile << "Initial Fitness: " << endl;
    printPopFits(MyFile, initFits);

    // Evolution
    for (int i = 0; i < numGens; ++i)
    {
    }

    return 0;
}

int Ep::printPopFits(ostream &outStrm, vector<double> &popFits)
{
    outStrm << "Fitness Values: ";
    int count = 0;
    bool first = true;
    for (double fit : popFits)
    {
        // This ensures commas after each fitness value other than the last
        if (!first)
        {
            outStrm << ", ";
        }
        outStrm << fit;
        if (fit > 150)
        {
            count++;
        }
        first = false;
    }
    outStrm << "\n";
    // outStrm << "Above 0.5: " << count << "\n";
    return 0;
}

// int Generational::genEvolver(int SDANumStates, int SDAOutputLen, int numGenerations, Topology T, ostream& MyFile) {
//     SDA *currentPop, *newPop, cp;
//     currentPop = new SDA[genPopSize];
//     newPop = new SDA[genPopSize];
//     genPopFits.reserve(genPopSize);
//     int modVal = numGenerations / 10;

//     // Step 1: initialize the population
//     for (int i = 0; i < genPopSize; ++i) {
//         currentPop[i] = SDA(SDANumStates, genSDANumChars, genSDAResponseLength, SDAOutputLen);;// place member into population
//         genPopFits.push_back(genCalcFitness(currentPop[i], T));// calculate new members fitness
//     }

//     MyFile << "Initial Pop Fitness Value: " << endl;
//     genPrintPopFits(MyFile, genPopFits); // print population fitness

//     // Step 2: Evolution
//     for (int gen = 0; gen < numGenerations; ++gen) {
//         genNewPopFits.clear();
//         genNewPopFits.reserve(genPopSize);

//         // Perform elitisim
//         vector<double> elite(elitism, DBL_MAX);// vector holding two most elite fitnesses
//         int idx1, idx2;// index locations of the two most elite SDA
//         for (int e = 0; e < genPopFits.size(); e++){// go through the genertions pop fitness
//             if(elite[0] > genPopFits[e]){
//                 elite[0] = genPopFits[e];// update elite 1 fitness value
//                 idx1 = e; // memner fitness better than elite 1
//             }
//             else if(elite[1] > genPopFits[e]){
//                 elite[1] = genPopFits[e];// update elite 2 fintess value
//                 idx2 = e;// member fitness better than elite 2 and not used by elite 1
//             }
//         }
//         // copy the elite SDA to the new population
//         newPop[0].copy(currentPop[idx1]);
//         newPop[1].copy(currentPop[idx2]);

//         // Store the fitness of the most elite members in the first two positions
//         genNewPopFits.push_back(genPopFits[idx1]);
//         genNewPopFits.push_back(genPopFits[idx2]);

//         // Generate the new population
//         genMatingEvent(currentPop, newPop, T);

//         // Replace current population with new population
//         for (int mem = 0; mem < genPopSize; mem++){
//             currentPop[mem] = newPop[mem];
//         }
//         if(gen % modVal == 0) genPrintPopFits(MyFile, genNewPopFits);// print every 10th generation
//     }
//     MyFile << "Final Fitness of Run: ";
//     genPrintPopFits(MyFile, genNewPopFits);
//     return 0;
// }