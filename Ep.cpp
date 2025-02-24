#include "Ep.h"

#include <algorithm>
#include <random>
// #include <unordered_set>

const int BOUTS_PER_INDIVIDUAL = 10;

/**
 * Calculates fitness based on Hamming distance between SDA output and target
 * sequence. Hamming distance counts positions where characters differ. Returns
 * raw Hamming distance (lower is better).
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
        return 9999999999.0; // SDA failed to generate output, return worst
                             // possible fitness
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

    // Add penalties for length differences (each extra/missing character counts
    // as a mismatch)
    differences += abs((int)(output.size() - sequence.size()));

    return differences; // Lower is better
}

vector<Individual> calculateRelativeFiteness(vector<Individual> &population,
                                             const vector<int> &target,
                                             mt19937 &gen)
{
    vector<Individual> individuals;
    individuals.reserve(population.size());
    for (Individual &indiv : population)
    {
        indiv.hammingFitness = hammingFitness(indiv.sda, target);
        indiv.boutWins = 0;
        individuals.push_back(indiv);
    }

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

    return individuals;
}

Ep::Ep(int SDANumStates, int SDAOutputLen, vector<int> &sequence, int numGens,
       ostream &MyFile, int numChars, int popSize, int boutSize, int seed,
       bool roulette)
    : rng(compute_seed(seed))
{
    this->numChars = numChars;
    this->popSize = popSize;
    this->boutSize = boutSize;
    this->roulette = roulette;
    currentPop.reserve(popSize);
    newPop.reserve(popSize * 2);

    initFits.reserve(popSize);
    currentFits.reserve(popSize);
    newGenFits.reserve(popSize * 2);
    // init population
    for (int i = 0; i < popSize; ++i)
    {
        SDA newSDA(SDANumStates, numChars, responseLength, SDAOutputLen);
        Individual newInd;
        newInd.hammingFitness = hammingFitness(newSDA, sequence);
        newInd.boutWins = 0;
        currentPop.push_back(newInd);
        initFits.push_back(newInd.hammingFitness);
    }
    newPop = currentPop;
    currentFits = initFits;
    Evolve(currentPop, sequence, numGens, MyFile, roulette);
}

/**
 * Destructor.
 */
Ep::~Ep() = default;

// double fitness(SDA &sda, vector<int> &sequence)
// {
//     return 0.0;
// }

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

/**
 * Selects 50% of individuals using roulette wheel selection based on relative
 * fitness. Ensures no duplicates in the selection.
 *
 * @param individuals Vector of evaluated individuals
 * @param seed Random seed for reproducibility (0 uses time-based seed)
 * @return Selected individuals (50% of original population)
 */
vector<Individual> selectByRouletteWheel(vector<Individual> &individuals,
                                         mt19937 &gen)
{
    int selectionSize = individuals.size() / 2;
    vector<Individual> selected;
    selected.reserve(selectionSize);

    uniform_real_distribution<> distrib(0.0, 1.0);

    // Track indices of selected individuals to avoid duplicates
    unordered_set<int> selectedIndices;

    while (selected.size() < selectionSize)
    {
        double randomValue = distrib(gen);
        double cumulativeProbability = 0.0;

        for (int i = 0; i < individuals.size(); i++)
        {
            cumulativeProbability += individuals[i].boutWins;

            if (randomValue <= cumulativeProbability &&
                selectedIndices.find(i) == selectedIndices.end())
            {
                selected.push_back(individuals[i]);
                selectedIndices.insert(i);
                break;
            }
        }

        // If we couldn't select a new individual, adjust the distribution
        if (selected.size() < selectionSize &&
            selectedIndices.size() >= individuals.size())
        {
            break; // Safety check in case we can't find enough unique individuals
        }
    }
    return selected;
}

/**
 * Selects the top 50% of individuals based on tournament wins.
 *
 * @param individuals Vector of evaluated individuals
 * @return Selected individuals (50% of original population)
 */
vector<Individual> selectByRank(vector<Individual> &individuals)
{
    // Select top 50%
    int selectionSize = individuals.size() / 2;
    vector<Individual> selected(individuals.begin(),
                                individuals.begin() + selectionSize);

    return selected;
}

int Ep::Evolve(vector<Individual> currentPop, const vector<int> &target,
               int numGens, ostream &MyFile, bool roulette)
{
    const int RUNS = 30;
    MyFile << "Initial Fitness: " << endl;
    printPopFits(MyFile, initFits);

    newGenFits.clear();
    newGenFits.reserve(popSize * 2);

    // Evolution
    for (int i = 0; i < numGens; ++i)
    {
        newPop.clear();
        newPop = currentPop;
        for (int j = 0; j < popSize; j++)
        {
            Individual ind = currentPop[j];
            Individual newInd = ind;
            newInd.sda.mutate(1);
            newPop.push_back(newInd);
        }
        newPop = calculateRelativeFiteness(newPop, target, rng);
        // Sort by tournament wins (descending order)
        sort(newPop.begin(), newPop.end(),
             [](const Individual &a, const Individual &b)
             {
                 return a.boutWins > b.boutWins;
             });
        if (roulette)
        {
            currentPop = selectByRouletteWheel(newPop, rng);
        }
        else
        {
            currentPop = selectByRank(newPop);
        }
        currentFits.clear();
        for (Individual &ind : currentPop)
        {
            currentFits.push_back(ind.hammingFitness);
        }
        if (i % 10 == 0)
        {
            MyFile << "Generation " << i << " Fitness: " << endl;
            printPopFits(MyFile, currentFits);
        }
    }

    return 0;
}
