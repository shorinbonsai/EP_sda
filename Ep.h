#ifndef EP_H
#define EP_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h>
#include "SDA.h"

struct Individual
{
    SDA sda;
    int hammingFitness;
    int boutWins;


    vector<Individual> calculateRelativeFiteness(vector<Individual> &population, const vector<int> &target, unsigned seed = 0);
};

class Ep
{

public:
    explicit Ep(int SDANumStates, int SDAOutputLen, vector<int> &sequence, int numGens, ostream &MyFile, int numChars, int popSize, int boutSize, int seed);
    ~Ep();
private:
    int Evolve(vector<Individual> currentPop, const vector<int> &target, int numGens, ostream &MyFile, unsigned seed = 0);
    double fitness(SDA &sda, vector<int> &sequence);
    int printPopFits(ostream &outStrm, vector<double> &popFits);
    int numChars;
    int popSize;
    int boutSize;
    int responseLength = 2;
    int seed;

    vector<double> initFits;
    vector<double> runningFits;
    vector<double> newGenFits;
    vector<Individual> currentPop;
    vector<Individual> newPop;
    // Individual *currentPop;
    // Individual *newPop;
};

#endif