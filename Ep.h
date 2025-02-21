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
    explicit Ep(int SDANumStates, int SDAOutputLen, vector<int> &sequence, int numGens, ostream &MyFile, int numChars, int popSize, int boutSize);
    ~Ep();
private:
    int Evolve(vector<SDA> &population, const vector<int> &target, int numGens, ostream &MyFile, unsigned seed = 0);
    double fitness(SDA &sda, vector<int> &sequence);
    int printPopFits(ostream &outStrm, vector<double> &popFits);
    int numChars;
    int popSize;
    int boutSize;
    int responseLength = 2;

    vector<double> initFits;
    vector<double> runningFits;
    vector<double> newGenFits;
    Individual *newPop;
};

#endif