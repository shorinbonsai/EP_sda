#ifndef EP_H
#define EP_H

#include <bits/stdc++.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "SDA.h"

using namespace std;

struct Individual {
  SDA sda;
  int hammingFitness;
  int boutWins;
};
vector<Individual> calculateRelativeFiteness(vector<Individual> &population,
                                             const vector<int> &target,
                                             mt19937 &gen);

class Ep {
 public:
  explicit Ep(int SDANumStates, int SDAOutputLen, vector<int> &sequence,
              int numGens, ostream &MyFile, int numChars, int popSize,
              int boutSize, int seed, bool roulette);
  ~Ep();

 private:
  int Evolve(vector<Individual> currentPop, const vector<int> &target,
             int numGens, ostream &MyFile, bool roulette);
  double fitness(SDA &sda, vector<int> &sequence);
  int printPopFits(ostream &outStrm, vector<double> &popFits);
  int numChars;
  int popSize;
  int boutSize;
  int responseLength = 2;
  int seed;
  bool roulette;
  mt19937 rng;

  vector<double> initFits;
  vector<double> currentFits;
  vector<double> newGenFits;
  vector<Individual> currentPop;
  vector<Individual> newPop;

  // Helper to compute seed (time-based if seed=0)
  static unsigned compute_seed(int seed_param) {
    if (seed_param == 0) {
      // Use system clock for time-based seed
      return static_cast<unsigned>(
          std::chrono::system_clock::now().time_since_epoch().count());
    } else {
      return static_cast<unsigned>(seed_param);
    }
  }
};

#endif