#ifndef EP_H
#define EP_H

#include <bits/stdc++.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "SDA.h"
#define BIGGER_BETTER (bool)false

using namespace std;

extern int populationBestIdx;
extern double populationBestFit;

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
  int Evolve(vector<Individual> currentPop, const vector<int> &target,
             int numGens, ostream &MyFile, bool roulette, int run);
  vector<int> sequence;
  int seqLen;
  vector<double> initFits;
  vector<double> currentFits;
  vector<double> newGenFits;
  vector<Individual> currentPop;
  vector<Individual> newPop;

 private:
  double fitness(SDA &sda, vector<int> &sequence);
  int printPopFits(ostream &outStrm, vector<double> &popFits);
  int numChars;
  int popSize;
  int boutSize;
  int responseLength = 2;
  int seed;
  bool roulette;

  mt19937 rng;

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

int printExpStatsHeader(ostream &outp);
template <class T>
vector<double> calcStats(vector<T> vals, bool biggerBetter);
double report(ostream &outp, int run, int rptNum, bool biggerBetter,
              vector<double> fits, const Ep &ep);
vector<int> seqToVector(const string &seq);
int intToChar(const vector<int> &from, vector<char> &to, const Ep &ep);

// Helper Class
class multiStream : public ostream {
 public:
  multiStream(ostream &os1, ostream &os2) : os1(os1), os2(os2) {}

  template <class T>
  multiStream &operator<<(const T &x) {
    os1 << x;
    os2 << x;
    return *this;
  }

 private:
  ostream &os1;
  ostream &os2;
};

template <class T1, class T2>
int printVector(T1 &outp, vector<T2> vec, const string &msg, const string &sep,
                bool newline) {
  outp << msg;
  for (auto val : vec) {
    outp << val << sep;
  }
  if (newline) outp << "\n";
  return 0;
}

template <class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs,
                      const string &msg, const string &sep, bool newline) {
  outp << msg;
  for (auto idx : idxs) {
    outp << vec[idx] << sep;
  }
  if (newline) outp << "\n";
  return 0;
}

template <class T1, class T2>
int printMatches(T1 &outp, const vector<T2> &test, const vector<T2> &goal,
                 bool newline) {
  int count = 0;
  for (int idx = 0; idx < goal.size(); idx++) {
    if (test[idx] == goal[idx]) {
      outp << "X";
      count++;
    } else {
      outp << " ";
    }
  }
  if (newline) outp << "\n";
  return count;
}

#endif