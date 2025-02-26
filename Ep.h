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

int printExpStatsHeader(ostream &outp) {
  outp << left << setw(5) << "Run";
  outp << left << setw(4) << "RI";
  outp << left << setw(10) << "Mean";
  outp << left << setw(12) << "95% CI";
  outp << left << setw(10) << "SD";
  outp << left << setw(8) << "Best";
  outp << left << setw(10) << "% Correct";
  outp << endl;
  return 0;
}

template <class T>
vector<double> calcStats(vector<T> vals, bool biggerBetter) {
  double sum = 0.0;
  double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);

  int val;
  for (int idx = 0; idx < vals.size(); ++idx) {
    val = vals[idx];
    sum += val;
    if ((biggerBetter && val > bestVal) || (!biggerBetter && val < bestVal)) {
      bestVal = val;
      populationBestIdx = idx;
      populationBestFit = bestVal;
    }
  }

  double mean = sum / (double)vals.size();
  double stdDevSum = 0.0;
  for (int val : vals) {
    stdDevSum += pow((double)val - mean, 2);
  }
  double stdDev = sqrt(stdDevSum / ((double)vals.size() - 1.0));
  double CI95 = 1.96 * (stdDev / sqrt(vals.size()));

  return {mean, stdDev, CI95, bestVal};  // {mean, stdDev, 95CI, best}
}

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

double report(ostream &outp, int run, int rptNum, bool biggerBetter,
              vector<double> fits, const Ep &ep) {
  vector<double> stats =
      calcStats<double>(fits, biggerBetter);  // {mean, stdDev, 95CI, best}
  multiStream printAndSave(cout, outp);

  printAndSave << left << setw(5) << run;
  printAndSave << left << setw(4) << rptNum;
  printAndSave << left << setw(10) << stats[0];
  printAndSave << left << setw(12) << stats[2];
  printAndSave << left << setw(10) << stats[1];
  printAndSave << left << setw(8) << stats[3];
  printAndSave << left << setw(8) << (stats[3] / ep.seqLen) * 100 << "%";
  printAndSave << "\n";
  return stats[3];
}

vector<int> seqToVector(const string &seq) {
  vector<int> sequence;
  for (char c : seq) {
    if (c == 'g' || c == 'G') {
      sequence.push_back(0);
    } else if (c == 'c' || c == 'C') {
      sequence.push_back(1);
    } else if (c == 'a' || c == 'A') {
      sequence.push_back(2);
    } else if (c == 't' || c == 'T') {
      sequence.push_back(3);
    }
  }
  return sequence;
}

int intToChar(const vector<int> &from, vector<char> &to, const Ep &ep) {
  for (int idx = 0; idx < ep.seqLen; ++idx) {
    switch (from[idx]) {
      case 0:
        to[idx] = 'G';
        break;
      case 1:
        to[idx] = 'C';
        break;
      case 2:
        to[idx] = 'A';
        break;
      case 3:
        to[idx] = 'T';
        break;
    }
  }
  return 0;
}

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