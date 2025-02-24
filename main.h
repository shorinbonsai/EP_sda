#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

#include "Ep.h"
#include "SDA.h"

using namespace std;

#define BIGGER_BETTER (bool)false

// Experiment Parameters
int popsize;
int numChars;
int sdaStates;
double *sdaProbs;
int seed;
int runs;
int maxGens;
int tournSize;
int seqNum;
double mutationRate;
int populationBestIdx;
double populationBestFit;  // current best population fitness
double prevBestFitness = 0;
int RICounter;  // Report Interval counter
int seqLen;

char pathToOut[150];


vector<int> goalSeq;
vector<int> testSeq;
vector<char> charSeq;

vector<vector<int>> getSequences(const string &pathToSeqs);
vector<int> seqToVector(const string &seq);
int intToChar(const vector<int> &from, vector<char> &to);
int getArgs(char *arguments[]);

vector<vector<int>> getSequences(const string &pathToSeqs) {
  string tmp;
  ifstream in(pathToSeqs);
  vector<vector<int>> rtn;
  for (int seqIdx = 0; seqIdx < 6; ++seqIdx) {
    if (in.is_open()) {
      getline(in, tmp);
      getline(in, tmp);
      getline(in, tmp);
      rtn.push_back(seqToVector(tmp));
      getline(in, tmp);
    }
  }
  in.close();
  return rtn;
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

int intToChar(const vector<int> &from, vector<char> &to) {
  for (int idx = 0; idx < seqLen; ++idx) {
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

/**
 * This method collects the command line arguments and places them in the
 * respective variable.
 *
 * @param arguments popsize, numChars, sdaStates, seed, runs, maxGens, seqNum,
 * tournSize
 *
 * @return
 */
int getArgs(char *arguments[]) {
  size_t pos;
  string arg;
  arg = arguments[1];  // popsize
  popsize = stoi(arg, &pos);
  arg = arguments[2];  // numChars
  numChars = stoi(arg, &pos);
  arg = arguments[3];  // sdaStates
  sdaStates = stoi(arg, &pos);
  arg = arguments[4];  // seed
  seed = stoi(arg, &pos);
  arg = arguments[5];  // runs
  runs = stoi(arg, &pos);
  arg = arguments[6];  // maxGens
  maxGens = stoi(arg, &pos);

  arg = arguments[7];  // seqNum
  seqNum = stoi(arg, &pos);
  arg = arguments[8];  // tournSize
  tournSize = stoi(arg, &pos);

  cout << "Arguments Captured!" << endl;
  return 0;
}

template<class T>
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

    double mean = sum / (double) vals.size();
    double stdDevSum = 0.0;
    for (int val: vals) {
        stdDevSum += pow((double) val - mean, 2);
    }
    double stdDev = sqrt(stdDevSum / ((double) vals.size() - 1.0));
    double CI95 = 1.96 * (stdDev / sqrt(vals.size()));

    return {mean, stdDev, CI95, bestVal}; // {mean, stdDev, 95CI, best}
}

// int initAlg(const string &pathToSeqs) {
//     srand48(time(nullptr)); // use system time as random number seed
//     srand48(seed);           // read the random number seed
//     //vector<vector<int>> sequences = getSequences(pathToSeqs);
//     goalSeq = getSequences(pathToSeqs)[seqNum];
//     seqLen = (int) goalSeq.size();
//     fits.reserve(popsize);

//     pop = new SDA[popsize];
//     for (int idx = 0; idx < popsize; ++idx) {
//         pop[idx] = SDA(sdaStates, numChars, 2, seqLen);
//     }

//     testSeq.reserve(seqLen);
//     charSeq.reserve(seqLen);
//     for (int idx = 0; idx < seqLen; ++idx) {
//         testSeq.push_back(-1);
//         charSeq.push_back('a');
//     }

//     populationBestIdx = -1;
//     populationBestFit = (BIGGER_BETTER ? 0.0 : MAXFLOAT);
//     return 0;
// }