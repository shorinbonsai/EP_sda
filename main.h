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
#define roulette (bool)false

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
// int populationBestIdx;
// double populationBestFit;  // current best population fitness
double prevBestFitness = 0;
int RICounter;  // Report Interval counter
// extern int seqLen;

char pathToOut[150];

vector<int> goalSeq;
vector<int> testSeq;
vector<char> charSeq;

vector<vector<int>> getSequences(const string &pathToSeqs);
vector<int> seqToVector(const string &seq);
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
