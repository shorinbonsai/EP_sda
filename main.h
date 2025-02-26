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

// int printExpStatsHeader(ostream &outp) {
//     outp << left << setw(5) << "Run";
//     outp << left << setw(4) << "RI";
//     outp << left << setw(10) << "Mean";
//     outp << left << setw(12) << "95% CI";
//     outp << left << setw(10) << "SD";
//     outp << left << setw(8) << "Best";
//     outp << left << setw(10) << "% Correct";
//     outp << endl;
//     return 0;
// }

// template <class T>
// vector<double> calcStats(vector<T> vals, bool biggerBetter)
// {
//     double sum = 0.0;
//     double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);

//     int val;
//     for (int idx = 0; idx < vals.size(); ++idx)
//     {
//         val = vals[idx];
//         sum += val;
//         if ((biggerBetter && val > bestVal) || (!biggerBetter && val <
//         bestVal))
//         {
//             bestVal = val;
//             populationBestIdx = idx;
//             populationBestFit = bestVal;
//         }
//     }

//     double mean = sum / (double)vals.size();
//     double stdDevSum = 0.0;
//     for (int val : vals)
//     {
//         stdDevSum += pow((double)val - mean, 2);
//     }
//     double stdDev = sqrt(stdDevSum / ((double)vals.size() - 1.0));
//     double CI95 = 1.96 * (stdDev / sqrt(vals.size()));

//     return {mean, stdDev, CI95, bestVal}; // {mean, stdDev, 95CI, best}
// }

// // Helper Class
// class multiStream : public ostream
// {
// public:
//     multiStream(ostream &os1, ostream &os2) : os1(os1), os2(os2) {}

//     template <class T>
//     multiStream &operator<<(const T &x)
//     {
//         os1 << x;
//         os2 << x;
//         return *this;
//     }

// private:
//     ostream &os1;
//     ostream &os2;
// };

// double report(ofstream &outp, int run, int rptNum, bool biggerBetter,
// vector<double> fits)
// {
//     vector<double> stats = calcStats<double>(fits, biggerBetter); // {mean,
//     stdDev, 95CI, best} multiStream printAndSave(cout, outp);

//     printAndSave << left << setw(5) << run;
//     printAndSave << left << setw(4) << rptNum;
//     printAndSave << left << setw(10) << stats[0];
//     printAndSave << left << setw(12) << stats[2];
//     printAndSave << left << setw(10) << stats[1];
//     printAndSave << left << setw(8) << stats[3];
//     printAndSave << left << setw(8) << (stats[3] / seqLen) * 100 << "%";
//     printAndSave << "\n";
//     return stats[3];
// }

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