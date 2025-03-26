#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>
namespace filesystem = std::filesystem;

#include "SDA.h"

using namespace std;

// System Parameters
#define REPORT_EVERY (int)10
#define TERM_CRIT 50
int CULLING_EVERY;
#define BIGGER_BETTER (bool)true
double MIN_GEN_RATIO = 0.5;
#define ROULETTE (bool)false
const int CULL_CAP = 80;

// Experiment Parameters
int popsize;
int numChars;
int sdaStates;
double *sdaProbs;
int seed;
int runs;
int maxGens;
int initNumTransMuts;
int initNumRespMuts;
int curNumTransMuts;
int curNumRespMuts;
int upBoundMuts;
int numMuts;
int dynamicMutOperator;
int tournSize;
int seqNum;
int crossoverOp;
double crossoverRate;
double mutationRate;
double cullingRate;
bool randomCulling;
int populationBestIdx;
double populationBestFit;  // current best population fitness
double prevBestFitness = 0;
int RICounter;  // Report Interval counter

vector<int> goalSeq;
vector<int> testSeq;
vector<char> charSeq;
int seqLen;

// SDA *pop;
vector<SDA> pop;
vector<SDA> doublePop;
vector<double> fits;
vector<double> doubleFits;
vector<double> relativeFits;
vector<double> doubleRelativeFits;

char pathToOut[150];

// Program Method Declarations:
int getArgs(char *arguments[]);
int initAlg(const string &pathToSeqs);
vector<vector<int>> getSequences(const string &pathToSeqs);
int cmdLineIntro(ostream &outp);
int makeReadMe(ostream &outp);
int initPop(int run);
double fitness(SDA &sda);
int printExpStatsHeader(ostream &outp);
double report(ofstream &outp, int run, int rptNum, bool biggerBetter);
template <class T>
vector<double> calcStats(vector<T> vals, bool biggerBetter);
int matingEvent(bool biggerBetter, double cullingRate);
vector<int> tournSelect(int size, bool decreasing);
bool compareFitness(int popIdx1, int popIdx2);
int culling(double percentage, bool biggerBetter);
int runReport(ostream &outp, bool biggerBetter);
template <class T1, class T2>
int printMatches(T1 &outp, const vector<T2> &test, const vector<T2> &goal,
                 bool newline);
int expReport(ostream &outp, vector<double> bestFits, SDA bestSDA,
              bool biggerBetter);
int calcRelativeFitness();
int selectByRank();
int selectByRoulette();
// Helper Method Declarations:
vector<int> seqToVector(const string &seq);
int intToChar(const vector<int> &from, vector<char> &to);
template <class T1, class T2>
int printVector(T1 &outp, vector<T2> vec, const string &msg, const string &sep,
                bool newline);
template <class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs,
                      const string &msg, const string &sep, bool newline);

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
