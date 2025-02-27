#include "main.h"

#include <filesystem>
namespace filesystem = std::filesystem;
#include <iostream>

int populationBestIdx;
double populationBestFit;
int expReport(ostream &outp, vector<double> bestFits, SDA bestSDA,
              bool biggerBetter, const Ep &ep);

int runReport(ostream &outp, bool biggerBetter, const Ep &ep) {
  auto maxIterator =
      minmax_element(ep.currentFits.begin(), ep.currentFits.end());
  int bestIdx;
  if (biggerBetter) {
    bestIdx = (int)distance(ep.currentFits.begin(), maxIterator.second);
  } else {
    bestIdx = (int)distance(ep.currentFits.begin(), maxIterator.first);
  }

  multiStream printAndSave(cout, outp);
  printAndSave << "The best fitness is " << ep.currentFits[bestIdx] << "\n";
  SDA sdaCopy = ep.currentPop[bestIdx].sda;
  sdaCopy.fillOutput(testSeq);
  printAndSave << left << setw(20) << "Best Match: ";
  intToChar(testSeq, charSeq, ep);
  printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
  printAndSave << left << setw(20) << "Matches: ";
  printMatches<multiStream, int>(printAndSave, testSeq, goalSeq, true);
  printAndSave << left << setw(20) << "Desired Sequence: ";
  intToChar(goalSeq, charSeq, ep);
  printVector<multiStream, char>(printAndSave, charSeq, "", "", true);

  outp << "SDA" << endl;
  sdaCopy.print(outp);

  vector<double> fitsCopy = ep.currentFits;
  sort(fitsCopy.begin(), fitsCopy.end());
  if (biggerBetter) {
    reverse(fitsCopy.begin(), fitsCopy.end());
  }
  printVector<multiStream, double>(printAndSave, fitsCopy,
                                   "Fitness Values: ", " ", true);
  return bestIdx;
}

int main(int argc, char *argv[]) {
  getArgs(argv);
  string pathToSeqs = "./Sequences.dat";
  string filename;
  ofstream runStats, expStats, readMe;

  vector<double> bests;
  bests.reserve(runs);
  double expBestFit = (BIGGER_BETTER ? 0 : MAXFLOAT);

  goalSeq = getSequences(pathToSeqs)[seqNum];
  seqLen = (int)goalSeq.size();
  SDA expBestSDA = SDA(sdaStates, numChars, 2, seqLen);

  sprintf(pathToOut, "./SEMOut/SEMatch on Seq%d with %dGens, %04dPS, %02dSt/",
          seqNum, maxGens, popsize, sdaStates);
  filesystem::create_directories(pathToOut);
  // mkdir(pathToOut, 0777);
  expStats.open(string(pathToOut) + "./exp.dat", ios::out);
  ofstream outFile(filename);
  int tmp;
  Ep newEp(sdaStates, seqLen, goalSeq, maxGens, outFile, numChars, popsize,
           tournSize, seed, roulette);
  for (int run = 1; run < runs; ++run) {
    char runNumStr[20];
    sprintf(runNumStr, "%02d", run);
    filename = string(pathToOut) + "run" + string(runNumStr) + ".dat";
    Ep newEp(sdaStates, seqLen, goalSeq, maxGens, outFile, numChars, popsize,
             tournSize, seed, roulette);

    runStats.open(filename, ios::out);
    printExpStatsHeader(cout);
    printExpStatsHeader(runStats);

    newEp.Evolve(newEp.currentPop, newEp.sequence, maxGens, runStats, roulette,
                 run);

    tmp = runReport(expStats, BIGGER_BETTER, newEp);
    if ((BIGGER_BETTER && newEp.currentFits[tmp] > expBestFit) ||
        (!BIGGER_BETTER && newEp.currentFits[tmp] < expBestFit)) {
      expBestFit = newEp.currentFits[tmp];
      expBestSDA = newEp.currentPop[tmp].sda;
    }
    bests.push_back(expBestFit);
    runStats.close();
    seed++;
  }

  ofstream best;
  best.open(string(pathToOut) + "./best.dat", ios::out);
  expReport(best, bests, expBestSDA, BIGGER_BETTER, newEp);
  best.close();
  expStats.close();
  return 0;
}

int expReport(ostream &outp, vector<double> bestFits, SDA bestSDA,
              bool biggerBetter, const Ep &ep) {
  vector<double> stats = calcStats<double>(bestFits, BIGGER_BETTER);
  multiStream printAndSave(cout, outp);
  printAndSave << "Experiment Report:" << "\n";
  printAndSave << "Best Run: " << populationBestIdx + 1 << "\n";
  printAndSave << "Best Fitness: " << bestFits[populationBestIdx] << " of "
               << seqLen << "\n";
  printAndSave << "Fitness 95% CI: " << stats[0] << " +- " << stats[2] << "\n";
  printAndSave << "\n";
  printAndSave << left << setw(20) << "Best Match: ";
  bestSDA.fillOutput(testSeq);
  intToChar(testSeq, charSeq, ep);
  printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
  printAndSave << left << setw(20) << "Matches: ";
  printMatches<multiStream, int>(printAndSave, testSeq, goalSeq, true);
  printAndSave << left << setw(20) << "Desired Sequence: ";
  intToChar(goalSeq, charSeq, ep);
  printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
  printAndSave << "\nSDA:\n";
  bestSDA.print(cout);
  bestSDA.print(outp);
  return 0;
}