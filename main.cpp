#include "main.h"

#include <iostream>

// Method Definitions
/**
 * This method collects the command line arguments and places them in the
 * respective variable.
 *
 * @param arguments popsize, numChars, sdaStates, seed, runs, maxGens, maxMuts,
 * seqNum, tournSize, crossoverOp, crossoverRate, mutationRate, cullingRate,
 * randomCulling
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
  arg = arguments[11];  // seqNum
  seqNum = stoi(arg, &pos);
  arg = arguments[12];  // tournSize
  tournSize = stoi(arg, &pos);
  cout << "Arguments Captured!" << endl;
  return 0;
}

int initAlg(const string &pathToSeqs) {
  srand48(time(nullptr));  // use system time as random number seed
  srand48(seed);           // read the random number seed
  // vector<vector<int>> sequences = getSequences(pathToSeqs);
  goalSeq = getSequences(pathToSeqs)[seqNum];
  seqLen = (int)goalSeq.size();
  fits.reserve(popsize);

  pop = new SDA[popsize];
  for (int idx = 0; idx < popsize; ++idx) {
    pop[idx] = SDA(sdaStates, numChars, 2, seqLen);
  }

  testSeq.reserve(seqLen);
  charSeq.reserve(seqLen);
  for (int idx = 0; idx < seqLen; ++idx) {
    testSeq.push_back(-1);
    charSeq.push_back('a');
  }

  populationBestIdx = -1;
  populationBestFit = (BIGGER_BETTER ? 0.0 : MAXFLOAT);
  return 0;
}

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

int cmdLineIntro(ostream &outp) {
  outp << "Amino Acid/DNA Sequence Matcher." << endl;
  outp << "Fitness Function: Naive." << endl;
  outp << "Matching Sequence (length " << seqLen << "): ";
  intToChar(goalSeq, charSeq);
  printVector<ostream, char>(outp, charSeq, "", "", true);
  outp << "See read.me for system and experiment parameters." << endl;
  outp << endl;
  return 0;
}

int makeReadMe(ostream &outp) {
  cmdLineIntro(outp);
  outp << "Experiment Parameters" << endl;
  outp << "Population Size: " << popsize << endl;
  outp << "Number of States: " << sdaStates << endl;
  outp << "Alphabet Size: " << numChars << endl;
  outp << "Max Generations: " << maxGens << endl;
  outp << "Report Every: " << REPORT_EVERY << " generations" << endl;
  outp << "Tournament Size: " << tournSize << endl;
  return 0;
}

int initPop(int run) {
  if (runs > 1) {
    cout << "Beginning Run " << run << " of " << runs << endl;
  }
  fits.clear();
  for (int idx = 0; idx < popsize; ++idx) {
    pop[idx].randomize();
    fits.push_back(fitness(pop[idx]));
  }
  cout << "Population Generated!" << endl;
  return 0;
}

double fitness(SDA &sda) {
  double val = 0.0;
  sda.fillOutput(testSeq);
  for (int i = 0; i < seqLen; ++i) {
    if (testSeq[i] == goalSeq[i]) {
      val += 1;
    }
  }
  return val;
}

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

double report(ofstream &outp, int run, int rptNum, bool biggerBetter) {
  vector<double> stats =
      calcStats<double>(fits, biggerBetter);  // {mean, stdDev, 95CI, best}
  multiStream printAndSave(cout, outp);

  printAndSave << left << setw(5) << run;
  printAndSave << left << setw(4) << rptNum;
  printAndSave << left << setw(10) << stats[0];
  printAndSave << left << setw(12) << stats[2];
  printAndSave << left << setw(10) << stats[1];
  printAndSave << left << setw(8) << stats[3];
  printAndSave << left << setw(8) << (stats[3] / seqLen) * 100 << "%";
  printAndSave << "\n";
  return stats[3];
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

int matingEvent(bool biggerBetter) {
  int numMuts;
  SDA child1, child2;

  vector<int> idxs = tournSelect(tournSize, biggerBetter);

  child1.copy(pop[idxs[0]]);
  child2.copy(pop[idxs[1]]);
  if (drand48() < crossoverRate) {
    if (crossoverOp == 0)
      child1.twoPointCrossover(child2);
    else if (crossoverOp == 1)
      child1.oneStateCrossover(child2);
  }

  if (drand48() < mutationRate) {
    if (dynamicMutOperator == 0) {
      child1.mutate(curNumTransMuts, curNumRespMuts);
      child2.mutate(curNumTransMuts, curNumRespMuts);
    }
  }

  pop[idxs.end()[-1]] = child1;
  pop[idxs.end()[-2]] = child2;

  fits[idxs.end()[-1]] = fitness(child1);
  fits[idxs.end()[-2]] = fitness(child2);
  return 0;
}

vector<int> tournSelect(int size, bool decreasing) {
  vector<int> tournIdxs;
  int idxToAdd;

  tournIdxs.reserve(size);
  if (size == popsize) {
    for (int idx = 0; idx < size; idx++) {
      tournIdxs.push_back(idx);
    }
  } else {
    do {
      idxToAdd = (int)lrand48() % popsize;
      if (count(tournIdxs.begin(), tournIdxs.end(), idxToAdd) == 0) {
        tournIdxs.push_back(idxToAdd);
      }
    } while (tournIdxs.size() < size);
  }

  sort(tournIdxs.begin(), tournIdxs.end(), compareFitness);
  if (decreasing) {
    reverse(tournIdxs.begin(), tournIdxs.end());
  }
  return tournIdxs;
}

bool compareFitness(int popIdx1, int popIdx2) {
  return fits[popIdx1] < fits[popIdx2];
}

int culling(double percentage, bool rndPick, bool biggerBetter) {
  if (percentage == 1) {  // Cull all but the best SDA
    pop[0].copy(pop[populationBestIdx]);
    fits[0] = fitness(pop[0]);
    for (int idx = 1; idx < popsize; ++idx) {
      pop[idx].randomize();
      fits[idx] = fitness(pop[idx]);
    }
    return 0;
  }

  // Otherwise, cull percentage% of the population
  int numKillings = (int)(popsize * percentage);
  vector<int> contestants;
  if (rndPick) {  // Cull random percentage% of the population
    contestants.reserve(numKillings);
    contestants = tournSelect(numKillings, !biggerBetter);
  } else {  // Cull worst percentage% of the population
    contestants.reserve(popsize);
    contestants = tournSelect(popsize, !biggerBetter);
  }

  int idxToCull;
  for (int cnt = 0; cnt < numKillings; cnt++) {
    idxToCull = contestants[cnt];
    pop[idxToCull].randomize();
    fits[idxToCull] = fitness(pop[idxToCull]);
  }
  return 0;
}

int runReport(ostream &outp, bool biggerBetter) {
  auto maxIterator = minmax_element(fits.begin(), fits.end());
  int bestIdx;
  if (biggerBetter) {
    bestIdx = (int)distance(fits.begin(), maxIterator.second);
  } else {
    bestIdx = (int)distance(fits.begin(), maxIterator.first);
  }

  multiStream printAndSave(cout, outp);
  printAndSave << "The best fitness is " << fits[bestIdx] << "\n";
  pop[bestIdx].fillOutput(testSeq);
  printAndSave << left << setw(20) << "Best Match: ";
  intToChar(testSeq, charSeq);
  printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
  printAndSave << left << setw(20) << "Matches: ";
  printMatches<multiStream, int>(printAndSave, testSeq, goalSeq, true);
  printAndSave << left << setw(20) << "Desired Sequence: ";
  intToChar(goalSeq, charSeq);
  printVector<multiStream, char>(printAndSave, charSeq, "", "", true);

  outp << "SDA" << endl;
  pop[bestIdx].printSDA(outp);

  vector<double> fitsCopy = fits;
  sort(fitsCopy.begin(), fitsCopy.end());
  if (biggerBetter) {
    reverse(fitsCopy.begin(), fitsCopy.end());
  }
  printVector<multiStream, double>(printAndSave, fitsCopy,
                                   "Fitness Values: ", " ", true);
  return bestIdx;
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

int expReport(ostream &outp, vector<double> bestFits, SDA bestSDA,
              bool biggerBetter) {
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
  intToChar(testSeq, charSeq);
  printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
  printAndSave << left << setw(20) << "Matches: ";
  printMatches<multiStream, int>(printAndSave, testSeq, goalSeq, true);
  printAndSave << left << setw(20) << "Desired Sequence: ";
  intToChar(goalSeq, charSeq);
  printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
  printAndSave << "\nSDA:\n";
  bestSDA.printSDA(cout);
  bestSDA.printSDA(outp);
  return 0;
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

/**
 * This method...
 *
 * Input variables:
 * 1.   Population size
 * 2.   Number of characters (in SDA)
 * 3.   Number of states
 * 4.   Seed
 * 5.   Number of runs
 * 6.   Maximum number of mating events
 * 7.   Default Number of Transition Mutations
 * 8.   Default Number of Response Mutations
 * 9.   Dynamic Mutation Operator? (0 -> static, >0 -> dynamic implementation to
 * use)
 * 10.  Upper Bound on Mutations
 * 11.  Sequence number
 * 12.  Tournament size
 * 13.  Crossover operator
 * 14.  Crossover Rate
 * 15.  Mutation Rate
 * 16.  Culling Rate
 * 17.  Random Culling
 * 18.  Culling Every _ Reporting Intervals
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
  getArgs(argv);
  string pathToSeqs = "./Sequences.dat";
  string filename;
  ofstream runStats, expStats, readMe;

  vector<double> bests;
  bests.reserve(runs);
  double expBestFit = (BIGGER_BETTER ? 0 : MAXFLOAT);

  initAlg(pathToSeqs);
  SDA expBestSDA = SDA(sdaStates, numChars, 2, seqLen);
  cmdLineIntro(cout);
  char dynamicMessage[20];
  sprintf(dynamicMessage, "%s%d",
          (dynamicMutOperator == 0 ? "Static" : "Dynamic"), dynamicMutOperator);
  sprintf(pathToOut,
          "./AAMOut/AAMatch on Seq%d with %.1fmilMMEs, %04dPS, %02dSt, "
          " %dTS/",
          seqNum, (double)maxGens / 1000000, popsize, sdaStates, tournSize);
  mkdir(pathToOut, 0777);
  expStats.open(string(pathToOut) + "./exp.dat", ios::out);
  readMe.open(string(pathToOut) + "./read.me", ios::out);
  makeReadMe(readMe);
  readMe.close();

  int tmp;
  for (int run = 1; run < runs + 1; ++run) {
    curNumTransMuts = initNumTransMuts;
    curNumRespMuts = initNumRespMuts;
    initPop(run);
    char runNumStr[20];
    sprintf(runNumStr, "%02d", run);
    filename = string(pathToOut) + "run" + string(runNumStr) + ".dat";
    runStats.open(filename, ios::out);
    printExpStatsHeader(cout);
    printExpStatsHeader(runStats);
    report(runStats, run, 0, BIGGER_BETTER);

    int gen = 1;
    int stallCount = 0;
    double best = (BIGGER_BETTER ? 0 : MAXFLOAT);
    while (gen <= maxGens &&
           (stallCount < TERM_CRIT || gen < MIN_GEN_RATIO * maxGens)) {
      matingEvent(BIGGER_BETTER);

      if (gen % REPORT_EVERY == 0) {
        tmp = report(runStats, run, (int)gen / (REPORT_EVERY), BIGGER_BETTER);
        if ((BIGGER_BETTER && tmp > best) || (!BIGGER_BETTER && tmp < best)) {
          best = tmp;
          stallCount = 0;
        } else {
          stallCount++;
        }
      }

      if (gen % (int)(CULLING_EVERY * REPORT_EVERY) == 0 &&
          stallCount < TERM_CRIT) {
        culling(cullingRate, randomCulling, BIGGER_BETTER);
      }
      gen++;
    }

    tmp = runReport(expStats, BIGGER_BETTER);
    if ((BIGGER_BETTER && fits[tmp] > expBestFit) ||
        (!BIGGER_BETTER && fits[tmp] < expBestFit)) {
      expBestFit = fits[tmp];
      expBestSDA.copy(pop[tmp]);
    }
    bests.push_back(fits[tmp]);
    runStats.close();
  }

  ofstream best;
  best.open(string(pathToOut) + "./best.dat", ios::out);
  expReport(best, bests, expBestSDA, BIGGER_BETTER);
  best.close();
  delete[] pop;
  cout << "Program Completed Successfully!" << endl;
  return 0;
}
