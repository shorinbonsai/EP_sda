#include "main.h"

#include <iostream>

// Method Definitions
/**
 * This method collects the command line arguments and places them in the
 * respective variable.
 *
 * @param arguments popsize, numChars, sdaStates, seed, runs, maxGens,
 * seqNum, tournSize, numMuts, cullingRate, CULLING_EVERY
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
  arg = arguments[9];  // numMuts
  numMuts = stoi(arg, &pos);
  arg = arguments[10];  // cullingRate
  cullingRate = stod(arg, &pos);
  arg = arguments[11];
  CULLING_EVERY = stoi(arg, &pos);
  // arg = arguments[12];
  // runarg = stoi(arg, &pos);
  cout << "Arguments Captured!" << endl;
  return 0;
}

int initAlg(const string &pathToSeqs) {
  srand48(seed);  // read the random number seed

  goalSeq = getSequences(pathToSeqs)[seqNum];
  seqLen = (int)goalSeq.size();
  fits.resize(popsize);
  doubleFits.resize(popsize * 2);
  relativeFits.resize(popsize);
  doubleRelativeFits.resize(popsize * 2);

  pop.resize(popsize);
  doublePop.reserve(popsize * 2);
  for (auto &sda : pop) {
    sda = SDA(sdaStates, numChars, 2, seqLen, 1.5 * sdaStates);
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
  outp << "Culling Every: " << CULLING_EVERY * REPORT_EVERY << " generations"
       << endl;
  outp << "Culling Winners: Worst " << (int)(cullingRate * 100)
       << "% of the population" << endl;
  if (numMuts == 0) {
    outp << "Number of Mutations: Random (1-3)" << endl;
  } else {
    outp << "Number of Mutations: " << numMuts << endl;
  }

  return 0;
}

int initPop(int run) {
  if (runs > 1) {
    cout << "Beginning Run " << run << " of " << runs << endl;
  }
  fits.clear();
  populationBestIdx = -1;
  populationBestFit = (BIGGER_BETTER ? 0.0 : MAXFLOAT);
  for (int idx = 0; idx < popsize; ++idx) {
    pop[idx].randomize();
    fits.push_back(fitness(pop[idx]));
  }
  for (int i = 0; i < popsize * 2; ++i) {
    noveltyFits.push_back(-1);
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

int calcNoveltyFit(int idx) {
  int val = 0;

  vector<int> idx_same_match_fit;
  idx_same_match_fit.reserve(popsize * 2);
  for (int i = 0; i < popsize * 2; ++i) {
    if (i != idx && doubleFits[i] == doubleFits[idx]) {
      idx_same_match_fit.push_back(i);
    }
  }

  vector<vector<int>> strings_same_match_fit;
  strings_same_match_fit.reserve(idx_same_match_fit.size());

  for (int i : idx_same_match_fit) {
    doublePop[i].fillOutput(testSeq);
    strings_same_match_fit.push_back(testSeq);
  }

  doublePop[idx].fillOutput(testSeq);
  for (int i = 0; i < goalSeq.size(); ++i) {
    if (testSeq[i] != goalSeq[i]) {
      for (vector<int> vec : strings_same_match_fit) {
        if (testSeq[i] == vec[i]) {
          val++;
        }
      }
    }
  }

  return val;
}


double diversify(SDA &sda, SDA &best) { return 0.0; }

int random_range(int min, int max) {
  long range = (long)max - min + 1;
  // Get a random number and scale it to the range
  return min + (int)(lrand48() % range);
}

int calcRelativeFitness() {
  doubleRelativeFits.clear();
  doubleRelativeFits.resize(popsize * 2, 0.0);
  for (int idx = 0; idx < popsize * 2; ++idx) {
    for (int i = 0; i < tournSize; ++i) {
      int opponentIdx = random_range(0, popsize * 2 - 1);
      while (opponentIdx == idx) {
        opponentIdx = random_range(0, popsize * 2 - 1);
      }
      if (BIGGER_BETTER && doubleFits[idx] > doubleFits[opponentIdx]) {
        doubleRelativeFits[idx] += 1.0;
      } else if (doubleFits[idx] == doubleFits[opponentIdx]) {
       if (noveltyFits[idx] < noveltyFits[opponentIdx]) {
         doubleRelativeFits[idx] += 1.0;
       } else if (noveltyFits[idx] == noveltyFits[opponentIdx]) {
         doubleRelativeFits[idx] += 0.5;
       }
      }

    }
  }

  return 0;
}

int printExpStatsHeader(ostream &outp) {
  outp << left << setw(5) << "Run";
  outp << left << setw(4) << "RI";
  outp << left << setw(10) << "Mean";
  outp << left << setw(12) << "95% CI";
  outp << left << setw(10) << "SD";
  outp << left << setw(8) << "Best";
  outp << left << setw(12) << "% Correct";
  outp << left << setw(10) << "AvgStates";
  outp << endl;
  return 0;
}

double report(ofstream &outp, int run, int rptNum, bool biggerBetter) {
  vector<double> stats =
      calcStats<double>(fits, biggerBetter);  // {mean, stdDev, 95CI, best}
  multiStream printAndSave(cout, outp);

  // Calculate average number of states
  double avgStates = 0.0;
  for (int idx = 0; idx < popsize; ++idx) {
    avgStates += pop[idx].getNumStates();
  }
  avgStates /= popsize;

  printAndSave << left << setw(5) << run;
  printAndSave << left << setw(4) << rptNum;
  printAndSave << left << setw(10) << stats[0];
  printAndSave << left << setw(12) << stats[2];
  printAndSave << left << setw(10) << stats[1];
  printAndSave << left << setw(8) << stats[3];
  // Format percentage with % as part of the field
  ostringstream percStream;
  percStream << fixed << setprecision(2) << (stats[3] / seqLen) * 100 << "%";
  printAndSave << left << setw(12) << percStream.str();
  printAndSave << left << setw(10) << avgStates;
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

int matingEvent(bool biggerBetter, double cullingRate, int generation,
                ofstream &tournStats) {
  SDA child1;
  doublePop.clear();
  doubleFits.clear();

  for (int idx = 0; idx < popsize; ++idx) {
    doublePop.push_back(pop[idx]);
    doubleFits.push_back(fits[idx]);
  }
  for (int idx = 0; idx < popsize; ++idx) {
    child1 = pop[idx];
    child1.mutate(numMuts);
    doublePop.push_back(child1);
    doubleFits.push_back(fitness(child1));
  }
  for (int i = 0; i < popsize * 2; ++i) {
    noveltyFits[i] = calcNoveltyFit(i);
  }

  calcRelativeFitness();

  if (ROULETTE) {
    return -1;
  } else {
    vector<int> indices(doublePop.size());
    iota(indices.begin(), indices.end(), 0);

    // Sort indices based on compareFitness
    sort(indices.begin(), indices.end(), compareFitness);
    // Reorder doublePop and doubleFits based on sorted indices
    vector<SDA> sortedDoublePop;
    vector<double> sortedDoubleFits;
    vector<double> sortedDoubleRelativeFits;
    vector<int> sortedNoveltyFits;
    for (int idx : indices) {
      sortedDoublePop.push_back(doublePop[idx]);
      sortedDoubleFits.push_back(doubleFits[idx]);
      sortedDoubleRelativeFits.push_back(doubleRelativeFits[idx]);
      sortedNoveltyFits.push_back(noveltyFits[idx]);
    }
    doublePop = sortedDoublePop;
    doubleFits = sortedDoubleFits;
    doubleRelativeFits = sortedDoubleRelativeFits;
    if (generation % 20 == 0 || generation == 1 || generation % 101 == 0) {
      tournStats << "Event at Generation " << generation << "\n";
      tournStats << "Relative Fitness Values: ";
      for (int i = 0; i < doubleRelativeFits.size(); ++i) {
        if (i != 0) tournStats << ", ";
        tournStats << doubleRelativeFits[i];
      }
      double avg = accumulate(doubleRelativeFits.begin(),
                              doubleRelativeFits.end(), 0.0) /
                   doubleRelativeFits.size();
      tournStats << "\nAverage: " << avg << "\n";
      for (int i = 0; i < doubleFits.size(); ++i) {
        if (i != 0) tournStats << ", ";
        tournStats << doubleFits[i];
      }
      tournStats << "\nNoveltyFits\n";
      for (int i = 0; i < sortedNoveltyFits.size(); ++i) {
        if (i != 0) tournStats << ", ";
        tournStats << sortedNoveltyFits[i];
      }
      tournStats << "\n\n";
    }

    selectByRank();

    if (cullingRate > 0.0) {

      culling(cullingRate, BIGGER_BETTER);
    }
  }

  return 0;
}

int selectByRoulette() {
  vector<double> tempRelativeFits = doubleRelativeFits;  // Temporary copy
  vector<bool> selected(doublePop.size(), false);  // Track selected individuals
  int selectedCount = 0;

  while (selectedCount < popsize) {
    // Calculate total fitness of remaining candidates
    double total =
        accumulate(tempRelativeFits.begin(), tempRelativeFits.end(), 0.0);

    if (total == 0.0) {
      // Assign equal probability if all remaining fitnesses are zero
      for (double &relFit : tempRelativeFits) relFit = 1.0;
      total = tempRelativeFits.size();
    }

    // Build cumulative probabilities
    vector<double> cumulative;
    cumulative.reserve(tempRelativeFits.size());
    double sum = 0.0;
    for (double relFit : tempRelativeFits) {
      sum += relFit / total;
      cumulative.push_back(sum);
    }

    // Spin the wheel
    double randVal = drand48();
    auto it = lower_bound(cumulative.begin(), cumulative.end(), randVal);
    int idx = it - cumulative.begin();

    // Handle edge case
    if (idx >= tempRelativeFits.size()) idx = tempRelativeFits.size() - 1;

    // Ensure unique selection
    if (!selected[idx]) {
      pop[selectedCount].copy(doublePop[idx]);
      fits[selectedCount] = doubleFits[idx];
      selected[idx] = true;
      tempRelativeFits[idx] = 0.0;  // Prevent reselection
      selectedCount++;
    }
  }
  return 0;
}

int selectByRank() {
  pop.resize(popsize);
  copy(doubleFits.begin(), doubleFits.begin() + popsize, fits.begin());
  copy(doublePop.begin(), doublePop.begin() + popsize, pop.begin());

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
  if (BIGGER_BETTER) {
    return doubleRelativeFits[popIdx1] > doubleRelativeFits[popIdx2];
  }
  return doubleRelativeFits[popIdx1] < doubleRelativeFits[popIdx2];
}

int culling(double percentage, bool biggerBetter) {
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
  contestants.reserve(popsize);
  contestants = tournSelect(popsize, !biggerBetter);

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
  int finalNumStates = pop[bestIdx].getNumStates();
  printAndSave << "Final Number of states: " << finalNumStates << "\n";
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

int main(int argc, char *argv[]) {
  getArgs(argv);
  string pathToSeqs = "./Sequences.dat";
  string filename;
  ofstream runStats, expStats, readMe, tournStats;

  vector<double> bests;
  bests.reserve(runs);
  double expBestFit = (BIGGER_BETTER ? 0 : MAXFLOAT);

  initAlg(pathToSeqs);
  SDA expBestSDA = SDA(sdaStates, numChars, 2, seqLen, 1.5 * sdaStates);
    cmdLineIntro(cout);
    char dynamicMessage[20];
    sprintf(pathToOut,
      "./parallel/SQMatch on Seq%d with %dGens, %04dPS, %02dSt, "
      " %dTS, %dMuts, %02d%%CuR, %dCE/",
      seqNum, maxGens, popsize, sdaStates, tournSize, numMuts,
      (int)(cullingRate * 100), CULLING_EVERY);
    filesystem::create_directories(pathToOut);
    expStats.open(string(pathToOut) + "./exp.dat", ios::out);
    readMe.open(string(pathToOut) + "./read.me", ios::out);
    makeReadMe(readMe);
    readMe.close();

    int tmp;
  // for (int run = runarg; run < runarg + 1; ++run) {
  for (int run = 1; run < runs + 1; ++run) {
    initPop(run);
    char runNumStr[20];
    sprintf(runNumStr, "%02d", run);
    filename = string(pathToOut) + "run" + string(runNumStr) + ".dat";
    runStats.open(filename, ios::out);
    tournStats.open(string(pathToOut) + "./tourn" + string(runNumStr) + ".dat",
                    ios::out);
    printExpStatsHeader(cout);
    printExpStatsHeader(runStats);
    report(runStats, run, 0, BIGGER_BETTER);

    int gen = 1;
    double best = (BIGGER_BETTER ? 0 : MAXFLOAT);
    while (gen <= maxGens) {
      int genRatio = (int)(100 * (gen / (double)maxGens));
      if (gen % (int)(CULLING_EVERY * REPORT_EVERY) == 0 && genRatio > 30 &&
          genRatio < 95) {
        cout << "Culling at generation " << gen << endl;
        matingEvent(BIGGER_BETTER, cullingRate, gen, tournStats);

      } else {
        matingEvent(BIGGER_BETTER, 0.0, gen, tournStats);
      }

      if (gen % REPORT_EVERY == 0) {
        tmp = report(runStats, run, (int)gen / (REPORT_EVERY), BIGGER_BETTER);
        if ((BIGGER_BETTER && tmp > best) || (!BIGGER_BETTER && tmp < best)) {
          best = tmp;
        }
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
    tournStats.close();
  }

  ofstream best;
  best.open(string(pathToOut) + "./best.dat", ios::out);
  expReport(best, bests, expBestSDA, BIGGER_BETTER);
  best.close();
  // delete[] pop;
  cout << "Program Completed Successfully!" << endl;
  return 0;
}
