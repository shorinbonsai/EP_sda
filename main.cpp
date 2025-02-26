#include "main.h"

#include <iostream>

int main(int argc, char *argv[]) {

    getArgs(argv);
    string pathToSeqs = "./Sequences.dat";
    string filename;
    ofstream runStats, expStats, readMe;

    vector<double> bests;
    bests.reserve(runs);
    double expBestFit = (BIGGER_BETTER ? 0 : MAXFLOAT);
    SDA expBestSDA = SDA(sdaStates, numChars, 2, seqLen);
    goalSeq = getSequences(pathToSeqs)[seqNum];
    seqLen = (int) goalSeq.size();

    sprintf(pathToOut, "./SEMOut/SEMatch on Seq%d with %dGens, %04dPS, %02dSt/", seqNum, maxGens,
            popsize, sdaStates);
    mkdir(pathToOut, 0777);
    expStats.open(string(pathToOut) + "./exp.dat", ios::out);

    ofstream outFile(filename);

    for (int run = 1; run < runs; ++run){
        char runNumStr[20];
        sprintf(runNumStr, "%02d", run);
        filename = string(pathToOut) + "run" + string(runNumStr) + ".dat";
        Ep newEp(sdaStates, seqLen, goalSeq, maxGens, outFile, numChars, popsize, tournSize, seed,
       false);

        char runNumStr[20];
        sprintf(runNumStr, "%02d", run);
        runStats.open(filename, ios::out);
        printExpStatsHeader(cout);
        printExpStatsHeader(runStats);

        seed++;
    }
}