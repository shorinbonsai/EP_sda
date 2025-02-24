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
}