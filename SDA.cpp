#include "SDA.h"

SDA::SDA(int numStates, int numChars, int maxRespLen, int outputLen,
         int initState, bool verbose, int maxStates) {
  initChar = -1;
  this->numStates = numStates;
  this->initNumStates = numStates;
  this->initState = initState;
  this->numChars = numChars;
  this->maxRespLen = maxRespLen;
  this->outputLen = outputLen;
  this->verbose = verbose;
  this->maxStates = maxStates;

  transitions.reserve(numStates);
  for (vector<int> v : transitions) {
    v.reserve(numChars);
  }

  responses.reserve(numStates);
  for (vector<vector<int>> vec1 : responses) {
    vec1.reserve(numChars);
    for (vector<int> vec2 : vec1) {
      vec1.reserve(maxRespLen);
    }
  }
  create();
  if (verbose) cout << "SDA made with " << numStates << " numStates." << endl;
}

SDA::SDA() : SDA(10, 2, 2, 1000) {}

SDA::SDA(const SDA &other)
    : maxStates(other.maxStates),
      initChar(other.initChar),
      numStates(other.numStates),
      initNumStates(other.initNumStates),
      initState(other.initState),
      curState(other.curState),
      numChars(other.numChars),
      maxRespLen(other.maxRespLen),
      outputLen(other.outputLen),
      verbose(other.verbose),
      transitions(other.transitions),
      responses(other.responses) {
  if (verbose) cout << "SDA Copied (via copy constructor)." << endl;
}

SDA::~SDA() = default;

int SDA::create() {
  initChar = (int)lrand48() % numChars;

  vector<int> oneState;
  oneState.reserve(numChars);
  for (int state = 0; state < numStates; ++state) {
    oneState.clear();
    for (int val = 0; val < numChars; ++val) {
      oneState.push_back((int)lrand48() % numStates);
    }
    transitions.push_back(oneState);
  }

  vector<int> oneResponse;
  oneResponse.reserve(maxRespLen);
  vector<vector<int>> oneStateResps;
  oneStateResps.reserve(numChars);
  int respSize;
  for (int state = 0; state < numStates; ++state) {
    oneStateResps.clear();
    for (int trans = 0; trans < numChars; ++trans) {
      oneResponse.clear();
      respSize = (int)lrand48() % maxRespLen + 1;
      for (int val = 0; val < respSize; ++val) {
        oneResponse.push_back((int)lrand48() % numChars);
      }
      oneStateResps.push_back(oneResponse);
    }
    responses.push_back(oneStateResps);
  }
  if (verbose) cout << "SDA initialized." << endl;
  return 0;
}

int SDA::setOutputLen(int newLen) {
  this->outputLen = newLen;
  return 0;
}

int SDA::randomize() {
  if (initChar < 0) {
    cout << "Error in SDA Class: randomize(): this SDA has not been "
            "initialized.";
    return -1;
  }

  initChar = (int)lrand48() % numChars;

  vector<int> oneResponse;
  oneResponse.reserve(maxRespLen);
  int respLen;
  numStates = initNumStates;

  for (int state = 0; state < numStates; ++state) {
    for (int trans = 0; trans < numChars; ++trans) {
      transitions[state][trans] = (int)lrand48() % numStates;
      oneResponse.clear();
      respLen = (int)lrand48() % maxRespLen + 1;
      for (int val = 0; val < respLen; ++val) {
        oneResponse.push_back((int)lrand48() % numChars);
      }
      responses[state][trans] = oneResponse;
    }
  }
  if (verbose) cout << "SDA Randomized." << endl;
  return 0;
}

int SDA::copy(SDA &other) {
  // Directly copy all member variables
  maxStates = other.maxStates;
  initChar = other.initChar;
  numStates = other.numStates;
  initState = other.initState;
  numChars = other.numChars;
  maxRespLen = other.maxRespLen;
  outputLen = other.outputLen;
  verbose = other.verbose;

  // Use vector assignments (automatic deep copy)
  transitions = other.transitions;
  responses = other.responses;
  if (verbose) cout << "SDA Copied." << endl;
  return 0;
}

int SDA::twoPointCrossover(SDA &other, int firstCP, int secondCP) {
  if (initChar < 0) {
    cout << "Error in SDA Class: twoPointCrossover(...): this SDA has not been "
            "initialized.";
    return -1;
  }
  if (other.initChar < 0) {
    cout << "Error in SDA Class: twoPointCrossover(...): other SDA has not "
            "been initialized.";
    return -1;
  }
  if (numStates != other.numStates) {
    cout << "Error in SDA Class: twoPointCrossover(...): the two SDAs have a "
            "different numStates.";
    return 1;
  }
  if (numChars != other.numChars) {
    cout << "Error in SDA Class: twoPointCrossover(...): the two SDAs have a "
            "different numChars.";
    return 1;
  }
  if (maxRespLen != other.maxRespLen) {
    cout << "Error in SDA Class: twoPointCrossover(...): the two SDAs have a "
            "different maxRespLen.";
    return 1;
  }

  int crossStart, crossEnd, swapInt;
  vector<int> swapVec;
  swapVec.reserve(numChars);

  if (firstCP == -1 && secondCP == -1) {
    do {
      crossStart = (int)lrand48() % (numStates + 1);
      crossEnd = (int)lrand48() % (numStates + 1);
      if (crossStart > crossEnd) {
        swapInt = crossStart;
        crossStart = crossEnd;
        crossEnd = swapInt;
      }
    } while (crossStart == crossEnd);
  } else {
    crossStart = firstCP;
    crossEnd = secondCP;
  }

  if (crossStart == 0) {
    swapInt = initChar;
    initChar = other.initChar;
    other.initChar = swapInt;
  }

  for (int state = crossStart; state < crossEnd; state++) {
    swapVec = transitions.at(state);
    transitions.at(state) = other.transitions.at(state);
    other.transitions.at(state) = swapVec;
    for (int trans = 0; trans < numChars; trans++) {
      swapVec = responses.at(state).at(trans);
      responses.at(state).at(trans) = other.responses.at(state).at(trans);
      other.responses.at(state).at(trans) = swapVec;
    }
  }
  return 0;
}

int SDA::oneStateCrossover(SDA &other) {
  if (initChar < 0) {
    cout << "Error in SDA Class: oneStateCrossover(...): this SDA has not been "
            "initialized.";
    return -1;
  }
  if (other.initChar < 0) {
    cout << "Error in SDA Class: oneStateCrossover(...): other SDA has not "
            "been initialized.";
    return -1;
  }
  if (numStates != other.numStates) {
    cout << "Error in SDA Class: oneStateCrossover(...): the two SDAs have a "
            "different numStates.";
    return 1;
  }
  if (numChars != other.numChars) {
    cout << "Error in SDA Class: oneStateCrossover(...): the two SDAs have a "
            "different numChars.";
    return 1;
  }
  if (maxRespLen != other.maxRespLen) {
    cout << "Error in SDA Class: oneStateCrossover(...): the two SDAs have a "
            "different maxRespLen.";
    return 1;
  }

  int crossState, swapInt;
  vector<int> swapVec;

  crossState = (int)lrand48() % numStates;
  if (crossState == 0) {
    swapInt = initChar;
    initChar = other.initChar;
    other.initChar = swapInt;
  }

  swapVec = transitions.at(crossState);
  transitions.at(crossState) = other.transitions.at(crossState);
  other.transitions.at(crossState) = swapVec;
  for (int trans = 0; trans < numChars; trans++) {
    swapVec = responses.at(crossState).at(trans);
    responses.at(crossState).at(trans) =
        other.responses.at(crossState).at(trans);
    other.responses.at(crossState).at(trans) = swapVec;
  }
  return 0;
}

// new function to add a state
int SDA::addState() {
  if (numStates >= maxStates) {
    if (verbose)
      cout << "Error: Cannot add state. Maximum states reached." << endl;
    return -1;
  }

  vector<int> new_trans(numChars);
  vector<vector<int>> new_resp(numChars);

  for (int c = 0; c < numChars; ++c) {
    new_trans[c] = lrand48() % numStates;  // Transition to existing states
    int resp_len = lrand48() % maxRespLen + 1;
    vector<int> resp(resp_len);
    for (int i = 0; i < resp_len; ++i) {
      resp[i] = lrand48() % numChars;
    }
    new_resp[c] = resp;
  }

  transitions.push_back(new_trans);
  responses.push_back(new_resp);
  numStates++;

  int connect_state = lrand48() % numStates;
  int connect_char = lrand48() % numChars;
  transitions[connect_state][connect_char] = numStates - 1;

  if (verbose) {
    cout << "Added state " << (numStates - 1) << ". Connected via state "
         << connect_state << ", char " << connect_char << endl;
  }
  return 0;
}

// new function to delete a state
int SDA::deleteState() {
  if (numStates <= 5) {
    if (verbose)
      cout << "Error: Cannot delete state. Minimum states reached." << endl;
    return -1;
  }

  int delete_state = lrand48() % numStates;

  // Update transitions pointing to delete_state to follow its transitions
  for (int s = 0; s < numStates; ++s) {
    for (int c = 0; c < numChars; ++c) {
      if (transitions[s][c] == delete_state) {
        transitions[s][c] = transitions[delete_state][c];
      }
    }
  }

  // Adjust transitions and handle any remaining references to delete_state
  for (int s = 0; s < numStates; ++s) {
    for (int c = 0; c < numChars; ++c) {
      if (transitions[s][c] == delete_state) {
        transitions[s][c] = 0;  // Replace with a valid state
      } else if (transitions[s][c] > delete_state) {
        transitions[s][c]--;
      }
    }
  }
  for (int s = 0; s < numStates; ++s) {
    for (int c = 0; c < numChars; ++c) {
      if (transitions[s][c] >=
          numStates) {  // After deletion, numStates is decremented
        transitions[s][c] = numStates - 1;
      }
    }
  }

  // Update initState if needed
  if (initState == delete_state) {
    initState = (delete_state == numStates - 1) ? 0 : delete_state + 1;
  } else if (initState > delete_state) {
    initState--;
  }

  // Remove the state
  transitions.erase(transitions.begin() + delete_state);
  responses.erase(responses.begin() + delete_state);
  numStates--;

  if (verbose) {
    cout << "Deleted state " << delete_state << ". New numStates: " << numStates
         << endl;
  }
  return 0;
}

int SDA::getNumStates() const { return numStates; }

int SDA::mutate(int numMuts) {
  if (initChar < 0) {
    cout << "Error in SDA Class: mutate(...): this SDA has not been "
            "initialized.";
    return -1;
  }
  // Randomly choose number of mutations if numMuts is 0 (1-3 mutations)
  int tmpMut = numMuts;
  if (numMuts == 0) {
    tmpMut = (int)lrand48() % 3 + 1;
  }

  int mutPt, respSize;
  vector<int> oneResponse;

  for (int mut = 0; mut < tmpMut; ++mut) {
    double randVal = drand48();
    if (randVal < 0.04) {  // 4% chance of mutating initial character
      initChar = (int)lrand48() % numChars;
      if (verbose) {
        cout << "Completed mutation on the SDA's initial character." << endl;
      }
      return 0;
    } else if (0.04 <= randVal && randVal < 0.12) {
      deleteState();
    } else if (0.12 <= randVal && randVal < 0.20) {
      addState();
    } else {
      mutPt = (int)lrand48() % numStates;
      int transNum = (int)lrand48() % numChars;

      if ((int)lrand48() % 2 == 0) {  // Mutate transition (50%)
        transitions.at(mutPt).at(transNum) = (int)lrand48() % numStates;
      } else {  // Mutate response (50%)
        oneResponse.clear();
        respSize = (int)lrand48() % maxRespLen + 1;
        for (int i = 0; i < respSize; ++i) {
          oneResponse.push_back((int)lrand48() % numChars);
        }
        responses.at(mutPt).at(transNum) = oneResponse;
      }
    }
  }
  return 0;
}

//!!!!!!!!!!!!!!!!!!!!!! Static version (does fixed amount of mutations for each
//! type)
// int SDA::mutate(int transMuts, int respMuts) {
//   if (initChar < 0) {
//     cout << "Error in SDA Class: mutate(...): this SDA has not been "
//             "initialized.";
//     return -1;
//   }

//   if (drand48() < 0.1) {  // 10% chance of mutating initial character
//     initChar = (int)lrand48() % numChars;
//   }

//   int mutPt, respSize;
//   vector<int> oneResponse;

//   for (int mut = 0; mut < transMuts; ++mut) {
//     mutPt = (int)lrand48() % numStates;
//     int transNum = (int)lrand48() % numChars;
//     transitions.at(mutPt).at(transNum) = (int)lrand48() % numStates;
//   }

//   for (int mut = 0; mut < respMuts; ++mut) {
//     mutPt = (int)lrand48() % numStates;
//     int transNum = (int)lrand48() % numChars;
//     oneResponse.clear();
//     respSize = (int)lrand48() % maxRespLen + 1;
//     for (int i = 0; i < respSize; ++i) {
//       oneResponse.push_back((int)lrand48() % numChars);
//     }
//     responses.at(mutPt).at(transNum) = oneResponse;
//   }
//   return 0;
// }

int SDA::fillOutput(vector<int> &output, bool printToo, ostream &outStream) {
  if (initChar < 0) {
    cout << "Error in SDA Class: fillOutput(...): this SDA has not been "
            "initialized.";
    return -1;
  }
  if (output.capacity() < outputLen) {
    cout << "Error in SDA Class: fillOutput(...): output vector capacity is "
         << output.capacity();
    cout << " but the outputLen is " << outputLen << "." << endl;
    return -1;
  }

  int headIdx = 0;
  int tailIdx = 0;
  curState = initState;
  output[headIdx++] = initChar;
  if (printToo) outStream << initChar;

  while (headIdx < outputLen) {
    for (int val : responses[curState][output[tailIdx]]) {
      if (headIdx < outputLen) {
        output[headIdx++] = val;
        if (printToo) outStream << " " << val;
      }
    }
    curState = transitions[curState][output[tailIdx++]];
  }
  if (printToo) outStream << endl;
  return 0;
}

vector<int> SDA::rtnOutput(bool printToo, ostream &outStream) {
  if (initChar < 0) {
    cout << "Error in SDA Class: rtnOutput(...): this SDA has not been "
            "initialized.";
    return {-1};
  }

  vector<int> output(outputLen);
  fillOutput(output, printToo, outStream);
  return output;
}

int SDA::printSDA(ostream &outStream = cout) {
  if (initChar < 0) {
    cout << "Error in SDA Class: printSDA(...): this SDA has not been "
            "initialized.";
    return -1;
  }

  outStream << initState << " <- " << initChar << endl;
  for (int state = 0; state < numStates; ++state) {
    for (int trans = 0; trans < numChars; ++trans) {
      outStream << state << " + " << trans << " -> "
                << transitions[state][trans] << " [";
      for (int vec : responses[state][trans]) {
        outStream << " " << vec;
      }
      outStream << " ]" << endl;
    }
  }
  if (verbose) cout << "SDA Printed." << endl;
  return 0;
}