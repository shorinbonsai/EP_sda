#ifndef SDA_H
#define SDA_H

#include <iostream>
#include <vector>

using namespace std;

class SDA {
public:
    SDA();
    SDA(SDA &other);
    // The equals sign after the arguments provides the default value for the parameters.
    // explicit SDA(int numStates, int numChars, int maxRespLen, int outputLen, int initState = 0, bool verbose = false);
    explicit SDA(int numStates, int numChars, int maxRespLen, int outputLen, int initState = 0, bool verbose = false, int maxStates = 100);
    ~SDA();

    int randomize();
    int copy(SDA &other);
    int crossover(SDA &other);
    int mutate(int numMuts);
    int print(ostream &to);
    bool operator!=(SDA &other);
    bool operator==(SDA &other);
    vector<vector<vector<int> > > getResponses();

    int setOutputLen(int newLen);
    int fillOutput(vector<int> &output, bool printToo, ostream &outStream);
    vector<int> rtnOutput(bool printToo, ostream &outStream);

    // New functions for adding and deleting states
    int addState();
    int deleteState();

private:
    int create();       // Randomly initialize the SDA by setting initChar and filling the transitions and response vectors
    int maxStates;      // Maximum allowed states

    int initChar{};     // The initial character to drive the first transition
    int numStates{};    // The number of states in the SDA
    int initState{};    // The initial state for the SDA
    int curState{};     // The current state of the SDA
    int numChars{};     // The number of characters in the SDA's alphabet
    int maxRespLen{};   // The maximum length of the response vector
    int outputLen{};    // The length of the response vector expected from the SDA
    bool verbose{};     // Should information be displayed to console?

    /**
     * The transition vectors for each state.  For example transitions[3][1] = 5 means the SDA will transition from
     * state 3 to state 5 if a 1 is read.
     */
    vector<vector<int> > transitions;
    /**
     * The response vectors for each transition.  For example responses[3][1] = {1,0} means the SDA will append 0 and 1
     * to output if a 1 is read at state 3.
     */
    vector<vector<vector<int> > > responses;
};

#endif // SDA_H
