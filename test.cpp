#include "SDA.h"
#include <iostream>
#include <cassert>
#include <sstream>
#include <string>

using namespace std;

void testAddState()
{
    cout << "Testing addState..." << endl;
    SDA sda(2, 2, 2, 100, 0, true, 5);
    int initialStates = sda.getResponses().size();
    assert(initialStates == 2);

    // Add a state and check the count
    int result = sda.addState();
    assert(result == 0);
    assert(sda.getResponses().size() == 3);

    // Add until reaching maxStates
    sda.addState();          // 4
    sda.addState();          // 5 (maxStates)
    result = sda.addState(); // Should fail
    assert(result == -1);
    assert(sda.getResponses().size() == 5);

    cout << "testAddState passed." << endl;
}

void testDeleteState()
{
    cout << "Testing deleteState..." << endl;
    SDA sda(3, 2, 2, 100, 0, true, 100);
    assert(sda.getResponses().size() == 3);

    sda.deleteState();
    assert(sda.getResponses().size() == 2);

    sda.deleteState();
    assert(sda.getResponses().size() == 1);

    int result = sda.deleteState(); // Should fail
    assert(result == -1);

    cout << "testDeleteState passed." << endl;
}

void testOutputAfterAddDelete()
{
    cout << "Testing output generation after add/delete..." << endl;
    SDA sda(2, 2, 2, 10, 0, true, 5);

    // Generate output before changes
    vector<int> output = sda.rtnOutput(false, cout);
    assert(output.size() == 10);

    // Add a state and generate output
    sda.addState();
    output = sda.rtnOutput(false, cout);
    assert(output.size() == 10);

    // Delete a state and generate output
    sda.deleteState();
    output = sda.rtnOutput(false, cout);
    assert(output.size() == 10);

    cout << "testOutputAfterAddDelete passed." << endl;
}

void testInitStateAfterDelete()
{
    cout << "Testing initState adjustment after delete..." << endl;
    srand48(0); // Seed for deterministic test

    // Create SDA with initState 2 (third state)
    SDA sda(3, 2, 2, 10, 2, false, 100);

    // Delete state 0 (seed ensures lrand48() % 3 = 0)
    srand48(0);
    sda.deleteState();

    // Capture print output to check initState
    stringstream ss;
    sda.print(ss);
    string line;
    getline(ss, line);
    size_t pos = line.find(" <- ");
    int initState = stoi(line.substr(0, pos));

    // Original initState was 2, after deleting state 0, it should be 1
    assert(initState == 1);

    cout << "testInitStateAfterDelete passed." << endl;
}

void testAddStateConnections()
{
    cout << "Testing connections after addState..." << endl;
    SDA sda(2, 2, 2, 10, 0, false, 5);
    sda.addState(); // New state index is 2

    // Check if any transition points to the new state (2)
    stringstream ss;
    sda.print(ss);
    string line;
    bool found = false;
    while (getline(ss, line))
    {
        if (line.find("-> 2 [") != string::npos)
        {
            found = true;
            break;
        }
    }
    assert(found);

    cout << "testAddStateConnections passed." << endl;
}

void testDeleteStateTransitions()
{
    cout << "Testing transitions after deleteState..." << endl;
    srand48(0); // Seed for deterministic test

    SDA sda(3, 2, 2, 10, 0, false, 100);
    srand48(0);
    sda.deleteState(); // Deletes state 0

    // Check all transitions are within valid state range
    stringstream ss;
    sda.print(ss);
    string line;
    while (getline(ss, line))
    {
        if (line.find("->") != string::npos)
        {
            size_t arrowPos = line.find("->");
            string statePart = line.substr(arrowPos + 3);
            size_t spacePos = statePart.find(' ');
            int targetState = stoi(statePart.substr(0, spacePos));
            assert(targetState >= 0 && targetState < 2); // New numStates is 2
        }
    }

    cout << "testDeleteStateTransitions passed." << endl;
}

int main()
{
    testAddState();
    testDeleteState();
    testOutputAfterAddDelete();
    testInitStateAfterDelete();
    testAddStateConnections();
    testDeleteStateTransitions();

    cout << "All tests passed!" << endl;
    return 0;
}