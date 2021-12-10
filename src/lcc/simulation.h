#ifndef LCC_SIMULATION_H
#define LCC_SIMULATION_H

#include <string>
#include <vector>
#include <map>

// #include "option.h"
// #include "node.h"
#include "context.h"
#include "datamanager.h"
#include "state.h"

class FDataset {
public:
    std::vector<std::string> Legend;
    std::vector<float> Data;

    void PrintLegend();
    void PrintData();

};

class FSimulation {
public:
    int N_SimSteps;

    void SetSimSteps(const int& SimSteps);

    void Init(FState& State, FDataset& Dataset, FDataManager& DataManager, const int& SimSteps);

    void Run(FState& State, FCompilerContext& Context, FDataset& Dataset);

};

#endif /* LCC_SIMULATION_H */
