#ifndef LCC_WRITER_H
#define LCC_WRITER_H

#include <string>
#include <vector>
#include <map>
#include <iterator>

#include <iostream>
#include <sstream>
#include <fstream>

#include "option.h"
#include "node.h"
#include "context.h"
#include "util.h"

using namespace std;

class FWriter {
public:
    FOption * pOption;
    FCompilerContext * pContext;

    std::string in;
    int N_MoleculesAllowed;
    std::string Name_Pseudo;
    float Float_Init = Numbers::GetFloatDefault(); // random initialized float
    int Int_Init = Numbers::GetIntDefault(); // random initialized int

    FWriter() {}

    void LinkOptionContext(FOption * InpOption, FCompilerContext * InpContext) { pOption = InpOption; pContext = InpContext; }
    void SetUpDefaultVariables(int InN_MoleculesAllowed, std::string InName_Pseudo, float InFloat_Default, int InInt_Default);

    void Print_InitializeStandardReaction(ofstream& ofs, std::string Type);
    void Print_SetUpStandardReaction(ofstream& ofs, std::string Type, std::vector<const FReaction *> ReactionSubList);
    void Print_InitializeEnzymeReaction(ofstream& ofs, std::string Type);
    void Print_SetUpEnzymeReaction(ofstream& ofs, std::string Type, std::vector<const FReaction *> ReactionSubList);
    void Print_InitializePolymeraseReaction(ofstream& ofs, const FPolymerase* Polymerase);
    void Print_SetUpPolymeraseReaction(ofstream& ofs, const FPolymerase* Polymerase, float Rate, std::string FreqBBFileName, std::string MaxLenFileName, int Idx_Pol, std::vector<int> Idx_Template, std::vector<int> Idx_TemplateSubset, std::vector<int> Idx_Target, std::vector<int> Idx_PolSub, std::vector<int> Idx_PolBB, int Threshold);
    void Print_InitiationReaction(ofstream& ofs, const FPolymerase* Polymerase);
    void Print_ElongationReaction(ofstream& ofs, const FPolymerase* Polymerase);
    void Print_TerminationReaction(ofstream& ofs, const FPolymerase* Polymerase);
    void Print_Initialize_SpatialSimulation(ofstream& ofs);
    void Print_SetUp_SpatialSimulation(ofstream& ofs);
    void SimIdx();
    void SimModule(int Sim_Steps, int Sim_Resolution);
    void SimVis2D();

};

#endif /* LCC_WRITER_H */
