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
    FOption Option;
    FCompilerContext Context;

    std::string in;
    int N_MoleculesAllowed;
    std::string Name_Pseudo;
    float Float_Init = Numbers::GetFloatDefault(); // random initialized float
    int Int_Init = Numbers::GetIntDefault(); // random initialized int

    FWriter() {}

    // Initialize
    void LinkOptionContext(FOption InOption, FCompilerContext InContext) { Option = InOption; Context = InContext; }
    void SetUpDefaultVariables(int InN_MoleculesAllowed, std::string InName_Pseudo, float InFloat_Default, int InInt_Default);

    // StandardReaction
    void Initialize_StandardReaction(ofstream& ofs, std::string Type, std::string NameSpace_Pathway="");
    void SetUp_StandardReaction(ofstream& ofs, std::string Type, std::vector<FReaction *> ReactionSubList, std::string NameSpace_Pathway="");

    // EnzymeReaction
    void Initialize_EnzymeReaction(ofstream& ofs, std::string Type, std::string NameSpace_Pathway="");
    void SetUp_EnzymeReaction(ofstream& ofs, std::string Type, std::vector<FReaction *> ReactionSubList, std::string NameSpace_Pathway="");

    // Polymerase
    void Initialize_PolymeraseReaction(ofstream& ofs, FPolymerase* Polymerase);
    void SetUp_PolymeraseReaction(ofstream& ofs, FPolymerase* Polymerase, float Rate, std::string FreqBBFileName, std::string MaxLenFileName, int Idx_Pol, std::vector<int> Idx_Template, std::vector<int> Idx_TemplateSubset, std::vector<int> Idx_Target, std::vector<int> Idx_PolSub, std::vector<int> Idx_PolBB, int Threshold);
    void Polymerase_InitiationReaction(ofstream& ofs, FPolymerase* Polymerase);
    void Polymerase_ElongationReaction(ofstream& ofs, FPolymerase* Polymerase);
    void Polymerase_TerminationReaction(ofstream& ofs, FPolymerase* Polymerase);

    // TransporterReaction
    void Initialize_TransporterReaction(ofstream& ofs, std::string Type);
    void SetUp_TransporterReaction(ofstream& ofs, std::string Type, std::vector<FReaction *> ReactionSubList);

    // Spatial Simulation
    void Initialize_SpatialSimulation(ofstream& ofs);
    void SetUp_SpatialSimulation(ofstream& ofs);

    // Simulation
    void SimIdx();
    void SimExecutor();
    void SimModule(int Sim_Steps, int Sim_Resolution);
    void SimVis2D();

};

#endif /* LCC_WRITER_H */
