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
    void Initialize_PolymeraseReaction_Matrix(ofstream& ofs, std::vector<std::vector<FMolecule *>> PolymeraseTypes);
    void Initialize_PolymeraseReaction_Index(ofstream& ofs, std::string Process);
    void SetUp_PolymeraseReaction_Matrix(ofstream& ofs, std::vector<std::vector<FMolecule *>> PolymeraseTypes);
    void SetUp_PolymeraseReaction_Index(ofstream& ofs, std::vector<FMolecule *> Polymerases, int Threshold);
    void Polymerase_InitiationReaction(ofstream& ofs, std::vector<FMolecule*> Polymerases);
    void Polymerase_ElongationReaction(ofstream& ofs, std::vector<FMolecule*> Polymerases);
    void Polymerase_TerminationReaction(ofstream& ofs, std::vector<FMolecule*> Polymerases);

    // Polymerase type-specific functions
    void Initialize_PolymeraseReaction_DNAP(ofstream& ofs, std::vector<FMolecule *> ListOfPolymerases);
    void Initialize_PolymeraseReaction_RNAP(ofstream& ofs, std::vector<FMolecule *> ListOfPolymerases);
    void Initialize_PolymeraseReaction_Ribosome(ofstream& ofs, std::vector<FMolecule *> ListOfPolymerases);
    void SetUp_Idx_mRNAInRNA(ofstream& ofs);
    void SetUp_PolymeraseReaction_DNAP(ofstream& ofs, std::vector<FMolecule*> Polymerase);
    void SetUp_PolymeraseReaction_RNAP(ofstream& ofs, std::vector<FMolecule*> Polymerase);
    void SetUp_PolymeraseReaction_Ribosome(ofstream& ofs, std::vector<FMolecule*> Polymerase);

    // TransporterReaction
    void Initialize_TransporterReaction(ofstream& ofs, std::string Type);
    void SetUp_TransporterReaction(ofstream& ofs, std::string Type, std::vector<FReaction *> ReactionSubList);

    // Spatial Simulation
    void Initialize_SpatialSimulation(ofstream& ofs);
    void Initialize_ChromosomeSimulation(ofstream& ofs);
    void SetUp_SpatialSimulation(ofstream& ofs);
    void SetUp_ChromosomeSimulation(ofstream& ofs);

    // Simulation
    void SimIdx();
    void SimExecutor();
    void SimModule(int Sim_Steps, int Sim_Resolution);
    void SimVis2D();
    void SimServer();

};

#endif /* LCC_WRITER_H */
