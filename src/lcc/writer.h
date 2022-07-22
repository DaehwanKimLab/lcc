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

    std::string Str_Empty = "";
    std::string Str_Comma = ", ";
    std::string Str_Equal = " = ";
    std::string Str_LRB = "(";
    std::string Str_RRB = ")";

    int N_MoleculesAllowed;
    std::string Name_Pseudo;
    float Float_Init = Numbers::GetFloatDefault(); // random initialized float
    int Int_Init = Numbers::GetIntDefault(); // random initialized int

    FWriter() {}

    // Initialize
    void LinkOptionContext(FOption InOption, FCompilerContext InContext) { Option = InOption; Context = InContext; }
    void SetUpDefaultVariables(int InN_MoleculesAllowed, std::string InName_Pseudo, float InFloat_Default, int InInt_Default);

    // Writing Utility
    void DisplayWriterFunctionName(ofstream& ofs, std::string WriterName);

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

    // Polymerase type-specific functions
    void Initialize_PolymeraseReaction_DNAP(ofstream& ofs, std::vector<FMolecule *> ListOfPolymerases);
    void Initialize_PolymeraseReaction_RNAP(ofstream& ofs, std::vector<FMolecule *> ListOfPolymerases);
    void Initialize_PolymeraseReaction_Ribosome(ofstream& ofs, std::vector<FMolecule *> ListOfPolymerases, std::string Name_mRNASubIdx);
    void SetUp_Idx_mRNAInRNA(ofstream& ofs, std::string Name_mRNASubIdx);
    void SetUp_PolymeraseReaction_DNAP(ofstream& ofs, std::vector<FMolecule*> Polymerase);
    void SetUp_PolymeraseReaction_RNAP(ofstream& ofs, std::vector<FMolecule*> Polymerase);
    void SetUp_PolymeraseReaction_Ribosome(ofstream& ofs, std::vector<FMolecule*> Polymerase, std::string Name_mRNASubIdx);
    void Polymerase_InitiationReaction_DNAP(ofstream& ofs, std::vector<FMolecule*> Polymerases);
    void Polymerase_InitiationReaction_RNAP(ofstream& ofs, std::vector<FMolecule*> Polymerases);
    void Polymerase_InitiationReaction_Ribosome(ofstream& ofs, std::vector<FMolecule*> Polymerases, std::string Name_mRNASubIdx);
    void Polymerase_ElongationReaction_DNAP(ofstream& ofs, std::vector<FMolecule*> Polymerases);
    void Polymerase_ElongationReaction_RNAP(ofstream& ofs, std::vector<FMolecule*> Polymerases);
    void Polymerase_ElongationReaction_Ribosome(ofstream& ofs, std::vector<FMolecule*> Polymerases, std::string Name_mRNASubIdx);
    void Polymerase_TerminationReaction_DNAP(ofstream& ofs, std::vector<FMolecule*> Polymerases);
    void Polymerase_TerminationReaction_RNAP(ofstream& ofs, std::vector<FMolecule*> Polymerases);
    void Polymerase_TerminationReaction_Ribosome(ofstream& ofs, std::vector<FMolecule*> Polymerases, std::string Name_mRNASubIdx);


    // TransporterReaction
    void Initialize_TransporterReaction(ofstream& ofs, std::string Type);
    void SetUp_TransporterReaction(ofstream& ofs, std::string Type, std::vector<FReaction *> ReactionSubList);

    // Spatial Simulation
    void Initialize_SpatialSimulation(ofstream& ofs, int Map_Width, int Map_Height);
    void Initialize_ChromosomeSimulation(ofstream& ofs);
    void SetUp_SpatialSimulation(ofstream& ofs);
    void SetUp_ChromosomeSimulation(ofstream& ofs);

    // SimServer Utility
    void GenerateVisObjects(std::ofstream& ofs, int indents, std::string ObjectFamilyName, std::string N_Objects_Str);
    void GenerateOrganizationTree(std::ofstream& ofs, std::string Node, std::vector<std::string> Leaf); // not implemented yet

    // Simulation
    void SimIdx();  // to store long arrays to call from SimModule
    void SimExecutor();
    void SimState(int Map_Width, int Map_Height);
    void SimData();
    void SimModule(int Sim_Steps, int Sim_Resolution);
    void SimVis2D(int Sim_Steps_SteadyState);
    void SimServer(int Sim_Steps_SteadyState);

};

#endif /* LCC_WRITER_H */
