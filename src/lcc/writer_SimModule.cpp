#include "writer.h"
#include <algorithm>

using namespace std;

void FWriter::SimModule(int Sim_Steps, int Sim_Resolution)
{
    std::cout << "Generating simulation module..." << std::endl;

    // write SimModule.py
    std::ofstream ofs(Option.SimModuleFile.c_str());
    std::string endl = "\n";

    // IMPORT
    ofs << "import os, sys" << endl;
    ofs << endl;

    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "print('\\n%%%%%% Run " << Option.SimExecutorFile << " to execute simulation. %%%%%%')" << endl;
    ofs << in+ "sys.exit()" << endl;
    ofs << endl;

    ofs << "import numpy as np" << endl;
    // ofs << "import tensorflow as tf" << endl;
    ofs << "from datetime import datetime" << endl;
    ofs << "import csv" << endl;
    ofs << "import SimFunctions as SimF" << endl;
    // ofs << "import SimIdx as idx" << endl;
    ofs << "from argparse import ArgumentParser" << endl;
    ofs << endl;

    //TODO: Take options from SimModule cmd line
    ofs << "# Temporary global variables" << endl;
    ofs << "N_SimSteps = " << Sim_Steps << endl;
    ofs << "SimStepTimeResolution = " << Sim_Resolution << endl;
    ofs << "DisplayCount = 0   # 0: Display Counts only, 1: Display both Counts & dCounts" << endl;
    ofs << "Unit = 'nM'   # nM or uM supported for now" << endl;
    ofs << endl;

    ofs << "def ListToStr(List, Delimiter):" << endl;
    ofs << in+ "List_Str = ''" << endl;
    ofs << in+ "if List:" << endl;
    ofs << in+ in+ "for item in List:" << endl;
    ofs << in+ in+ in+ "List_Str += item + '_'" << endl;
    ofs << in+ "return List_Str" << endl;
    ofs << endl;

    // class FState
    ofs << "class FState:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "self.Vol = 0" << endl;
    ofs << endl;

    ofs << in+ in+ "# State Arrays" << endl;

    ofs << in+ in+ "# Temporary legend attribute" << endl;
    ofs << in+ in+ "self.Mol_Names = list()" << endl;
    ofs << in+ in+ "self.Mol_Name2Idx = dict()" << endl;
    ofs << in+ in+ "self.Legends = list()" << endl;
    ofs << endl;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto ObjLoc = Context.GetSubList_LocationList("Compartment");

    int MatrixSize = 1;
    if (!ObjLoc.empty()) {
        MatrixSize = 0;
        std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
        for (auto& UniqueName : ObjUniqueNames) {
            int Count = int(Context.GetInitialCountByName_CountList(UniqueName));
            MatrixSize += Count;
        }
    }
    ofs << in+ in+ "self.Count_All = np.zeros([" << MatrixSize << ", " << Context.MoleculeList.size() << "])" << endl;
    ofs << in+ in+ "self.dCount_All = np.zeros([" << MatrixSize << ", " << Context.MoleculeList.size() << "])" << endl;
    ofs << endl;

    if (!Context.LocationList.empty()) {
        Initialize_SpatialSimulation(ofs);
    }

    // for standard reactions
    std::vector<std::string> ReactionTypes;
    ReactionTypes.emplace_back("Standard_Unregulated");
    ReactionTypes.emplace_back("Standard_Inhibition_Allosteric");
    ReactionTypes.emplace_back("Standard_Activation_Allosteric");

    std::vector<std::string> StandardReactionTypes = ReactionTypes; // to reuse later

    for (auto& Type : ReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Initialize_StandardReaction(ofs, Type);
        }
    }

    // for Enz reactions
    ReactionTypes.clear();
    ReactionTypes.push_back("Enz_Standard_Unregulated");
    ReactionTypes.push_back("Enz_Standard_Inhibition_Allosteric");
    ReactionTypes.push_back("Enz_Standard_Activation_Allosteric");
    ReactionTypes.push_back("Enz_MichaelisMenten_Unregulated");
    ReactionTypes.push_back("Enz_MichaelisMenten_Inhibition_Allosteric");
    ReactionTypes.push_back("Enz_MichaelisMenten_Inhibition_Competitive");
    ReactionTypes.push_back("Enz_MichaelisMenten_Activation_Allosteric");

    std::vector<std::string> EnzReactionTypes = ReactionTypes; // to reuse later

    for (auto& Type : ReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Initialize_EnzymeReaction(ofs, Type);
        }
    }

    // for polymerase reactions (Template-based)
    std::vector<FPolymerase *> PolymeraseList = Context.GetList_Polymerase_MoleculeList();
//    std::vector<std::string> PolymeraseNames = Context.GetNames_PolymeraseList(PolymeraseList);
    std::vector<FPolymeraseReaction *> PolymeraseReactionList = Context.GetList_Polymerase_ReactionList();

    if (!PolymeraseList.empty()) {
        for (auto& Polymerase : PolymeraseList) {
            Initialize_PolymeraseReaction(ofs, Polymerase);
        }
    }

    // for Transporter reactions
    ReactionTypes.clear();
    ReactionTypes.push_back("Transporter_Unregulated");
    ReactionTypes.push_back("Transporter_Inhibition_Allosteric");
    ReactionTypes.push_back("Transporter_Inhibition_Competitive");
    ReactionTypes.push_back("Transporter_Activation_Allosteric");

    std::vector<std::string> TransporterReactionTypes = ReactionTypes; // to reuse later

    for (auto& Type : ReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Initialize_TransporterReaction(ofs, Type);
        }
    }

    ofs << in+ "def Initialize(self):" << endl;
    ofs << in+ in+ "self.Vol = 1" << endl;
    ofs << endl;

    // Print SetUp_SpatialReaction for all spatial simulation (to be updated)
    if (!Context.LocationList.empty()) {
        SetUp_SpatialSimulation(ofs);
    }

    // Print SetUp_StandardReaction for each Reaction Type
    for (auto& Type : StandardReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            SetUp_StandardReaction(ofs, Type, ReactionSubList);
        }
    }

    // Print SetUp_EnzymeReaction for each Reaction Type
    for (auto& Type : EnzReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            SetUp_EnzymeReaction(ofs, Type, ReactionSubList);
        }
    }

    if (!PolymeraseList.empty()) {
        for (auto& Polymerase : PolymeraseList) {
            std::string PolymeraseName = Polymerase->Name;
            float Rate = Polymerase->Rate;
            int Idx_Pol = Context.GetIdxByName_MoleculeList(PolymeraseName);
            std::vector<int> Idx_Template;
            std::vector<int> Idx_Target;
            std::vector<int> Idx_TemplateSubset;
            std::vector<int> Idx_PolSub = Context.GetIdx_PolymeraseReactionSubstrate_ByPolymeraseName_MoleculeList(PolymeraseName);
            std::vector<int> Idx_PolBB;
            std::vector<std::string> BuildingBlockNames;
            std::string MaxLenFileName;
            std::string FreqBBFileName;
            int Threshold;

            for (auto& reaction : Context.GetList_Polymerase_ReactionList()){
                if (reaction->Name == PolymeraseName){
                    for (auto& BuildingBlock : reaction->BuildingBlocks){
                        BuildingBlockNames.push_back(BuildingBlock);
                    }
                    Idx_PolBB = Context.GetIdxByStrList_MoleculeList(BuildingBlockNames);
                }
            }

            std::vector<std::vector<int>> StoichMatrix_PolymeraseReaction = Context.GetStoichiometryMatrix_PolymeraseReaction(PolymeraseReactionList);

            if (Polymerase->Process == "DNAReplication") {
                MaxLenFileName = "./Database/Len_ChromosomesInGenome.npy";
                FreqBBFileName = "./Database/Freq_NTsInChromosomesInGenome.npy";
                Idx_Template   = Context.GetIdxListFromMoleculeList("Chromosome");
                Idx_Target     = Context.GetIdxListFromMoleculeList("Chromosome");
                Threshold      = 50;  // The number of polymerases required for a functional unit to initiate the polymerase reaction

                // Idx_TemplateSubset
                for (int i = 0; i < Idx_Template.size(); i++) {
                    Idx_TemplateSubset.push_back(i);
                }
            }

            else if (Polymerase->Process == "RNATranscription") {
                MaxLenFileName = "./Database/Len_RNAs.npy";
                FreqBBFileName = "./Database/Freq_NTsInRNAs.npy";
                Idx_Template   = Context.GetIdxListFromMoleculeList("Gene");
                Idx_Target     = Context.GetIdxListFromMoleculeList("RNA"); // update to link gene to RNA matching when importing and generating this list
                Threshold      = 1;

                // Idx_TemplateSubset
                for (int i = 0; i < Idx_Template.size(); i++) {
                    Idx_TemplateSubset.push_back(i);
                }
            }

            else if (Polymerase->Process == "ProteinTranslation") {
                MaxLenFileName = "./Database/Len_Proteins.npy";
                FreqBBFileName = "./Database/Freq_AAsInProteins.npy";
                Idx_Template   = Context.GetIdxListFromMoleculeList("mRNA");
                Idx_Target     = Context.GetIdxListFromMoleculeList("Protein"); // update to link mRNA to protein matching
                Threshold      = 1;

                Idx_TemplateSubset = Context.GetIdxOfStrListFromStrList(Context.GetNameListFromMoleculeList("mRNA"), Context.GetNameListFromMoleculeList("RNA"));
            }

            else { // add exception handling
            }

            // debugging
            std::cout << Polymerase->Process;
            if (Polymerase->Process == "DNAReplication") {
                std::cout << "\t";
            }
            std::cout << "\t | Idx_Template.size(): " << Idx_Template.size() << "\t Idx_Target.size(): " << Idx_Target.size() << endl;
            Utils::Assertion((Idx_Template.size() == Idx_Target.size()), "# of Target and Target do not match for Polymerase: " + PolymeraseName);

            SetUp_PolymeraseReaction(ofs, Polymerase, Rate, FreqBBFileName, MaxLenFileName, Idx_Pol, Idx_Template, Idx_TemplateSubset, Idx_Target, Idx_PolSub, Idx_PolBB, Threshold);
        }
    }

    // Print SetUp_TransporterReaction for each Reaction Type
    for (auto& Type : TransporterReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            SetUp_TransporterReaction(ofs, Type, ReactionSubList);
        }
    }

    // for legends
    std::vector<std::string> MolNames = Context.GetNames_MoleculeList();
    ofs << in+ in+ "self.Mol_Names = [" << Utils::JoinStr2Str(MolNames) << "]" << endl;
    ofs << in+ in+ "self.Mol_Name2Idx = {";
    int i = 0;
    for (auto& MolName : MolNames) {
        ofs << "'" << MolName << "' : np.array([" << i << "]), ";
        i++;
    } ofs << "}" << endl;
    ofs << in+ in+ "self.Legends = ['SimStep', 'Vol', " << Utils::JoinStr2Str(MolNames) << "]" << endl;
    ofs << endl;

    // Initialize all molecule counts
    std::vector<int> Idx_Mol = Context.GetIdxListFromMoleculeList("Molecule");
    std::vector<float> InitialCount_Molecules;
    std::vector<float> MolarityFactor_Molecules;

    for (auto& molecule : Context.MoleculeList) {
        float Count = Context.GetInitialCountByName_CountList(molecule->Name);
        float MolarityFactor = Context.GetMolarityFactorByName_CountList(molecule->Name);

        // TODO: special consideration for 'L' and 'qL' for chemotaxis
        if ((molecule->Name == "L") || (molecule->Name == "qL")) { Count = 0; }

        InitialCount_Molecules.push_back(Count);
        MolarityFactor_Molecules.push_back(MolarityFactor);
    }

    ofs << in+ in+ "Idx_Mol = np.array([";
    if (ObjLoc.empty()) {
        ofs << Utils::JoinInt2Str_Idx(Idx_Mol);
    } else {
        for (auto& objLoc : ObjLoc) {
            ofs << "[" << Utils::JoinInt2Str_Idx(Idx_Mol) << "], ";
        }
    }
    ofs << "], ndmin=2)" << endl;

    ofs << in+ in+ "Count_Mol = np.array([";
    if (ObjLoc.empty()) {
        ofs << Utils::JoinFloat2Str(InitialCount_Molecules);
    } else {
        for (auto& objLoc : ObjLoc) {
            ofs << "[" << Utils::JoinFloat2Str(InitialCount_Molecules) << "], ";
        }
    }
    ofs << "])" << endl;

    ofs << in+ in+ "MolarityFactor_Mol = np.array([" << Utils::JoinFloat2Str(MolarityFactor_Molecules) << "])" << endl;
    ofs << in+ in+ "MolarityFactor_Mol = np.where(MolarityFactor_Mol == 1, self.Vol, 1)" << endl;
    ofs << in+ in+ "Count_Mol *= MolarityFactor_Mol" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, Idx_Mol, Count_Mol, axis=1)" << endl;
    ofs << endl;

    ofs << in+ "def GetMolNames(self):" << endl;
    ofs << in+ in+ "return self.Mol_Names" << endl;
    ofs << endl;

    ofs << in+ "def GetDistNames(self):" << endl;
    ofs << in+ in+ "return self.Dist_Names" << endl;
    ofs << endl;

    ofs << in+ "def GetPosNames(self):" << endl;
    ofs << in+ in+ "return self.Pos_Names" << endl;
    ofs << endl;

    ofs << in+ "def ExportLegend(self):" << endl;
    ofs << in+ in+ "return self.Legends" << endl;
    ofs << endl;

    ofs << in+ "def ExportData(self, Time):" << endl;
    ofs << in+ in+ "Data = np.asmatrix(np.zeros(2 + " << MolNames.size() << "))" << endl;
    i = 0;
    int i_SimStep = i + 1;
    ofs << in+ in+ "Data[:, " << i << ":" << i_SimStep << "] = Time" << endl;

    int i_Vol = i_SimStep + 1;
    ofs << in+ in+ "Data[:, " << i_SimStep << ":" << i_Vol << "] = self.Vol" << endl;

    int i_Count_Mol = i_Vol + MolNames.size();

    ofs << in+ in+ "Data[:, " << i_Vol << ":" << i_Count_Mol << "] = self.Count_All";
    if (!ObjLoc.empty()) { ofs << "[0]"; } ofs << endl;

    ofs << in+ in+ "return Data" << endl;
    ofs << endl;

    // class FDataset
    ofs << "class FDataset:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "self.Legend = list()" << endl;
    ofs << in+ in+ "self.Data = None" << endl;
    ofs << endl;

    ofs << in+ "def PrintLegend(self):" << endl;
    ofs << in+ in+ "print(self.Legend)" << endl;
    ofs << endl;

    ofs << in+ "def PrintData(self):" << endl;
    ofs << in+ in+ "print(self.Data)" << endl;
    ofs << endl;

    ofs << in+ "def ExportDataAsList(self):" << endl;
    ofs << in+ in+ "return np.array(self.Data).flatten().tolist()" << endl;
    ofs << endl;

    // class FSimulation
    ofs << "class FSimulation:" << endl;
    ofs << in+ "def __init__(self, InState, InDataset, InDataManager):" << endl;
    ofs << in+ in+ "self.N_SimSteps = 0" << endl;
    ofs << in+ in+ "self.SimStep = 0" << endl;
    ofs << in+ in+ "self.SimTimeResolutionPerSecond = 0" << endl;
    ofs << endl;
    ofs << in+ in+ "self.State = InState" << endl;
    ofs << in+ in+ "self.Dataset = InDataset" << endl;
    ofs << in+ in+ "self.DataManager = InDataManager" << endl;
    ofs << endl;


    // Distribution to Count
    // Update spatial distribution to a specific coord values --> to LocationList method

    if (!MolLoc.empty()) {
        for (auto& molLoc : MolLoc) {
            ofs << in+ in+ "self.Idx_DistToCoord_" << molLoc->Name << " = None" << endl;
        }
    }

    if (!MolLoc.empty()) {
        for (auto& molLoc : MolLoc) {
            ofs << in+ in+ "self.Idx_CoordToDist_" << molLoc->Name << " = None" << endl;
        }
    }

    std::vector<FCount *> Counts_Restore = Context.GetSubList_CountList("Restore");
    for (auto& Count : Counts_Restore) {
        ofs << in+ in+ "self.Idx_Restore_" << Count->Name << " = None" << endl;
    }

    std::vector<FCount *> Counts_Event = Context.GetSubList_CountList("Event");
    for (auto& Count : Counts_Event) {
        ofs << in+ in+ "self.Idx_Event_" << Count->Name << " = None" << endl;
    }

    ofs << in+ in+ "# Debugging" << endl;
    ofs << in+ in+ "self.Debug_Idx_Molecules = list()" << endl;
    ofs << in+ in+ "self.UnitTxt = ''" << endl;
    ofs << in+ in+ "self.Unit = 0" << endl;
    ofs << endl;

    ofs << in+ in+ "# Debugging_Spatial" << endl;
    ofs << in+ in+ "self.Debug_Idx_Pos = list()" << endl;
    ofs << in+ in+ "self.Debug_Idx_Dist = list()" << endl;
    ofs << endl;

    ofs << in+ "def Initialize(self, InN_SimSteps=1000, InTimeResolution=100):" << endl;
    ofs << in+ in+ "print('Simulation Initialized...')" << endl;
    ofs << in+ in+ "self.N_SimSteps = np.array([InN_SimSteps])" << endl;
    ofs << in+ in+ "self.SimTimeResolutionPerSecond = InTimeResolution" << endl;
    ofs << endl;

    // Distribution To Count
    if (!MolLoc.empty()) {
        for (i = 0; i < MolLoc.size(); i++) {
            int Idx = Context.GetIdxByName_MoleculeList(MolLoc[i]->Name);
            ofs << in+ in+ "self.Idx_DistToCoord_" << MolLoc[i]->Name << " = np.array([" << std::to_string(Idx) << "])" << endl;
        }
    }

    if (!MolLoc.empty()) {
        for (i = 0; i < MolLoc.size(); i++) {
            ofs << in+ in+ "self.Idx_CoordToDist_" << MolLoc[i]->Name << " = np.array([" << i << "])" << endl;
        }
    }

    // Restore --> to CountList method
    for (auto& count : Counts_Restore) {
        int Idx = Context.GetIdxByName_MoleculeList(count->Name);
        ofs << in+ in+ "self.Idx_Restore_" << count->Name << " = np.array([" << std::to_string(Idx) << "], ndmin=2)" << endl;
    }
    ofs << endl;

    // Event --> to CountList method
    for (auto& count : Counts_Event) {
        int Idx = Context.GetIdxByName_MoleculeList(count->Name);
        ofs << in+ in+ "self.Idx_Event_" << count->Name << " = np.array([" << std::to_string(Idx) << "], ndmin=2)" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.State.Initialize()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Legend Export" << endl;
    ofs << in+ in+ "self.Dataset.Legend = self.State.ExportLegend()" << endl;
    ofs << in+ in+ "self.DataManager.SetLegend(self.Dataset.Legend)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Data Export" << endl;
    ofs << in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Debugging" << endl;
    ofs << in+ in+ "self.Debug_SetIdxMoleculesToTrack()" << endl;
    ofs << in+ in+ "self.Debug_SetIdxDistAndPosToTrack()" << endl;
    ofs << in+ in+ "self.Debug_SetUnit(Unit)" << endl;
    ofs << endl;

    // regular simloop
    ofs << in+ "def SimLoop_WithSpatialSimulation(self):" << endl;
    bool bDebug_SimFlow = false;
//    bDebug_SimFlow = true;

    ofs << in+ in+ "self.IncrementSimStep()" << endl;

    if (!Option.bDebug) { ofs << "# ";}
    ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
    if (!Option.bDebug) { ofs << "# "; }
    ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
    if (!Option.bDebug) { ofs << "# "; }
    ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Run Reactions" << endl;
    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after NonSpatialSimulation')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }

//    ofs << in+ in+ "# Restore Substrate Count for Sustained Substrate InTransporter" << endl;
//    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
//    ofs << endl;
//
//    if (bDebug_SimFlow) {
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "print('@ after RestoreMoleculeCount')" << endl;
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
//        ofs << endl;
//    }

    ofs << in+ in+ "# Update Substrate Count" << endl;
    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after UpdateCounts')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }
//    ofs << in+ in+ "# temporary: Update Threshold" << endl;
//    ofs << in+ in+ "self.UpdateThreshold()" << endl;
//    ofs << endl;


    ofs << in+ in+ "# Update Spatially Distributed Molecules from Count (special treatment on 'qL' for now by addition)" << endl;
    ofs << in+ in+ "self.CountToDistribution()" << endl;
    ofs << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after CountToDistribution')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }

    ofs << in+ in+ "# Run Spatial Simulation" << endl;
    ofs << in+ in+ "self.SpatialSimulation()" << endl; // TODO: Take delta, instead of updating directly then move up before UpdateCounts
    ofs << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after SpatialSimulation')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }

//    ofs << in+ in+ "# Restore Substrate Count for Sustained Substrate InTransporter" << endl;
//    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
//    ofs << endl;
//
//    if (bDebug_SimFlow) {
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "print('@ after RestoreMoleculeCount')" << endl;
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
//        ofs << endl;
//    }

    ofs << in+ in+ "# Update Spatially Distributed Molecules On Count (special treatment on 'L' and 'qL' for now)" << endl;
    ofs << in+ in+ "self.DistributionToCount()" << endl;
    ofs << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after DistributionToCount')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }

    ofs << in+ "def SimLoop_WithoutSpatialSimulation(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

    if (!Option.bDebug) { ofs << "# ";}
    ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
    if (!Option.bDebug) { ofs << "# "; }
    ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
    if (!Option.bDebug) { ofs << "# "; }
    ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
    ofs << endl;

    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;

    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << endl;

    // simloop for receptivity
    ofs << in+ "def SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

    if (!Option.bDebug) { ofs << "# ";}
    ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
    if (!Option.bDebug) { ofs << "# ";}
    ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
    ofs << endl;

    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after NonSpatialSimulation')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }

//    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
//
//    if (bDebug_SimFlow) {
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "print('@ after RestoreMoleculeCount')" << endl;
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
//        ofs << endl;
//    }

    ofs << in+ in+ "self.UpdateCounts()" << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after UpdateCounts')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }

//    ofs << in+ in+ "self.CountToDistribution()" << endl;
//
//    if (bDebug_SimFlow) {
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "print('@ after CountToDistribution')" << endl;
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
//        if (!Option.bDebug) { ofs << "# "; }
//        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
//        ofs << endl;
//    }

    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after RestoreMoleculeCount')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }

    ofs << in+ in+ "self.DistributionToCount()" << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after DistributionToCount')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }
    ofs << endl;

    ofs << in+ "def Run(self, Spatial=0):" << endl;
    ofs << in+ in+ "if Spatial:" << endl;
    ofs << in+ in+ in+ "self.Run_WithSpatialSimulation()" << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "self.Run_WithoutSpatialSimulation()" << endl;
    ofs << endl;

    ofs << in+ "def Run_WithSpatialSimulation(self):" << endl;
    ofs << in+ in+ "print('Simulation Run Begins...')" << endl;
    ofs << endl;

    ofs << in+ in+ "while self.SimStep < self.N_SimSteps:" << endl;
    ofs << in+ in+ in+ "self.SimLoop_WithSpatialSimulation()" << endl;

    ofs << in+ in+ in+ "# Trigger Event on Substrate Count" << endl;
    ofs << in+ in+ in+ "self.TriggerEventMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ "def Run_WithoutSpatialSimulation(self):" << endl;
    ofs << in+ in+ "print('Simulation Run_WithoutSpatialSimulation Begins...')" << endl;
    ofs << endl;

    ofs << in+ in+ "while self.SimStep < self.N_SimSteps:" << endl;
    ofs << in+ in+ in+ "self.SimLoop_WithoutSpatialSimulation()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Trigger Event on Substrate Count" << endl;
    ofs << in+ in+ in+ "self.TriggerEventMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ in+ "print('Simulation Run_WithoutSpatialSimulation Completed')" << endl;
    ofs << endl;

    // TODO: scale up exporting
    ofs << in+ "def ExportData(self):" << endl;
//    if (Context.LocationList.empty()) {
        ofs << in+ in+ "self.Dataset.Data = self.State.ExportData(self.SimStep/self.SimTimeResolutionPerSecond)" << endl;
        ofs << in+ in+ "self.DataManager.Add(self.Dataset.Data)" << endl;
//    } else { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "def SaveState(self, FileName):" << endl;
    ofs << in+ in+ "self.DataManager.SaveToNpyFile(self.GetCount_All(), FileName)" << endl;
    ofs << endl;

    ofs << in+ "def ApplySimTimeResolution(self, Rate):" << endl;
    ofs << in+ in+ "return Rate / self.SimTimeResolutionPerSecond" << endl;
    ofs << endl;

    // Distribution to Count

    // get coord for 'e'
    // get 'L' count from distribution by 'e' coord
    // get idx for 'L'
    // Update count at idx


    // TODO: make it conditional for specially treated molecules in the system
    ofs << in+ "def DistributionToCount(self):" << endl;
    bool PassSwitch = true;
    if (!MolLoc.empty() & !ObjLoc.empty()) {
        for (auto& molLoc : MolLoc) {
            std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
            for (auto& UniqueName : ObjUniqueNames) {
                if ((molLoc->Name == "L") || (molLoc->Name == "qL")) {
                    ofs << in+ in+ "Count = self.GetCountFromDistributionByNamesOfDistAndPos('" << molLoc->Name << "', " << "'" << UniqueName << "')" << endl;
                    ofs << in+ in+ "self.State.Count_All[:, self.Idx_DistToCoord_" << molLoc->Name
                        << "] = Count.reshape(-1, 1)" << endl;
                    PassSwitch = false;
                }
            }
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "def CountToDistribution(self):" << endl;
    PassSwitch = true;
    if (!MolLoc.empty() & !ObjLoc.empty()) {
        for (auto& molLoc : MolLoc) {
            // TODO: Hardcoding to add qL count to distribution
            if (molLoc->Name == "qL") {
                std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
                for (auto& UniqueName : ObjUniqueNames) {
                    ofs << in+ in+ "Count = self.State.Count_All[:, self.Idx_DistToCoord_" << molLoc->Name << "].reshape(-1)" << endl;
                    ofs << in+ in+ "self.AddToDist(self.Idx_CoordToDist_" << molLoc->Name << ", self.State.Pos_X, self.State.Pos_Y, Count)" << endl;
                    PassSwitch = false;
                }
            }
        }
    }
    if (PassSwitch) {
        ofs << in+ in+ "pass" << endl; // here
    }
    ofs << endl;

    // Restore
    ofs << in+ "def RestoreMoleculeCount(self):" << endl;
    if (!Counts_Restore.empty()){
        for (auto& count : Counts_Restore) {
            std::string Amount;
            float MolarityFactor = Context.GetMolarityFactorByName_CountList(count->Name);

            if (MolarityFactor) { Amount = Utils::SciFloat2Str(count->Amount) + " * self.State.Vol"; }
            else                { Amount = Utils::SciFloat2Str(count->Amount); }

            // TODO: Special treatment for 'qL' as if qL[:] = 0;
            if (count->Name == "qL") {
                Amount = "0";
            }

            ofs << in + in + "np.put_along_axis(self.State.Count_All, self.Idx_Restore_" << count->Name
                << ".reshape(1, -1), " << Amount << ", axis=1)" << endl;
        }
    } else { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    // Event
    ofs << in+ "def TriggerEventMoleculeCount(self):" << endl;
    if (!Counts_Event.empty()){
        ofs << in+ in+ "Time = self.SimStep / self.SimTimeResolutionPerSecond" << endl;
        for (auto& count : Counts_Event) {
            std::string Amount;
            float MolarityFactor = Context.GetMolarityFactorByName_CountList(count->Name);

            if (MolarityFactor) { Amount = Utils::SciFloat2Str(count->Amount) + " * self.State.Vol"; }
            else                { Amount = Utils::SciFloat2Str(count->Amount); }

            ofs << in+ in+ "if (Time >= " << Utils::SciFloat2Str(count->Begin) << ") & (Time < " << Utils::SciFloat2Str(count->End) << "):" << endl;
            ofs << in+ in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Event_" << count->Name << ", " << Amount << ", axis=1)" << endl;

        }
    } else { ofs << in+ in+ "pass" << endl; }
    ofs << endl;


    if (!Context.ThresholdList.empty()) {
        // temporary code
        ofs << in+ "def GetNewThresholdValues(self):" << endl;
        float ThresholdFactor = 0.999;
        ofs << in+ in+ "return self.GetCount_All()[:, self.State.Idx_Count_Threshold].transpose() * " << Utils::SciFloat2Str(ThresholdFactor) << endl;
        ofs << endl;

        ofs << in+ "def SetThreshold(self, bSteadyState):" << endl;
        ofs << in+ in+ "NewThreshold = self.GetNewThresholdValues()" << endl;
        ofs << in+ in+ "self.State.Pos_Threshold = np.where(np.logical_and(self.State.Pos_Threshold == 0, bSteadyState), NewThreshold, self.State.Pos_Threshold)" << endl;
        ofs << endl;

        ofs << in+ "def PrintThreshold(self, MoleculeNames_Str):" << endl;
        ofs << in+ in+ "self.Debug_PrintSimStepTime(end='\\n')" << endl;
        ofs << in+ in+ "print('Thresholds: ' + MoleculeNames_Str + '(' + self.UnitTxt + ')')" << endl;
        ofs << in+ in+ "print(self.Debug_ApplyUnit(self.State.Pos_Threshold.transpose()))" << endl;
        ofs << endl;

//        ofs << in+ "def UpdateThreshold(self):" << endl;
//        float HomeostasisFactor = 1e-7;
//        ofs << in+ in+ "Count_Now = self.GetCount_All()[:, self.State.Idx_Count_Threshold].transpose()" << endl;
//        ofs << in+ in+ "bSteadyState = abs(Count_Now - self.Count_Prev) / Count_Now < " << Utils::SciFloat2Str(HomeostasisFactor) << endl;
//        ofs << in+ in+ "Idx_SteadyState = np.where(bSteadyState)" << endl;
//        ofs << in+ in+ "self.State.Pos_Threshold[Idx_SteadyState] = Count_Now[Idx_SteadyState]" << endl;
//    //    ofs << in+ in+ "if np.any(Idx_SteadyState):" << endl;
//    //    ofs << in+ in+ in+ "print('[Homeostasis] Steady state has been updated.')" << endl;
//        ofs << in+ in+ "self.State.Count_Prev = Count_Now" << endl;
//        ofs << endl;
    }

    if (!Context.LocationList.empty()) {
        ofs << in+ "# Spatial Simulation related routines" << endl;
        ofs << in+ "def SpatialSimulation(self):" << endl;
//        if (!TransporterReactionTypes.empty()) {
//            ofs << in+ in+ "self.TransporterReactions()" << endl;
//        }
//
//            ofs << endl;
//            ofs << in+ in+ "#self.Debug_PrintSimStepTime()" << endl;
//            ofs << in+ in+ "#print()" << endl;
//            ofs << in+ in+ "#self.Debug_PrintCountsAndDistributions() # Temporary placement. Update with implementing delta for spatial simulation" << endl;
//            ofs << endl;
        if (!Context.MotilityList.empty()) {
            ofs << in+ in+ "self.SpatialLocation()" << endl;
        }

        // diffuse either before or after (nonspatial simulation + spatial location simulation)
        ofs << in+ in+ "self.SpatialDiffusion()" << endl;
        ofs << endl;
    }

    ofs << in+ "def SpatialDiffusion(self):" << endl;
    PassSwitch = true;
    if (!MolLoc.empty()) {
        int N_Dist = Context.GetNames_LocationList("Molecule").size();
        // TODO: update to 3d array
        for (i = 0; i < N_Dist; i++) {


//            if (MolLoc[i]->Name == "L") { continue; }
//            else {
                ofs << in+ in+ "self.State.Dist_All[" << i << "] = SimF.DiffuseDistribution_4Cell(self.State.Dist_All[" << i << "])" << endl;
//                PassSwitch = false;
//            }


            // Future implementation
            // ofs << in+ in+ "self.State.Dist_All[" << i << "] = SimF.DiffuseDistribution_8Cell(self.State.Dist_All[" << i << "])" << endl;
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    // Print StandardReaction for each Reaction Type
    ofs << in+ "def TransporterReactions(self):" << endl;
    PassSwitch = true;
    for (auto& Type : TransporterReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            ofs << in+ in+ "self.TransporterReaction_" << Type << "()" << endl;
            PassSwitch = false;
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    // TransporterReaction function
    for (auto& Type : TransporterReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            bool bMolaritySys = Context.CheckMolarityFactorTrueForAny_CountList();

            std::string AmountTextStr;
            if (bMolaritySys) { AmountTextStr = "Conc_"; }
            else              { AmountTextStr = "Count_"; }

            ofs << in+ "def TransporterReaction_" << Type << "(self):" << endl;

//            if (Option.bDebug) {
                ofs << endl;
                ofs << in+ in+ "# Debugging tool" << endl;
                ofs << in+ in+ "# self.Debug_PrintNames(PUT_IDX_HERE)" << endl;
                ofs << endl;
//            }

            std::string Idx_Cargo = "self.State.Idx_Cargo_" + Type;
            std::string Idx_Dist_Cargo = "self.State.Idx_Dist_Cargo_" + Type;
            std::string Idx_Transporter = "self.State.Idx_Transporter_" + Type;

            std::string Amount_Cargo_Inside = AmountTextStr + "Cargo_Inside";
            std::string Amount_Cargo_Outside = AmountTextStr + "Cargo_Outside";
            std::string Amount_Transporter = AmountTextStr + "Transporter";

            ofs << in+ in+ Amount_Transporter << " = ";
            std::string Line = "self.State.Count_All[:, " + Idx_Transporter + "]";
            if (bMolaritySys) { Line = "SimF.CountToConc(" + Line + ", self.State.Vol)"; }
            ofs << Line << endl;

            // Get Concentrations
            ofs << in+ in+ Amount_Cargo_Inside << " = ";
            Line = "self.State.Count_All[:, " + Idx_Cargo + "]";
            if (bMolaritySys) { Line = "SimF.CountToConc(" + Line + ", self.State.Vol)"; }
            ofs << Line << endl;

            ofs << in+ in+ Amount_Cargo_Outside << " = ";
            Line = "self.GetCountFromDistribution(" + Idx_Dist_Cargo + ", self.State.Pos_X, self.State.Pos_Y).transpose()";
            if (bMolaritySys) { Line = "SimF.CountToConc(" + Line + ", self.State.Vol)"; }
            ofs << Line << endl;

//            // Regulators
//            if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
//                ofs << in+ in+ AmountTextStr << "Regulator = ";
//                if (bMolaritySys) {
//                    ofs << "SimF.CountToConc(self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "], self.State.Vol)" << endl;
//                } else {
//                    ofs << "self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "]" << endl;
//                }
//            }

            // Calculate Rate
            ofs << in+ in+ "Rate = SimF.Eqn_" << Type << "(" <<
            Amount_Transporter << ", " <<
            Amount_Cargo_Inside << ", " <<
            Amount_Cargo_Outside << ", " <<
            "self.State.Const_ki_" << Type << ", " <<
            "self.State.Const_ko_" << Type << ")" << endl;

            //            // Regulators
//            if (Typing[substrate] != StandardReactionTypes[0]) {
//                ofs << AmountTextStr << "Regulator, ";
//            }
//            // Reaction constants
//            ofs << "self.State.Const_k_" << substrate << "_" << Type;
//            if (Typing[substrate] != StandardReactionTypes[0]) {
//                ofs << ", self.State.Const_K_" << Type << ", ";
//                ofs << "self.State.Const_n_" << Type;
//            }
//            ofs << ")" << endl;

            // Apply Time Resolution
            ofs << in+ in+ "Rate = self.ApplySimTimeResolution(Rate)" << endl;

            if (bMolaritySys) {
                // Apply stoichiometry
                ofs << in+ in+ "dConc_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
                // Convert to counts
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
            } else {
                // Apply stoichiometry
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
            }
            // Apply delta counts for molecules in the stoichiometry matrix
            ofs << in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
            ofs << in+ in+ "self.AddToDist(" << Idx_Dist_Cargo << ", self.State.Pos_X, self.State.Pos_Y, -dCount_Mol_InStoichMatrix.transpose())" << endl;
            ofs << endl;
        }
    }

    ofs << in+ "def SpatialLocation(self):" << endl;
    ofs << in+ in+ "Evaluations = self.EvaluateChemotaxisThreshold()" << endl;
    ofs << in+ in+ "self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle = SimF.BacterialChemotaxis(Evaluations, self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle)" << endl;
    ofs << in+ in+ "self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle = SimF.CorrectOutOfBounds(self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle, self.State.Dimension_X, self.State.Dimension_Y)" << endl;
    ofs << endl;

    ofs << in+ "def NonSpatialSimulation(self):" << endl;
    PassSwitch = true;
    if (!StandardReactionTypes.empty()) {
        for (auto& Type: StandardReactionTypes) {
            std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
            if (!ReactionSubList.empty()) {
                ofs << in+ in+ "self.StandardReactions()" << endl;
                ofs << endl;
                PassSwitch = false;
                break;
            }
        }
    }

    if (!EnzReactionTypes.empty()) {
        for (auto& Type: EnzReactionTypes) {
            std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
            if (!ReactionSubList.empty()) {
                ofs << in+ in+ "self.EnzymaticReactions()" << endl;
                ofs << endl;
                PassSwitch = false;
                break;
            }
        }
    }

    // TODO: encapsulate this part with each polymerase to allow more process-specific customization
    if (!PolymeraseList.empty()){
        ofs << in+ in+ "self.Polymerase_InitiationReactions()" << endl;
        ofs << in+ in+ "self.Polymerase_ElongationReactions()" << endl;
        ofs << in+ in+ "self.Polymerase_TerminationReactions()" << endl;
        PassSwitch = false;
    }

    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "# Biochemical Reaction related routines" << endl;
    // StandardReaction function
    for (auto& Type : StandardReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            // Typing set up
            std::map<std::string, std::string> Typing;
            std::vector<std::string> SubstrateTypes { "Reactant", "Product" };
            Typing.emplace("Reactant", Type);
            Typing.emplace("Product", StandardReactionTypes[0]); // default type

            bool bMolaritySys = Context.CheckMolarityFactorTrueForAny_CountList();
            std::string AmountTextStr;
            if (bMolaritySys) { AmountTextStr = "Conc_"; }
            else              { AmountTextStr = "Count_"; }

            ofs << in+ "def StandardReaction_" << Type << "(self):" << endl;

//            if (Option.bDebug) {
                ofs << endl;
                ofs << in+ in+ "# Debugging tool" << endl;
                ofs << in+ in+ "# self.Debug_PrintNames(PUT_IDX_HERE)" << endl;
                ofs << endl;
//            }

            // Get Concentrations
            // Reactants and products
            for (auto& substrate : SubstrateTypes) {
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << in+ in+ AmountTextStr << substrate << "_" << i << " = ";
                    if (bMolaritySys) {
                        ofs << "SimF.CountToConc(self.State.Count_All[:, self.State.Idx_" << substrate << "_" << i << "_" << Type << "], self.State.Vol)" << endl;
                    } else {
                        ofs << "self.State.Count_All[:, self.State.Idx_" << substrate << "_" << i << "_" << Type << "]" << endl;
                    }
                }
            }

            // Regulators
            if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                ofs << in+ in+ AmountTextStr << "Regulator = ";
                if (bMolaritySys) {
                    ofs << "SimF.CountToConc(self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "], self.State.Vol)" << endl;
                } else {
                    ofs << "self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "]" << endl;
                }
            }

            // Calculate Rate
            for (auto& substrate : SubstrateTypes) {
                ofs << in+ in+ "Rate_" << substrate << " = SimF.Eqn_" << Typing[substrate] << "_" << N_MoleculesAllowed << "(";

                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << AmountTextStr << substrate << "_" << i << ", ";
                }

                // Regulators
                if (Typing[substrate] != StandardReactionTypes[0]) {
                    ofs << AmountTextStr << "Regulator, ";
                }
                // Reaction constants
                ofs << "self.State.Const_k_" << substrate << "_" << Type;
                if (Typing[substrate] != StandardReactionTypes[0]) {
                    ofs << ", self.State.Const_K_" << Type << ", ";
                    ofs << "self.State.Const_n_" << Type;
                }
                ofs << ")" << endl;

                // Apply Time Resolution
                ofs << in+ in+ "Rate_" << substrate << " = self.ApplySimTimeResolution(Rate_" << substrate << ")" << endl;

                // compare conc
                ofs << in+ in+ "Rate_" << substrate << " = SimF.CheckRateAndConc_" << N_MoleculesAllowed << "(Rate_" << substrate;
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << ", " << AmountTextStr << substrate << "_" << i;
                }
                ofs << ")" << endl;
            }

            // Tally Rates
            ofs << in+ in+ "Rate = Rate_Reactant - Rate_Product" << endl;

            if (bMolaritySys) {
                // Apply stoichiometry
                ofs << in+ in+ "dConc_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
                // Convert to counts
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
            } else {
                // Apply stoichiometry
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
            }
            // Apply delta counts for molecules in the stoichiometry matrix
            ofs << in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
            ofs << endl;
        }
    }

    // EnzReaction function // TODO: Add bMolaritySys option as above
    for (auto& Type : EnzReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {

            // Typing set up
            std::map<std::string, std::string> Typing;
            std::vector<std::string> SubstrateTypes { "Reactant", "Product" };
            Typing.emplace("Reactant", Type);
            Typing.emplace("Product", EnzReactionTypes[0]); // default type

            ofs << in+ "def EnzymaticReaction_" << Type << "(self):" << endl;

//            if (Option.bDebug) {
                ofs << endl;
                ofs << in+ in+ "# Debugging tool" << endl;
                ofs << in+ in+ "# self.Debug_PrintNames(PUT_IDX_HERE)" << endl;
                ofs << endl;
//            }

            // Get Concentrations
            // Enzyme
            ofs << in+ in+ "Conc_Enz = SimF.CountToConc(self.State.Count_All[:, self.State.Idx_Enz_" << Type << "], self.State.Vol)" << endl;
            // Reactants and products or EnzSubstrate
            if (Type.find("Enz_Standard") != std::string::npos) {
                for (auto& substrate : SubstrateTypes) {
                    for (int i = 0; i < N_MoleculesAllowed; i++) {
                        ofs << in+ in+ "Conc_" << substrate << "_" << i << " = SimF.CountToConc(self.State.Count_All[:, self.State.Idx_" << substrate << "_" << i << "_" << Type << "], self.State.Vol)" << endl;
                    }
                }
            } else if (Type.find("Enz_MichaelisMenten") != std::string::npos) {
                ofs << in+ in+ "Conc_EnzSub = SimF.CountToConc(self.State.Count_All[:, self.State.Idx_EnzSub_" << Type << "], self.State.Vol)" << endl;
            } else {
                Utils::Assertion(false, "Unsupported Enz Reaction Type: " + Type);
            }

            // Regulators
            if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                ofs << in+ in+ "Conc_Regulator = SimF.CountToConc(self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "], self.State.Vol)" << endl;
            }

            // Calculate Rate
            if (Type.find("Enz_Standard") != std::string::npos) {
                for (auto& substrate : SubstrateTypes) {
                    ofs << in+ in+ "Rate_" << substrate << " = SimF.Eqn_" << Typing[substrate] << "_" << N_MoleculesAllowed << "(Conc_Enz, ";

                    for (int i = 0; i < N_MoleculesAllowed; i++) {
                        ofs << "Conc_" << substrate << "_" << i << ", ";
                    }

                    // Regulators
                    if (Typing[substrate] != EnzReactionTypes[0]) {
                        ofs << "Conc_Regulator, ";
                    }
                    // Reaction constants
                    ofs << "self.State.Const_k_" << substrate << "_" << Type;
                    if (Typing[substrate] != EnzReactionTypes[0]) {
                        ofs << ", self.State.Const_K_" << Type << ", ";
                        ofs << "self.State.Const_n_" << Type;
                    }
                    ofs << ")" << endl;

                    // Apply Time Resolution
                    ofs << in+ in+ "Rate_" << substrate << " = self.ApplySimTimeResolution(Rate_" << substrate << ")" << endl;

                    // compare conc
                    ofs << in+ in+ "Rate_" << substrate << " = SimF.CheckRateAndConc_" << N_MoleculesAllowed << "(Rate_" << substrate;
                    for (int i = 0; i < N_MoleculesAllowed; i++) {
                        ofs << ", Conc_" << substrate << "_" << i;
                    }
                    ofs << ")" << endl;
                }

                // Tally Rates
                ofs << in+ in+ "Rate = Rate_Reactant - Rate_Product" << endl;

            } else if (Type.find("Enz_MichaelisMenten") != std::string::npos) {
                ofs << in+ in+ "Rate = SimF.Eqn_" << Type << "(Conc_Enz, Conc_EnzSub";

                // Regulators
                if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                    ofs << ", Conc_Regulator";
                }

                // Reaction constants
                ofs << ", self.State.Const_kcat_" << Type << ", self.State.Const_KM_" << Type;

                if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                    ofs << ", self.State.Const_K_" << Type;
                    if (Type.find("Allosteric") != std::string::npos) {
                        ofs << ", self.State.Const_n_" << Type;
                    }
                }
                ofs << ")" << endl;

                // Apply TimeResolution to the rate
                ofs << in+ in+ "Rate = self.ApplySimTimeResolution(Rate)" << endl;
                // compare conc for michaelis menten?
            }
            // Apply stoichiometry
            ofs << in+ in+ "dConc_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
            // Convert to counts
            ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
            // Apply delta counts for molecules in the stoichiometry matrix
            ofs << in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
            ofs << endl;
        }
    }

    // Print StandardReaction for each Reaction Type
    ofs << in+ "def StandardReactions(self):" << endl;
    PassSwitch = true;
    for (auto& Type : StandardReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            ofs << in+ in+ "self.StandardReaction_" << Type << "()" << endl;
            PassSwitch = false;
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;


    // Print EnzymeReaction for each Reaction Type
    ofs << in+ "def EnzymaticReactions(self):" << endl;
    PassSwitch = true;
    for (auto& Type : EnzReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            ofs << in+ in+ "self.EnzymaticReaction_" << Type << "()" << endl;
            PassSwitch = false;
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }  
    ofs << endl;

    if (!PolymeraseList.empty()) {

        ofs << in+ "def Polymerase_InitiationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {

            Polymerase_InitiationReaction(ofs, Polymerase);
        }


        ofs << in+ "def Polymerase_ElongationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {

            Polymerase_ElongationReaction(ofs, Polymerase);
        }

        ofs << in+ "def Polymerase_TerminationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {

            Polymerase_TerminationReaction(ofs, Polymerase);
        }
    }

//    if (EnzymeList.empty() & PolymeraseList.empty()) {
//            ofs << in+ in+ "pass" << endl;
//            ofs << endl;
//    }

    ofs << in+ "# Useful routines" << endl;

    ofs << in+ "def GetSimStep(self):" << endl;
    ofs << in+ in+ "return self.SimStep" << endl;
    ofs << endl;

    ofs << in+ "def GetSimTime(self):" << endl;
    ofs << in+ in+ "return self.SimStep / self.SimTimeResolutionPerSecond" << endl;
    ofs << endl;

    ofs << in+ "def IncrementSimStep(self):" << endl;
    ofs << in+ in+ "self.SimStep += 1" << endl;
    ofs << endl;

    ofs << in+ "def GetCount(self, Idx):" << endl;
    ofs << in+ in+ "return self.State.Count_All[:, Idx]" << endl;
    ofs << endl;

    ofs << in+ "def GetCount_All(self):" << endl;
    ofs << in+ in+ "return np.copy(self.State.Count_All)" << endl;
    ofs << endl;

    ofs << in+ "def GetConcentration(self, Idx):" << endl;
    ofs << in+ in+ "return SimF.CountToConc(self.State.Count_All[:, Idx], self.State.Vol)" << endl;
    ofs << endl;

    ofs << in+ "def GetDistribution(self, Idx):" << endl;
    ofs << in+ in+ "return self.State.Dist_All[Idx]" << endl;
    ofs << endl;

    ofs << in+ "def GetDistWidth(self):" << endl;
    ofs << in+ in+ "return self.State.Dimension_X" << endl;
    ofs << endl;

    ofs << in+ "def GetDistHeight(self):" << endl;
    ofs << in+ in+ "return self.State.Dimension_Y" << endl;
    ofs << endl;

    ofs << in+ "def AddTodCount(self, Idx, Values):" << endl;
    ofs << in+ in+ "dCountToAdd = np.zeros_like(self.State.dCount_All)" << endl;
    ofs << in+ in+ "np.put_along_axis(dCountToAdd, Idx, Values, axis=1)" << endl;
    ofs << in+ in+ "dCount_All_New = self.State.dCount_All + dCountToAdd" << endl;
    ofs << in+ in+ "ZeroTest = dCount_All_New + self.State.Count_All" << endl;
    ofs << in+ in+ "self.State.dCount_All = np.where(ZeroTest < 0, dCount_All_New - ZeroTest, dCount_All_New)" << endl;
    ofs << endl;

    ofs << in+ "def UpdateCounts(self):" << endl;
    ofs << in+ in+ "self.State.Count_All += self.State.dCount_All" << endl;
    ofs << in+ in+ "self.CleardCounts()" << endl;
    ofs << endl;

    ofs << in+ "def CleardCounts(self):" << endl;
    ofs << in+ in+ "self.State.dCount_All = np.zeros_like(self.State.dCount_All)" << endl;
    ofs << endl;

    ofs << in+ "def AddToDist(self, Idx, X, Y, Value):" << endl;
    ofs << in+ in+ "self.State.Dist_All[Idx, X.astype(int), Y.astype(int)] += Values" << endl;

//    ofs << in+ in+ "X, Y = self.Rescale(X, Y)" << endl;
//    ofs << in+ in+ "self.State.Dist_All[Idx] = SimF.BilinearExtrapolation(self.State.Dist_All[Idx], X, Y, Value)" << endl;
//    ofs << in+ in+ "self.State.Dist_All[Idx, X.astype(int), Y.astype(int)] += Value" << endl;
    ofs << endl;

    // TODO: Use SimIdx in the future
    ofs << in+ "# Temporary routines" << endl;

    ofs << in+ "def GetMolIdx(self, Name):" << endl;
    ofs << in+ in+ "return self.State.Mol_Name2Idx[Name]" << endl;
    ofs << endl;

    ofs << in+ "def GetCountByName(self, Name):" << endl;
    ofs << in+ in+ "return self.GetCount(self.GetMolIdx(Name))" << endl;
    ofs << endl;

    ofs << in+ "def GetConcentrationByName(self, Name):" << endl;
    ofs << in+ in+ "return self.GetConcentration(self.GetMolIdx(Name))" << endl;
    ofs << endl;

    ofs << in+ "def GetThreshold(self, Pos_Idx):" << endl;
    ofs << in+ in+ "return self.State.Pos_Threshold[:, Pos_Idx]" << endl;
    ofs << endl;

    ofs << in+ "def GetPosIdx(self, Name):" << endl;
    ofs << in+ in+ "return self.State.Pos_Name2Idx[Name]" << endl;
    ofs << endl;

    ofs << in+ "def GetPositionXY(self, Idx):" << endl;
    ofs << in+ in+ "return np.take(self.State.Pos_X, Idx), np.take(self.State.Pos_Y, Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetPositionXYByName(self, Name):" << endl;
    ofs << in+ in+ "Idx = self.GetPosIdx(Name)" << endl;
    ofs << in+ in+ "return self.GetPositionXY(Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetPositionXYAngle(self, Idx):" << endl;
    ofs << in+ in+ "return np.take(self.State.Pos_X, Idx), np.take(self.State.Pos_Y, Idx), np.take(self.State.Pos_Angle, Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetPositionXYAngleByName(self, Name):" << endl;
    ofs << in+ in+ "Idx = self.GetPosIdx(Name)" << endl;
    ofs << in+ in+ "return self.GetPositionXYAngle(Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetDistIdx(self, Name):" << endl;
    ofs << in+ in+ "return self.State.GetDistNames().index(Name)" << endl;
    ofs << endl;

    ofs << in+ "def GetDistributionByName(self, Name):" << endl;
    ofs << in+ in+ "Idx = self.GetDistIdx(Name)" << endl;
    ofs << in+ in+ "return self.GetDistribution(Idx)" << endl;
    ofs << endl;

    ofs << in+ "def Rescale(self, X, Y):" << endl;
    ofs << in+ in+ "return X / self.State.FoldReduction, Y / self.State.FoldReduction " << endl;
    ofs << endl;


    ofs << in+ "def GetCountFromDistribution(self, Dist_Idx, X, Y):" << endl;
//    ofs << in+ in+ "return SimF.BilinearInterpolation(self.State.Dist_All[Dist_Idx], X, Y)" << endl;
    ofs << in+ in+ "return self.State.Dist_All[Dist_Idx, X.astype(int), Y.astype(int)]" << endl;
    ofs << endl;

    ofs << in+ "def GetCountFromDistributionByNamesOfDistAndPos(self, NameOfDist, NameOfPos):" << endl;
    // temporary code
    ofs << in+ in+ "X, Y = self.GetPositionXYByName(NameOfPos)" << endl;
    ofs << in+ in+ "return self.GetCountFromDistributionByNameOfDistAndXY(NameOfDist, X, Y)" << endl;
    ofs << endl;

    ofs << in+ "def GetCountFromDistributionByNameOfDistAndXY(self, NameOfDist, X, Y):" << endl;
    // temporary code
    ofs << in+ in+ "Dist_Idx = self.GetDistIdx(NameOfDist)" << endl;
//    ofs << in+ in+ "X, Y = self.Rescale(X, Y)" << endl;
    ofs << in+ in+ "return self.GetCountFromDistribution(Dist_Idx, X, Y)" << endl;
    ofs << endl;

    ofs << in+ "def ApplyCountToDistribution(self, Dist, X, Y, Count):" << endl;
    ofs << in+ in+ "Dist[X.astype(int), Y.astype(int)] += Count" << endl;
    ofs << endl;

    ofs << in+ "def ApplyCountToDistributionByNameAndPos(self, NameOfDist, NameOfPos):" << endl;
    ofs << in+ in+ "Count = self.GetCountByName" << endl;
    ofs << endl;


    ofs << in+ "# Temporary routines" << endl;

    ofs << in+ "def OverElongationCorrection(self, Len_Elongated, Max):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ "Len_Over = np.where(Len_Elongated > Max, Len_Elongated - Max, 0)" << endl;
    ofs << in+ in+ "return Len_Elongated - Len_Over" << endl;
    ofs << endl;

    ofs << in+ "def BuildingBlockConsumption(self, Freq, N_Elongated_PerSpecies):" << endl;
    ofs << in+ in+ "Raw = SimF.DetermineAmountOfBuildingBlocks(Freq, N_Elongated_PerSpecies)" << endl;
    ofs << in+ in+ "Rounded = np.around(Raw)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Discrepancy handling" << endl;
    ofs << in+ in+ "N_Elongated = np.sum(N_Elongated_PerSpecies)" << endl;
    ofs << in+ in+ "Discrepancy = np.sum(Rounded) - N_Elongated" << endl;

    ofs << in+ in+ "NUniq_BuildingBlocks = Freq.shape[1]" << endl;
    ofs << in+ in+ "Sets, Remainder = np.divmod(Discrepancy, NUniq_BuildingBlocks)" << endl;
    ofs << in+ in+ "return Rounded + np.ones(NUniq_BuildingBlocks) * np.int32(Sets) + np.concatenate((np.ones(np.int32(np.round(Remainder))), np.zeros(np.int32(np.around(NUniq_BuildingBlocks - Remainder)))))" << endl;
    ofs << endl;

    ofs << in+ "# Polymerase Reaction related" << endl;
    ofs << in+ "def Initiation(self, Len_Template, Len_Target, Idx_Pol, Idx_Template, Idx_TemplateSubset, Weight, PolThreshold):" << endl;
    ofs << in+ in+ "# Get available, active polymerase count - TO BE UPDATED with more regulatory algorithms" << endl;
    ofs << in+ in+ "Count_Pol = self.GetCount(Idx_Pol)" << endl;
    ofs << in+ in+ "Count_Pol_Active = np.floor_divide(Count_Pol, 2).astype(int)" << endl;
    ofs << in+ in+ "Count_Pol_Occupied = np.count_nonzero(np.where(Len_Target != -1, 1, 0)) * PolThreshold" << endl;
    ofs << in+ in+ "Count_Pol_Avail = Count_Pol_Active - Count_Pol_Occupied" << endl;
    ofs << in+ in+ "Count_Pol_FunctionalUnit = np.floor_divide(Count_Pol_Avail, PolThreshold)" << endl;
    ofs << in+ in+ "Count_Pol_Avail = np.where(Count_Pol_FunctionalUnit > 0, Count_Pol_FunctionalUnit, 0)[0, 0]" << endl;
    ofs << endl;

    ofs << in+ in+ "# Get final initiation weight by applying initiation site count" << endl;
    ofs << in+ in+ "Count_Template_Complete = self.GetCount(Idx_Template)" << endl;
    ofs << in+ in+ "Count_Template_Nascent = np.count_nonzero(np.where(Len_Template != -1, 1, 0), axis=0)   # Assumption: each nascent template has one highly efficient initiation site" << endl;
    ofs << in+ in+ "Count_TemplateSubset_Nascent = np.take(Count_Template_Nascent, Idx_TemplateSubset)   # May need to update np.take with fancy indexing [:, idx]" << endl;
    ofs << in+ in+ "Count_InitiationSite = Count_Template_Complete + Count_TemplateSubset_Nascent" << endl;
    ofs << in+ in+ "Weight_Initiation = Count_InitiationSite * Weight " << endl;
    ofs << endl;

    ofs << in+ in+ "# Get randomly selected target indices" << endl;
    ofs << in+ in+ "Idx_Selected = SimF.PickRandomIdx(Count_Pol_Avail, Idx_Template, Weight_Initiation)" << endl;
    ofs << in+ in+ "Len_Target_Initiated = SimF.InsertZeroIntoNegOneElementInLenMatrix(Len_Target, Idx_Selected)" << endl;
    ofs << in+ in+ "# Export Data" << endl;
    ofs << in+ in+ "# N_Initiated" << endl;
    ofs << endl;

    ofs << in+ in+ "return Len_Target_Initiated" << endl;
    ofs << endl;

    ofs << in+ "def Elongation(self, Len, Max, Rate, Weight, Freq, Idx_PolSub, Idx_BB):" << endl;
    ofs << in+ in+ "NUniq_BuildingBlocks = Freq.shape[1]" << endl;
    ofs << in+ in+ "NUniq_Species = Freq.shape[0]" << endl;
    ofs << endl;

//    ofs << in+ in+ "dLength = np.matmul(SMatrix,Rate)
    ofs << in+ in+ "dLength = self.ApplySimTimeResolution(Rate)   # this is not necessarily true based on the reaction input" << endl;
    ofs << in+ in+ "Len_Elongated = np.where(Len >= 0, Len + dLength, Len)" << endl;
    ofs << endl;

    ofs << in+ in+ "Len_Trimmed = self.OverElongationCorrection(Len_Elongated, Max)" << endl;

    ofs << in+ in+ "N_Elongated_PerSpecies = np.asmatrix(np.sum(Len_Trimmed - Len, axis=0))   # This step loses shape for some reason, hence apply matrix again" << endl;
    ofs << in+ in+ "N_Elongated = np.sum(N_Elongated_PerSpecies)" << endl;
    ofs << endl;

    ofs << in+ in+ "Consumed_BB = self.BuildingBlockConsumption(Freq, N_Elongated_PerSpecies)" << endl;

    ofs << in+ in+ "# Update dCount for BuildingBlocks" << endl;
    ofs << in+ in+ "self.AddTodCount(Idx_BB, -Consumed_BB)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Update dCount for Polymerase Reaction Substrates (To be updated by the reaction matrix form" << endl;
    ofs << in+ in+ "self.AddTodCount(Idx_PolSub, N_Elongated)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Export Data" << endl;
    ofs << in+ in+ "# N_Elongated" << endl;

    ofs << in+ in+ "return Len_Trimmed" << endl;
    ofs << endl;

    ofs << in+ "def Termination(self, Len, Max, Idx_Target):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ "Bool_Completed = (Len == Max)" << endl;
    ofs << in+ in+ "N_Completed_PerSpecies = np.sum(Bool_Completed, axis=0)" << endl;
    ofs << in+ in+ "N_Completed = np.sum(N_Completed_PerSpecies)" << endl;
    ofs << in+ in+ "Len_Completed = np.where(Bool_Completed, -1, Len)" << endl;

    ofs << in+ in+ "# Update dCount for BuildingBlocks" << endl;
    ofs << in+ in+ "self.AddTodCount(Idx_Target, N_Completed_PerSpecies)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Export Data" << endl;
    ofs << in+ in+ "# N_Completed" << endl;

    ofs << in+ in+ "return Len_Completed" << endl;
    ofs << endl;

    ofs << in+ "def EvaluateChemotaxisThreshold(self):" << endl;
    ofs << in+ in+ "ThresholdedMolecules = self.GetCount(self.State.Idx_Count_Threshold).transpose()" << endl;
    ofs << in+ in+ "return np.array(ThresholdedMolecules) < self.State.Pos_Threshold" << endl;
    ofs << endl;

    ofs << in+ "# Debugging tools" << endl;
    ofs << in+ "def Debug_PrintNames(self, IdxArray):" << endl;
    ofs << in+ in+ "ListOfNames = list()" << endl;
    ofs << in+ in+ "if IdxArray.dim == 1:" << endl;
    ofs << in+ in+ in+ "IdxArray = IdxArray.tolist()" << endl;
    ofs << in+ in+ "if IdxArray.dim == 2:" << endl;
    ofs << in+ in+ in+ "IdxArray = IdxArray.tolist()[0]" << endl;
    ofs << in+ in+ "for Idx in IdxArray.tolist()[0]:" << endl;
    ofs << in+ in+ in+ "ListOfNames.append(self.State.GetMolNames()[Idx])" << endl;
    ofs << in+ in+ "print(ListOfNames)" << endl;
    ofs << endl;

    ofs << in+ "def Debug_SetIdxMoleculesToTrack(self):" << endl;
    ofs << in+ in+ "# Add a list of molecules to track for debugging every simulation step" << endl;
    ofs << in+ in+ "Debug_Names_Molecules = []" << endl; // TODO: take input from command line
    if (MolLoc.size() == 1) {
        ofs << in+ in   + "# Debug_Names_Molecules = ['L', 'Am', ]" << endl;
    } else if (MolLoc.size() == 2) {
        ofs << in+ in+ "# Debug_Names_Molecules = ['L', 'qL', 'Am', 'qAm', 'qA', 'qAL', 'qAmL']" << endl;
    }
    ofs << in+ in+ "# Debug_Names_Molecules = ['Am', 'AmL', 'L']" << endl;
    ofs << in+ in+ "# Debug_Names_Molecules = ['Am', 'AmL', 'L', 'qL', 'pc_qL']" << endl;
    ofs << endl;
    ofs << in+ in+ "if Debug_Names_Molecules == []:" << endl;
    ofs << in+ in+ in+ "Debug_Names_Molecules = self.State.GetMolNames()" << endl;
    ofs << endl;
    ofs << in+ in+ "for Name in Debug_Names_Molecules:" << endl;
    ofs << in+ in+ in+ "self.Debug_Idx_Molecules.append(self.State.GetMolNames().index(Name))" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCounts(self, Switch):" << endl;
    ofs << in+ in+ "for Idx_Pos in self.Debug_Idx_Pos:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintCountsForSinglePosition(Switch, Idx_Pos)" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCountsForSinglePosition(self, Switch, Idx_Pos=0):" << endl;
    ofs << in+ in+ "for Idx_Mol in self.Debug_Idx_Molecules:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintCount(Idx_Mol, Idx_Pos)" << endl;
    ofs << in+ in+ "print()" << endl;
    ofs << in+ in+ "if Switch:" << endl;
    ofs << in+ in+ in+ "for Idx_Mol in self.Debug_Idx_Molecules:" << endl;
    ofs << in+ in+ in+ in+ "self.Debug_PrintdCount(Idx_Mol, Idx_Pos)" << endl;
    ofs << in+ in+ in+ "print()" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintSimStepTime(self, end='\\t| '):" << endl;
    ofs << in+ in+ "Time = self.GetSimTime()" << endl;
    ofs << in+ in+ "print(self.SimStep, '(', round(Time, 3), 's)', end=end)" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCount(self, Idx_Mol, Idx_Pos=0):" << endl;
    ofs << in+ in+ "print(' ' + self.State.GetMolNames()[Idx_Mol], end=': ')" << endl;
    ofs << in+ in+ "print('{:010e}'.format(self.Debug_ApplyUnit(self.State.Count_All[Idx_Pos][Idx_Mol])), self.UnitTxt, end=' | ')" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintdCount(self, Idx_Mol, Idx_Pos=0):" << endl;
    ofs << in+ in+ "print('d' + self.State.GetMolNames()[Idx_Mol], end=': ')" << endl;
    ofs << in+ in+ "print('{:010e}'.format(self.Debug_ApplyUnit(self.State.dCount_All[Idx_Pos][Idx_Mol])), self.UnitTxt, end=' | ')" << endl;
    ofs << endl;

    ofs << in+ "def Debug_SetUnit(self, Input):" << endl;
    ofs << in+ in+ "self.UnitTxt = Input" << endl;
    ofs << in+ in+ "if Unit == 'nM':" << endl;
    ofs << in+ in+ in+ "self.Unit = 1e-9 * " << Numbers::GetAvogadroStr() << endl;
    ofs << in+ in+ "elif Unit == 'uM':" << endl;
    ofs << in+ in+ in+ "self.Unit = 1e-6 * " << Numbers::GetAvogadroStr() << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "self.Unit = 1" << endl;
    ofs << endl;

    ofs << in+ "def Debug_ApplyUnit(self, Value):" << endl;
    ofs << in+ in+ "return Value / self.Unit" << endl;
    ofs << endl;

    // Debugging for spatial simulation
    ofs << in+ "def Debug_SetIdxDistAndPosToTrack(self):" << endl;
    if (!Context.LocationList.empty()) {
        ofs << in+ in+ "# Add a list of distributions at and around selected positions to track for debugging every simulation step" << endl;
        ofs << in+ in+ "Debug_Names_Positions = []" << endl; // TODO: take input from command line
        ofs << in+ in+ "if Debug_Names_Positions == []:" << endl;
        ofs << in+ in+ in+ "self.Debug_Idx_Pos = list(range(len(self.State.GetPosNames())))" << endl;
        ofs << in+ in+ "else:   # does not properly cover all cases due to redundant names allowed through for loops" << endl;
        ofs << in+ in+ in+ "for Name in Debug_Names_Positions:" << endl;
        ofs << in+ in+ in+ in+ "self.Debug_Idx_Pos.append(self.State.GetPosNames().index(Name))" << endl;
        ofs << endl;
        ofs << in+ in+ "Debug_Names_Distributions = []" << endl; // TODO: take input from command line
        ofs << in+ in+ "if Debug_Names_Distributions == []:" << endl;
        ofs << in+ in+ in+ "Debug_Names_Distributions = self.State.GetDistNames()" << endl;
        ofs << in+ in+ "for Name in Debug_Names_Distributions:" << endl;
        ofs << in+ in+ in+ "self.Debug_Idx_Dist.append(self.State.GetDistNames().index(Name))" << endl;
    } else { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "def Debug_PrintDistributions(self):" << endl;
    ofs << in+ in+ "for Idx_Pos in self.Debug_Idx_Pos:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintDistributionsForSinglePosition(Idx_Pos)" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCountsAndDistributions(self):" << endl;
    ofs << in+ in+ "for Idx_Pos in self.Debug_Idx_Pos:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintDistributionsForSinglePosition(Idx_Pos)" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintCountsForSinglePosition(DisplayCount, Idx_Pos)" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintDistributionsForSinglePosition(self, Idx_Pos):" << endl;
    ofs << in+ in+ "X = self.State.Pos_X[Idx_Pos].astype(int)" << endl;
    ofs << in+ in+ "Y = self.State.Pos_Y[Idx_Pos].astype(int)" << endl;
    ofs << in+ in+ "print('Position #' + str(Idx_Pos) + ' Name: ' + self.State.GetPosNames()[Idx_Pos] + ' @ ({}, {})'.format(X, Y))" << endl;
    ofs << in+ in+ "DistNames = ''" << endl;
    ofs << in+ in+ "DistValues = dict()" << endl;
    ofs << in+ in+ "DistValues[0] = ''" << endl;
    ofs << in+ in+ "DistValues[1] = ''" << endl;
    ofs << in+ in+ "DistValues[2] = ''" << endl;
    ofs << in+ in+ "for Idx_Dist in self.Debug_Idx_Dist:" << endl;
    ofs << in+ in+ in+ "DistNames, DistValues = self.Debug_PrintDistribution(DistNames, DistValues, Idx_Dist, X, Y)" << endl;
    ofs << in+ in+ "print(DistNames, '\\n', DistValues[0], '\\n', DistValues[1], '\\n', DistValues[2])" << endl;
    ofs << endl;

    // When no reduction and interpolation
    ofs << in+ "def Debug_PrintDistribution(self, DistNames, DistValues, Idx_Dist, X, Y):" << endl;
    ofs << in+ in+ "DistNames += self.State.GetDistNames()[Idx_Dist] + ' (' + self.UnitTxt + ')' + '\\t\\t\\t\\t\\t\\t\\t\\t'" << endl;
    ofs << in+ in+ "DistValues[0] += self.Debug_GetArrayInString(self.Debug_ApplyUnit(self.State.Dist_All[Idx_Dist][X-1:X+2, Y-1])) + '\\t'" << endl;
    ofs << in+ in+ "DistValues[1] += self.Debug_GetArrayInString(self.Debug_ApplyUnit(self.State.Dist_All[Idx_Dist][X-1:X+2, Y]))   + '\\t'" << endl;
    ofs << in+ in+ "DistValues[2] += self.Debug_GetArrayInString(self.Debug_ApplyUnit(self.State.Dist_All[Idx_Dist][X-1:X+2, Y+1])) + '\\t'" << endl;
    ofs << in+ in+ "return DistNames, DistValues" << endl;
    ofs << endl;

    ofs << in+ "def Debug_GetArrayInString(self, Array_1D):" << endl;
    ofs << in+ in+ "Array_1D_Str = ''" << endl;
    ofs << in+ in+ "for Idx in range(Array_1D.shape[0]):" << endl;
    ofs << in+ in+ in+ "Array_1D_Str += SimF.SciFloat(Array_1D[Idx]) + ' '" << endl;
    ofs << in+ in+ "return Array_1D_Str" << endl;
    ofs << endl;



    // class FDataManager
    ofs << "class FDataManager:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "# SimOut" << endl;
    ofs << in+ in+ "self.Legend = list()" << endl;
    ofs << in+ in+ "self.DataBuffer = list()" << endl;
    ofs << endl;

    ofs << in+ "def SetLegend(self, InLegend):" << endl;
    ofs << in+ in+ "self.Legend = InLegend" << endl;
    ofs << endl;

    ofs << in+ "def Add(self, InData):" << endl;
    ofs << in+ in+ "self.DataBuffer.append(InData)" << endl;
    ofs << endl;

    ofs << in+ "def SaveToTsvFile(self, Legend, Data, InFileName):" << endl;
    ofs << in+ in+ "with open(InFileName, 'w', newline='', encoding='utf-8') as OutFile:" << endl;
    ofs << in+ in+ in+ "TsvWriter = csv.writer(OutFile, delimiter='\\t')" << endl;
    ofs << in+ in+ in+ "if Legend:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter.writerow(Legend)" << endl;
    ofs << in+ in+ in+ "for Row in Data:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter.writerow(np.array(Row).flatten().tolist())" << endl;
    ofs << endl;

    ofs << in+ "def SaveToNpyFile(self, Data, InFileName):" << endl;
    ofs << in+ in+ "np.save('%s.npy' % InFileName, Data)" << endl;
    ofs << endl;

    ofs << in+ "def SaveCountToFile(self, InFileName):" << endl;
    ofs << in+ in+ "self.SaveToTsvFile(self.Legend, self.DataBuffer, InFileName)" << endl;
    ofs << endl;

    if (Context.LocationList.empty()) {
        // BODY
        ofs << "def main():   # add verbose" << endl;
        ofs << endl;

        // Instantiate Objects
        ofs << in+ "State = FState()" << endl;
        ofs << in+ "Data = FDataset()" << endl;
        ofs << in+ "DataManager = FDataManager()" << endl;
        ofs << in+ "Simulation = FSimulation(State, Data, DataManager)" << endl;
        ofs << endl;

        // Simulation
        ofs << in+ "Simulation.Initialize(N_SimSteps, SimStepTimeResolution)" << endl;
        ofs << in+ "Simulation.Run(Spatial=0) # 0: WithoutSpatialSimulation, 1: WithSpatialSimulation" << endl;
        ofs << endl;
        ofs << in+ "DataManager.SaveCountToFile('" << Option.SimResultFile.c_str() << "')" << endl;
        ofs << endl;
    }
        // MAIN
//        ofs << "if __name__ == '__main__':" << endl;
//        ofs << in+ "parser = ArgumentParser()" << endl;
//        ofs << in+ "parser.add_argument('--save-fig', dest='save_fname', type=str, help='Save figure to file')" << endl;
//        ofs << in+ "args = parser.parse_args()" << endl;
//        ofs << in+ "main()" << endl;
//        ofs << endl;

    std::cout << "  Simulation module has been generated: ";
}
