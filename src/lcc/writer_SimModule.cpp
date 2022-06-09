#include "writer.h"
#include <algorithm>

using namespace std;

void FWriter::SimModule(int Sim_Steps, int Sim_Resolution, int Map_Width, int Map_Height)
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

    if (Option.bDebug) {
        ofs << "import plot" << endl;
    }

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
    auto Compartments = Context.GetSubList_ContainerList("Compartment");

    int MatrixSize = 1;
    if (!Compartments.empty()) {
        MatrixSize = Context.GetCounts_ContainerList("Compartment");
    }
    ofs << in+ in+ "self.Count_All = np.zeros([" << MatrixSize << ", " << Context.MoleculeList.size() << "])" << endl;
    ofs << in+ in+ "self.dCount_All = np.zeros([" << MatrixSize << ", " << Context.MoleculeList.size() << "])" << endl;
    ofs << endl;

    if (!Context.LocationList.empty()) {
        Initialize_SpatialSimulation(ofs, Map_Width, Map_Height);
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
    std::vector<FMolecule *> PolymeraseList = Context.GetSubList_MoleculeList("Polymerase");
    std::vector<FMolecule *> DNAPs = Context.GetSubList_MoleculeList("DNAP");
    std::vector<FMolecule *> RNAPs = Context.GetSubList_MoleculeList("RNAP");
    std::vector<FMolecule *> Ribosomes = Context.GetSubList_MoleculeList("Ribosome");

    if (!PolymeraseList.empty()) {
        if      (!DNAPs.empty())        {}
        else if (!RNAPs.empty())        {}
        else if (!Ribosomes.empty())    {}
        else                            { Utils::Assertion(false, "ERROR: Polymerase List is empty."); }
    }

    std::vector<std::vector<FMolecule *>> PolymeraseTypes = {DNAPs, RNAPs, Ribosomes};
    Initialize_PolymeraseReaction_Matrix(ofs, PolymeraseTypes);

    if (!DNAPs.empty())             { Initialize_PolymeraseReaction_DNAP(ofs, DNAPs); }
    else if (!RNAPs.empty())        { Initialize_PolymeraseReaction_RNAP(ofs, RNAPs); }
    else if (!Ribosomes.empty())    { Initialize_PolymeraseReaction_Ribosome(ofs, Ribosomes); }

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

    // Temporary database from tsv
    if (Context.CheckForEcoli()) {
        ofs << in+ in+ "DatabaseFileName = r'./Database/genes.tsv'" << endl;
        ofs << in+ in+ "Database = self.OpenTSVDatabase(DatabaseFileName)" << endl;
        ofs << endl;
    }

    // for legends
    std::vector<std::string> MolNames = Context.GetNameListByType_MoleculeList("All");
    ofs << in+ in+ "self.Mol_Names = [" << Utils::JoinStr2Str(MolNames) << "]" << endl;
    ofs << in+ in+ "self.Mol_Name2Idx = {";
    for (int i = 0; i < MolNames.size(); i++) {
        ofs << "'" << MolNames[i] << "' : np.array([" << i << "]), ";
    } ofs << "}" << endl;
    ofs << in+ in+ "self.Legends = ['SimStep', 'Vol', " << Utils::JoinStr2Str(MolNames) << "]" << endl;
    ofs << endl;

    // Initialize all molecule counts
    std::vector<int> Idx_Mol = Context.GetIdxListByType_MoleculeList("All");
    std::vector<float> InitialCount_Molecules;
    std::vector<float> MolarityFactor_Molecules;
    std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");

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
        ofs << Utils::JoinInt2Str_Idx(Idx_Mol) << "], ndmin=2)" << endl;
    }
    else {
        for (auto& ObjName : ObjUniqueNames) {
            int Count = int(Context.GetInitialCountByName_CountList(ObjName));
            for (int i = 0; i < (Count); i++) {
                ofs << "[" << Utils::JoinInt2Str_Idx(Idx_Mol) << "], ";
            }
        }
        ofs << "])" << endl;
    }

    ofs << in+ in+ "Count_Mol = np.array([";
    if (ObjLoc.empty()) {
        ofs << Utils::JoinFloat2Str(InitialCount_Molecules);
    }
    else {
        for (auto& ObjName : ObjUniqueNames) {
            int Count = int(Context.GetInitialCountByName_CountList(ObjName));
            for (int i = 0; i < (Count); i++) {
                ofs << "[" << Utils::JoinFloat2Str(InitialCount_Molecules) << "], ";
            }
        }
    }
    ofs << "])" << endl;

    ofs << in+ in+ "MolarityFactor_Mol = np.array([" << Utils::JoinFloat2Str(MolarityFactor_Molecules) << "])" << endl;
    ofs << in+ in+ "MolarityFactor_Mol = np.where(MolarityFactor_Mol == 1, self.Vol, 1)" << endl;
    ofs << in+ in+ "Count_Mol *= MolarityFactor_Mol" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, Idx_Mol, Count_Mol, axis=1)" << endl;
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

    SetUp_PolymeraseReaction_Matrix(ofs, PolymeraseTypes);

    if (!DNAPs.empty())             { SetUp_PolymeraseReaction_DNAP(ofs, DNAPs); }
    else if (!RNAPs.empty())        { SetUp_PolymeraseReaction_RNAP(ofs, RNAPs); }
    else if (!Ribosomes.empty())    { SetUp_PolymeraseReaction_Ribosome(ofs, Ribosomes); }

//
//        //    std::vector<std::string> PolymeraseNames = Context.GetNames_PolymeraseList(PolymeraseList);
//        std::vector<FPolymeraseReaction *> PolymeraseReactionList = Context.GetList_Polymerase_ReactionList();
//
//        std::vector<std::vector<int>> StoichMatrix_PolymeraseReaction = Context.GetStoichiometryMatrix_PolymeraseReaction(PolymeraseReactionList);
//
//
//
//        SetUp_PolymeraseReaction(ofs, Polymerase, Rate, FreqBBFileName, MaxLenFileName, Idx_Pol, Idx_Template, Idx_TemplateSubset, Idx_Target, Idx_PolSub, Idx_PolBB, Threshold);
//    }

    // Print SetUp_TransporterReaction for each Reaction Type
    for (auto& Type : TransporterReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            SetUp_TransporterReaction(ofs, Type, ReactionSubList);
        }
    }


    ofs << in+ "def GetMolNames(self):" << endl;
    ofs << in+ in+ "return self.Mol_Names" << endl;
    ofs << endl;

    ofs << in+ "def GetDistNames(self):" << endl;
    if (!MolLoc.empty() || !ObjLoc.empty())  { ofs << in+ in+ "return self.Dist_Names" << endl; }
    else                                    { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "def GetPosNames(self):" << endl;
    if (!MolLoc.empty() || !ObjLoc.empty())  { ofs << in+ in+ "return self.Pos_Names" << endl; }
    else                                    { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "def ExportLegend(self):" << endl;
    ofs << in+ in+ "return self.Legends" << endl;
    ofs << endl;

    ofs << in+ "def ExportData(self, Time):" << endl;
    ofs << in+ in+ "Data = np.asmatrix(np.zeros(2 + " << MolNames.size() << "))" << endl;
    int Initial = 0;
    int i_SimStep = Initial + 1;
    ofs << in+ in+ "Data[:, " << Initial << ":" << i_SimStep << "] = Time" << endl;

    int i_Vol = i_SimStep + 1;
    ofs << in+ in+ "Data[:, " << i_SimStep << ":" << i_Vol << "] = self.Vol" << endl;

    int i_Count_Mol = i_Vol + MolNames.size();

    ofs << in+ in+ "Data[:, " << i_Vol << ":" << i_Count_Mol << "] = self.Count_All";
    if (!Compartments.empty()) { ofs << "[0]"; } ofs << endl;

    ofs << in+ in+ "return Data" << endl;
    ofs << endl;

    ofs << in+ "# Temporary database routines" << endl;
    ofs << endl;

    ofs << in+ "def LoadFASTADatabase(self, db_fname):" << endl;
    ofs << in+ in+ "db = list()" << endl;
    ofs << in+ in+ "with open(db_fname) as fp:" << endl;
    ofs << in+ in+ in+ "sequences = fp.read().split('>')[1:]" << endl;
    ofs << in+ in+ in+ "for i, sequence in enumerate(sequences):" << endl;
    ofs << in+ in+ in+ in+ "list_of_rows = sequence.split('\\n')" << endl;
    ofs << in+ in+ in+ in+ "db.append(''.join(list_of_rows[1:]))" << endl;
    ofs << in+ in+ "return db" << endl;
    ofs << endl;

    ofs << in+ "def OpenFASTADatabase(self, db_fname):" << endl;
    ofs << in+ in+ "db = self.LoadFASTADatabase(db_fname)" << endl;
    ofs << in+ in+ "return db" << endl;
    ofs << endl;

    ofs << in+ "def LoadTSVDatabase(self, db_fname):" << endl;
    ofs << in+ in+ "db = None" << endl;
    ofs << in+ in+ "with open(db_fname) as fp:" << endl;
    ofs << in+ in+ in+ "csv_reader = csv.reader(fp, delimiter='\\t')" << endl;
    ofs << in+ in+ in+ "list_of_rows = list(csv_reader)" << endl;
    ofs << in+ in+ in+ "db = list_of_rows[1:]" << endl;
    ofs << in+ in+ "return db" << endl;
    ofs << endl;

    ofs << in+ "def OpenTSVDatabase(self, db_fname):" << endl;
    ofs << in+ in+ "db = self.LoadTSVDatabase(db_fname)" << endl;
    ofs << in+ in+ "Database_Gene = self.ParseGenes(db)" << endl;
    ofs << in+ in+ "return Database_Gene" << endl;
    ofs << endl;

    ofs << in+ "def ParseGenes(self, db_genes):" << endl;
    ofs << in+ in+ "db = dict()" << endl;
    ofs << in+ in+ "NUniq_Genes = len(db_genes)" << endl;
    ofs << in+ in+ "db['Symbol'] = list()" << endl;
    ofs << in+ in+ "db['Length'] = np.zeros(NUniq_Genes)" << endl;
    ofs << in+ in+ "db['Coord'] = np.zeros(NUniq_Genes)" << endl;
    ofs << in+ in+ "db['Dir'] = np.zeros(NUniq_Genes)" << endl;
    ofs << in+ in+ "db['Seq'] = list()" << endl;
    ofs << endl;
    ofs << in+ in+ "Dir = dict()" << endl;
    ofs << in+ in+ "Dir['+'] = 1" << endl;
    ofs << in+ in+ "Dir['-'] = -1" << endl;
    ofs << endl;
    ofs << in+ in+ "for i, Value in enumerate(db_genes):" << endl;
    ofs << in+ in+ in+ "Length, Name, Seq, RNAID, Coordinate, Direction, Symbol, Type, GeneID, MonomerID = Value" << endl;
    ofs << in+ in+ in+ "db['Symbol'].append(Symbol)" << endl;
    ofs << in+ in+ in+ "db['Length'][i] = (len(Seq))" << endl;
    ofs << in+ in+ in+ "db['Coord'][i] = int(Coordinate)" << endl;
    ofs << in+ in+ in+ "db['Dir'][i] = Dir[Direction]" << endl;
    ofs << in+ in+ in+ "db['Seq'].append(Seq)" << endl;
    ofs << endl;
    ofs << in+ in+ "return db" << endl;
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
        for  (int i = 0; i < MolLoc.size(); i++) {
            int Idx = Context.GetIdxByName_MoleculeList(MolLoc[i]->Name);
            ofs << in+ in+ "self.Idx_DistToCoord_" << MolLoc[i]->Name << " = np.array([" << std::to_string(Idx) << "])" << endl;
        }
    }

    if (!MolLoc.empty()) {
        for  (int i = 0; i < MolLoc.size(); i++) {
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

    ofs << in+ in+ "# Restore Substrate Count for Sustained Substrate InTransporter" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << endl;

    if (bDebug_SimFlow) {
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "print('@ after RestoreMoleculeCount')" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        if (!Option.bDebug) { ofs << "# "; }
        ofs << in+ in+ "self.Debug_PrintDistributions()" << endl;
        ofs << endl;
    }

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


    // temporary for qL efflux
    if (!Context.ReactionList.empty()) {
        for (auto& reaction : Context.ReactionList) {
            if (reaction->CheckIfProduct("qL")) {
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
                break;
            }
        }
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

    ofs << in+ in+ "self.CellDivision()" << endl;
    ofs << endl;

    ofs << in+ "def SimLoop_WithoutSpatialSimulation(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

    if (!Option.bDebug) { ofs << "# ";}
    ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
    if (!Option.bDebug) { ofs << "# "; }
    ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
    ofs << endl;

    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;

    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ in+ "self.CellDivision()" << endl;
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

    ofs << in+ in+ "self.CellDivision()" << endl;
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

            ofs << in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Restore_" << count->Name
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
        for  (int i = 0; i < N_Dist; i++) {

            if (MolLoc[i]->Name == "L") {
                // skip glucose diffusion
                continue;
            } else if (MolLoc[i]->Name == "qL") {
                // fast diffusion applied for qL
                ofs << in+ in+ "self.State.Dist_All[" << i << "] = SimF.DiffuseDistribution_FAST(self.State.Dist_All[" << i << "])" << endl;
                PassSwitch = false;
            } else {
                ofs << in+ in+ "self.State.Dist_All[" << i << "] = SimF.DiffuseDistribution_4Cell(self.State.Dist_All[" << i << "])" << endl;
                PassSwitch = false;
            }

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

    if (!Context.MotilityList.empty()) {
        ofs << in+ "def SpatialLocation(self):" << endl;
        ofs << in+ in+ "Evaluations = self.EvaluateChemotaxisThreshold()" << endl;
        ofs << in+ in+ "self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle = SimF.BacterialChemotaxis(Evaluations, self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle)" << endl;
        ofs << in+ in+ "self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle = SimF.CorrectOutOfBounds(self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle, self.State.Dimension_X, self.State.Dimension_Y)" << endl;
        ofs << endl;
    }

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
        ofs << in+ in+ "self.PolymeraseReactions()" << endl;
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
        ofs << in+ "def PolymeraseReactions(self):" << endl;
        if      (!DNAPs.empty())        { ofs << in+ in+ "self.Replication()" << endl; }
        else if (!RNAPs.empty())        { ofs << in+ in+ "self.Transcription()" << endl; }
        else if (!Ribosomes.empty())    { ofs << in+ in+ "self.Translation()" << endl; }
        ofs << endl;

//        ofs << in+ in+ "self.Polymerase_InitiationReactions()" << endl;
//        ofs << in+ in+ "self.Polymerase_ElongationReactions()" << endl;
//        ofs << in+ in+ "self.Polymerase_TerminationReactions()" << endl;
//

        if (!DNAPs.empty()) {
            ofs << in+ "def Replication(self):" << endl;
            Polymerase_InitiationReaction_DNAP(ofs, DNAPs);
            Polymerase_ElongationReaction_DNAP(ofs, DNAPs);
            Polymerase_TerminationReaction_DNAP(ofs, DNAPs);
            ofs << endl;

        } else if (!RNAPs.empty()) {
            ofs << in+ "def Transcription(self):" << endl;
            Polymerase_InitiationReaction_RNAP(ofs, RNAPs);
            Polymerase_ElongationReaction_RNAP(ofs, RNAPs);
            Polymerase_TerminationReaction_RNAP(ofs, RNAPs);
            ofs << endl;

        } else if (!Ribosomes.empty()) {
            ofs << in+ "def Translation(self):" << endl;
            Polymerase_InitiationReaction_Ribosome(ofs, Ribosomes);
            Polymerase_ElongationReaction_Ribosome(ofs, Ribosomes);
            Polymerase_TerminationReaction_Ribosome(ofs, Ribosomes);
            ofs << endl;
        }
    }

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
    ofs << in+ in+ "self.State.Dist_All[Idx, X.astype(int), Y.astype(int)] += Value" << endl;

//    ofs << in+ in+ "X, Y = self.Rescale(X, Y)" << endl;
//    ofs << in+ in+ "self.State.Dist_All[Idx] = SimF.BilinearExtrapolation(self.State.Dist_All[Idx], X, Y, Value)" << endl;
//    ofs << in+ in+ "self.State.Dist_All[Idx, X.astype(int), Y.astype(int)] += Value" << endl;
    ofs << endl;

    ofs << in+ "# External Simulation Control routines" << endl;
    ofs << endl;

    ofs << in+ "def Receptivity(self, N_SimulationsToPass=50):" << endl;
    ofs << in+ in+ "for _ in range(N_SimulationsToPass):" << endl;
    ofs << in+ in+ in+ "self.SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution()" << endl;
    ofs << in+ in+ in+ "self.ExportData()" << endl;
    ofs << in+ in+ "self.SimLoop_WithSpatialSimulation()" << endl;
    ofs << in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    // TODO: Use SimIdx in the future
    ofs << in+ "# Temporary routines" << endl;
    ofs << endl;

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

    ofs << in+ "def GetReplicationCompletionRate(self, Name):" << endl;
    if (!DNAPs.empty()) {
        ofs << in+ in+ "Idx = self.GetPosIdx(Name)" << endl;
        ofs << in+ in+ "ReplicationCompletion = self.State.Len_NascentChromosome / self.State.MaxLen_NascentChromosome" << endl;
        ofs << in+ in+ "ReplicationCompletion = np.where(ReplicationCompletion < 0, 0, ReplicationCompletion)" << endl;
        ofs << in+ in+ "return ReplicationCompletion[Idx]" << endl;
    } else {
        ofs << in+ in+ "return 0" << endl;
    }
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

    ofs << in+ "# Central dogma routines" << endl;
    ofs << endl;

    ofs << in+ "def OverElongationCorrection(self, Len_Elongated, Max):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ "Len_Over = np.where(Len_Elongated > Max, Len_Elongated - Max, 0)" << endl;
    ofs << in+ in+ "return Len_Elongated - Len_Over" << endl;
    ofs << endl;

    ofs << in+ "def OverElongationCorrection_Transcription(self, Pos_Pol_Elongated, Max, Dir_Pol):" << endl;
    ofs << in+ in+ "Delta = Pos_Pol_Elongated * Dir_Pol - Max * Dir_Pol" << endl;
    ofs << in+ in+ "Pos_Pol_OverElongated = np.where(Delta > 0, Delta, 0)" << endl;
    ofs << in+ in+ "return Pos_Pol_Elongated - (Pos_Pol_OverElongated * Dir_Pol)" << endl;
    ofs << endl;

    ofs << in+ "def BuildingBlockConsumption(self, Freq, N_Elongated_PerSpecies):" << endl;
    ofs << in+ in+ "Raw = np.array(SimF.DetermineAmountOfBuildingBlocks(Freq, N_Elongated_PerSpecies))" << endl;
    ofs << in+ in+ "Rounded = np.around(Raw)" << endl;
    ofs << in+ in+ "return Rounded" << endl;
    ofs << endl;

    // TODO: Discrepancy handling needs to be reimplemented
//    ofs << in+ in+ "# Discrepancy handling" << endl;
//    ofs << in+ in+ "N_Elongated = np.sum(N_Elongated_PerSpecies, axis=1)" << endl;
//    ofs << in+ in+ "Discrepancy = np.sum(Rounded, axis=1) - N_Elongated" << endl;
//
//    ofs << in+ in+ "NUniq_BuildingBlocks = Freq.shape[1]" << endl;
//    ofs << in+ in+ "Sets, Remainder = np.divmod(Discrepancy, NUniq_BuildingBlocks)" << endl;
//    ofs << in+ in+ "Sets = np.array(np.asmatrix(Sets).transpose() * np.asmatrix(np.ones(NUniq_BuildingBlocks))).astype(int)" << endl;
//    ofs << in+ in+ "Remainder = np.concatenate((np.ones(Remainder.astype(int)), np.zeros((NUniq_BuildingBlocks - Remainder).astype(int))))" << endl;
//    ofs << in+ in+ "return Rounded + Sets + Remainder" << endl;
//    ofs << endl;

    ofs << in+ "# Polymerase Reaction related" << endl;

    if (!DNAPs.empty()) {
        ofs << in+ "def Replication_GetAvailablePolymerases(self, Len_Target, Idx_Pol, PolThreshold):" << endl;
        ofs << in+ in+ "# Get available, active polymerase count - TO BE UPDATED with more regulatory algorithms" << endl;
        ofs << in+ in+ "Count_Pol = self.GetCount(Idx_Pol)" << endl;
        ofs << in+ in+ "Count_Pol_Active = np.floor_divide(Count_Pol, 2).astype(int)" << endl;
        ofs << in+ in+ "Count_Pol_Occupied = np.where(Len_Target != -1, 1, 0) * PolThreshold" << endl;
        ofs << in+ in+ "Count_Pol_Avail = Count_Pol_Active - Count_Pol_Occupied" << endl;
        ofs << in+ in+ "Count_Pol_FunctionalUnit = np.floor_divide(Count_Pol_Avail, PolThreshold)" << endl;
        ofs << in+ in+ "Count_Pol_Avail = np.where(Count_Pol_FunctionalUnit > 0, Count_Pol_FunctionalUnit, 0)" << endl;
        ofs << in+ in+ "return Count_Pol_Avail" << endl;
        ofs << endl;

        ofs << in+ "def Replication_GetInitiationSites(self, Idx_Template, Len_Template):" << endl;
        ofs << in+ in+ "# Get final initiation weight by applying initiation site count" << endl;
        ofs << in+ in+ "Count_Template_Complete = self.GetCount(Idx_Template)" << endl;
        ofs << in+ in+ "Count_Template_Nascent = np.where(Len_Template != -1, 1, 0)" << endl;
        ofs << in+ in+ "Count_InitiationSite = Count_Template_Complete + Count_Template_Nascent" << endl;
        ofs << in+ in+ "return Count_InitiationSite" << endl;
        ofs << endl;

        // TODO: Update with more mechanistic algorithm later
        ofs << in+ "def Replication_Initiation(self, Len_Template, Len_Target, Idx_Pol, Idx_Template, Idx_TemplateSubset, Weight, PolThreshold):" << endl;
        ofs << in+ in+ "Count_Pol_Avail = self.Replication_GetAvailablePolymerases(Len_Target, Idx_Pol, PolThreshold)" << endl;
        ofs << in+ in+ "Count_InitiationSite = self.Replication_GetInitiationSites(Idx_Template, Len_Template)" << endl;

        // Only compatible with a single nascent chromosomes
        ofs << in+ in+ "Count_ToBind = np.where(np.logical_and(Count_Pol_Avail > 0, Count_InitiationSite > 0), np.minimum(Count_Pol_Avail, Count_InitiationSite), 0)" << endl;
        ofs << in+ in+ "Len_Target_Initiated = Len_Target + Count_ToBind" << endl;
        ofs << in+ in+ "return Len_Target_Initiated" << endl;
        ofs << endl;
        
        ofs << in+ "def Replication_Elongation(self, Len, Max, Rate, Freq, Idx_PolSub, Idx_BB):" << endl;
        ofs << in+ in+ "NUniq_BuildingBlocks = Freq.shape[1]" << endl;
        ofs << in+ in+ "NUniq_Species = Freq.shape[0]" << endl;
        ofs << endl;

    //    ofs << in+ in+ "dLength = np.matmul(SMatrix,Rate)
        ofs << in+ in+ "dLength = self.ApplySimTimeResolution(Rate)   # this is not necessarily true based on the reaction input" << endl;
        ofs << in+ in+ "Len_Elongated = np.where(Len >= 0, Len + dLength, Len)" << endl;
        ofs << in+ in+ "Len_Trimmed = self.OverElongationCorrection(Len_Elongated, Max)" << endl;
        ofs << in+ in+ "N_Elongated = np.array(np.sum(Len_Trimmed - Len, axis=1), ndmin=2).transpose()" << endl;
        ofs << endl;

        ofs << in+ in+ "Consumed_BB = self.BuildingBlockConsumption(Freq, N_Elongated)" << endl;
        ofs << in+ in+ "# Update dCount for BuildingBlocks" << endl;
        ofs << in+ in+ "self.AddTodCount(Idx_BB, -Consumed_BB)" << endl;
        ofs << endl;

        // TODO: Update this with matrix calculation form
        ofs << in+ in+ "# Update dCount for Polymerase Reaction Substrates" << endl;
        ofs << in+ in+ "self.AddTodCount(Idx_PolSub, N_Elongated)" << endl;
        ofs << endl;
        ofs << in+ in+ "return Len_Trimmed" << endl;
        ofs << endl;

        ofs << in+ "def Replication_Termination(self, Len, Max, Idx_Target, Idx_Pol, PolThreshold):   # Some polymerization process may not have max" << endl;
        ofs << in+ in+ "Bool_Completed = (Len == Max)" << endl;
        ofs << in+ in+ "N_Completed = np.array(np.sum(Bool_Completed, axis=1), ndmin=2).transpose()" << endl;
        ofs << in+ in+ "Len_Completed = np.where(Bool_Completed, -1, Len)" << endl;

        ofs << in+ in+ "self.AddTodCount(Idx_Target, N_Completed)" << endl;
        ofs << in+ in+ "self.AddTodCount(np.array(Idx_Pol, ndmin=2), N_Completed * PolThreshold)" << endl;
        ofs << endl;

        ofs << in+ in+ "# Export Data" << endl;
        ofs << in+ in+ "# N_Completed" << endl;

        ofs << in+ in+ "return Len_Completed" << endl;
        ofs << endl;

    }

    if (!RNAPs.empty()) {
        ofs << in+ "def Transcription_GetAvailablePolymerases(self, Count_NascentTarget, Idx_Pol, PolThreshold):" << endl;
        ofs << in+ in+ "# Get available, active polymerase count - TO BE UPDATED with more regulatory algorithms" << endl;
        ofs << in+ in+ "Count_Pol = self.GetCount(Idx_Pol)" << endl;
        ofs << in+ in+ "Count_Reg = self.GetCount(Idx_Pol)" << endl;  // TODO: Implement Sigma factors. Temporary same count as pol.
        ofs << in+ in+ "Count_Pol_Reg_Complex = np.where(Count_Pol > Count_Reg, Count_Reg, Count_Pol)" << endl;
        ofs << in+ in+ "Count_Pol_Active = np.floor_divide(Count_Pol_Reg_Complex, 2).astype(int)" << endl;
        ofs << in+ in+ "Count_Pol_Occupied = np.reshape(np.sum(Count_NascentTarget, axis=1), [-1, 1]) * PolThreshold" << endl;
        ofs << in+ in+ "Count_Pol_Avail = Count_Pol_Active - Count_Pol_Occupied" << endl;
        ofs << in+ in+ "Count_Pol_FunctionalUnit = np.floor_divide(Count_Pol_Avail, PolThreshold)" << endl;
        ofs << in+ in+ "Count_Pol_Avail = np.where(Count_Pol_FunctionalUnit > 0, Count_Pol_FunctionalUnit, 0)" << endl;
        ofs << in+ in+ "return Count_Pol_Avail" << endl;
        ofs << endl;

        ofs << in+ "def Transcription_GetInitiationSites(self, Idx_Template, Count_Template_Nascent):" << endl;
        ofs << in+ in+ "# Get final initiation weight by applying initiation site count" << endl;
        ofs << in+ in+ "Count_Template_Complete = self.GetCount(Idx_Template)" << endl;
        ofs << in+ in+ "Count_InitiationSitePerSite = Count_Template_Complete + Count_Template_Nascent" << endl;
        ofs << in+ in+ "return Count_InitiationSitePerSite" << endl;
        ofs << endl;

        // TODO: Update with more mechanistic algorithm later
        ofs << in+ "def Transcription_Initiation(self, Pos_Pol, Pos_Pol_End, Pos_Pol_Template, Dir_Pol, Pos_Template_Start, Pos_Template_End, Dir_Template, Count_NascentTemplate, Count_NascentTarget, Idx_Pol, Idx_Template, Weight, PolThreshold):" << endl;
        ofs << in+ in+ "# Get Available Polymerase complex count" << endl;
        ofs << in+ in+ "Count_Pol_Avail = self.Transcription_GetAvailablePolymerases(Count_NascentTarget, Idx_Pol, PolThreshold)" << endl;
        ofs << endl;

        ofs << in+ in+ "# Get Count and Weight of transcription initiation sites" << endl;
        ofs << in+ in+ "Count_InitiationSites_PerSite = self.Transcription_GetInitiationSites(Idx_Template, Count_NascentTemplate)" << endl;
        ofs << in+ in+ "Count_InitiationSites_Total = np.reshape(np.sum(Count_InitiationSites_PerSite, axis=1), [-1, 1])" << endl;
        ofs << in+ in+ "Weight_InitiationSite = Count_InitiationSites_PerSite * Weight" << endl;
        ofs << in+ in+ "Count_ToBind = np.logical_and(Count_Pol_Avail > 0, Count_InitiationSites_Total > 0)" << endl;
        ofs << in+ in+ "Count_ToBindThisStep = np.where(Count_ToBind, 1, 0)   # Only allowing one binding at a time to use random function in matrix operation" << endl;
        ofs << endl;

        // TODO: Update to work on all rows
        // TODO: Incorporate Weighted system that works with arrays
        ofs << in+ in+ "# Get New polymerase binding sites." << endl;
        ofs << in+ in+ "Idx_MolToBind = np.reshape(np.random.randint(0, high=Weight_InitiationSite.shape[1], size=Weight_InitiationSite.shape[0]), [-1, 1])" << endl;
        ofs << endl;

        ofs << in+ in+ "# Apply New Binding to Count_NascentTarget" << endl;
        ofs << in+ in+ "Count_NascentTarget_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Count_NascentTarget, Idx_MolToBind, Count_ToBindThisStep) " << endl;
        ofs << endl;

        ofs << in+ in+ "# Apply New Binding at start position to Pos_Pol" << endl;
        ofs << in+ in+ "Idx_EmptyElement = np.reshape(SimF.GetIdxOfEmptyElement(Pos_Pol), [-1, 1]) " << endl;
        ofs << endl;

        if (Option.bDebug) {
            ofs << in+ in+ "print('Count_ToBind\\n', Count_ToBind.T)" << endl;
            ofs << in+ in+ "print('Count_ToBindThisStep\\n', Count_ToBindThisStep.T)" << endl;
            ofs << in+ in+ "print('Idx_MolToBind\\n', Idx_MolToBind)" << endl;
            ofs << in+ in+ "print('Idx_EmptyElement\\n', Idx_EmptyElement.T)" << endl;
            ofs << in+ in+ "print('Count_NascentTarget_Initiated\\n', Count_NascentTarget_Initiated)" << endl;
            ofs << endl;
        }

        ofs << in+ in+ "Pos_Template_Start_Initiated = Pos_Template_Start[Idx_MolToBind]" << endl;
        ofs << in+ in+ "Pos_Pol_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Pos_Pol, Idx_EmptyElement, Pos_Template_Start_Initiated * Count_ToBindThisStep + Count_ToBindThisStep) # Additional Count_ToBindThisStep to offset -1 value in the unbound element convention" << endl;
        ofs << endl;
        ofs << in+ in+ "Pos_Template_End_Initiated = Pos_Template_End[Idx_MolToBind]" << endl;
        ofs << in+ in+ "Pos_Pol_End_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Pos_Pol_End, Idx_EmptyElement, Pos_Template_End_Initiated * Count_ToBindThisStep) " << endl;
        ofs << endl;
        ofs << in+ in+ "Pos_Pol_Template_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Pos_Pol_Template, Idx_EmptyElement, Idx_MolToBind * Count_ToBindThisStep + Count_ToBindThisStep) # Additional Count_ToBindThisStep to offset -1 value in the unbound element convention" << endl;
        ofs << endl;
        ofs << in+ in+ "Dir_Template_Initiated = Dir_Template[Idx_MolToBind]" << endl;
        ofs << in+ in+ "Dir_Pol_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Dir_Pol, Idx_EmptyElement, Dir_Template_Initiated * Count_ToBindThisStep) " << endl;
        ofs << endl;

        if (Option.bDebug) {
            ofs << in+ in+ "print('Pos_Pol_Initiated\\n', Pos_Pol_Initiated)" << endl;
            ofs << in+ in+ "print('Pos_Pol_End_Initiated\\n', Pos_Pol_End_Initiated)" << endl;
            ofs << in+ in+ "print('Pos_Pol_Template_Initiated\\n', Pos_Pol_Template_Initiated)" << endl;
            ofs << in+ in+ "print('Dir_Pol_Initiated\\n', Dir_Pol_Initiated)" << endl;
        }

        ofs << in+ in+ "return Pos_Pol_Initiated, Pos_Pol_End_Initiated, Pos_Pol_Template_Initiated, Dir_Pol_Initiated, Count_NascentTarget_Initiated" << endl;
        ofs << endl;
        
        ofs << in+ "def Transcription_Elongation(self, Pos_Pol, Pos_Pol_End, Dir_Pol, Count_Nascent_Target, Rate, Freq_BB, Idx_PolSub, Idx_PolBB):" << endl;
        //    ofs << in+ in+ "dLength = np.matmul(SMatrix,Rate)
        ofs << in+ in+ "dLength = self.ApplySimTimeResolution(Rate)   # this is not necessarily true based on the reaction input" << endl;
        ofs << in+ in+ "Pos_Pol_Elongated = np.where(Pos_Pol >= 0, Pos_Pol + dLength * Dir_Pol, Pos_Pol)" << endl;
        ofs << in+ in+ "Pos_Pol_Trimmed = self.OverElongationCorrection_Transcription(Pos_Pol_Elongated, Pos_Pol_End, Dir_Pol)" << endl;
        ofs << endl;

        ofs << in+ in+ "# N_Elongated_PerPol = Pos_Pol_Trimmed - Pos_Pol * Dir_Pol" << endl;
        ofs << in+ in+ "N_Elongated_Total = np.array(np.sum(Pos_Pol_Trimmed - Pos_Pol, axis=1), ndmin=2).transpose()" << endl;
        ofs << endl;

        ofs << in+ in+ "Consumed_BB = self.BuildingBlockConsumption(Freq_BB, Count_Nascent_Target)" << endl;
        //ofs << in+ in+ "Consumed_BB = self.BuildingBlockConsumption(Freq_BB, N_Elongated_PerPol)" << endl;
        ofs << in+ in+ "# Update dCount for BuildingBlocks" << endl;
        ofs << in+ in+ "self.AddTodCount(Idx_PolBB, -Consumed_BB)" << endl;
        ofs << endl;

        // TODO: Update this with matrix calculation form
        ofs << in+ in+ "# Update dCount for Polymerase Reaction Substrates" << endl;
        ofs << in+ in+ "self.AddTodCount(Idx_PolSub, N_Elongated_Total)" << endl;
        ofs << endl;
        ofs << in+ in+ "return Pos_Pol_Trimmed" << endl;
        ofs << endl;
        
        ofs << in+ "def Transcription_Termination(self, Pos_Pol, Pos_Pol_End, Pos_Pol_Template, Dir_Pol, Count_NascentTarget, Idx_Target):" << endl;
        ofs << in+ in+ "# May be taken out to Sim Function" << endl;
        ofs << endl;

        ofs << in+ in+ "Idx_Pol_Completed = np.where(Pos_Pol == Pos_Pol_End)" << endl;
        ofs << in+ in+ "Idx_Template_Completed = Pos_Pol_Template[Idx_Pol_Completed]" << endl;
        ofs << endl;

        ofs << in+ in+ "Pos_Pol_Completed = SimF.ReplaceValueInArrayOnlyAtIdx(Pos_Pol, Idx_Pol_Completed, -1)   # Pol Position on Chromosome" << endl;
        ofs << in+ in+ "Pos_Pol_End_Completed = SimF.ReplaceValueInArrayOnlyAtIdx(Pos_Pol_End, Idx_Pol_Completed, 0)   # End chromosome position of the gene to which Pol is bound to" << endl;
        ofs << in+ in+ "Pos_Pol_Template_Completed = SimF.ReplaceValueInArrayOnlyAtIdx(Pos_Pol_Template, Idx_Pol_Completed, -1)   # Genes to which Pol is bound to" << endl;
        ofs << in+ in+ "Dir_Pol_Completed = SimF.ReplaceValueInArrayOnlyAtIdx(Dir_Pol, Idx_Pol_Completed, 0)   # Direction of Gene to which Pol is bound to" << endl;
        ofs << endl;

        if (Option.bDebug) {
            ofs << in+ in+ "print('Idx_Pol_Completed', Idx_Pol_Completed)" << endl;
            ofs << in+ in+ "print('Idx_Target_Completed', Idx_Template_Completed)" << endl;
            ofs << in+ in+ "print('Pos_Pol_Completed\\n', Pos_Pol_Completed)" << endl;
            ofs << in+ in+ "print('Pos_Pol_End_Completed\\n', Pos_Pol_End_Completed)" << endl;
            ofs << in+ in+ "print('Pos_Pol_Template_Completed\\n', Pos_Pol_Template_Completed)" << endl;
            ofs << in+ in+ "print('Dir_Pol_Completed\\n', Dir_Pol_Completed)" << endl;
            ofs << endl;
        }

        ofs << in+ in+ "Idx_Pol_Completed_Rows = Idx_Pol_Completed[0]" << endl;
        ofs << in+ in+ "Count_NascentTarget_Completed = SimF.AddValueToArrayOnlyAtIdx(Count_NascentTarget, (Idx_Pol_Completed_Rows, Idx_Template_Completed), -1)" << endl;
        if (Option.bDebug) {
            ofs << in+ in+ "print('Count_NascentTarget_Completed\\n', Count_NascentTarget_Completed)" << endl;
        }
        ofs << endl;

        ofs << in+ in+ "Idx_Mol_Completed = Idx_Target[Idx_Template_Completed]" << endl;
        if (Option.bDebug) {
            ofs << in+ in+ "print('Idx_Mol_Completed\\n', Idx_Mol_Completed)" << endl;
        }
        ofs << endl;

        ofs << in+ in+ "self.State.dCount_All = SimF.AddValueToArrayOnlyAtIdx(self.State.dCount_All, (Idx_Pol_Completed_Rows, Idx_Mol_Completed), 1)" << endl;
        if (Option.bDebug) {
            ofs << in+ in+ "print('self.State.dCount_All After\\n', self.State.dCount_All)" << endl;
            
        }
        ofs << endl;

        ofs << in+ in+ "return Pos_Pol_Completed, Pos_Pol_End_Completed, Pos_Pol_Template_Completed, Dir_Pol_Completed, Count_NascentTarget_Completed" << endl;
        ofs << endl;
    }

    if (!Ribosomes.empty()) {
        ofs << in+ "def Translation_GetAvailablePolymerases(self, Count_NascentTarget, Idx_Pol, PolThreshold):" << endl;
        ofs << in+ in+ "# Get available, active polymerase count - TO BE UPDATED with more regulatory algorithms" << endl;
        ofs << in+ in+ "Count_Pol = self.GetCount(Idx_Pol)" << endl;
        ofs << in+ in+ "Count_Reg = self.GetCount(Idx_Pol)" << endl;  // TODO: Implement regulation
        ofs << in+ in+ "Count_Pol_Reg_Complex = np.where(Count_Pol > Count_Reg, Count_Reg, Count_Pol)" << endl;
        ofs << in+ in+ "Count_Pol_Active = np.floor_divide(Count_Pol_Reg_Complex, 2).astype(int)" << endl;
        ofs << in+ in+ "Count_Pol_Occupied = np.reshape(np.sum(Count_NascentTarget, axis=1), [-1, 1]) * PolThreshold" << endl;
        ofs << in+ in+ "Count_Pol_Avail = Count_Pol_Active - Count_Pol_Occupied" << endl;
        ofs << in+ in+ "Count_Pol_FunctionalUnit = np.floor_divide(Count_Pol_Avail, PolThreshold)" << endl;
        ofs << in+ in+ "Count_Pol_Avail = np.where(Count_Pol_FunctionalUnit > 0, Count_Pol_FunctionalUnit, 0)" << endl;
        ofs << in+ in+ "return Count_Pol_Avail" << endl;
        ofs << endl;

        ofs << in+ "def Translation_GetInitiationSites(self, Idx_Template, Count_Template_Nascent):" << endl;
        ofs << in+ in+ "# Get final initiation weight by applying initiation site count" << endl;
        ofs << in+ in+ "Count_Template_Complete = self.GetCount(Idx_Template)" << endl;
        ofs << in+ in+ "Count_InitiationSitePerSite = Count_Template_Complete + Count_Template_Nascent" << endl;
        ofs << in+ in+ "return Count_InitiationSitePerSite" << endl;
        ofs << endl;

        // TODO: Update with more mechanistic algorithm later
        ofs << in+ "def Translation_Initiation(self, Pos_Pol, Pos_Pol_End, Pos_Pol_Template, Dir_Pol, Pos_Template_Start, Pos_Template_End, Dir_Template, Count_NascentTemplate, Count_NascentTarget, Idx_Pol, Idx_Template, Weight, PolThreshold):" << endl;
        ofs << in+ in+ "# Get Available Polymerase complex count" << endl;
        ofs << in+ in+ "Count_Pol_Avail = self.Translation_GetAvailablePolymerases(Count_NascentTarget, Idx_Pol, PolThreshold)" << endl;
        ofs << endl;

        ofs << in+ in+ "# Get Count and Weight of Translation initiation sites" << endl;
        ofs << in+ in+ "Count_InitiationSites_PerSite = self.Translation_GetInitiationSites(Idx_Template, Count_NascentTemplate)" << endl;
        ofs << in+ in+ "Count_InitiationSites_Total = np.reshape(np.sum(Count_InitiationSites_PerSite, axis=1), [-1, 1])" << endl;
        ofs << in+ in+ "Weight_InitiationSite = Count_InitiationSites_PerSite * Weight" << endl;
        ofs << in+ in+ "Count_ToBind = np.logical_and(Count_Pol_Avail > 0, Count_InitiationSites_Total > 0)" << endl;
        ofs << in+ in+ "Count_ToBindThisStep = np.where(Count_ToBind, 1, 0)   # Only allowing one binding at a time to use random function in matrix operation" << endl;
        ofs << endl;

        // TODO: Update to work on all rows
        // TODO: Incorporate Weighted system that works with arrays
        ofs << in+ in+ "# Get New polymerase binding sites." << endl;
        ofs << in+ in+ "Idx_MolToBind = np.reshape(np.random.randint(0, high=Weight_InitiationSite.shape[1], size=Weight_InitiationSite.shape[0]), [-1, 1])" << endl;
        ofs << endl;

        ofs << in+ in+ "# Apply New Binding to Count_NascentTarget" << endl;
        ofs << in+ in+ "Count_NascentTarget_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Count_NascentTarget, Idx_MolToBind, Count_ToBindThisStep) " << endl;
        ofs << endl;

        ofs << in+ in+ "# Apply New Binding at start position to Pos_Pol" << endl;
        ofs << in+ in+ "Idx_EmptyElement = np.reshape(SimF.GetIdxOfEmptyElement(Pos_Pol), [-1, 1]) " << endl;
        ofs << endl;

        if (Option.bDebug) {
            ofs << in+ in+ "print('Count_ToBind\\n', Count_ToBind.T)" << endl;
            ofs << in+ in+ "print('Count_ToBindThisStep\\n', Count_ToBindThisStep.T)" << endl;
            ofs << in+ in+ "print('Idx_MolToBind\\n', Idx_MolToBind)" << endl;
            ofs << in+ in+ "print('Idx_EmptyElement\\n', Idx_EmptyElement.T)" << endl;
            ofs << in+ in+ "print('Count_NascentTarget_Initiated\\n', Count_NascentTarget_Initiated)" << endl;
            ofs << endl;
        }

        ofs << in+ in+ "Pos_Template_Start_Initiated = Pos_Template_Start[Idx_MolToBind]" << endl;
        ofs << in+ in+ "Pos_Pol_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Pos_Pol, Idx_EmptyElement, Pos_Template_Start_Initiated * Count_ToBindThisStep + Count_ToBindThisStep) # Additional Count_ToBindThisStep to offset -1 value in the unbound element convention" << endl;
        ofs << endl;
        ofs << in+ in+ "Pos_Template_End_Initiated = Pos_Template_End[Idx_MolToBind]" << endl;
        ofs << in+ in+ "Pos_Pol_End_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Pos_Pol_End, Idx_EmptyElement, Pos_Template_End_Initiated * Count_ToBindThisStep) " << endl;
        ofs << endl;
        ofs << in+ in+ "Pos_Pol_Template_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Pos_Pol_Template, Idx_EmptyElement, Idx_MolToBind * Count_ToBindThisStep + Count_ToBindThisStep) # Additional Count_ToBindThisStep to offset -1 value in the unbound element convention" << endl;
        ofs << endl;
        ofs << in+ in+ "Dir_Template_Initiated = Dir_Template[Idx_MolToBind]" << endl;
        ofs << in+ in+ "Dir_Pol_Initiated = SimF.AddValueToArrayAllRowsAtIdx(Dir_Pol, Idx_EmptyElement, Dir_Template_Initiated * Count_ToBindThisStep) " << endl;
        ofs << endl;

        if (Option.bDebug) {
            ofs << in+ in+ "print('Pos_Pol_Initiated\\n', Pos_Pol_Initiated)" << endl;
            ofs << in+ in+ "print('Pos_Pol_End_Initiated\\n', Pos_Pol_End_Initiated)" << endl;
            ofs << in+ in+ "print('Pos_Pol_Template_Initiated\\n', Pos_Pol_Template_Initiated)" << endl;
            ofs << in+ in+ "print('Dir_Pol_Initiated\\n', Dir_Pol_Initiated)" << endl;
        }

        ofs << in+ in+ "return Pos_Pol_Initiated, Pos_Pol_End_Initiated, Pos_Pol_Template_Initiated, Dir_Pol_Initiated, Count_NascentTarget_Initiated" << endl;
        ofs << endl;
        
        ofs << in+ "def Translation_Elongation(self, Pos_Pol, Pos_Pol_End, Dir_Pol, Count_Nascent_Target, Rate, Freq_BB, Idx_PolSub, Idx_PolBB):" << endl;
        //    ofs << in+ in+ "dLength = np.matmul(SMatrix,Rate)
        ofs << in+ in+ "dLength = self.ApplySimTimeResolution(Rate)   # this is not necessarily true based on the reaction input" << endl;
        ofs << in+ in+ "Pos_Pol_Elongated = np.where(Pos_Pol >= 0, Pos_Pol + dLength * Dir_Pol, Pos_Pol)" << endl;
        ofs << in+ in+ "Pos_Pol_Trimmed = self.OverElongationCorrection_Translation(Pos_Pol_Elongated, Pos_Pol_End, Dir_Pol)" << endl;
        ofs << endl;

        ofs << in+ in+ "# N_Elongated_PerPol = Pos_Pol_Trimmed - Pos_Pol * Dir_Pol" << endl;
        ofs << in+ in+ "N_Elongated_Total = np.array(np.sum(Pos_Pol_Trimmed - Pos_Pol, axis=1), ndmin=2).transpose()" << endl;
        ofs << endl;

        ofs << in+ in+ "Consumed_BB = self.BuildingBlockConsumption(Freq_BB, Count_Nascent_Target)" << endl;
        //ofs << in+ in+ "Consumed_BB = self.BuildingBlockConsumption(Freq_BB, N_Elongated_PerPol)" << endl;
        ofs << in+ in+ "# Update dCount for BuildingBlocks" << endl;
        ofs << in+ in+ "self.AddTodCount(Idx_PolBB, -Consumed_BB)" << endl;
        ofs << endl;

        // TODO: Update this with matrix calculation form
        ofs << in+ in+ "# Update dCount for Polymerase Reaction Substrates" << endl;
        ofs << in+ in+ "self.AddTodCount(Idx_PolSub, N_Elongated_Total)" << endl;
        ofs << endl;
        ofs << in+ in+ "return Pos_Pol_Trimmed" << endl;
        ofs << endl;
        
        ofs << in+ "def Translation_Termination(self, Pos_Pol, Pos_Pol_End, Pos_Pol_Template, Dir_Pol, Count_NascentTarget, Idx_Target):" << endl;
        ofs << in+ in+ "# May be taken out to Sim Function" << endl;
        ofs << endl;

        ofs << in+ in+ "Idx_Pol_Completed = np.where(Pos_Pol == Pos_Pol_End)" << endl;
        ofs << in+ in+ "Idx_Template_Completed = Pos_Pol_Template[Idx_Pol_Completed]" << endl;
        ofs << endl;

        ofs << in+ in+ "Pos_Pol_Completed = SimF.ReplaceValueInArrayOnlyAtIdx(Pos_Pol, Idx_Pol_Completed, -1)   # Pol Position on Chromosome" << endl;
        ofs << in+ in+ "Pos_Pol_End_Completed = SimF.ReplaceValueInArrayOnlyAtIdx(Pos_Pol_End, Idx_Pol_Completed, 0)   # End chromosome position of the gene to which Pol is bound to" << endl;
        ofs << in+ in+ "Pos_Pol_Template_Completed = SimF.ReplaceValueInArrayOnlyAtIdx(Pos_Pol_Template, Idx_Pol_Completed, -1)   # Genes to which Pol is bound to" << endl;
        ofs << in+ in+ "Dir_Pol_Completed = SimF.ReplaceValueInArrayOnlyAtIdx(Dir_Pol, Idx_Pol_Completed, 0)   # Direction of Gene to which Pol is bound to" << endl;
        ofs << endl;

        if (Option.bDebug) {
            ofs << in+ in+ "print('Idx_Pol_Completed', Idx_Pol_Completed)" << endl;
            ofs << in+ in+ "print('Idx_Target_Completed', Idx_Template_Completed)" << endl;
            ofs << in+ in+ "print('Pos_Pol_Completed\\n', Pos_Pol_Completed)" << endl;
            ofs << in+ in+ "print('Pos_Pol_End_Completed\\n', Pos_Pol_End_Completed)" << endl;
            ofs << in+ in+ "print('Pos_Pol_Template_Completed\\n', Pos_Pol_Template_Completed)" << endl;
            ofs << in+ in+ "print('Dir_Pol_Completed\\n', Dir_Pol_Completed)" << endl;
            ofs << endl;
        }

        ofs << in+ in+ "Idx_Pol_Completed_Rows = Idx_Pol_Completed[0]" << endl;
        ofs << in+ in+ "Count_NascentTarget_Completed = SimF.AddValueToArrayOnlyAtIdx(Count_NascentTarget, (Idx_Pol_Completed_Rows, Idx_Template_Completed), -1)" << endl;
        if (Option.bDebug) {
            ofs << in+ in+ "print('Count_NascentTarget_Completed\\n', Count_NascentTarget_Completed)" << endl;
        }
        ofs << endl;

        ofs << in+ in+ "Idx_Mol_Completed = Idx_Target[Idx_Template_Completed]" << endl;
        if (Option.bDebug) {
            ofs << in+ in+ "print('Idx_Mol_Completed\\n', Idx_Mol_Completed)" << endl;
        }
        ofs << endl;

        ofs << in+ in+ "self.State.dCount_All = SimF.AddValueToArrayOnlyAtIdx(self.State.dCount_All, (Idx_Pol_Completed_Rows, Idx_Mol_Completed), 1)" << endl;
        if (Option.bDebug) {
            ofs << in+ in+ "print('self.State.dCount_All After\\n', self.State.dCount_All)" << endl;
            
        }
        ofs << endl;

        ofs << in+ in+ "return Pos_Pol_Completed, Pos_Pol_End_Completed, Pos_Pol_Template_Completed, Dir_Pol_Completed, Count_NascentTarget_Completed" << endl;
        ofs << endl;
    }

    ofs << in+ "def EvaluateChemotaxisThreshold(self):" << endl;
    ofs << in+ in+ "ThresholdedMolecules = self.GetCount(self.State.Idx_Count_Threshold).transpose()" << endl;
    ofs << in+ in+ "return np.array(ThresholdedMolecules) < self.State.Pos_Threshold" << endl;
    ofs << endl;
    
    ofs << in+ "def CellDivision(self):" << endl;
    if (!RNAPs.empty()) {
        ofs << in+ in+ "# Temporary rate ramping for viewing" << endl;
        ofs << in+ in+ "self.State.Rate_Transcription = np.array([1e4])" << endl;

    }

    if (!DNAPs.empty()) {
        ofs << in+ in+ "bDividingCells = np.count_nonzero(np.where(self.GetCount(self.State.Idx_Template_Replication) >= 2, 1, 0), axis=1)" << endl;
//        ofs << in+ in+ "bDividingCells[0] = 1" << endl;
        ofs << in+ in+ "Idx_DividingCells = np.where(bDividingCells > 0)" << endl;
        // One liner:
        //        ofs << in+ in+ "Idx_DividingCells = np.where(np.count_nonzero(np.where(self.GetCount(self.State.Idx_Template_Replication) >= 2), axis=1) > 0)" << endl;
        ofs << endl;

        ofs << in+ in+ "# Temporary rate ramping for viewing" << endl;
        ofs << in+ in+ "self.State.Rate_Replication = np.array([1e6])" << endl;

        if (Option.bDebug) {
            ofs << in+ in+ "#Debugging" << endl;
            //ofs << in+ in+ "self.State.Rate_Replication = np.array([1e6])" << endl;

            ofs << in+ in+ "self.Debug_PrintSimStepTime(end='\\n')" << endl;
            ofs << in+ in+ "print('Chromosomes:     ', self.GetCount(self.State.Idx_Template_Replication).transpose())" << endl;
            ofs << in+ in+ "print('Dividing T/F:    ', bDividingCells.transpose())" << endl;
            ofs << in+ in+ "print('PolCounts:       ', self.GetCount(self.State.Idx_Pol_Replication).transpose())" << endl;
            ofs << in+ in+ "print('ReplicationRate: ', self.GetReplicationCompletionRate('E').transpose())" << endl;
            if (!Context.ThresholdList.empty()) {
                ofs << in+ in+ "print('Pos_Threshold:   ', self.State.Pos_Threshold)" << endl;
            }
            ofs << endl;
        }


        ofs << in+ in+ "DividingCell_Count_All = self.State.Count_All[Idx_DividingCells]" << endl;
        ofs << in+ in+ "Distributed_Count_All = DividingCell_Count_All / 2" << endl;
        ofs << in+ in+ "self.State.Count_All = np.vstack([self.State.Count_All, Distributed_Count_All])" << endl;
        ofs << in+ in+ "self.State.Count_All[Idx_DividingCells] = Distributed_Count_All" << endl;
        ofs << endl;
        ofs << in+ in+ "DividingCell_dCount_All = self.State.dCount_All[Idx_DividingCells]" << endl;
        ofs << in+ in+ "Distributed_dCount_All = DividingCell_dCount_All / 2" << endl;
        ofs << in+ in+ "self.State.dCount_All = np.vstack([self.State.dCount_All, Distributed_dCount_All])" << endl;
        ofs << in+ in+ "self.State.dCount_All[Idx_DividingCells] = Distributed_dCount_All" << endl;
        ofs << endl;
        // TODO: to be updated for concurrent replication mode
        ofs << in+ in+ "DividingCell_Len_NascentChromosome = self.State.Len_NascentChromosome[Idx_DividingCells]" << endl;
        ofs << in+ in+ "self.State.Len_NascentChromosome = np.vstack([self.State.Len_NascentChromosome, DividingCell_Len_NascentChromosome])" << endl;
        ofs << endl;
        if (!ObjLoc.empty()) {
            ofs << in+ in+ "BodyLength = 20" << endl;
            ofs << in+ in+ "dX =  np.cos(self.State.Pos_Angle) * -BodyLength" << endl;
            ofs << in+ in+ "dY = -np.sin(self.State.Pos_Angle) * -BodyLength" << endl;
            
            ofs << in+ in+ "Pos_X = self.State.Pos_X[Idx_DividingCells]" << endl;
//            ofs << in+ in+ "MotherCell_Pos_X = Pos_X - dX[Idx_DividingCells] /2" << endl;
//            ofs << in+ in+ "DaughterCell_Pos_X = Pos_X + dX[Idx_DividingCells] /2" << endl;
            ofs << in+ in+ "DaughterCell_Pos_X = Pos_X + dX[Idx_DividingCells] * 1.2" << endl;
//            ofs << in+ in+ "self.State.Pos_X[Idx_DividingCells] = MotherCell_Pos_X" << endl;
            ofs << in+ in+ "self.State.Pos_X = np.concatenate([self.State.Pos_X, DaughterCell_Pos_X])" << endl;
            ofs << endl;
            ofs << in+ in+ "Pos_Y = self.State.Pos_Y[Idx_DividingCells]" << endl;
//            ofs << in+ in+ "MotherCell_Pos_Y = Pos_Y - dY[Idx_DividingCells] /2" << endl;
//            ofs << in+ in+ "DaughterCell_Pos_Y = Pos_Y + dY[Idx_DividingCells] /2" << endl;
            ofs << in+ in+ "DaughterCell_Pos_Y = Pos_Y + dY[Idx_DividingCells] * 1.2" << endl;
//            ofs << in+ in+ "self.State.Pos_Y[Idx_DividingCells] = MotherCell_Pos_Y" << endl;
            ofs << in+ in+ "self.State.Pos_Y = np.concatenate([self.State.Pos_Y, DaughterCell_Pos_Y])" << endl;
            ofs << endl;
            ofs << in+ in+ "Pos_Angle = self.State.Pos_Angle[Idx_DividingCells]" << endl;
            ofs << in+ in+ "self.State.Pos_Angle = np.concatenate([self.State.Pos_Angle, Pos_Angle])" << endl;
//            ofs << in+ in+ "MotherCell_Pos_Angle = Pos_Angle - dX" << endl;
//            ofs << in+ in+ "DaughterCell_Pos_Angle = Pos_Angle + dX" << endl;
//            ofs << in+ in+ "self.State.Pos_Angle[Idx_DividingCells] = MotherCell_Pos_Angle" << endl;
//            ofs << in+ in+ "self.State.Pos_Angle = np.concatenate([self.State.Pos_Angle, DaughterCell_Pos_Angle])" << endl;
            ofs << endl;
            ofs << in+ in+ "self.State.Pos_Name2Idx['E'] = np.array(range(self.State.Pos_X.shape[0]))" << endl;
            ofs << endl;
//            ofs << in+ in+ "DividingCell_Pos_X = self.State.Pos_X[Idx_DividingCells]" << endl;
//            ofs << in+ in+ "self.State.Pos_X = np.vstack([self.State.Pos_X, DividingCell_Pos_X])" << endl;

        }
        if (!Context.ThresholdList.empty()) {
            ofs << in+ in+ "DividingCell_Pos_Threshold = self.State.Pos_Threshold[:, Idx_DividingCells][0]" << endl;
            ofs << in+ in+ "self.State.Pos_Threshold = np.concatenate([self.State.Pos_Threshold, DividingCell_Pos_Threshold], axis=1)" << endl;
            ofs << endl;
        }
        

//        self.Pos_Names = ['E', ]
//        self.Idx_Pos_E = np.array([0, 1, ])
//        self.Pos_Name2Idx['E'] = self.Idx_Pos_E



//        ofs << in+ in+ "DividingCell_Pos = self.State.Pos"
    } else {
        ofs << in+ in+ "pass" << endl;
    } ofs << endl;


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
    } else { ofs << in+ in+ "self.Debug_Idx_Pos = [0]" << endl; }
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
