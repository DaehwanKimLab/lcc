#include "writer.h"
#include <algorithm>

using namespace std;

// function that returns a threshold value

// run while i < X (500 seconds?)

// prev - now / now < 0.0001

void FWriter::SimThreshold(FMolecule * ThresholdMol)
{
    std::string Name_ThresholdMol = "Threshold_" + ThresholdMol->Name;
    std::string Name_SimThreshold = "SimThreshold_" + ThresholdMol->Name;
    std::string FileName_SimThreshold = Name_SimThreshold + ".py";

    auto Pathway = Context.GetAssociatedPathway_ReactionList(ThresholdMol);

    auto ReactionList_Threshold = Context.GetSubList_ReactionList("All", Pathway);
    auto MoleculeList_Threshold = Context.GetSubListByReactionList_MoleculeList(ReactionList_Threshold);

    std::cout << "Generating " << Name_SimThreshold << "..." << std::endl;

    // write SimThreshold_ThresholdMol.py
    std::ofstream ofs(FileName_SimThreshold);
    std::string endl = "\n";

    // IMPORT
    ofs << "import os, sys" << endl;
    ofs << "import numpy as np" << endl;
    ofs << "import SimFunctions as SimF" << endl;
    // ofs << "import SimIdx as idx" << endl;
    ofs << "from argparse import ArgumentParser" << endl;
    ofs << endl;

    //TODO: Take options from simmodule cmd line
    ofs << "# Temporary global variables" << endl;
    ofs << "N_SimSteps = " << 30000 << endl;
    ofs << "SimStepTimeResolution = " << 100 << endl;
    ofs << "DisplayCount = 0   # 0: Display Counts only, 1: Display both Counts & dCounts" << endl;
    ofs << "Unit = 'nM'   # nM or uM supported for now" << endl;
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

    ofs << in+ in+ "self.Count_All = np.zeros([1, " << MoleculeList_Threshold.size() << "])" << endl;
    ofs << in+ in+ "self.dCount_All = np.zeros([1, " << MoleculeList_Threshold.size() << "])" << endl;
    ofs << endl;

    // for standard reactions
    std::vector<std::string> ReactionTypes;
    ReactionTypes.emplace_back("Standard_Unregulated");
    ReactionTypes.emplace_back("Standard_Inhibition_Allosteric");
    ReactionTypes.emplace_back("Standard_Activation_Allosteric");

    std::vector<std::string> StandardReactionTypes = ReactionTypes; // to reuse later

    for (auto& Type : ReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Initialize_StandardReaction(ofs, Type, Pathway);
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
            Initialize_EnzymeReaction(ofs, Type, Pathway);
        }
    }

    // for polymerase reactions (Template-based)
    std::vector<FPolymerase *> PolymeraseList = Context.GetList_Polymerase_MoleculeList();
//    std::vector<std::string> PolymeraseNames = Context.GetNames_PolymeraseList(PolymeraseList);
    std::vector<FPolymeraseReaction *> PolymeraseReactionList = Context.GetList_Polymerase_ReactionList();

    if (!PolymeraseList.empty()) {
        for (auto& Polymerase : PolymeraseList) {
            Initialize_PolymeraseReaction(ofs, Polymerase); // TODO: Update polymerase reactions
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
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
        if (!ReactionSubList.empty()) {
            Initialize_TransporterReaction(ofs, Type);
        }
    }

    ofs << in+ "def Initialize(self):" << endl;
    ofs << in+ in+ "self.Vol = 1" << endl;
    ofs << endl;

    // Print SetUp_StandardReaction for each Reaction Type
    for (auto& Type : StandardReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
        if (!ReactionSubList.empty()) {
            SetUp_StandardReaction(ofs, Type, ReactionSubList);
        }
    }

    // Print SetUp_EnzymeReaction for each Reaction Type
    for (auto& Type : EnzReactionTypes) {
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
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
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
        if (!ReactionSubList.empty()) {
            SetUp_TransporterReaction(ofs, Type, ReactionSubList);
        }
    }

    // for legends

    std::vector<std::string> MolNames;

    for (auto Molecule : MoleculeList_Threshold){
        if (std::find(MolNames.begin(), MolNames.end(), Molecule->Name) == MolNames.end()) {
            MolNames.push_back(Molecule->Name);
        }
    }

    ofs << in+ in+ "self.Mol_Names = [" << Utils::JoinStr2Str(MolNames) << "]" << endl;
    ofs << in+ in+ "self.Mol_Name2Idx = {";
    for (int i = 0; i < MolNames.size(); i++) {
        ofs << "'" << MolNames[i] << "' : np.array([" << i << "]), ";
    } ofs << "}" << endl;
    ofs << in+ in+ "self.Legends = ['SimStep', 'Vol', " << Utils::JoinStr2Str(MolNames) << "]" << endl;
    ofs << endl;

    // Initialize all molecule counts

    std::vector<int> Idx_Mol;

    for (int i = 0; i < MoleculeList_Threshold.size(); i++) {
        if (Utils::is_class_of<FMolecule, FMolecule>(MoleculeList_Threshold[i])) {
            Idx_Mol.push_back(i);
        }
    }

    std::vector<float> InitialCount_Molecules;
    std::vector<float> MolarityFactor_Molecules;

    // This may be leaky if redundant molecules are modified outside the current pathway
    for (auto& molecule : MoleculeList_Threshold) {
        float Count = Context.GetInitialCountByName_CountList(molecule->Name);
        float MolarityFactor = Context.GetMolarityFactorByName_CountList(molecule->Name);
        InitialCount_Molecules.push_back(Count);
        MolarityFactor_Molecules.push_back(MolarityFactor);
    }

    ofs << in+ in+ "Idx_Mol = np.array([[" << Utils::JoinInt2Str_Idx(Idx_Mol) << "]])" << endl;
    ofs << in+ in+ "Count_Mol = np.array([[" << Utils::JoinFloat2Str(InitialCount_Molecules) << "]])" << endl;

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


    int i = 0;
    int i_SimStep = i + 1;
    ofs << in+ in+ "Data[:, " << i << ":" << i_SimStep << "] = Time" << endl;

    int i_Vol = i_SimStep + 1;
    ofs << in+ in+ "Data[:, " << i_SimStep << ":" << i_Vol << "] = self.Vol" << endl;

    int i_Count_Mol = i_Vol + MolNames.size();
    ofs << in+ in+ "Data[:, " << i_Vol << ":" << i_Count_Mol << "] = self.Count_All" << endl;

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

    // Restore --> to CountList method
    std::vector<std::string> Names;
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;
        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }
        if (count->End == -1) { // indicator for fixed
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                ofs << in+ in+ "self.Idx_Restore_" << Name << " = None" << endl;
                Names.push_back(Name);
            }
        }
    }

    // Event --> to CountList method
    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;
        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }
        if ((count->End >= 0) & (count->Begin != count->End)) {
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                ofs << in+ in+ "self.Idx_Event_" << Name << " = None" << endl;
                Names.push_back(Name);
            }
        }
    }

    // Threshold-related
    ofs << in+ in+ "self.ThresholdValue = 0" << endl;
    ofs << in+ in+ "self.Count_Prev = 0" << endl;
    ofs << in+ in+ "self.Idx_Count_Homeostasis = None" << endl;
    ofs << endl;

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

    // Restore --> to CountList method
    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;
        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }
        if (count->End == -1) { // indicator for fixed
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                int Idx = Context.GetIdxByName_MoleculeList(Name);
                ofs << in+ in+ "self.Idx_Restore_" << Name << " = np.array([" << std::to_string(Idx) << "])" << endl;
                Names.push_back(Name);
            }
        }
    }

    // Event --> to CountList method
    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;
        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }
        if ((count->End >= 0) & (count->Begin != count->End)) { //
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                int Idx = Context.GetIdxByName_MoleculeList(Name);
                ofs << in+ in+ "self.Idx_Event_" << Name << " = np.array([" << std::to_string(Idx) << "], ndmin=2)" << endl;
                Names.push_back(Name);
            }
        }
    }
    ofs << endl;

    ofs << in+ in+ "self.State.Initialize()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Threshold-related" << endl;
    ofs << in+ in+ "self.Idx_Count_Homeostasis = self.GetMolIdx(" << ThresholdMol->Name << ")" << endl;
    ofs << endl;

    ofs << in+ in+ "# Legend Export" << endl;
    ofs << in+ in+ "self.Dataset.Legend = self.State.ExportLegend()" << endl;
    ofs << in+ in+ "self.DataManager.SetLegend(self.Dataset.Legend)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Data Export" << endl;
    ofs << in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Debugging" << endl;
    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#"; }
    ofs << "self.Debug_SetIdxMoleculesToTrack()" << endl;
    if (!Context.LocationList.empty()) {
        ofs << in+ in;
//        if (!Option.bDebug) { ofs << "#"; }
        ofs << "self.Debug_SetIdxDistAndPosToTrack()" << endl;
    }
    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#"; }
    ofs << "self.Debug_SetUnit(Unit)" << endl;
    ofs << endl;


    ofs << in+ "def SimLoop_WithoutSpatialSimulation(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

//    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#";}
//    ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
//    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#";}
//    ofs << "self.Debug_PrintCounts(DisplayCount)" << endl;
//    ofs << endl;

    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;

    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ "def Run_FindThreshold(self):" << endl;
    ofs << in+ in+ "print('Simulation Run_FindThreshold(" << ThresholdMol->Name << ") Begins...')" << endl;
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

    ofs << in+ in+ "print('Simulation Run_FindThreshold(" << ThresholdMol->Name << ")  Completed')" << endl;
    ofs << endl;

    // TODO: scale up exporting
    ofs << in+ "def ExportData(self):" << endl;
    if (Context.LocationList.empty()) {
        ofs << in+ in+ "self.Dataset.Data = self.State.ExportData(self.SimStep/self.SimTimeResolutionPerSecond)" << endl;
        ofs << in+ in+ "self.DataManager.Add(self.Dataset.Data)" << endl;
    } else { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "def SaveState(self, FileName):" << endl;
    ofs << in+ in+ "self.DataManager.SaveToNpyFile(self.GetCount_All(), FileName)" << endl;
    ofs << endl;

    ofs << in+ "def ApplySimTimeResolution(self, Rate):" << endl;
    ofs << in+ in+ "return Rate / self.SimTimeResolutionPerSecond" << endl;
    ofs << endl;

    // Restore
    ofs << in+ "def RestoreMoleculeCount(self):" << endl;

    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;

        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }

        if (count->End == -1) { // indicator for fixed
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                std::string Amount;
                float MolarityFactor = Context.GetMolarityFactorByName_CountList(Name);

                if (MolarityFactor)      { Amount = Utils::SciFloat2Str(count->Amount) + " * self.State.Vol"; }
                else                     { Amount = Utils::SciFloat2Str(count->Amount); }

                ofs << in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Restore_" << Name << ".reshape(1, -1), " << Amount << ", axis=1)" << endl;
                Names.push_back(Name);
            }
        }
    }

    if (Names.empty()) {
        ofs << in+ in+ "pass" << endl;
    }

    // Event
    ofs << in+ "def TriggerEventMoleculeCount(self):" << endl;
    ofs << in+ in+ "Time = self.SimStep / self.SimTimeResolutionPerSecond" << endl;
    bool ElseSwitch = false;

    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;

        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }

        std::string Amount;
        float MolarityFactor = Context.GetMolarityFactorByName_CountList(Name);

        if (MolarityFactor) { Amount = Utils::SciFloat2Str(count->Amount) + " * self.State.Vol"; }
        else                { Amount = Utils::SciFloat2Str(count->Amount); }

        if ((count->End >= 0) & (count->Begin != count->End)) { //
            ofs << in+ in+ "if (Time >= " << Utils::SciFloat2Str(count->Begin) << ") & (Time < " << Utils::SciFloat2Str(count->End) << "):" << endl;
            ofs << in+ in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Event_" << Name << ", " << Amount << ", axis=1)" << endl;
            Names.push_back(Name);
        }
    }
    ofs << endl;

    if (!Context.ThresholdList.empty()) {
        ofs << in+ "def FindThresholds(self):" << endl;
        for (auto ThresholdMol : Context.ThresholdList) {
//            SimThreshold(ThresholdMol);
        }

        // temporary code
        ofs << in+ "def GetNewThresholdValues(self):" << endl;
        float ThresholdFactor = 0.999;
        ofs << in+ in+ "return self.GetCount_All()[:, self.Idx_Count_Homeostasis].transpose() * " << Utils::SciFloat2Str(ThresholdFactor) << endl;
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

        ofs << in+ "def UpdateThreshold(self):" << endl;
        float HomeostasisFactor = 1e-7;
        ofs << in+ in+ "Count_Now = self.GetCount_All()[:, self.Idx_Count_Homeostasis].transpose()" << endl;
        ofs << in+ in+ "bSteadyState = abs(Count_Now - self.Count_Prev) / Count_Now < " << Utils::SciFloat2Str(HomeostasisFactor) << endl;
        ofs << in+ in+ "Idx_SteadyState = np.where(bSteadyState)" << endl;
        ofs << in+ in+ "self.State.Pos_Threshold[Idx_SteadyState] = Count_Now[Idx_SteadyState]" << endl;
        //    ofs << in+ in+ "if np.any(Idx_SteadyState):" << endl;
        //    ofs << in+ in+ in+ "print('[Homeostasis] Steady state has been updated.')" << endl;
        ofs << in+ in+ "self.Count_Prev = Count_Now" << endl;
        ofs << endl;

    }

    float HomeostasisFactor = 1e-7;
    ofs << in+ "def Homeostasis(self, MoleculeNames=None):" << endl;
    // TODO: Get homeostasis info from FPathway,
    ofs << in+ in+ "bNotHomeostasis = True" << endl;
    ofs << in+ in+ "self.SetIdxForHomeostasis(MoleculeNames)" << endl;
    ofs << in+ in+ "MoleculeNames_Str = ListToStr(MoleculeNames, '_')" << endl;

    ofs << in+ in+ "# Temporary homeostasis save file" << endl;
    ofs << in+ in+ "FileName_Homeostasis = 'SimSave_Homeostasis_%s" << int(Numbers::RandomNumber(0, 1000)) << "' % MoleculeNames_Str" << endl;
    ofs << endl;
    ofs << in+ in+ "if os.path.exists('%s.npy' % FileName_Homeostasis):" << endl;
    ofs << in+ in+ in+ "print('[Homeostasis] Loading previous steady state...')" << endl;
    ofs << in+ in+ in+ "self.State.Count_All = np.load('%s.npy' % FileName_Homeostasis)" << endl;
    ofs << in+ in+ in+ "print('[Homeostasis] Steady state has been achieved!')" << endl;
    ofs << in+ in+ in+ "if MoleculeNames:" << endl;
    ofs << in+ in+ in+ in+ "self.SetThreshold(True)" << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "print('[Homeostasis] Running simulation to achieve steady state...')" << endl;

    ofs << in+ in+ in+ "while (bNotHomeostasis):" << endl;
    ofs << in+ in+ in+ in+ "self.SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution()" << endl;
    ofs << in+ in+ in+ in+ "Count_Now = self.GetCount_All()[:, self.Idx_Count_Homeostasis].transpose()" << endl;
    ofs << in+ in+ in+ in+ "bSteadyState = abs(Count_Now - self.Count_Prev) / Count_Now < " << Utils::SciFloat2Str(HomeostasisFactor) << endl;
    ofs << in+ in+ in+ in+ "if MoleculeNames:" << endl;
    ofs << in+ in+ in+ in+ in+ "self.SetThreshold(bSteadyState)" << endl;
    ofs << in+ in+ in+ in+ "if np.all(Count_Now > 0) and np.all(bSteadyState):" << endl;
    ofs << in+ in+ in+ in+ in+ "if MoleculeNames:" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "self.PrintThreshold(MoleculeNames_Str)" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "print('Current: ')" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
    ofs << in+ in+ in+ in+ in+ "bNotHomeostasis = False" << endl;
    ofs << in+ in+ in+ in+ in+ "print('[Homeostasis] Steady state has been achieved.')" << endl;
    ofs << in+ in+ in+ in+ in+ "self.SaveState(FileName_Homeostasis)" << endl;
    ofs << in+ in+ in+ in+ in+ "print('[Homeostasis] Steady state has been saved: %s' % FileName_Homeostasis)" << endl;
    ofs << in+ in+ in+ in+ "else:" << endl;
    ofs << in+ in+ in+ in+ in+ "self.Count_Prev = Count_Now" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ in+ "# Trigger Event on Substrate Count" << endl;
    ofs << in+ in+ in+ in+ in+ "self.TriggerEventMoleculeCount()" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ in+ in+ "self.ExportData()" << endl;
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

    ofs << in+ "def NonSpatialSimulation(self):" << endl;
    bool PassSwitch = true;
    if (!StandardReactionTypes.empty()) {
        for (auto &Type: StandardReactionTypes) {
            std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
            if (!ReactionSubList.empty()) {
                ofs << in+ in+ "self.StandardReactions()" << endl;
                ofs << endl;
                PassSwitch = false;
                break;
            }
        }
    }

    if (!EnzReactionTypes.empty()) {
        for (auto &Type: EnzReactionTypes) {
            std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
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
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
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
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
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
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
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
        std::vector<FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type, Pathway);
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
    ofs << in+ in+ "return SimF.CountToConc(self.State.Count_All[:, Idx])" << endl;
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

    // TODO: Use SimIdx in the future
    ofs << in+ "# Temporary routines" << endl;

    ofs << in+ "def GetMolIdx(self, Name):" << endl;
    ofs << in+ in+ "return self.State.Mol_Name2Idx[Name]" << endl;
    ofs << endl;

    ofs << in+ "def GetCountByName(self, Name):" << endl;
    ofs << in+ in+ "return self.GetCount(self.GetMolIdx(Name))" << endl;
    ofs << endl;

    ofs << in+ "def GetConcentrationByName(self, Name):" << endl;
    ofs << in+ in+ "return self.GetConcentration(self.GetMolIdx(Name), self.State.Vol)" << endl;
    ofs << endl;

    ofs << in+ "def GetThreshold(self, Pos_Idx):" << endl;
    ofs << in+ in+ "return self.State.Pos_Threshold[:, Pos_Idx]" << endl;
    ofs << endl;

    ofs << in+ "def GetPosIdx(self, Name):" << endl;
    ofs << in+ in+ "return self.State.Pos_Name2Idx[Name]" << endl;
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
    ofs << in+ in+ "#Debug_Names_Molecules = ['Am', 'AmL', 'L']" << endl;
    ofs << in+ in+ "#Debug_Names_Molecules = ['Am', 'AmL', 'L', 'qL', 'pc_qL']" << endl;
    ofs << endl;
    ofs << in+ in+ "if Debug_Names_Molecules == []:" << endl;
    ofs << in+ in+ in+ "Debug_Names_Molecules = self.State.GetMolNames()" << endl;
    ofs << endl;
    ofs << in+ in+ "for Name in Debug_Names_Molecules:" << endl;
    ofs << in+ in+ in+ "self.Debug_Idx_Molecules.append(self.State.GetMolNames().index(Name))" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCounts(self, Switch):" << endl;
    ofs << in+ in+ "self.Debug_PrintCountsForSinglePosition(Switch)" << endl;
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
    ofs << in+ in+ "print(self.SimStep, '(', round(Time,3), 's)', end=end)" << endl;
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
    ofs << in+ "Simulation.Run_FindThreshold()" << endl;
    ofs << in+ "return Simulation.ThresholdValue" << endl;
    ofs << endl;

    // API to get the threshold value
    ofs << "def DetermineThreshold():" << endl;
    ofs << in+ "return main()" << endl;
    ofs << endl;

    // MAIN
    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "parser = ArgumentParser()" << endl;
    ofs << in+ "parser.add_argument('--save-fig', dest='save_fname', type=str, help='Save figure to file')" << endl;
    ofs << in+ "args = parser.parse_args()" << endl;
    ofs << in+ "main()" << endl;
    ofs << endl;

    std::cout << "  Threshold search program has been generated: " << FileName_SimThreshold;
}
