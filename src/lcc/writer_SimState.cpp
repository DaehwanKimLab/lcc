#include "writer.h"
#include <algorithm>

using namespace std;

void FWriter::SimState(int Map_Width, int Map_Height)
{
    std::cout << "Generating simulation state..." << std::endl;

    // write SimState.py
    std::ofstream ofs(Option.SimStateFile.c_str());
    std::string endl = "\n";

    // IMPORT
    ofs << "import os, sys" << endl;
    ofs << endl;

    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "print('\\n%%%%%% Run " << Option.SimExecutorFile << " to execute simulation. %%%%%%')" << endl;
    ofs << in+ "sys.exit()" << endl;
    ofs << endl;

    ofs << "import numpy as np" << endl;
    ofs << "import csv" << endl;
    ofs << "import SimFunctions as SimF" << endl;
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
    ofs << in+ in+ "self.VolInit = None" << endl;
    ofs << in+ in+ "self.Vol = None" << endl;
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

    std::vector<std::vector<FMolecule *>> PolymeraseTypes = {DNAPs, RNAPs, Ribosomes};
    std::string Name_mRNASubIdx = "Idx_mRNAInRNA";

    // Conditional initialization is deactivated for unconditional initialization of variables used to communicate to visualization

    //if (!PolymeraseList.empty()) {
    //    Initialize_PolymeraseReaction_Matrix(ofs, PolymeraseTypes);

    //    if (!DNAPs.empty())        { Initialize_PolymeraseReaction_DNAP(ofs, DNAPs); }
    //    if (!RNAPs.empty())        { Initialize_PolymeraseReaction_RNAP(ofs, RNAPs); }
    //    if (!Ribosomes.empty())    { Initialize_PolymeraseReaction_Ribosome(ofs, Ribosomes, Name_mRNASubIdx); }
    //} 

    // Unconditional initialization of polymerase variables 
    Initialize_PolymeraseReaction_Matrix(ofs, PolymeraseTypes);
    Initialize_PolymeraseReaction_DNAP(ofs, DNAPs);
    Initialize_PolymeraseReaction_RNAP(ofs, RNAPs);
    Initialize_PolymeraseReaction_Ribosome(ofs, Ribosomes, Name_mRNASubIdx);



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

    float VolInit = 1;

    ofs << in+ "def Initialize(self):" << endl;
    ofs << in+ in+ "self.VolInit = " << VolInit << endl;
    ofs << in+ in+ "self.Vol = np.full([" << MatrixSize << ", 1], " << VolInit << ")" << endl;
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
    ofs << in+ in+ "MolarityFactor_Mol = np.where(MolarityFactor_Mol == 1, self.Vol[0], 1)" << endl;
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

    if (!PolymeraseList.empty()) {
        SetUp_PolymeraseReaction_Matrix(ofs, PolymeraseTypes);

        if (!DNAPs.empty()) { SetUp_PolymeraseReaction_DNAP(ofs, DNAPs); }
        if (!RNAPs.empty()) { SetUp_PolymeraseReaction_RNAP(ofs, RNAPs); }
        if (!Ribosomes.empty()) { SetUp_PolymeraseReaction_Ribosome(ofs, Ribosomes, Name_mRNASubIdx); }
    }   
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
    ofs << in+ in+ "Data[:, " << i_SimStep << ":" << i_Vol << "] = self.Vol[0]" << endl;

    int i_Count_Mol = i_Vol + MolNames.size();

    ofs << in+ in+ "Data[:, " << i_Vol << ":" << i_Count_Mol << "] = self.Count_All[0]" << endl;

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
    ofs << in+ in+ "Gene2UniprotID = dict()" << endl;
    ofs << in+ in+ "db_UniprotID = self.LoadTSVDatabase(r'./Database/UniprotID_Ecoli.tsv')" << endl;
    ofs << in+ in+ "for i, Value in enumerate(db_UniprotID):" << endl;
    ofs << in+ in+ in+ "Entry, Reviewed, Entry_Name, Protein_names, Gene_Names, Organism, Length = Value" << endl;
    ofs << in+ in+ in+ "Gene2UniprotID[Gene_Names.split(' ')[0]] = Entry" << endl;
    ofs << endl;

    ofs << in+ in+ "Gene2PDBID = dict()" << endl;
    ofs << in+ in+ "db_PDBID = self.LoadTSVDatabase(r'./Database/Gene2PDBID.tsv')" << endl;
    ofs << in+ in+ "for i, Value in enumerate(db_PDBID):" << endl;
    ofs << in+ in+ in+ "Symbol, PDBID = Value" << endl;
    //ofs << in+ in+ in+ "if not PDBID:" << endl;
    //ofs << in+ in+ in+ in+ "PDBID = ''" << endl;
    ofs << in+ in+ in+ "Gene2PDBID[Symbol] = PDBID" << endl;
    ofs << endl;

    ofs << in+ in+ "db = dict()" << endl;
    ofs << in+ in+ "NUniq_Genes = len(db_genes)" << endl;
    ofs << in+ in+ "db['Symbol'] = list()" << endl;
    ofs << in+ in+ "db['Length'] = np.zeros(NUniq_Genes)" << endl;
    ofs << in+ in+ "db['Coord'] = np.zeros(NUniq_Genes)" << endl;
    ofs << in+ in+ "db['Dir'] = np.zeros(NUniq_Genes)" << endl;
    ofs << in+ in+ "db['Seq'] = list()" << endl;
    ofs << in+ in+ "db['UniprotID'] = list()" << endl;
    ofs << in+ in+ "db['PDBID'] = list()" << endl;
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
    ofs << in+ in+ in+ "UniprotID = ''" << endl;
    ofs << in+ in+ in+ "if Symbol in Gene2UniprotID:" << endl;
    ofs << in+ in+ in+ in+ "UniprotID = Gene2UniprotID[Symbol]" << endl;
    ofs << in+ in+ in+ "db['UniprotID'].append(UniprotID)" << endl;
    ofs << endl;
    ofs << in+ in+ in+ "PDBID = ''" << endl;
    ofs << in+ in+ in+ "if Symbol in Gene2PDBID:" << endl;
    ofs << in+ in+ in+ in+ "PDBID = Gene2PDBID[Symbol]" << endl;
    ofs << in+ in+ in+ "db['PDBID'].append(PDBID)" << endl;
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

    std::cout << "  Simulation state has been generated: ";
}
