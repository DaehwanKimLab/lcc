#include <iostream>
#include <queue>
#include <cassert>
#include <algorithm>

#include <unistd.h>
#include <getopt.h>

#include "node.h"
#include "option.h"
#include "context.h"
#include "simulation.h"
#include "datamanager.h"
#include "util.h"

extern NBlock* ProgramBlock;
extern int yyparse();
extern int yylex_destroy();
extern FILE* yyin;

using namespace std;

FOption Option;
FCompilerContext Context;
FSimulation Simulation;
FDataset Dataset;
FState State;
FDataManager DataManager;

const char *VersionString = "1.0.0";
using namespace std;

void DumpNBlock(const NBlock* InProgramBlock) 
{
    if (!InProgramBlock) {
        return;
    }
    for (const auto& stmt : InProgramBlock->Statements) {
        if (stmt) {
            stmt->Print(std::cout);
            std::cout << std::endl;
        }
    }
}

void TraversalNode(NBlock* InProgramBlock)
{
    ostream& os = std::cout;
    FTraversalContext tc(std::cerr);
    tc.Queue.push(InProgramBlock);

    os << endl << "## TraversalNode ##" << endl;

    while(!tc.Queue.empty()) {
        const NNode* node = tc.Queue.front(); tc.Queue.pop();

        if (Utils::is_class_of<NProteinDeclaration, NNode>(node)) {
            auto Protein = dynamic_cast<const NProteinDeclaration *>(node);
            os << "Protein Id: " << Protein->Id.Name << endl;
            // Protein->Print(os);

            auto& Id = Protein->Id;
	    
            string Name = Id.Name;
            string Substrate = Context.QueryEnzymeTable(Name, "Substrate");
            float kcat = std::stof(Context.QueryEnzymeTable(Name, "kcat"));
            float kM = std::stof(Context.QueryEnzymeTable(Name, "kM"));

            os << "  Enzyme Query Results: " << Substrate << ", " << kcat << ", " << kM << endl;

            FEnzyme Enzyme(Name, Substrate, kcat, kM);

            auto& OverallReaction = Protein->OverallReaction;
            // os << "  OverallReaction:" << endl;
            map<string, int> Stoichiometry;

            int Coefficient;
            for (const auto& reactant : OverallReaction.Reactants) {
                Coefficient = -1; // update when coeff is fully implemented in parser
                os << "    Reactants: " << "(" << Coefficient << ")" << reactant->Name << ", " << endl;
                Stoichiometry[reactant->Name]= Coefficient;
            }

            for (const auto& product : OverallReaction.Products) {
                Coefficient = 1; // update when coeff is fully implemented in parser
                os << "    Products: " << "(" << Coefficient << ")" << product->Name << ", " << endl;
                Stoichiometry[product->Name]= Coefficient;
            }

            FEnzymaticReaction EnzymaticReaction(Name, Stoichiometry, Name);

//            if (Protein->Block) {
//                auto& Block = Protein->Block;
//                for (auto& stmt: Block->Statements) {
//                    os << "  "; stmt->Print(os);
//                }
//
//            }
//            Context.ProteinList.emplace_back(*Protein);

            Context.EnzymeList.emplace_back(Enzyme);
            Context.EnzymaticReactionList.emplace_back(EnzymaticReaction);
 
        } else if (Utils::is_class_of<NPathwayDeclaration, NNode>(node)) {
            auto Pathway = dynamic_cast<const NPathwayDeclaration *>(node);
            os << "Pathway: " << Pathway->Id.Name << endl;

            string Name = Pathway->Id.Name;
            vector<string> Sequence;

            os << "  Enzymes: ";
            if (Pathway->PathwayChainReaction) {
                auto& PathwayChainReaction = Pathway->PathwayChainReaction;
                auto& Exprs = PathwayChainReaction->Exprs;
                for (auto& expr: Exprs) {
//                    os << "  "; expr->Print(os);
                    auto& Identifiers = expr->Identifiers;
                    for (auto& Id: Identifiers) {
                        os << Id.Name << ", ";
                        Sequence.push_back(Id.Name);
                    }  
                }
            }
            os << endl;

            FPathway Pathway_New(Name, Sequence);
            Context.PathwayList.emplace_back(Pathway_New);

        } else if (Utils::is_class_of<NOrganismDeclaration, NNode>(node)) {
            auto Organism = dynamic_cast<const NOrganismDeclaration *>(node);
            os << "Organism: " << Organism->Id.Name << endl;
            os << "  " << Organism->Description << endl;

        } else if (Utils::is_class_of<NExperimentDeclaration, NNode>(node)) {
            auto Experiment = dynamic_cast<const NExperimentDeclaration *>(node);
            os << "Experiment: " << Experiment->Id.Name << endl;
            os << "  " << Experiment->Description << endl;
            if (Experiment->Block) {
                for(const auto& stmt : Experiment->Block->Statements) {
                    os << "  "; stmt->Print(os); os << endl;
                }
            }
        } else if (Utils::is_class_of<NUsingStatement, NNode>(node)) {
            auto UsingStmt = dynamic_cast<const NUsingStatement *>(node);
            os << "Using(" << UsingStmt->Type << "): " << UsingStmt->Id.Name << endl;
            Context.UsingModuleList.emplace_back(UsingStmt->Id.Name);
        }
#if 0
else if (Utils::is_class_of<NIdentifier, NNode>(node)) {
            auto Identifier = dynamic_cast<const NIdentifier *>(node);
            os << "Identifier: " << Identifier->Name  << endl;
            Context.IdentifierList.emplace_back(Identifier->Name);
        }
#endif

        node->Visit(tc);
    }
}

void ScanNodes(const NBlock* InProgramBlock)
{
    if (!InProgramBlock) {
        return;
    }
}


// Utils for string expression
std::string JoinStr2Str(std::vector<std::string> StringList)
{
   std::string JoinedStr;
   for (auto& Str : StringList){
       JoinedStr += "'" + Str + "'" + ", ";
   }
   return JoinedStr;
}

std::string JoinInt2Str(std::vector<int> IntList) 
{
   std::string JoinedStr;
   for (auto& Int : IntList){
       JoinedStr += std::to_string(Int) + ", ";
   }
   return JoinedStr;
}

std::string JoinInt2Str_Idx(std::vector<int> IntList) 
{
   std::string JoinedStr;
   for (auto& Int : IntList){
       JoinedStr += std::to_string(Int) + ", ";
   }
   return JoinedStr;
}

std::string JoinFloat2Str(std::vector<float> FloatList) 
{
   std::string JoinedStr;
   for (auto& Float : FloatList){
       JoinedStr += std::to_string(Float) + ", ";
   }
   return JoinedStr;
}

std::string Matrix2Str(std::vector<std::vector<int>> Matrix) 
{
    std::string MatrixStr;
    int N_Rows;
    int N_Columns;

    N_Rows = Matrix.size();
    for (auto& Row : Matrix) {
        MatrixStr += "[";
        N_Columns = Row.size();
        for (auto& Item : Row) {
            MatrixStr += std::to_string(Item) + ", ";
        }
        MatrixStr += "], ";
    }
    return MatrixStr;
}

void WriteSimModule(int TestInt)
{
    // write simulation.py
    std::ofstream ofs(Option.SimModuleFile.c_str());
    std::string endl = "\n";
    std::string in = "    ";

    // IMPORT
    ofs << "import os, sys" << endl;
    ofs << "import numpy as np" << endl;
    // ofs << "import tensorflow as tf" << endl;
    ofs << "from datetime import datetime" << endl;
    ofs << "import csv" << endl;
    ofs << endl;

    // C++ to Python Data conversion
    ofs << "PathwayList = list()" << endl;
    for (auto& pathway : Context.PathwayList) {
       ofs << "PathwayList.append({";
       ofs << "'Name' : '" << pathway.Name << "', ";
       ofs << "'Sequence' : [";
       for (auto& enzyme : pathway.Sequence) {
           ofs << "'" << enzyme << "', ";
       }
       ofs << "]" << "})" << endl;
    }
    ofs << endl;

    ofs << "EnzymeList = list()" << endl;
    for (auto& enzyme : Context.EnzymeList) {
        ofs << "EnzymeList.append({";
        ofs << "'Name' : '" << enzyme.Name << "', ";
        ofs << "'Substrate' : '" << enzyme.Substrate << "', ";
        ofs << "'kcat' : " << enzyme.kcat << ", ";
        ofs << "'kM' : " << enzyme.kM << ", })" << endl;
    }
    ofs << endl;

    // BODY
    ofs << "def main():   # add verbose" << endl;
    ofs << endl;

    // user input
    ofs << in+ "N_SimSteps = 100" << endl;
    

    ofs << endl;

    // utilities
    ofs << in+ "def ConcToCount(Conc_Molecule, Volume):" << endl;
    ofs << in+ in+ "return Conc_Molecule * Volume" << endl;
    ofs << endl;

    ofs << in+ "def CountToConc(Count_Molecule, Volume):" << endl;
    ofs << in+ in+ "return Count_Molecule / Volume" << endl;
    ofs << endl;

    ofs << in+ "def MichaelisMentenEqn(Conc_Enzyme, Conc_Substrate, kcat, kM):" << endl;
    ofs << in+ in+ "# Return Rate of Reaction" << endl;
    ofs << in+ in+ "return (kcat * Conc_Enzyme * Conc_Substrate) / (Conc_Substrate + kM)" << endl;

    ofs << in+ "def MichaelisMentenEqn_Array(Conc_Enzyme, Conc_Substrate, kcat, kM):" << endl;
    ofs << in+ in+ "# Return Rate of Reaction" << endl;
    ofs << in+ in+ "return np.multiply(np.multiply(kcat, Conc_Enzyme), Conc_Substrate) / (Conc_Substrate + kM)" << endl;

    ofs << endl;

    // class FState 
    ofs << in+ "class FState:" << endl;
    ofs << in+ in+ "def __init__(self):" << endl;
    ofs << in+ in+ in+ "self.Vol = 0" << endl;
    ofs << endl;

    std::vector<std::string> EnzymeNames = Context.GetNames_EnzymeList();
    std::vector<std::string> SubstrateNames = Context.GetSubstrateNames_EnzymaticReactionList();
    std::vector<float> kcats = Context.Getkcats_EnzymeList();
    std::vector<float> kMs = Context.GetkMs_EnzymeList();
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix();
    std::vector<int> Idx_EnzSubInAllSub = Context.GetEnzSubstrateIdxFromAllSubstrates();

    ofs << in+ in+ in+ "# State Arrays" << endl;
    ofs << in+ in+ in+ "self.Count_Enz = np.zeros(" << EnzymeNames.size() << ")" << endl;
    ofs << in+ in+ in+ "self.Count_Sub = np.zeros(" << SubstrateNames.size() << ")" << endl;
    ofs << in+ in+ in+ "self.dCount_Enz = np.zeros(len(self.Count_Enz))" << endl;
    ofs << in+ in+ in+ "self.dCount_Sub = np.zeros(len(self.Count_Sub))" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# K Constant Arrays" << endl;
    ofs << in+ in+ in+ "self.Const_kcats = np.zeros(len(self.Count_Enz))" << endl;
    ofs << in+ in+ in+ "self.Const_kMs = np.zeros(len(self.Count_Enz))" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Stoichiometry Matrix" << endl;
    ofs << in+ in+ in+ "self.Const_StoichMatrix = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Indices" << endl;
    ofs << in+ in+ in+ "self.Idx_EnzSubInAllSub = 0" << endl;

//    ofs << in+ in+ "def SetEnzCount(self, Name, Count):" << endl;
//    ofs << in+ in+ in+ "self.Enz2Count[Name] = Count" << endl;
//    ofs << endl;
//    ofs << in+ in+ "def SetSubCount(self, Name, Count):" << endl;
//    ofs << in+ in+ in+ "self.Sub2Count[Name] = Count" << endl;
//    ofs << endl;
//    ofs << in+ in+ "def GetEnzCount(self, Name):" << endl;
//    ofs << in+ in+ in+ "return self.Enz2Count[Name]" << endl;
//    ofs << endl;
//    ofs << in+ in+ "def GetSubCount(self, Name):" << endl;
//    ofs << in+ in+ in+ "return self.Sub2Count[Name]" << endl;
//    ofs << endl;
//    ofs << in+ in+ "def GetEnzConc(self, Name):" << endl;
//    ofs << in+ in+ in+ "return self.Enz2Count[Name] / self.Vol" << endl;
//    ofs << endl;
//    ofs << in+ in+ "def GetSubConc(self, Name):" << endl;
//    ofs << in+ in+ in+ "return self.Sub2Count[Name] / self.Vol" << endl;
//    ofs << endl;


    ofs << in+ in+ "def Initialize(self):" << endl;
    ofs << in+ in+ in+ "self.Vol = np.asmatrix([1])" << endl;
    ofs << in+ in+ in+ "self.Count_Enz = np.asmatrix(np.random.randint(10, high=100, size=len(self.Count_Enz)))" << endl;
    ofs << in+ in+ in+ "self.Count_Sub = np.asmatrix(np.random.randint(100, high=1000, size=len(self.Count_Sub)))" << endl;
    // ofs << in+ in+ in+ "self.Counts_Enz = np.matrix(np.array(5, size=len(self.Counts_Enz)))" << endl;
    // ofs << in+ in+ in+ "self.Counts_Sub = np.matrix(np.array(500, size=len(self.Counts_Sub)))" << endl;
    ofs << in+ in+ in+ "self.Const_kcats = np.asmatrix([" << JoinFloat2Str(kcats) << "])" << endl;
    ofs << in+ in+ in+ "self.Const_kMs = np.asmatrix([" << JoinFloat2Str(kMs) << "])" << endl;
    ofs << in+ in+ in+ "self.Const_StoichMatrix = np.asmatrix([" << Matrix2Str(StoichMatrix) << "])" << endl;

    ofs << in+ in+ in+ "self.Idx_EnzSubInAllSub = np.asmatrix([" << JoinInt2Str_Idx(Idx_EnzSubInAllSub) << "])" << endl;
    ofs << endl;

    ofs << in+ in+ "def ExportLegend(self):" << endl;
    ofs << in+ in+ in+ "return ['SimStep', 'Vol', " << JoinStr2Str(EnzymeNames) << JoinStr2Str(SubstrateNames) << "]" << endl;
    ofs << endl;

    ofs << in+ in+ "def ExportData(self, SimStep):" << endl;
    ofs << in+ in+ in+ "Data = np.asmatrix(np.zeros(2 + " << EnzymeNames.size() << " + " << SubstrateNames.size() << "))" << endl;
    int i = 0;
    int i_SimStep = i + 1;
    ofs << in+ in+ in+ "Data[0, " << i << ":" << i_SimStep << "] = SimStep" << endl;
    int i_Vol = i_SimStep + 1;
    ofs << in+ in+ in+ "Data[0, " << i_SimStep << ":" << i_Vol << "] = self.Vol" << endl;
    int i_Count_Enz = i_Vol + EnzymeNames.size();
    ofs << in+ in+ in+ "Data[0, " << i_Vol << ":" << i_Count_Enz << "] = self.Count_Enz" << endl;
    int i_Count_Sub = i_Count_Enz + SubstrateNames.size();
    ofs << in+ in+ in+ "Data[0, " << i_Count_Enz << ":" << i_Count_Sub << "] = self.Count_Sub" << endl;
    ofs << in+ in+ in+ "return Data" << endl;
    ofs << endl;
    
    // class FDataset
    ofs << in+ "class FDataset:" << endl;
    ofs << in+ in+ "def __init__(self):" << endl;
    ofs << in+ in+ in+ "self.Legend = list()" << endl;
    ofs << in+ in+ in+ "self.Data = 0" << endl;
    ofs << endl;

    ofs << in+ in+ "def PrintLegend(self):" << endl;
    ofs << in+ in+ in+ "print(self.Legend)" << endl;
    ofs << endl;

    ofs << in+ in+ "def PrintData(self):" << endl;
    ofs << in+ in+ in+ "print(self.Data)" << endl;
    ofs << endl;

    // class FSimulation
    ofs << in+ "class FSimulation:" << endl;
    ofs << in+ in+ "def __init__(self, InState, InDataset, InDM):" << endl;
    ofs << in+ in+ in+ "self.N_SimSteps = 0" << endl;
    ofs << in+ in+ in+ "self.SimStep = 0" << endl;
    ofs << in+ in+ in+ "self.State = InState" << endl;
    ofs << in+ in+ in+ "self.Dataset = InDataset" << endl;
    ofs << in+ in+ in+ "self.DM = InDM" << endl;
    ofs << endl;

    ofs << in+ in+ "def Initialize(self, InN_SimSteps):" << endl;
    ofs << in+ in+ in+ "print('Simulation Initialized...')" << endl;
    ofs << in+ in+ in+ "self.N_SimSteps = np.asmatrix([InN_SimSteps])" << endl;
    ofs << in+ in+ in+ "self.State.Initialize()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Legend Export" << endl;
    ofs << in+ in+ in+ "self.Dataset.Legend = self.State.ExportLegend()" << endl;
    ofs << in+ in+ in+ "self.DM.SetLegend(self.Dataset.Legend)" << endl;

    ofs << in+ in+ in+ "# Data Export" << endl;
    ofs << in+ in+ in+ "self.Dataset.Data = self.State.ExportData(self.SimStep)" << endl;
    ofs << in+ in+ in+ "self.DM.Add(self.Dataset.Data)" << endl;
    ofs << endl;

    ofs << in+ in+ "def Run(self):" << endl;
    ofs << in+ in+ in+ "print('Simulation Run Begins...')" << endl;
    ofs << in+ in+ in+ "while self.SimStep < self.N_SimSteps:" << endl;
    ofs << in+ in+ in+ in+ "self.SimStep += 1" << endl;

    ofs << in+ in+ in+ in+ "# Run Enzymatic Reactions" << endl;

    ofs << in+ in+ in+ in+ "Conc_Enz = self.State.Count_Enz / self.State.Vol" << endl;
    ofs << in+ in+ in+ in+ "Conc_EnzSub = np.take(self.State.Count_Sub, self.State.Idx_EnzSubInAllSub) / self.State.Vol" << endl;
    ofs << in+ in+ in+ in+ "Rate = MichaelisMentenEqn_Array(Conc_Enz, Conc_EnzSub, self.State.Const_kcats, self.State.Const_kMs)" << endl;
    ofs << in+ in+ in+ in+ "self.State.dCount_Sub = np.transpose(np.matmul(np.transpose(self.State.Const_StoichMatrix), np.transpose(Rate)))" << endl;

//    ofs << in+ in+ in+ in+ "# Display" << endl;
//    ofs << in+ in+ in+ in+ "EnzName = '" << EnzymeName << "' + '\\t|'" << endl;
//    ofs << in+ in+ in+ in+ "EnzStr = 'EnzConc: ' + str(EnzConc) + '\\t|'" << endl;
//    ofs << in+ in+ in+ in+ "SubStr = 'SubConc: ' + str(SubConc) + '\\t|'" << endl;
//    ofs << in+ in+ in+ in+ "RateStr = 'Rate: ' + str(Rate)" << endl;
//    ofs << in+ in+ in+ in+ "print(EnzName + EnzStr + SubStr + RateStr)" << endl;
// 
//    for (auto& EnzymaticReaction : Context.EnzymaticReactionList) {
//        if ((EnzymaticReaction.Enzyme == EnzymeName) & EnzymaticReaction.CheckIfReactant(Substrate)){
//            for (std::pair<std::string, int> Stoich : EnzymaticReaction.Stoichiometry) {
//                ofs << in+ in+ in+ in+ "Coeff = " << Stoich.second << endl;
//                ofs << in+ in+ in+ in+ "SubDelConc = Rate * Coeff" << endl;
//                ofs << in+ in+ in+ in+ "self.State.Sub2DelCount['" << Stoich.first << "'] = ConcToCount(SubDelConc, self.State.Vol)" << endl;
//                ofs << endl;
//            }
//        }
//    }
                           
    ofs << in+ in+ in+ in+ "# Update Substrate Count" << endl;
    ofs << in+ in+ in+ in+ "self.State.Count_Sub = self.State.Count_Sub + self.State.dCount_Sub" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ in+ "self.Dataset.Data = self.State.ExportData(self.SimStep)" << endl;
    ofs << in+ in+ in+ in+ "self.DM.Add(self.Dataset.Data)" << endl;
    

    ofs << endl;

    // class FDataManager
    ofs << in+ "class FDataManager:" << endl;
    ofs << in+ in+ "def __init__(self):" << endl;
    ofs << in+ in+ in+ "self.Legend = list()" << endl;
    ofs << in+ in+ in+ "self.DataBuffer = list()" << endl;
    ofs << endl;
    ofs << in+ in+ "def SetLegend(self, InLegend):" << endl;
    ofs << in+ in+ in+ "self.Legend = InLegend" << endl;
    ofs << endl;
    ofs << in+ in+ "def Add(self, InData):" << endl;
    ofs << in+ in+ in+ "self.DataBuffer.append(InData)" << endl;
    ofs << endl;
    ofs << in+ in+ "def SaveToFile(self, InFileName):" << endl;
    ofs << in+ in+ in+ "with open(InFileName, 'w') as OutFile:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter = csv.writer(OutFile, delimiter='\\t')" << endl;
    ofs << in+ in+ in+ in+ "if self.Legend:" << endl;
    ofs << in+ in+ in+ in+ in+ "TsvWriter.writerow(self.Legend)" << endl;
    ofs << in+ in+ in+ in+ "for Row in self.DataBuffer:" << endl;
    ofs << in+ in+ in+ in+ in+ "TsvWriter.writerow(Row)" << endl;
    ofs << endl;

    // Instantiate Objects
    ofs << in+ "State = FState()" << endl;
    ofs << in+ "Data = FDataset()" << endl;
    ofs << in+ "DM = FDataManager()" << endl;
    ofs << in+ "Sim = FSimulation(State, Data, DM)" << endl;
    ofs << endl;

    // Simulation Module
    ofs << in+ "Sim.Initialize(N_SimSteps)" << endl;
    ofs << in+ "Sim.Run()" << endl;
    ofs << endl;
    ofs << in+ "DM.SaveToFile('" << Option.SimResultFile.c_str() << "')" << endl;
    ofs << endl;

    // MAIN
    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "main()" << endl;
    ofs << endl;

}


int main(int argc, char *argv[])
{
    if (Option.Parse(argc, argv)) {
        Option.Usage(argv[0]);
        return -1;
    }
    if (Option.bVersion) {
        Option.ShowVersion(argv[0]);
        return 0;
    }
    if (Option.bShowHelp || Option.SourceFiles.empty()) {
        Option.Usage(argv[0]);
        return 0;
    }
    if (Option.Verbose) {
        Option.Dump();
    }

    ostream& os = std::cout;

    // Load genes.tsv
    if (!Option.bParseOnly)
    {
        std::cout<< endl << "## Loading Database ##" << std::endl;
        Context.Init(Option);
        vector<string> Keys;
        Keys.emplace_back("symbol");
        Keys.emplace_back("name");
        Keys.emplace_back("rnaId");
        Context.GeneTable.Dump(Keys);

		Keys.clear();
        Keys.emplace_back("reaction id");
        Keys.emplace_back("stoichiometry");
		Context.ReactionTable.Dump(Keys);

                Keys.clear();
        Keys.emplace_back("EnzymeName");
        Keys.emplace_back("Substrate");
        Keys.emplace_back("kcat");
        Keys.emplace_back("kM");
                Context.EnzymeTable.Dump(Keys);
    }

    for (const auto& SourceFile: Option.SourceFiles) {
        ProgramBlock = nullptr;
        yyin = fopen(SourceFile.c_str(), "r");
        if (yyin == NULL) {
            std::cerr << "Can't open file: " << SourceFile << std::endl;
            continue;
        }

        yyparse();
        fclose(yyin);
        yylex_destroy();

        if (Option.bDebug) {
            std::cout << ProgramBlock << std::endl;
            DumpNBlock(ProgramBlock);
        }

        TraversalNode(ProgramBlock);

        Context.PrintLists(os);

        delete ProgramBlock;
    }


    if (Option.bDebug) {
        for(const auto& protein : Context.ProteinList) {
            cout << protein.GetName() << endl;
        }
        for(const auto& name : Context.UsingModuleList) {
            cout << name << endl;
        }
    }

    if (Option.bSimPython) {
        cout << endl << "## Simulation_Python ##" << endl;
        
        int TestInt = 1;
        WriteSimModule(TestInt);

        cout << "Simulation_Python module has been generated: ";
        cout << Option.SimModuleFile << endl;
    }

    if (Option.bSimCpp) {
        // temporary C++ simulation code
        cout << endl << "## Simulation_C++ ##" << endl;
  
        Simulation.Init(State, Dataset, DataManager, 100);
        Simulation.Run(State, Context, Dataset);
    }
            if (!Option.SimResultFile.empty()) {
		DataManager.SaveToFile(Option.SimResultFile.c_str());

                
 	    }       

    return 0;
}

