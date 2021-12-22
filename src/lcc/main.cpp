#include <iostream>
#include <queue>
#include <cassert>
#include <algorithm>
#include <new>
#include <string>
#include <vector>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#include <getopt.h>
#endif

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

        // This is inteded for NEnzymeDeclaration, to be fixed later on.
        if (Utils::is_class_of<NProteinDeclaration, NNode>(node)) {
            auto N_Enzyme = dynamic_cast<const NProteinDeclaration *>(node);
            os << "Enzyme Id: " << N_Enzyme->Id.Name << endl;
            // Enzyme->Print(os);

            auto& Id = N_Enzyme->Id;	    
	    
            string Name = Id.Name;
            string Substrate = Context.QueryTable(Name, "Substrate", Context.EnzymeTable);
            float kcat = std::stof(Context.QueryTable(Name, "kcat", Context.EnzymeTable));
            float kM = std::stof(Context.QueryTable(Name, "kM", Context.EnzymeTable));

            FEnzyme * Enzyme = new FEnzyme(Name, Substrate, kcat, kM);
            Enzyme->Print(os);
            Context.AddToMoleculeList(Enzyme);

            auto& OverallReaction = N_Enzyme->OverallReaction;
            // os << "  OverallReaction:" << endl;
            map<string, int> Stoichiometry;
			string Location = OverallReaction.Location.Name;

            int Coefficient;
            for (const auto& reactant : OverallReaction.Reactants) {
                Coefficient = -1; // update when coeff is fully implemented in parser
                os << "    Reactants: " << "(" << Coefficient << ")" << reactant->Name << ", " << endl;
                Stoichiometry[reactant->Name]= Coefficient;

                FSmallMolecule * Molecule = new FSmallMolecule(reactant->Name);
                Molecule->Print(os);
                Context.AddToMoleculeList(Molecule);
            }

            for (const auto& product : OverallReaction.Products) {
                Coefficient = 1; // update when coeff is fully implemented in parser
                os << "    Products: " << "(" << Coefficient << ")" << product->Name << ", " << endl;
                Stoichiometry[product->Name]= Coefficient;

                FSmallMolecule * Molecule = new FSmallMolecule(product->Name);
                Molecule->Print(os);
                Context.AddToMoleculeList(Molecule);
            }

			if (!Location.empty()) {
                os << "    Location: " << Location << endl;
			}

            FEnzymaticReaction *EnzymaticReaction = new FEnzymaticReaction(Name, Stoichiometry, Name);
            EnzymaticReaction->Print(os);
            Context.AddToReactionList(EnzymaticReaction);

//            if (N_Enzyme->Block) {
//                auto& Block = N_Enzyme->Block;
//                for (auto& stmt: Block->Statements) {
//                    os << "  "; stmt->Print(os);
//                }
//
//            }
 
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

            FPathway Pathway_New(Name, Sequence); // Fixme
            Context.PathwayList.emplace_back(Pathway_New);

        } else if (Utils::is_class_of<NPolymeraseDeclaration, NNode>(node)) {
            auto N_Polymerase = dynamic_cast<const NPolymeraseDeclaration *>(node);
            os << "Polymerase Id: " << N_Polymerase->Id.Name << endl;
            // N_Polymerase->Print(os);

            auto& Id = N_Polymerase->Id;
	    
            string Name = Id.Name;
            string Substrate = Context.QueryTable(Name, "Substrate", Context.PolymeraseTable);
            float Rate = std::stof(Context.QueryTable(Name, "Rate", Context.PolymeraseTable));

            FPolymerase * Polymerase = new FPolymerase(Name, Substrate, Rate);
            Polymerase->Print(os);
            Context.AddToMoleculeList(Polymerase);


#if 1
            for (const shared_ptr<NStatement>& stmt: N_Polymerase->Statements) {
                if (Utils::is_class_of<NElongationStatement, NStatement>(stmt.get())) {
                    const shared_ptr<NElongationStatement> elongstmt = dynamic_pointer_cast<NElongationStatement>(stmt);
                    os << "---This is an elongation statement of the polymerase stmt---" << endl;
                    elongstmt->Print(os);
                    os << "-----------------" << endl;
                } else if (Utils::is_class_of<NInitiationStatement, NStatement>(stmt.get())) {
                    const shared_ptr<NInitiationStatement> initstmt = dynamic_pointer_cast<NInitiationStatement>(stmt);
                    os << "---This is an initiation statement of the polymerase stmt---" << endl;
                    initstmt->Print(os);
                    os << "-----------------" << endl;
                } else if (Utils::is_class_of<NTerminationStatement, NStatement>(stmt.get())) {
                    const shared_ptr<NTerminationStatement> termstmt = dynamic_pointer_cast<NTerminationStatement>(stmt);
                    os << "---This is a termination statement of the polymerase stmt---" << endl;
                    termstmt->Print(os);
                    os << "-----------------" << endl;
                }
            }



#endif

//            int i = 0;
//
//            for (const auto stmt : N_Polymerase->Statements) {
//                stmt->Print(os);
//                NStatement Statement = *stmt;
//
//                if (i == 1) {
//                    auto Elongation = static_cast<NElongationStatement *>(&Statement);
//                    NReaction ElongationReaction = Elongation->Reaction;
//                    ElongationReaction.Print(os);
//                    os << "AAAAAAAAAAA" << endl;
//
////                if (Utils::is_class_of<const NElongationStatement, const NStatement>(&Statement)) {
////                    auto Elongation = static_cast<const NElongationStatement *>(&Statement);
////                    NReaction ElongationReaction = Elongation->Reaction;
//
//                    os << "11111" << endl; 
//
//                    os << "  Elongation:"; ElongationReaction.Print(os);
//                    map<string, int> Stoichiometry;
//        			string Location = ElongationReaction.Location.Name;
//                    int Coefficient;
//                    std::vector<std::string> BuildingBlocks;
//
//                    for (const auto& reactant : ElongationReaction.Reactants) {
//                        Coefficient = -1; // update when coeff is fully implemented in parser
//                        os << "    Reactants: " << "(" << Coefficient << ")" << reactant->Name << ", " << endl;
//                        if ((reactant->Name == "dna_{n}") | (reactant->Name == "rna_{n}") | (reactant->Name == "peptide_{n}")) {
//                            continue;
//                        } else if (reactant->Name == "dnt") {
//                            BuildingBlocks = {"dATP", "dCTP", "dGTP", "dUTP"};
//                            continue;
//                        } else if (reactant->Name == "nt") {
//                            BuildingBlocks = {"ATP", "CTP", "GTP", "UTP"};
//                            continue;
//                        } else if (reactant->Name == "nt") {
//                            BuildingBlocks = {"ALPHA-ALANINE", "ARG", "ASN", "L-ASPARTATE", "CYS", "GLT", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "L-SELENOCYSTEINE", "VAL"};
//                            continue;
//                        }
//        
//                        Stoichiometry[reactant->Name]= Coefficient;
//        
//                        FSmallMolecule * Molecule = new FSmallMolecule(reactant->Name);
//                        Molecule->Print(os);
//                        Context.AddToMoleculeList(Molecule);
//                    }
//
//                    os << "2222" << endl;
//                    for (const auto& product : ElongationReaction.Products) {
//                        Coefficient = 1; // update when coeff is fully implemented in parser
//                        os << "    Products: " << "(" << Coefficient << ")" << product->Name << ", " << endl;
//                        if (product->Name == "rna_{n+1}") {
//                            continue;
//                        }
//                        Stoichiometry[product->Name]= Coefficient;
//
//                        FSmallMolecule * Molecule = new FSmallMolecule(product->Name);
//                        Molecule->Print(os);
//                        Context.AddToMoleculeList(Molecule);
//                    }
//        
//        			if (!Location.empty()) {
//                        os << "    Location: " << Location << endl;
//        			}
//
//                    os << "3333" << endl;
//                    for (auto& BuildingBlock : BuildingBlocks) {
//                        FSmallMolecule * Molecule = new FSmallMolecule(BuildingBlock);
//                        Molecule->Print(os);
//                        Context.AddToMoleculeList(Molecule);
//                    }
//
//                    FPolymeraseReaction *PolymeraseReaction = new FPolymeraseReaction(Name, Stoichiometry, Name, BuildingBlocks);
//                    PolymeraseReaction->Print(os);
//                    Context.AddToReactionList(PolymeraseReaction);
//                } // if
//            i++;
//            }





        } else if (Utils::is_class_of<NElongationStatement, NNode>(node)) {
            auto N_Elongation = dynamic_cast<const NElongationStatement *>(node);
            std::string Name = "rnap";

            auto& ElongationReaction = N_Elongation->Reaction;

            os << "  Elongation:"; ElongationReaction.Print(os);
            map<string, int> Stoichiometry;
			string Location = ElongationReaction.Location.Name;
            int Coefficient;
            std::vector<std::string> BuildingBlocks;

            for (const auto& reactant : ElongationReaction.Reactants) {
                Coefficient = -1; // update when coeff is fully implemented in parser
                os << "    Reactants: " << "(" << Coefficient << ")" << reactant->Name << ", " << endl;
                if ((reactant->Name == "dna_{n}") | (reactant->Name == "rna_{n}") | (reactant->Name == "peptide_{n}")) {
                    continue;
                } else if (reactant->Name == "dnt") {
                    BuildingBlocks = {"dATP", "dCTP", "dGTP", "dUTP"};
                    continue;
                } else if (reactant->Name == "nt") {
                    BuildingBlocks = {"ATP", "CTP", "GTP", "UTP"};
                    continue;
                } else if (reactant->Name == "nt") {
                    BuildingBlocks = {"ALPHA-ALANINE", "ARG", "ASN", "L-ASPARTATE", "CYS", "GLT", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "L-SELENOCYSTEINE", "VAL"};
                    continue;
                }

                Stoichiometry[reactant->Name]= Coefficient;

                FSmallMolecule * Molecule = new FSmallMolecule(reactant->Name);
                Molecule->Print(os);
                Context.AddToMoleculeList(Molecule);
            }

            for (const auto& product : ElongationReaction.Products) {
                Coefficient = 1; // update when coeff is fully implemented in parser
                os << "    Products: " << "(" << Coefficient << ")" << product->Name << ", " << endl;
                if (product->Name == "rna_{n+1}") {
                    continue;
                }
                Stoichiometry[product->Name]= Coefficient;

                FSmallMolecule * Molecule = new FSmallMolecule(product->Name);
                Molecule->Print(os);
                Context.AddToMoleculeList(Molecule);
            }

			if (!Location.empty()) {
                os << "    Location: " << Location << endl;
			}

            for (auto& BuildingBlock : BuildingBlocks) {
                os << "BuildingBlock" << endl;
                FSmallMolecule * Molecule = new FSmallMolecule(BuildingBlock);
                Molecule->Print(os);
                Context.AddToMoleculeList(Molecule);
            }

            FPolymeraseReaction *PolymeraseReaction = new FPolymeraseReaction(Name, Stoichiometry, Name, BuildingBlocks);
            PolymeraseReaction->Print(os);
            Context.AddToReactionList(PolymeraseReaction);
    
       






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

void WriteSimModule()
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
    ofs << in+ in+ "return (kcat * Conc_Enzyme * Conc_Substrate) / (Conc_Substrate + kM)" << endl;
    ofs << endl;

    ofs << in+ "def MichaelisMentenEqn_Array(Conc_Enzyme, Conc_Substrate, kcat, kM):" << endl;
    ofs << in+ in+ "return np.multiply(np.multiply(kcat, Conc_Enzyme), Conc_Substrate) / (Conc_Substrate + kM)" << endl;
    ofs << endl;

    // Elementary simulation functions
    ofs << in+ "def DetermineAmountOfBuildingBlocks(Freq, Rate):" << endl;
    ofs << in+ in+ "return np.transpose(np.multiply(np.transpose(Freq), np.transpose(Rate)))" << endl;
    ofs << endl;

    // class FState 
    ofs << in+ "class FState:" << endl;
    ofs << in+ in+ "def __init__(self):" << endl;
    ofs << in+ in+ in+ "self.Vol = 0" << endl;
    ofs << endl;

    // for enzyme reactions
    std::vector<const FEnzyme *> EnzymeList = Context.GetList_Enzyme_MoleculeList();

    std::vector<std::string> EnzNames = Context.GetNames_EnzymeList(EnzymeList);
    std::cout << endl;

    std::vector<const FEnzymaticReaction *> EnzymaticReactionList = Context.GetList_Enzymatic_ReactionList();

//    std::vector<const FSmallMolecule *> SMolList = Context.GetList_SmallMolecule_MoleculeList();

//    std::vector<std::string> SMolNames = Context.GetNames_SmallMoleculeList(SMolList);
//    std::cout << in+ in+ in+ "# SMolNames: " << JoinStr2Str(SMolNames) << endl;
//    os << endl;

    std::vector<float> kcats = Context.Getkcats_EnzymeList(EnzymeList);
    std::vector<float> kMs = Context.GetkMs_EnzymeList(EnzymeList);
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix_EnzymaticReaction(EnzymaticReactionList);

    std::vector<int> Idx_Enz = Context.GetIdx_Enzyme_MoleculeList();
    std::vector<int> Idx_EnzSub = Context.GetIdx_EnzymeSubstrate_MoleculeList();
    std::vector<int> Idx_SMol = Context.GetIdx_SmallMolecule_MoleculeList();
  
    //if (Utils::is_class_of<NProteinDeclaration, NNode>(node)) {
    //    auto Protein = dynamic_cast<const NProteinDeclaration *>(node);

    ofs << in+ in+ in+ "# State Arrays" << endl;
    ofs << in+ in+ in+ "self.Count_All = np.zeros([1, " << Context.MoleculeList.size() << "])" << endl;
    ofs << in+ in+ in+ "self.dCount_All = np.zeros([1, " << Context.MoleculeList.size() << "])" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Idx_Enz = 0" << endl;
    ofs << in+ in+ in+ "self.Idx_EnzSub = 0" << endl;
    ofs << in+ in+ in+ "self.Idx_SMol = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# K Constant Arrays" << endl;
    ofs << in+ in+ in+ "self.Const_kcats = 0" << endl;
    ofs << in+ in+ in+ "self.Const_kMs = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Stoichiometry Matrix" << endl;
    ofs << in+ in+ in+ "self.Const_StoichMatrix = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Indices" << endl;
    ofs << in+ in+ in+ "self.Idx_EnzSubInAllSub = 0" << endl;

    ofs << in+ in+ in+ "# Frequency Matrix" << endl;
    ofs << in+ in+ in+ "self.Freq_ReplicatingStrands = 0" << endl;
    ofs << in+ in+ in+ "self.Freq_RNAs = 0" << endl;
    ofs << in+ in+ in+ "self.Freq_Proteins = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.NUniq_ReplicatingStrands = 0" << endl;
    ofs << in+ in+ in+ "self.NUniq_RNAs = 0" << endl;
    ofs << in+ in+ in+ "self.NUniq_Proteins = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Len_ReplicatingStrands = 0" << endl;
    ofs << in+ in+ in+ "self.Len_NascentRNAs = 0" << endl;
    ofs << in+ in+ in+ "self.Len_NascentProteins = 0" << endl;

    ofs << in+ in+ in+ "self.Max_ReplicatingStrands = 0" << endl;
    ofs << in+ in+ in+ "self.Max_NascentRNAs = 0" << endl;
    ofs << in+ in+ in+ "self.Max_NascentProteins = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Rate_Replication = 0" << endl;
    ofs << in+ in+ in+ "self.Rate_Transcription = 0" << endl;
    ofs << in+ in+ in+ "self.Rate_Translation = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Weight_Replication = 0" << endl;
    ofs << in+ in+ in+ "self.Weight_Transcription = 0" << endl;
    ofs << in+ in+ in+ "self.Weight_Translation = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Count_Sub_Replication = 0" << endl;
    ofs << in+ in+ in+ "self.Count_Sub_Transcription = 0" << endl;
    ofs << in+ in+ in+ "self.Count_Sub_Translation = 0" << endl;
    ofs << endl;

//    ofs << in+ in+ in+ "self.Count_BuildingBlocks_Replication = 0" << endl;
//    ofs << in+ in+ in+ "self.Count_BuildingBlocks_Transcription = 0" << endl;
//    ofs << in+ in+ in+ "self.Count_BuildingBlocks_Translation = 0" << endl;
//    ofs << endl;

    ofs << in+ in+ in+ "self.dCount_Sub_Replication = 0" << endl;
    ofs << in+ in+ in+ "self.dCount_Sub_Transcription = 0" << endl;
    ofs << in+ in+ in+ "self.dCount_Sub_Translation = 0" << endl;
    ofs << endl;

//    ofs << in+ in+ in+ "self.dCount_BuildingBlocks_Replication = 0" << endl;
//    ofs << in+ in+ in+ "self.dCount_BuildingBlocks_Transcription = 0" << endl;
//    ofs << in+ in+ in+ "self.dCount_BuildingBlocks_Translation = 0" << endl;
//    ofs << endl;


    ofs << in+ in+ "def Initialize(self):" << endl;
    ofs << in+ in+ in+ "self.Vol = np.asmatrix([1])" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Enzyme Reactions (small molecules at ~mM range)" << endl;
    // ofs << in+ in+ in+ "self.Counts_Enz = np.matrix(np.array(5, size=len(self.Counts_Enz)))" << endl;
    // ofs << in+ in+ in+ "self.Counts_Sub = np.matrix(np.array(500, size=len(self.Counts_Sub)))" << endl;
    ofs << in+ in+ in+ "self.Const_kcats = np.asmatrix([" << JoinFloat2Str(kcats) << "])" << endl;
    ofs << in+ in+ in+ "self.Const_kMs = np.asmatrix([" << JoinFloat2Str(kMs) << "])" << endl;
    ofs << in+ in+ in+ "self.Const_StoichMatrix = np.asmatrix([" << Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Idx_Enz = np.asmatrix([" << JoinInt2Str_Idx(Idx_Enz) << "])" << endl;
    ofs << in+ in+ in+ "self.Idx_EnzSub = np.asmatrix([" << JoinInt2Str_Idx(Idx_EnzSub) << "])" << endl;
    ofs << in+ in+ in+ "self.Idx_SMol = np.asmatrix([" << JoinInt2Str_Idx(Idx_SMol) << "])" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "Count_Enz = np.asmatrix(np.random.randint(10, high=100, size=self.Idx_Enz.size))" << endl;
    ofs << in+ in+ in+ "Count_SMol = np.asmatrix(np.random.randint(100, high=1000, size=self.Idx_SMol.size))" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Enz, Count_Enz, axis=1)" << endl;
    ofs << in+ in+ in+ "np.put_along_axis(self.Count_All, self.Idx_SMol, Count_SMol, axis=1)" << endl;

    ofs << endl;


    ofs << in+ in+ in+ "# Elongation Reaction" << endl;
//    ofs << in+ in+ in+ "self.Freq_Chromosomes = " << STR() Context.GetFreqMatrixForChromosomes << endl;
//    ofs << in+ in+ in+ "self.Freq_RNAs = " << STR() Context.GetFreqMatrixForRNAs << endl;
//    ofs << in+ in+ in+ "self.Freq_Proteins = " << STR() Context.GetFreqMatrixForProteins.npy')" << endl;
    ofs << in+ in+ in+ "self.Freq_ReplicationStrands = np.asmatrix(np.load(r'./Database/Freq_NTsInChromosomes.npy'))" << endl;
    ofs << in+ in+ in+ "self.Freq_RNAs = np.asmatrix(np.load(r'./Database/Freq_NTsInRNAs.npy'))" << endl;
    ofs << in+ in+ in+ "self.Freq_Proteins = np.asmatrix(np.load(r'./Database/Freq_AAsInProteins.npy'))" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.NUniq_ReplicatingStrands = len(self.Freq_ReplicationStrands)" << endl;
    ofs << in+ in+ in+ "self.NUniq_RNAs = len(self.Freq_RNAs)" << endl;
    ofs << in+ in+ in+ "self.NUniq_Proteins = len(self.Freq_Proteins)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Len_ReplicatingStrands = np.asmatrix(np.full([4, len(self.Freq_ReplicationStrands)], 0))" << endl;
    ofs << in+ in+ in+ "self.Len_NascentRNAs = np.asmatrix(np.full([10, len(self.Freq_RNAs)], 0))" << endl;
    ofs << in+ in+ in+ "self.Len_NascentProteins = np.asmatrix(np.full([10, len(self.Freq_Proteins)], 0))" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Max_ReplicatingStrands = np.asmatrix(np.full(len(self.Freq_ReplicationStrands), np.load(r'./Database/Len_ChromosomesInGenome.npy') / 2))" << endl;
    ofs << in+ in+ in+ "self.Max_NascentRNAs = np.asmatrix(np.load(r'./Database/Len_RNAs.npy'))" << endl;
    ofs << in+ in+ in+ "self.Max_NascentProteins = np.asmatrix(np.load(r'./Database/Len_Proteins.npy'))" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Rate_Replication = 1000" << endl;
    ofs << in+ in+ in+ "self.Rate_Transcription = 60" << endl;
    ofs << in+ in+ in+ "self.Rate_Translation = 20" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Weight_Replication = 1" << endl;
    ofs << in+ in+ in+ "self.Weight_Transcription = 1" << endl;
    ofs << in+ in+ in+ "self.Weight_Translation = 1" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Count_Sub_Replication = np.asmatrix(np.full(len(self.Freq_ReplicationStrands), 100000))" << endl;
    ofs << in+ in+ in+ "self.Count_Sub_Transcription = np.asmatrix(np.full(len(self.Freq_RNAs), 100000))" << endl;
    ofs << in+ in+ in+ "self.Count_Sub_Translation = np.asmatrix(np.full(len(self.Freq_Proteins), 100000))" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Count_BuildingBlocks_Replication = np.random.randint(400000, high=500000, size=self.Freq_ReplicationStrands.shape[1])" << endl;
    ofs << in+ in+ in+ "self.Count_BuildingBlocks_Transcription = np.random.randint(400000, high=500000, size=self.Freq_RNAs.shape[1])" << endl;
    ofs << in+ in+ in+ "self.Count_BuildingBlocks_Translation = np.random.randint(400000, high=500000, size=self.Freq_Proteins.shape[1])" << endl;
    ofs << endl; 


    ofs << in+ in+ "def ExportLegend(self):" << endl;
    // for legends
    std::vector<std::string> MolNames = Context.GetNames_MoleculeList();
    // std::cout << "Legend_MoleNames: " << JoinStr2Str(MolNames) << endl;
    ofs << in+ in+ in+ "return ['SimStep', 'Vol', " << JoinStr2Str(MolNames) << "]" << endl;
    ofs << endl;

    ofs << in+ in+ "def ExportData(self, SimStep):" << endl;
    ofs << in+ in+ in+ "Data = np.asmatrix(np.zeros(2 + " << MolNames.size() << "))" << endl;
    int i = 0;
    int i_SimStep = i + 1;
    ofs << in+ in+ in+ "Data[0, " << i << ":" << i_SimStep << "] = SimStep" << endl;

    int i_Vol = i_SimStep + 1;
    ofs << in+ in+ in+ "Data[0, " << i_SimStep << ":" << i_Vol << "] = self.Vol" << endl;

    int i_Count_Mol = i_Vol + MolNames.size();
    ofs << in+ in+ in+ "Data[0, " << i_Vol << ":" << i_Count_Mol << "] = self.Count_All" << endl;

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
    ofs << in+ in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ in+ "def Run(self):" << endl;
    ofs << in+ in+ in+ "print('Simulation Run Begins...')" << endl;
    ofs << endl;
    ofs << in+ in+ in+ "while self.SimStep < self.N_SimSteps:" << endl;
    ofs << in+ in+ in+ in+ "self.SimStep += 1" << endl;

    ofs << in+ in+ in+ in+ "# Reset Substrate Count" << endl;
    ofs << in+ in+ in+ in+ "self.State.dCount_All = np.zeros_like(self.State.dCount_All)" << endl;
    
    ofs << in+ in+ in+ in+ "# Run Reactions" << endl;
    ofs << in+ in+ in+ in+ "self.EnzymaticReactions()" << endl;
    ofs << in+ in+ in+ in+ "self.ElongationReactions()" << endl;
                          
    ofs << in+ in+ in+ in+ "# Update Substrate Count" << endl;
    ofs << in+ in+ in+ in+ "self.State.Count_All += self.State.dCount_All" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "print('Simulation Run Completed')" << endl;
    ofs << endl;

    ofs << in+ in+ "def ExportData(self):" << endl;    
    ofs << in+ in+ in+ in+ "self.Dataset.Data = self.State.ExportData(self.SimStep)" << endl;
    ofs << in+ in+ in+ in+ "self.DM.Add(self.Dataset.Data)" << endl;
    ofs << endl;

    ofs << in+ in+ "def EnzymaticReactions(self):" << endl;
    ofs << in+ in+ in+ "Conc_Enz = np.take(self.State.Count_All, self.State.Idx_Enz) / self.State.Vol" << endl;
    ofs << in+ in+ in+ "Conc_EnzSub = np.take(self.State.Count_All, self.State.Idx_EnzSub) / self.State.Vol" << endl;
    ofs << in+ in+ in+ "Rate = MichaelisMentenEqn_Array(Conc_Enz, Conc_EnzSub, self.State.Const_kcats, self.State.Const_kMs)" << endl;
    // Update with mole indexes from EnzReactions
    ofs << in+ in+ in+ "dCount_SMol = np.zeros_like(self.State.dCount_All)" << endl;
    ofs << in+ in+ in+ "np.put_along_axis(dCount_SMol, self.State.Idx_SMol, np.transpose(np.matmul(np.transpose(self.State.Const_StoichMatrix), np.transpose(Rate))), axis=1)" << endl;
    ofs << in+ in+ in+ "self.State.dCount_All += dCount_SMol" << endl;
    ofs << endl;

    ofs << in+ in+ "def ElongationReactions(self):" << endl;
    ofs << in+ in+ in+ "# Replication" << endl;
    ofs << in+ in+ in+ "self.State.Len_ReplicatingStrands = self.Elongation(self.State.Len_ReplicatingStrands, self.State.Max_ReplicatingStrands, self.State.Rate_Replication, self.State.Weight_Replication, self.State.Freq_ReplicationStrands)" << endl;
    ofs << in+ in+ in+ "# Transcription" << endl;
    ofs << in+ in+ in+ "self.State.Len_NascentRNAs = self.Elongation(self.State.Len_NascentRNAs, self.State.Max_NascentRNAs, self.State.Rate_Transcription, self.State.Weight_Transcription, self.State.Freq_RNAs)" << endl;
    ofs << in+ in+ in+ "# Translation" << endl;
    ofs << in+ in+ in+ "self.State.Len_NascentProteins = self.Elongation(self.State.Len_NascentProteins, self.State.Max_NascentProteins, self.State.Rate_Translation, self.State.Weight_Translation, self.State.Freq_Proteins)" << endl;
    ofs << in+ in+ in+ "" << endl;
    ofs << endl;
//    ofs << in+ in+ in+ "NUniq_RNAs = len(self.Freq_RNAs)" << endl;
//    ofs << endl;


    ofs << in+ in+ "def Elongation(self, Len, Max, Rate, Weight, Freq):   # should additionally take buildingblock indices for count deduction" << endl;
    ofs << in+ in+ in+ "NUniq_BuildingBlocks = Freq.shape[1]" << endl;
    ofs << in+ in+ in+ "NUniq_Species = Freq.shape[0]" << endl;

//    ofs << in+ in+ in+ "Count_BuildingBlocks = self.State.Count_BuildingBlocks_Transcription" << endl;
//    ofs << in+ in+ in+ "dCount_BuildingBlocks = self.State.dCount_BuildingBlocks_Transcription" << endl;
    ofs << endl;

//    ofs << in+ in+ in+ "dLength = np.matmul(SMatrix,Rate )
    ofs << in+ in+ in+ "dLength = Rate   # this is not necessarily true based on the reaction input" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# TODO: Initiation " << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# def Elongation(self, ):" << endl;
    ofs << in+ in+ in+ "Len_Elongated = np.where(Len >= 0, Len + dLength, -1)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# def OverElongation Correction(self, ):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ in+ "Len_Over = np.where(Len_Elongated > Max, Len_Elongated - Max, 0)" << endl;
    ofs << in+ in+ in+ "N_OverElongated = np.sum(Len_Over, axis=0)" << endl;
    ofs << in+ in+ in+ "Len_Trimmed = Len_Elongated - Len_Over" << endl;
    ofs << in+ in+ in+ "N_Elongated_PerSpecies = np.sum(Len_Trimmed - Len, axis=0)" << endl;
    ofs << in+ in+ in+ "N_Elongated = np.sum(N_Elongated_PerSpecies)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# def Building block Consumption(self, ):" << endl;
    ofs << in+ in+ in+ "Raw = np.transpose(np.matmul(np.transpose(Freq), np.transpose(N_Elongated_PerSpecies)))" << endl;
    ofs << in+ in+ in+ "Rounded = np.around(Raw)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Discrepancy handling" << endl;
    ofs << in+ in+ in+ "Discrepancy = np.sum(Rounded) - N_Elongated" << endl;
    ofs << in+ in+ in+ "Sets, Remainder = np.divmod(Discrepancy, NUniq_BuildingBlocks)" << endl;
    ofs << in+ in+ in+ "Revised = Rounded + np.multiply(np.ones(NUniq_BuildingBlocks), np.int32(Sets)) + np.concatenate((np.ones(np.int32(Remainder)), np.zeros(np.int32(NUniq_BuildingBlocks - Remainder))))" << endl;
    ofs << endl;
    
    ofs << in+ in+ in+ "# Update dCount" << endl;
    ofs << in+ in+ in+ "# self.State.dCount_BuildingBlocks_Transcription = Revised" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# def Elongation Completion(self, ):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ in+ "N_Completed = np.sum(np.where(Len_Trimmed == Max, 1, 0))" << endl;
    ofs << in+ in+ in+ "Len_Completed =  np.where(Len_Trimmed == Max, -1, Len_Trimmed)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Export Data" << endl;
    ofs << in+ in+ in+ "# Rate" << endl;
    ofs << in+ in+ in+ "# N_Initiated" << endl;
    ofs << in+ in+ in+ "# N_Completed" << endl;
    
    ofs << in+ in+ in+ "return Len_Completed" << endl;
   
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
    ofs << in+ in+ in+ "with open(InFileName, 'w', newline='', encoding='utf-8') as OutFile:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter = csv.writer(OutFile, delimiter='\\t')" << endl;
    ofs << in+ in+ in+ in+ "if self.Legend:" << endl;
    ofs << in+ in+ in+ in+ in+ "TsvWriter.writerow(self.Legend)" << endl;
    ofs << in+ in+ in+ in+ "for Row in self.DataBuffer:" << endl;
    ofs << in+ in+ in+ in+ in+ "TsvWriter.writerow(np.array(Row).flatten().tolist())" << endl;
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
        Keys.emplace_back("Name");
        Keys.emplace_back("Substrate");
        Keys.emplace_back("kcat");
        Keys.emplace_back("kM");
                os << "# EnzymeTable #" << endl;
                Context.EnzymeTable.Dump(Keys);

                Keys.clear();
        Keys.emplace_back("Name");
        Keys.emplace_back("Substrate");
        Keys.emplace_back("Rate");
                os << "# PolymeraseTable #" << endl;
                Context.PolymeraseTable.Dump(Keys);
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
    
        for(const auto& name : Context.UsingModuleList) {
            cout << name << endl;
        }
    }

    if (Option.bSimPython) {
        // cout << "## Simulation_Python ##" << endl;
     
        WriteSimModule();

        cout << "Simulation_Python module has been generated: ";
        cout << Option.SimModuleFile << endl;
    }

//    if (Option.bSimCpp) {
//        // temporary C++ simulation code
//        cout << endl << "## Simulation_C++ ##" << endl;
//  
//        Simulation.Init(State, Dataset, DataManager, 100);
//        Simulation.Run(State, Context, Dataset);
//    }
            if (!Option.SimResultFile.empty()) {
		DataManager.SaveToFile(Option.SimResultFile.c_str());

                
 	    }       

    return 0;
}
