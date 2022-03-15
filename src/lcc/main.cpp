#include <iostream>
#include <queue>
#include <cassert>
#include <algorithm>
#include <new>
#include <string>
#include <vector>
#include <random>
#include <cctype>
#include <locale>
#include <stdio.h>
#include <string.h>

#ifdef _MSC_VER
#include <io.h>
//#include <sqlite3.h>

#else
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#endif

#include "node.h"
#include "option.h"
#include "context.h"
#include "simulation.h"
#include "datamanager.h"
#include "util.h"

#include "lpp.y.hpp"

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

ostream& os = std::cout;

// additional global variables
std::string in = "    ";
int N_MoleculesAllowed = 1; // Set number of molecules accepted for reactants and products
std::string Name_Pseudo = "Pseudo";
float Float_Init = Numbers::GetFloatDefault(); // random initialized float
int Int_Init = Numbers::GetIntDefault(); // random initialized int

// temporary simulation control parameters
int Sim_Steps = 1000;
int Sim_Resolution = 100;

const char *VersionString = "1.0.0";

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

float RandomNumber(float Min=0.0, float Max=1.0)
{
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<float> distr(Min, Max);

    return distr(eng);
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
       JoinedStr += Utils::SciFloat2Str(Float) + ", ";
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

void Update_N_MoleculesAllowed(int N_Molecules) {
    if (N_MoleculesAllowed < N_Molecules) {
        N_MoleculesAllowed = N_Molecules;
    }
}

// Node parsing routines
std::vector<std::pair<std::string, int>> GetStoichFromReaction(const NReaction* Reaction, bool bProductIsReaction=false) {

    map<string, int> Stoichiometry;
		string Location = Reaction->Location.Name;

    int N_Molecules = 0;
    std::string Name;
    int Coefficient;
    for (const auto& reactant : Reaction->Reactants) {
        Name = reactant->Id.Name;
        Coefficient = -reactant->Coeff; // update when coeff is fully implemented in parser
        std::cout << "    Reactant: " << "(" << Coefficient << ") " << Name << ", " << std::endl;
        Stoichiometry[Name] = Coefficient;
        N_Molecules++;
    }
    Update_N_MoleculesAllowed(N_Molecules);

    N_Molecules = 0;
    for (const auto& product : Reaction->Products) {
        Name = product->Id.Name;
        Coefficient = product->Coeff; // update when coeff is fully implemented in parser
        std::cout << "    Product: " << "(" << Coefficient << ") " << Name << ", " << std::endl;
        if (Stoichiometry.count(Name) > 0) {
            Coefficient += Stoichiometry[Name];
            std::cout << "    Updated Stoichiometry: " << "(" << Coefficient << ") " << Name << ", " << std::endl;
        }
        Stoichiometry[Name] = Coefficient;
        N_Molecules++;
    }
    Update_N_MoleculesAllowed(N_Molecules);

		if (!Location.empty()) {
        std::cout << "    Location: " << Location << endl;
		}

    // convert stoichiometry map to vector of pairs and add new molecules to the system
    std::vector<std::pair<std::string, int>> Stoichiometry_Ordered;

    for (auto& stoich : Stoichiometry) {
        std::pair<std::string, int> Stoich(stoich.first, stoich.second);
        Stoichiometry_Ordered.push_back(Stoich);

        // Do not add new molecule if the product is a reaction name
        if ((stoich.second > 0) & (bProductIsReaction)) {
            continue;
        }
        FMolecule * NewMolecule = new FMolecule(stoich.first);
        if (Option.bDebug) { NewMolecule->Print(os); }
        Context.AddToMoleculeList(NewMolecule);
    } 

    return Stoichiometry_Ordered;
}
            
void AddEnzReaction(std::string ReactionName, const NReaction* Reaction, std::string EnzymeName)
{
    float k = Float_Init;
    float krev = Float_Init;

    // properties
    const auto& propertylist = Reaction->Property;
    for (auto& property :propertylist) {
        auto& Key = property->Key;
        // auto Value = property->Value->Evaluate();
        const auto Value_Exp = dynamic_pointer_cast<const NConstantExpression>(property->Value);
        auto Value = Value_Exp->EvaluateValueAndPrefix();

        if (Key == "k") {
            k = std::stof(Value);
        } else if (Key == "krev") {
            krev = std::stof(Value);
        }
    }

    std::vector<std::pair<std::string, int>> Stoichiometry = GetStoichFromReaction(Reaction);
    string Location = Reaction->Location.Name;

    // for Enz_Standard Reaction
    // Fill in presumably irreversible reaction kinetic values 
    if      ((k != Float_Init) & (krev == Float_Init))   { krev = 0; }
    else if ((k == Float_Init) & (krev != Float_Init))   { k = 0; }

    // add new enzymatic reaction to the system
    if ((k >= 0) & (krev >= 0)) {
        FEnz_StandardReaction *NewReaction = new FEnz_StandardReaction(ReactionName, Stoichiometry, EnzymeName, k, krev);
        if (Option.bDebug) { NewReaction->Print(os); }
        Context.AddToReactionList(NewReaction);

    } else {
        FEnzymaticReaction *NewReaction = new FEnzymaticReaction(ReactionName, Stoichiometry, EnzymeName);
        if (Option.bDebug) { NewReaction->Print(os); }
        Context.AddToReactionList(NewReaction);
    }
}

std::pair<std::string, std::vector<float>> GetEnzKinetics(std::string EnzymeName, const NReaction* Reaction)
{
    std::pair<std::string, std::vector<float>> SubConstPair;
    std::string Substrate;
    float kcat = Float_Init;
    float KM = Float_Init;                

    // properties
    const auto& propertylist = Reaction->Property;
    for (auto& property :propertylist) {
        auto& Key = property->Key;

        if (Key == "Substrate") {
            auto Value = property->Value->Evaluate();
            Substrate = std::stof(Value);

        } else if ((Key == "kcat") || (Key == "kCat")) {
            const auto Value_Exp = dynamic_pointer_cast<const NConstantExpression>(property->Value);
            auto Value = Value_Exp->EvaluateValueAndPrefix();
            kcat = std::stof(Value);

        } else if ((Key == "KM") || (Key == "kM") || (Key == "km")) {
            const auto Value_Exp = dynamic_pointer_cast<const NConstantExpression>(property->Value);
            auto Value = Value_Exp->EvaluateValueAndPrefix();
            KM = std::stof(Value);
        }
    }
     

    // if Substrate not defined by user input, search the database            
    if (Substrate.empty()) {
        Substrate = Context.QueryTable(EnzymeName, "Substrate", Context.EnzymeTable);
        if (!Substrate.empty()) {
            os << "  Substrate imported from database: " << Substrate << endl;
        }
    }

    // import constants from database
    if (kcat == Float_Init) {
        string kcat_Database = Context.QueryTable(EnzymeName, "kcat", Context.EnzymeTable);
        if (!kcat_Database.empty()) {
            kcat = std::stof(kcat_Database);
            os << "  kcat imported from database: " << kcat_Database << endl;
        }
    }
    if (KM == Float_Init) {
        string KM_Database = Context.QueryTable(EnzymeName, "KM", Context.EnzymeTable);
        if (!KM_Database.empty()) {
            KM = std::stof(KM_Database);
            os << "  KM imported from database: " << KM_Database << endl;
        }
    }

    // Use first reactant as substrate for MichaelisMenten if still not assigned
    // This may not always work with Michaelis Menten without database. Excluding common molecules will improve a chance.
    if (Substrate.empty()) {
        for (const auto& reactant : Reaction->Reactants) {
            if (reactant->Id.Name == EnzymeName) {
                continue;
            }
            Substrate = reactant->Id.Name;
            break;
        }
    }

    Utils::Assertion(((kcat >= 0) & (KM >= 0)), "Kinetic constants are not properly extracted: " + EnzymeName);
 
    std::vector<float> KConstants{ kcat, KM };
    std::cout << EnzymeName << ", " << Substrate << ", " << kcat << ", " << KM << std::endl;
    SubConstPair = {Substrate, KConstants};
    
    return SubConstPair;
}

std::vector<float> ParseLocationFromVariable(const NVariableExpression * Variable)
{

}

// PseudoMolecule (placeholder to support current version of matrix operation)
void AddPseudoMolecule()
{
    // int Idx = 0
    FMolecule * NewMolecule = new FMolecule(Name_Pseudo);   // old count input system
    if (Option.bDebug) { NewMolecule->Print(os); }
    Context.AddToMoleculeList(NewMolecule);

    float Amount = 1;
    std::vector<float> Range = {0, -1, 0};
    bool bMolarity = false; // bMolarity gets reverted after traversal and data import if no molecule in the system has molarity unit

    FCount * NewCount = new FCount(Name_Pseudo, Amount, Range, bMolarity);
    if (Option.bDebug) { NewCount->Print(os); }
    Context.AddToCountList(NewCount);
}

void TraversalNode(NBlock* InProgramBlock)
{
    FTraversalContext tc(std::cerr);
    tc.Queue.push(InProgramBlock);

    // std::locale loc;

    os << endl << "## TraversalNode ##" << endl;
    AddPseudoMolecule();

    while(!tc.Queue.empty()) {
        const NNode* node = tc.Queue.front(); tc.Queue.pop();

        if (Utils::is_class_of<NReactionDeclaration, NNode>(node)) {
            auto N_Reaction = dynamic_cast<const NReactionDeclaration *>(node);
            os << "Reaction Id: " << N_Reaction->Id.Name << endl;
            // Reaction->Print(os);
        
            auto& Id = N_Reaction->Id;	    

            // Reaction Information to extract
            string Name = Id.Name;
        
            enum ReactionType {
        	Standard = 0,
        	Regulatory = 1,
            };
        
            ReactionType Type;
        
            std::vector<std::pair<std::string, int>> Stoichiometry;
        
            // for Standard reactions
            float k1 = Float_Init;
            float k2 = Float_Init;
        
            // for Regulatory reactions
            float K = Float_Init;
            float n = Float_Init; // TODO: Update how to indicate competitiveness. Temporarily, use n=-1 to indicate competitive mode for now
            string Effect;
            string Mode = "Allosteric"; // default setting

            // parse overall reaction
            auto& Reaction = N_Reaction->OverallReaction;
            // os << "  Reaction:" << endl;
        
            auto& bEffect = Reaction->bEffect;
            bool bProductIsReaction;
        
            const auto& propertylist = Reaction->Property;
            for (auto& property :propertylist) {
        	auto& Key = property->Key;
        	// auto Value = property->Value->Evaluate();
                const auto Value_Exp = dynamic_pointer_cast<const NConstantExpression>(property->Value);
                auto Value = Value_Exp->EvaluateValueAndPrefix();        

        	if (Key == "k") {
        	    k1 = std::stof(Value);
        	    Type = Standard;
        	} else if (Key == "krev") {
        	    k2 = std::stof(Value);
        	    Type = Standard;
        	} else if ((Key == "Ki") || (Key == "ki") || (Key == "Ka") || (Key == "ka") || (Key == "K")) {
        	    K = std::stof(Value);
        	    Type = Regulatory;
        	} else if (Key == "n") {
        	    n = std::stof(Value);
        					    // temporary code for competitive mode of inhibition regulation
        					    if (n == -1) {
        						Mode = "Competitive";
        	    }
        	} else {
        	    //                    os << "Unsupported reaction parameter: '" << property->Key << "' for the protein '" << Name << "'" << endl;
        	}
            }
        
            // Effect
            if ((!bEffect) & (K != Float_Init))       { Effect = "Inhibition"; } 
            else if ((bEffect) & (K != Float_Init))   { Effect = "Activation"; } 
        
            // Fill in presumably irreversible reaction kinetic values 
            if ((k1 != Float_Init) & (k2 == Float_Init)) { k2 = 0; }
            if ((k1 == Float_Init) & (k2 != Float_Init)) { k1 = 0; }
            if ((K != Float_Init) & (n == Float_Init))   { n = 1; }
        
            // Stoichiometry
            if      (Type == 0) { Stoichiometry = GetStoichFromReaction(Reaction, bProductIsReaction=false); } 
            else if (Type == 1) { Stoichiometry = GetStoichFromReaction(Reaction, bProductIsReaction=true); } // do not add the product to the molecule list
        
            // add new reaction to the system
            if (Type == 0) {
        	FStandardReaction *NewReaction = new FStandardReaction(Name, Stoichiometry, k1, k2);
        	if (Option.bDebug) { NewReaction->Print(os); }
        	Context.AddToReactionList(NewReaction);
        
            } else if (Type == 1) {
        	FRegulatoryReaction *NewReaction = new FRegulatoryReaction(Name, Stoichiometry, K, n, Effect, Mode);
            if (Option.bDebug) { NewReaction->Print(os); }
        	Context.AddToReactionList(NewReaction);
            }
 
        // This is intended for NEnzymeDeclaration, to be fixed later on.
        } else if (Utils::is_class_of<NProteinDeclaration, NNode>(node)) {
            auto N_Enzyme = dynamic_cast<const NProteinDeclaration *>(node);           
            os << "Enzyme Id: " << N_Enzyme->Id.Name << endl;
            // Enzyme->Print(os);
  
            auto& EnzymeName = N_Enzyme->Id.Name;	    
            auto& Reaction = N_Enzyme->OverallReaction;
            std::vector<std::pair<std::string, std::vector<float>>> Kinetics;

            if (Reaction) {   

                // parse overall reaction.
                // Note: Using enzyme name as reaction name for the overall reaction
                AddEnzReaction(EnzymeName, Reaction, EnzymeName);

                // Extract Substrate and MichaelisMenten Reaction Parameters
                bool bGetEnzKinetics = false;
                const auto& propertylist = Reaction->Property;
                for (auto& property :propertylist) {
                    auto& Key = property->Key;
                    if ((Key == "Substrate") || (Key == "kcat") || (Key == "kCat") || (Key == "KM") || (Key == "kM") || (Key == "km")) {
                        bGetEnzKinetics = true;
                        break;
                    }
                }
                if (bGetEnzKinetics) {
                    std::pair<std::string, std::vector<float>> SubConstPair = GetEnzKinetics(EnzymeName, Reaction);
                    Kinetics.push_back(SubConstPair);
                }
            }
  
            // TODO: if the block contains subreactions, priortize subreactions over main reaction for simulation?


            if (N_Enzyme->Block) {

                auto& Block = N_Enzyme->Block;
                int i_reaction = 0;
                for (auto& stmt: Block->Statements) {
                    os << "  "; stmt->Print(os);

                    // WITHOUT overall reaction: enzyme takes care of multiple reactions
                    if (!Reaction & (Utils::is_class_of<NReaction, NNode>(stmt.get()))) {

                        const auto Reaction = &*dynamic_pointer_cast<NReaction>(stmt);
                        os << "Reaction Id: " << Reaction->Id.Name << endl;
                        // Reaction->Print(os);

                        if (Reaction) {
                            // parse overall reaction.
                            std::string ReactionName = EnzymeName + "_" + Reaction->Id.Name;
                            AddEnzReaction(ReactionName, Reaction, EnzymeName);
            
                            // Extract Substrate and MichaelisMenten Reaction Parameters
                            bool bGetEnzKinetics = false;
                            const auto& propertylist = Reaction->Property;
                            for (auto& property :propertylist) {
                                auto& Key = property->Key;
                                if ((Key == "Substrate") || (Key == "kcat") || (Key == "kCat") || (Key == "KM") || (Key == "kM") || (Key == "km")) {
                                    bGetEnzKinetics = true;
                                    break;
                                }
                            }
                            if (bGetEnzKinetics) {
                                std::pair<std::string, std::vector<float>> SubConstPair = GetEnzKinetics(EnzymeName, Reaction);
                                Kinetics.push_back(SubConstPair);
                            }
                        }
                    i_reaction++;
                    }
                } // closing for stmt loop
                
            } // closing if block

            // add new enzyme to the system
            if (Kinetics.empty()) { 
                FEnzyme * NewEnzyme = new FEnzyme(EnzymeName);
                if (Option.bDebug) { NewEnzyme->Print(os); }
                Context.AddToMoleculeList(NewEnzyme);
            } 
            else { 
                FEnzyme * NewEnzyme = new FEnzyme(EnzymeName, Kinetics);
                if (Option.bDebug) { NewEnzyme->Print(os); }
                Context.AddToMoleculeList(NewEnzyme);
            }

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
            string Template = Context.QueryTable(Name, "Template", Context.PolymeraseTable);
            string Target = Context.QueryTable(Name, "Target", Context.PolymeraseTable);
            string Process = Context.QueryTable(Name, "Process", Context.PolymeraseTable);
            float Rate = std::stof(Context.QueryTable(Name, "Rate", Context.PolymeraseTable));

            FPolymerase * NewPolymerase = new FPolymerase(Name, Template, Target, Process, Rate);
            if (Option.bDebug) { NewPolymerase->Print(os); }
            Context.AddToMoleculeList(NewPolymerase);


#if 1
            for (const shared_ptr<NStatement>& stmt: N_Polymerase->Statements) {
                if (Utils::is_class_of<NElongationStatement, NStatement>(stmt.get())) {
                    const shared_ptr<NElongationStatement> elongstmt = dynamic_pointer_cast<NElongationStatement>(stmt);
                    // os << "---This is an elongation statement of the polymerase stmt---" << endl;
                    // elongstmt->Print(os);

                    NReaction ElongationReaction = elongstmt->Reaction;
                    os << "  Elongation:";
                    // ElongationReaction.Print(os);

                    os << "-----------------" << endl;
                } else if (Utils::is_class_of<NInitiationStatement, NStatement>(stmt.get())) {
                    const shared_ptr<NInitiationStatement> initstmt = dynamic_pointer_cast<NInitiationStatement>(stmt);
                    // os << "---This is an initiation statement of the polymerase stmt---" << endl;
                    // initstmt->Print(os);
                    // os << "-----------------" << endl;
                } else if (Utils::is_class_of<NTerminationStatement, NStatement>(stmt.get())) {
                    const shared_ptr<NTerminationStatement> termstmt = dynamic_pointer_cast<NTerminationStatement>(stmt);
                    // os << "---This is a termination statement of the polymerase stmt---" << endl;
                    // termstmt->Print(os);
                    // os << "-----------------" << endl;
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
            //                    NReaction EPickRandomNumberWithWeight>Reaction;
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
            //                            BuildingBlocks = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLT", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "SEL", "VAL"};
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




            // Temporary PolymeraseReactionCode
        } else if (Utils::is_class_of<NElongationStatement, NNode>(node)) {
            auto N_Elongation = dynamic_cast<const NElongationStatement *>(node);
            std::string Name;

            auto& ElongationReaction = N_Elongation->Reaction;

            os << "  Polymerase Reaction | Elongation:"; ElongationReaction.Print(os);
            std::vector<std::pair<std::string, int>> Stoichiometry;
            string Location = ElongationReaction.Location.Name;
            int Coefficient;
            std::vector<std::string> BuildingBlocks;

            for (const auto& reactant : ElongationReaction.Reactants) {
                const std::string& ReactantName = reactant->Id.Name;
                Coefficient = -reactant->Coeff;
                os << "    Reactant: " << "(" << Coefficient << ")" << ReactantName << ", " << endl;
                if ((ReactantName == "dna_{n}") | (ReactantName == "rna_{n}") | (ReactantName == "peptide_{n}")) {
                    continue;
                } else if (ReactantName == "dnt") {
                    Name = "pol1";
                    BuildingBlocks = {"dATP", "dCTP", "dGTP", "dUTP"};
                    continue;
                } else if (ReactantName == "nt") {
                    Name = "rnap";
                    BuildingBlocks = {"ATP", "CTP", "GTP", "UTP"};
                    continue;
                } else if (ReactantName == "aa") {
                    Name = "r1";
                    BuildingBlocks = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLT", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "SEL", "VAL"};
                    continue;
                }
                std::pair<std::string, int> Stoich(ReactantName, Coefficient);
                Stoichiometry.push_back(Stoich);

                FMolecule * NewMolecule = new FMolecule(ReactantName);
                if (Option.bDebug) { NewMolecule->Print(os); }
                Context.AddToMoleculeList(NewMolecule);
            }

            for (const auto& product : ElongationReaction.Products) {
                const std::string ProductName = product->Id.Name;
                Coefficient = product->Coeff;
                os << "    Product: " << "(" << Coefficient << ")" << ProductName << ", " << endl;
                if ((ProductName == "dna_{n+1}") | (ProductName == "rna_{n+1}") | (ProductName == "peptide_{n+1}")) {
                    continue;
                }
                std::pair<std::string, int> Stoich(ProductName, Coefficient);
                Stoichiometry.push_back(Stoich);

                FMolecule * NewMolecule = new FMolecule(ProductName);
                if (Option.bDebug) { NewMolecule->Print(os); }
                Context.AddToMoleculeList(NewMolecule);
            }

            if (!Location.empty()) {
                os << "    Location: " << Location << endl;
            }

            if (!BuildingBlocks.empty()){
                os << "    BuildingBlocks: [";
            }
            for (auto& BuildingBlock : BuildingBlocks) {
                FMolecule * NewMolecule = new FMolecule(BuildingBlock);
                os << NewMolecule->Name << ", ";
//                if (Option.bDebug) { NewMolecule->Print(os); }
                Context.AddToMoleculeList(NewMolecule);
            }
            os << "]" << endl;
            FPolymeraseReaction *NewReaction = new FPolymeraseReaction(Name, Stoichiometry, Name, BuildingBlocks);
            if (Option.bDebug) { NewReaction->Print(os); }
            Context.AddToReactionList(NewReaction);

        } else if (Utils::is_class_of<NOrganismDeclaration, NNode>(node)) {
            auto Organism = dynamic_cast<const NOrganismDeclaration *>(node);
            os << "Organism: " << Organism->Id.Name << endl;
            os << "  " << Organism->Description << endl;

            if (Utils::is_class_of<NEcoliStatement, NNode>(node)) {
                auto Ecoli = dynamic_cast<const NEcoliStatement *>(node);
                os << "Ecoli: " << Ecoli->Id.Name << endl;

                std::string OrganismSpace = "Ecoli";

                if (Ecoli->Block) {
                    for(const auto& stmt : Ecoli->Block->Statements) {
                        os << "  "; stmt->Print(os); os << endl;
                        
                    }
                }

                FOrganism * NewOrganism = new FOrganism(Ecoli->Id.Name, "Ecoli");
//                if (Option.bDebug) { NewOrganism->Print(os); }
                Context.AddToContainerList(NewOrganism);

            } else if (Organism->Id.Name == "ecoli") {

                if (Organism->Description == "E. coli K-12 MG1655") {
                    std::string Strain = "K-12 MG1655";
    
                    int ChromosomeSize = 4641652;
                    os << "Chromosome_I: " << std::to_string(ChromosomeSize) << "bp" << endl;
                    FChromosome * NewChromosome = new FChromosome("ChI", ChromosomeSize);
//                    if (Option.bDebug) { NewChromosome->Print(os); }
                    Context.AddToMoleculeList(NewChromosome);                
    
                    int i;
                    int i_cap = 5000;
    
                    i = 0;
                    os << ">> Genes being imported... : ";
                    for (auto& record : Context.GeneTable.Records) {
                        // os << record["symbol"] << ", ";
                        FGene * NewGene = new FGene(record["id"], record["symbol"]);
//                        if (Option.bDebug) { NewGene->Print(os); }
                        Context.AddToMoleculeList(NewGene);
    
                        i++;
                        // temporary capping
                        if (i == i_cap) {
                            os << "Gene importing is capped at " <<  std::to_string(i_cap) << endl;
                            break;
                        }
                    } os << "done" << endl;
                    // os << endl;
    
                    i = 0;
                    os << ">> RNAs being imported... : ";
                    for (auto& record : Context.RNATable.Records) {
                        // os << record["id"] << ", ";
                        FRNA * NewRNA = new FRNA(record["id"], record["type"]);
//                        if (Option.bDebug) { NewRNA->Print(os); }
                        Context.AddToMoleculeList(NewRNA);
    
                        i++;
                        // temporary capping
                        if (i == i_cap) {
                            os << "RNA importing is capped at " <<  std::to_string(i_cap) << endl;
                            break;
                        }
                    } os << "done" << endl;
                    // os << endl;
    
                    i = 0;
                    os << ">> Proteins being imported... : ";
                    for (auto& record : Context.ProteinTable.Records) {
                        // os << record["id"] << ", ";
                        FProtein * NewProtein = new FProtein(record["id"]);
//                        if (Option.bDebug) { NewProtein->Print(os); }
                        Context.AddToMoleculeList(NewProtein);
    
                        // temporary capping
                        i++;
                        if (i == i_cap) {
                            os << "Protein importing is capped at " <<  std::to_string(i_cap) << endl;
                            break;
                        }
                    } os << "done" << endl;
                    // os << endl;

                    FOrganism * NewOrganism = new FOrganism(Organism->Id.Name, Organism->Id.Name, Strain);
//                    if (Option.bDebug) { NewOrganism->Print(os); }
                    Context.AddToContainerList(NewOrganism);
    
                } else { 

                    FOrganism * NewOrganism = new FOrganism(Organism->Id.Name, Organism->Id.Name);
//                    if (Option.bDebug) { NewOrganism->Print(os); }
                    Context.AddToContainerList(NewOrganism);
                }
            }

        } else if (Utils::is_class_of<NAExpression, NNode>(node)) {
            auto AExpression = dynamic_cast<const NAExpression *>(node);
//            auto VarExp = dynamic_pointer_cast<const NVariableExpression *>(AExpression->OpA);
//            auto VarAssigned = dynamic_pointer_cast<const NExpression *>(AExpression->OpB); 

            if (Option.bDebug) {
                AExpression->Print(os);
                os << endl;
            }

            if (AExpression->Oper == T_ASSIGN) {

                // info to seek
                std::string Name;
                float Amount = Float_Init;
                std::vector<float> Range, Location;
                bool bMolarity = false;

                // get count from OpB
                if (Utils::is_class_of<const NConstantExpression, const NExpression>(AExpression->OpB.get())) {
                    const auto VarAssigned = dynamic_pointer_cast<const NConstantExpression>(AExpression->OpB); 

                    Amount = std::stof(VarAssigned->Evaluate());
                    bMolarity = VarAssigned->Molarity();

                } else if (Utils::is_class_of<const NVariableExpression, const NExpression>(AExpression->OpB.get())) {
                    const auto VarAssigned = dynamic_pointer_cast<const NVariableExpression>(AExpression->OpB); 
                    // TODO: Evaluate may not work yet
//                    Amount = std::stof(VarAssigned->Evaluate());
//                    bMolarity = VarAssigned->Molarity();
                }

                // get else from OpA
                if (Utils::is_class_of<const NVariableExpression, const NExpression>(AExpression->OpA.get())) {
                    const auto VarExp = dynamic_pointer_cast<const NVariableExpression>(AExpression->OpA);

                    // parsing VarExp (may be made into a function in the future)

                    Name = VarExp->GetName();

                    if (Utils::is_class_of<const NFunctionCallExpression, const NExpression>(VarExp->Variable.get())) {
                        const auto FCExp = dynamic_pointer_cast<const NFunctionCallExpression>(VarExp->Variable);

                        Location = FCExp->GetParameters("Location");
                    }

                    // VarExp->Index, for time info
                    if (VarExp->Index) {
                        if (Utils::is_class_of<const NRangeExpression, const NExpression>(VarExp->Index.get())) {
                            const auto RangeExp = dynamic_pointer_cast<const NRangeExpression>(VarExp->Index);

                            Range = RangeExp->GetBeginEndStep();
                        }
                    }

                } else if (Utils::is_class_of<const NFunctionCallExpression, const NExpression>(AExpression->OpA.get())) {
                    const auto FCExp = dynamic_pointer_cast<const NFunctionCallExpression>(AExpression->OpA);

                    Name = FCExp->GetName();
                    Location = FCExp->GetParameters("Location");
                }

                if (Option.bDebug) {
                    os << "@ AExpression parsing result" << endl;
                    os << "Name : " << Name << ", ";
                    os << "Amount: " << Utils::SciFloat2Str(Amount) << ", ";
                    os << "Range: ";
                    if (!Range.empty()) { os << "[" << JoinFloat2Str(Range) << "], "; }
                    else                { os << "None to [0], "; }
                    os << "bMolarity: "; if (bMolarity) { os << "true, ";  }
                                         else           { os << "false, "; }
                    os << "Location: ";
                    if (!Location.empty()) { os << "(" << JoinFloat2Str(Location) << ")"; }
                    else                   { os << "None"; }
                    os << endl;
                }

                // No range provided assumes [0] as input by default
                if (Range.empty()) {
                    Range = {0, 0, 0};
                }
                
					// temporary simulation control system
					if (Name == "SimSteps") {
					    Sim_Steps = static_cast<int>(Amount);
					    os << "# Temp Sim Control: SimSteps = " << Sim_Steps << endl;
					    continue;
					
					} else if (Name == "SimRes") {
					    Sim_Resolution = static_cast<int>(Amount);
					    os << "# Temp Sim Control: SimResolution = " << Sim_Resolution << endl;
					    continue;
					}

                // Add to Count and location
                FCount * NewCount = new FCount(Name, Amount, Range, bMolarity);
                if (Option.bDebug) { NewCount->Print(os); }
                Context.AddToCountList(NewCount);

                if (!Location.empty()) {
                    FLocation * NewLocation = new FLocation(Name, Location);
                    if (Option.bDebug) { NewLocation->Print(os); }
                    Context.AddToLocationList(NewLocation);
                }

            } // end of AExpression

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

void ImportCountFromDatabase()
{
    // user input
    std::vector<std::string> Names;
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;
        if (count->Begin == 0) {
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                Names.push_back(Name);
            }
        }
    }

    // TODO: Read in data
    // Find Molecule name in 'Names' vector to check if defined by user input
    // if not found in 'Names', make new FCount object and add to Context.CountList
    // TODO: Take out the var range parsing code to use it here again

// old database input for TCA
//        if (count.Initial < 0) {
//            for (auto& Pathway : Context.PathwayList) {
//                if (Pathway.Name == "TCA") {
//                    if (!Context.QueryTable(mol->Name, "Count", Context.InitialCountTable_TCA).empty()) {
//                        count.Initial = std::stof(Context.QueryTable(mol->Name, "Count", Context.InitialCountTable_TCA));
//                        std::cout << "InitialCount Imported | Molecule: " << mol->Name << "\t| Count: " << count.Initial << endl;
//                    }
//                }
//            }
//        }
    
}

std::string GetRegType(std::string Type)
{
    std::string RegType;

    if      (Type.find("Inhibition_Allosteric") != string::npos)   { RegType = "Regulatory_Inhibition_Allosteric"; }
    else if (Type.find("Inhibition_Competitive") != string::npos)  { RegType = "Regulatory_Inhibition_Competitive"; }
    else if (Type.find("Activation_Allosteric") != string::npos)   { RegType = "Regulatory_Activation_Allosteric"; }
    else if (Type.find("Unregulated") != string::npos)             { RegType = "Unregulated"; }

    return RegType;
}
void Print_InitializeStandardReaction(ofstream& ofs, std::string Type)
{
    ofs << in+ in+ "# " << Type << endl;

    // Standard vs. MichaelisMenten
    ofs << in+ in+ "self.Const_k_Reactant_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Const_k_Product_" << Type << " = None" << endl;
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = None" << endl;
    }
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Product_" << i << "_" << Type << " = None" << endl;
    }
    ofs << endl;

    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Const_n_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = None" << endl;
    }
    ofs << endl;
     
    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = None" << endl;
    ofs << endl;
}

void Print_SetUpStandardReaction(ofstream& ofs, std::string Type, std::vector<const FReaction *> ReactionSubList)
{
    // Types allowed: "Standard_Unregulated", "Standard_Inhibition", "Standard_Activation"
    // TODO: Multiplexing for regulation
    std::string RegType = GetRegType(Type);

    // for standard reactions 
    std::vector<float> k;
    std::vector<float> krev;

    // common variables
    int Idx_Pseudo = Context.GetIdxByName_MoleculeList(Name_Pseudo);
    std::vector<std::vector<int>> Idx_Reactants(N_MoleculesAllowed);
    std::vector<std::vector<int>> Idx_Products(N_MoleculesAllowed);

    // for regulatory mechanisms
    std::vector<float> K;
    std::vector<float> n;
    std::vector<int> Idx_Regulator; // reactant of the regulatory reaction

    // for stoichiometry matrix
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(Type); // TODO: Add new type
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(Type); // TODO: Add new type

    // relevant reaction lists
    std::vector<const FReaction *> RegulatoryReactionSubList = Context.GetSubList_ReactionList(RegType);

    for (auto Reaction : ReactionSubList) {
        const auto& reaction = dynamic_cast<const FStandardReaction *>(Reaction);
        
        k.push_back(reaction->k);
        krev.push_back(reaction->krev);

        std::vector<int> Idx_Reactant;
        std::vector<int> Idx_Product;

        for (auto& stoich : reaction->Stoichiometry) {
            int Idx = Context.GetIdxByName_MoleculeList(stoich.first);
            if (stoich.second <= 0) {
                Idx_Reactant.push_back(Idx);
            } else if (stoich.second >= 0) {
                Idx_Product.push_back(Idx);
            }
        }
        // Determines the number of substrates to handle (pseudomolecule is implemented to fill in blank spots for now)
        while (Idx_Reactant.size() != N_MoleculesAllowed) {
            Idx_Reactant.push_back(Idx_Pseudo);
        }
        while (Idx_Product.size() != N_MoleculesAllowed) {
            Idx_Product.push_back(Idx_Pseudo);
        }

        for (int i = 0; i < N_MoleculesAllowed; i++) {
            // std::cout << i << " : " << Idx_Reactant[i] << ", " << Idx_Product[i] << endl;
            Idx_Reactants[i].push_back(Idx_Reactant[i]);
            Idx_Products[i].push_back(Idx_Product[i]);
        }

        // Check Idx lengths
        Utils::Assertion((Idx_Reactants.size() == Idx_Products.size()), "# of parsed reactants and products do not match" );

        if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
            for (auto &Reg: RegulatoryReactionSubList) {
                const auto& reg = dynamic_cast<const FRegulatoryReaction *>(Reg);

                bool Import = false;
                std::string reactant_reg;
                std::string product_reg;
                for (auto &stoich: reg->Stoichiometry) {
                    if (stoich.second < 0) { reactant_reg = stoich.first; }
                    else if (stoich.second > 0) { product_reg = stoich.first; }
                    // import only if the product of the regulatory reaction targets the current standard reaction
                    if (stoich.first == reaction->Name) {
                        Import = true;
                    }
                }
                if (Import) {
                    Idx_Regulator.push_back(Context.GetIdxByName_MoleculeList(reactant_reg));

                    K.push_back(reg->K);
                    n.push_back(reg->n);
                    break;
                }
            }
        }
    }
 
    ofs << in+ in+ "# " << Type << endl;

    ofs << in+ in+ "self.Const_k_Reactant_" << Type << " = np.array([" << JoinFloat2Str(k) << "])" << endl;
    ofs << in+ in+ "self.Const_k_Product_" << Type << " = np.array([" << JoinFloat2Str(krev) << "])" << endl;
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Reactants[i]) << "])" << endl;
    }
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Product_" << i << "_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Products[i]) << "])" << endl;
    }
 
    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << Type << " = np.array([" << JoinFloat2Str(K) << "])" << endl;
        ofs << in+ in+ "self.Const_n_" << Type << " = np.array([" << JoinFloat2Str(n) << "])" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Regulator) << "])" << endl;
    }
    ofs << endl;
    
    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = np.asmatrix([" << Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;
}

void Print_InitializeEnzymeReaction(ofstream& ofs, std::string Type)
{
    ofs << in+ in+ "# " << Type << endl;

    ofs << in+ in+ "self.Idx_Enz_" << Type << " = None" << endl;

    // Standard vs. MichaelisMenten
    if (Type.find("Enz_Standard") != string::npos) {
        ofs << in+ in+ "self.Const_k_Reactant_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Const_k_Product_" << Type << " = None" << endl;
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = None" << endl;
        }
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Product_" << i << "_" << Type << " = None" << endl;
        }
 
    } else if (Type.find("Enz_MichaelisMenten") != string::npos) {
        ofs << in+ in+ "self.Const_kcat_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Const_KM_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Idx_EnzSub_" << Type << " = None" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Const_n_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = None" << endl;
    }
    ofs << endl;
    
    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = None" << endl;
    ofs << endl;
}

void Print_SetUpEnzymeReaction(ofstream& ofs, std::string Type, std::vector<const FReaction *> ReactionSubList) // to be changed with reaction list
{
    // TODO: Multiplexing for regulation
    std::string RegType = GetRegType(Type);

    // for standard reactions 
    std::vector<float> k;
    std::vector<float> krev;

    // for michaelis menten kinetics
    std::vector<float> kcats;
    std::vector<float> KMs;
    std::vector<int> Idx_EnzSub;

    // common variables
    int Idx_Pseudo = Context.GetIdxByName_MoleculeList(Name_Pseudo);
    std::vector<int> Idx_Enz; // for Enzyme where En is not included in the reaction
    std::vector<std::vector<int>> Idx_Reactants(N_MoleculesAllowed);
    std::vector<std::vector<int>> Idx_Products(N_MoleculesAllowed);

    // for regulatory mechanisms
    std::vector<float> K;
    std::vector<float> n;
    std::vector<int> Idx_Regulator;

    // for stoichiometry matrix
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(Type);
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(Type);

    // relevant reaction lists
    std::vector<const FReaction *> RegulatoryReactionSubList = Context.GetSubList_ReactionList(RegType);

    for (auto& Reaction : ReactionSubList) {
        const auto& reaction = dynamic_cast<const FEnzymaticReaction *>(Reaction);
        const auto& enzyme = Context.GetEnzyme_EnzymeList(reaction->Enzyme);

        Idx_Enz.push_back(Context.GetIdxByName_MoleculeList(enzyme->Name));

        // Enz_Standard
        if ((reaction->Type >= 10) & (reaction->Type < 20)) {
            const auto& reaction = dynamic_cast<const FEnz_StandardReaction *>(Reaction);
	    
            k.push_back(reaction->k);
            krev.push_back(reaction->krev);

            std::vector<int> Idx_Reactant;
            std::vector<int> Idx_Product;
    
            for (auto& stoich : reaction->Stoichiometry) {
                int Idx = Context.GetIdxByName_MoleculeList(stoich.first);
                if (stoich.second <= 0) {
                    Idx_Reactant.push_back(Idx);
                } else if (stoich.second >= 0) {
                    Idx_Product.push_back(Idx);
                }
            }
            // Determines the number of substrates to handle (pseudomolecule is implemented to fill in blank spots for now)
            while (Idx_Reactant.size() != N_MoleculesAllowed) {
                Idx_Reactant.push_back(Idx_Pseudo);
            }
            while (Idx_Product.size() != N_MoleculesAllowed) {
                Idx_Product.push_back(Idx_Pseudo);
            }
    
            for (int i = 0; i < N_MoleculesAllowed; i++) {
                // std::cout << i << " : " << Idx_Reactant[i] << ", " << Idx_Product[i] << endl;
                Idx_Reactants[i].push_back(Idx_Reactant[i]);
                Idx_Products[i].push_back(Idx_Product[i]);
            }

        // Enz_Michaelis Menten
        } else if ((reaction->Type >= 20) & (reaction->Type < 30)) {
            bool bSuccess = false;
            for (auto& kinetics : enzyme->Kinetics) {             
                // get substrate, kcat, KM when enzyme's substrate is found as a reactant of the reaction 
        
                if (Reaction->CheckIfReactant(kinetics.first)) {
                    Idx_EnzSub.push_back(Context.GetIdxByName_MoleculeList(kinetics.first));
                    auto& constants = kinetics.second;

                    kcats.push_back(constants[0]);
                    KMs.push_back(constants[1]);
                    bSuccess = true;
                    break;
                }
            }
            Utils::Assertion (bSuccess, "Unable to load kinetic constants for Michaelis Menten Reaction");
        }

        if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
            for (auto &Reg: RegulatoryReactionSubList) {
                const auto &reg = dynamic_cast<const FRegulatoryReaction *>(Reg);

                bool Import = false;
                std::string reactant_reg;
                std::string product_reg;
                for (auto &stoich: reg->Stoichiometry) {
                    if (stoich.second < 0) { reactant_reg = stoich.first; }
                    else if (stoich.second > 0) { product_reg = stoich.first; }
                    // import only if the product of the regulatory reaction targets the current standard reaction
                    if (stoich.first == reaction->Name) {
                        Import = true;
                    }
                }
                if (Import) {
                    Idx_Regulator.push_back(Context.GetIdxByName_MoleculeList(reactant_reg));
                    K.push_back(reg->K);
                    n.push_back(reg->n);
                    break;
                }
            }
        }
    }
 
    // Check Idx lengths
    Utils::Assertion((Idx_Reactants.size() == Idx_Products.size()), "# of parsed reactants and products do not match");

    ofs << in+ in+ "# " << Type << endl;
    ofs << in+ in+ "self.Idx_Enz_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Enz) << "])" << endl;

    // Standard vs. MichaelisMenten
    if (Type.find("Enz_Standard") != string::npos) {
        ofs << in+ in+ "self.Const_k_Reactant_" << Type << " = np.array([" << JoinFloat2Str(k) << "])" << endl;
        ofs << in+ in+ "self.Const_k_Product_" << Type << " = np.array([" << JoinFloat2Str(krev) << "])" << endl;
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Reactants[i]) << "])" << endl;
        }
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Product_" << i << "_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Products[i]) << "])" << endl;
        }
 
    } else if (Type.find("Enz_MichaelisMenten") != string::npos) {
        ofs << in+ in+ "self.Const_kcat_" << Type << " = np.array([" << JoinFloat2Str(kcats) << "])" << endl;
        ofs << in+ in+ "self.Const_KM_" << Type << " = np.array([" << JoinFloat2Str(KMs) << "])" << endl;
        ofs << in+ in+ "self.Idx_EnzSub_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_EnzSub) << "])" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << Type << " = np.array([" << JoinFloat2Str(K) << "])" << endl;
        ofs << in+ in+ "self.Const_n_" << Type << " = np.array([" << JoinFloat2Str(n) << "])" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Regulator) << "])" << endl;
    }
    ofs << endl;
    
    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = np.asmatrix([" << Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;
}

void Print_InitializePolymeraseReaction(ofstream& ofs, const FPolymerase* Polymerase)
{
    ofs << in+ in+ "self.Freq_BB_" << Polymerase->Target << "s = None" << endl;

    // Initiation
    ofs << in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.Idx_Template_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.Idx_TemplateSubset" << Polymerase->Process << " = None" << endl; // local indexing within the template population for mRNA in RNA for protein translation
    		ofs << in+ in+ "self.Weight_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.Pol_Threshold_" << Polymerase->Process << " = None" << endl;
    ofs << endl;

    // Check if template exists as a target of any polymerases in the system
    bool TemplateExists = false;
//    auto& PolymeraseList = Context.GetList_PolymeraseMoleculeList
    for (auto& PolymeraseInSystem : Context.GetList_Polymerase_MoleculeList()) {
        if (PolymeraseInSystem->Target == Polymerase->Template) {
        TemplateExists = true;
        break;
        }
    } 

    if (!TemplateExists) {
        ofs << in+ in+ "self.Len_Nascent" << Polymerase->Template << "s = None" << endl;
    }

    // Elongation
    ofs << in+ in+ "# " << Polymerase->Process << " Initialize ElongationReaction" << endl;
    ofs << in+ in+ "self.Rate_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.MaxLen_Nascent" << Polymerase->Target << "s = None" << endl;
    ofs << in+ in+ "self.Len_Nascent" << Polymerase->Target << "s = None" << endl;
    ofs << endl;

    ofs << in+ in+ "self.Idx_PolSub_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.Idx_PolBB_" << Polymerase->Process << " = None" << endl;
    ofs << endl;
}

void Print_SetUpPolymeraseReaction(ofstream& ofs, const FPolymerase* Polymerase, float Rate, std::string FreqBBFileName, std::string MaxLenFileName, int Idx_Pol, std::vector<int> Idx_Template, std::vector<int> Idx_TemplateSubset, std::vector<int> Idx_Target, std::vector<int> Idx_PolSub, std::vector<int> Idx_PolBB, int Threshold)
{
    ofs << in+ in+ "self.Freq_BB_" << Polymerase->Target << "s = np.asmatrix(np.load(r'" << FreqBBFileName << "'))" << endl;

    // Initiation
    ofs << in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = np.asmatrix([" << std::to_string(Idx_Pol) << "])" << endl;
    ofs << in+ in+ "self.Idx_Template_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Template) << "])" << endl;
    ofs << in+ in+ "self.Idx_TemplateSubset_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_TemplateSubset) << "])" << endl;
    ofs << in+ in+ "self.Idx_Target_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Target) << "])" << endl;
    		ofs << in+ in+ "self.Weight_" << Polymerase->Process << " = np.array([" << "1" << "])" << endl;
    ofs << endl;

    // Check if template exists as a target of any polymerases in the system
    bool TemplateExists = false;
    for (auto& PolymeraseInSystem : Context.GetList_Polymerase_MoleculeList()){
        if (PolymeraseInSystem->Target == Polymerase->Template){
        TemplateExists = true;
        break;
        }
    } 

    if (!TemplateExists) {
        ofs << in+ in+ "self.Len_Nascent" << Polymerase->Template << "s = np.asmatrix(np.full([10, self.Freq_BB_" << Polymerase->Target << "s.shape[0]], -1))" << endl;
        ofs << endl;
    }

    int Template_Min, Template_Max, Target_Min, Target_Max;

    if      (Polymerase->Process == "DNAReplication") 	        {Template_Min = 1; Template_Max = 2; Target_Min = 1; Target_Max = 2;}
    else if (Polymerase->Process == "RNATranscription") 	{Template_Min = 1; Template_Max = 2; Target_Min = 0; Target_Max = 1;}
    else if (Polymerase->Process == "ProteinTranslation") 	{Template_Min = 0; Template_Max = 1; Target_Min = 0; Target_Max = 1;}
    else    				                        {Template_Min = 0; Template_Max = 1; Target_Min = 0; Target_Max = 1;} // improve exception handling
        

    ofs << in+ in+ "self.Pol_Threshold_" << Polymerase->Process << " = " << std::to_string(Threshold) << endl;
    ofs << endl;
    ofs << in+ in+ "Count_Template_" << Polymerase->Process << " = np.random.randint(" << std::to_string(Template_Min) << ", high=" << std::to_string(Template_Max) << ", size=self.Idx_Template_" << Polymerase->Process << ".shape[1])" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Template_" << Polymerase->Process << ", Count_Template_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl; 

    ofs << in+ in+ "Count_Target_" << Polymerase->Process << " = np.random.randint(" << std::to_string(Target_Min) << ", high=" << std::to_string(Target_Max) << ", size=self.Idx_Target_" << Polymerase->Process << ".shape[1])" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Target_" << Polymerase->Process << ", Count_Target_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl; 


    // Elongation
    ofs << in+ in+ "# " << Polymerase->Process << " SetUp ElongationReaction" << endl;
    ofs << in+ in+ "self.Rate_" << Polymerase->Process << " = " << std::to_string(int(Rate)) << endl;
    ofs << in+ in+ "self.MaxLen_Nascent" << Polymerase->Target << "s = np.asmatrix(np.load(r'" << MaxLenFileName << "'))" << endl;
    ofs << in+ in+ "self.Len_Nascent" << Polymerase->Target << "s = np.asmatrix(np.full([10, self.Freq_BB_" << Polymerase->Target << "s.shape[0]], -1))" << endl;
    ofs << endl;

    ofs << in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = np.asmatrix([" << std::to_string(Context.GetIdxByName_MoleculeList(Polymerase->Name)) << "])" << endl;
    ofs << in+ in+ "self.Idx_PolSub_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_PolSub) << "])" << endl;
    ofs << in+ in+ "self.Idx_PolBB_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_PolBB) << "])" << endl;
    ofs << endl;

    ofs << in+ in+ "Count_Pol_" << Polymerase->Process << " = np.asmatrix(np.full(len(self.Idx_Pol_" << Polymerase->Process << "), 100))" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Pol_" << Polymerase->Process << ", Count_Pol_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl;

    ofs << in+ in+ "Count_PolSub_" << Polymerase->Process << " = np.asmatrix(np.full(len(self.Idx_PolSub_" << Polymerase->Process << "), 1000))" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_PolSub_" << Polymerase->Process << ", Count_PolSub_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl;

    ofs << in+ in+ "Count_PolBB_" << Polymerase->Process << " = np.random.randint(4000, high=5000, size=self.Freq_BB_" << Polymerase->Target << "s.shape[1])" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_PolBB_" << Polymerase->Process << ", Count_PolBB_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl; 

}

void Print_InitiationReaction(ofstream& ofs, const FPolymerase* Polymerase)
{
    ofs << in+ in+ "# " << Polymerase->Process << endl;
    ofs << in+ in+ "self.State.Len_Nascent" << Polymerase->Target << "s = self.Initiation(";
                       ofs << "self.State.Len_Nascent" << Polymerase->Template << "s, ";
                       ofs << "self.State.Len_Nascent" << Polymerase->Target << "s, ";
                       ofs << "self.State.Idx_Pol_" << Polymerase->Process << ", ";
                       ofs << "self.State.Idx_Template_" << Polymerase->Process << ", ";
                       ofs << "self.State.Idx_TemplateSubset_" << Polymerase->Process << ", ";
                       ofs << "self.State.Weight_" << Polymerase->Process << ", ";
                       ofs << "self.State.Pol_Threshold_" << Polymerase->Process << ") " << endl;
    ofs << endl;
}

void Print_ElongationReaction(ofstream& ofs, const FPolymerase* Polymerase)
{
    ofs << in+ in+ "# " << Polymerase->Process << endl;
    ofs << in+ in+ "self.State.Len_Nascent" << Polymerase->Target << "s = self.Elongation(";
                       ofs << "self.State.Len_Nascent" << Polymerase->Target << "s, ";
                       ofs << "self.State.MaxLen_Nascent" << Polymerase->Target << "s, ";
                       ofs << "self.State.Rate_" << Polymerase->Process << ", ";
                       ofs << "1, ";
                       ofs << "self.State.Freq_BB_" << Polymerase->Target << "s, ";
                       ofs << "self.State.Idx_PolSub_" << Polymerase->Process << ", ";
                       ofs << "self.State.Idx_PolBB_" << Polymerase->Process << ") " << endl;
    ofs << endl;
}

void Print_TerminationReaction(ofstream& ofs, const FPolymerase* Polymerase)
{
    ofs << in+ in+ "# " << Polymerase->Process << endl;
    ofs << in+ in+ "self.State.Len_Nascent" << Polymerase->Target << "s = self.Termination(";
                       ofs << "self.State.Len_Nascent" << Polymerase->Target << "s, ";
                       ofs << "self.State.MaxLen_Nascent" << Polymerase->Target << "s, ";
                       ofs << "self.State.Idx_Target_" << Polymerase->Process << ")" << endl;
    ofs << endl;
}

void Print_Initialize_SpatialSimulation(ofstream& ofs)
{
    ofs << in+ in+ "# Spatial Simulation" << endl;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto ObjLoc = Context.GetSubList_LocationList("Compartment");
    auto Organisms = Context.OrganismList;

    ofs << in+ in+ "self.Dist_Names = list()" << endl;
    // TODO: update to 3d array
    ofs << in+ in+ "self.Dist_All = list()" << endl;
    for (auto& location : MolLoc) {
        ofs << in+ in+ "self.Idx_Dist_" << location->Name << " = None" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Pos_Names = list()" << endl;
    ofs << in+ in+ "self.Pos_Name2Idx = dict()" << endl;

    std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
    for (auto& UniqueName : ObjUniqueNames) {
        ofs << in+ in+ "self.Idx_Pos_" << UniqueName << " = None" << endl;
    }
    ofs << in+ in+ "self.Pos_X = None" << endl;
    ofs << in+ in+ "self.Pos_Y = None" << endl;
    ofs << in+ in+ "self.Pos_Angle = None" << endl;
    ofs << in+ in+ "self.Pos_Threshold = None" << endl;
    ofs << endl;
}

void Print_SetUp_SpatialSimulation(ofstream& ofs)
{
    int Map_Width = 1200;
    int Map_Height = 800;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto ObjLoc = Context.GetSubList_LocationList("Compartment");

    for (auto& location : MolLoc) {
        ofs << in+ in+ "self.Idx_Dist_" << location->Name << " = None" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Dist_Names = [" << JoinStr2Str(Context.GetNames_LocationList("Molecule")) << "]" << endl;
    int i = 0;
    for (auto location : MolLoc) {
        auto Coord = location->Coord;
        auto Amount = Context.GetInitialCountByName_CountList(location->Name);
        ofs << in+ in+ "self.Idx_Dist_" << location->Name << " = np.asmatrix([" << i << "])" << endl;
        ofs << in+ in+ "Dist = sim.InitializeDistribution(" << Map_Width << ", " << Map_Height << ", " << Coord[0] << ", " << Coord[1] << ", " << Amount << ")" << endl;
        // TODO: update to 3d array
        ofs << in+ in+ "self.Dist_All.append(Dist)" << endl;
        i++;
    }
    ofs << endl;
    
    ofs << in+ in+ "self.Pos_Names = [" << JoinStr2Str(Context.GetNames_LocationList("Compartment")) << "]" << endl;

    i = 0;
    std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
    for (auto& UniqueName : ObjUniqueNames) {
        int Count = int(Context.GetInitialCountByName_CountList(UniqueName));
        ofs << in+ in+ "self.Idx_Pos_" << UniqueName << " = np.asmatrix([";
        for (int j = i; j < (Count); j++) {
            ofs << j << ", ";
        }
        ofs << "])" << endl;
        ofs << in+ in+ "self.Pos_Name2Idx['" << UniqueName << "'] = self.Idx_Pos_" << UniqueName << endl;
        i++;
    }

    ofs << in+ in+ "# Currently support X, Y, Angle, Threshold" << endl;
    ofs << in+ in+  "self.Pos_X = np.array([";
    for (auto location: ObjLoc) {
        auto Coord = location->Coord;
        ofs << Coord[0] << ", ";
    }
    ofs << "]) " << endl;

    ofs << in+ in+  "self.Pos_Y = np.array([";
    for (auto location: ObjLoc) {
        auto Coord = location->Coord;
        ofs << Coord[1] << ", ";
    }
    ofs << "]) " << endl;

    ofs << in+ in+  "self.Pos_Angle = np.array([";
    for (auto location: ObjLoc) { ofs << "0.0, ";
    }
    ofs << "]) " << endl;

    ofs << in+ in+  "self.Pos_Threshold = np.asmatrix([";
    for (auto location: ObjLoc) {
        ofs << "0.0, ";
//        ofs << Numbers::MultiplyByAvogadro(0.983405e-9) << ", ";
    }
    ofs << "]) " << endl;
    ofs << endl;
}

void WriteSimIdx() {}

void WriteSimModule()
{
    // write SimModule.py
    std::cout << std::endl << "Generating SimModule..." << std::endl;

    std::ofstream ofs(Option.SimModuleFile.c_str());
    std::string endl = "\n";

    // IMPORT
    ofs << "import os, sys" << endl;
    ofs << "import numpy as np" << endl;
    // ofs << "import tensorflow as tf" << endl;
    ofs << "from datetime import datetime" << endl;
    ofs << "import csv" << endl;
    ofs << "import SimFunctions as sim" << endl;
    // ofs << "import SimIdx as idx" << endl;
    ofs << "import plot" << endl;
    ofs << endl;

    //TODO: Take options from SimModule cmd line
    ofs << "N_SimSteps = " << Sim_Steps << endl;
    ofs << "SimStepTimeResolution = " << Sim_Resolution << endl;
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

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto ObjLoc = Context.GetSubList_LocationList("Compartment");

    int MatrixSize = 1;
    if (!ObjLoc.empty()) { 
        MatrixSize = ObjLoc.size();
    }
    ofs << in+ in+ "self.Count_All = np.zeros([" << MatrixSize << ", " << Context.MoleculeList.size() << "])" << endl;
    ofs << in+ in+ "self.dCount_All = np.zeros([" << MatrixSize << ", " << Context.MoleculeList.size() << "])" << endl;
    ofs << endl;

    if (!Context.LocationList.empty()) {
        Print_Initialize_SpatialSimulation(ofs);
    }

    // for standard reactions
    std::vector<std::string> ReactionTypes {"Standard_Unregulated", "Standard_Inhibition_Allosteric", "Standard_Activation_Allosteric"};

    std::vector<std::string> StandardReactionTypes = ReactionTypes; // to reuse later

    for (auto& Type : ReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Print_InitializeStandardReaction(ofs, Type);
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
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Print_InitializeEnzymeReaction(ofs, Type);
        }
    }

    // for polymerase reactions (Template-based)
    std::vector<const FPolymerase *> PolymeraseList = Context.GetList_Polymerase_MoleculeList();
//    std::vector<std::string> PolymeraseNames = Context.GetNames_PolymeraseList(PolymeraseList);
    std::vector<const FPolymeraseReaction *> PolymeraseReactionList = Context.GetList_Polymerase_ReactionList();

    if (!PolymeraseList.empty()) {
        for (auto& Polymerase : PolymeraseList) {
            Print_InitializePolymeraseReaction(ofs, Polymerase);
        }
    }

    ofs << in+ "def Initialize(self):" << endl;
    ofs << in+ in+ "self.Vol = 1" << endl;
    ofs << endl;

    // Print SetUp_SpatialReaction for all spatial simulation (to be updated)
    if (!Context.LocationList.empty()) {
        Print_SetUp_SpatialSimulation(ofs);
    }

    // Print SetUpStandardReaction for each Reaction Type 
    for (auto& Type : StandardReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Print_SetUpStandardReaction(ofs, Type, ReactionSubList);
        }
    }

    // Print SetUpEnzymeReaction for each Reaction Type 
    for (auto& Type : EnzReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Print_SetUpEnzymeReaction(ofs, Type, ReactionSubList);
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
 
            Print_SetUpPolymeraseReaction(ofs, Polymerase, Rate, FreqBBFileName, MaxLenFileName, Idx_Pol, Idx_Template, Idx_TemplateSubset, Idx_Target, Idx_PolSub, Idx_PolBB, Threshold);
        }
    }

    // for legends
    std::vector<std::string> MolNames = Context.GetNames_MoleculeList();
    ofs << in+ in+ "self.Mol_Names = [" << JoinStr2Str(MolNames) << "]" << endl;
    int i = 0;
    for (auto MolName : MolNames) {
        ofs << in+ in+ "self.Mol_Name2Idx['" << MolName << "'] = np.asmatrix([" << i << "])" <<endl;
        i++;
    }
    ofs << in+ in+ "self.Legends = ['SimStep', 'Vol', " << JoinStr2Str(MolNames) << "]" << endl;
    ofs << endl;

    // Initialize all molecule counts
    std::vector<int> Idx_Mol = Context.GetIdxListFromMoleculeList("Molecule");
    std::vector<float> InitialCount_Molecules;
    std::vector<float> MolarityFactor_Molecules;

    for (auto& molecule : Context.MoleculeList) {
        float Count = Context.GetInitialCountByName_CountList(molecule->Name);
        float MolarityFactor = Context.GetMolarityFactorByName_CountList(molecule->Name); 
        InitialCount_Molecules.push_back(Count);
        MolarityFactor_Molecules.push_back(MolarityFactor);
    }

    
    ofs << in+ in+ "Idx_Mol = np.asmatrix([";
    if (ObjLoc.empty()) {
        ofs << JoinInt2Str_Idx(Idx_Mol);
    } else {
        for (auto objLoc : ObjLoc) {
            ofs << "[" << JoinInt2Str_Idx(Idx_Mol) << "], ";
        }
    }
    ofs << "])" << endl;

    ofs << in+ in+ "Count_Mol = np.array([";
    if (ObjLoc.empty()) {
        ofs << JoinFloat2Str(InitialCount_Molecules);
    } else {
        for (auto objLoc : ObjLoc) {
            ofs << "[" << JoinFloat2Str(InitialCount_Molecules) << "], ";
        }
    }
    ofs << "])" << endl;

    ofs << in+ in+ "MolarityFactor_Mol = np.array([" << JoinFloat2Str(MolarityFactor_Molecules) << "])" << endl;
    ofs << in+ in+ "MolarityFactor_Mol = np.where(MolarityFactor_Mol == 1, self.Vol, 1)" << endl;
    ofs << in+ in+ "Count_Mol *= MolarityFactor_Mol" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, Idx_Mol, Count_Mol, axis=1)" << endl;
    ofs << endl;

    ofs << in+ "def GetMolNames(self):" << endl;
    ofs << in+ in+ "return self.Mol_Names" << endl;
    ofs << endl;

    ofs << in+ "def GetDistributionNames(self):" << endl;
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
        for (auto molLoc : MolLoc) {
            ofs << in+ in+ "self.Idx_DistToCoord_" << molLoc->Name << " = None" << endl;
        }
    }

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

    // Homeostasis --> to SimModule cmd option
    std::vector<std::string> HomeostasisList = {"Am"};// TODO: take user input from cmd line
    for (auto& name : HomeostasisList) {
        int Idx = Context.GetIdxByName_MoleculeList(name);
        ofs << in+ in+ "self.Idx_Count_Homeostasis_" << name << " = None" << endl;
        ofs << in+ in+ "self.Idx_Pos_Homeostasis_" << name << " = None" << endl;
        ofs << in+ in+ "self.Homeostasis_Prev_" << name << " = None" << endl;
    }
    ofs << endl;

//    old 
//    // Restore
//    for (auto& molecule : Context.MoleculeList) {
//        auto& count = molecule->Count;
//        if (count.Fixed) {
//            ofs << in+ in+ "self.Idx_Restore_" << molecule->Name << " = None" << endl;
//        }
//    }
//
//    // Event
//    for (auto& count : Context.CountList) {
//        if (Begin != 0) {
//            }
//        }
//    }

//    if (!Context.PathwayList.empty()){
//        for (auto& Pathway : Context.PathwayList) {
//            if (Pathway.Name == "TCA") {
//                ofs << in+ in+ "self.Idx_Restore_" << Pathway.Name << " = None" << endl;
//            }
//        }
//    }

    ofs << in+ in+ "# Debugging" << endl;
    ofs << in+ in+ "self.Debug_Idx_Molecules = list()" << endl;
    ofs << in+ in+ "self.UnitTxt = ''" << endl;
    ofs << in+ in+ "self.Unit = 0" << endl;
    ofs << endl;

    ofs << in+ "def Initialize(self, InN_SimSteps=1000, InTimeResolution=100):" << endl;
    ofs << in+ in+ "print('Simulation Initialized...')" << endl;
    ofs << in+ in+ "self.N_SimSteps = np.asmatrix([InN_SimSteps])" << endl;
    ofs << in+ in+ "self.SimTimeResolutionPerSecond = InTimeResolution" << endl;
    ofs << endl;

    // Distribution To Count
    i = 0;
    if (!MolLoc.empty()) {
        for (auto molLoc : MolLoc) {
            int Idx = Context.GetIdxByName_MoleculeList(molLoc->Name);
            ofs << in+ in+ "self.Idx_DistToCoord_" << molLoc->Name << " = np.array([" << std::to_string(Idx) << "])" << endl;
        i++;
        }
    }

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
                ofs << in+ in+ "self.Idx_Restore_" << Name << " = np.asmatrix([" << std::to_string(Idx) << "])" << endl;
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
                ofs << in+ in+ "self.Idx_Event_" << Name << " = np.asmatrix([" << std::to_string(Idx) << "])" << endl;
                Names.push_back(Name);
            }
        }
    }
    ofs << endl;

    // Homeostasis
    for (auto& name : HomeostasisList) {
        int Idx = Context.GetIdxByName_MoleculeList(name);
        ofs << in+ in+ "self.Idx_Count_Homeostasis_" << name << " = np.array([" << std::to_string(Idx) << "])" << endl;
        Idx = 0; // Get index of "E" ihe location list? Use MolName to connect to "E"
        ofs << in+ in+ "self.Idx_Pos_Homeostasis_" << name << " = np.array([" << std::to_string(Idx) << "])" << endl;
        ofs << in+ in+ "self.Homeostasis_Prev_" << name << " = 0" << endl;
    }
    ofs << endl;

//    if (!Context.PathwayList.empty()){
//        for (auto& Pathway : Context.PathwayList) {
//            if (Pathway.Name == "TCA") {
//                std::string MoleculeToRestore;
//                int Idx;
//
//                MoleculeToRestore = "acetyl-CoA";
//                Idx = Context.GetIdxByName_MoleculeList(MoleculeToRestore);
//                ofs << in+ in+ "self.Idx_Restore_" << Pathway.Name << " = np.asmatrix([" << Idx << "]) # " << MoleculeToRestore << endl;
//
//	                MoleculeToRestore.clear();
//                MoleculeToRestore = "malate";
//                Idx = Context.GetIdxByName_MoleculeList(MoleculeToRestore);
//                ofs << in+ in+ "self.Idx_Restore_" << Pathway.Name << " = np.asmatrix([" << Idx << "]) # " << MoleculeToRestore << endl;
//            }
//        }
//    }
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
    ofs << in+ in+ "self.Debug_SetUnit(Unit)" << endl;
    ofs << endl;

    // regular simloop
    ofs << in+ "def SimLoop_WithSpatialSimulation(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

    ofs << in+ in+ "# Run Spatial Simulation" << endl;
    ofs << in+ in+ "self.SpatialSimulation()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Run Reactions" << endl;
    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;

    if (Option.bDebug) {
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        ofs << endl;
    }

    ofs << in+ in+ "# Update Substrate Count" << endl;
    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Update Spatially Distributed Molecules On Count" << endl;
    ofs << in+ in+ "self.DistributionToCount()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Restore Substrate Count for Sustained Substrate Influx" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ "def SimLoop_WithoutSpatialSimulation(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;
    if (Option.bDebug) {
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        ofs << endl;
    }
    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << endl;

    // simloop for homeostasis
    ofs << in+ "def SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;
    if (Option.bDebug) {
        ofs << in+ in+ "self.Debug_PrintCounts(DisplayCount)" << endl;
        ofs << endl;
    }
    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << in+ in+ "self.DistributionToCount()" << endl;
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
    if (Context.LocationList.empty()) {
        ofs << in+ in+ "self.Dataset.Data = self.State.ExportData(self.SimStep/self.SimTimeResolutionPerSecond)" << endl;
        ofs << in+ in+ "self.DataManager.Add(self.Dataset.Data)" << endl;
    } else {
        ofs << in+ in+ "pass" << endl;
    }
    ofs << endl;

    ofs << in+ "def ApplySimTimeResolution(self, Rate):" << endl;
    ofs << in+ in+ "return Rate / self.SimTimeResolutionPerSecond" << endl;
    ofs << endl;

    // Distribution to Count

    // get coord for 'e'
    // get 'L' count from distribution by 'e' coord
    // get idx for 'L'
    // Update count at idx


    // here
    ofs << in+ "def DistributionToCount(self):" << endl;
    if (!MolLoc.empty() & !ObjLoc.empty()) {
        for (auto molLoc : MolLoc) {
            std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
            for (auto UniqueName : ObjUniqueNames) {
                ofs << in + in + "Count = self.GetCountFromDistributionByNameAndPos('" << molLoc->Name << "', " << "'" << UniqueName << "')" << endl;
                ofs << in + in + "self.State.Count_All[:, self.Idx_DistToCoord_" << molLoc->Name
                    << "] = Count.transpose()" << endl;
            }
        }
    } else {
        ofs << in+ in+ "pass" << endl;
    }
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

                ofs << in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Restore_" << Name << ", " << Amount << ", axis=1)" << endl;
                Names.push_back(Name);
            }
        }
    }
   
    if (Names.empty()) {
        ofs << in+ in+ "pass" << endl;
    }


//    if (!Context.PathwayList.empty()){
//        for (auto& Pathway : Context.PathwayList) {
//            if (Pathway.Name == "TCA") {
//                std::string MoleculeToRestore = "acetyl-CoA";
//                int Count = 446331; 
//                ofs << in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Restore_" << Pathway.Name << ", " << std::to_string(Count) << ", axis=1)   # " << MoleculeToRestore << endl;
//                Pass = false;
//            }
//        }
//    }
    ofs << endl;

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

    ofs << in + "def Homeostasis(self):" << endl;
    // TODO: Get homeostasis info from FPathway,
    float HomeostasisFactor = 1e-7;
    float ThresholdFactor = 0.99999;
    ofs << in+ in+ "print('Running simulation to achieve homeostasis for : ";
    for (auto& name : HomeostasisList) {
        ofs << name << ", ";
    }
    ofs << "')" << endl;

    for (auto& name : HomeostasisList) {
        ofs << in + in + "bNotHomeostasis_" << name << " = True" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "while (";
    for (auto& name : HomeostasisList) {
        ofs << "bNotHomeostasis_" << name;
        if (name != HomeostasisList.back()) {
            ofs << " and ";
        }
    }
    ofs << "):" << endl;

    ofs << in+ in+ in+ "self.SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution()" << endl;
    ofs << endl;

    for (auto& name : HomeostasisList) {
        int Idx = Context.GetIdxByName_MoleculeList(name);
        std::string Now = "Homeostasis_Now_" + name;
        std::string Prev = "self.Homeostasis_Prev_" + name;
        std::string Idx_Count = "self.Idx_Count_Homeostasis_" + name;
        std::string Idx_Pos = "self.Idx_Pos_Homeostasis_" + name;
        std::string Threshold = "self.State.Pos_Threshold[" + Idx_Pos + "]";

        ofs << in+ in+ in+ Now << " = self.GetCountByName('" << name << "')" << endl;
        ofs << in+ in+ in+ "if np.all(" << Now << " > 0) and np.all(abs(" << Now <<  " - " << Prev << ") / " << Now << " < " << Utils::SciFloat2Str(HomeostasisFactor) << "):" << endl;
        ofs << in+ in+ in+ in+ "bNotHomeostasis_" << name << " = False" << endl;
        ofs << in+ in+ in+ in+ "self.State.Pos_Threshold[" << Idx_Pos << ", :] = (" << Now << " * " << Utils::SciFloat2Str(ThresholdFactor) << ").transpose()"<< endl;
        ofs << in+ in+ in+ in+ "print('[Homeostasis] achieved for      : " << name << " @', self.Debug_ApplyUnit(" << Now << ").transpose(), self.UnitTxt)" << endl;
        ofs << in+ in+ in+ in+ "print('[Homeostasis] threshold set for : " << name << " @', self.Debug_ApplyUnit(" << Now << " * " << Utils::SciFloat2Str(ThresholdFactor) << ").transpose(), self.UnitTxt)" << endl;
        //ofs << in+ in+ in+ in+ "print('Homeostasis achieved for      : " << name << " @ {:.010f}'.format(self.Debug_ApplyUnit(" << Now << ")), self.UnitTxt)" << endl;
        //ofs << in+ in+ in+ in+ "print('Homeostasis threshold set for : " << name << " @ {:.010f}'.format(self.Debug_ApplyUnit(" << Now << " * " << Utils::SciFloat2Str(ThresholdFactor) << ")), self.UnitTxt)" << endl;
        ofs << in+ in+ in+ "self.Homeostasis_Prev_" << name << " = Homeostasis_Now_" << name << endl;
    }
    ofs << endl;
    //    print("[Homeostasis {:06d}] Glucose:{:.6f}{} Am:{:.6f}{}".format(self.SimCount, GlucoseLvl / Unit / NA, UnitTxt, self.Am / Unit / NA, UnitTxt))

    ofs << in+ in+ in+ "# Trigger Event on Substrate Count" << endl;
    ofs << in+ in+ in+ "self.TriggerEventMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ "# Spatial Simulation related routines" << endl;
    if (!Context.LocationList.empty()) {
        ofs << in + "def SpatialSimulation(self):" << endl;
        ofs << in + in + "self.SpatialDiffusion()" << endl;
        ofs << in + in + "self.SpatialLocation()" << endl;
        ofs << endl;
    }

    ofs << in+ "def SpatialDiffusion(self):" << endl;
    if (!MolLoc.empty()) {
        int N_Dist = Context.GetNames_LocationList("Molecule").size();
        // TODO: update to 3d array
        for (int i = 0; i < N_Dist; i++) {
            ofs << in+ in+ "self.State.Dist_All[" << i << "] = sim.DiffuseDistribution(self.State.Dist_All[" << i << "])" << endl;
        }
    } else {
        ofs << in+ in+ "pass" << endl;
    }
    ofs << endl;

    ofs << in+ "def SpatialLocation(self):" << endl;
    ofs << in+ in+ "HomeostasisMolecule = self.GetCount(self.Idx_Count_Homeostasis_Am)" << endl; // TODO:set up dynamic indexing
    ofs << in+ in+ "self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle = sim.BacterialChemotaxis(np.array(HomeostasisMolecule), self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle, self.State.Pos_Threshold)" << endl;
    ofs << endl;

    ofs << in+ "def NonSpatialSimulation(self):" << endl;
    if (!StandardReactionTypes.empty()) {
        for (auto &Type: StandardReactionTypes) {
            std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
            if (!ReactionSubList.empty()) {
                ofs << in + in + "self.StandardReactions()" << endl;
                ofs << endl;
                break;
            }
        }
    }

    if (!EnzReactionTypes.empty()) {
        for (auto &Type: EnzReactionTypes) {
            std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
            if (!ReactionSubList.empty()) {
                ofs << in + in + "self.EnzymaticReactions()" << endl;
                ofs << endl;
                break;
            }
        }
    }

    // TODO: encapsulate this part with each polymerase to allow more process-specific customization
    if (!PolymeraseList.empty()){
        ofs << in+ in+ "self.InitiationReactions()" << endl;
        ofs << in+ in+ "self.ElongationReactions()" << endl;
        ofs << in+ in+ "self.TerminationReactions()" << endl;
        ofs << endl;
    }

    ofs << in+ "# Biochemical Reaction related routines" << endl;
    // StandardReaction function
    for (auto& Type : StandardReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
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
            // Get Concentrations
            // Reactants and products
            for (auto& substrate : SubstrateTypes) {
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                        ofs << in+ in+ AmountTextStr << substrate << "_" << i << " = ";
                    if (bMolaritySys) {
                        ofs << "sim.CountToConc(self.State.Count_All[:, self.State.Idx_" << substrate << "_" << i << "_" << Type << "], self.State.Vol)" << endl;
                    } else {
                        ofs << "self.State.Count_All[:, self.State.Idx_" << substrate << "_" << i << "_" << Type << "]" << endl;
                    }
                } 
            }

            // Regulators 
            if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                ofs << in+ in+ AmountTextStr << "Regulator = ";
                if (bMolaritySys) {
                    ofs << "sim.CountToConc(self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "], self.State.Vol)" << endl;
                } else {
                    ofs << "self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "]" << endl;
                }
            } 

            // Calculate Rate
            for (auto& substrate : SubstrateTypes) {
                ofs << in+ in+ "Rate_" << substrate << " = sim.Eqn_" << Typing[substrate] << "_" << N_MoleculesAllowed << "(";

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
                ofs << in+ in+ "Rate_" << substrate << " = sim.CheckRateAndConc_" << N_MoleculesAllowed << "(Rate_" << substrate;
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << ", " << AmountTextStr << substrate << "_" << i;
                }
                ofs << ")" << endl;
            }

            // Tally Rates
            ofs << in+ in+ "Rate = Rate_Reactant - Rate_Product" << endl;

            if (bMolaritySys) {
                // Apply stoichiometry
                ofs << in+ in+ "dConc_Mol_InStoichMatrix = sim.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
                // Convert to counts
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = sim.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
            } else {
                // Apply stoichiometry
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = sim.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
            }
            // Apply delta counts for molecules in the stoichiometry matrix
            ofs << in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
            ofs << endl;
        }
    }

    // EnzReaction function // TODO: Add bMolaritySys option as above
    for (auto& Type : EnzReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {

            // Typing set up
            std::map<std::string, std::string> Typing;
            std::vector<std::string> SubstrateTypes { "Reactant", "Product" };
            Typing.emplace("Reactant", Type);
            Typing.emplace("Product", EnzReactionTypes[0]); // default type

            ofs << in+ "def EnzymaticReaction_" << Type << "(self):" << endl;
            // Get Concentrations
            // Enzyme
             ofs << in+ in+ "Conc_Enz = sim.CountToConc(self.State.Count_All[:, self.State.Idx_Enz_" << Type << "], self.State.Vol)" << endl;
            // Reactants and products or EnzSubstrate
            if (Type.find("Enz_Standard") != std::string::npos) {
                for (auto& substrate : SubstrateTypes) {
                    for (int i = 0; i < N_MoleculesAllowed; i++) {
                        ofs << in+ in+ "Conc_" << substrate << "_" << i << " = sim.CountToConc(self.State.Count_All[:, self.State.Idx_" << substrate << "_" << i << "_" << Type << "], self.State.Vol)" << endl;
                    } 
                }
            } else if (Type.find("Enz_MichaelisMenten") != std::string::npos) {
                ofs << in+ in+ "Conc_EnzSub = sim.CountToConc(self.State.Count_All[:, self.State.Idx_EnzSub_" << Type << "], self.State.Vol)" << endl;
            } else {
                Utils::Assertion(false, "Unsupported Enz Reaction Type: " + Type);
            }

            // Regulators 
            if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                ofs << in+ in+ "Conc_Regulator = sim.CountToConc(self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "], self.State.Vol)" << endl;
            }

            // Calculate Rate
            if (Type.find("Enz_Standard") != std::string::npos) {
                for (auto& substrate : SubstrateTypes) {
                    ofs << in+ in+ "Rate_" << substrate << " = sim.Eqn_" << Typing[substrate] << "_" << N_MoleculesAllowed << "(Conc_Enz, ";
    
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
                    ofs << in+ in+ "Rate_" << substrate << " = sim.CheckRateAndConc_" << N_MoleculesAllowed << "(Rate_" << substrate;
                    for (int i = 0; i < N_MoleculesAllowed; i++) {
                        ofs << ", Conc_" << substrate << "_" << i;
                    }
                    ofs << ")" << endl;
                }
    
                // Tally Rates
                ofs << in+ in+ "Rate = Rate_Reactant - Rate_Product" << endl;

            } else if (Type.find("Enz_MichaelisMenten") != std::string::npos) {
                ofs << in+ in+ "Rate = sim.Eqn_" << Type << "(Conc_Enz, Conc_EnzSub";
    
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
            ofs << in+ in+ "dConc_Mol_InStoichMatrix = sim.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
            // Convert to counts
            ofs << in+ in+ "dCount_Mol_InStoichMatrix = sim.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
            // Apply delta counts for molecules in the stoichiometry matrix
            ofs << in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
            ofs << endl;
        }
    }

    // Print StandardReaction for each Reaction Type 
    ofs << in+ "def StandardReactions(self):" << endl;
    bool PassSwitch = true;
    for (auto& Type : StandardReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            ofs << in+ in+ "self.StandardReaction_" << Type << "()" << endl;
            PassSwitch = false;
        }
    } 
    if (PassSwitch) {
        ofs << in+ in+ "pass" << endl;
    }
    ofs << endl;
    

    // Print EnzymeReaction for each Reaction Type 
    ofs << in+ "def EnzymaticReactions(self):" << endl;
    PassSwitch = true;
    for (auto& Type : EnzReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            ofs << in+ in+ "self.EnzymaticReaction_" << Type << "()" << endl;
            PassSwitch = false;
        }
    } 
    if (PassSwitch) {
        ofs << in+ in+ "pass" << endl;
    }
    ofs << endl;

    if (!PolymeraseList.empty()) {

        ofs << in+ "def InitiationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {
            
            Print_InitiationReaction(ofs, Polymerase);
        }


        ofs << in+ "def ElongationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {
            
            Print_ElongationReaction(ofs, Polymerase);
        }

        ofs << in+ "def TerminationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {
            
            Print_TerminationReaction(ofs, Polymerase);
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

    ofs << in+ "def GetConcentration(self, Idx):" << endl;
    ofs << in+ in+ "return sim.CountToConc(self.State.Count_All[:, Idx])" << endl;
    ofs << endl;

    ofs << in+ "def GetDistribution(self, Idx):" << endl;
    ofs << in+ in+ "return self.State.Dist_All[Idx]" << endl;
    ofs << endl;

    ofs << in+ "def AddTodCount(self, Idx, Values):" << endl;
    ofs << in+ in+ "dCountToAdd = np.zeros_like(self.State.dCount_All)" << endl;
    ofs << in+ in+ "np.put_along_axis(dCountToAdd, Idx, Values, axis=1)" << endl;
    ofs << in+ in+ "dCount_All_New = self.State.dCount_All + dCountToAdd" << endl;
    ofs << in+ in+ "ZeroTest = dCount_All_New + self.State.Count_All" << endl;
    ofs << in+ in+ "self.State.dCount_All =  np.where(ZeroTest < 0, dCount_All_New - ZeroTest, dCount_All_New)" << endl;
    ofs << endl;

    ofs << in+ "def UpdateCounts(self):" << endl;
    ofs << in+ in+ "self.State.Count_All += self.State.dCount_All" << endl;
    ofs << in+ in+ "self.CleardCounts()" << endl;

    ofs << in+ "def CleardCounts(self):" << endl;
    ofs << in+ in+ "self.State.dCount_All = np.zeros_like(self.State.dCount_All)" << endl;

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
    ofs << in+ in+ "return self.State.GetDistributionNames().index(Name)" << endl;
    ofs << endl;

    ofs << in+ "def GetDistributionByName(self, Name):" << endl;
    ofs << in+ in+ "Idx = self.GetDistIdx(Name)" << endl;
    ofs << in+ in+ "return self.GetDistribution(Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetCountFromDistribution(self, Dist, X, Y):" << endl;
    ofs << in+ in+ "return Dist[X.astype(int), Y.astype(int)]" << endl;
//    ofs << in+ in+ "return np.take(Dist, ([X,Y]))" << endl;
    ofs << endl;

    ofs << in+ "def GetCountFromDistributionByNameAndPos(self, NameOfDist, NameOfPos):" << endl;
    // temporary code
    ofs << in+ in+ "X, Y = self.GetPositionXYByName(NameOfPos)" << endl;
    ofs << in+ in+ "Dist = self.GetDistributionByName(NameOfDist)" << endl;
    ofs << in+ in+ "return self.GetCountFromDistribution(Dist, X, Y)" << endl;
//    ofs << in+ in+ "return np.take(Dist, ([X,Y]))" << endl;
    ofs << endl;

    ofs << in+ "# Temporary routines" << endl;

    ofs << in+ "def OverElongationCorrection(self, Len_Elongated, Max):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ "Len_Over = np.where(Len_Elongated > Max, Len_Elongated - Max, 0)" << endl;
    ofs << in+ in+ "return Len_Elongated - Len_Over" << endl;
    ofs << endl;

    ofs << in+ "def BuildingBlockConsumption(self, Freq, N_Elongated_PerSpecies):" << endl;
    ofs << in+ in+ "Raw = sim.DetermineAmountOfBuildingBlocks(Freq, N_Elongated_PerSpecies)" << endl;
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
    ofs << in+ in+ "Idx_Selected = sim.PickRandomIdx(Count_Pol_Avail, Idx_Template, Weight_Initiation)" << endl;
    ofs << in+ in+ "Len_Target_Initiated = sim.InsertZeroIntoNegOneElementInLenMatrix(Len_Target, Idx_Selected)" << endl;
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

    ofs << in+ "def Debug_SetIdxMoleculesToTrack(self):" << endl;
    ofs << in+ in+ "# Add a list of molecules to track for debugging every simulation step" << endl;
    ofs << in+ in+ "Debug_Names_Molecules = []" << endl; // TODO: take input from command line
    ofs << endl;

    ofs << in+ in+ "if Debug_Names_Molecules:" << endl;
    ofs << in+ in+ in+ "for Name in Debug_Names_Molecules:" << endl;
    ofs << in+ in+ in+ in+ "self.Debug_Idx_Molecules.append(self.State.GetMolNames().index(Name))" << endl;
    ofs << in+ in+ in+ "self.Debug_AllCountsSwitch = False" << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "self.Debug_Idx_Molecules = list(range(len(self.State.GetMolNames())))" << endl;
    ofs << in+ in+ in+ "self.Debug_AllCountsSwitch = True" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCounts(self, Switch):" << endl;
    ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
    ofs << in+ in+ "for Idx in self.Debug_Idx_Molecules:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintCount(Idx)" << endl;
    ofs << in+ in+ "print()" << endl;
    ofs << in+ in+ "if Switch:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintSimStepTime()" << endl;
    ofs << in+ in+ in+ "for Idx in self.Debug_Idx_Molecules:" << endl;
    ofs << in+ in+ in+ in+ "self.Debug_PrintdCount(Idx)" << endl;
    ofs << in+ in+ in+ "print()" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintSimStepTime(self):" << endl;
    ofs << in+ in+ "Time = self.GetSimTime()" << endl;
    ofs << in+ in+ "print(self.SimStep, '(', round(Time,3), 's)', end='\t| ')" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCount(self, Idx):" << endl;
    ofs << in+ in+ "print(' ' + self.State.GetMolNames()[Idx], end=': ')" << endl;
    ofs << in+ in+ "print('{:010e}'.format(self.Debug_ApplyUnit(self.State.Count_All[0][Idx])), self.UnitTxt, end=' | ')" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintdCount(self, Idx):" << endl;
    ofs << in+ in+ "print('d' + self.State.GetMolNames()[Idx], end=': ')" << endl;
    ofs << in+ in+ "print('{:010e}'.format(self.Debug_ApplyUnit(self.State.dCount_All[0][Idx])), self.UnitTxt, end=' | ')" << endl;
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


    // class FDataManager
    ofs << "class FDataManager:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "self.Legend = list()" << endl;
    ofs << in+ in+ "self.DataBuffer = list()" << endl;
    ofs << endl;

    ofs << in+ "def SetLegend(self, InLegend):" << endl;
    ofs << in+ in+ "self.Legend = InLegend" << endl;
    ofs << endl;

    ofs << in+ "def Add(self, InData):" << endl;
    ofs << in+ in+ "self.DataBuffer.append(InData)" << endl;
    ofs << endl;

    ofs << in+ "def SaveToFile(self, InFileName):" << endl;
    ofs << in+ in+ "with open(InFileName, 'w', newline='', encoding='utf-8') as OutFile:" << endl;
    ofs << in+ in+ in+ "TsvWriter = csv.writer(OutFile, delimiter='\\t')" << endl;
    ofs << in+ in+ in+ "if self.Legend:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter.writerow(self.Legend)" << endl;
    ofs << in+ in+ in+ "for Row in self.DataBuffer:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter.writerow(np.array(Row).flatten().tolist())" << endl;
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

    // Simulation Module
    ofs << in+ "Simulation.Initialize(N_SimSteps, SimStepTimeResolution)" << endl;
    ofs << in+ "Simulation.Run(Spatial=0) # 0: WithoutSpatialSimulation, 1: WithSpatialSimulation" << endl;
    ofs << endl;
    ofs << in+ "DataManager.SaveToFile('" << Option.SimResultFile.c_str() << "')" << endl;
    ofs << endl;

    // MAIN
    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "main()" << endl;
    if (!Option.bRunInOmVisim) {
        ofs << in + "plot.main()" << endl; // temporary for convenience
    }
    ofs << endl;

    cout << "Simulation_Python module has been generated: ";
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

    // Load genes.tsv
    if (!Option.bParseOnly)
    {
        std::cout<< endl << "## Loading Database ##" << std::endl;
        Context.Init(Option);
        if (Option.Verbose) {
            vector<string> Keys;
            Keys.emplace_back("symbol");
            Keys.emplace_back("name");
            Keys.emplace_back("id");
            Keys.emplace_back("rnaId");
            Context.GeneTable.Dump(Keys);

            Keys.clear();   
            Keys.emplace_back("id");
            Keys.emplace_back("name");
            Keys.emplace_back("type");
            Keys.emplace_back("location");
            Keys.emplace_back("geneId");
            Keys.emplace_back("monomerId");
            Context.RNATable.Dump(Keys);

            Keys.clear();
            Keys.emplace_back("id");
            Keys.emplace_back("name");
            Keys.emplace_back("location");
            Keys.emplace_back("geneId");
            Keys.emplace_back("rnaId");
            Context.ProteinTable.Dump(Keys);

            Keys.clear();
            Keys.emplace_back("reaction id");
            Keys.emplace_back("stoichiometry");
            Context.ReactionTable.Dump(Keys);

            Keys.clear();
            Keys.emplace_back("Name");
            Keys.emplace_back("Substrate");
            Keys.emplace_back("kcat");
            Keys.emplace_back("KM");
            Keys.emplace_back("Inhibitor");
            Keys.emplace_back("Ki");
            os << "# EnzymeTable #" << endl;
            Context.EnzymeTable.Dump(Keys);

            Keys.clear();
            Keys.emplace_back("Name");
            Keys.emplace_back("Template");
            Keys.emplace_back("Target");
            Keys.emplace_back("Process");
            Keys.emplace_back("Rate");
            os << "# PolymeraseTable #" << endl;
            Context.PolymeraseTable.Dump(Keys);

            Keys.clear();
            Keys.emplace_back("Name");
            Keys.emplace_back("Count");
            os << "# InitialCountTable_TCA #" << endl;
            Context.InitialCountTable_TCA.Dump(Keys);
        }
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

        if (!Option.bParseOnly) {
            TraversalNode(ProgramBlock);
            Context.Organize();

            if (Option.bDebug) {
                Context.PrintLists(os);

                // initial conditions
                Context.PrintInitialCounts(os);
                Context.PrintInitialLocations(os);
            }
        }

        delete ProgramBlock;
    }

    if (Option.bDebug) {
    
        for(const auto& name : Context.UsingModuleList) {
            cout << name << endl;
        }
    }

    if (Option.bParseOnly) {
        return 0;
    }

    if (Option.bSimPython) {

//        WriteSimIdx();
        WriteSimModule();

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
