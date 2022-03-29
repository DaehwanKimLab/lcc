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
#include "writer.h"
#include "util.h"

#include "lpp.y.hpp"

extern NBlock* ProgramBlock;
extern int yyparse();
extern int yylex_destroy();
extern FILE* yyin;

using namespace std;

FOption Option, *pOption;
FCompilerContext Context, *pContext;
FWriter Writer;

ostream& os = std::cout;

// additional global variables

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

void AddReaction(std::string ReactionName, const NReaction* Reaction)
{
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
        FStandardReaction *NewReaction = new FStandardReaction(ReactionName, Stoichiometry, k1, k2);
        if (Option.bDebug) { NewReaction->Print(os); }
        Context.AddToReactionList(NewReaction);

    } else if (Type == 1) {
        FRegulatoryReaction *NewReaction = new FRegulatoryReaction(ReactionName, Stoichiometry, K, n, Effect, Mode);
        if (Option.bDebug) { NewReaction->Print(os); }
        Context.AddToReactionList(NewReaction);
    }
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
            kcat = Value_Exp->EvaluateValueAndPrefixInFloat();

        } else if ((Key == "KM") || (Key == "kM") || (Key == "km")) {
            const auto Value_Exp = dynamic_pointer_cast<const NConstantExpression>(property->Value);
            KM = Value_Exp->EvaluateValueAndPrefixInFloat();
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

void ParseCountLocation_AExpression(const NAExpression *AExpression) {
    // info to seek
    std::string Name;
    float Amount = Float_Init;
    std::vector<float> Range, Location;
    bool bMolarity = false;

    // get count from OpB
    if (Utils::is_class_of<const NConstantExpression, const NExpression>(AExpression->OpB.get())) {
        const auto VarAssigned = dynamic_pointer_cast<const NConstantExpression>(AExpression->OpB);

        Amount = VarAssigned->EvaluateInFloat();
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

    } else if (Utils::is_class_of<const NFunctionCallExpression, const NExpression>(
            AExpression->OpA.get())) {
        const auto FCExp = dynamic_pointer_cast<const NFunctionCallExpression>(AExpression->OpA);

        Name = FCExp->GetName();
        Location = FCExp->GetParameters("Location");
    }

    if (Option.bDebug) {
        os << "@ AExpression parsing result" << endl;
        os << "Name : " << Name << ", ";
        os << "Amount: " << Utils::SciFloat2Str(Amount) << ", ";
        os << "Range: ";
        if (!Range.empty()) { os << "[" << Utils::JoinFloat2Str(Range) << "], "; }
        else { os << "None to [0], "; }
        os << "bMolarity: ";
        if (bMolarity) { os << "true, "; }
        else { os << "false, "; }
        os << "Location: ";
        if (!Location.empty()) { os << "(" << Utils::JoinFloat2Str(Location) << ")"; }
        else { os << "None"; }
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
    } else if (Name == "SimRes") {
        Sim_Resolution = static_cast<int>(Amount);
        os << "# Temp Sim Control: SimResolution = " << Sim_Resolution << endl;
    } else {
        // Add to Count and location
        FCount *NewCount = new FCount(Name, Amount, Range, bMolarity);
        if (Option.bDebug) { NewCount->Print(os); }
        Context.AddToCountList(NewCount);

        if (!Location.empty()) {
            FLocation *NewLocation = new FLocation(Name, Location);
            if (Option.bDebug) { NewLocation->Print(os); }
            Context.AddToLocationList(NewLocation);
        }
    }
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

void TraversalNode_Core(NNode * node)
{
    if (Utils::is_class_of<NReactionDeclaration, NNode>(node)) {
        auto N_Reaction = dynamic_cast<const NReactionDeclaration *>(node);
        os << "Reaction Id: " << N_Reaction->Id.Name << endl;
        // Reaction->Print(os);

        auto& Id = N_Reaction->Id;

        // Reaction Information to extract
        string Name = Id.Name;

        auto& Reaction = N_Reaction->OverallReaction;
        AddReaction(Name, Reaction);

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
        auto N_Pathway = dynamic_cast<const NPathwayDeclaration *>(node);
        // os << "Pathway: " << N_Pathway->Id.Name << endl;

        string Name = N_Pathway->Id.Name;
        vector<string> Sequence;

        if (N_Pathway->PathwayChainReaction) {
            auto& PathwayChainReaction = N_Pathway->PathwayChainReaction;
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

        if (N_Pathway->OverallReaction) {
            auto Reaction = N_Pathway->OverallReaction;
            // TODO: if the block contains subreactions, priortize subreactions over main reaction for simulation?

            os << "Reaction Id: " << Reaction->Id.Name << endl;
            // Reaction->Print(os);

            if (Reaction) {
                // parse overall reaction.
                std::string ReactionName = Name + "_" + Reaction->Id.Name;
                AddReaction(ReactionName, Reaction);
            }
        }

        if (N_Pathway->Block) {

            auto& Block = N_Pathway->Block;
            int i_reaction = 0;
            for (auto& stmt: Block->Statements) {
                os << "  "; stmt->Print(os);

                // WITHOUT overall reaction: enzyme takes care of multiple reactions
                if (!N_Pathway->OverallReaction & (Utils::is_class_of<NStatement, NNode>(stmt.get()))) {
                    TraversalNode_Core(stmt.get());
                }
            } // closing for stmt loop

        } // closing if block

        FPathway * NewPathway = new FPathway(Name); // Fixme
        if (Option.bDebug) { NewPathway->Print(os); }
        Context.AddToPathwayList(NewPathway);

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

    } else if (Utils::is_class_of<NAExpression, NNode>(node)) {
        auto AExpression = dynamic_cast<const NAExpression *>(node);
//            auto VarExp = dynamic_pointer_cast<const NVariableExpression *>(AExpression->OpA);
//            auto VarAssigned = dynamic_pointer_cast<const NExpression *>(AExpression->OpB); 

        if (Option.bDebug) {
            AExpression->Print(os);
            os << endl;
        }

        if (AExpression->Oper == T_ASSIGN) {
            ParseCountLocation_AExpression(AExpression);
        } // end of AExpression

    } else if (Utils::is_class_of<NLoopStatement, NNode>(node)) {
        auto LoopStatement = dynamic_cast<const NLoopStatement *>(node);

        if (Option.bDebug) {
            LoopStatement->Print(os);
            os << endl;
        }

        // info to extract
        std::string ControlVar_Name;
        int ControlVar_Start = Int_Init;
        int ControlVar_Final = Int_Init;
        int ControlVar_Increment = Int_Init;

        // Get ControlVar_Name from InitStatement (AExpression)
        if (!LoopStatement->InitStatements.empty()) {
            for (const auto& InitStatement : LoopStatement->InitStatements) {
                if (Utils::is_class_of<const NExpressionStatement, const NStatement>(InitStatement.get())) {
                    auto NExpStmt = dynamic_pointer_cast<const NExpressionStatement>(InitStatement);
                    if (Utils::is_class_of<const NAExpression, const NExpression>(NExpStmt->Expression.get())) {
                        auto AExpression = dynamic_pointer_cast<const NAExpression>(NExpStmt->Expression);

                        if (AExpression->Oper == T_ASSIGN) {
                            // get ControlVar_Name from OpA
                            if (Utils::is_class_of<const NVariableExpression, const NExpression>(
                                    AExpression->OpA.get())) {
                                const auto VarExp = dynamic_pointer_cast<const NVariableExpression>(
                                        AExpression->OpA);

                                ControlVar_Name = VarExp->GetName();
                            }

                            // get ControlVar_Start from OpB
                            if (Utils::is_class_of<const NConstantExpression, const NExpression>(
                                    AExpression->OpB.get())) {
                                const auto VarAssigned = dynamic_pointer_cast<const NConstantExpression>(
                                        AExpression->OpB);

                                ControlVar_Start = floor(VarAssigned->EvaluateInFloat());

                            } else if (Utils::is_class_of<const NVariableExpression, const NExpression>(
                                    AExpression->OpB.get())) {
                                const auto VarAssigned = dynamic_pointer_cast<const NVariableExpression>(
                                        AExpression->OpB);
                                // TODO: Evaluate may not work yet
                                //                    Amount = std::stoi(VarAssigned->Evaluate());
                            }
                        }
                    }
                }
            }
        } // end of InitStatement

        // Get ControlVar_Start from InitStatement (AExpression)
        if (LoopStatement->CondExpression) {
            if (Utils::is_class_of<const NAExpression, const NExpression>(LoopStatement->CondExpression.get())) {
                auto AExpression = dynamic_pointer_cast<const NAExpression>(LoopStatement->CondExpression);

                if ((AExpression->Oper == T_LT) || (AExpression->Oper == T_LE)) {
                    // get Var_Name from OpA to confirm the Control Var
                    if (Utils::is_class_of<const NVariableExpression, const NExpression>(AExpression->OpA.get())) {
                        const auto VarExp = dynamic_pointer_cast<const NVariableExpression>(AExpression->OpA);

                        Utils::Assertion(VarExp->GetName() == ControlVar_Name, "Inconsistent control variable in the loop expression");
                    }

                    // get ControlVar_Start from OpB
                    if (Utils::is_class_of<const NConstantExpression, const NExpression>(AExpression->OpB.get())) {
                        const auto VarAssigned = dynamic_pointer_cast<const NConstantExpression>(AExpression->OpB);

                        ControlVar_Final = floor(VarAssigned->EvaluateInFloat());

                    } else if (Utils::is_class_of<const NVariableExpression, const NExpression>(AExpression->OpB.get())) {
                        const auto VarAssigned = dynamic_pointer_cast<const NVariableExpression>(AExpression->OpB);
                        // TODO: Evaluate may not work yet
    //                    Amount = std::stoi(VarAssigned->Evaluate());
                    }
                }
                // consider inclusivity
                if (AExpression->Oper == T_LE) {
                    ControlVar_Final ++;
                }
            }
        } // end of CondExpression

        // Get ControlVar_Increment from LoopExpression (AExpression within AExpression)
        if (LoopStatement->LoopExpression) {
            if (Utils::is_class_of<const NAExpression, const NExpression>(LoopStatement->LoopExpression.get())) {
                auto AExpression = dynamic_pointer_cast<const NAExpression>(LoopStatement->LoopExpression);

                if (AExpression->Oper == T_ASSIGN) {
                    // get Var_Name from OpA to confirm the Control Var
                    if (Utils::is_class_of<const NVariableExpression, const NExpression>(AExpression->OpA.get())) {
                        const auto VarExp = dynamic_pointer_cast<const NVariableExpression>(AExpression->OpA);

                        Utils::Assertion(VarExp->GetName() == ControlVar_Name, "Inconsistent control variable in the loop expression");
                    }
                    // get ControlVar_Start from OpB (recursive)
                    if (Utils::is_class_of<const NAExpression, const NExpression>(AExpression->OpB.get())) {
                        auto AExp = dynamic_pointer_cast<const NAExpression>(AExpression->OpB);

                        if (AExp->Oper == T_PLUS) {
                            // get Var_Name from OpA to confirm the Control Var
                            if (Utils::is_class_of<const NVariableExpression, const NExpression>(AExp->OpA.get())) {
                                const auto VarExp = dynamic_pointer_cast<const NVariableExpression>(AExp->OpA);

                                Utils::Assertion(VarExp->GetName() == ControlVar_Name, "Inconsistent control variable in the loop expression");
                            }

                            // get ControlVar_Start from OpB
                            if (Utils::is_class_of<const NConstantExpression, const NExpression>(AExp->OpB.get())) {
                                const auto VarAssigned = dynamic_pointer_cast<const NConstantExpression>(AExp->OpB);

                                ControlVar_Increment = floor(VarAssigned->EvaluateInFloat());
                            }
                        }
                    }
                }
            }

        } // end of LoopExpression

        if (Option.bDebug) {
            os << "@ LoopExpression parsing result" << endl;
            os << "Control Variable Name : " << ControlVar_Name      << ", ";
            os << "Start: "                  << ControlVar_Start     << ", ";
            os << "Final: "                  << ControlVar_Final     << ", ";
            os << "Increment: "              << ControlVar_Increment << ", ";
        }
// Body
        if (LoopStatement->Body) {
// {Block: {Expression: {Arithmetic: {OpCode: 325, OpA: {FunctionCall: {Name: {Variable: {Name: Id: E}}, Args: [FunctionCall: {Name: {Variable: {Name: Id: rand}}, Args: [Constant: 400, Constant: 600, ]}, FunctionCall: {Name: {Variable: {Name: Id: rand}},
// Args: [Constant: 400, Constant: 600, ]}, FunctionCall: {Name: {Variable: {Name: Id: rand}},  Args: [Constant: 400, Constant: 600, ]}, ]}}, OpB: {Constant: 1}}}, }}

//            int N_Loop = floor((ControlVar_Final - ControlVar_Start) / ControlVar_Increment);
            for (int i = ControlVar_Start; i < ControlVar_Final - ControlVar_Start; i = i + ControlVar_Increment) {
                for (auto &statement: LoopStatement->Body->Statements) {
                    if (Utils::is_class_of<NExpressionStatement, NStatement>(statement.get())) {
                        auto ExpStmt = dynamic_cast<NExpressionStatement *>(statement.get());
                        auto AExpression = dynamic_cast<const NAExpression *>(ExpStmt->Expression.get());

                        if (AExpression->Oper == T_ASSIGN) {
                            ParseCountLocation_AExpression(AExpression);

                        }
                    }
                }

            }
        }

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

    } else if (Utils::is_class_of<NDeclarationStatement, NNode>(node)) {
        auto DeclStmt = dynamic_cast<const NDeclarationStatement *>(node);

        // temporary Transporter syntax parsing
        if (DeclStmt->Type == "transporter") {

            // extract the following info
            std::string Name_Reaction, Name_Transporter;
            std::vector<std::pair<std::string, int>> Stoichiometry;
//            float D = Float_Init;
            float ki = Float_Init;
            float ko = Float_Init;

            // parsing
            Name_Reaction = DeclStmt->Id.Name;
            Name_Transporter = DeclStmt->Id.Name;

            std::string MolName; // for Stoich

            if (DeclStmt->Initializer) {
                auto Inits = dynamic_cast<const NInitializerExpression *>(DeclStmt->Initializer.get());
                for (const auto Exp : Inits->ExpressionList) {
                    if (Utils::is_class_of<NVariableExpression, NExpression>(Exp.get())) {
                        auto VarExp = dynamic_pointer_cast<NVariableExpression>(Exp);
                        MolName = VarExp->GetName();
                    } else if (Utils::is_class_of<NAExpression, NExpression>(Exp.get())) {
                        auto AExp = dynamic_pointer_cast<NAExpression>(Exp);
                        if (AExp->Oper == T_ASSIGN) {
                            std::string Key;
                            float Value;
                            if (Utils::is_class_of<const NVariableExpression, const NExpression>(AExp->OpA.get())) {
                                const auto VarExp = dynamic_pointer_cast<const NVariableExpression>(AExp->OpA);
                                Key = VarExp->GetName();
                                Utils::Assertion((Key == "ki" || Key == "ko"), "Transporter syntax currently takes 'ki or ko' as a parameter key");
                            }
                            if (Utils::is_class_of<const NConstantExpression, const NExpression>(AExp->OpB.get())) {
                                const auto ConstExp = dynamic_pointer_cast<const NConstantExpression>(AExp->OpB);
                                Value = ConstExp->EvaluateInFloat();
                            }
                            // add info
                            if      (Key == "ki") { ki = Value; }
                            else if (Key == "ko") { ko = Value; }
                        }
                    }
                }
            }

            // stoichiometry set up
            std::pair<std::string, int> Stoich(MolName, -1);
            Stoichiometry.push_back(Stoich);

            // default 0 for unfilled rates
            if      (ki == Float_Init) { ki = 0; }
            else if (ko == Float_Init) { ko = 0; }

//            if (Option.bDebug) {
                os << "@ Transporter Expression parsing result" << endl;
                os << "Transporter Reaction Name : " << Name_Reaction    << ", ";
                os << "Target Molecule: "            << MolName          << ", ";
                os << "Transporter: "                << Name_Transporter << ", ";
                os << "ki: "                         << ki << ", ";
                os << "ko: "                         << ko << endl;
//            }

            // add Transporter object
            FTransporterReaction * NewReaction = new FTransporterReaction(Name_Reaction, Stoichiometry, Name_Transporter);
            if (Option.bDebug) { NewReaction->Print(os); }
            Context.AddToReactionList(NewReaction);

            FTransporter * NewTransporter = new FTransporter(Name_Transporter, ki, ko);
            if (Option.bDebug) { NewTransporter->Print(os); }
            Context.AddToMoleculeList(NewTransporter);
        }
    }

#if 0
    else if (Utils::is_class_of<NIdentifier, NNode>(node)) {
        auto Identifier = dynamic_cast<const NIdentifier *>(node);
        os << "Identifier: " << Identifier->Name  << endl;
        Context.IdentifierList.emplace_back(Identifier->Name);
    }
#endif
}

void TraversalNode(NBlock* InProgramBlock)
{
    FTraversalContext tc(std::cerr);
    tc.Queue.push(InProgramBlock);

    // std::locale loc;

    os << endl << "## TraversalNode_Code ##" << endl;
    AddPseudoMolecule();

    while(!tc.Queue.empty()) {
        const NNode* ConstNode = tc.Queue.front(); tc.Queue.pop();
        NNode* node = const_cast<NNode *>(ConstNode);

        if (Utils::is_class_of<NContainerStatement, NNode>(node)) {
            auto ContainerStmt = dynamic_cast<NContainerStatement *>(node);

            ContainerStmt->Print(os);

            FContainer *NewContainer = new FContainer(ContainerStmt->Id.Name);
            //                if (Option.bDebug) { NewOrganism->Print(os); }
            Context.AddToContainerList(NewContainer);

            if (ContainerStmt->Body) {
                for (auto &statement: ContainerStmt->Body->Statements) {
                    os << "  ";
                    statement->Print(os);

                    if (Utils::is_class_of<NOrganismDeclaration, NNode>(statement.get())) {
                        auto Organism = dynamic_cast<NOrganismDeclaration *>(statement.get());
                        os << "Organism: " << Organism->Id.Name << endl;
                        os << "  " << Organism->Description << endl;

                        if (Utils::is_class_of<NEcoliStatement, NNode>(statement.get())) {
                            auto Ecoli = dynamic_cast<NEcoliStatement *>(statement.get());
                            os << "Ecoli: " << Ecoli->Id.Name << endl;

                            std::string OrganismSpace = "Ecoli";

                            if (Ecoli->Block) {
                                for (auto &stmt: Ecoli->Block->Statements) {
                                    os << "  ";
                                    stmt->Print(os);
                                    os << endl;

                                    if (Utils::is_class_of<NExpressionStatement, NNode>(stmt.get())) {
                                        auto ExpStmt = dynamic_cast<NExpressionStatement *>(stmt.get());

                                        TraversalNode_Core(ExpStmt->Expression.get());
//
//                                        if (Utils::is_class_of<NAExpression, NNode>(ExpStmt->Expression.get())) {
//                                            auto AExpression = dynamic_cast<NAExpression *>(ExpStmt->Expression.get());
//
//                                            if (Option.bDebug) {
//                                                AExpression->Print(os);
//                                                os << endl;
//                                            }
//
//                                            if (AExpression->Oper == T_ASSIGN) {
//                                                ParseCountLocation_AExpression(AExpression);
//                                            } // end of AExpression
//                                        }
//                                    } else {
//                                        TraversalNode_Core(stmt.get());
                                    } else {
                                        TraversalNode_Core(stmt.get());
                                    }
                                }
                            }

                            FOrganism *NewOrganism = new FOrganism(Ecoli->Id.Name, "Ecoli");
                            //                if (Option.bDebug) { NewOrganism->Print(os); }
                            Context.AddToContainerList(NewOrganism);

                        } else if (Organism->Id.Name == "ecoli") {

                            if (Organism->Description == "E. coli K-12 MG1655") {
                                std::string Strain = "K-12 MG1655";

                                int ChromosomeSize = 4641652;
                                os << "Chromosome_I: " << std::to_string(ChromosomeSize) << "bp" << endl;
                                FChromosome *NewChromosome = new FChromosome("ChI", ChromosomeSize);
                                //                    if (Option.bDebug) { NewChromosome->Print(os); }
                                Context.AddToMoleculeList(NewChromosome);

                                int i;
                                int i_cap = 5000;

                                i = 0;
                                os << ">> Genes being imported... : ";
                                for (auto &record: Context.GeneTable.Records) {
                                    // os << record["symbol"] << ", ";
                                    FGene *NewGene = new FGene(record["id"], record["symbol"]);
                                    //                        if (Option.bDebug) { NewGene->Print(os); }
                                    Context.AddToMoleculeList(NewGene);

                                    i++;
                                    // temporary capping
                                    if (i == i_cap) {
                                        os << "Gene importing is capped at " << std::to_string(i_cap) << endl;
                                        break;
                                    }
                                }
                                os << "done" << endl;
                                // os << endl;

                                i = 0;
                                os << ">> RNAs being imported... : ";
                                for (auto &record: Context.RNATable.Records) {
                                    // os << record["id"] << ", ";
                                    FRNA *NewRNA = new FRNA(record["id"], record["type"]);
                                    //                        if (Option.bDebug) { NewRNA->Print(os); }
                                    Context.AddToMoleculeList(NewRNA);

                                    i++;
                                    // temporary capping
                                    if (i == i_cap) {
                                        os << "RNA importing is capped at " << std::to_string(i_cap) << endl;
                                        break;
                                    }
                                }
                                os << "done" << endl;
                                // os << endl;

                                i = 0;
                                os << ">> Proteins being imported... : ";
                                for (auto &record: Context.ProteinTable.Records) {
                                    // os << record["id"] << ", ";
                                    FProtein *NewProtein = new FProtein(record["id"]);
                                    //                        if (Option.bDebug) { NewProtein->Print(os); }
                                    Context.AddToMoleculeList(NewProtein);

                                    // temporary capping
                                    i++;
                                    if (i == i_cap) {
                                        os << "Protein importing is capped at " << std::to_string(i_cap) << endl;
                                        break;
                                    }
                                }
                                os << "done" << endl;
                                // os << endl;

                                FOrganism *NewOrganism = new FOrganism(Organism->Id.Name, Organism->Id.Name, Strain);
                                //                    if (Option.bDebug) { NewOrganism->Print(os); }
                                Context.AddToContainerList(NewOrganism);

                            } else {

                                FOrganism *NewOrganism = new FOrganism(Organism->Id.Name, Organism->Id.Name);
                                //                    if (Option.bDebug) { NewOrganism->Print(os); }
                                Context.AddToContainerList(NewOrganism);
                            }
                        }
                    } else if (Utils::is_class_of<NExpressionStatement, NNode>(statement.get())) {
                        auto ExpStmt = dynamic_cast<NExpressionStatement *>(statement.get());
                        if ((Utils::is_class_of<NAExpression, NNode>(ExpStmt->Expression.get())) ||
                            (Utils::is_class_of<NLoopStatement, NNode>(ExpStmt->Expression.get()))) {

                            TraversalNode_Core(ExpStmt->Expression.get());

                        }
                    } else {
                        TraversalNode_Core(statement.get());
                    }
                }
            }
        } else {
            TraversalNode_Core(node);
        }

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

//            if (Option.bDebug) {
                Context.PrintLists(os);

                // initial conditions
                Context.PrintInitialCounts(os);
                Context.PrintInitialLocations(os);
//            }
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

        Writer.LinkOptionContext(Option, Context);
        Writer.SetUpDefaultVariables(N_MoleculesAllowed, Name_Pseudo, Float_Init, Int_Init);

//        WriteSimIdx();
        Writer.SimModule(Sim_Steps, Sim_Resolution);
        cout << Option.SimModuleFile << std::endl;

        if (!Context.LocationList.empty()) {
            Writer.SimVis2D();
            cout << Option.SimVis2DFile << std::endl;
        }

    }

//    if (Option.bSimCpp) {
//        // temporary C++ simulation code
//        cout << endl << "## Simulation_C++ ##" << endl;
//  
//        Simulation.Init(State, Dataset, DataManager, 100);
//        Simulation.Run(State, Context, Dataset);
//    }

//    if (!Option.SimResultFile.empty()) {
//        DataManager.SaveToFile(Option.SimResultFile.c_str());
//    }

    return 0;
}
