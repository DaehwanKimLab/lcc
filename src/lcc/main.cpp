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

std::string in = "    ";
int N_MoleculesAllowed = 3; // Set number of molecules accepted for reactants and products

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

std::vector<std::pair<std::string, int>> GetStoichFromOverallReaction(const NReaction& OverallReaction) {

    map<string, int> Stoichiometry;
		string Location = OverallReaction.Location.Name;

    std::string Name;
    int Coefficient;
    for (const auto& reactant : OverallReaction.Reactants) {
        Name = reactant->Id.Name;
        Coefficient = -reactant->Coeff; // update when coeff is fully implemented in parser
        std::cout << "    Reactants: " << "(" << Coefficient << ") " << Name << ", " << std::endl;
        Stoichiometry[Name] = Coefficient;

    }

    for (const auto& product : OverallReaction.Products) {
        Name = product->Id.Name;
        Coefficient = product->Coeff; // update when coeff is fully implemented in parser
        std::cout << "    Products: " << "(" << Coefficient << ") " << Name << ", " << std::endl;
        if (Stoichiometry.count(Name) > 0) {
            Coefficient += Stoichiometry[Name];
            std::cout << "    Updated Stoichiometry: " << "(" << Coefficient << ") " << Name << ", " << std::endl;
        }
        Stoichiometry[Name] = Coefficient;
    }

		if (!Location.empty()) {
        std::cout << "    Location: " << Location << endl;
		}

    // convert stoichiometry map to vector of pairs and add new molecules to the system
    std::vector<std::pair<std::string, int>> Stoichiometry_Ordered;

    for (auto& stoich : Stoichiometry) {
        std::pair<std::string, int> Stoich(stoich.first, stoich.second);
        Stoichiometry_Ordered.push_back(Stoich);
 
        FMolecule * Molecule = new FMolecule(Name);
        // Molecule->Print(os);
        Context.AddToMoleculeList(Molecule);
    } 

    return Stoichiometry_Ordered;
}

void TraversalNode(NBlock* InProgramBlock)
{
    ostream& os = std::cout;
    FTraversalContext tc(std::cerr);
    tc.Queue.push(InProgramBlock);
    float Float_Init = -0.09876723; // random initialized float
    int Int_Init = -128; // random initialized int 
    // std::locale loc;

    os << endl << "## TraversalNode ##" << endl;

    // Pseudomolecule (placeholder to support current version of matrix operation)
    // Idx = 0
    FMolecule * Molecule = new FMolecule("Pseudo", 1, true);
    // Molecule->Print(os);
    Context.AddToMoleculeList(Molecule);

    while(!tc.Queue.empty()) {
        const NNode* node = tc.Queue.front(); tc.Queue.pop();

        if (Utils::is_class_of<NReactionDeclaration, NNode>(node)) {
            auto N_Reaction = dynamic_cast<const NReactionDeclaration *>(node);
            os << "Reaction Id: " << N_Reaction->Id.Name << endl;
            // Reaction->Print(os);

            auto& Id = N_Reaction->Id;	    
            auto& OverallReaction = N_Reaction->OverallReaction;
            // os << "  OverallReaction:" << endl;

            // Reaction Information
            string Name = Id.Name;
            std::vector<std::pair<std::string, int>> Stoichiometry = GetStoichFromOverallReaction(OverallReaction);

            float k1 = Float_Init;
            float k2 = Float_Init;
            float K = Float_Init;
            float n = Float_Init;
            string Effect;

            auto& bEffect = OverallReaction.bEffect;
            const auto& propertylist = OverallReaction.Property;
            for (auto& property :propertylist) {
                auto& Key = property->Key;
		auto& Value = property->Value;

                if (Key == "k") {
                    k1 = std::stof(Value);
                } else if (Key == "krev") {
                    k2 = std::stof(Value);
                } else if ((Key == "Ki") || (Key == "ki") || (Key == "Ka") || (Key == "ka") || (Key == "K")) {
                    K = std::stof(Value);
                } else if (Key == "n") {
                    n = std::stof(Value);
                } else {
//                    os << "Unsupported reaction parameter: '" << property->Key << "' for the protein '" << Name << "'" << endl;
                }
            }

            // Effect
            if ((!bEffect) & (K != Float_Init)) {
                Effect = "Inhibition";
            } else if ((bEffect) & (K != Float_Init)) {
                Effect = "Activation";
            } 

            // Fill in presumably irreversible reaction kinetic values 
            if ((k1 != Float_Init) & (k2 == Float_Init)) {
                k2 = 0;
            }

            if ((k1 == Float_Init) & (k2 != Float_Init)) {
                k1 = 0;
            }

            if ((K != Float_Init) & (n == Float_Init)) {
                n = 1;
            }

            // add new reaction to the system
            if ((k1 >= 0) & (k2 >= 0)) {
                FStandardReaction *Reaction = new FStandardReaction(Name, Stoichiometry, k1, k2);
                Reaction->Print(os);
                Context.AddToReactionList(Reaction);
            } else if (K >= 0) {
                FRegulatoryReaction *Reaction = new FRegulatoryReaction(Name, Stoichiometry, K, n, Effect);
                Reaction->Print(os);
                Context.AddToReactionList(Reaction);
            }

        }

        // This is inteded for NEnzymeDeclaration, to be fixed later on.
        enum EnzReactionType {
            Standard = 0,
            Standard_Inhibition_Allosteric,
            Standard_Inhibition_Competitive,
            Standard_Activation_Allosteric,
            Standard_Activation_Competitive,
            MichaelisMenten,
            MichaelisMenten_Inhibition_Allosteric,
            MichaelisMenten_Inhibition_Competitive,
            MichaelisMenten_Activation_Allosteric,
            MichaelisMenten_Activation_Competitive
        };

        if (Utils::is_class_of<NProteinDeclaration, NNode>(node)) {
            auto N_Enzyme = dynamic_cast<const NProteinDeclaration *>(node);
            os << "Enzyme Id: " << N_Enzyme->Id.Name << endl;
            // Enzyme->Print(os);

            auto& Id = N_Enzyme->Id;	    
            auto& OverallReaction = N_Enzyme->OverallReaction;
            // os << "  OverallReaction:" << endl;

            // Enzyme Information
            string Name = Id.Name;

            EnzReactionType Type = Standard;

            string Substrate;
            float k1 = Float_Init;
            float k2 = Float_Init;

            string Regulator; // Inhibitor or Activator will be assigned as a Regulator
            string Inhibitor;
            string Activator;
            string Mode; // "Allosteric" or "Competitive"
            float K = Float_Init; // Reduced regulatory kinetic constant
            float n = Float_Init; // Hill's coefficient

            float InitialCount = Float_Init;
            bool Fixed = false;
            vector<pair<pair<float, float>, float>> Ranges;

            //parse properties
            const auto& propertylist = OverallReaction.Property;
            for (auto& property :propertylist) {
                auto& Key = property->Key;
		auto& Value = property->Value;

                if (Key == "k") {
                    k1 = std::stof(Value);
                    Type = Standard;
                } else if (Key == "krev") {
                    k2 = std::stof(Value);
                    Type = Standard; 
                } else if ((Key == "kcat") || (Key == "kCat")) {
                    k1 = std::stof(Value);
                    Type = MichaelisMenten;
                } else if ((Key == "KM") || (Key == "kM") || (Key == "km")) {
                    k2 = std::stof(Value);
                    Type = MichaelisMenten;

                } else if ((Key == "inhibitor") || (Key == "Inhibitor")) {
                    Inhibitor = Value; // TODO: improve coding style
                    Regulator = Value;
                } else if ((Key == "activator") || (Key == "Activator")) {
                    Activator = Value; // TODO: improve coding style
                    Regulator = Value;
                } else if ((Key == "mode") || (Key == "Mode")) {
                    Mode = Value;
                    if (Type == Standard) {
                        if (!Inhibitor.empty() & Activator.empty()){
                            if ((Mode == "allosteric") || ("Allosteric")) 	{ Type = Standard_Inhibition_Allosteric; }
                       else if ((Mode == "competitive") || ("Competitive"))	{ Type = Standard_Inhibition_Competitive; }
                        } else if (Inhibitor.empty() & !Activator.empty()){
                            if ((Mode == "allosteric") || ("Allosteric")) 	{ Type = Standard_Activation_Allosteric; }
                       else if ((Mode == "competitive") || ("Competitive")) 	{ Type = Standard_Activation_Allosteric; }
                        }
                    } else if (Type == MichaelisMenten) {
                        if (!Inhibitor.empty() & Activator.empty()){
                            if ((Mode == "allosteric") || ("Allosteric")) 	{ Type = MichaelisMenten_Inhibition_Allosteric; }
                       else if ((Mode == "competitive") || ("Competitive")) 	{ Type = MichaelisMenten_Inhibition_Competitive; }
                        } else if (Inhibitor.empty() & !Activator.empty()){
                            if ((Mode == "allosteric") || ("Allosteric")) 	{ Type = MichaelisMenten_Activation_Allosteric; }
                       else if ((Mode == "competitive") || ("Competitive")) 	{ Type = MichaelisMenten_Activation_Competitive; }
                        }
                    }
                } else if ((Key == "Ki") || (Key == "ki") || (Key == "Ka") || (Key == "ka")) {
                    K = std::stof(Value);
                } else if (Key == "n") {
                    n = std::stof(Value);

                } else if (Key == Name) {
                    // temporary parsing code

                    // Value: [0]=0, currently replaced with _0__0
                    // Value: [1:3]=2, _1to3__2

                    
                    std::vector<std::string> Value_Parsed;
                    Utils::tokenize(Value, "_", Value_Parsed);
                    string MolName = Value_Parsed[0];
                    string Range = Value_Parsed[1];
                    float Count = std::stof(Value_Parsed[2]);
                    os << "\tValue_Parsed | MolName: " << MolName << "Range: " << Range << ", Count: " << Count << endl;

                    string Delimiter = "to";
                    if (Range.find(Delimiter) != string::npos) {
                        std::vector<std::string> Range_Parsed;
                        Utils::tokenize(Range, Delimiter, Range_Parsed);
                        float Range_Begin = std::stof(Range_Parsed[0]);
                        float Range_End = std::stof(Range_Parsed[1]);
                        os << "\tRange_Parsed | Range_Begin: " << Range_Parsed[0] << ", Range_End: " << Range_Parsed[1] << endl;
                        std::pair<float, float> Range(Range_Begin, Range_End);
                        std::pair<std::pair<float, float>, float> Range_Count(Range, Count);
                        Ranges.push_back(Range_Count);

                    } else {
                        InitialCount = Count;
                        os << "\tInitial_Count: " << InitialCount << endl;
                    }

                } else if (Key == "Fixed") {
                    if (Value == Name) {
                        Fixed = true;
                    }
                
                } else {
//                    os << "Unsupported reaction parameter: '" << property->Key << "' for the protein '" << Name << "'" << endl;
                }
            }

            // if not defined by user input, search the database            
            if (Substrate.empty()) {
                Substrate = Context.QueryTable(Name, "Substrate", Context.EnzymeTable);
                if (!Substrate.empty()) {
                    os << "  Substrate imported from database: " << Substrate << endl;
                }
            }
 
            if (k1 == Float_Init) {
                string kcat_Database = Context.QueryTable(Name, "kcat", Context.EnzymeTable);
                if (!kcat_Database.empty()) {
                    k1 = std::stof(kcat_Database);
                    os << "  kcat imported from database: " << kcat_Database << endl;
                }
            }
            if (k2 == Float_Init) {
                string KM_Database = Context.QueryTable(Name, "KM", Context.EnzymeTable);
                if (!KM_Database.empty()) {
                    k2 = std::stof(KM_Database);
                    os << "  KM imported from database: " << KM_Database << endl;
                }
            }

            // Fill in presumably irreversible reaction kinetic values 
            if ((k1 != Float_Init) & (k2 == Float_Init)) {
                k2 = 0;
            }

            if ((k1 == Float_Init) & (k2 != Float_Init)) {
                k1 = 0;
            }
          
//            if (Inhibitor.empty()) {
//                Inhibitor = Context.QueryTable(Name, "Inhibitor", Context.EnzymeTable);
//                if (!Inhibitor.empty()) {
//                    os << "  Inhibitor imported from database: " << Substrate << endl;
//                }
//            }
//
//            if (Activator.empty()) {
//                Activator = Context.QueryTable(Name, "Activator", Context.EnzymeTable);
//                if (!Inhibitor.empty()) {
//                    os << "  Activator imported from database: " << Substrate << endl;
//                }
//            }
//
//            if (Mode.empty()) {
//                Mode = Context.QueryTable(Name, "Mode", Context.EnzymeTable);
//                if (!Mode.empty()) {
//                    os << "  Mode imported from database: " << Substrate << endl;
//                }
//            }
//
//            if (!Inhibitor.empty() and (K == Float_Init)) {
//                string Ki_Database = Context.QueryTable(Name, "Ki", Context.EnzymeTable); // * RandomNumber();
//                if (!Ki_Database.empty()) {
//                    K = std::stof(Ki_Database);
//                    os << "  Ki imported from database: " << Ki_Database << endl;
//                }
//            } else if (!Activator.empty() and (K == Float_Init)) {
//                string Ka_Database = Context.QueryTable(Name, "Ka", Context.EnzymeTable); // * RandomNumber();
//                if (!Ka_Database.empty()) {
//                    K = std::stof(Ka_Database);
//                    os << "  Ka imported from database: " << Ka_Database << endl;
//                }
//
//            }
//
//            if (n == Float_Init) {
//                string n_Database = Context.QueryTable(Name, "n", Context.EnzymeTable); // * RandomNumber();
//                if (!n_Database.empty()) {
//                    n = std::stof(n_Database);
//                    os << "  n imported from database: " << n_Database << endl;
//                }
//            }


            std::vector<std::pair<std::string, int>> Stoichiometry;
			string Location = OverallReaction.Location.Name;

            int Coefficient;
            for (const auto& reactant : OverallReaction.Reactants) {
                if (reactant->Id.Name == Name) {
                    break;
                }

                // This may not always work with Michaelis Menten without database. Excluding common molecules will improve a chance.
                if (Substrate.empty()) {
                    Substrate = reactant->Id.Name;
                }

                Coefficient = -reactant->Coeff; // update when coeff is fully implemented in parser
                os << "    Reactants: " << "(" << Coefficient << ") " << reactant->Id.Name << ", " << endl;

                // if also found in products, set the first stoichiometry coefficient as the sum of their coefficients.
                for (const auto& product : OverallReaction.Products) {
                    if (product->Id.Name == reactant->Id.Name) {
                        Coefficient = 0; // TODO: use sum of the coefficients when coefficient node is implemented.
                    }
                } 
                std::pair<std::string, int> Stoich(reactant->Id.Name, Coefficient);
                Stoichiometry.push_back(Stoich);

                float InitialCount = Float_Init;
                bool Fixed = false;
                vector<pair<pair<float, float>, float>> Ranges;

                const auto& propertylist = OverallReaction.Property;
                for (auto& property :propertylist) {
                    if (property->Key == reactant->Id.Name) {
                        // temporary parsing code
    
                        // Value: [0]=0, currently replaced with _0__0
                        // Value: [1:3]=2, _1to3__2
    
                        string Value = property->Value;
                        std::vector<std::string> Value_Parsed;
                        Utils::tokenize(Value, "_", Value_Parsed);
                        string MolName = Value_Parsed[0];
                        string Range = Value_Parsed[1];
                        float Count = std::stof(Value_Parsed[2]);
                        os << "\tValue_Parsed | MolName: " << MolName << "Range: " << Range << ", Count: " << Count << endl;
    
                        string Delimiter = "to";
                        if (Range.find(Delimiter) != string::npos) {
                            std::vector<std::string> Range_Parsed;
                            Utils::tokenize(Range, Delimiter, Range_Parsed);
                            float Range_Begin = std::stof(Range_Parsed[0]);
                            float Range_End = std::stof(Range_Parsed[1]);
                            os << "\tRange_Parsed | Range_Begin: " << Range_Parsed[0] << ", Range_End: " << Range_Parsed[1] << endl;
                            std::pair<float, float> Range(Range_Begin, Range_End);
                            std::pair<std::pair<float, float>, float> Range_Count(Range, Count);
                            Ranges.push_back(Range_Count);
    
                        } else {
                            InitialCount = Count;
                            os << "\tInitial_Count: " << InitialCount << endl;
                        }
                    } else if ((property->Key == "Fixed") & (property->Value == reactant->Id.Name)) {
                        Fixed = true;
                    }
                }

                // Add a new small molecule (small molecule is an assumption made on enzyme, but it may not always be true)
                if (Ranges.empty()) {
                    FSmallMolecule * Molecule = new FSmallMolecule(reactant->Id.Name, InitialCount, Fixed);
                    // Molecule->Print(os);
                    Context.AddToMoleculeList(Molecule);
                } else {
                    FSmallMolecule * Molecule = new FSmallMolecule(reactant->Id.Name, InitialCount, Ranges); 
                    // Molecule->Print(os);
                    Context.AddToMoleculeList(Molecule);
                }

            }

            for (const auto& product : OverallReaction.Products) {
                const std::string& ProductName = product->Id.Name;
                if (ProductName == Name) {
                    break;
                }

                Coefficient = product->Coeff; // update when coeff is fully implemented in parser
                os << "    Products: " << "(" << Coefficient << ") " << ProductName << ", " << endl;

                bool SkipStoich = false;
                for (const auto& reactant : OverallReaction.Reactants) {
                    if (ProductName == reactant->Id.Name) {
                        SkipStoich = true;
                    }
                } 
                std::pair<std::string, int> Stoich(ProductName, Coefficient);
                if (!SkipStoich) {
                    Stoichiometry.push_back(Stoich);
                }

                float InitialCount = Float_Init;
                bool Fixed = false;

                const auto& propertylist = OverallReaction.Property;
                for (auto& property :propertylist) {
                    if (property->Key == ProductName) {
                        // temporary parsing code
    
                        // Value: [0]=0, currently replaced with _0__0
                        // Value: [1:3]=2, _1to3__2
    
                        string Value = property->Value;
                        std::vector<std::string> Value_Parsed;
                        Utils::tokenize(Value, "_", Value_Parsed);
                        string MolName = Value_Parsed[0];
                        string Range = Value_Parsed[1];
                        float Count = std::stof(Value_Parsed[2]);
                        os << "\tValue_Parsed | MolName: " << MolName << "Range: " << Range << ", Count: " << Count << endl;
    
                        string Delimiter = "to";
                        if (Range.find(Delimiter) != string::npos) {
                            std::vector<std::string> Range_Parsed;
                            Utils::tokenize(Range, Delimiter, Range_Parsed);
                            float Range_Begin = std::stof(Range_Parsed[0]);
                            float Range_End = std::stof(Range_Parsed[1]);
                            os << "\tRange_Parsed | Range_Begin: " << Range_Parsed[0] << ", Range_End: " << Range_Parsed[1] << endl;
                            std::pair<float, float> Range(Range_Begin, Range_End);
                            std::pair<std::pair<float, float>, float> Range_Count(Range, Count);
                            Ranges.push_back(Range_Count);
    
                        } else {
                            InitialCount = Count;
                            os << "\tInitial_Count: " << InitialCount << endl;
                        }
                    } else if ((property->Key == "Fixed") & (property->Value == ProductName)) {
                        Fixed = true;
                    }
                }

                // Add a new small molecule (small molecule is an assumption made on enzyme, but it may not always be true)
                if (Ranges.empty()) {
                    FSmallMolecule * Molecule = new FSmallMolecule(ProductName, InitialCount, Fixed);
                    // Molecule->Print(os);
                    Context.AddToMoleculeList(Molecule);
                } else {
                    FSmallMolecule * Molecule = new FSmallMolecule(ProductName, InitialCount, Ranges); 
                    // Molecule->Print(os);
                    Context.AddToMoleculeList(Molecule);
                }
            }

			if (!Location.empty()) {
                os << "    Location: " << Location << endl;
			}

            // add new enzymatic reaction to the system
            FEnzymaticReaction *EnzymaticReaction = new FEnzymaticReaction(Name, Stoichiometry, Name);
            EnzymaticReaction->Print_Stoichiometry(os);
            Context.AddToReactionList(EnzymaticReaction);

            // add new enzyme to the system
            if (Ranges.empty()) {
                if ((Type == Standard) || (Type == MichaelisMenten)) {
                    FEnzyme * Enzyme = new FEnzyme(Type, Name, Substrate, k1, k2, InitialCount, Fixed);
                    Enzyme->Print(os);
                    Context.AddToMoleculeList(Enzyme);
                } else {
                    FEnzyme * Enzyme = new FEnzyme(Type, Name, Substrate, k1, k2, Regulator, Mode, K, n, InitialCount, Fixed);
                    Enzyme->Print(os);
                    Context.AddToMoleculeList(Enzyme);
                }
            } else {
                if ((Type == Standard) || (Type == MichaelisMenten)) {
                    FEnzyme * Enzyme = new FEnzyme(Type, Name, Substrate, k1, k2, InitialCount, Ranges);
                    Enzyme->Print(os);
                    Context.AddToMoleculeList(Enzyme);
                } else {
                    FEnzyme * Enzyme = new FEnzyme(Type, Name, Substrate, k1, k2, Regulator, Mode, K, n, InitialCount, Ranges);
                    Enzyme->Print(os);
                    Context.AddToMoleculeList(Enzyme);
                }
            }
            // TODO: if the block contains subreactions, priortize subreactions over main reaction for simulation?
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
            string Template = Context.QueryTable(Name, "Template", Context.PolymeraseTable);
            string Target = Context.QueryTable(Name, "Target", Context.PolymeraseTable);
            string Process = Context.QueryTable(Name, "Process", Context.PolymeraseTable);
            float Rate = std::stof(Context.QueryTable(Name, "Rate", Context.PolymeraseTable));

            FPolymerase * Polymerase = new FPolymerase(Name, Template, Target, Process, Rate);
            Polymerase->Print(os);
            Context.AddToMoleculeList(Polymerase);


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
                os << "    Reactants: " << "(" << Coefficient << ")" << ReactantName << ", " << endl;
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

                FSmallMolecule * Molecule = new FSmallMolecule(ReactantName);             
                // Molecule->Print(os);
                Context.AddToMoleculeList(Molecule);
            }

            for (const auto& product : ElongationReaction.Products) {
                const std::string ProductName = product->Id.Name;
                Coefficient = product->Coeff;
                os << "    Products: " << "(" << Coefficient << ")" << ProductName << ", " << endl;
                if ((ProductName == "dna_{n+1}") | (ProductName == "rna_{n+1}") | (ProductName == "peptide_{n+1}")) {
                    continue;
                }
                std::pair<std::string, int> Stoich(ProductName, Coefficient);
                Stoichiometry.push_back(Stoich);

                FSmallMolecule * Molecule = new FSmallMolecule(ProductName);
                // Molecule->Print(os);
                Context.AddToMoleculeList(Molecule);
            }

			if (!Location.empty()) {
                os << "    Location: " << Location << endl;
			}

            if (!BuildingBlocks.empty()){
                os << "    BuildingBlocks: [";
            }
            for (auto& BuildingBlock : BuildingBlocks) {
                float InitialCount_BuildingBlock = 50000;
                FSmallMolecule * Molecule = new FSmallMolecule(BuildingBlock, InitialCount_BuildingBlock);
                os << Molecule->Name << ", ";
                // Molecule->Print(os);
                Context.AddToMoleculeList(Molecule);
            }
            os << "]" << endl;
            FPolymeraseReaction *PolymeraseReaction = new FPolymeraseReaction(Name, Stoichiometry, Name, BuildingBlocks);
            // PolymeraseReaction->Print(os);
            Context.AddToReactionList(PolymeraseReaction);
    
       






        } else if (Utils::is_class_of<NOrganismDeclaration, NNode>(node)) {
            auto Organism = dynamic_cast<const NOrganismDeclaration *>(node);
            os << "Organism: " << Organism->Id.Name << endl;
            os << "  " << Organism->Description << endl;

            if ((Organism->Id.Name == "ecoli") | (Organism->Description == "E. coli K-12 MG1655")) {

                int ChromosomeSize = 4641652;
                os << "Chromosome_I: " << std::to_string(ChromosomeSize) << "bp" << endl;
                FChromosome * Chromosome = new FChromosome("ChI", ChromosomeSize);
                Context.AddToMoleculeList(Chromosome);                

                int i;
                int i_cap = 5000;

                i = 0;
                os << ">> Genes being imported... : ";
                for (auto& record : Context.GeneTable.Records) {
                    // os << record["symbol"] << ", ";
                    FGene * Gene = new FGene(record["id"], record["symbol"]);
                    Context.AddToMoleculeList(Gene);

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
                    FRNA * RNA = new FRNA(record["id"], record["type"]);
                    Context.AddToMoleculeList(RNA);

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
                    FProtein * Protein = new FProtein(record["id"]);
                    Context.AddToMoleculeList(Protein);

                    // temporary capping
                    i++;
                    if (i == i_cap) {
                        os << "Protein importing is capped at " <<  std::to_string(i_cap) << endl;
                        break;
                    }
                } os << "done" << endl;
                // os << endl;
                
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

void CheckInitialCounts()
{
    std::vector<int> Idx_Mol = Context.GetIdxListFromMoleculeList("Molecule");
    for (auto& mol : Context.MoleculeList) {
        auto& count = mol->Count;
        if (count.Initial < 0) {
            for (auto& Pathway : Context.PathwayList) {
                if (Pathway.Name == "TCA") {
                    if (!Context.QueryTable(mol->Name, "Count", Context.InitialCountTable_TCA).empty()) {
                        count.Initial = std::stof(Context.QueryTable(mol->Name, "Count", Context.InitialCountTable_TCA));
                        std::cout << "InitialCount Imported | Molecule: " << mol->Name << "\t| Count: " << count.Initial << endl;
                    }
                }
            }
        }
        if (count.Initial < 0) {
            count.Initial = 0;
        }
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

void Print_InitializeStandardReaction(ofstream& ofs, std::string Type)
{
    ofs << in+ in+ in+ "# " << Type << endl;

    // Standard vs. MichaelisMenten
    ofs << in+ in+ in+ "self.Const_k_" << Type << " = None" << endl;
    ofs << in+ in+ in+ "self.Const_krev_" << Type << " = None" << endl;
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = None" << endl;
        ofs << in+ in+ in+ "self.Idx_Product_" << i << "_" << Type << " = None" << endl;
    }
    ofs << endl;

    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ in+ "self.Const_K_" << Type << " = None" << endl;
        ofs << in+ in+ in+ "self.Const_n_" << Type << " = None" << endl;
        ofs << in+ in+ in+ "self.Idx_Regulator_" << Type << " = None" << endl;
    }
    ofs << endl;
     
    ofs << in+ in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = None" << endl;
    ofs << in+ in+ in+ "self.Const_StoichMatrix_" << Type << " = None" << endl;
    ofs << endl;
}

void Print_SetUpStandardReaction(ofstream& ofs, std::string Type)
{
    // Types allowed: "Standard_Unregulated", "Standard_Inhibited", "Standard_Activated"
    // TODO: Multiplexing for regulation
    std::string RegType;
    if      (Type == "Standard_Inhibited") { RegType = "Inhibition"; }
    else if (Type == "Standard_Activated") { RegType = "Activation"; }

    // for standard reactions 
    std::vector<float> k;
    std::vector<float> krev;

    // common variables
    int Idx_Pseudo = Context.GetIdxByName_MoleculeList("Pseudo");
    std::vector<std::vector<int>> Idx_Reactants(N_MoleculesAllowed);
    std::vector<std::vector<int>> Idx_Products(N_MoleculesAllowed);

    // for regulatory mechanisms
    std::vector<float> K;
    std::vector<float> n;
    std::vector<int> Idx_Regulator; // reactant of the regulatory reaction

    // for stoichiometry matrix
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(Type); // TODO: Add new type
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(Type); // TODO: Add new type

    std::vector<const FStandardReaction *> StandardReactionSubList = Context.GetList_Standard_ReactionList(Type);
    std::vector<const FRegulatoryReaction *> RegulatoryReactionSubList = Context.GetList_Regulatory_ReactionList(RegType);

    for (auto& reaction : StandardReactionSubList) {
        
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
    void assert (Idx_Reactants.size() == Idx_Products.size());

        for (auto& reg : RegulatoryReactionSubList) {

            bool Import = false;
            std::string reactant_reg;
            std::string product_reg;
            for (auto& stoich : reg->Stoichiometry) {
                if      (stoich.second < 0) { reactant_reg = stoich.first; }
                else if (stoich.second > 0) { product_reg  = stoich.first; }               
                // import only if the product of the regulatory reaction targets the current standard reaction
                if (stoich.first == reaction->Name) { Import = true; }
            }
            if (Import) {
                Idx_Regulator.push_back(Context.GetIdxByName_MoleculeList(reactant_reg));
                K.push_back(reg->K);
                n.push_back(reg->n);
            }
        }
    }
 
    ofs << in+ in+ in+ "# " << Type << endl;

    ofs << in+ in+ in+ "self.Const_k_" << Type << " = np.array([" << JoinFloat2Str(k) << "])" << endl;
    ofs << in+ in+ in+ "self.Const_krev_" << Type << " = np.array([" << JoinFloat2Str(krev) << "])" << endl;
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Reactants[i]) << "])" << endl;
        ofs << in+ in+ in+ "self.Idx_Product_" << i << "_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Products[i]) << "])" << endl;
    }
 
    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ in+ "self.Const_K_" << Type << " = np.array([" << JoinFloat2Str(K) << "])" << endl;
        ofs << in+ in+ in+ "self.Const_n_" << Type << " = np.array([" << JoinFloat2Str(n) << "])" << endl;
        ofs << in+ in+ in+ "self.Idx_Regulator_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Regulator) << "])" << endl;
    }
    ofs << endl;
    
    ofs << in+ in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ in+ "self.Const_StoichMatrix_" << Type << " = np.asmatrix([" << Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;
}

void Print_InitializeEnzymeReaction(ofstream& ofs, std::string Type)
{
    ofs << in+ in+ in+ "# " << Type << endl;

    ofs << in+ in+ in+ "self.Idx_Enz_" << Type << " = None" << endl;

    // Standard vs. MichaelisMenten
    if (Type.find("Standard") != string::npos) {
        ofs << in+ in+ in+ "self.Const_k_" << Type << " = None" << endl;
        ofs << in+ in+ in+ "self.Const_krev_" << Type << " = None" << endl;
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = None" << endl;
            ofs << in+ in+ in+ "self.Idx_Product_" << i << "_" << Type << " = None" << endl;
        }
 
    } else if (Type.find("MichaelisMenten") != string::npos) {
        ofs << in+ in+ in+ "self.Const_kcat_" << Type << " = None" << endl;
        ofs << in+ in+ in+ "self.Const_KM_" << Type << " = None" << endl;
        ofs << in+ in+ in+ "self.Idx_EnzSub_" << Type << " = None" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ in+ "self.Const_K_" << Type << " = None" << endl;
        ofs << in+ in+ in+ "self.Const_n_" << Type << " = None" << endl;
        if (Type.find("Inhibition") != string::npos) { 
            ofs << in+ in+ in+ "self.Idx_Inhibitor_" << Type << " = None" << endl;
        }
        if (Type.find("Activation") != string::npos) {
            ofs << in+ in+ in+ "self.Idx_Activator_" << Type << " = None" << endl;
        }
    }
    ofs << endl;
    
    ofs << in+ in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = None" << endl;
    ofs << in+ in+ in+ "self.Const_StoichMatrix_" << Type << " = None" << endl;
    ofs << endl;
}

void Print_SetUpEnzymeReaction(ofstream& ofs, std::string Type) // to be changed with reaction list
{
    // Useful copy-paste conditions
//    // Standard vs. MichaelisMenten
//    if ((Type_Enum >= 1) || (Type_Enum <= 5)) {
//    } else if ((Type_Enum >= 6) || (Type_Enum <= 10)) {
//    }

//    // Inhibition vs. Activation
//    if ((Type_Enum == 2) || (Type_Enum == 3) || (Type_Enum == 7) || (Type_Enum == 8)) {
//    } else if ((Type_Enum == 4) || (Type_Enum == 5) || (Type_Enum == 9) || (Type_Enum == 10)) {
//    }

    // for standard reactions 
    std::vector<float> k;
    std::vector<float> krev;


    // for michaelis menten kinetics
    std::vector<float> kcats;
    std::vector<float> KMs;
    std::vector<int> Idx_EnzSub;

    // common variables
    int Idx_Pseudo = Context.GetIdxByName_MoleculeList("Pseudo");
    std::vector<int> Idx_Enz; // for Enzyme where En is not included in the reaction
    std::vector<std::vector<int>> Idx_Reactants(N_MoleculesAllowed);
    std::vector<std::vector<int>> Idx_Products(N_MoleculesAllowed);

    // for regulatory mechanisms
    std::vector<float> K;
    std::vector<float> n;
    std::vector<int> Idx_Inhibitor;
    std::vector<int> Idx_Activator;

    // for stoichiometry matrix
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(Type);
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(Type);

    // loop through enzymelist and push back if k >= 0
    std::vector<const FEnzymaticReaction *> EnzymaticReactionList = Context.GetList_Enzymatic_ReactionList();
    std::vector<const FEnzyme *> EnzymeList = Context.GetList_Enzyme_MoleculeList();
    std::vector<const FEnzyme *> EnzymeList_Sub = Context.GetSubList_EnzymeList(EnzymeList, Type);

    for (auto& enzyme : EnzymeList_Sub) {
        Idx_Enz.push_back(Context.GetIdxByName_MoleculeList(enzyme->Name));

        if (Type.find("Standard") != string::npos) {

            k.push_back(enzyme->k);
            krev.push_back(enzyme->krev);

        } else if (Type.find("MichaelisMenten") != string::npos) {

            kcats.push_back(enzyme->kcat);
            KMs.push_back(enzyme->KM);
            Idx_EnzSub.push_back(Context.GetIdxByName_MoleculeList(enzyme->Substrate));
        }

        // TODO: There may be more than one Regulator for the enzyme

        if (Type.find("Inhibition") != string::npos) {

            K.push_back(enzyme->Ki);
            n.push_back(enzyme->n_i);
            Idx_Inhibitor.push_back(Context.GetIdxByName_MoleculeList(enzyme->Inhibitor));

        } else if (Type.find("Activation") != string::npos) {

            K.push_back(enzyme->Ka);
            n.push_back(enzyme->n_a);
            Idx_Activator.push_back(Context.GetIdxByName_MoleculeList(enzyme->Activator));

        }

        for (auto& reaction : EnzymaticReactionList) {
            if (enzyme->Name == reaction->Name) {
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
            }
        } 
    }
 
    // Check Idx lengths
    void assert (Idx_Reactants.size() == Idx_Products.size());

    ofs << in+ in+ in+ "# " << Type << endl;
    ofs << in+ in+ in+ "self.Idx_Enz_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Enz) << "])" << endl;

    // Standard vs. MichaelisMenten
    if (Type.find("Standard") != string::npos) {
        ofs << in+ in+ in+ "self.Const_k_" << Type << " = np.array([" << JoinFloat2Str(k) << "])" << endl;
        ofs << in+ in+ in+ "self.Const_krev_" << Type << " = np.array([" << JoinFloat2Str(krev) << "])" << endl;
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Reactants[i]) << "])" << endl;
            ofs << in+ in+ in+ "self.Idx_Product_" << i << "_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Products[i]) << "])" << endl;
        }
 
    } else if (Type.find("MichaelisMenten") != string::npos) {
        ofs << in+ in+ in+ "self.Const_kcat_" << Type << " = np.array([" << JoinFloat2Str(kcats) << "])" << endl;
        ofs << in+ in+ in+ "self.Const_KM_" << Type << " = np.array([" << JoinFloat2Str(KMs) << "])" << endl;
        ofs << in+ in+ in+ "self.Idx_EnzSub_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_EnzSub) << "])" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ in+ "self.Const_K_" << Type << " = np.array([" << JoinFloat2Str(K) << "])" << endl;
        ofs << in+ in+ in+ "self.Const_n_" << Type << " = np.array([" << JoinFloat2Str(n) << "])" << endl;
        if (Type.find("Inhibition") != string::npos) { 
            ofs << in+ in+ in+ "self.Idx_Inhibitor_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Inhibitor) << "])" << endl;
        }
        if (Type.find("Activation") != string::npos) {
            ofs << in+ in+ in+ "self.Idx_Activator_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Activator) << "])" << endl;
        }
    }
    ofs << endl;
    
    ofs << in+ in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ in+ "self.Const_StoichMatrix_" << Type << " = np.asmatrix([" << Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;
}

void Print_InitializePolymeraseReaction(ofstream& ofs, const FPolymerase* Polymerase)
{
    ofs << in+ in+ in+ "self.Freq_BB_" << Polymerase->Target << "s = None" << endl;

    // Initiation
    ofs << in+ in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ in+ "self.Idx_Template_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ in+ "self.Idx_TemplateSubset" << Polymerase->Process << " = None" << endl; // local indexing within the template population for mRNA in RNA for protein translation
    		ofs << in+ in+ in+ "self.Weight_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ in+ "self.Pol_Threshold_" << Polymerase->Process << " = None" << endl;
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
        ofs << in+ in+ in+ "self.Len_Nascent" << Polymerase->Template << "s = None" << endl;
    }

    // Elongation
    ofs << in+ in+ in+ "# " << Polymerase->Process << " Initialize ElongationReaction" << endl;
    ofs << in+ in+ in+ "self.Rate_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ in+ "self.MaxLen_Nascent" << Polymerase->Target << "s = None" << endl;
    ofs << in+ in+ in+ "self.Len_Nascent" << Polymerase->Target << "s = None" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Idx_PolSub_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ in+ "self.Idx_PolBB_" << Polymerase->Process << " = None" << endl; 
    ofs << endl;
}

void Print_SetUpPolymeraseReaction(ofstream& ofs, const FPolymerase* Polymerase, float Rate, std::string FreqBBFileName, std::string MaxLenFileName, int Idx_Pol, std::vector<int> Idx_Template, std::vector<int> Idx_TemplateSubset, std::vector<int> Idx_Target, std::vector<int> Idx_PolSub, std::vector<int> Idx_PolBB, int Threshold)
{
    ofs << in+ in+ in+ "self.Freq_BB_" << Polymerase->Target << "s = np.asmatrix(np.load(r'" << FreqBBFileName << "'))" << endl;

    // Initiation
    ofs << in+ in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = np.asmatrix([" << std::to_string(Idx_Pol) << "])" << endl;
    ofs << in+ in+ in+ "self.Idx_Template_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Template) << "])" << endl;
    ofs << in+ in+ in+ "self.Idx_TemplateSubset_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_TemplateSubset) << "])" << endl;
    ofs << in+ in+ in+ "self.Idx_Target_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Target) << "])" << endl;
    		ofs << in+ in+ in+ "self.Weight_" << Polymerase->Process << " = np.array([" << "1" << "])" << endl;
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
        ofs << in+ in+ in+ "self.Len_Nascent" << Polymerase->Template << "s = np.asmatrix(np.full([10, self.Freq_BB_" << Polymerase->Target << "s.shape[0]], -1))" << endl;
        ofs << endl;
    }

    int Template_Min, Template_Max, Target_Min, Target_Max;

    if      (Polymerase->Process == "DNAReplication") 	        {Template_Min = 1; Template_Max = 2; Target_Min = 1; Target_Max = 2;}
    else if (Polymerase->Process == "RNATranscription") 	{Template_Min = 1; Template_Max = 2; Target_Min = 0; Target_Max = 1;}
    else if (Polymerase->Process == "ProteinTranslation") 	{Template_Min = 0; Template_Max = 1; Target_Min = 0; Target_Max = 1;}
    else    				                        {Template_Min = 0; Template_Max = 1; Target_Min = 0; Target_Max = 1;} // improve exception handling
        

    ofs << in+ in+ in+ "self.Pol_Threshold_" << Polymerase->Process << " = " << std::to_string(Threshold) << endl;
    ofs << endl;
    ofs << in+ in+ in+ "Count_Template_" << Polymerase->Process << " = np.random.randint(" << std::to_string(Template_Min) << ", high=" << std::to_string(Template_Max) << ", size=self.Idx_Template_" << Polymerase->Process << ".shape[1])" << endl;
    ofs << in+ in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Template_" << Polymerase->Process << ", Count_Template_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl; 

    ofs << in+ in+ in+ "Count_Target_" << Polymerase->Process << " = np.random.randint(" << std::to_string(Target_Min) << ", high=" << std::to_string(Target_Max) << ", size=self.Idx_Target_" << Polymerase->Process << ".shape[1])" << endl;
    ofs << in+ in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Target_" << Polymerase->Process << ", Count_Target_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl; 


    // Elongation
    ofs << in+ in+ in+ "# " << Polymerase->Process << " SetUp ElongationReaction" << endl;
    ofs << in+ in+ in+ "self.Rate_" << Polymerase->Process << " = " << std::to_string(int(Rate)) << endl;
    ofs << in+ in+ in+ "self.MaxLen_Nascent" << Polymerase->Target << "s = np.asmatrix(np.load(r'" << MaxLenFileName << "'))" << endl;
    ofs << in+ in+ in+ "self.Len_Nascent" << Polymerase->Target << "s = np.asmatrix(np.full([10, self.Freq_BB_" << Polymerase->Target << "s.shape[0]], -1))" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = np.asmatrix([" << std::to_string(Context.GetIdxByName_MoleculeList(Polymerase->Name)) << "])" << endl;
    ofs << in+ in+ in+ "self.Idx_PolSub_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_PolSub) << "])" << endl;
    ofs << in+ in+ in+ "self.Idx_PolBB_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_PolBB) << "])" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "Count_Pol_" << Polymerase->Process << " = np.asmatrix(np.full(len(self.Idx_Pol_" << Polymerase->Process << "), 100))" << endl;
    ofs << in+ in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Pol_" << Polymerase->Process << ", Count_Pol_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "Count_PolSub_" << Polymerase->Process << " = np.asmatrix(np.full(len(self.Idx_PolSub_" << Polymerase->Process << "), 1000))" << endl;
    ofs << in+ in+ in+ "np.put_along_axis(self.Count_All, self.Idx_PolSub_" << Polymerase->Process << ", Count_PolSub_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "Count_PolBB_" << Polymerase->Process << " = np.random.randint(4000, high=5000, size=self.Freq_BB_" << Polymerase->Target << "s.shape[1])" << endl;
    ofs << in+ in+ in+ "np.put_along_axis(self.Count_All, self.Idx_PolBB_" << Polymerase->Process << ", Count_PolBB_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl; 

}

void Print_InitiationReaction(ofstream& ofs, const FPolymerase* Polymerase)
{
    ofs << in+ in+ in+ "# " << Polymerase->Process << endl;
    ofs << in+ in+ in+ "self.State.Len_Nascent" << Polymerase->Target << "s = self.Initiation(";
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
    ofs << in+ in+ in+ "# " << Polymerase->Process << endl;
    ofs << in+ in+ in+ "self.State.Len_Nascent" << Polymerase->Target << "s = self.Elongation(";
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
    ofs << in+ in+ in+ "# " << Polymerase->Process << endl;
    ofs << in+ in+ in+ "self.State.Len_Nascent" << Polymerase->Target << "s = self.Termination(";
                       ofs << "self.State.Len_Nascent" << Polymerase->Target << "s, ";
                       ofs << "self.State.MaxLen_Nascent" << Polymerase->Target << "s, ";
                       ofs << "self.State.Idx_Target_" << Polymerase->Process << ")" << endl;
    ofs << endl;
}


void WriteSimModule()
{
    // write simulation.py
    std::ofstream ofs(Option.SimModuleFile.c_str());
    std::string endl = "\n";

    // IMPORT
    ofs << "import os, sys" << endl;
    ofs << "import numpy as np" << endl;
    // ofs << "import tensorflow as tf" << endl;
    ofs << "from datetime import datetime" << endl;
    ofs << "import csv" << endl;
    ofs << "import SimulationFunctions as sim" << endl;
    ofs << "import plot" << endl;
    ofs << endl;

    // BODY
    ofs << "def main():   # add verbose" << endl;
    ofs << endl;

    // user input
    ofs << in+ "N_SimSteps = 1000" << endl;
    ofs << in+ "SimStepTimeResolution = 100" << endl;
    ofs << endl;

    // class FState 
    ofs << in+ "class FState:" << endl;
    ofs << in+ in+ "def __init__(self):" << endl;
    ofs << in+ in+ in+ "self.Vol = 0" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# State Arrays" << endl;
    ofs << in+ in+ in+ "self.Count_All = np.zeros([1, " << Context.MoleculeList.size() << "])" << endl;
    ofs << in+ in+ in+ "self.dCount_All = np.zeros([1, " << Context.MoleculeList.size() << "])" << endl;
    ofs << endl;

    // for standard reactions
    std::vector<const FStandardReaction *> StandardReactionList = Context.GetList_Standard_ReactionList("Standard_All");
    std::vector<std::string> StandardReactionTypes {"Standard_Unregulated", "Standard_Inhibited", "Standard_Activated"};

    if (!StandardReactionList.empty()) {      
        for (auto& Type : StandardReactionTypes) {
            std::vector<const FStandardReaction *> StandardReactionSubList = Context.GetList_Standard_ReactionList(Type);
            if (!StandardReactionSubList.empty()) {
                Print_InitializeStandardReaction(ofs, Type);
            }
        }
    }

    // for enzyme reactions
    std::vector<const FEnzyme *> EnzymeList = Context.GetList_Enzyme_MoleculeList();
    std::vector<std::string> EnzReactionTypes;

    // Get all existing enzymatic reaction types
    if (!EnzymeList.empty()) {
        for (auto& enzyme : EnzymeList) {
            if ((enzyme->k >= 0) & (enzyme->KM < 0) & (enzyme->Mode.empty())) {
                EnzReactionTypes.push_back("Standard");
            } else if ((enzyme->k >= 0) & (enzyme->KM < 0) & (!enzyme->Inhibitor.empty()) & (enzyme->Activator.empty()) & (enzyme->Mode == "Allosteric")) {
                EnzReactionTypes.push_back("Standard_Inhibition_Allosteric");
            } else if ((enzyme->k >= 0) & (enzyme->KM < 0) & (enzyme->Inhibitor.empty()) & (!enzyme->Activator.empty()) & (enzyme->Mode == "Allosteric")) {
                EnzReactionTypes.push_back("Standard_Activation_Allosteric");
            } else if ((enzyme->KM >= 0) & (enzyme->k < 0) & (enzyme->Mode.empty())) {
                EnzReactionTypes.push_back("MichaelisMenten");
            } // add more types
        }
        // remove duplicates in the reaction types
        std::sort(EnzReactionTypes.begin(), EnzReactionTypes.end());
        vector<std::string>::iterator iter;
        iter = std::unique(EnzReactionTypes.begin(), EnzReactionTypes.begin() + EnzReactionTypes.size());
        EnzReactionTypes.resize(std::distance(EnzReactionTypes.begin(), iter));
    }

    // Print InitializeEnzymeReaction for each Reaction Type
    if (!EnzymeList.empty()) {
        for (auto& Type : EnzReactionTypes) {
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





    ofs << in+ in+ "def Initialize(self):" << endl;
    ofs << in+ in+ in+ "self.Vol = 1" << endl;
    ofs << endl;

    // Print SetUpEnzymeReaction for each Reaction Type 
    if (!StandardReactionList.empty()) {      
        for (auto& Type : StandardReactionTypes) {
            std::vector<const FStandardReaction *> StandardReactionSubList = Context.GetList_Standard_ReactionList(Type);
            if (!StandardReactionSubList.empty()) {
                Print_SetUpStandardReaction(ofs, Type);
            }
        }
    }

    if (!EnzymeList.empty()) {
        for (auto& Type : EnzReactionTypes) {
            Print_SetUpEnzymeReaction(ofs, Type);
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
            std::cout << "\t | Idx_Template.size(): " << Idx_Template.size() << "\t Idx_Taret.size(): " << Idx_Target.size() << endl;
            assert (Idx_Template.size() == Idx_Target.size());   
 
            Print_SetUpPolymeraseReaction(ofs, Polymerase, Rate, FreqBBFileName, MaxLenFileName, Idx_Pol, Idx_Template, Idx_TemplateSubset, Idx_Target, Idx_PolSub, Idx_PolBB, Threshold);
        }
    }

    // Initialize all molecule counts
    std::vector<int> Idx_Mol = Context.GetIdxListFromMoleculeList("Molecule");
    std::vector<float> InitialCount_Molecules;
    for (auto& mol : Context.MoleculeList) {
        auto& count = mol->Count;
        float InitialCount = count.Initial;
        InitialCount_Molecules.push_back(InitialCount);
    }
    ofs << in+ in+ in+ "Idx_Mol = np.asmatrix([" << JoinInt2Str_Idx(Idx_Mol) << "])" << endl;
    ofs << in+ in+ in+ "Count_Mol = np.asmatrix([" << JoinFloat2Str(InitialCount_Molecules) << "])" << endl;
    ofs << in+ in+ in+ "np.put_along_axis(self.Count_All, Idx_Mol, Count_Mol, axis=1)" << endl;
    ofs << endl;

    ofs << in+ in+ "def ExportLegend(self):" << endl;
    // for legends
    std::vector<std::string> MolNames = Context.GetNames_MoleculeList();
    // std::cout << "Legend_MoleNames: " << JoinStr2Str(MolNames) << endl;
    ofs << in+ in+ in+ "return ['SimStep', 'Vol', " << JoinStr2Str(MolNames) << "]" << endl;
    ofs << endl;

    ofs << in+ in+ "def ExportData(self, Time):" << endl;
    ofs << in+ in+ in+ "Data = np.asmatrix(np.zeros(2 + " << MolNames.size() << "))" << endl;
    int i = 0;
    int i_SimStep = i + 1;
    ofs << in+ in+ in+ "Data[0, " << i << ":" << i_SimStep << "] = Time" << endl;

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
    ofs << in+ in+ in+ "self.Data = None" << endl;
    ofs << endl;

    ofs << in+ in+ "def PrintLegend(self):" << endl;
    ofs << in+ in+ in+ "print(self.Legend)" << endl;
    ofs << endl;

    ofs << in+ in+ "def PrintData(self):" << endl;
    ofs << in+ in+ in+ "print(self.Data)" << endl;
    ofs << endl;

    // class FSimulation
    ofs << in+ "class FSimulation:" << endl;
    ofs << in+ in+ "def __init__(self, InState, InDataset, InDataManager):" << endl;
    ofs << in+ in+ in+ "self.N_SimSteps = 0" << endl;
    ofs << in+ in+ in+ "self.SimStep = 0" << endl;
    ofs << in+ in+ in+ "self.SimTimeResolutionPerSecond = 0" << endl;
    ofs << in+ in+ in+ "self.State = InState" << endl;
    ofs << in+ in+ in+ "self.Dataset = InDataset" << endl;
    ofs << in+ in+ in+ "self.DataManager = InDataManager" << endl;

    // Restore
    for (auto& molecule : Context.MoleculeList) {
        auto& count = molecule->Count;
        if (count.Fixed) {
            ofs << in+ in+ in+ "self.Idx_Restore_" << molecule->Name << " = None" << endl;
        }
    }

    // Event
    for (auto& molecule : Context.MoleculeList) {
        auto& count = molecule->Count;
        for (auto& range : count.Range) { 
            if ((range.first.first >= 0) & (range.first.second >= 0)) {
                ofs << in+ in+ in+ "self.Idx_Event_" << molecule->Name << " = None" << endl;
                break;
            }
        }
    }

//    if (!Context.PathwayList.empty()){
//        for (auto& Pathway : Context.PathwayList) {
//            if (Pathway.Name == "TCA") {
//                ofs << in+ in+ in+ "self.Idx_Restore_" << Pathway.Name << " = None" << endl;
//            }
//        }
//    }
    ofs << endl;

    ofs << in+ in+ "def Initialize(self, InN_SimSteps, InTimeResolution):" << endl;
    ofs << in+ in+ in+ "print('Simulation Initialized...')" << endl;
    ofs << in+ in+ in+ "self.N_SimSteps = np.asmatrix([InN_SimSteps])" << endl;
    ofs << in+ in+ in+ "self.SimTimeResolutionPerSecond = InTimeResolution" << endl;
    ofs << endl;

    // Restore
    for (auto& molecule : Context.MoleculeList) {
        auto& count = molecule->Count;
        if (count.Fixed) {
            int Idx = Context.GetIdxByName_MoleculeList(molecule->Name);
            ofs << in+ in+ in+ "self.Idx_Restore_" << molecule->Name << " = np.asmatrix([" << std::to_string(Idx) << "])" << endl;
        }
    }

    // Event
    for (auto& molecule : Context.MoleculeList) {
        auto& count = molecule->Count;
        for (auto& range : count.Range) { 
            if ((range.first.first >= 0) & (range.first.second >= 0)) {
                int Idx = Context.GetIdxByName_MoleculeList(molecule->Name);
                ofs << in+ in+ in+ "self.Idx_Event_" << molecule->Name << " = np.asmatrix([" << std::to_string(Idx) << "])" << endl;
                break;
            }
        }
    }


//    if (!Context.PathwayList.empty()){
//        for (auto& Pathway : Context.PathwayList) {
//            if (Pathway.Name == "TCA") {
//                std::string MoleculeToRestore;
//                int Idx;
//
//                MoleculeToRestore = "acetyl-CoA";
//                Idx = Context.GetIdxByName_MoleculeList(MoleculeToRestore);
//                ofs << in+ in+ in+ "self.Idx_Restore_" << Pathway.Name << " = np.asmatrix([" << Idx << "]) # " << MoleculeToRestore << endl;
//
//	                MoleculeToRestore.clear();
//                MoleculeToRestore = "malate";
//                Idx = Context.GetIdxByName_MoleculeList(MoleculeToRestore);
//                ofs << in+ in+ in+ "self.Idx_Restore_" << Pathway.Name << " = np.asmatrix([" << Idx << "]) # " << MoleculeToRestore << endl;
//            }
//        }
//    }
    ofs << in+ in+ in+ "self.State.Initialize()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Legend Export" << endl;
    ofs << in+ in+ in+ "self.Dataset.Legend = self.State.ExportLegend()" << endl;
    ofs << in+ in+ in+ "self.DataManager.SetLegend(self.Dataset.Legend)" << endl;

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

    if (!StandardReactionList.empty()){
        ofs << in+ in+ in+ in+ "self.StandardReactions()" << endl;
    }

    if (!EnzymeList.empty()){
        ofs << in+ in+ in+ in+ "self.EnzymaticReactions()" << endl;
    }

    // TODO: encapsulate this part with each polymerase to allow more process-specific customization
    if (!PolymeraseList.empty()){
        ofs << in+ in+ in+ in+ "self.InitiationReactions()" << endl;
        ofs << in+ in+ in+ in+ "self.ElongationReactions()" << endl;
        ofs << in+ in+ in+ in+ "self.TerminationReactions()" << endl;
    }
                          
    ofs << in+ in+ in+ in+ "# Update Substrate Count" << endl;
    ofs << in+ in+ in+ in+ "self.State.Count_All += self.State.dCount_All" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "# Restore Substrate Count for Sustained Substrate Influx" << endl;
    ofs << in+ in+ in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "# Trigger Event on Substrate Count" << endl;
    ofs << in+ in+ in+ in+ "self.TriggerEventMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "print('Simulation Run Completed')" << endl;
    ofs << endl;

    ofs << in+ in+ "def ExportData(self):" << endl;    
    ofs << in+ in+ in+ "self.Dataset.Data = self.State.ExportData(self.SimStep/self.SimTimeResolutionPerSecond)" << endl;
    ofs << in+ in+ in+ "self.DataManager.Add(self.Dataset.Data)" << endl;
    ofs << endl;

    ofs << in+ in+ "def ApplySimTimeResolution(self, Rate):" << endl;    
    ofs << in+ in+ in+ "return Rate / self.SimTimeResolutionPerSecond" << endl;
    ofs << endl;

    // Restore
    ofs << in+ in+ "def RestoreMoleculeCount(self):" << endl;    
    bool Pass = true;

    for (auto& molecule : Context.MoleculeList) {
        auto& count = molecule->Count;
        if (count.Fixed) {
            ofs << in+ in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Restore_" << molecule->Name << ", ";
            if (molecule->Name == "Pseudo") {
                ofs << "self.State.Vol";
            } else {
                ofs << std::to_string(count.Initial);
            }
            ofs << ", axis=1)" << endl;
            Pass = false;
        }
    }
//    if (!Context.PathwayList.empty()){
//        for (auto& Pathway : Context.PathwayList) {
//            if (Pathway.Name == "TCA") {
//                std::string MoleculeToRestore = "acetyl-CoA";
//                int Count = 446331; 
//                ofs << in+ in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Restore_" << Pathway.Name << ", " << std::to_string(Count) << ", axis=1)   # " << MoleculeToRestore << endl;
//                Pass = false;
//            }
//        }
//    }
    if (Pass) {
        ofs << in+ in+ in+ "pass" << endl;
    }
    ofs << endl;

    // Event
    ofs << in+ in+ "def TriggerEventMoleculeCount(self):" << endl;
    ofs << in+ in+ in+ "Time = self.SimStep / self.SimTimeResolutionPerSecond" << endl;
    bool ElseSwitch = false;

    for (auto& molecule : Context.MoleculeList) {
        auto& count = molecule->Count;
        int i = 0;
        for (auto& range : count.Range) { 
            if ((range.first.first >= 0) & (range.first.second >= 0)) {
                ofs << in+ in+ in;
                if (i != 0) {
                    ofs << "el";
                }
                ofs << "if (Time >= " << std::to_string(range.first.first) << ") & (Time < " << std::to_string(range.first.second) << "):" << endl;
                ofs << in+ in+ in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Event_" << molecule->Name << ", " << std::to_string(range.second) << ", axis=1)" << endl;
                ElseSwitch = true;
                i++;
            }
        }
        if (ElseSwitch) {
            ofs << in+ in+ in+ "else:" << endl;
            ofs << in+ in+ in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Event_" << molecule->Name << ", " << "0"  << ", axis=1)" << endl;
            ElseSwitch = false;  
        }
    }
    ofs << endl;

    ofs << in+ in+ "# Biochemical Reaction related routines" << endl;
// here
    if (!StandardReactionList.empty()) {      
        for (auto& Type : StandardReactionTypes) {
            std::vector<const FStandardReaction *> StandardReactionSubList = Context.GetList_Standard_ReactionList(Type);
            if (!StandardReactionSubList.empty()) {

                ofs << in+ in+ "def StandardReaction_" << Type << "(self):" << endl;
                // Get Concentrations
                // Reactants and products or EnzSubstrate
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << in+ in+ in+ "Conc_Reactant_" << i << " = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Reactant_" << i << "_" << Type << "), self.State.Vol)" << endl;
                    ofs << in+ in+ in+ "Conc_Product_" << i << " = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Product_" << i << "_" << Type << "), self.State.Vol)" << endl;
                } 
    
                // Regulators 
                if ((Type.find("Inhibited") != std::string::npos) || (Type.find("Activated") != std::string::npos)) { 
                    ofs << in+ in+ in+ "Conc_Regulator = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Regulator_" << Type << "), self.State.Vol)" << endl;
                } 
                // Calculate Rate
                // Eqn type and Enz
                ofs << in+ in+ in+ "Rate = sim.Eqn_" << Type << "(";
    
                // Reactants and products or EnzSubstrate
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << "Conc_Reactant_" << i << ", "; 
                }
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << "Conc_Product_" << i << ", "; 
                }
    
                // Regulators
                if ((Type.find("Inhibited") != std::string::npos) || (Type.find("Activated") != std::string::npos)) { 
                    ofs << in+ in+ in+ "Conc_Regulator, ";
                }
                // Reaction constants
                ofs << "self.State.Const_k_" << Type << ", self.State.Const_krev_" << Type;
                if ((Type.find("Inhibited") != std::string::npos) || (Type.find("Activated") != std::string::npos)) { 
                    ofs << ", self.State.Const_K_" << Type << ", ";             
                    ofs << "self.State.Const_n_" << Type;
                }
                ofs << ")" << endl;
    
                // Apply TimeResolution to the rate
                ofs << in+ in+ in+ "Rate = self.ApplySimTimeResolution(Rate)" << endl;
    
                // Apply stoichiometry 
                ofs << in+ in+ in+ "dConc_Mol_InStoichMatrix = sim.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
                // Convert to counts
                ofs << in+ in+ in+ "dCount_Mol_InStoichMatrix = sim.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
                // Apply delta counts for molecules in the stoichiometry matrix
                ofs << in+ in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
                ofs << endl;
            }
        }
    }

    if (!EnzymeList.empty()) {
        for (auto& Type : EnzReactionTypes) {
            ofs << in+ in+ "def EnzymaticReaction_" << Type << "(self):" << endl;
            // Get Concentrations
            // Enzyme
             ofs << in+ in+ in+ "Conc_Enz = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Enz_" << Type << "), self.State.Vol)" << endl;
            // Reactants and products or EnzSubstrate
            if (Type.find("Standard") != std::string::npos) {
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << in+ in+ in+ "Conc_Reactant_" << i << " = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Reactant_" << i << "_" << Type << "), self.State.Vol)" << endl;
                    ofs << in+ in+ in+ "Conc_Product_" << i << " = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Product_" << i << "_" << Type << "), self.State.Vol)" << endl;
                } 
            } else if (Type.find("MichaelisMenten") != std::string::npos) {
                ofs << in+ in+ in+ "Conc_EnzSub = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_EnzSub_" << Type << "), self.State.Vol)" << endl;
            }

            // Regulators 
            if (Type.find("Inhibition") != std::string::npos) { 
                ofs << in+ in+ in+ "Conc_Inhibitor = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Inhibitor_" << Type << "), self.State.Vol)" << endl;
            } else if (Type.find("Activation") != std::string::npos) {
                ofs << in+ in+ in+ "Conc_Activator = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Activator_" << Type << "), self.State.Vol)" << endl;
            }

            // Calculate Rate
            // Eqn type and Enz
            ofs << in+ in+ in+ "Rate = sim.Eqn_Enz_" << Type << "(Conc_Enz";

            // Reactants and products or EnzSubstrate
            if (Type.find("Standard") != std::string::npos) {
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << ", Conc_Reactant_" << i; 
                }
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << ", Conc_Product_" << i; 
                }
            } else if (Type.find("MichaelisMenten") != std::string::npos) {
                ofs << ", Conc_EnzSub";
            }

            // Regulators
            if (Type.find("Inhibition") != std::string::npos) { 
                ofs << in+ in+ in+ ", Conc_Inhibitor";
            } else if (Type.find("Activation") != std::string::npos) {
                ofs << in+ in+ in+ ", Conc_Activator";
            }
            // Reaction constants
            if (Type.find("Standard") != std::string::npos) {
                ofs << ", self.State.Const_k_" << Type << ", self.State.Const_krev_" << Type;
            } else if (Type.find("MichaelisMenten") != std::string::npos) {
                ofs << ", self.State.Const_kcat_" << Type << ", self.State.Const_KM_" << Type; 
            }
            if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                ofs << ", self.State.Const_K_" << Type;             
                if (Type.find("Allosteric") != std::string::npos) {
                    ofs << ", self.State.Const_n_" << Type;
                }
            }
            ofs << ")" << endl;

            // Apply TimeResolution to the rate
            ofs << in+ in+ in+ "Rate = self.ApplySimTimeResolution(Rate)" << endl;

            // Apply stoichiometry 
            ofs << in+ in+ in+ "dConc_Mol_InStoichMatrix = sim.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
            // Convert to counts
            ofs << in+ in+ in+ "dCount_Mol_InStoichMatrix = sim.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
            // Apply delta counts for molecules in the stoichiometry matrix
            ofs << in+ in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
            ofs << endl;
        }
    }

    ofs << in+ in+ "def StandardReactions(self):" << endl;
    // Print SetUpEnzymeReaction for each Reaction Type 
    if (!StandardReactionList.empty()) {
        for (auto& Type : StandardReactionTypes) {
            ofs << in+ in+ in+ "self.StandardReaction_" << Type << "()" << endl;
        }
    } 
    ofs << endl;
    

    ofs << in+ in+ "def EnzymaticReactions(self):" << endl;
    // Print SetUpEnzymeReaction for each Reaction Type 
    if (!EnzymeList.empty()) {
        for (auto& Type : EnzReactionTypes) {
            ofs << in+ in+ in+ "self.EnzymaticReaction_" << Type << "()" << endl;
        }
    } 
    ofs << endl;

    if (!PolymeraseList.empty()) {

        ofs << in+ in+ "def InitiationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {
            
            Print_InitiationReaction(ofs, Polymerase);
        }


        ofs << in+ in+ "def ElongationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {
            
            Print_ElongationReaction(ofs, Polymerase);
        }

        ofs << in+ in+ "def TerminationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {
            
            Print_TerminationReaction(ofs, Polymerase);
        }
    } 

    if (EnzymeList.empty() & PolymeraseList.empty()) {
            ofs << in+ in+ in+ "pass" << endl;
            ofs << endl;
    }
    
    ofs << in+ in+ "# Useful routines" << endl;
    ofs << in+ in+ "def GetCount(self, Idx):" << endl;
    ofs << in+ in+ in+ "return np.take(self.State.Count_All, Idx)" << endl;
    ofs << endl;

    ofs << in+ in+ "def AddTodCount(self, Idx, Values):" << endl;
    ofs << in+ in+ in+ "dCountToAdd = np.zeros_like(self.State.dCount_All)" << endl;
    ofs << in+ in+ in+ "np.put_along_axis(dCountToAdd, Idx, Values, axis=1)" << endl;
    ofs << in+ in+ in+ "dCount_All_New = self.State.dCount_All + dCountToAdd" << endl;
    ofs << in+ in+ in+ "ZeroTest = dCount_All_New + self.State.Count_All" << endl;
    ofs << in+ in+ in+ "self.State.dCount_All =  np.where(ZeroTest < 0, dCount_All_New - ZeroTest, dCount_All_New)" << endl;
    ofs << endl;

    ofs << in+ in+ "def OverElongationCorrection(self, Len_Elongated, Max):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ in+ "Len_Over = np.where(Len_Elongated > Max, Len_Elongated - Max, 0)" << endl;
    ofs << in+ in+ in+ "return Len_Elongated - Len_Over" << endl;
    ofs << endl;

    ofs << in+ in+ "def BuildingBlockConsumption(self, Freq, N_Elongated_PerSpecies):" << endl;
    ofs << in+ in+ in+ "Raw = sim.DetermineAmountOfBuildingBlocks(Freq, N_Elongated_PerSpecies)" << endl;
    ofs << in+ in+ in+ "Rounded = np.around(Raw)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Discrepancy handling" << endl;
    ofs << in+ in+ in+ "N_Elongated = np.sum(N_Elongated_PerSpecies)" << endl;
    ofs << in+ in+ in+ "Discrepancy = np.sum(Rounded) - N_Elongated" << endl;

    ofs << in+ in+ in+ "NUniq_BuildingBlocks = Freq.shape[1]" << endl;
    ofs << in+ in+ in+ "Sets, Remainder = np.divmod(Discrepancy, NUniq_BuildingBlocks)" << endl;
    ofs << in+ in+ in+ "return Rounded + np.ones(NUniq_BuildingBlocks) * np.int32(Sets) + np.concatenate((np.ones(np.int32(np.round(Remainder))), np.zeros(np.int32(np.around(NUniq_BuildingBlocks - Remainder)))))" << endl;
    ofs << endl;

    ofs << in+ in+ "# Polymerase Reaction related" << endl;
    ofs << in+ in+ "def Initiation(self, Len_Template, Len_Target, Idx_Pol, Idx_Template, Idx_TemplateSubset, Weight, PolThreshold):" << endl;
    ofs << in+ in+ in+ "# Get available, active polymerase count - TO BE UPDATED with more regulatory algorithms" << endl;
    ofs << in+ in+ in+ "Count_Pol = self.GetCount(Idx_Pol)" << endl;
    ofs << in+ in+ in+ "Count_Pol_Active = np.floor_divide(Count_Pol, 2).astype(int)" << endl;
    ofs << in+ in+ in+ "Count_Pol_Occupied = np.count_nonzero(np.where(Len_Target != -1, 1, 0)) * PolThreshold" << endl;
    ofs << in+ in+ in+ "Count_Pol_Avail = Count_Pol_Active - Count_Pol_Occupied" << endl;
    ofs << in+ in+ in+ "Count_Pol_FunctionalUnit = np.floor_divide(Count_Pol_Avail, PolThreshold)" << endl;
    ofs << in+ in+ in+ "Count_Pol_Avail = np.where(Count_Pol_FunctionalUnit > 0, Count_Pol_FunctionalUnit, 0)[0, 0]" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Get final initiation weight by applying initiation site count" << endl;
    ofs << in+ in+ in+ "Count_Template_Complete = self.GetCount(Idx_Template)" << endl;
    ofs << in+ in+ in+ "Count_Template_Nascent = np.count_nonzero(np.where(Len_Template != -1, 1, 0), axis=0)   # Assumption: each nascent template has one highly efficient initiation site" << endl;
    ofs << in+ in+ in+ "Count_TemplateSubset_Nascent = np.take(Count_Template_Nascent, Idx_TemplateSubset)" << endl;
    ofs << in+ in+ in+ "Count_InitiationSite = Count_Template_Complete + Count_TemplateSubset_Nascent" << endl;
    ofs << in+ in+ in+ "Weight_Initiation = Count_InitiationSite * Weight " << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Get randomly selected target indices" << endl;
    ofs << in+ in+ in+ "Idx_Selected = sim.PickRandomIdx(Count_Pol_Avail, Idx_Template, Weight_Initiation)" << endl;
    ofs << in+ in+ in+ "Len_Target_Initiated = sim.InsertZeroIntoNegOneElementInLenMatrix(Len_Target, Idx_Selected)" << endl;
    ofs << in+ in+ in+ "# Export Data" << endl;
    ofs << in+ in+ in+ "# N_Initiated" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "return Len_Target_Initiated" << endl;
    ofs << endl;
 
    ofs << in+ in+ "def Elongation(self, Len, Max, Rate, Weight, Freq, Idx_PolSub, Idx_BB):" << endl;
    ofs << in+ in+ in+ "NUniq_BuildingBlocks = Freq.shape[1]" << endl;
    ofs << in+ in+ in+ "NUniq_Species = Freq.shape[0]" << endl;
    ofs << endl;

//    ofs << in+ in+ in+ "dLength = np.matmul(SMatrix,Rate)
    ofs << in+ in+ in+ "dLength = self.ApplySimTimeResolution(Rate)   # this is not necessarily true based on the reaction input" << endl;
    ofs << in+ in+ in+ "Len_Elongated = np.where(Len >= 0, Len + dLength, Len)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "Len_Trimmed = self.OverElongationCorrection(Len_Elongated, Max)" << endl;

    ofs << in+ in+ in+ "N_Elongated_PerSpecies = np.asmatrix(np.sum(Len_Trimmed - Len, axis=0))   # This step loses shape for some reason, hence apply matrix again" << endl;
    ofs << in+ in+ in+ "N_Elongated = np.sum(N_Elongated_PerSpecies)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "Consumed_BB = self.BuildingBlockConsumption(Freq, N_Elongated_PerSpecies)" << endl;
   
    ofs << in+ in+ in+ "# Update dCount for BuildingBlocks" << endl;
    ofs << in+ in+ in+ "self.AddTodCount(Idx_BB, -Consumed_BB)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Update dCount for Polymerase Reaction Substrates (To be updated by the reaction matrix form" << endl;
    ofs << in+ in+ in+ "self.AddTodCount(Idx_PolSub, N_Elongated)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Export Data" << endl;
    ofs << in+ in+ in+ "# N_Elongated" << endl;
    
    ofs << in+ in+ in+ "return Len_Trimmed" << endl;
    ofs << endl; 

    ofs << in+ in+ "def Termination(self, Len, Max, Idx_Target):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ in+ "Bool_Completed = (Len == Max)" << endl;
    ofs << in+ in+ in+ "N_Completed_PerSpecies = np.sum(Bool_Completed, axis=0)" << endl;
    ofs << in+ in+ in+ "N_Completed = np.sum(N_Completed_PerSpecies)" << endl;
    ofs << in+ in+ in+ "Len_Completed = np.where(Bool_Completed, -1, Len)" << endl;

    ofs << in+ in+ in+ "# Update dCount for BuildingBlocks" << endl;
    ofs << in+ in+ in+ "self.AddTodCount(Idx_Target, N_Completed_PerSpecies)" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Export Data" << endl;
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
    ofs << in+ "DataManager = FDataManager()" << endl;
    ofs << in+ "Simulation = FSimulation(State, Data, DataManager)" << endl;
    ofs << endl;

    // Simulation Module
    ofs << in+ "Simulation.Initialize(N_SimSteps, SimStepTimeResolution)" << endl;
    ofs << in+ "Simulation.Run()" << endl;
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

    ostream& os = std::cout;

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
//            Context.PrintInitialCounts(os);
            CheckInitialCounts();

            if (Option.Verbose) {
                Context.PrintLists(os);
                //Context.PrintInitialCounts(os);
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
        // cout << "## Simulation_Python ##" << endl;
     
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
