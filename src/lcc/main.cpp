#include <iostream>
#include <queue>
#include <cassert>

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
void DumpNBlock(const NBlock* InProgramBlock) {
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

    if (!Option.bParseOnly) {
        string UsingModuleFilename;
        if (!Option.OutputPrefix.empty()) {
            UsingModuleFilename = Option.OutputPrefix + "/";
            if (!Utils::CreatePaths(UsingModuleFilename.c_str())) {
                cerr << "Can't create directory" << endl;
            }
        }
        UsingModuleFilename += "process_module.tsv";

        Context.SaveUsingModuleList(UsingModuleFilename.c_str());
    }

    // temporary simulation code
    cout << endl << "## Simulation ##" << endl;

    Simulation.Init(State, Dataset, 20);
    Simulation.Run(State, Context, Dataset);

	if (!Option.SimResultFile.empty()) {
		DataManager.SaveToFile(Option.SimResultFile.c_str());
	}
    return 0;
}
