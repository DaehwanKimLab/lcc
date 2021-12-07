#include <iostream>
#include <queue>
#include <cassert>

#include <unistd.h>
#include <getopt.h>

#include "node.h"
#include "option.h"
#include "context.h"
#include "util.h"

extern NBlock* ProgramBlock;
extern int yyparse();
extern int yylex_destroy();
extern FILE* yyin;

using namespace std;

FOption Option;
FCompilerContext Context;


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

class Database {


};

class NCompilerData {
public:
    std::vector<string> EnzymeList;
    std::vector<string> SubstrateList;

    const float CellVolume = 0.005; // 7e-16;
    const float MolWeight = 100.0;
    const float kcat = 5.0;
    const float kM = 0.004;
    int Count_Enzymes = 50;
    int Count_Substrates = 100000;
//    std::vector<float> MolWeight;
//    std::vector<float> kcat;
//    std::vector<float> kM;
//    std::vector<int> Count_Enzymes;
//    std::vector<int> Count_Substrates;
    
    NCompilerData() {}
    
    void Print() const {
        cout << "CompilerData(" << std::endl;

        cout << "  EnzymeList(" << std::endl;
        for (auto& enzymeName: EnzymeList) {
            cout << "  " << enzymeName << ", ";
        } cout << "  )" << std::endl; 

        cout << "  SubstrateList(" << std::endl;
        for (auto& substrateName: SubstrateList) {
            cout << "  " << substrateName << ", ";
        } cout << "  )" << std::endl; 
    }
   
    void Print_SystemState() const {
        cout << "System State(" << std::endl;
        cout << "  CellVolume:" << CellVolume << std::endl;
        cout << "  MolWeight:" << MolWeight << std::endl;
        cout << "  kcat:" << kcat << std::endl;
        cout << "  kM:" << kM << std::endl;
        cout << "  Count_Enzymes:" << Count_Enzymes << std::endl;
        cout << "  Count_Substrates:" << Count_Substrates << std::endl;
        cout << "  )" << std::endl; 
    }
};

float ConcentrationToCount(float Conc_Molecule, float Volume) {
    return Conc_Molecule * Volume;
};

float CountToConcentration(int Count_Molecule, float Volume) {
    return Count_Molecule / Volume;
};

float MichaelisMentenEqn(int Count_Enzyme, int Count_Substrate, float CellVolume ,float kcat, float kM) {
    float Conc_Enzyme = CountToConcentration(Count_Enzyme, CellVolume);
    float Conc_Substrate = CountToConcentration(Count_Substrate, CellVolume);

    float Rate = (kcat * Conc_Enzyme * Conc_Substrate) / (Conc_Substrate + kM);
    cout << Rate << std::endl;
    return Rate;
}; 

//void AddToList(vector<string> List, string Item) {
//    if (std::find(std::begin(List), std::end(List), Item) != std::end(List)) {
//        IdList.push_back(Item);
//    } else {
//        return; 
//
//    }  
//};

void TraversalNode(NBlock* InProgramBlock)
{
    ostream& os = std::cout;
    FTraversalContext tc(std::cerr);
    tc.Queue.push(InProgramBlock);

    while(!tc.Queue.empty()) {
        const NNode* node = tc.Queue.front(); tc.Queue.pop();

        if (Utils::is_class_of<NProteinDeclaration, NNode>(node)) {
            auto Protein = dynamic_cast<const NProteinDeclaration *>(node);
            os << "Protein Id: " << Protein->Id.Name << endl;
            // Protein->Print(os);

            auto& Id = Protein->Id;
	    os << "  Id: " << Id.Name << endl;
	    
            string Name = Id.Name;
            string Substrate = Context.QueryEnzymeTable(Name, "Substrate");
            float kcat = std::stof(Context.QueryEnzymeTable(Name, "kcat"));
            float kM = std::stof(Context.QueryEnzymeTable(Name, "kM"));

            os << "Query Results: " << Substrate << ", " << kcat << ", " << kM << endl;

            FEnzyme Enzyme(Name, Substrate, kcat, kM);

            auto& OverallReaction = Protein->OverallReaction;
            // os << "  OverallReaction:" << endl;

            for (const auto& reactant : OverallReaction.Reactants) {
                os << "  Reactant: " << reactant->Name << ", " << endl;
                // CompilerData->SubstrateList.push_back(reactant->Name);
            }

//            os << "  Products: " << endl;
//            os << "    ";

//            for (const auto& product : OverallReaction.Products) {
//                SubstrateList.push_back(product);
//                for (auto& substrate: SubstrateList) {
//                    os << substrate->Name << ", ";
//                }

//                product->Print(os); os << ", ";
//            }
//            os << std::endl;

//            os << "-----" << endl;
//
//            if (Protein->Block) {
//                auto& Block = Protein->Block;
//                for (auto& stmt: Block->Statements) {
//                    os << "  "; stmt->Print(os);
//                }
//
//            }
//            Context.ProteinList.emplace_back(*Protein);
            Context.EnzymeList.emplace_back(Enzyme);
 
        } else if (Utils::is_class_of<NPathwayDeclaration, NNode>(node)) {
            auto Pathway = dynamic_cast<const NPathwayDeclaration *>(node);
            os << "Pathway: " << Pathway->Id.Name << endl;

            string Name = Pathway->Id.Name;
            vector<string> Sequence;

            if (Pathway->PathwayChainReaction) {
                auto& PathwayChainReaction = Pathway->PathwayChainReaction;
                auto& Exprs = PathwayChainReaction->Exprs;
                for (auto& expr: Exprs) {
//                    os << "  "; expr->Print(os);
                    auto& Identifiers = expr->Identifiers;
                    for (auto& Id: Identifiers) {
                        os << "  Enzyme: " << Id.Name << endl;
                        Sequence.push_back(Id.Name);
                    }  
                }
            }

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

    // Load genes.tsv
    if (!Option.bParseOnly)
    {
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

    // NCompilerData CompilerData;

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

        // organize context database (protein list, pathway list)
        TraversalNode(ProgramBlock);

        delete ProgramBlock;

        // CompilerData.Print_SystemState();
        // CompilerData.Print();
        
        // MichaelisMentenEqn(CompilerData.Count_Enzymes, CompilerData.Count_Substrates, CompilerData.CellVolume, CompilerData.kcat, CompilerData.kM);
        

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

    return 0;
}
