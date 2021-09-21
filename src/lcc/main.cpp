#include <iostream>
#include <queue>
#include <cassert>

#include <unistd.h>
#include <getopt.h>

#include "node.h"
#include "option.h"
#include "context.h"

extern NBlock* ProgramBlock;
extern int yyparse();
extern int yylex_destroy();
extern FILE* yyin;

using namespace std;

FOption Option;
FCompilerContext Context;


class FExperiment {
public:
    std::string Name;
    std::string Description;

};

class FOrganism {
public:
    std::string Name;
    std::string Description;

};

class FMolecule {
public:
    std::string Name;
    std::string Id;

};

class FPathway {
public:

};

class FProtein {
public:

};

class FReaction {
public:

};

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


template <typename Derived>
static bool is_class_of(const NNode *Node) {
    Derived* DerivedNode = dynamic_cast<Derived *>(const_cast<NNode *>(Node));
    return DerivedNode != nullptr;
}

void TraversalNode(NBlock* InProgramBlock)
{
    ostream& os = std::cout;
    FTraversalContext Context(std::cerr);
    Context.Queue.push(InProgramBlock);

    while(!Context.Queue.empty()) {
        const NNode* node = Context.Queue.front(); Context.Queue.pop();

        if (is_class_of<NProteinDeclaration>(node)) {
            auto* Protein = dynamic_cast<const NProteinDeclaration *>(node);
            os << "Protein: " << Protein->Id.Name << endl;

        } else if (is_class_of<NPathwayDeclaration>(node)) {
            auto* Pathway = dynamic_cast<const NPathwayDeclaration *>(node);
            os << "Pathway: " << Pathway->Id.Name << endl;

        } else if (is_class_of<NOrganismDeclaration>(node)) {
            auto* Organism = dynamic_cast<const NOrganismDeclaration *>(node);
            os << "Organism: " << Organism->Id.Name << endl;
            os << "  " << Organism->Description << endl;
        } else if (is_class_of<NExperimentDeclaration>(node)) {
            auto* Experiment = dynamic_cast<const NExperimentDeclaration *>(node);
            os << "Experiment: " << Experiment->Id.Name << endl;
            os << "  " << Experiment->Description << endl;
            if (Experiment->Block) {
                for(const auto& stmt : Experiment->Block->Statements) {
                    os << "  "; stmt->Print(os); os << endl;
                }
            }
        }

        node->Visit(Context);
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
#if 0
    if (argc < 2) {
        std::cerr << argv[0] << " LPPSourceFile" << std::endl;
        return -1;
    }
#endif

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

        delete ProgramBlock;
    }



    return 0;
}
