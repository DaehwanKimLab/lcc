#include <iostream>
#include <queue>
#include "node.h"

extern NBlock* ProgramBlock;
extern int yyparse();
extern int yylex_destroy();
extern FILE* yyin;


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
    FTraversalContext Context(std::cerr);
    Context.Queue.push(InProgramBlock);

    while(!Context.Queue.empty()) {
        const NNode* node = Context.Queue.front(); Context.Queue.pop();

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
    if (argc < 2) {
        std::cerr << argv[0] << " LPPSourceFile" << std::endl;
        return -1;
    }

    yyin = fopen(argv[1], "r");
    yyparse();
    fclose(yyin);
    yylex_destroy();

//    std::cout << ProgramBlock << std::endl;

//    DumpNBlock(ProgramBlock);

    TraversalNode(ProgramBlock);

    delete ProgramBlock;
    return 0;
}
