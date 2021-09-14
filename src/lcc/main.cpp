#include <iostream>
#include "node.h"

extern NBlock* ProgramBlock;
extern int yyparse();
extern int yylex_destroy();
extern FILE* yyin;

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

    std::cout << ProgramBlock << std::endl;

    DumpNBlock(ProgramBlock);


    delete ProgramBlock;
    return 0;
}
