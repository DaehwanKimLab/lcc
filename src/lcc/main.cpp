#include <iostream>
#include "node.h"

extern NBlock* ProgramBlock;
extern int yyparse();

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
    yyparse();

    std::cout << ProgramBlock << std::endl;

    DumpNBlock(ProgramBlock);

    return 0;
}
