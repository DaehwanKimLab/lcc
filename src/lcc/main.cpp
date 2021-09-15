#include <iostream>
#include <queue>
#include <cassert>

#include <unistd.h>
#include <getopt.h>

#include "node.h"

extern NBlock* ProgramBlock;
extern int yyparse();
extern int yylex_destroy();
extern FILE* yyin;

enum {
    ARG_START = 256,
    ARG_HELP,
    ARG_VERBOSE,
    ARG_DEBUG,
    ARG_LIBRARY_PATH,
    ARG_DATA_PATH,
    ARG_OUTPUT_PREFIX,
};


class FOption {
public:
    FOption() {
        Reset();
    }

    int Parse(int argc, char *argv[])
    {
        static const char *ShortOptions = "hvL:i:o:";
        static const struct option LongOptions[] = {
                {"help", no_argument, NULL, 'h'},
                {"verbose", no_argument, NULL, 'v'},
                {"datadir", required_argument, NULL, 'i'},
                {"output", required_argument, NULL, 'o'},
                {"debug", no_argument, NULL, ARG_DEBUG},
                {NULL, 0, NULL, 0},
        };

        int opt;

        while ((opt = getopt_long(argc, argv, ShortOptions, LongOptions, NULL)) != -1) {

            switch(opt) {
                case 'h':
                case ARG_HELP:
                    bShowHelp = true;
                    return 0;

                case 'v':
                case ARG_VERBOSE:
                    Verbose++;
                    break;

                case 'i':
                case ARG_DATA_PATH:
                    DataPaths.emplace_back(optarg);
                    break;

                case 'L':
                    LibraryPaths.emplace_back(optarg);
                    break;

                case 'o':
                case ARG_OUTPUT_PREFIX:
                    OutputPrefix.assign(optarg);
                    break;

                case ARG_DEBUG:
                    bDebug = true;
                    break;

                case '?':
                default:
                    return -1;
            }
        }

        assert(!bShowHelp);

        for (int i = optind; i < argc; i++) {
            SourceFiles.emplace_back(argv[i]);
        }

        return 0;
    }

    void Reset() {
        DataPaths.clear();
        LibraryPaths.clear();
        SourceFiles.clear();

        OutputPrefix.clear();

        bDebug = false;
        Verbose = 0;
        bShowHelp = false;
    };


    void Usage()
    {
        std::cerr << "Usage" << std::endl;
    }

    void Dump()
    {
        std::cerr << "bDebug: " << bDebug << std::endl;
        std::cerr << "Verbose: " << Verbose << std::endl;
        std::cerr << "bShowHelp: " << bShowHelp << std::endl;
        std::cerr << "OutputPrefix: " << OutputPrefix << std::endl;

        std::cerr << "DataPaths: " << std::endl;
        for(const auto& p : DataPaths) {
            std::cerr << "\t" << p << std::endl;
        }

        std::cerr << "LibraryPaths: " << std::endl;
        for(const auto& p : LibraryPaths) {
            std::cerr << "\t" << p << std::endl;
        }

        std::cerr << "SourceFiles: " << std::endl;
        for(const auto& p : SourceFiles) {
            std::cerr << "\t" << p << std::endl;
        }

    }

public:
    std::vector<std::string> DataPaths;
    std::vector<std::string> LibraryPaths;
    std::vector<std::string> SourceFiles;

    std::string OutputPrefix;

    bool bDebug;
    int Verbose;
    bool bShowHelp;
};

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

    FOption Option;
    if (Option.Parse(argc, argv)) {
        Option.Usage();
        return -1;
    }
    if (Option.bShowHelp) {
        Option.Usage();
        return 0;
    }

    if (Option.Verbose) {
        Option.Dump();
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
