#include <iostream>
#include <cassert>
#include <getopt.h>

#include "option.h"

enum {
    ARG_START = 256,
    ARG_HELP,
    ARG_VERSION,
    ARG_VERBOSE,
    ARG_DEBUG,
    ARG_PARSEONLY,
    ARG_LIBRARY_PATH,
    ARG_DATA_PATH,
    ARG_OUTPUT_PREFIX,
	ARG_SIMOUT,
    ARG_SIMMODULE,
};

int FOption::Parse(int argc, char *argv[])
{
    static const char *ShortOptions = "hvL:i:o:";
    static const struct option LongOptions[] = {
            {"help", no_argument, NULL, 'h'},
            {"version", no_argument, NULL, ARG_VERSION},
            {"verbose", no_argument, NULL, 'v'},
            {"datadir", required_argument, NULL, 'i'},
            {"output", required_argument, NULL, 'o'},
            {"debug", no_argument, NULL, ARG_DEBUG},
            {"parse-only", no_argument, NULL, ARG_PARSEONLY},
		            {"simout", required_argument, NULL, ARG_SIMOUT},
            {"simmodule", required_argument, NULL, ARG_SIMMODULE},
            {NULL, 0, NULL, 0},
    };

    int opt;

    while ((opt = getopt_long(argc, argv, ShortOptions, LongOptions, NULL)) != -1) {

        switch(opt) {
            case 'h':
            case ARG_HELP:
                bShowHelp = true;
                return 0;

            case ARG_VERSION:
                bVersion = true;
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

            case ARG_PARSEONLY:
                bParseOnly = true;
                break;

			case ARG_SIMOUT:
				SimResultFile = std::string(optarg);
				break;
			
            case ARG_SIMMODULE:
                SimModuleFile = std::string(optarg);
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

void FOption::Reset() {
    DataPaths.clear();
    LibraryPaths.clear();
    SourceFiles.clear();

    OutputPrefix.clear();

    bDebug = false;
    bVersion = false;
    Verbose = 0;
    bShowHelp = false;
    bParseOnly = false;
}


static std::string GetFilename(const char *Fullname)
{
    std::string s(Fullname);
    return s.substr(s.find_last_of("/") + 1);
}

extern const char *VersionString;
void FOption::ShowVersion(const char *argv0)
{
    std::ostream& os = std::cout;
    std::string AppName = GetFilename(argv0);

    os << AppName << " (" << VersionString << ")" << std::endl;
}

void FOption::Usage(const char *argv0)
{
    std::ostream& os = std::cout;
    std::string AppName = GetFilename(argv0);

    os << "Usage: " << AppName << " [options] file..." << std::endl;
    os << "Options:" << std::endl;
    os << "  " << "-h, --help          Display this help" << std::endl;
    os << "  " << "--version           Display the version of the compiler" << std::endl;
    os << "  " << "-v, --verbose       Show more messages" << std::endl;
    os << "  " << "-L <dir>            Add <dir> to search path for libraries" << std::endl;
    os << "  " << "-i <dir>, --datadir=<dir>" << std::endl;
    os << "  " << "                    Add <dir> to search path for database" << std::endl;
    os << "  " << "-o <name>, --output=<name>" << std::endl;
    os << "  " << "                    Use <name> as output file prefix" << std::endl;
    os << "  " << "--debug             Enable debug mode" << std::endl;
    os << "  " << "--parse-only        Check syntax" << std::endl;
    os << "  " << "--simmodule           Simulation Module File" << std::endl;
}

void FOption::Dump()
{
    std::cerr << "bDebug: " << bDebug << std::endl;
    std::cerr << "Verbose: " << Verbose << std::endl;
    std::cerr << "bShowHelp: " << bShowHelp << std::endl;
    std::cerr << "bParseOnly: " << bParseOnly << std::endl;
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
