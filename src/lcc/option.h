#ifndef LCC_OPTION_H
#define LCC_OPTION_H

#include <vector>
#include <string>

class FOption {
public:
    FOption() {
        Reset();
    }

    int Parse(int argc, char *argv[]);
    void Reset();
    void Usage(const char *argv0);
    void ShowVersion(const char *argv0);
    void Dump();


public:
    std::vector<std::string> DataPaths;
    std::vector<std::string> LibraryPaths;
    std::vector<std::string> SourceFiles;

    std::string OutputPrefix;

    bool bDebug;
    bool bVersion;
    int Verbose;
    bool bShowHelp;
};

#endif //LCC_OPTION_H
