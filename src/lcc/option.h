#ifndef LCC_OPTION_H
#define LCC_OPTION_H

#include <vector>
#include <string>

class FOption {
public:
    FOption() {
        Reset();
        PreferredSetting();

    }

    int Parse(int argc, char *argv[]);
    void Reset();
    void Usage(const char *argv0);
    void ShowVersion(const char *argv0);
    void Dump();
    void PreferredSetting(); 


public:
    std::vector<std::string> DataPaths;
    std::vector<std::string> LibraryPaths;
    std::vector<std::string> SourceFiles;

    std::string OutputPrefix;

    std::string SimResultFile = "SimOut.tsv";
    std::string SimModuleFile = "SimModule.py";

    bool bSimCpp;
    bool bSimPython = true;

    bool bDebug;
    bool bVersion;
    int Verbose;
    bool bShowHelp;
    bool bParseOnly;
    bool bRunInOmVisim;
};

#endif //LCC_OPTION_H
