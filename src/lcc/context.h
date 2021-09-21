#ifndef LCC_CONTEXT_H
#define LCC_CONTEXT_H

#include <string>
#include <vector>
#include <map>

#include "option.h"
#include "node.h"

typedef std::map<std::string, std::string> FTableRecord;

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
    FProtein(const NProteinDeclaration& InProteinDecl) : ProteinDecl(InProteinDecl) {};

    const std::string& GetName() const {
        return ProteinDecl.Id.Name;
    }

private:
    const NProteinDeclaration ProteinDecl;
};

class FReaction {
public:

};

class FModule {
public:
};

class FTable {
public:
	std::vector<FTableRecord> Records;

	void LoadFromTSV(const char *Filename);
	void Dump();
	void Dump(const std::vector<std::string>& InKeys);

};


class FCompilerContext {
public:

    void Init(const FOption& InOption);

    FTable GeneTable;
    FTable ReactionTable;


    std::vector<std::string> UsingModuleList;
    std::vector<FProtein> ProteinList;
    std::vector<std::string> PathwayList;



    void SaveUsingModuleList(const char *Filename);
};

#endif /* LCC_CONTEXT_H */
