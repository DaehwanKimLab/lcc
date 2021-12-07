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
    std::string Name;
    std::vector<std::string> Sequence;

    FPathway(const std::string& InName, std::vector<std::string>& InSequence)
        : Name(InName), Sequence(InSequence) {};
};

class FProtein {
public:
//    std::string Name;

    FProtein(const NProteinDeclaration& InProteinDecl) : ProteinDecl(InProteinDecl) {};
//    FProtein(std::string& InName) : Name(InName) {};

    const std::string& GetName() const {
        return ProteinDecl.Id.Name;
    }

private:
    const NProteinDeclaration ProteinDecl;
};

class FEnzyme { //: public FProtein {
public:
    std::string Name; // remove after resolving inheritance issue
    std::string Substrate;
    float kcat;
    float kM;

    FEnzyme(const std::string& InName, const std::string& InSubstrate, const float& Inkcat, const float& InkM)
        : Name(InName), Substrate(InSubstrate), kcat(Inkcat), kM(InkM) {}; // , FProtein(InName) {};

};

class FReaction {
public:
    std::string ReactionName; // currently equivalent to EnzymeName

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
    FTable EnzymeTable;


    std::vector<std::string> UsingModuleList;
    std::vector<FProtein> ProteinList;
    std::vector<FEnzyme> EnzymeList;
    std::vector<FPathway> PathwayList;
    std::vector<std::string> IdentifierList;

    const std::string QueryEnzymeTable(const std::string& Name, const std::string& Property){
        for (auto& record : EnzymeTable.Records) {
            if (Name == record["EnzymeName"]) {
                return record[Property];
            }          
        }
        std::cout << "No such enzyme found in the database: " << Name << std::endl;
    }

    const std::string QueryReactionTable(const std::string& Name, const std::string& Property){
        for (auto& record : ReactionTable.Records) {
            if (Name == record["EnzymeName"]) {
                return record[Property];
            }          
        }
        std::cout << "No such reaction found in the database: " << Name << std::endl;
    }

    void SaveUsingModuleList(const char *Filename);
};

#endif /* LCC_CONTEXT_H */
