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
    std::string Name;

    FProtein() {}

    FProtein(const std::string& InName) : Name(InName) {}

    const std::string GetName() const {
        return Name;
    }
    
};

class FEnzyme : public FProtein {
public:
    std::string Substrate;
    float kcat;
    float kM;

    FEnzyme() {}

    FEnzyme(const std::string& InName, const std::string& InSubstrate, const float& Inkcat, const float& InkM)
        : Substrate(InSubstrate), kcat(Inkcat), kM(InkM), FProtein(InName) {}

};

class FReaction {
public:
    std::string Name; // Name is equivalent to EnzymeName for now
    std::map<std::string, int> Stoichiometry;

    FReaction(const std::string& InName, const std::map<std::string, int>& InStoichiometry) 
        : Name(InName), Stoichiometry(InStoichiometry) {}

    bool CheckIfReactant(const std::string& Name) {
        return (Stoichiometry[Name] < 0);
    }
    bool CheckIfProduct(const std::string& Name) {
        return (Stoichiometry[Name] > 0);
    }
};

class FEnzymaticReaction : public FReaction {
public:
    std::string Enzyme;

    FEnzymaticReaction(const std::string& InName, const std::map<std::string, int>& InStoichiometry, const std::string& InEnzyme)
        : Enzyme(InEnzyme), FReaction(InName, InStoichiometry) {}
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
    std::vector<FEnzymaticReaction> EnzymaticReactionList;
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

//    void PrintList(std::vector List) {
//        for (auto& item : List){
//            std::cout << item.Name << ", ";
//        };
//        std::cout << std::endl;
//    }

    void PrintLists(std::ostream& os) {
        os << "Compiler Context Lists:" << std::endl;
        
        os << "  PathwayList: " << std::endl << "  " << "  ";
        for (auto& item : PathwayList){
            os << item.Name << ", ";
        };
        os << std::endl;

        os << "  EnzymeList: " << std::endl << "  " << "  ";
        for (auto& item : EnzymeList){
            os << item.Name << ", ";
        };
        os << std::endl;

        os << "  EnzymaticReactionList: " << std::endl << "  " << "  ";
        for (auto& item : EnzymaticReactionList){
            os << item.Name << ", ";
        };
        os << std::endl;      
    }

    void SaveUsingModuleList(const char *Filename);
};

#endif /* LCC_CONTEXT_H */
