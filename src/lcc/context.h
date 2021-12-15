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
    FTable ProteinTable;
    FTable EnzymeTable;
    FTable PathwayTable;

    std::vector<std::string> UsingModuleList;
    std::vector<FProtein> ProteinList;
    std::vector<FEnzyme> EnzymeList;
    std::vector<FPathway> PathwayList;
    std::vector<FEnzymaticReaction> EnzymaticReactionList;
    std::vector<std::string> IdentifierList;

    void PrintLists(std::ostream& os);
    void SaveUsingModuleList(const char *Filename);

    const std::string QueryEnzymeTable(const std::string& Name, const std::string& Property);
    const std::string QueryReactionTable(const std::string& Name, const std::string& Property);

    std::vector<std::string> GetNames_EnzymeList();
    std::vector<std::string> GetSubstrates_EnzymeList();
    std::vector<float> Getkcats_EnzymeList();
    std::vector<float> GetkMs_EnzymeList();

    std::vector<std::string> GetNames_ReactionList();
    std::vector<std::string> GetSubstrateNames_ReactionList();
    std::vector<std::string> GetReactantNames_ReactionList();
    std::vector<std::string> GetProductNames_ReactionList();

    std::vector<std::string> GetNames_EnzymaticReactionList();
    std::vector<std::string> GetSubstrateNames_EnzymaticReactionList();
    std::vector<std::string> GetReactantNames_EnzymaticReactionList();
    std::vector<std::string> GetProductNames_EnzymaticReactionList();
    std::vector<std::string> GetEnzymeNames_EnzymaticReactionList();

    std::vector<std::string> GetNames_PathwayList();
    std::vector<std::string> GetSequences_PathwayList();

    std::vector<std::vector<int>> GetStoichiometryMatrix();

    std::vector<int> GetIdxListFromList(std::vector<std::string> InputList, std::vector<std::string> RefList);
    std::vector<int> GetEnzSubstrateIdxFromAllSubstrates();

    std::vector<float> GetFreqMatrixForChromosomes();
    std::vector<float> GetFreqMatrixForRNAs();
    std::vector<float> GetFreqMatrixForProteins();
};

#endif /* LCC_CONTEXT_H */
