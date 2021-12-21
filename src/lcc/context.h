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

    FMolecule() {}
    virtual ~FMolecule() {}

    FMolecule(const std::string& InName) : Name(InName), Id(InName) {}

    const std::string GetName() const {
        return Name;
    }
    const std::string GetId() const {
        return Id;
    }

    void Print(std::ostream& os) {
        os << "  Molecule Id: " << Id << std::endl;
    }

};

class FPathway {
public:
    std::string Name;
    std::vector<std::string> Sequence;

    FPathway(const std::string& InName, std::vector<std::string>& InSequence)
        : Name(InName), Sequence(InSequence) {};
};

class FSmallMolecule : public FMolecule {
public:
    FSmallMolecule() {}

    FSmallMolecule(const std::string& InName) : FMolecule(InName) {}

    void Print(std::ostream& os) {
        os << "  SmallMolecule Id: " << Name << std::endl;
    }
};
 
class FProtein : public FMolecule {
public:
    FProtein() {}

    FProtein(const std::string& InName) : FMolecule(InName) {}

};

class FEnzyme : public FProtein {
public:
    std::string Substrate;
    float kcat;
    float kM;

    FEnzyme() {}

    FEnzyme(const std::string& InName, const std::string& InSubstrate, const float& Inkcat, const float& InkM)
        : Substrate(InSubstrate), kcat(Inkcat), kM(InkM), FProtein(InName) {}

    void Print(std::ostream& os) {
        os << "  Enzyme Id: " << Name << " | Substrate: " << Substrate << "\tkcat:  " << std::to_string(kcat) << "\tkM: " << std::to_string(kM) << std::endl;
    }
};

class FComplex : public FMolecule {
public:
    FComplex() {}

    FComplex(const std::string& InName) : FMolecule(InName) {}
};

class FPolymerase : public FMolecule{
public:
    std::string Substrate;
    float Rate;

    FPolymerase() {}

    FPolymerase(const std::string& InName, const std::string& InSubstrate, const float& InRate) 
        : Substrate(InSubstrate), Rate(InRate), FMolecule(InName) {}

    void Print(std::ostream& os) {
        os << "  Polymerase Id: " << Name << "\tSubstrate: " << Substrate << "\tRate:  " << std::to_string(Rate) << std::endl;
    }
};

class FReaction {
public:
    std::string Name; // Name is equivalent to EnzymeName for now
    std::map<std::string, int> Stoichiometry;

    virtual ~FReaction() {}

    FReaction(const std::string& InName, const std::map<std::string, int>& InStoichiometry) 
        : Name(InName), Stoichiometry(InStoichiometry) {}

    bool CheckIfReactant(const std::string& Name) {
        return (Stoichiometry[Name] < 0);
    }
    bool CheckIfProduct(const std::string& Name) {
        return (Stoichiometry[Name] > 0);
    }

    void Print(std::ostream& os) {
        os << "  Reaction Id: " << Name << " | ";
        for (auto& Stoich : Stoichiometry) {
            os << "[" << Stoich.first << ", " << Stoich.second << "], ";
        }
        os << std::endl;
    }
};

class FEnzymaticReaction : public FReaction {
public:
    std::string Enzyme;

    FEnzymaticReaction(const std::string& InName, const std::map<std::string, int>& InStoichiometry, const std::string& InEnzyme)
        : Enzyme(InEnzyme), FReaction(InName, InStoichiometry) {}
};

class FPolymeraseReaction : public FReaction {
public:
    std::string Polymerase;
    std::vector<std::string> BuildingBlocks;
    // std::map<std::string, int> Stoichiometry;

    FPolymeraseReaction(const std::string& InName, const std::map<std::string, int>& InStoichiometry, const std::string& InPolymerase, const std::vector<std::string>& InBuildingBlocks)
        : Polymerase(InPolymerase), BuildingBlocks(InBuildingBlocks), FReaction(InName, InStoichiometry) {}
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
    FTable PolymeraseTable;
    FTable PathwayTable;

    std::vector<std::string> UsingModuleList;
    std::vector<FMolecule*> MoleculeList;
    std::vector<FReaction*> ReactionList;
    std::vector<FPathway> PathwayList;
    std::vector<std::string> IdentifierList;

    void PrintLists(std::ostream& os);
    void SaveUsingModuleList(const char *Filename);
    void AddToMoleculeList(FMolecule *NewMolecule);
    void AddToReactionList(FReaction *NewReaction);

    const std::string QueryTable(const std::string& Name, const std::string& Property, FTable Table);

    std::vector<std::string> GetNames_MoleculeList();

    std::vector<std::string> GetNames_EnzymeList(std::vector<const FEnzyme *> EnzymeList);
    std::vector<std::string> GetSubstrateNames_EnzymeList(std::vector<const FEnzyme *> EnzymeList);
    std::vector<float> Getkcats_EnzymeList(std::vector<const FEnzyme *> EnzymeList);
    std::vector<float> GetkMs_EnzymeList(std::vector<const FEnzyme *> EnzymeList);

    std::vector<std::string> GetNames_ReactionList();
    std::vector<std::string> GetSubstrateNames_ReactionList();
    std::vector<std::string> GetReactantNames_ReactionList();
    std::vector<std::string> GetProductNames_ReactionList();

    std::vector<std::string> GetNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);
    std::vector<std::string> GetSubstrateNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);
    std::vector<std::string> GetReactantNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);
    std::vector<std::string> GetProductNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);
    std::vector<std::string> GetEnzymeNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);

    std::vector<std::string> GetNames_PathwayList();
    std::vector<std::string> GetSequences_PathwayList();

    std::vector<std::vector<int>> GetStoichiometryMatrix_EnzymaticReaction(std::vector<const FEnzymaticReaction *> EnzymaticReactionList);

    int GetIdxByName_MoleculeList(std::string InputName);

    std::vector<const FEnzyme *> GetList_Enzyme_MoleculeList();
    std::vector<const FSmallMolecule *> GetList_SmallMolecule_MoleculeList();
    std::vector<const FEnzymaticReaction *> GetList_Enzymatic_ReactionList();
    
    std::vector<int> GetIdx_Enzyme_MoleculeList();
    std::vector<int> GetIdx_EnzymeSubstrate_MoleculeList();
    std::vector<int> GetIdx_SmallMolecule_MoleculeList();
    

    std::vector<int> GetIdxListFromList(std::vector<std::string> InputList, std::vector<std::string> RefList);
//    std::vector<int> GetEnzSubstrateIdxFromAllSubstrates();

    std::vector<float> GetFreqMatrixForChromosomes();
    std::vector<float> GetFreqMatrixForRNAs();
    std::vector<float> GetFreqMatrixForProteins();
};

#endif /* LCC_CONTEXT_H */
