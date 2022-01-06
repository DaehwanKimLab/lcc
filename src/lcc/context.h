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
    int InitialCount = -1;
    bool Fixed = false;

    FMolecule() {}
    virtual ~FMolecule() {}

    FMolecule(const std::string& InName) : Name(InName), Id(InName) {}

    FMolecule(const std::string& InName, const std::string& InId) : Name(InName), Id(InId) {}

    FMolecule(const std::string& InName, const int InInitialCount) : Name(InName), Id(InName), InitialCount(InInitialCount), Fixed(false) {}

    FMolecule(const std::string& InName, const int InInitialCount, const bool InFixed) : Name(InName), Id(InName), InitialCount(InInitialCount), Fixed(InFixed) {}


    const std::string GetName() const {
        return Name;
    }
    const std::string GetId() const {
        return Id;
    }

    void Print(std::ostream& os) {
        os << "  Molecule Id: " << Id << std::endl;
        if (InitialCount < -1) {
            os << "    InitialCount: " << InitialCount;
        }
        if (Fixed) {
            os << "(Fixed)" << std::endl;
        }
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

    FSmallMolecule(const std::string& InName, const int InInitialCount) : FMolecule(InName, InInitialCount) {}

    FSmallMolecule(const std::string& InName, const int InInitialCount, const bool InFixed) : FMolecule(InName, InInitialCount, InFixed) {}

    void Print(std::ostream& os) {
        os << "  Small Molecule Id: " << Id;
        if (InitialCount >= 0 ) {
            os << "\t| InitialCount: " << InitialCount;
        }
        if (Fixed) {
            os << " (Fixed)";
        }
        os << std::endl;
    }
};

class FPolymer : public FMolecule {
public:
    std::map<std::string, float> Composition;

    FPolymer() {}

    FPolymer(const std::string& InName) : FMolecule(InName) {}

    FPolymer(const std::string& InName, std::map<std::string, float> InComposition) 
        : Composition(InComposition), FMolecule(InName) {}
 
};

class FPolymer_Static : public FPolymer {
public:
    int Size;

    FPolymer_Static() {}

    FPolymer_Static(const std::string& InName) 
        : FPolymer(InName) {}

    FPolymer_Static(const std::string& InName, int InSize) 
        : Size(InSize), FPolymer(InName) {}

    FPolymer_Static(const std::string& InName, const std::map<std::string, float> InComposition) 
        : FPolymer(InName, InComposition) {}

    FPolymer_Static(const std::string& InName, const std::map<std::string, float> InComposition, int InSize) 
        : Size(InSize), FPolymer(InName, InComposition) {}
};

class FChromosome : public FPolymer_Static {   
public:
    std::string Symbol; // roman numeral? such as ChI, ChII?
    
    FChromosome() {}

    FChromosome(const std::string& InName) 
        : FPolymer_Static(InName) {}

    FChromosome(const std::string& InName, const std::string& InSymbol) 
        : Symbol(InSymbol), FPolymer_Static(InName) {}

    FChromosome(const std::string& InName, int InSize) 
        : FPolymer_Static(InName, InSize) {}

    FChromosome(const std::string& InName, const std::map<std::string, float> InComposition, int InSize, const std::string& InSymbol) 
        : Symbol(InSymbol), FPolymer_Static(InName, InComposition, InSize) {}

    void Print(std::ostream& os) {
        os << "  Chromosome Id: " << Name << std::endl;
    }
};

class FGene : public FPolymer_Static {   // to be revisted for its categorization. Make a new FEATURE class?
public:
    std::string Symbol;

    FGene() {}

    FGene(const std::string& InName) 
        : FPolymer_Static(InName) {}

    FGene(const std::string& InName, const std::string& InSymbol) 
        : Symbol(InSymbol), FPolymer_Static(InName) {}

    FGene(const std::string& InName, const std::map<std::string, float> InComposition, int InSize, const std::string& InSymbol) 
        : Symbol(InSymbol), FPolymer_Static(InName, InComposition, InSize) {}

    void Print(std::ostream& os) {
        os << "  Gene Id: " << Name << std::endl;
    }
};

class FRNA : public FPolymer_Static {
public:
    std::string Type;

    FRNA() {}

    FRNA(const std::string& InName) 
        : FPolymer_Static(InName) {}

    FRNA(const std::string& InName, const std::string& InType) 
        : Type(InType), FPolymer_Static(InName) {}

    FRNA(const std::string& InName, const std::map<std::string, float> InComposition, int InSize, const std::string& InType) 
        : Type(InType), FPolymer_Static(InName, InComposition, InSize) {}

    void Print(std::ostream& os) {
        os << "  RNA Id: " << Name << std::endl;
    }
};
 
class FProtein : public FPolymer_Static {
public:
    std::string Symbol;

    FProtein() {}  

    FProtein(const std::string& InName) 
        : FPolymer_Static(InName) {}

    FProtein(const std::string& InName, const std::string& InSymbol) 
        : Symbol(InSymbol), FPolymer_Static(InName) {}

    FProtein(const std::string& InName, const std::map<std::string, float> InComposition, int InSize, const std::string& InSymbol) 
        : Symbol(InSymbol), FPolymer_Static(InName, InComposition, InSize) {}

    void Print(std::ostream& os) {
        os << "  Protein Id: " << Name << std::endl;
    }
};

class FEnzyme : public FMolecule { // update to FProtein when ID system is set up
public:
    std::string Substrate;
    float kcat = -1;
    float kM = -1;

    std::string Inhibitor;
    float ki = -1;

    // standard chemical reaction rates for approximation
    float k = -1;
    float krev = -1;

    FEnzyme() {}

    // UPDATE TO FProtein(InName) when ID System is set up
// without InitialCount
    FEnzyme(const std::string& InName, const std::string& InSubstrate, const float& Inkcat, const float& InkM)
        : Substrate(InSubstrate), kcat(Inkcat), kM(InkM), FMolecule(InName) {}

    FEnzyme(const std::string& InName, const std::string& InSubstrate, const float& Inkcat, const float& InkM, const std::string& InInhibitor, const float& Inki)
        : Substrate(InSubstrate), kcat(Inkcat), kM(InkM),  Inhibitor(InInhibitor), ki(Inki), FMolecule(InName) {}

    // standard chemical reaction case may have Product specified
    FEnzyme(const std::string& InName, const float& Ink, const float& Inkrev)
        : k(Ink), krev(Inkrev), FMolecule(InName) {}

// with InitialCount
    FEnzyme(const std::string& InName, const std::string& InSubstrate, const float& Inkcat, const float& InkM, const int& InInitialCount, const bool Fixed)
        : Substrate(InSubstrate), kcat(Inkcat), kM(InkM), FMolecule(InName, InInitialCount, Fixed) {}

    FEnzyme(const std::string& InName, const std::string& InSubstrate, const float& Inkcat, const float& InkM, const std::string& InInhibitor, const float& Inki, const int& InInitialCount, const bool Fixed)
        : Substrate(InSubstrate), kcat(Inkcat), kM(InkM),  Inhibitor(InInhibitor), ki(Inki), FMolecule(InName, InInitialCount, Fixed) {}

    // standard chemical reaction case may have Product specified
    FEnzyme(const std::string& InName, const float& Ink, const float& Inkrev, const int& InInitialCount, const bool Fixed)
        : k(Ink), krev(Inkrev), FMolecule(InName, InInitialCount, Fixed) {}

    void Print(std::ostream& os) {
        os << "  Enzyme Id: " << Name; 
        if (kM > 0) {
            os << "\t| Substrate: " << Substrate << "\t| kcat:  " << std::to_string(kcat) << "\t| kM: " << std::to_string(kM);
        }
        if (!Inhibitor.empty()) {
            os << "\t| Inhibitor: " << Inhibitor << "\t| ki:  " << std::to_string(ki);
        }        
        if (k > 0) {
            os << "\t| k: " << std::to_string(k) << "\t| krev:  " << std::to_string(krev);
        }
        if (InitialCount > 0) {
            os << "\t| InitialCount: " << InitialCount;
        }
        if (Fixed) {
            os << " (Fixed)";
        }
        os << std::endl;
    }
};

class FPolymerase : public FMolecule{ // update to FProtein when ID system is set up with SQL // Note: specific to template-dependent polymerase class.
public:
    std::string Template;
    std::string Target;
    std::string Process; // temporary
    float Rate;

    FPolymerase() {}

    FPolymerase(const std::string& InName, const std::string& InTemplate, const std::string& InTarget, const float& InRate) 
        : Template(InTemplate), Target(InTarget), Rate(InRate), FMolecule(InName) {}

    FPolymerase(const std::string& InName, const std::string& InTemplate, const std::string& InTarget, const std::string& InProcess, const float& InRate) 
//        : Template(InTemplate), Target(InTarget), Process(InProcess), Rate(InRate), FProtein(InName) {}
        : Template(InTemplate), Target(InTarget), Process(InProcess), Rate(InRate), FMolecule(InName) {}

    void Print(std::ostream& os) {
        os << "  Polymerase Id: " << Name << "\tTemplate: " << Template << "\tTarget: " << Target << "\tRate:  " << std::to_string(Rate) << std::endl;
    }
};

class FComplex : public FMolecule {
public:
    FComplex() {}

    FComplex(const std::string& InName) : FMolecule(InName) {}
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
    FTable RNATable;
    FTable ProteinTable;
	    FTable ReactionTable;
    FTable EnzymeTable;
    FTable PolymeraseTable;
    FTable PathwayTable;
	    FTable InitialCountTable_TCA;

    std::vector<std::string> 	UsingModuleList;
    std::vector<FMolecule*> 	MoleculeList;
    std::vector<FGene*> 	GeneList;
    std::vector<FRNA*>		RNAList;
    std::vector<FProtein*>	ProteinList;
    std::vector<FReaction*> 	ReactionList;
    std::vector<FPathway> 	PathwayList;
    std::vector<std::string> 	IdentifierList;

    void PrintLists(std::ostream& os);
    void SaveUsingModuleList(const char *Filename);
    void AddToMoleculeList(FMolecule *NewMolecule);
    void AddToReactionList(FReaction *NewReaction);

    const std::string QueryTable(const std::string& Name, const std::string& Property, FTable Table);

    std::vector<std::string> GetNames_MoleculeList();

    // TODO: merge the below methods to return to have vector<string> and vector<float> functions, taking key strings. 
    std::vector<std::string> GetNames_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<std::string> GetSubstrateNames_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<float> Getkcats_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<float> GetkMs_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<std::string> GetInhibitorNames_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<float> Getkis_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<float> Getks_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<float> Getkrevs_EnzymeList(std::vector<const FEnzyme *>);

    std::vector<std::string> GetNames_PolymeraseList(std::vector<const FPolymerase *> PolymeraseList);
    std::vector<float> GetRates_PolymeraseList(std::vector<const FPolymerase *> PolymeraseList);

    std::vector<std::string> GetNames_ReactionList();
    std::vector<std::string> GetSubstrateNames_ReactionList();
    std::vector<std::string> GetReactantNames_ReactionList();
    std::vector<std::string> GetProductNames_ReactionList();

    // update by removing input data
    std::vector<std::string> GetNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);
    std::vector<std::string> GetSubstrateNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);
    std::vector<std::string> GetReactantNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);
    std::vector<std::string> GetProductNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);
    std::vector<std::string> GetEnzymeNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *>);

    std::vector<std::string> GetNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *>);
    std::vector<std::string> GetSubstrateNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *>);
    std::vector<std::string> GetReactantNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *>);
    std::vector<std::string> GetProductNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *>);
    std::vector<std::string> GetBuildingBlockNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *>);

    std::vector<std::string> GetNames_PathwayList();
    std::vector<std::string> GetSequences_PathwayList();

    // Stoichiometry Matrix-related
    std::vector<int> GetIdxForStoichiometryMatrix(std::string);
    std::vector<std::vector<int>> GetStoichiometryMatrix(std::string);
    std::vector<std::vector<int>> GetStoichiometryMatrix_PolymeraseReaction(std::vector<const FPolymeraseReaction *>);



    int GetIdxByName_MoleculeList(std::string InputName);
    int GetInitialCountByName_MoleculeList(std::string InputName);
    std::vector<int> GetIdxByStrList_MoleculeList(std::vector<std::string>);
    std::vector<int> GetInitialCountByStrList_MoleculeList(std::vector<std::string>);

    std::vector<const FGene *> GetList_Gene_MoleculeList();
    std::vector<const FRNA *> GetList_RNA_MoleculeList();
    std::vector<const FProtein *> GetList_Protein_MoleculeList();
    std::vector<const FEnzyme *> GetList_Enzyme_MoleculeList();
    std::vector<const FPolymerase *> GetList_Polymerase_MoleculeList();
    std::vector<const FSmallMolecule *> GetList_SmallMolecule_MoleculeList();
    std::vector<const FEnzymaticReaction *> GetList_Enzymatic_ReactionList();
    std::vector<const FPolymeraseReaction *> GetList_Polymerase_ReactionList();
    
    std::vector<int> GetIdxListFromMoleculeList(std::string FClassName);
    std::vector<std::string> GetNameListFromMoleculeList(std::string FClassName);

    std::vector<int> GetIdx_EnzymeSubstrate_MoleculeList();
//    std::vector<int> GetIdx_PolymeraseSubstrate_MoleculeList();

    std::vector<int> GetIdx_PolymeraseReactionSubstrate_ByPolymeraseName_MoleculeList(std::string);

    std::vector<int> GetIdxOfStrListFromStrList(std::vector<std::string> InputList, std::vector<std::string> RefList);
//    std::vector<int> GetEnzSubstrateIdxFromAllSubstrates();

    std::vector<float> GetFreqMatrixForChromosomes();
    std::vector<float> GetFreqMatrixForRNAs();
    std::vector<float> GetFreqMatrixForProteins();
};

#endif /* LCC_CONTEXT_H */
