#ifndef LCC_CONTEXT_H
#define LCC_CONTEXT_H

#include <string>
#include <vector>
#include <map>

#include "option.h"
#include "node.h"
#include "util.h"

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

class FCount {
public:
    // new 
    std::string Name;
    float Amount;   
    float Begin;
    float End;
    float Step;

    FCount() {}

    // new 
    FCount(std::string InName, float InAmount, float InBegin, float InEnd, float InStep) : Name(InName), Amount(InAmount), Begin(InBegin), End(InEnd), Step(InStep) {}

    void Print(std::ostream& os) {
        os << "Name : " << Name << 
         "\t| Amount: " << Utils::SciFloat2Str(Amount) << 
         "\t| Begin: " << Utils::SciFloat2Str(Begin) << 
         "\t| End: " << Utils::SciFloat2Str(End) << 
         "\t| Step: " << Utils::SciFloat2Str(Step) << std::endl;
    }
};

class FMolecule {
public:
    std::string Name;
    std::string Id;

    FMolecule() {}
    virtual ~FMolecule() {}

    FMolecule(const std::string& InName) : Name(InName), Id(InName) {}

    FMolecule(const std::string& InName, const std::string& InId) : Name(InName), Id(InId) {}

    const std::string GetName() const {
        return Name;
    }
    const std::string GetId() const {
        return Id;
    }

    void Print(std::ostream& os) {
        os << "  Molecule Id: " << Id << std::endl;
        os << std::endl;
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
        os << "  Small Molecule Id: " << Id;
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
    int Type; // 0 : Enz_Standard
              // 1 : Enz_MichaelisMenten

    std::string Substrate;
    float kcat = -1;
    float KM = -1;

    // standard chemical reaction rates
    float k = -1;
    float krev = -1;

    std::string Mode;

    FEnzyme() {}

    // UPDATE TO FProtein(InName) when ID System is set up
    FEnzyme(int InType, const std::string& InName, const std::string& InSubstrate, const float& Ink1, const float& Ink2)
        : Type(InType), Substrate(InSubstrate), FMolecule(InName) {
            if (Type == 0) {
                k = Ink1;	krev = Ink2;
            } else if (Type == 1) {
                kcat = Ink1; 	KM = Ink2;           
            }       
        }

    void Print(std::ostream& os) {
        os << "  Enzyme Id: " << Name << "\t| Substrate: " << Substrate; 
        if (Type == 0) {
            os << "\t| k: " << Utils::SciFloat2Str(k) << "\t| krev:  " << Utils::SciFloat2Str(krev);
        }
        if (Type == 1) {
            os << "\t| kcat:  " << Utils::SciFloat2Str(kcat) << "\t| KM: " << Utils::SciFloat2Str(KM);
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
        os << "  Polymerase Id: " << Name << "\tTemplate: " << Template << "\tTarget: " << Target << "\tRate:  " << Utils::SciFloat2Str(Rate) << std::endl;
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
    std::vector<std::pair<std::string, int>> Stoichiometry;

    int Type = -1; // default

    virtual ~FReaction() {}

    FReaction(const std::string& InName, const std::vector<std::pair<std::string, int>>& InStoichiometry) 
        : Name(InName), Stoichiometry(InStoichiometry) {}

    bool CheckIfReactant(const std::string& Query) {
        for (auto& stoich : Stoichiometry) {
           if (stoich.first == Query) {
               if (stoich.second < 0) {
                   return true;
               } else {
                   return false;
               }
           }
        }
        std::cout << "ERROR: Unable to find the molecule '" << Query << "'in the reaction of '" << Name << "'" << std::endl;
        return false;
    }
    bool CheckIfProduct(const std::string& Query) {
        for (auto& stoich : Stoichiometry) {
           if (stoich.first == Query) {
               if (stoich.second > 0) {
                   return true;
               } else {
                   return false;
               }
           }
        }
        std::cout << "ERROR: Unable to find the molecule '" << Query << "'in the reaction of '" << Name << "'" << std::endl;
        return false;
    }

    void Print_Stoichiometry(std::ostream& os) {
        os << "  Reaction Id: " << Name << " | ";
        for (auto& Stoich : Stoichiometry) {
            os << "[" << Stoich.first << ", " << Stoich.second << "], ";
        }
        os << std::endl;
    }
};

class FStandardReaction : public FReaction {
public:
    float k;
    float krev;

    FStandardReaction(const std::string& InName, const std::vector<std::pair<std::string, int>>& InStoichiometry, const float Ink, const float Inkrev)
        : k(Ink), krev(Inkrev), FReaction(InName, InStoichiometry) {}

    void Print(std::ostream& os) {
        os << "  [Standard Reaction]" << std::endl;
        Print_Stoichiometry(os);
        os << "  k: " << k << "\t| krev: " << krev << std::endl;
        os << std::endl;
    }
};

class FRegulatoryReaction : public FReaction {
public:
    float K;
    float n; // if Allosteric
    std::string Effect; // Activation or Inhibition
    std::string Mode; // Allosteric by default, competitive when n=-1 (temporary)

    FRegulatoryReaction(const std::string& InName, const std::vector<std::pair<std::string, int>>& InStoichiometry, const float InK, const float Inn, const std::string& InEffect, const std::string& InMode)
        : K(InK), n(Inn), Effect(InEffect), Mode(InMode), FReaction(InName, InStoichiometry) {}

    void Print(std::ostream& os) {
        os << "  [Regulatory Reaction]" << std::endl;
        Print_Stoichiometry(os);
        if 	(Effect == "Inhibition") 	{ os << "   Ki: "; } 
        else if (Effect == "Activation") 	{ os << "   Ka: "; }
        os << K;
        if 	(n > 0) 			{ os << "\t| n: " << n; }
        os << std::endl;
    }
};

class FEnzymaticReaction : public FReaction {
public:
    std::string Enzyme;

    FEnzymaticReaction(const std::string& InName, const std::vector<std::pair<std::string, int>>& InStoichiometry, const std::string& InEnzyme)
        : Enzyme(InEnzyme), FReaction(InName, InStoichiometry) {}
};

class FPolymeraseReaction : public FReaction {
public:
    std::string Polymerase;
    std::vector<std::string> BuildingBlocks; // may be moved to FPolymer
    // std::map<std::string, int> Stoichiometry;

    FPolymeraseReaction(const std::string& InName, const std::vector<std::pair<std::string, int>>& InStoichiometry, const std::string& InPolymerase, const std::vector<std::string>& InBuildingBlocks)
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
    std::vector<FEnzyme*>       EnzymeList;
    std::vector<FReaction*> 	ReactionList;
    std::vector<FPathway> 	PathwayList;
    std::vector<std::string> 	IdentifierList;
    std::vector<FCount*>        CountList;

    void PrintLists(std::ostream& os);
    void PrintInitialCounts(std::ostream& os);
    void SaveUsingModuleList(const char *Filename);
    void AddToMoleculeList(FMolecule *NewMolecule);
    void AddToReactionList(FReaction *NewReaction);

    // Compiler organization
    void Organize();
    void AssignReactionType(FReaction *Reaction, int Regulation);
    void AssignReactionTypesForReactionList();
    void MakeListsFromMoleculeList();

    const std::string QueryTable(const std::string& Name, const std::string& Property, FTable Table);

    std::vector<std::string> GetNames_MoleculeList();

    // TODO: merge the below methods to return to have vector<string> and vector<float> functions, taking key strings. 
    const FEnzyme * GetEnzyme_EnzymeList(std::string Name);
    std::vector<std::string> GetNames_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<std::string> GetSubstrateNames_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<float> Getkcats_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<float> GetKMs_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<float> Getks_EnzymeList(std::vector<const FEnzyme *>);
    std::vector<float> Getkrevs_EnzymeList(std::vector<const FEnzyme *>);

    // new merged methods for EnzymeList
    float GetFloatAttributeByName_EnzymeList(std::vector<const FEnzyme *> EnzymeList, std::string EnzymeName, std::string Attribute);
    std::string GetStringAttributeByName_EnzymeList(std::vector<const FEnzyme *> EnzymeList, std::string EnzymeName, std::string Attribute);

    std::vector<std::string> GetNames_PolymeraseList(std::vector<const FPolymerase *> PolymeraseList);
    std::vector<float> GetRates_PolymeraseList(std::vector<const FPolymerase *> PolymeraseList);

    std::vector<std::string> GetNames_ReactionList(std::string);
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
    std::vector<int> AddUniqueSubstrateIdxToIdxList(const FReaction *, std::vector<int>);
    std::vector<int> GetIdxForStoichiometryMatrix(std::string);
    std::vector<int> GetCoefficientArray(const FReaction *, std::vector<int>);
    std::vector<std::vector<int>> GetStoichiometryMatrix(std::string);
    std::vector<std::vector<int>> GetStoichiometryMatrix_PolymeraseReaction(std::vector<const FPolymeraseReaction *>);



    int GetIdxByName_MoleculeList(std::string InputName);
    std::vector<int> GetIdxByStrList_MoleculeList(std::vector<std::string>);
    // std::vector<int> GetInitialCountByStrList_MoleculeList(std::vector<std::string>);

    float GetInitialCountByName_CountList(std::string InputName);

    std::vector<const FGene *> GetList_Gene_MoleculeList();
    std::vector<const FRNA *> GetList_RNA_MoleculeList();
    std::vector<const FProtein *> GetList_Protein_MoleculeList();
    std::vector<const FEnzyme *> GetList_Enzyme_MoleculeList();
    std::vector<const FPolymerase *> GetList_Polymerase_MoleculeList();
    std::vector<const FSmallMolecule *> GetList_SmallMolecule_MoleculeList();

    std::vector<const FEnzymaticReaction *> GetList_Enzymatic_ReactionList();
    std::vector<const FPolymeraseReaction *> GetList_Polymerase_ReactionList();

    std::vector<const FReaction *> GetSubList_ReactionList(std::string); // useful for stoichiometry
    std::vector<const FStandardReaction *> GetList_Standard_ReactionList(std::string Type);
    std::vector<const FRegulatoryReaction *> GetList_Regulatory_ReactionList(std::string Type);
    
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
