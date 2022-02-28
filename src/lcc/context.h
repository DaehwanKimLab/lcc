#ifndef LCC_CONTEXT_H
#define LCC_CONTEXT_H

#include <string>
#include <vector>
#include <map>

#include "option.h"
#include "node.h"
#include "util.h"

typedef std::map<std::string, std::string> FTableRecord;

class FContainerSpace {
public:
    std::vector<std::string> Names;

    FContainerSpace() {}

    void Print(std::ostream& os) {
        os << "Compartments: ";
        if (!Names.empty()) {
            for (auto name: Names) {
                os << name << ", ";
            }
            os << std::endl;
        } else {
            os << "None" << std::endl;
        }
    }
};

class FExperiment {
public:
    std::string Name;
    std::string Description;

    FExperiment() {}
};

class FLocation {
public:
    std::string Name;
    std::vector<float> Coord_Abs;
    std::vector<float> Coord_Rel; // when all values < 1 
    float Angle; // to be expanded 

    FLocation() {}

    FLocation(const std::string InName, const std::vector<float> InCoord) : Name(InName) {
    // if all values > 1, then make it absolute
        bool Abs = false;
        for (auto& coord : InCoord) {
            if (coord > 1) {
                Abs = true;
            }
        }
        if (Abs) {
            for (auto& coord : InCoord) {
                Coord_Abs.push_back(coord);
            }
        } else     {
            for (auto& coord : InCoord) {
                Coord_Rel.push_back(coord);
            }
        }
    }

};

class FComposition {
public:
    std::string Name;
    std::vector<std::pair<std::string, float>> Constituent_Surface; // material name and percentage: if single name provided, assume 100%
    std::vector<std::pair<std::string, float>> Constituent_Core; // material name and percentage: if single name provided, assume 100%
    int Pattern_Surface = 0; // 0: well-mixed, 1: linear, 2: radial, etc.
    int Pattern_Core = 0;

    FComposition() {}

    FComposition(const std::string InName) : Name(InName) {}

    FComposition(const std::string InName, const std::string InSurface, const std::string InCore) : Name(InName) {
        std::pair <std::string, float> S(InSurface, 100); 
        Constituent_Surface.push_back(S);
        std::pair <std::string, float> C(InCore, 100 );
        Constituent_Core.push_back(C); 
    }

    FComposition(const std::string InName, const std::vector<std::pair<std::string, float>> InSurface, const std::vector<std::pair<std::string, float>> InCore) 
        : Name(InName), Constituent_Surface(InSurface), Constituent_Core(InCore) {}
};

class FContainer {
public:
    std::string Name;

    FContainer() {}
    virtual ~FContainer() {}

    FContainer(const std::string InName) : Name(InName) {}

};

class FCompartment : public FContainer { // spatial container
public:
    std::string Shape; // rod(bacilli), sphere(cocci), spiral, filamentous, box, comma(vibrio), amoeba
    std::vector<float> Dimension; // width, height, depth: if single provided, assume isotropy, like a diameter

    float pH = -1; // -1 as default
    float Temp = -1; // -1 as default
    
    FCompartment() {}

    FCompartment(const std::string InName) 
        : Shape("Sphere"), FContainer(InName) {} // sphere as default

    FCompartment(const std::string InName, const std::string InShape) 
        : Shape(InShape), FContainer(InName) {}

    FCompartment(const std::string InName, const std::string InShape, const float InDimension) 
        : Shape(InShape), Dimension{InDimension}, FContainer(InName) {}

    FCompartment(const std::string InName, const std::string InShape, const std::vector<float> InDimension) 
        : Shape(InShape), Dimension(InDimension), FContainer(InName) {}

};

class FNonBiologicalCompartment : public FCompartment {
public:
    std::string Description;
};

class FBiologicalCompartment : public FCompartment {
public:
    std::string Type;

    FBiologicalCompartment() {}

    FBiologicalCompartment(const std::string InName) : FCompartment(InName) {} 

    FBiologicalCompartment(const std::string InName, const std::string InType) : Type(InType), FCompartment(InName) {} 

};

class FOrganism : public FBiologicalCompartment {
public:
    std::string Species; // special type set aside for organism
    std::string Strain;
 
    FOrganism() {}

    FOrganism(const std::string InName) : FBiologicalCompartment(InName) {} 

    FOrganism(const std::string InName, const std::string InSpecies) : Species(InSpecies), FBiologicalCompartment(InName) {}

    FOrganism(const std::string InName, const std::string InSpecies, const std::string InStrain) : Species(InSpecies), Strain(InStrain), FBiologicalCompartment(InName) {}
};

class FOrganSystem : public FBiologicalCompartment {
public:
    std::string Type;
};

class FOrgan : public FBiologicalCompartment {
public:
    std::string Type;

};

class FTissue : public FBiologicalCompartment {
public:

};

class FCell : public FBiologicalCompartment {
public:

};

class FOrganelle : public FBiologicalCompartment {
public:

};


class FProkaryote : public FOrganism, public FCell { // Note multiple inheritance
public:
    bool Gram;

    FProkaryote() {}

    FProkaryote(const std::string InName) : FOrganism(InName) {} 

    FProkaryote(const std::string InName, const std::string InSpecies) : FOrganism(InName, InSpecies) {} 
    
};

class FEukaryote : public FOrganism {
public:
    bool Multicellularity;
    
};

class FCount {
public:
    std::string Name;
    float Amount;   
    float Begin;
    float End;
    float Step;

    bool bMolarity;

    FCount() {}

    FCount(std::string InName, float InAmount, float InBegin, float InEnd, float InStep, bool InbMolarity) : Name(InName), Amount(InAmount), Begin(InBegin), End(InEnd), Step(InStep), bMolarity(InbMolarity) {}

    void Print(std::ostream& os) {
        os << "Name : " << Name << 
         "\t| Amount: " << Utils::SciFloat2Str(Amount) << 
         "\t| Begin: " << Utils::SciFloat2Str(Begin) << 
         "\t| End: " << Utils::SciFloat2Str(End) <<
         "\t| Step: " << Utils::SciFloat2Str(Step);
         if (bMolarity) { os << "\t| (Molarity)"; }
         os << std::endl;
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
        os << "[Molecule] Id: " << Id << std::endl;
        os << std::endl;
    }

};

class FPathway {
public:
    std::string Name;
    std::vector<std::string> Sequence;

    FPathway(const std::string& InName, std::vector<std::string>& InSequence)
        : Name(InName), Sequence(InSequence) {}
};

class FSmallMolecule : public FMolecule {
public:
    FSmallMolecule() {}

    FSmallMolecule(const std::string& InName) : FMolecule(InName) {}

    void Print(std::ostream& os) {
        os << "[Small Molecule] Id: " << Id;
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
        os << "[Gene] Id: " << Name << std::endl;
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
        os << "[RNA] Id: " << Name << std::endl;
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
        os << "[Protein] Id: " << Name << std::endl;
    }
};

class FEnzyme : public FMolecule { // update to FProtein when ID system is set up
public:
    // Kinetics is a vector of Substrate:(kcat, KM...) data. Maybe extended by adding temperature, pH, etc.
    std::vector<std::pair<std::string, std::vector<float>>> Kinetics;

    FEnzyme() {}

    // UPDATE TO FProtein(InName) when ID System is set up
    FEnzyme(const std::string& InName, const std::vector<std::pair<std::string, std::vector<float>>> InKinetics)
        : Kinetics(InKinetics), FMolecule(InName) {}

    FEnzyme(const std::string& InName)
        : FMolecule(InName) {}

    void Print(std::ostream& os) {
        os << "[Enzyme] Id: " << Name; 
        if (!Kinetics.empty()) {
            for (auto kinetics : Kinetics) {
                os << "\t| Substrate:  " << kinetics.first;
                int i = 0;
                for (auto Float: kinetics.second) {
                    if (i == 0) {
                        os << ", kcat: " << Utils::SciFloat2Str(Float);
                    } else if (i == 1) {
                        os << ", KM: " << Utils::SciFloat2Str(Float);
                    }
                    i++;
                }
            }
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
        os << "[Polymerase] Id: " << Name << "\tTemplate: " << Template << "\tTarget: " << Target << "\tRate:  " << Utils::SciFloat2Str(Rate) << std::endl;
    }
};

class FComplex : public FMolecule {
public:
    FComplex() {}

    FComplex(const std::string& InName) : FMolecule(InName) {}
};

class FReaction {
public:
    std::string Name; 
    std::vector<std::pair<std::string, int>> Stoichiometry;

    int Type = -1; // default

    virtual ~FReaction() {}

    FReaction(const std::string& InName, const std::vector<std::pair<std::string, int>>& InStoichiometry) 
        : Name(InName), Stoichiometry(InStoichiometry) {}

    bool CheckIfReactant(const std::string& Query) const {
        for (auto& stoich : Stoichiometry) {
           if (stoich.first == Query) {
               if (stoich.second < 0) {
// std::cout << "Found the molecule '" << Query << "'in the reaction of '" << Name << "'" << std::endl;
                   return true;
               }
           }
        }
// std::cout << "Unable to find the molecule '" << Query << "'in the reaction of '" << Name << "'" << std::endl;
        return false;
    }
    bool CheckIfProduct(const std::string& Query) const {
        for (auto& stoich : Stoichiometry) {
           if (stoich.first == Query) {
               if (stoich.second > 0) {
// std::cout << "Found the molecule '" << Query << "'in the reaction of '" << Name << "'" << std::endl;
                   return true;
               }
           }
        }
// std::cout << "Unable to find the molecule '" << Query << "'in the reaction of '" << Name << "'" << std::endl;
        return false;
    }

    void Print_IdStoichiometry(std::ostream& os) {
        os << "  Reaction Id: " << Name << " | ";
        for (auto& Stoich : Stoichiometry) {
            os << "[" << Stoich.first << ", " << Stoich.second << "], ";
        }
    }
};

class FStandardReaction : public FReaction {
public:
    float k;
    float krev;

    FStandardReaction(const std::string& InName, const std::vector<std::pair<std::string, int>>& InStoichiometry, const float Ink, const float Inkrev)
        : k(Ink), krev(Inkrev), FReaction(InName, InStoichiometry) {}

    void Print(std::ostream& os) {
        os << "[Standard Reaction]" ;
        Print_IdStoichiometry(os);
        os << "  k: " << k << ", krev: " << krev << std::endl;
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
        os << "[Regulatory Reaction]" ;
        Print_IdStoichiometry(os);
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

    void Print(std::ostream& os) {
        os << "[Enzymatic Reaction]" ;
        Print_IdStoichiometry(os);
        os << std::endl;
    }
};

class FEnz_StandardReaction : public FEnzymaticReaction { // may have multiple inheritance (FStandardReaction) but it may complicate the logic for FStandardReaction
public:
    float k;
    float krev;

    FEnz_StandardReaction(const std::string& InName, const std::vector<std::pair<std::string, int>>& InStoichiometry, const std::string& InEnzyme, const float Ink, const float Inkrev)
        : k(Ink), krev(Inkrev), FEnzymaticReaction(InName, InStoichiometry, InEnzyme) {}

    void Print(std::ostream& os) {
        os << "[Enz_Standard Reaction]" ;
        Print_IdStoichiometry(os);
        os << "  k: " << k << ", krev: " << krev << std::endl;
    }

};

class FPolymeraseReaction : public FReaction {
public:
    std::string Polymerase;
    std::vector<std::string> BuildingBlocks; // may be moved to FPolymer
    // std::map<std::string, int> Stoichiometry;

    FPolymeraseReaction(const std::string& InName, const std::vector<std::pair<std::string, int>>& InStoichiometry, const std::string& InPolymerase, const std::vector<std::string>& InBuildingBlocks)
        : Polymerase(InPolymerase), BuildingBlocks(InBuildingBlocks), FReaction(InName, InStoichiometry) {}

    void Print(std::ostream& os) {
        os << "[Polymerase Reaction]" ;
        Print_IdStoichiometry(os);
        os << "  Polymerase: " << Polymerase << "\t| BuildingBlocks: ";
        for (auto& buildingblock : BuildingBlocks) {
            os << buildingblock << ", ";
        } 
        os << std::endl;
    }
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

    std::vector<std::string>    UsingModuleList;
    std::vector<FMolecule*>     MoleculeList;
    std::vector<FGene*>         GeneList;
    std::vector<FRNA*>	        RNAList;
    std::vector<FProtein*>      ProteinList;
    std::vector<FEnzyme*>       EnzymeList;
    std::vector<FReaction*>     ReactionList;
    std::vector<FPathway>       PathwayList;
    std::vector<std::string>    IdentifierList;
    std::vector<FContainer*>    ContainerList;
    std::vector<FOrganism*>     OrganismList;
    std::vector<FCount*>        CountList;
    std::vector<FLocation*>     LocationList;
    std::vector<FComposition*>  CompositionList;

    void PrintLists(std::ostream& os);
    void PrintInitialCounts(std::ostream& os);
    void SaveUsingModuleList(const char *Filename);

    void AddToMoleculeList(FMolecule *NewMolecule);
    void AddToReactionList(FReaction *NewReaction);
    void AddToContainerList(FContainer *NewContainer);
    void AddToCountList(FCount *NewCount);
    void AddToLocationList(FLocation *NewLocation);
//    void AddToCompositionList(FComposition *NewComposition);

    // Compiler organization
    void Organize();
    void AssignReactionType(FReaction *Reaction, int Regulation);
    void AssignReactionTypesForReactionList();
    void MakeListsFromMoleculeList();
    void MakeListsFromContainerList();

    // Tables
    const std::string QueryTable(const std::string& Name, const std::string& Property, FTable Table);

    // EnzymeList
    const FEnzyme * GetEnzyme_EnzymeList(std::string Name);
//    float GetFloatAttributeByName_EnzymeList(std::vector<const FEnzyme *> EnzymeList, std::string EnzymeName, std::string Attribute);
//    std::string GetStringAttributeByName_EnzymeList(std::vector<const FEnzyme *> EnzymeList, std::string EnzymeName, std::string Attribute);

    // PolymeraseList - to be updated
    std::vector<std::string> GetNames_PolymeraseList(std::vector<const FPolymerase *> PolymeraseList);
    std::vector<float> GetRates_PolymeraseList(std::vector<const FPolymerase *> PolymeraseList);

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

    // CountList
    float GetInitialCountByName_CountList(std::string InputName);
    bool GetMolarityFactorByName_CountList(std::string InputName);
    bool CheckMolarityFactorTrueForAny_CountList();
    void RevertMolarity_CountList(std::string InName); 
    void AdjustMolarity_PseudoMolecule();

    // MoleculeList
    std::vector<std::string> GetNames_MoleculeList();
    int GetIdxByName_MoleculeList(std::string InputName);
    std::vector<int> GetIdxByStrList_MoleculeList(std::vector<std::string>);
    // std::vector<int> GetInitialCountByStrList_MoleculeList(std::vector<std::string>);

    // useful?
    std::vector<int> GetIdxListFromMoleculeList(std::string FClassName);
    std::vector<std::string> GetNameListFromMoleculeList(std::string FClassName);
    std::vector<int> GetIdx_PolymeraseReactionSubstrate_ByPolymeraseName_MoleculeList(std::string);

    std::vector<const FGene *> GetList_Gene_MoleculeList();
    std::vector<const FRNA *> GetList_RNA_MoleculeList();
    std::vector<const FProtein *> GetList_Protein_MoleculeList();
    std::vector<const FEnzyme *> GetList_Enzyme_MoleculeList();
    std::vector<const FPolymerase *> GetList_Polymerase_MoleculeList();
    std::vector<const FSmallMolecule *> GetList_SmallMolecule_MoleculeList();



    // ReactionList
    std::vector<std::string> GetNames_ReactionList(std::string);
    std::vector<const FReaction *> GetSubList_ReactionList(std::string); // useful for stoichiometry
    std::vector<const FStandardReaction *> GetList_Standard_ReactionList(std::string Type);
    std::vector<const FRegulatoryReaction *> GetList_Regulatory_ReactionList(std::string Type);
    std::vector<std::string> GetNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *> EnzymaticReactionList);

    // useful?
    std::vector<const FEnzymaticReaction *> GetList_Enzymatic_ReactionList();
    std::vector<const FPolymeraseReaction *> GetList_Polymerase_ReactionList();
    

    std::vector<int> GetIdxOfStrListFromStrList(std::vector<std::string> InputList, std::vector<std::string> RefList);
//    std::vector<int> GetEnzSubstrateIdxFromAllSubstrates();

    std::vector<float> GetFreqMatrixForChromosomes();
    std::vector<float> GetFreqMatrixForRNAs();
    std::vector<float> GetFreqMatrixForProteins();
};


#endif /* LCC_CONTEXT_H */
