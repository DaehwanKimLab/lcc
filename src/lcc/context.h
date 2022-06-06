#ifndef LCC_CONTEXT_H
#define LCC_CONTEXT_H

#include <string>
#include <vector>
#include <map>

#include "option.h"
#include "node.h"
#include "util.h"
#include "bioinfo.h"

typedef std::map<std::string, std::string> FTableRecord;

class FCount {
public:
    std::string Name;
    float Amount;
    float Begin;
    float End;
    float Step;

    bool bMolarity;

    FCount() {}

    FCount(std::string InName, float InAmount) : Name(InName), Amount(InAmount), Begin(0), End(0), Step(0), bMolarity(false) {}

    FCount(std::string InName, float InAmount, std::vector<float> InRange, bool InbMolarity) : Name(InName), Amount(InAmount), bMolarity(InbMolarity) {
        Utils::Assertion((InRange.size() == 3), "Range requires three float values");
        Begin = InRange[0];
        End   = InRange[1];
        Step  = InRange[2];
    }

    void Print(std::ostream& os) {
        os << "[Count] "
              "Name : " << Name <<
           "\t| Amount: " << Utils::SciFloat2Str(Amount) <<
           "\t| Begin: " << Utils::SciFloat2Str(Begin) <<
           "\t| End: " << Utils::SciFloat2Str(End) <<
           "\t| Step: " << Utils::SciFloat2Str(Step);
        if (bMolarity) { os << "\t| (Molarity)"; }
        os << std::endl;
    }
};

class FContainerSpace {
public:
    std::vector<std::string> Names;

    FContainerSpace() {}

    void Print(std::ostream& os) {
        os << "Compartments: ";
        if (!Names.empty()) {
            for (auto& name: Names) {
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
    std::vector<float> Coord;
    float Angle; // initial angle not connected to simulation yet
    FCount * Count;

    FLocation() {}

    FLocation(std::string InName, std::vector<float> InCoord) : Name(InName), Coord(InCoord), Angle(0) {}

    FLocation(std::string InName, std::vector<float> InCoord, FCount* InCount) : Name(InName), Coord(InCoord), Angle(0), Count(InCount) {}

    void Print(std::ostream& os) {
        os << "[Location] ";

        os << Name << " : ";

        os << "(";
        for (auto& coord : Coord) {
            os << Utils::SciFloat2Str(coord) << ", ";
        } os << "), ";

        os << "(Angle: " << Utils::SciFloat2Str(Angle) << ")";
        os << std::endl;

        os << "(Count: " << Count->Amount << ")";
        os << std::endl;
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

    FComposition(std::string InName) : Name(InName) {}

    FComposition(std::string InName, std::string InSurface, std::string InCore) : Name(InName) {
        std::pair <std::string, float> S(InSurface, 100);
        Constituent_Surface.push_back(S);
        std::pair <std::string, float> C(InCore, 100 );
        Constituent_Core.push_back(C);
    }

    FComposition(std::string InName, std::vector<std::pair<std::string, float>> InSurface, std::vector<std::pair<std::string, float>> InCore)
        : Name(InName), Constituent_Surface(InSurface), Constituent_Core(InCore) {}
};

class FContainer {
public:
    std::string Name;

    FContainer() {}
    virtual ~FContainer() {}

    FContainer(std::string InName) : Name(InName) {}

};

class FCompartment : public FContainer { // spatial container
public:
    std::string Shape; // rod(bacilli), sphere(cocci), spiral, filamentous, box, comma(vibrio), amoeba
    std::vector<float> Dimension; // width, height, depth: if single provided, assume isotropy, like a diameter

    float pH = -1; // -1 as default
    float Temp = -1; // -1 as default

    FCompartment() {}

    FCompartment(std::string InName)
        : Shape("Sphere"), FContainer(InName) {} // sphere as default

    FCompartment(std::string InName, std::string InShape)
        : Shape(InShape), FContainer(InName) {}

    FCompartment(std::string InName, std::string InShape, float InDimension)
        : Shape(InShape), Dimension{InDimension}, FContainer(InName) {}

    FCompartment(std::string InName, std::string InShape, std::vector<float> InDimension)
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

    FBiologicalCompartment(std::string InName) : FCompartment(InName) {}

    FBiologicalCompartment(std::string InName, std::string InType) : Type(InType), FCompartment(InName) {}

};

class FOrganism : public FBiologicalCompartment {
public:
    std::string Species; // special type set aside for organism
    std::string Strain;

    FOrganism() {}

    FOrganism(std::string InName) : FBiologicalCompartment(InName) {}

    FOrganism(std::string InName, std::string InSpecies) : Species(InSpecies), FBiologicalCompartment(InName) {}

    FOrganism(std::string InName, std::string InSpecies, std::string InStrain) : Species(InSpecies), Strain(InStrain), FBiologicalCompartment(InName) {}
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

    FProkaryote(std::string InName) : FOrganism(InName) {}

    FProkaryote(std::string InName, std::string InSpecies) : FOrganism(InName, InSpecies) {}

};

class FEukaryote : public FOrganism {
public:
    bool Multicellularity;

};

class FMolecule {
public:
    std::string Name;
    std::string Id;
    
    float Threshold = -1;

    FMolecule() {}
    virtual ~FMolecule() {}

    FMolecule(std::string InName) : Name(InName), Id(InName) {}

    FMolecule(std::string InName, std::string InId) : Name(InName), Id(InId) {}

    std::string GetName() {
        return Name;
    }
    std::string GetId() {
        return Id;
    }
    void SetThreshold(float InThreshold) {
        Threshold = InThreshold;
    }

    void Print(std::ostream& os) {
        os << "[Molecule] Id: " << Id << std::endl;
    }

};

class FPathway {
public:
    std::string Name;
    std::vector<std::string> Sequence;

    FPathway(std::string InName)
        : Name(InName) {}
    FPathway(std::string InName, std::vector<std::string>& InSequence)
        : Name(InName), Sequence(InSequence) {}

    void Print(std::ostream& os) {
        os << "[Pathway] Id: " << Name;
        if (!Sequence.empty()) {
            os << ", Sequence: [" << Utils::JoinStr2Str(Sequence) << "]";
        }
        os << std::endl;
    }
};

class FSmallMolecule : public FMolecule {
public:
    FSmallMolecule() {}

    FSmallMolecule(std::string InName) : FMolecule(InName) {}

    void Print(std::ostream& os) {
        os << "[Small Molecule] Id: " << Id;
        os << std::endl;
    }
};

class FPolymer : public FMolecule {
public:
    std::vector<std::pair<std::string, float>> Composition; // Name of monomers and their proportion

    FPolymer() {}

    FPolymer(std::string InName) : FMolecule(InName) {}

    FPolymer(std::string InName, std::vector<std::pair<std::string, float>> InComposition)
        : Composition(InComposition), FMolecule(InName) {}

    FPolymer(std::string InName, std::vector<std::string> InCompositionNameOnly)
        : FMolecule(InName) {
        for (auto &CompositionName: InCompositionNameOnly) {
            Composition.push_back(std::pair<std::string, float> (CompositionName, 1 / InCompositionNameOnly.size()));
        }
    }

    void Print_Composition(std::ostream& os) {
        for (auto& composition : Composition){
            os << composition.first << " : " << composition.second << ", ";
        }
    }
};

class FPolymer_Static : public FPolymer {
public:
    int Size;

    FPolymer_Static() {}

    FPolymer_Static(std::string InName)
        : FPolymer(InName) {}

    FPolymer_Static(std::string InName, std::vector<std::string> InCompositionNameOnly)
            : FPolymer(InName, InCompositionNameOnly) {}

    FPolymer_Static(std::string InName, int InSize, std::vector<std::pair<std::string, float>> InComposition)
        : Size(InSize), FPolymer(InName, InComposition) {}

    FPolymer_Static(std::string InName, int InSize, std::vector<std::string> InCompositionNameOnly)
        : Size(InSize), FPolymer(InName, InCompositionNameOnly) {}

    void SetSize(int InSize) {
        Size = InSize;
    }
};

class FPolymer_TemplateBased : public FPolymer_Static {
public:

    FMolecule* Template;   // Added in context organization

    FPolymer_TemplateBased() {}

    FPolymer_TemplateBased(std::string InName)
        : FPolymer_Static(InName) {}

    FPolymer_TemplateBased(std::string InName, std::vector<std::string> InCompositionNameOnly)
        : FPolymer_Static(InName, InCompositionNameOnly) {}

    FPolymer_TemplateBased(std::string InName, int InSize, std::vector<std::pair<std::string, float>> InComposition)
        : FPolymer_Static(InName, InSize, InComposition) {}

    FPolymer_TemplateBased(std::string InName, int InSize, std::vector<std::string> InCompositionNameOnly)
        : FPolymer_Static(InName, InSize, InCompositionNameOnly) {}

    void SetTemplate(FMolecule* Molecule) {
        Template = Molecule;
    }
};

class FGeneticMaterial : public FPolymer_TemplateBased {
public:
    // std::string Symbol;
    std::string Sequence;
    std::vector<std::string> BuildingBlocks;

    // Either enter (size + composition) || sequence
    // Sequence can be used to calculate size and composition)

    FGeneticMaterial() {}

    FGeneticMaterial(std::string InName)
            : FPolymer_TemplateBased(InName) {}

    FGeneticMaterial(std::string InName, std::vector<std::string> InBuildingBlocks)
            : BuildingBlocks(InBuildingBlocks), FPolymer_TemplateBased(InName, InBuildingBlocks) {} // Convenient default example for genes would be FGene(A, 1000, NT);

    FGeneticMaterial(std::string InName, int InSize, std::vector<std::pair<std::string, float>> InComposition)
        : FPolymer_TemplateBased(InName, InSize, InComposition) {}

    FGeneticMaterial(std::string InName, int InSize, std::vector<std::string> InBuildingBlocks)
        : BuildingBlocks(InBuildingBlocks), FPolymer_TemplateBased(InName, InSize, InBuildingBlocks) {} // Convenient default example for genes would be FGene(A, 1000, NT);

    FGeneticMaterial(std::string InName, std::string InSequence, std::vector<std::string> InBuildingBlocks)
        : Sequence(InSequence), BuildingBlocks(InBuildingBlocks), FPolymer_TemplateBased(InName) { GetCompositionFromSequence(InSequence, InBuildingBlocks); }

    void GetCompositionFromSequence(std::string InSequence, std::vector<std::string> InBuildingBlocks) {
        Size = InSequence.size();
        for (auto& buildingblock : InBuildingBlocks) {
            auto Char = BioInfo::GetBuildingBlockAbbr(buildingblock);
            int Count = std::count(InSequence.begin(), InSequence.end(), Char);
            std::string Char_Str(1, Char);
            std::pair<std::string, int> CharNCount(Char_Str, Count);
            Composition.push_back(CharNCount);
        }
    }
    
    void UpdateInfoWithSequence(std::string InSequence) {
        Sequence = InSequence;
        GetCompositionFromSequence(Sequence, BuildingBlocks);
    }

};

class FChromosome : public FGeneticMaterial {
public:
    FChromosome() {}

    FChromosome(std::string InName)
        : FGeneticMaterial(InName, BioInfo::GetBuildingBlocks("dNT")) {}

    FChromosome(std::string InName, int InSize)
        : FGeneticMaterial(InName, InSize, BioInfo::GetBuildingBlocks("dNT")) {}

    FChromosome(std::string InName, std::string InSequence)
        : FGeneticMaterial(InName, InSequence, BioInfo::GetBuildingBlocks("dNT")) {}

    void Print(std::ostream& os) {
        os << "[Chromosome] Id: " << Name << " | ";
        Print_Composition(os);
        os << std::endl;
    }
};

class FGene : public FGeneticMaterial {   // to be revisted for its categorization. Make a new FEATURE class?
public:
    std::vector<std::string> Promoters;

    FGene() {}

    FGene(std::string InName)
        : FGeneticMaterial(InName, 1000, BioInfo::GetBuildingBlocks("dNT")) {}

    FGene(std::string InName, int InSize)
        : FGeneticMaterial(InName, InSize, BioInfo::GetBuildingBlocks("dNT")) {} // Convenient default example for genes would be FGene(A_g, 1000, NT);

    FGene(std::string InName, std::string InSequence)
        : FGeneticMaterial(InName, InSequence, BioInfo::GetBuildingBlocks("dNT")) {}

    void Print(std::ostream& os) {
        os << "[Gene] Id: " << Name << " | ";
        Print_Composition(os);
        os << std::endl;
    }
};

class FRNA : public FGeneticMaterial {
public:

    FRNA() {}

    FRNA(std::string InName)
    : FGeneticMaterial(InName, 1000, BioInfo::GetBuildingBlocks("NT")) {}

    FRNA(std::string InName, int InSize)
    : FGeneticMaterial(InName, InSize, BioInfo::GetBuildingBlocks("NT")) {} // Convenient default example for genes would be FRNA(A_r, 1000, NT);

    FRNA(std::string InName, std::string InSequence)
    : FGeneticMaterial(InName, InSequence, BioInfo::GetBuildingBlocks("NT")) {}

    void Print(std::ostream& os) {
        os << "[RNA] Id: " << Name << " | ";
        Print_Composition(os);
        os << std::endl;
    }
};

class FmRNA : public FRNA {
public:

    FmRNA() {}

    FmRNA(std::string InName)
    : FRNA(InName) {}

    FmRNA(std::string InName, int InSize)
    : FRNA(InName, InSize) {}

    FmRNA(std::string InName, std::string InSequence)
    : FRNA(InName, InSequence) {}

    void Print(std::ostream& os) {
        os << "[mRNA] Id: " << Name << " | ";
        Print_Composition(os);
        os << std::endl;
    }
};

class FrRNA : public FRNA {
public:

    FrRNA() {}

    FrRNA(std::string InName)
    : FRNA(InName) {}

    FrRNA(std::string InName, int InSize)
    : FRNA(InName, InSize) {}

    FrRNA(std::string InName, std::string InSequence)
    : FRNA(InName, InSequence) {}

    void Print(std::ostream& os) {
        os << "[rRNA] Id: " << Name << " | ";
        Print_Composition(os);
        os << std::endl;
    }
};

class FtRNA : public FRNA {
public:

    FtRNA() {}

    FtRNA(std::string InName)
    : FRNA(InName) {}

    FtRNA(std::string InName, int InSize)
    : FRNA(InName, InSize) {}

    FtRNA(std::string InName, std::string InSequence)
    : FRNA(InName, InSequence) {}

    void Print(std::ostream& os) {
        os << "[tRNA] Id: " << Name << " | ";
        Print_Composition(os);
        os << std::endl;
    }
};

class FmiscRNA : public FRNA {
public:

    FmiscRNA() {}

    FmiscRNA(std::string InName)
    : FRNA(InName) {}

    FmiscRNA(std::string InName, int InSize)
    : FRNA(InName, InSize) {}

    FmiscRNA(std::string InName, std::string InSequence)
    : FRNA(InName, InSequence) {}

    void Print(std::ostream& os) {
        os << "[miscRNA] Id: " << Name << " | ";
        Print_Composition(os);
        os << std::endl;
    }
};

class FProtein : public FGeneticMaterial {
public:

    FProtein() {}

    FProtein(std::string InName)
    : FGeneticMaterial(InName, 333, BioInfo::GetBuildingBlocks("AA")) {}

    FProtein(std::string InName, int InSize)
    : FGeneticMaterial(InName, InSize, BioInfo::GetBuildingBlocks("AA")) {} // Convenient default example for genes would be FProtein(A_p, 333, AA);

    FProtein(std::string InName, std::string InSequence)
    : FGeneticMaterial(InName, InSequence, BioInfo::GetBuildingBlocks("AA")) {}

    void Print(std::ostream& os) {
        os << "[Protein] Id: " << Name << " | ";
        Print_Composition(os);
        os << std::endl;
    }
};

class FEnzyme : public FProtein {
public:
    // Kinetics is a vector of Substrate:(kcat, KM...) data. Maybe extended by adding temperature, pH, etc.
    std::vector<std::pair<std::string, std::vector<float>>> Kinetics;

    FEnzyme() {}

    // UPDATE TO FProtein(InName) when ID System is set up
    FEnzyme(std::string InName, std::vector<std::pair<std::string, std::vector<float>>> InKinetics)
        : Kinetics(InKinetics), FProtein(InName) {}

    FEnzyme(std::string InName)
        : FProtein(InName) {}

    void Print(std::ostream& os) {
        os << "[Enzyme] Id: " << Name;
        if (!Kinetics.empty()) {
            for (auto& kinetics : Kinetics) {
                os << "\t| Substrate:  " << kinetics.first;
                int i = 0;
                for (auto& Float: kinetics.second) {
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

class FPolymerase : public FProtein { // Note: specific to template-dependent polymerase class.
public:
    std::string InitiationSite;
    std::string TerminationSite;
    std::string TemplateClass;
    std::string TargetClass;
    std::string Process;
    float Rate;

    FPolymerase() {}

    FPolymerase(std::string InName, std::string InInitiationSite, float& InRate, std::string InTerminationSite)
        : InitiationSite(InInitiationSite), Rate(InRate), TerminationSite(InTerminationSite), FProtein(InName) {}

    FPolymerase(std::string InName, std::string InInitiationSite, float& InRate, std::string InTerminationSite, std::string InTemplateClass, std::string InTargetClass, std::string InProcess)
        : InitiationSite(InInitiationSite), Rate(InRate), TerminationSite(InTerminationSite), TemplateClass(InTemplateClass), TargetClass(InTargetClass), Process(InProcess), FProtein(InName) {}
//
//    FPolymerase(std::string InName, std::string InType, std::string InTemplate, std::string InTarget, std::string InProcess, float& InRate)
//        : Template(InTemplate), Type(InType), Target(InTarget), Process(InProcess), Rate(InRate), FProtein(InName) {}

    void Print(std::ostream& os) {
        os << "[Polymerase] Id: " << Name << "\tTemplate: " << Template << "\tRate:  " << Utils::SciFloat2Str(Rate) << std::endl;
    }
};

class FDNAP : public FPolymerase {
public:
    FDNAP() {}

    FDNAP(std::string InName, std::string InInitiationSite, float& InRate, std::string InTerminationSite)
        : FPolymerase(InName, InInitiationSite, InRate, InTerminationSite, "Chromosome", "Chromosome", "Replication") {}

    void Print(std::ostream& os) {
        os << "[DNAP] Id: " << Name << "\tTemplate: " << Template << "\tRate:  " << Utils::SciFloat2Str(Rate) << std::endl;
        os << "\tTemplateClass: " << TemplateClass << "\tTargetClass: " << TargetClass << std::endl;
    }
};

class FRNAP : public FPolymerase {
public:
    FRNAP() {}

    FRNAP(std::string InName, std::string InInitiationSite, float& InRate, std::string InTerminationSite)
        : FPolymerase(InName, InInitiationSite, InRate, InTerminationSite, "Gene", "RNA", "Transcription") {}

    void Print(std::ostream& os) {
        os << "[RNAP] Id: " << Name << "\tTemplate: " << Template << "\tRate:  " << Utils::SciFloat2Str(Rate) << std::endl;
        os << "\tTemplateClass: " << TemplateClass << "\tTargetClass: " << TargetClass << std::endl;
    }
};

class FRibosome : public FPolymerase {
public:
    FRibosome() {}

    FRibosome(std::string InName, std::string InInitiationSite, float& InRate, std::string InTerminationSite)
        : FPolymerase(InName, InInitiationSite, InRate, InTerminationSite, "mRNA", "Protein", "Translation") {}

    void Print(std::ostream& os) {
        os << "[Ribosome] Id: " << Name << "\tTemplate: " << Template << "\tRate:  " << Utils::SciFloat2Str(Rate) << std::endl;
        os << "\tTemplateClass: " << TemplateClass << "\tTargetClass: " << TargetClass << std::endl;
    }
};

class FTransporter : public FProtein {
public:
    float ki;
    float ko;

    FTransporter(std::string InName, float Inki, float Inko)
        : ki(Inki), ko(Inko), FProtein(InName) {}

    void Print(std::ostream& os) {
        os << "[Transporter] Id: " << Name << "\t| ki:  " << Utils::SciFloat2Str(ki) << "\t| ko:  " << Utils::SciFloat2Str(ko) << std::endl;
    }
};

class FComplex : public FPolymer_Static {
public:
    FComplex() {}

    FComplex(std::string InName) : FPolymer_Static(InName) {}
};

class FReaction {
public:
    std::string Name;
    std::vector<std::pair<std::string, int>> Stoichiometry;

    int Type = -1; // default

    // temporary pathway namespace
    std::string Pathway;

    virtual ~FReaction() {}

    FReaction(std::string InName, std::vector<std::pair<std::string, int>>& InStoichiometry)
        : Name(InName), Stoichiometry(InStoichiometry) {}

    bool CheckIfReactant(std::string Query) {
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
    bool CheckIfProduct(std::string Query) {
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

    virtual void Print(std::ostream& os) {
        os << "[Reaction]";
        Print_IdStoichiometry(os);
        os << std::endl;
    }

    void Print_IdStoichiometry(std::ostream& os) {
        os << "  Reaction Id: " << Name << " | ";
        for (auto& Stoich : Stoichiometry) {
            os << "[" << Stoich.first << ", " << Stoich.second << "], ";
        }
    }

    void AddPathway(std::string InPathway) {
        Pathway = InPathway;
        std::cout << Name << " in " << Pathway << std::endl;
    }

};

class FStandardReaction : public FReaction {
public:
    float k;
    float krev;

    FStandardReaction(std::string InName, std::vector<std::pair<std::string, int>>& InStoichiometry, float Ink, float Inkrev)
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

    FRegulatoryReaction(std::string InName, std::vector<std::pair<std::string, int>>& InStoichiometry, float InK, float Inn, std::string InEffect, std::string InMode)
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

    FEnzymaticReaction(std::string InName, std::vector<std::pair<std::string, int>>& InStoichiometry, std::string InEnzyme)
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

    FEnz_StandardReaction(std::string InName, std::vector<std::pair<std::string, int>>& InStoichiometry, std::string InEnzyme, float Ink, float Inkrev)
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

    FPolymeraseReaction(std::string InName, std::vector<std::pair<std::string, int>>& InStoichiometry, std::string InPolymerase, std::vector<std::string> InBuildingBlocks)
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

class FTransporterReaction : public FReaction {
public:
    float D;
    std::string Transporter;

    FTransporterReaction(std::string InName, std::vector<std::pair<std::string, int>>& InStoichiometry, std::string InTransporter)
        : Transporter(InTransporter), FReaction(InName, InStoichiometry) {}

    FTransporterReaction(std::string InName, std::vector<std::pair<std::string, int>>& InStoichiometry, float InD)
        : D(InD), FReaction(InName, InStoichiometry) {}

    FTransporterReaction(std::string InName, std::vector<std::pair<std::string, int>>& InStoichiometry, float InD, std::string InTransporter)
        : Transporter(InTransporter), D(InD), FReaction(InName, InStoichiometry) {}

    void Print(std::ostream& os) {
        os << "[Transporter Reaction]" ;
        Print_IdStoichiometry(os);
        os << "  D: " << D << ", ";
    }
};

class FMotility {
public:
    std::string Name;

    FMotility() {}
    virtual ~FMotility() {}

    FMotility(std::string InName) : Name(InName) {}

    virtual void Print(std::ostream& os) {
        os << "[Motility]" ;
        os << "  Name: " << Name;
    }
};

class FSwim : public FMotility {
public:
    std::vector<std::pair<std::string, float>> Thresholds;

    FSwim() {}

    FSwim(std::string InName)
            : FMotility(InName) {}

    FSwim(std::string InName, std::vector<std::pair<std::string, float>> InThresholds)
            : Thresholds(InThresholds), FMotility(InName) {}

    void Print(std::ostream& os) override {
        os << "[Swim]" ;
        os << " Name: " << Name;
        if (!Thresholds.empty()) {
            os << ", Thresholds: ";
            for (int i = 0; i < Thresholds.size(); i++) {
                os << "(Thresholded Molecule: " << Thresholds[i].first << ", Threshold Value: " << Thresholds[i].second << "), ";
            }
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

    void Init(FOption& InOption);

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
//    std::vector<FGene*>         GeneList;
//    std::vector<FRNA*>	        RNAList;
//    std::vector<FProtein*>      ProteinList;
//    std::vector<FEnzyme*>       EnzymeList;
    std::vector<FReaction*>     ReactionList;
    std::vector<FPathway*>      PathwayList;
    std::vector<std::string>    IdentifierList;
    std::vector<FContainer*>    ContainerList;
//    std::vector<FOrganism*>     OrganismList;
    std::vector<FCount*>        CountList;
    std::vector<FLocation*>     LocationList;
    std::vector<FComposition*>  CompositionList;
    std::vector<FMotility*>     MotilityList;
            std::vector<std::pair<std::string, float>>     ThresholdList; // temporary, list of mol name and threshold value pairs


    void PrintLists(std::ostream& os);
    void PrintInitialCounts(std::ostream& os);
    void PrintInitialLocations(std::ostream& os);
    void SaveUsingModuleList(const char *Filename);

    // Base class list addition
    void AddToMoleculeList(FMolecule *NewMolecule);
    void AddToReactionList(FReaction *NewReaction);
    void AddToPathwayList(FPathway *NewPathway);
    void AddToContainerList(FContainer *NewContainer);
    void AddToCountList(FCount *NewCount);
        void AddToCountList(std::string Name, float Count);
    void AddToLocationList(FLocation *NewLocation);
    void AddToMotilityList(FMotility *NewMotility);
//    void AddToCompositionList(FComposition *NewComposition);

    // Compiler organization
    void Organize();
    void ApplyDefaultGeneticInformationProcessingOnMoleculeList();
    void AssignReactionType(FReaction *Reaction, int Regulation);
    void AssignReactionTypesForReactionList();
    void MakeListsFromMoleculeList();
    void MakeListsFromContainerList();
    void MakePathwayLists();
                void MakeListsFromMotilityList();

    bool CheckForEcoli();
    void DefaultSetUp_Ecoli();

    // Tables
    std::string QueryTable(std::string Name, std::string Property, FTable Table);

    // EnzymeList
//    FEnzyme * GetEnzyme_EnzymeList(std::string Name);
//    float GetFloatAttributeByName_EnzymeList(std::vector<FEnzyme *> EnzymeList, std::string EnzymeName, std::string Attribute);
//    std::string GetStringAttributeByName_EnzymeList(std::vector<FEnzyme *> EnzymeList, std::string EnzymeName, std::string Attribute);

    // PolymeraseList - to be updated
//    std::vector<std::string> GetNames_PolymeraseList(std::vector<FPolymerase *> PolymeraseList);
//    std::vector<float> GetRates_PolymeraseList(std::vector<FPolymerase *> PolymeraseList);

    std::vector<std::string> GetNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *>);
    std::vector<std::string> GetSubstrateNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *>);
    std::vector<std::string> GetReactantNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *>);
    std::vector<std::string> GetProductNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *>);
    std::vector<std::string> GetBuildingBlockNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *>);

    std::vector<std::string> GetNames_PathwayList();
    std::vector<std::string> GetSequences_PathwayList();

    // Compiler Utility
    FGeneticMaterial * GenerateChromosome(std::string MolName, int Count, int Size); // default
    FGeneticMaterial * GenerateChromosome(std::string MolName, int Count, std::string Sequence); // default

    FGeneticMaterial * GenerateCounterpart_Gene(std::string MolName, std::string Sequence, int Count); // default
    FGeneticMaterial * GenerateCounterpart_RNA(std::string MolName, int Count, std::string RNAType); // default
    FGeneticMaterial * GenerateCounterpart_Protein(std::string MolName, int Count); // default


    // Stoichiometry Matrix-related
    std::vector<int> AddUniqueSubstrateIdxToIdxList(FReaction *, std::vector<int>);
    std::vector<int> GetUniqueSubstrateIdx_ReactionList(std::vector<FReaction *> ReactionList);
    std::vector<int> GetCoefficientArray(FReaction *, std::vector<int>);
    std::vector<std::vector<int>> GetStoichiometryMatrix(std::vector<FReaction *> ReactionList);
//    std::vector<std::vector<int>> GetStoichiometryMatrix_PolymeraseReaction(std::vector<FPolymeraseReaction *>);

    // TODO: Update all XXList functions to take the XXList itself to enable more flexible list manipulation

    // LocationList
    std::vector<float> GetLocationByName_LocationList(std::string MolName);
    std::vector<FLocation *> GetSubList_LocationList(std::string Type);
    std::vector<std::string> GetNames_LocationList(std::string Type);
    std::vector<std::string> GetUniqueNames_LocationList(std::string Type);
    std::vector<FContainer*> FCompilerContext::GetUniqueContainers_LocationList(std::string Type);
    int GetCounts_LocationList(std::string Type);

    // CountList
    std::vector<std::string> GetNames_CountList(std::string Type);
    std::vector<FCount *> GetSubList_CountList(std::string Type);
    float GetInitialCountByName_CountList(std::string InputName);
    bool GetMolarityFactorByName_CountList(std::string InputName);
    bool CheckMolarityFactorTrueForAny_CountList();
    void RevertMolarity_CountList(std::string InName); 
    void AdjustMolarity_PseudoMolecule();

    // MoleculeList
    FMolecule * GetMolecule_MoleculeList(std::string Name);
    std::vector<FMolecule *> GetSubList_MoleculeList(std::string Type);
    std::vector<std::string> GetNameListByType_MoleculeList(std::string Type);
    std::vector<std::string> GetNameList_MoleculeList(std::vector<FMolecule *> ListOfMolecules);
    int GetIdxByName_MoleculeList(std::string InputName);
    int GetIdx_MoleculeList(FMolecule * Molecule);
    std::vector<int> GetIdxByStrList_MoleculeList(std::vector<std::string>);
    std::vector<int> GetIdxListByType_MoleculeList(std::string Type);
    std::vector<int> GetIdxList_MoleculeList(std::vector<FMolecule *> ListOfMolecules);
    std::vector<int> GetLocalIdxList_MoleculeList(std::vector<FMolecule *> SubListOfMolecules, std::vector<FMolecule *> ListOfMolecules);

//    std::vector<int> GetIdx_PolymeraseReactionSubstrate_ByPolymeraseName_MoleculeList(std::string);
    bool CheckDuplicates_List(std::string ListType, std::string Name);

    // ContainerList
    std::vector<FContainer *> GetSubList_ContainerList(std::string Type);
    std::vector<std::string> GetNames_ContainerList(std::string Type);
    int GetCounts_ContainerList(std::string Type);

    // ReactionList
    std::vector<std::string> GetNames_ReactionList(std::string Type);
    std::vector<FReaction *> GetSubList_ReactionList(std::string Type, std::string NameSpace_Pathway=""); // useful for stoichiometry
    std::vector<FStandardReaction *> GetList_Standard_ReactionList(std::string Type);
    std::vector<FRegulatoryReaction *> GetList_Regulatory_ReactionList(std::string Type);
    std::vector<std::string> GetNames_EnzymaticReactionList(std::vector<FEnzymaticReaction *> EnzymaticReactionList);

    FReaction * GetReactionByPolymeraseName_ReactionList(std::string InPolymeraseName);


    // TODO: Update polymerase reaction writing system
    std::vector<FPolymeraseReaction *> GetList_Polymerase_ReactionList();

//    std::vector<int> GetEnzSubstrateIdxFromAllSubstrates();

    std::vector<float> GetFreqMatrixForChromosomes();
    std::vector<float> GetFreqMatrixForRNAs();
    std::vector<float> GetFreqMatrixForProteins();

    // temporary namespace-related functions
    std::vector<FReaction *> FilterByPathway_ReactionList(std::vector<FReaction *> ListOfReactions, std::string PathwayName);
    std::string GetAssociatedPathway_ReactionList(FMolecule* Molecule);
    std::vector<FMolecule *> GetSubListByReactionList_MoleculeList(std::vector<FReaction *> ListOfReactions);
};


#endif /* LCC_CONTEXT_H */
