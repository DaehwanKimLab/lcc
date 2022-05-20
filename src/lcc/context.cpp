#include <iostream>
#include <fstream>
//#include <map>
#include <algorithm>
//#include <iterator>
#include <vector>

#include "context.h"
#include "option.h"
#include "number.h"
#include "util.h"

using namespace std;

extern FOption Option;

void FTable::LoadFromTSV(const char *Filename)
{
	if (Option.bDebug) {
		cerr << Filename << endl;
	}

	ifstream fp(Filename);
	string buf;

	vector<string> Headers;
	vector<string> Fields;

	// parsing header
	if (!getline(fp, buf)) {
		cerr << "Can't read header" << endl;
		return;
	}
    Utils::tokenize(buf, "\t", Headers);

	for(int i = 0; i < Headers.size(); i++) {
		Headers[i] = Utils::strip(Headers[i], "\"");

		if (Option.bDebug) {
			cerr << Headers[i] << endl;
		}
	}


	while(getline(fp, buf)) {
		Fields.clear();
        Utils::tokenize(buf, "\t", Fields);

		if (Fields.size() != Headers.size()) {
			cerr << "Wrong line: " << buf << endl;
			continue;
		}

		// Alloc new Record
		FTableRecord Record;

		for (int i = 0; i < Fields.size(); i++) {
			Record[Headers[i]] = Utils::strip(Fields[i], "\"");
		}
		Records.push_back(Record);
	}


	fp.close();
return;

}

void FTable::Dump()
{
	int index = 0;
	for(const auto& Record: Records) {
		cout << index++ << endl;
		for(const auto& it : Record) {
			cout << "  " << it.first << " ==> " << it.second << endl;
		}
	}
}

void FTable::Dump(const vector<string>& InKeys)
{
	int index = 0;
	for(auto& Record: Records) {
		cout << index++;
		for(const auto& key: InKeys) {
			cout << "\t" << "[" << key << ", " << Record[key] << "]";
		}
		cout << endl;
	}
}

void FCompilerContext::Init(FOption& InOption)
{
    if (InOption.DataPaths.size() > 0) {
        GeneTable.LoadFromTSV((InOption.DataPaths[0] + "/Database/genes.tsv").c_str());
        RNATable.LoadFromTSV((InOption.DataPaths[0] + "/Database/rnas.tsv").c_str());
        ProteinTable.LoadFromTSV((InOption.DataPaths[0] + "/Database/proteins.tsv").c_str());
		ReactionTable.LoadFromTSV((InOption.DataPaths[0] + "/Database/reactions.tsv").c_str());
        EnzymeTable.LoadFromTSV((InOption.DataPaths[0] + "/Database/EnzymeDatabase.txt").c_str());
        PolymeraseTable.LoadFromTSV((InOption.DataPaths[0] + "/Database/PolymeraseDatabase.txt").c_str());
        InitialCountTable_TCA.LoadFromTSV((InOption.DataPaths[0] + "/Database/InitialCount_TCA_Modified.txt").c_str());
    }
}

void FCompilerContext::PrintLists(std::ostream& os)
{
    os << std::endl << "## Compiler Data Entries ##" << std::endl;
    if (!MoleculeList.empty()) {
        os << "  MoleculeList: " << std::endl << "  " << "  ";
        for (int i = 0; i < MoleculeList.size(); i++) {
            os << "[" << i << "] " << MoleculeList[i]->Name << ", ";
        } os << std::endl;
    }

    if (!PathwayList.empty()) {
        os << "  PathwayList: " << std::endl << "  " << "  ";
        for (int i = 0; i < PathwayList.size(); i++) {
            os << "[" << i << "] " << PathwayList[i]->Name << ", ";
        } os << std::endl;
    }

    if (!ReactionList.empty()) {
        os << "  ReactionList: " << std::endl << "  " << "  ";
        for (int i = 0; i < ReactionList.size(); i++) {
            os << "[" << i << "] " << ReactionList[i]->Name << ", ";
        } os << std::endl;
    }

    if (!ContainerList.empty()) {
        os << "  ContainerList: " << std::endl << "  " << "  ";
        for (int i = 0; i < ContainerList.size(); i++) {
            os << "[" << i << "] " << ContainerList[i]->Name << ", ";
        } os << std::endl;
    }

    if (!CountList.empty()) {
        os << "  CountList: " << std::endl << "  " << "  ";
        for (int i = 0; i < CountList.size(); i++) {
            os << "[" << i << "] " << CountList[i]->Name << ", ";
        } os << std::endl;
    }

    if (!LocationList.empty()) {
        os << "  LocationList: " << std::endl << "  " << "  ";
        for (int i = 0; i < LocationList.size(); i++) {
            os << "[" << i << "] " << LocationList[i]->Name << ", ";
        } os << std::endl;
    }
}

void FCompilerContext::SaveUsingModuleList(const char* Filename)
{
    ofstream fp(Filename);

    if (!fp) {
        cerr << "Can't open file: " << Filename << endl;
        return;
    }

    vector<string> Headers;
    Headers.push_back("ProcessName");

    // Print Headers
    for(int i = 0; i < Headers.size() - 1; i++) {
        const auto& header = Headers[i];
        fp << "\"" << header << "\"" << "\t";
    }
    fp << "\"" << Headers.back() << "\"" << endl;

    // Items
    for(const auto& ModuleName: UsingModuleList) {
        fp << "\"" << ModuleName << "\"" << endl;
    }

}

void FCompilerContext::AddToMoleculeList(FMolecule *NewMolecule)
{
    bool Addition = true;
    // std::cout<< "Checking if " << NewMolecule->Name << "Exists in MoleculeList" << std::endl;
    for (auto& molecule : MoleculeList) {
        if (molecule->Name == NewMolecule->Name) {
            Addition = false;
            std::cout << "Redundant molecule name found in MoleculeList :" + NewMolecule->Name << std::endl;

            // overwrite if higher molecule class (TODO::cover more molecule classes)
            if ((!Utils::is_class_of<FEnzyme, FMolecule>(molecule)) &
               (  Utils::is_class_of<FEnzyme, FMolecule>(NewMolecule))) {
                molecule = NewMolecule;
                std::cout << "\tOverwritten in MoleculeList :" + NewMolecule->Name << std::endl;
            }

            break;
        }
    }
    if (Addition) {
        MoleculeList.push_back(NewMolecule);
    }
//        std::cout << "Molecule "<< NewMolecule->Name << " has been added to the system" << std::endl;
}

void FCompilerContext::AddToReactionList(FReaction *NewReaction)
{
    bool Addition = true;
    // std::cout<< "Checking if " << NewReaction->Name << "Exists in ReactionList" << std::endl;
    for (auto& reaction : ReactionList) {
        if (reaction->Name == NewReaction->Name) {
            Addition = false;
            std::cout << "Redundant molecule name found in ReactionList :" + NewReaction->Name << std::endl;
            break;
        }
    }
    if (Addition) {
        ReactionList.push_back(NewReaction);
    }
//        std::cout << "Reaction "<< NewReaction->Name << " has been added to the system" << std::endl;
}

void FCompilerContext::AddToPathwayList(FPathway *NewPathway)
{
    bool Addition = true;
    // std::cout<< "Checking if " << NewReaction->Name << "Exists in ReactionList" << std::endl;
    for (auto& reaction : ReactionList) {
        if (reaction->Name == NewPathway->Name) {
            Addition = false;
            std::cout << "Redundant pathway name found in NewPathway :" + NewPathway->Name << std::endl;
            break;
        }
    }
    if (Addition) {
        PathwayList.push_back(NewPathway);
    }
//        std::cout << "Pathway "<< NewPathway->Name << " has been added to the system" << std::endl;
}

void FCompilerContext::AddToContainerList(FContainer *NewContainer)
{
    bool Addition = true;
    // std::cout<< "Checking if " << NewContainer->Name << "Exists in ContainerList" << std::endl;
    for (auto& container : ContainerList) {
        if (container->Name == NewContainer->Name) {
            Addition = false;
            std::cout << "Redundant molecule name found in ContainerList :" + NewContainer->Name << std::endl;
            break;
        }
    }
    if (Addition) {
        ContainerList.push_back(NewContainer);
    }
//        std::cout << "Container "<< NewContainer->Name << " has been added to the system" << std::endl;
}

void FCompilerContext::AddToCountList(FCount *NewCount)
{   // special: this one does not check redundant entries
    CountList.push_back(NewCount);
//        std::cout << "Count "<< NewCount->Name << " has been added to the system" << std::endl;
}

void FCompilerContext::AddToLocationList(FLocation *NewLocation)
{   // special: this one does not check redundant entries
    LocationList.push_back(NewLocation);
//        std::cout << "Location "<< NewLocation->Name << " has been added to the system" << std::endl;
}

void FCompilerContext::AddToMotilityList(FMotility *NewMotility)
{   // special: this one does not check redundant entries
    MotilityList.push_back(NewMotility);
//        std::cout << "Motility "<< NewMotility->Name << " has been added to the system" << std::endl;
}

void FCompilerContext::Organize()
{
    std::cout<< std::endl << "## Organizing Compiler Data ## " << std::endl;


    if (!GetSubList_MoleculeList("Polymerase").empty()) {
        ApplyDefaultGeneticInformationProcessingOnMoleculeList();
    }
//    MakeListsFromMoleculeList();
//    MakeListsFromContainerList();
                MakeListsFromMotilityList(); // temporary

    AssignReactionTypesForReactionList();
    AdjustMolarity_PseudoMolecule();
}

void FCompilerContext::ApplyDefaultGeneticInformationProcessingOnMoleculeList() {
    std::cout << std::endl << "  Applying Genetic Information Processing on MoleculeList..." << std::endl;

    std::vector<FMolecule *> ListOfNewMolecules;

    bool bDNAP = !GetSubList_MoleculeList("DNAP").empty();
    bool bRNAP = !GetSubList_MoleculeList("RNAP").empty();
    bool bRibosome = !GetSubList_MoleculeList("Ribosome").empty();

    // Generate chromosome if not made by the user
    if ((bDNAP || bRNAP) & GetSubList_MoleculeList("Chromosome").empty()) {
        int ChromosomeSize = -1;

        // E coli specific chromosome
        auto Organisms = GetSubList_ContainerList("Organism");
        for (auto& organism : Organisms) {
            if (dynamic_cast<FOrganism *>(organism)->Species == "Ecoli") {
                //std::string Strain == "K-12 MG1655";
                ChromosomeSize = 4641652;
                FGeneticMaterial* Chromosome = GenerateChromosome("Ch_I", 1, ChromosomeSize);
                Chromosome->SetTemplate(Chromosome->Name);
                ListOfNewMolecules.push_back(Chromosome);
                break;
            }
        }

        if (ChromosomeSize < 0) {
            ChromosomeSize = 500000;
            FGeneticMaterial* Chromosome = GenerateChromosome("Chromosome", 1, ChromosomeSize);
            Chromosome->SetTemplate(Chromosome->Name);
            ListOfNewMolecules.push_back(Chromosome);
        }
    }

    //
    if (bRNAP & bRibosome) {
        // add tRNA, rRNA, misc RNA
        for (auto &molecule: MoleculeList) {
            if (Utils::is_class_of<FProtein, FMolecule>(molecule)) {
                auto Protein = dynamic_cast<FProtein *>(molecule);

                FGeneticMaterial* mRNA = GenerateCounterpart_RNA(molecule->Name, 0, "mRNA");
                FGeneticMaterial* Gene = GenerateCounterpart_Gene(molecule->Name, 1);
                Protein->SetTemplate(mRNA->Name);
                mRNA->SetTemplate(Gene->Name);
                ListOfNewMolecules.push_back(Gene);
                ListOfNewMolecules.push_back(mRNA);

            } else if (Utils::is_class_of<FRNA, FMolecule>(molecule)) {
                auto RNA = dynamic_cast<FRNA *>(molecule);

                FGeneticMaterial* Gene = GenerateCounterpart_Gene(molecule->Name, 1);
                RNA->SetTemplate(Gene->Name);
                ListOfNewMolecules.push_back(Gene);
            }
        }
    } else if (bRNAP & !bRibosome) {
        for (auto &molecule: MoleculeList) {
            if (Utils::is_class_of<FRNA, FMolecule>(molecule)) {
                auto RNA = dynamic_cast<FRNA *>(molecule);

                FGeneticMaterial* Gene = GenerateCounterpart_Gene(molecule->Name, 1);
                RNA->SetTemplate(Gene->Name);
                ListOfNewMolecules.push_back(Gene);
            }
        }
    } else if (!bRNAP & bRibosome) {
        for (auto &molecule: MoleculeList) {
            if (Utils::is_class_of<FProtein, FMolecule>(molecule)) {
                auto Protein = dynamic_cast<FProtein *>(molecule);

                FGeneticMaterial * mRNA = GenerateCounterpart_RNA(molecule->Name, 1, "mRNA");
                Protein->SetTemplate(mRNA->Name);
                ListOfNewMolecules.push_back(mRNA);
            }
        }
    }

    if (!ListOfNewMolecules.empty()) {
        for (auto& molecule : ListOfNewMolecules) {
            if (Option.bDebug) {
                std::string TextOutputTag = "[Counterpart Molecule Added] ";
                molecule->Print(std::cout);
            }
            AddToMoleculeList(molecule);
        }
    }
}

FGeneticMaterial * FCompilerContext::GenerateChromosome(std::string MolName, int Count, int Size)
{
    FChromosome * NewChromosome = new FChromosome(MolName, Size);

    FCount * NewCount = new FCount(MolName, Count);
    AddToCountList(NewCount);

    return NewChromosome;
}

FGeneticMaterial * FCompilerContext::GenerateCounterpart_Gene(std::string MolName, int Count)
{
    std::string NewName = MolName + "_Gene";
    FGene * NewGene = new FGene(NewName);

    FCount * NewCount = new FCount(NewName, Count);
    AddToCountList(NewCount);

    return NewGene;
}

FGeneticMaterial *  FCompilerContext::GenerateCounterpart_RNA(std::string MolName, int Count, std::string RNAType)
{
    std::string NewName = MolName + "_" + RNAType;
    FRNA * NewRNA;

    if (RNAType == "mRNA") {
        NewRNA = new FmRNA(NewName);
        FCount * NewCount = new FCount(NewName, Count);
    } else if (RNAType == "rRNA") {
        NewRNA = new FrRNA(NewName);
        FCount * NewCount = new FCount(NewName, Count);
    } else if (RNAType == "tRNA") {
        NewRNA = new FtRNA(NewName);
        FCount * NewCount = new FCount(NewName, Count);
    } else if (RNAType == "miscRNA") {
        NewRNA = new FmiscRNA(NewName);
        FCount * NewCount = new FCount(NewName, Count);
    } else {
        NewRNA = new FRNA(NewName);
        FCount * NewCount = new FCount(NewName, Count);
    }

    return NewRNA;
}

void FCompilerContext::MakeListsFromMotilityList()
{
    for (auto& motility : MotilityList) {
        if (Utils::is_class_of<FSwim, FMotility>(motility)) {
            auto swim = dynamic_cast<FSwim *>(motility);
            for (auto threshold : swim->Thresholds) {
                ThresholdList.push_back(threshold);
            }
        }
    }
}

enum ReactionTypeAssignment {
    //Standard Reactions
    Standard_Unregulated = 0,
    Standard_Inhibition_Allosteric = 1,
    Standard_Activation_Allosteric = 2,

    // Enzymatic Reactions:Standard
    Enz_Standard_Unregulated = 10,
    Enz_Standard_Inhibition_Allosteric = 11,
    Enz_Standard_Activation_Allosteric = 12,

    // Enzymatic Reactions:MichaelisMenten
    Enz_MichaelisMenten_Unregulated = 20,
    Enz_MichaelisMenten_Inhibition_Allosteric = 21,
    Enz_MichaelisMenten_Inhibition_Competitive = 22,
    Enz_MichaelisMenten_Activation_Allosteric = 23,

    // Transporter Reactions
    Transporter_Unregulated = 40,
    Transporter_Inhibition_Allosteric = 41,
    Transporter_Inhibition_Competitive = 41,
    Transporter_Activation_Allosteric = 43,

    // Regulatory Reactions
    Regulatory_Inhibition_Allosteric = 100,
    Regulatory_Inhibition_Competitive = 101,
    Regulatory_Activation_Allosteric = 102,

    // Polymerase Reactions
    Polymerase_TemplateBased_DNAP = 200,
    Polymerase_TemplateBased_RNAP = 201,
    Polymerase_TemplateBased_Ribosome = 202,

    Unassigned = -1,
    // Polymerase reactions to add

};

void FCompilerContext::AssignReactionType(FReaction *Reaction, int Regulation)
{
    ReactionTypeAssignment Type = Unassigned;

    if (Utils::is_class_of<FEnzymaticReaction, FReaction>(Reaction)) {
        auto reaction = dynamic_cast<FEnzymaticReaction *>(Reaction);
//        auto Enzyme = GetEnzyme_EnzymeList(reaction->Enzyme);

        if (Utils::is_class_of<FEnz_StandardReaction, FReaction>(Reaction)) {
            if      (Regulation == 0) { Type = Enz_Standard_Unregulated; }
            else if (Regulation == 1) { Type = Enz_Standard_Inhibition_Allosteric; }
            else if (Regulation == 3) { Type = Enz_Standard_Activation_Allosteric; }
        } else {
            if      (Regulation == 0) { Type = Enz_MichaelisMenten_Unregulated; }
            else if (Regulation == 1) { Type = Enz_MichaelisMenten_Inhibition_Allosteric; }
            else if (Regulation == 2) { Type = Enz_MichaelisMenten_Inhibition_Competitive; }
            else if (Regulation == 3) { Type = Enz_MichaelisMenten_Activation_Allosteric; }
        }

    } else if (Utils::is_class_of<FRegulatoryReaction, FReaction>(Reaction)) {
        auto reaction = dynamic_cast<FRegulatoryReaction *>(Reaction);
        auto& Effect = reaction->Effect;
        auto& n = reaction->n;

        if      ((Effect == "Inhibition") & (n >= 0))  { Type = Regulatory_Inhibition_Allosteric; }
        else if ((Effect == "Inhibition") & (n == -1)) { Type = Regulatory_Inhibition_Competitive; }
        else if ((Effect == "Activation") & (n >= 0))  { Type = Regulatory_Activation_Allosteric; }

    // Some reactions may have similar form to FStandardReaction, therefore more specific reactions are assigned before FStandardReaction.
    } else if ((Utils::is_class_of<FStandardReaction, FReaction>(Reaction)) & !(Utils::is_class_of<FEnzymaticReaction, FReaction>(Reaction))) {

        if      (Regulation == 0) { Type = Standard_Unregulated; }
        else if (Regulation == 1) { Type = Standard_Inhibition_Allosteric; }
        else if (Regulation == 3) { Type = Standard_Activation_Allosteric; }

    } else if (Utils::is_class_of<FTransporterReaction, FReaction>(Reaction)) {

        if      (Regulation == 0) { Type = Transporter_Unregulated; }
        else if (Regulation == 1) { Type = Transporter_Inhibition_Allosteric; }
        else if (Regulation == 2) { Type = Transporter_Inhibition_Competitive; }
        else if (Regulation == 3) { Type = Transporter_Activation_Allosteric; }

    } else if (Utils::is_class_of<FPolymeraseReaction, FReaction>(Reaction)) {

        // Template-based polymerase reactions

        auto PolymeraseReaction = dynamic_cast<FPolymeraseReaction *>(Reaction);
        auto Polymerase = dynamic_cast<FPolymerase *>(GetMolecule_MoleculeList(PolymeraseReaction->Polymerase));

        if      (Utils::is_class_of<FDNAP, FPolymerase>(Polymerase))        { Type = Polymerase_TemplateBased_DNAP; }
        else if (Utils::is_class_of<FRNAP, FPolymerase>(Polymerase))        { Type = Polymerase_TemplateBased_RNAP; }
        else if (Utils::is_class_of<FRibosome, FPolymerase>(Polymerase))    { Type = Polymerase_TemplateBased_Ribosome; }

    } else {
        Utils::Assertion(false, "ReactionTypeAssignment Failed: " + Reaction->Name);
    }

    Reaction->Type = Type;
    std::cout << " | Reaction Type: " << Type << std::endl;

}

void FCompilerContext::AssignReactionTypesForReactionList()
{
    std::cout << endl << "  Assigning Reaction Types..." << endl;
    int i = 0;
    // check if regulated
    int reg;

    // names of the reactions regulated by regulatory reactions
    std::vector<std::string> Inhibition_Allosteric = GetNames_ReactionList("Regulatory_Inhibition_Allosteric");
    std::vector<std::string> Inhibition_Competitive = GetNames_ReactionList("Regulatory_Inhibition_Competitive");
    std::vector<std::string> Activation_Allosteric = GetNames_ReactionList("Regulatory_Activation_Allosteric");

    for (auto& reaction : ReactionList) {

        // skip if assigned already
        if (reaction->Type >= 0) {
            continue;
            }

        // Unregulated
        if ((std::find(Inhibition_Allosteric.begin(), Inhibition_Allosteric.end(), reaction->Name) == Inhibition_Allosteric.end()) &
            (std::find(Inhibition_Competitive.begin(), Inhibition_Competitive.end(), reaction->Name) == Inhibition_Competitive.end()) &
            (std::find(Activation_Allosteric.begin(), Activation_Allosteric.end(), reaction->Name) == Activation_Allosteric.end())) {
            reg = 0;

        // Inhibition_Allosteric
        } else if (std::find(Inhibition_Allosteric.begin(), Inhibition_Allosteric.end(), reaction->Name) != Inhibition_Allosteric.end()) {
            reg = 1;

        // Inhibition_Competitive
        } else if (std::find(Inhibition_Competitive.begin(), Inhibition_Competitive.end(), reaction->Name) != Inhibition_Competitive.end()) {
            reg = 2;

        // Activation_Allosteric
        } else if (std::find(Activation_Allosteric.begin(), Activation_Allosteric.end(), reaction->Name) != Activation_Allosteric.end()) {
            reg = 3;
        }

        std::cout << "    [" << i << "] " << reaction->Name << "\t| Regulation mechanism :" << reg;
        AssignReactionType(reaction, reg);

        i++;
    }
}

void FCompilerContext::MakeListsFromMoleculeList()
{
//    std::cout << "Making Context Lists from Molecule List..." << std::endl;
//    for (auto& molecule : MoleculeList) {
//        if (Utils::is_class_of<FEnzyme, FMolecule>(molecule)) {
//            auto enzyme = dynamic_cast<FEnzyme *>(molecule);
//            EnzymeList.push_back(enzyme);
//// add more as needed
////        } else if {
//        }
//    }
}

void FCompilerContext::MakeListsFromContainerList()
{
//    std::cout << "Making Context Lists from Container List..." << std::endl;
//    for (auto& container : ContainerList) {
//        if (Utils::is_class_of<FOrganism, FContainer>(container)) {
//            auto organism = dynamic_cast<FOrganism *>(container);
//            OrganismList.push_back(organism);
// add more as needed
//        } else if {
//        }
//    }
}

void FCompilerContext::PrintInitialLocations(std::ostream& os) // TODO: NEEDS UPDATE
{
    os << std::endl << "## Compiler Initial Locations ##" << std::endl;
    if (!LocationList.empty()) {
        for (int i = 0; i < LocationList.size(); i++) {
            os << "  [" << i << "] ";
            os << LocationList[i]->Name << " : (";
            for (auto& coord : LocationList[i]->Coord) {
                os << Utils::SciFloat2Str(coord) << ", ";
            } os << ")" << std::endl;
        } os << std::endl;
    } else {
        os << "  None" << std::endl;
    }
}

void FCompilerContext::PrintInitialCounts(std::ostream& os)
{
    os << std::endl << "## Compiler Initial Counts: Molecules ##" << std::endl;
    // currently restricted to the molecules that are registered on the molecule list
    if (!MoleculeList.empty()) {
        for (int i = 0; i < MoleculeList.size(); i++) {
            float Count = GetInitialCountByName_CountList(MoleculeList[i]->Name);
            os << "  [" << i << "] " << MoleculeList[i]->Name << " : " << Count << std::endl;
        }
    } else {
        os << "  None" << std::endl;
    }
    os << std::endl << "## Compiler Initial Counts: Containers ##" << std::endl;
    // currently restricted to the molecules that are registered on the molecule list
    if (!ContainerList.empty()) {
        for (int i = 0; i < ContainerList.size(); i++) {
            float Count = GetInitialCountByName_CountList(ContainerList[i]->Name);
            os << "  [" << i << "] " << ContainerList[i]->Name << " : " << Count << std::endl;
        }
    } else {
        os << "  None" << std::endl;
    }

}

std::string FCompilerContext::QueryTable(std::string Name, std::string Property, FTable Table)
{
    for (auto& record : Table.Records) {
        if (record["Name"] == Name) {
            return record[Property];
        }          
    }
//    std::cout << "Record not found in the database. Query: " << Name << " | Property Keyword: " << Property << std::endl;
    return std::string();
}

std::vector<std::string> FCompilerContext::GetNameListByType_MoleculeList(std::string Type)
{
    return GetNameList_MoleculeList(GetSubList_MoleculeList(Type));
}


std::vector<FReaction *> FCompilerContext::GetSubList_ReactionList(std::string Type, std::string NameSpace_Pathway)
{
    int ReactionType = Unassigned;

    //Standard Reactions
  
    if      (Type == "Standard_Unregulated")             { ReactionType = Standard_Unregulated; }
    else if (Type == "Standard_Inhibition_Allosteric")   { ReactionType = Standard_Inhibition_Allosteric; }
    else if (Type == "Standard_Activation_Allosteric")   { ReactionType = Standard_Activation_Allosteric; }

    // Enzymatic Reactions:Standard
    else if (Type == "Enz_Standard_Unregulated")            { ReactionType = Enz_Standard_Unregulated; }
    else if (Type == "Enz_Standard_Inhibition_Allosteric")  { ReactionType = Enz_Standard_Inhibition_Allosteric; }
    else if (Type == "Enz_Standard_Activation_Allosteric")  { ReactionType = Enz_Standard_Activation_Allosteric; }

    // Enzymatic Reactions:MichaelisMenten
    else if (Type == "Enz_MichaelisMenten_Unregulated")            { ReactionType = Enz_MichaelisMenten_Unregulated; }
    else if (Type == "Enz_MichaelisMenten_Inhibition_Allosteric")  { ReactionType = Enz_MichaelisMenten_Inhibition_Allosteric; }
    else if (Type == "Enz_MichaelisMenten_Inhibition_Competitive") { ReactionType = Enz_MichaelisMenten_Inhibition_Competitive; }
    else if (Type == "Enz_MichaelisMenten_Activation_Allosteric")  { ReactionType = Enz_MichaelisMenten_Activation_Allosteric; }

    // Transporter Reactions
    else if (Type == "Transporter_Unregulated")            { ReactionType = Transporter_Unregulated; }
    else if (Type == "Transporter_Inhibition_Allosteric")  { ReactionType = Transporter_Inhibition_Allosteric; }
    else if (Type == "Transporter_Inhibition_Competitive") { ReactionType = Transporter_Inhibition_Competitive; }
    else if (Type == "Transporter_Activation_Allosteric")  { ReactionType = Transporter_Activation_Allosteric; }

    // Regulatory Reactions
    else if (Type == "Regulatory_Inhibition_Allosteric")    { ReactionType = Regulatory_Inhibition_Allosteric; }
    else if (Type == "Regulatory_Inhibition_Competitive")   { ReactionType = Regulatory_Inhibition_Competitive; }
    else if (Type == "Regulatory_Activation_Allosteric")    { ReactionType = Regulatory_Activation_Allosteric; }

    // Polymerase Reactions
    else if (Type == "Polymerase_TemplateBased_DNAP")       { ReactionType = Polymerase_TemplateBased_DNAP; }
    else if (Type == "Polymerase_TemplateBased_RNAP")       { ReactionType = Polymerase_TemplateBased_RNAP; }
    else if (Type == "Polymerase_TemplateBased_Ribosome")   { ReactionType = Polymerase_TemplateBased_Ribosome; }

    else {
        if (Option.bDebug) {
            std::cout << "Unable to find the reaction type for requested Reaction Type: " + Type << std:: endl;
        }
    }

//    Utils::Assertion(ReactionType >= 0, "Type to ReactionType Error. Type: " + Type + " | Assigned ReactionType: " + std::to_string(ReactionType));

    std::vector<FReaction *> SubList;

    for (auto& reaction : ReactionList) {
        if (Type == "All") {
            SubList.push_back(reaction);
        } else if (reaction->Type == ReactionType) {
            SubList.push_back(reaction);
        }
    }

    if (NameSpace_Pathway != "") {
        SubList = FilterByPathway_ReactionList(SubList, NameSpace_Pathway);
    }

    return SubList;
}

std::vector<FStandardReaction *> FCompilerContext::GetList_Standard_ReactionList(std::string ReactionType)
{
    std::vector<FStandardReaction *> SubList;
    std::vector<FReaction *> ReactionSubList = GetSubList_ReactionList(ReactionType);
    

    for (FReaction* Reaction: ReactionSubList) {
        auto reaction = dynamic_cast<FStandardReaction *>(Reaction);
        SubList.push_back(reaction);
    }

    return SubList;
}

std::vector<FRegulatoryReaction *> FCompilerContext::GetList_Regulatory_ReactionList(std::string ReactionType)
{
    // Reaction Type is either "Inhibition" or "Activation"

    std::vector<FRegulatoryReaction *> SubList;

    for (FReaction* Reaction: ReactionList) {
        if (Utils::is_class_of<FRegulatoryReaction, FReaction>(Reaction)) {
            auto reaction = dynamic_cast<FRegulatoryReaction *>(Reaction);
            if (reaction->Effect == ReactionType) {
                SubList.push_back(reaction);
            }
        }
    }

    return SubList;
}

FMolecule * FCompilerContext::GetMolecule_MoleculeList(std::string Name)
{
    for (auto& molecule : MoleculeList) {
        if (molecule->Name == Name) {
            return molecule;
        }
    }

    Utils::Assertion (false, "Unable to find enzyme in the Context.MoleculeList: " + Name);
    return nullptr;
}

std::vector<std::string> FCompilerContext::GetNames_EnzymaticReactionList(std::vector<FEnzymaticReaction *> EnzymaticReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymaticReactionList){
        for (auto& Stoich : item->Stoichiometry){
            if (std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) {
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

//std::vector<std::string> FCompilerContext::GetNames_PolymeraseList()
//{
//    std::vector<FMolecule *> MoleculeSubList = GetSubList_MoleculeList("Polymerase")
//
//    std::vector<std::string> StrList;
//    for (auto& item : PolymeraseList){
//        StrList.push_back(item->Name);
//    }
//    return StrList;
//}

//std::vector<float> FCompilerContext::GetRates_PolymeraseList(std::vector<FPolymerase *> PolymeraseList)
//{
//    std::vector<float> floatList;
//    for (auto& item : PolymeraseList){
//        floatList.push_back(item->Rate);
//    }
//    return floatList;
//}

std::vector<FPolymeraseReaction *> FCompilerContext::GetList_Polymerase_ReactionList()
{
    std::vector<FPolymeraseReaction *> SubList;

    for (FReaction* Reaction: ReactionList) {
        if (Utils::is_class_of<FPolymeraseReaction, FReaction>(Reaction)) {
            auto PolymeraseReaction = dynamic_cast<FPolymeraseReaction* >(Reaction);
            SubList.push_back(PolymeraseReaction);
        }
    }
    return SubList;
}

std::vector<std::string> FCompilerContext::GetNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *> PolymeraseReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : PolymeraseReactionList){
        for (auto& Stoich : item->Stoichiometry){
            if (std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) {
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetSubstrateNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *> PolymeraseReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : PolymeraseReactionList){
        for (auto& Stoich : item->Stoichiometry){
            if (std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) {
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetReactantNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *> PolymeraseReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : PolymeraseReactionList){
        for (auto& Stoich : item->Stoichiometry){
            if ((std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) & (Stoich.second < 0)){
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetProductNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *> PolymeraseReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : PolymeraseReactionList){
        for (auto& Stoich : item->Stoichiometry){
            if ((std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) & (Stoich.second > 0)){
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetBuildingBlockNames_PolymeraseReactionList(std::vector<FPolymeraseReaction *> PolymeraseReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : PolymeraseReactionList){
        for (auto& BuildingBlock : item->BuildingBlocks){
            StrList.push_back(BuildingBlock);
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetNames_ContainerList(std::string Type)
{
    std::vector<FContainer *> SubList = GetSubList_ContainerList(Type);
    std::vector<std::string> StrList;

    for (auto& container : SubList) {
        StrList.push_back(container->Name);
    }

    return StrList;
}

std::vector<FContainer *> FCompilerContext::GetSubList_ContainerList(std::string Type)
{
    std::vector<FContainer *> SubList;

    if (Type == "Compartment") {
        for (auto& container : ContainerList) {
            if (Utils::is_class_of<FCompartment, FContainer>(container)) {
                SubList.push_back(container);
            }
        }
    } else if (Type == "Organism") {
        for (auto& container : ContainerList) {
            if (Utils::is_class_of<FOrganism, FContainer>(container)) {
                SubList.push_back(container);
            }
        }
    } else if (Type == "All") {
        for (auto& container : ContainerList) {
            SubList.push_back(container);
        }
    }
    
 
    return SubList;
}

std::vector<std::string> FCompilerContext::GetNames_ReactionList(std::string Type)
{
    std::vector<std::string> StrList;

//    if (Type == "All") {
//        for (auto& reaction : ReactionList){
//            StrList.push_back(reaction->Name);       
//        return StrList;
//        }
//
//    std::vector<std::string> ReactionNames_All = GetNames_ReactionList("All"); // recursive

    if (Type.rfind("Regulatory", 0) == 0) {
        for (auto& Reaction : ReactionList){
            if (Utils::is_class_of<FRegulatoryReaction, FReaction>(Reaction)) {
                auto reaction = dynamic_cast<FRegulatoryReaction *>(Reaction);

                if ((Type.find("Inhibition") != std::string::npos) & (Type.find("Allosteric") != std::string::npos)) {
                    if ((reaction->Effect == "Inhibition") & (reaction->Mode == "Allosteric")) {
                        for (auto& stoich : reaction->Stoichiometry) {
                            if (stoich.second > 0) {
                                StrList.push_back(stoich.first);
                            }
                        }
                    }
                } else if ((Type.find("Inhibition") != std::string::npos) & (Type.find("Competitive") != std::string::npos)) {
                    if ((reaction->Effect == "Inhibition") & (reaction->Mode == "Competitive")) {
                        for (auto& stoich : reaction->Stoichiometry) {
                            if (stoich.second > 0) {
                                StrList.push_back(stoich.first);
                            }
                        }
                    }
                } else if ((Type.find("Activation") != std::string::npos) & (Type.find("Allosteric") != std::string::npos)) {
                    if ((reaction->Effect == "Activation") & (reaction->Mode == "Allosteric")) {
                        for (auto& stoich : reaction->Stoichiometry) {
                            if (stoich.second > 0) {
                               StrList.push_back(stoich.first);
                            }
                        }
                    }
                }
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetNames_PathwayList()
{
    std::vector<std::string> StrList;
    for (auto& item : PathwayList){
        StrList.push_back(item->Name);
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetSequences_PathwayList()
{
    std::vector<std::string> StrList;
    for (auto& item : PathwayList){
        for (auto& subitem : item->Sequence){
            if (std::find(StrList.begin(), StrList.end(), subitem) == StrList.end()) {
                StrList.push_back(subitem);
            }
        }
    }
    return StrList;
}

std::vector<int> FCompilerContext::GetCoefficientArray(FReaction* Reaction, std::vector<int> Idx_Substrates)
{
    // StoichiometryMatrix Routine
    std::vector<int> CoeffArray(Idx_Substrates.size(), 0);

    for (auto& stoich : Reaction->Stoichiometry) {
    		// std::cout << "enzymaticReaction-> Stoichiometry for loop. Now working on: " << stoich.first << endl;
        int Coeff = stoich.second;
        int MolIdx = GetIdxByName_MoleculeList(stoich.first);
    		// std::cout << "MolIdx: " << std::to_string(MolIdx) << endl;
        // get local index from Substrate list
        int Idx_Local = 0;
        for (auto& Idx_Substrate : Idx_Substrates) {
    		// std::cout << "Current Idx_Substrate: " << std::to_string(Idx_Substrate) << "\t| Idx_Substrate" << endl;
            if (Idx_Substrate == MolIdx) {
                CoeffArray[Idx_Local] = Coeff;
    		// std::cout << "Enz_Standard | SubstrateName: " << SubstrateName << ", MolIdx: " << MolIdx << ", Idx_Local: " << Idx_Local << ", Coeff: " << Coeff << endl;
            }
            Idx_Local++;
        }
    		// std::cout << "add" << endl;
    }
    return CoeffArray; 
}

//FEnzyme * FCompilerContext::GetEnzyme_EnzymeList(std::string Name)
//{
//    for (auto& enzyme : EnzymeList) {
//        if (enzyme->Name == Name) {
//            return enzyme;
//        }
//    }
//
//    Utils::Assertion (false, "Unable to find enzyme in the Context.EnzymeList: " + Name);
//    return nullptr;
//}

// float FCompilerContext::GetFloatAttributeByName_EnzymeList(std::string Name, std::string Attribute)
// {
//     float Value = -1;
// 
//     for (auto& enzyme : EnzymeList){
//         if (enzyme->Name == Name) {
//             if (Attribute == "kcat")   { Value = enzyme->kcat; }
//             else if (Attribute == "KM")     { Value = enzyme->KM; }
//             break;
//         }
//     }
// 
//     if (Value != -1) { return Value; }
//     else             { std::cout << "Undefined enzyme name: " << Name << std::endl; }
// }

// std::string FCompilerContext::GetStringAttributeByName_EnzymeList(std::vector<FEnzyme *> EnzymeList, std::string Name, std::string Attribute)
// {
//     std::string Value;
// 
//     for (auto& enzyme : EnzymeList){
//         if (enzyme->Name == Name) {
//             if      (Attribute == "Substrate")      { Value = enzyme->Substrate; } 
//             // TODO: add more after updating FEnzyme class
//             break;
//         }
//     }
// 
//     if (Value.empty()) { 
//         std::cout << "Undefined enzyme name or attribute: " << Name << ", " << Attribute << std::endl; 
//     }
//     assert (!Value.empty());
//     return Value;
// }

std::vector<std::vector<int>> FCompilerContext::GetStoichiometryMatrix(std::vector<FReaction *> ListOfReactions)
{
    // Types: "Standard_Unregulated", "Standard_Inhibition", "Standard_Activation",
    //        "Enz_Standard_Unregulated", "Enz_Standard_Inhibition", "Enz_Standard_Activation",
    //        "Enz_MichaelisMenten_Unregulated", "Enz_MichaelisMenten_Inhibition_Allosteric", "Enz_MichaelisMenten_Inhibition_Competitive", "Enz_MichaelisMenten_Inhibition_Competitive"

    std::vector<std::vector<int>> StoichMatrix;

    // for matrix generation
    std::vector<int> Idx_Substrates = GetUniqueSubstrateIdx_ReactionList(ListOfReactions);

    for (auto& reaction : ListOfReactions) {
        StoichMatrix.push_back(GetCoefficientArray(reaction, Idx_Substrates));
    }

    return StoichMatrix;
}

std::vector<int> FCompilerContext::AddUniqueSubstrateIdxToIdxList(FReaction* Reaction, std::vector<int>IdxList)
{
    for (auto& stoich : Reaction->Stoichiometry) {
        int MolIdx = GetIdxByName_MoleculeList(stoich.first);
        if (std::find(IdxList.begin(), IdxList.end(), MolIdx) == IdxList.end()) {
//std::cout << "New mol idx added for : " << stoich.first << endl;
            IdxList.push_back(MolIdx);
        }
    }
    return IdxList;
}

std::vector<int> FCompilerContext::GetUniqueSubstrateIdx_ReactionList(std::vector<FReaction *> ListOfReactions)
{
    // Types: "Standard_Unregulated", "Standard_Inhibition", "Standard_Activation",
    //        "Enz_Standard_Unregulated", "Enz_Standard_Inhibition", "Enz_Standard_Activation",
    //        "Enz_MichaelisMenten_Unregulated", "Enz_MichaelisMenten_Inhibition_Allosteric", "Enz_MichaelisMenten_Inhibition_Competitive", "Enz_MichaelisMenten_Inhibition_Competitive"

    std::vector<int> IdxList;

    for (auto& reaction : ListOfReactions) {
        IdxList = AddUniqueSubstrateIdxToIdxList(reaction, IdxList);
    }

    return IdxList;
}
//
//std::vector<std::vector<int>> FCompilerContext::GetStoichiometryMatrix_PolymeraseReaction(std::vector<FPolymeraseReaction *> PolymeraseReactionList)
//{
//    std::vector<std::vector<int>> StoichMatrix;
//    std::vector<FMolecule *> SMolList = GetSubList_MoleculeList("SmallMolecule");
//
//    for (auto& PolymeraseReaction : PolymeraseReactionList){
//        std::vector<int> CoeffArray(SMolList.size(), 0); // replace with substrate index for the reaction
//
//        for (auto& Stoich : PolymeraseReaction->Stoichiometry){
//            std::string SubstrateName = Stoich.first;
//            int Coeff = Stoich.second;
//	    int Index = 0;
//
//            for (auto& molecule : SMolList){
//                // std::cout << "Searching from the List: " << Substrate << " | " << "Index: " << Index << endl;
//                if (molecule->Name == SubstrateName){
//                    break;
//                }
//                Index++;
//            }
//            if (Index >= MoleculeList.size()) {
//	std::cout << "Substrate index searching in MoleculeList failed: " << SubstrateName << endl;
//            } else {
//            // std::cout << "Substrate Searching from the List: " << SubstrateName << " | " << "Index: " << Index << endl;
//            }
//
//            CoeffArray[Index] = Coeff;
//        }
//        StoichMatrix.push_back(CoeffArray);
//    }
//    return StoichMatrix;
//}

int FCompilerContext::GetIdxByName_MoleculeList(std::string InputName)
{
    if (InputName.empty()) {
        std::cout << "ERROR: Empty Name to search!" << InputName << std::endl;
        // exit(3);
    }

    int Index = 0;
    for (auto& molecule : MoleculeList) {
        if (molecule->Name == InputName) {
            return Index;
        }
        Index++; 
    }
    if (Option.bDebug) {
        std::cout << "** Unable to find index in MoleculeList : " << InputName << std::endl;
    }
    return 0;
}

int FCompilerContext::GetIdx_MoleculeList(FMolecule * Molecule)
{
    for (int i = 0; i < MoleculeList.size(); i++) {
        if (MoleculeList[i]->Name == Molecule->Name) {
            return i;
        }
    }

    if (Option.bDebug) {
        std::cout << "%% Unable to find index in MoleculeList : " << Molecule->Name << std::endl;
    }

    return 0;
}

std::vector<float> FCompilerContext::GetLocationByName_LocationList(std::string MolName)
{
    std::vector<float> coord;
    for (auto& location : LocationList) {
        if (location->Name == MolName) {
            return location->Coord;
        }
    }
    if (coord.empty()) {
        coord = {0, 0, 0};
        std::cout << "No location info found: " << MolName << std::endl;
    }

    return coord;
}

bool FCompilerContext::CheckDuplicates_List(std::string Type, std::string Name)
{
    bool Duplicate = false;
    std::vector<std::string> Names;

    if      (Type == "CountList")     { Names = GetNames_CountList("All"); }
    else if (Type == "LocationList")  { Names = GetNames_LocationList("All"); }
    else if (Type == "ContainerList") { Names = GetNames_ContainerList("All"); }
    
    if (std::find(Names.begin(), Names.end(), Name) != Names.end()) {
        Duplicate = true;
    }

    return Duplicate;
}

std::vector<FLocation *> FCompilerContext::GetSubList_LocationList(std::string Type)
{
    std::vector<FLocation *> SubList;
    std::vector<std::string> Names;

    if      (Type == "Molecule")    { Names = GetNameListByType_MoleculeList(Type); }
    else if (Type == "Compartment") { Names = GetNames_ContainerList(Type); }
    else if (Type == "Organism")    { Names = GetNames_ContainerList(Type); }
    else if (Type == "All")         {}

    if (!Names.empty()) {
        for (auto& item : LocationList) {
            if (std::find(Names.begin(), Names.end(), item->Name) != Names.end()) {
                SubList.push_back(item);
            }
        }
    } else {
        for (auto& item : LocationList) {
            SubList.push_back(item);
        }
    }

    return SubList;
}

std::vector<std::string> FCompilerContext::GetNames_LocationList(std::string Type)
{
    std::vector<FLocation *> SubList = GetSubList_LocationList(Type);
    std::vector<std::string> StrList; 

    for (auto& location : SubList) {
        StrList.push_back(location->Name);
    }

    return StrList;
}

std::vector<std::string> FCompilerContext::GetUniqueNames_LocationList(std::string Type)
{
    std::vector<FLocation *> SubList = GetSubList_LocationList(Type);
    std::vector<std::string> StrList; 

    for (auto& location : SubList) {
        if (std::find(StrList.begin(), StrList.end(), location->Name) == StrList.end()) {
            StrList.push_back(location->Name);
        }
    }

    return StrList;
}

int FCompilerContext::GetCounts_LocationList(std::string Type)
{
    std::vector<FLocation *> SubList = GetSubList_LocationList(Type);
    int Sum = 0;

    for (auto& location : SubList) {
        Sum += int(location->Count->Amount);
    }

    return Sum;
}

int FCompilerContext::GetCounts_ContainerList(std::string Type)
{
    std::vector<FContainer *> SubList = GetSubList_ContainerList(Type);
    int Sum = 0;

    for (auto& container : SubList) {
        for (auto& count : CountList) {
            if (container->Name == count->Name) {
                Sum += int(count->Amount);
            }
        }
    }

    return Sum;
}

std::vector<std::string> FCompilerContext::GetNames_CountList(std::string Type)
{
    std::vector<FCount *> SubList = GetSubList_CountList(Type);
    std::vector<std::string> StrList;

    for (auto& item : SubList){
        StrList.push_back(item->Name);
    }
    return StrList;
}

std::vector<FCount *> FCompilerContext::GetSubList_CountList(std::string Type)
{
    std::vector<FCount *> SubList;
    std::vector<std::string> Names;

    if     ((Type == "Molecule")
         || (Type == "Restore")
         || (Type == "Event"))     { Names = GetNameListByType_MoleculeList("All"); }

    else if (Type == "Compartment") { Names = GetNames_ContainerList(Type); }
    else if (Type == "Organism")    { Names = GetNames_ContainerList(Type); }

    else if (Type == "All")         {}

    if (!Names.empty()) {
        for (auto& item : CountList) {
            if (std::find(Names.begin(), Names.end(), item->Name) != Names.end()) {
                SubList.push_back(item);
            }
        }
    } else {
        for (auto& item : CountList) {
            SubList.push_back(item);
        }
    }

    std::vector<FCount *> Filtered;
    if (Type == "Restore") {
        for (auto& count : SubList) {
            if ((count->End == -1) || (count->Name == "qL")) {
                Filtered.push_back(count);
            }
        }
        SubList = Filtered;
    } else if (Type == "Event") {
        for (auto& count : SubList) {
            if ((count->End >= 0) & (count->Begin != count->End)) {
                Filtered.push_back(count);
            }
        }
        SubList = Filtered;
    }

    return SubList;
}

float FCompilerContext::GetInitialCountByName_CountList(std::string MolName)
{
    // returns 0 if there is no initial count set
    float sum = 0;
    for (auto& count : CountList) {
        if ((count->Name == MolName) & (count->Begin == 0)) {
            sum += count->Amount;
        }
    }
    return sum;
}

bool FCompilerContext::GetMolarityFactorByName_CountList(std::string MolName)
{
    for (auto& count : CountList) {
        if (count->Name == MolName) {
            return count->bMolarity;
        }
    }
    return false;
}

bool FCompilerContext::CheckMolarityFactorTrueForAny_CountList()
{
    for (auto& count : CountList) {
        if (count->bMolarity) {
            return true; // note: PseudoMolecule Molarity is set to false by default.
        }
    }
    return false;
}

void FCompilerContext::RevertMolarity_CountList(std::string InName)
{
    for (auto& count : CountList) {
        if (count->Name == InName) {
            // adjust amount by switching unit
            if (count->bMolarity) { count->Amount /= Numbers::GetAvogadro(); }
            else                  { count->Amount *= Numbers::GetAvogadro(); }
            // revert molarity
            count->bMolarity = !(count->bMolarity);
        }
    }
}

void FCompilerContext::AdjustMolarity_PseudoMolecule()
{
    if (CheckMolarityFactorTrueForAny_CountList()) {
        RevertMolarity_CountList("Pseudo"); // Name for PseudoMolecule
    }
}

std::vector<FMolecule *> FCompilerContext::GetSubList_MoleculeList(std::string Type)
{
    std::vector<FMolecule *> SubList;

    if ((Type == "All") || (Type == "Molecule")) {
        SubList = MoleculeList;
    } else if (Type == "Chromosome") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FChromosome, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "Gene") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FGene, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "RNA") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FRNA, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "mRNA") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FmRNA, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "Protein") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FProtein, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "Enzyme") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FEnzyme, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "Polymerase") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FPolymerase, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "DNAP") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FDNAP, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "RNAP") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FRNAP, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "Ribosome") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FRibosome, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    } else if (Type == "SmallMolecule") {
        for (FMolecule *molecule: MoleculeList) {
            if (Utils::is_class_of<FSmallMolecule, FMolecule>(molecule)) {
                SubList.push_back(molecule);
            }
        }
    }

    return SubList;
}

std::vector<int> FCompilerContext::GetIdxListByType_MoleculeList(std::string Type)
{
    std::vector<int> IdxList;
    std::vector<FMolecule *> SubList = GetSubList_MoleculeList(Type);

    for (FMolecule* molecule :SubList) {
        IdxList.push_back(GetIdxByName_MoleculeList(molecule->Name));
    }

    return IdxList;
}

std::vector<int> FCompilerContext::GetIdxList_MoleculeList(std::vector<FMolecule *> ListOfMolecules)
{
    std::vector<int> IdxList;

    for (FMolecule* molecule :ListOfMolecules) {
        IdxList.push_back(GetIdx_MoleculeList(molecule));
    }

    return IdxList;
}

std::vector<std::string> FCompilerContext::GetNameList_MoleculeList(std::vector<FMolecule *> ListOfMolecules)
{
    std::vector<std::string> StrList;

    for (FMolecule* molecule :ListOfMolecules) {
        StrList.push_back(molecule->Name);
    }

    return StrList;
}

std::vector<int> FCompilerContext::GetLocalIdxList_MoleculeList(std::vector<FMolecule *> SubListOfMolecules, std::vector<FMolecule *> ListOfMolecules)
{
    std::vector<int> IdxList;

    for (FMolecule* molecule : SubListOfMolecules) {
        for (int i = 0; i < ListOfMolecules.size(); i++) {
            if (ListOfMolecules[i]->Name == molecule->Name) {
                IdxList.push_back(i);
            }
        }
    }

    return IdxList;
}

FReaction * FCompilerContext::GetReactionByPolymeraseName_ReactionList(std::string InPolymeraseName)
{
    FReaction * PolymeraseReaction;
    std::vector<FReaction *> ListOfPolymeraseReactions;

    auto Reactions_DNAP = GetSubList_ReactionList("Polymerase_TemplateBased_DNAP");
    auto Reactions_RNAP = GetSubList_ReactionList("Polymerase_TemplateBased_RNAP");
    auto Reactions_Ribosome = GetSubList_ReactionList("Polymerase_TemplateBased_Ribosome");

    ListOfPolymeraseReactions.insert(ListOfPolymeraseReactions.end(), Reactions_DNAP.begin(), Reactions_DNAP.end());
    ListOfPolymeraseReactions.insert(ListOfPolymeraseReactions.end(), Reactions_RNAP.begin(), Reactions_RNAP.end());
    ListOfPolymeraseReactions.insert(ListOfPolymeraseReactions.end(), Reactions_Ribosome.begin(), Reactions_Ribosome.end());

    for (auto& reaction : ListOfPolymeraseReactions){
        auto Reaction = dynamic_cast<FPolymeraseReaction *>(reaction);
        if (Reaction->Polymerase == InPolymeraseName) {
            PolymeraseReaction = Reaction;
        }
    }

    return PolymeraseReaction;
}

//std::vector<int> FCompilerContext::GetIdx_PolymeraseReactionSubstrate_ByPolymeraseName_MoleculeList(std::string InPolymeraseName)
//{
//    std::vector<int> IndexArray;
//    int Index;
//
//    std::vector<FPolymeraseReaction *> PolymeraseReactionList = GetList_Polymerase_ReactionList();
//
//    for (auto& PolymeraseReaction : PolymeraseReactionList){
//        if (PolymeraseReaction->Polymerase == InPolymeraseName){
//            for (auto& stoich : PolymeraseReaction->Stoichiometry) {
//                Index = GetIdxByName_MoleculeList(stoich.first);
//                IndexArray.push_back(Index);
//            }
//        }
//    }
//    return IndexArray;
//}

std::vector<int> FCompilerContext::GetIdxByStrList_MoleculeList(std::vector<std::string> StrList)
{ 
    std::vector<int> IndexArray;
    int Index;

    for (std::string Item : StrList){
        Index = GetIdxByName_MoleculeList(Item);
        IndexArray.push_back(Index);
    }
    return IndexArray;
}

std::vector<float> FCompilerContext::GetFreqMatrixForChromosomes()
{
	return std::vector<float>();
}

std::vector<float> FCompilerContext::GetFreqMatrixForRNAs()
{
	return std::vector<float>();
}

std::vector<float> FCompilerContext::GetFreqMatrixForProteins()
{
	return std::vector<float>();
}

// temporary namespace functions
std::vector<FReaction *> FCompilerContext::FilterByPathway_ReactionList(std::vector<FReaction *> ListOfReactions, std::string PathwayName)
{
    std::vector<FReaction *> SubList;

    for (auto& reaction : ListOfReactions) {
        if (reaction->Pathway == PathwayName) {
            SubList.push_back(reaction);
        }
    }

    return SubList;
}

std::string FCompilerContext::GetAssociatedPathway_ReactionList(FMolecule* Molecule)
{
    std::string Pathway;
    for (auto& reaction : ReactionList) {
        if ((reaction->CheckIfProduct(Molecule->Name)) || reaction->CheckIfProduct(Molecule->Name)) {
            Pathway = reaction->Pathway;
        }
    }
    Utils::Assertion(!Pathway.empty(), "Associated Pathway is not found for Molecule: " + Molecule->Name);

    return Pathway;
}

std::vector<FMolecule *> FCompilerContext::GetSubListByReactionList_MoleculeList(std::vector<FReaction *> ListOfReactions)
{
    std::vector<FMolecule *> Molecules;

    for (auto& reaction : ListOfReactions) {
        for (auto& stoich: reaction->Stoichiometry) {
            auto& MolName = stoich.first;
            for (auto& Molecule : MoleculeList) {
                if (Molecule->Name == MolName) {
                    if (std::find(Molecules.begin(), Molecules.end(), Molecule) == Molecules.end()) {
                        Molecules.push_back(Molecule);
//                        std::cout << "MoleculeSubList | New molecule added: " << Molecule->Name << std::endl;
                        break;
                    } else {
//                        std::cout << "MoleculeSubList | Redundant molecule skipped: " << Molecule->Name << std::endl;
                    }
                }
            }
        }
    }
    return Molecules;
}
