#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <iterator>
#include <vector>

#include "context.h"
#include "option.h"

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

void FCompilerContext::Init(const FOption& InOption)
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
    os << std::endl << "## Compiler Context Lists ##" << std::endl;
    int i = 0;
    if (!MoleculeList.empty()) {
        os << "  MoleculeList: " << std::endl << "  " << "  ";
        for (auto& item : MoleculeList){
            os << "[" << i << "] " << item->Name << ", ";
            i++;
        }
        os << std::endl;
    }

    i = 0;
    if (!PathwayList.empty()) { 
        os << "  PathwayList: " << std::endl << "  " << "  ";
        for (auto& item : PathwayList){
            os << "[" << i << "] " << item.Name << ", ";
            i++;
        }
        os << std::endl;
    }
//    os << "  EnzymeList: " << std::endl << "  " << "  ";
//    for (auto& item : EnzymeList){
//        os << item.Name << ", ";
//    };
//    os << std::endl;
    i = 0;
    if (!ReactionList.empty()) {
        os << "  ReactionList: " << std::endl << "  " << "  ";
        for (auto& item : ReactionList){
            os << "[" << i << "] " << item->Name << ", ";
            i++;
        }
        os << std::endl;
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
    for (auto& molecule :MoleculeList) {
        if (molecule->Name == NewMolecule->Name){
            Addition = false;
//            std::cout << "Redundant molecule " << NewMolecule->Name << " Found in MoleculeList" << std::endl;
            break;
        }
    }
    if (Addition) {
        MoleculeList.push_back(NewMolecule);
//        std::cout << "Molecule "<< NewMolecule->Name << " has been added to the system" << std::endl;
    }
}

void FCompilerContext::AddToReactionList(FReaction *NewReaction)
{
    bool Addition = true;
    // std::cout<< "Checking if " << NewMolecule->Name << "Exists in MoleculeList" << std::endl;
    for (auto& reaction : ReactionList) {
        if (reaction->Name == NewReaction->Name){
            Addition = false;
//            std::cout << "Redundant reaction (" << NewReaction->Name << " Found in ReactionList" << std::endl;
            break;
        }
    }
    if (Addition) {
        ReactionList.push_back(NewReaction);
    }
}

void FCompilerContext::Organize() 
{
    MakeListsFromMoleculeList();
    // Dependency Note: AssignReactionTypes uses Lists made from above
    AssignReactionTypesForReactionList();
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

    // Regulatory Reactions
    Regulatory_Inhibition_Allosteric = 100,
    Regulatory_Inhibition_Competitive = 101,
    Regulatory_Activation_Allosteric = 102,

    Unassigned = -1,
    // Polymerase reactions to add

};

void FCompilerContext::AssignReactionType(FReaction *Reaction, int Regulation)
{
    ReactionTypeAssignment Type;
    if (Utils::is_class_of<FStandardReaction, FReaction>(Reaction)) {

        if      (Regulation == 0) { Type = Standard_Unregulated; }
        else if (Regulation == 1) { Type = Standard_Inhibition_Allosteric; }
        else if (Regulation == 3) { Type = Standard_Activation_Allosteric; }

    } else if (Utils::is_class_of<FEnzymaticReaction, FReaction>(Reaction)) {
        auto reaction = dynamic_cast<FEnzymaticReaction *>(Reaction);
        auto Enzyme = GetEnzyme_EnzymeList(reaction->Enzyme);
       
        if (Enzyme->Type == 0) {
            if      (Regulation == 0) { Type = Enz_Standard_Unregulated; }
            else if (Regulation == 1) { Type = Enz_Standard_Inhibition_Allosteric; }
            else if (Regulation == 3) { Type = Enz_Standard_Activation_Allosteric; }
        }
        else if (Enzyme->Type == 1) {
            if      (Regulation == 0) { Type = Enz_MichaelisMenten_Unregulated; }
            else if (Regulation == 1) { Type = Enz_MichaelisMenten_Inhibition_Allosteric; }
            else if (Regulation == 2) { Type = Enz_MichaelisMenten_Inhibition_Competitive; }
            else if (Regulation == 3) { Type = Enz_MichaelisMenten_Activation_Allosteric; }
        }
    } else if (Utils::is_class_of<FRegulatoryReaction, FReaction>(Reaction)) {
        auto reaction = dynamic_cast<FRegulatoryReaction *>(Reaction);
        auto Effect = reaction->Effect;
        auto n = reaction->n;

	if      ((Effect == "Inhibition") & (n >= 0))  { Type = Regulatory_Inhibition_Allosteric; }
	else if ((Effect == "Inhibition") & (n == -1)) { Type = Regulatory_Inhibition_Competitive; }
	else if ((Effect == "Activation") & (n >= 0))  { Type = Regulatory_Activation_Allosteric; }
    } else {
        Type = Unassigned; // keep the default value
    }
    std::cout << " | Reaction Type: " << Type << std::endl;
    Reaction->Type = Type;
}

void FCompilerContext::AssignReactionTypesForReactionList()
{
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

        std::cout << reaction->Name << " | Regulation mechanism :" << reg;
        AssignReactionType(reaction, reg);
    }
}

void FCompilerContext::MakeListsFromMoleculeList()
{
    for (auto& molecule : MoleculeList) {
        if (Utils::is_class_of<FEnzyme, FMolecule>(molecule)) {
        auto enzyme = dynamic_cast<FEnzyme *>(molecule);
            EnzymeList.push_back(enzyme);
// add more as needed
//        } else if {
        }
    }
}


void FCompilerContext::PrintInitialCounts(std::ostream& os) 
{
    os << std::endl << "## Compiler Initial Counts ##" << std::endl;
    if (!MoleculeList.empty()) {
        for (auto& molecule : MoleculeList){
            float Count = GetInitialCountByName_CountList(molecule->Name);
            os << molecule->Name << " : " << Count << "\t| ";
        }
        os << std::endl;
    }
}

const std::string FCompilerContext::QueryTable(const std::string& Name, const std::string& Property, FTable Table)
{
    for (auto& record : Table.Records) {
        if (record["Name"] == Name) {
            return record[Property];
        }          
    }
//    std::cout << "Record not found in the database. Query: " << Name << " | Property Keyword: " << Property << std::endl;
    return std::string();
}

std::vector<std::string> FCompilerContext::GetNames_MoleculeList() 
{
    std::vector<std::string> StrList;
    for (auto& item : MoleculeList){
        StrList.push_back(item->Name);
    }
    return StrList;
}

const FEnzyme* FCompilerContext::GetEnzyme_EnzymeList(std::string Name) 
{
    for (auto& enzyme : EnzymeList){
        if (enzyme->Name == Name) {
            return enzyme;
        }
    }
}

std::vector<std::string> FCompilerContext::GetNames_EnzymeList(std::vector<const FEnzyme *> EnzymeList) 
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymeList){
        StrList.push_back(item->Name);
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetSubstrateNames_EnzymeList(std::vector<const FEnzyme *> EnzymeList) 
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymeList){
        StrList.push_back(item->Substrate);
    }
    return StrList;
}

std::vector<float> FCompilerContext::Getkcats_EnzymeList(std::vector<const FEnzyme *> EnzymeList)
{
    std::vector<float> floatList;
    for (auto& item : EnzymeList){
        floatList.push_back(item->kcat);
    }
    return floatList;
}

std::vector<float> FCompilerContext::GetKMs_EnzymeList(std::vector<const FEnzyme *> EnzymeList)
{
    std::vector<float> floatList;
    for (auto& item : EnzymeList){
        floatList.push_back(item->KM);
    }
    return floatList;
}

std::vector<float> FCompilerContext::Getks_EnzymeList(std::vector<const FEnzyme *> EnzymeList) 
{
    std::vector<float> floatList;
    for (auto& item : EnzymeList){
        floatList.push_back(item->k);
    }
    return floatList;
}

std::vector<float> FCompilerContext::Getkrevs_EnzymeList(std::vector<const FEnzyme *> EnzymeList) 
{
    std::vector<float> floatList;
    for (auto& item : EnzymeList){
        floatList.push_back(item->krev);
    }
    return floatList;
}

std::vector<const FReaction *> FCompilerContext::GetSubList_ReactionList(std::string Type)
{
    int ReactionType;

    //Standard Reactions
  
    if      (Type == "Standard_Unregulated")             { ReactionType = 0; }
    else if (Type == "Standard_Inhibition_Allosteric")   { ReactionType = 1; }
    else if (Type == "Standard_Activation_Allosteric")   { ReactionType = 2; }

    // Enzymatic Reactions:Standard
    else if (Type == "Enz_Standard_Unregulated")            { ReactionType = 10; }
    else if (Type == "Enz_Standard_Inhibition_Allosteric")  { ReactionType = 11; }
    else if (Type == "Enz_Standard_Activation_Allosteric")  { ReactionType = 12; }

    // Enzymatic Reactions:MichaelisMenten
    else if (Type == "Enz_MichaelisMenten_Unregulated")            { ReactionType = 20; }
    else if (Type == "Enz_MichaelisMenten_Inhibition_Allosteric")  { ReactionType = 21; }
    else if (Type == "Enz_MichaelisMenten_Inhibition_Competitive") { ReactionType = 22; }
    else if (Type == "Enz_MichaelisMenten_Activation_Allosteric")  { ReactionType = 23; }

    // Regulatory Reactions
    else if (Type == "Regulatory_Inhibition_Allosteric")    { ReactionType = 100; }
    else if (Type == "Regulatory_Inhibition_Competitive")   { ReactionType = 101; }
    else if (Type == "Regulatory_Activation_Allosteric")    { ReactionType = 102; }

    else {
//        std::cout << "Unable to find the reaction type for requested Reaction Type: " + Type << std:: endl;
    }

    std::vector<const FReaction *> SubList;

    for (auto& reaction : ReactionList) {
        if (reaction->Type == ReactionType) {
            SubList.push_back(reaction);
        }
    }

    return SubList;
}

std::vector<const FStandardReaction *> FCompilerContext::GetList_Standard_ReactionList(std::string ReactionType)
{
    std::vector<const FStandardReaction *> SubList;
    std::vector<const FReaction *> ReactionSubList = GetSubList_ReactionList(ReactionType);
    

    for (const FReaction* Reaction: ReactionSubList) {
        auto reaction = dynamic_cast<const FStandardReaction *>(Reaction);
        SubList.push_back(reaction);
    }

    return SubList;
}

std::vector<const FRegulatoryReaction *> FCompilerContext::GetList_Regulatory_ReactionList(std::string ReactionType)
{
    // Reaction Type is either "Inhibition" or "Activation"

    std::vector<const FRegulatoryReaction *> SubList;

    for (const FReaction* Reaction: ReactionList) {
        if (Utils::is_class_of<FRegulatoryReaction, FReaction>(Reaction)) {
            auto reaction = dynamic_cast<const FRegulatoryReaction *>(Reaction);
            if (reaction->Effect == ReactionType) {
                SubList.push_back(reaction);
            }
        }
    }

    return SubList;
}

std::vector<const FEnzymaticReaction *> FCompilerContext::GetList_Enzymatic_ReactionList()
{
    std::vector<const FEnzymaticReaction *> SubList;

    for (const FReaction* Reaction: ReactionList) {
        if (Utils::is_class_of<FEnzymaticReaction, FReaction>(Reaction)) {
            auto EnzymaticReaction = dynamic_cast<const FEnzymaticReaction* >(Reaction);
            SubList.push_back(EnzymaticReaction);
        }
    }
    return SubList;
}

std::vector<std::string> FCompilerContext::GetNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *> EnzymaticReactionList)
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

std::vector<std::string> FCompilerContext::GetSubstrateNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *> EnzymaticReactionList)
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

std::vector<std::string> FCompilerContext::GetReactantNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *> EnzymaticReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymaticReactionList){
        for (auto& Stoich : item->Stoichiometry){
            if ((std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) & (Stoich.second <= 0)){
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetProductNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *> EnzymaticReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymaticReactionList){
        for (auto& Stoich : item->Stoichiometry){
            if ((std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) & (Stoich.second >= 0)){
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetEnzymeNames_EnzymaticReactionList(std::vector<const FEnzymaticReaction *> EnzymaticReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymaticReactionList){
        StrList.push_back(item->Enzyme);
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetNames_PolymeraseList(std::vector<const FPolymerase *> PolymeraseList) 
{
    std::vector<std::string> StrList;
    for (auto& item : PolymeraseList){
        StrList.push_back(item->Name);
    }
    return StrList;
}

std::vector<float> FCompilerContext::GetRates_PolymeraseList(std::vector<const FPolymerase *> PolymeraseList)
{
    std::vector<float> floatList;
    for (auto& item : PolymeraseList){
        floatList.push_back(item->Rate);
    }
    return floatList;
}

std::vector<const FPolymeraseReaction *> FCompilerContext::GetList_Polymerase_ReactionList()
{
    std::vector<const FPolymeraseReaction *> SubList;

    for (const FReaction* Reaction: ReactionList) {
        if (Utils::is_class_of<FPolymeraseReaction, FReaction>(Reaction)) {
            auto PolymeraseReaction = dynamic_cast<const FPolymeraseReaction* >(Reaction);
            SubList.push_back(PolymeraseReaction);
        }
    }
    return SubList;
}

std::vector<std::string> FCompilerContext::GetNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *> PolymeraseReactionList)
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

std::vector<std::string> FCompilerContext::GetSubstrateNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *> PolymeraseReactionList)
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

std::vector<std::string> FCompilerContext::GetReactantNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *> PolymeraseReactionList)
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

std::vector<std::string> FCompilerContext::GetProductNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *> PolymeraseReactionList)
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

std::vector<std::string> FCompilerContext::GetBuildingBlockNames_PolymeraseReactionList(std::vector<const FPolymeraseReaction *> PolymeraseReactionList)
{
    std::vector<std::string> StrList;
    for (auto& item : PolymeraseReactionList){
        for (auto& BuildingBlock : item->BuildingBlocks){
            StrList.push_back(BuildingBlock);
        }
    }
    return StrList;
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
                auto reaction = dynamic_cast<const FRegulatoryReaction *>(Reaction);

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
        StrList.push_back(item.Name);
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetSequences_PathwayList()
{
    std::vector<std::string> StrList;
    for (auto& item : PathwayList){
        for (auto& subitem : item.Sequence){
            if (std::find(StrList.begin(), StrList.end(), subitem) == StrList.end()) {
                StrList.push_back(subitem);
            }
        }
    }
    return StrList;
}

std::vector<int> FCompilerContext::GetCoefficientArray(const FReaction* Reaction, std::vector<int> Idx_Substrates)
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

float FCompilerContext::GetFloatAttributeByName_EnzymeList(std::vector<const FEnzyme *> EnzymeList, std::string Name, std::string Attribute)
{
    float Value = -1;

    for (auto& enzyme : EnzymeList){
        if (enzyme->Name == Name) {
            if      (Attribute == "k")      { Value = enzyme->k; }
            else if (Attribute == "krev")   { Value = enzyme->krev; }
            else if (Attribute == "kcat")   { Value = enzyme->kcat; }
            else if (Attribute == "KM")     { Value = enzyme->KM; }
            break;
        }
    }

    if (Value != -1) { return Value; }
    else             { std::cout << "Undefined enzyme name: " << Name << std::endl; }
}

std::string FCompilerContext::GetStringAttributeByName_EnzymeList(std::vector<const FEnzyme *> EnzymeList, std::string Name, std::string Attribute)
{
    std::string Value;

    for (auto& enzyme : EnzymeList){
        if (enzyme->Name == Name) {
            if      (Attribute == "Substrate")      { Value = enzyme->Substrate; } 
            // TODO: add more after updating FEnzyme class
            break;
        }
    }

    if (Value.empty()) { 
        std::cout << "Undefined enzyme name or attribute: " << Name << ", " << Attribute << std::endl; 
    }
    assert (!Value.empty());
    return Value;
}

std::vector<std::vector<int>> FCompilerContext::GetStoichiometryMatrix(std::string Type)
{
    // Types: "Standard_Unregulated", "Standard_Inhibition", "Standard_Activation",
    //        "Enz_Standard_Unregulated", "Enz_Standard_Inhibition", "Enz_Standard_Activation",
    //        "Enz_MichaelisMenten_Unregulated", "Enz_MichaelisMenten_Inhibition_Allosteric", "Enz_MichaelisMenten_Inhibition_Competitive", "Enz_MichaelisMenten_Inhibition_Competitive"

    std::vector<std::vector<int>> StoichMatrix;

    // for matrix generation
    std::vector<int> Idx_Substrates = GetIdxForStoichiometryMatrix(Type);
    std::vector<const FReaction *> ReactionList = GetSubList_ReactionList(Type);

    for (auto& reaction : ReactionList) {
        StoichMatrix.push_back(GetCoefficientArray(reaction, Idx_Substrates));
    }

    return StoichMatrix;
}

std::vector<int> FCompilerContext::AddUniqueSubstrateIdxToIdxList(const FReaction* Reaction, std::vector<int>IdxList)
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

std::vector<int> FCompilerContext::GetIdxForStoichiometryMatrix(std::string Type)
{
    // Types: "Standard_Unregulated", "Standard_Inhibition", "Standard_Activation",
    //        "Enz_Standard_Unregulated", "Enz_Standard_Inhibition", "Enz_Standard_Activation",
    //        "Enz_MichaelisMenten_Unregulated", "Enz_MichaelisMenten_Inhibition_Allosteric", "Enz_MichaelisMenten_Inhibition_Competitive", "Enz_MichaelisMenten_Inhibition_Competitive"

    std::vector<int> IdxList;
    std::vector<const FReaction *> ReactionList = GetSubList_ReactionList(Type);

    for (auto& reaction : ReactionList) {
        IdxList = AddUniqueSubstrateIdxToIdxList(reaction, IdxList);
    }

    return IdxList;
}

std::vector<std::vector<int>> FCompilerContext::GetStoichiometryMatrix_PolymeraseReaction(std::vector<const FPolymeraseReaction *> PolymeraseReactionList)
{
    std::vector<std::vector<int>> StoichMatrix;
    std::vector<const FSmallMolecule *> SMolList = GetList_SmallMolecule_MoleculeList();

    for (auto& PolymeraseReaction : PolymeraseReactionList){
        std::vector<int> CoeffArray(SMolList.size(), 0); // replace with substrate index for the reaction

        for (auto& Stoich : PolymeraseReaction->Stoichiometry){
            std::string SubstrateName = Stoich.first;
            int Coeff = Stoich.second;
	    int Index = 0;

            for (auto& molecule : SMolList){
                // std::cout << "Searching from the List: " << Substrate << " | " << "Index: " << Index << endl;
                if (molecule->Name == SubstrateName){
                    break;
                }
                Index++;
            }
            if (Index >= MoleculeList.size()) {
	std::cout << "Substrate index searching in MoleculeList failed: " << SubstrateName << endl;
            } else {
            // std::cout << "Substrate Searching from the List: " << SubstrateName << " | " << "Index: " << Index << endl;
            } 

            CoeffArray[Index] = Coeff;
        }
        StoichMatrix.push_back(CoeffArray);         
    }
    return StoichMatrix;
}

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
    std::cout << "Unable to find index in MoleculeList : " << InputName << std::endl;
}

float FCompilerContext::GetInitialCountByName_CountList(std::string MolName)
{
// returns 0 if there is no initial count set
    for (auto& count : CountList) {
        if ((count->Name == MolName) & (count->Begin == 0)) {
            return count->Amount;
            break;
        }
    }
    return 0;
}

std::vector<const FGene *> FCompilerContext::GetList_Gene_MoleculeList()
{
    std::vector<const FGene *> SubList;
    
    for (auto* molecule :MoleculeList) {
        if (Utils::is_class_of<FGene, FMolecule>(molecule)) {
            auto Item = dynamic_cast<const FGene *>(molecule);
            SubList.push_back(Item);
        }
    }
    return SubList;
}

std::vector<const FProtein *> FCompilerContext::GetList_Protein_MoleculeList()
{
    std::vector<const FProtein *> SubList;
    
    for (auto* molecule :MoleculeList) {
        if (Utils::is_class_of<FProtein, FMolecule>(molecule)) {
            auto Item = dynamic_cast<const FProtein *>(molecule);
            SubList.push_back(Item);
        }
    }
    return SubList;
}

std::vector<const FEnzyme *> FCompilerContext::GetList_Enzyme_MoleculeList()
{
    std::vector<const FEnzyme *> SubList;
    
    for (const FMolecule* molecule :MoleculeList) {
        if (Utils::is_class_of<FEnzyme, FMolecule>(molecule)) {
            auto Enzyme = dynamic_cast<const FEnzyme *>(molecule);
            SubList.push_back(Enzyme);
        }
    }
    return SubList;
}

std::vector<const FSmallMolecule *> FCompilerContext::GetList_SmallMolecule_MoleculeList()
{
    std::vector<const FSmallMolecule *> SubList;
    
    for (const FMolecule* molecule :MoleculeList) {
        if (Utils::is_class_of<FSmallMolecule, FMolecule>(molecule)) {
            auto SmallMolecule = dynamic_cast<const FSmallMolecule *>(molecule);
            SubList.push_back(SmallMolecule);
        }
    }
    return SubList;
}

std::vector<const FPolymerase *> FCompilerContext::GetList_Polymerase_MoleculeList()
{
    std::vector<const FPolymerase *> SubList;
    
    for (const FMolecule* molecule :MoleculeList) {
        if (Utils::is_class_of<FPolymerase, FMolecule>(molecule)) {
            auto Polymerase = dynamic_cast<const FPolymerase *>(molecule);
            SubList.push_back(Polymerase);
        }
    }
    return SubList;
}

std::vector<int> FCompilerContext::GetIdxListFromMoleculeList(std::string FClassName)
{ 
    std::vector<int> IndexArray;
    int Index = 0;

//    // update this code after studying template syntax
//    template <typename Derived, typename Base> 
//    static bool is_class_of(const Base*Node) {
//        Derived* DerivedNode = dynamic_cast<Derived *>(const_cast<Base *>(Node));
//        return DerivedNode != nullptr;
//    }

    if (FClassName == "Chromosome") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FChromosome, FMolecule>(molecule)) {
                IndexArray.push_back(Index);
            }
            Index++;
        }
    } else if (FClassName == "Gene") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FGene, FMolecule>(molecule)) {
                IndexArray.push_back(Index);
            }
            Index++;
        }
    } else if (FClassName == "RNA") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FRNA, FMolecule>(molecule)) {
                IndexArray.push_back(Index);
            }
            Index++;
        }
    } else if (FClassName == "Protein") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FProtein, FMolecule>(molecule)) {
                IndexArray.push_back(Index);
            }
            Index++;
        }
    } else if (FClassName == "Enzyme") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FEnzyme, FMolecule>(molecule)) {
                IndexArray.push_back(Index);
            }
            Index++;
        }
    } else if (FClassName == "Molecule") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FMolecule, FMolecule>(molecule)) {
                IndexArray.push_back(Index);
            }
            Index++;
        }
    } else if (FClassName == "SmallMolecule") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FSmallMolecule, FMolecule>(molecule)) {
                IndexArray.push_back(Index);
            }
            Index++;
        }
    } else if (FClassName == "Polymerase") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FPolymerase, FMolecule>(molecule)) {
                IndexArray.push_back(Index);
            }
            Index++;
        }
    } else if (FClassName == "mRNA") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FRNA, FMolecule>(molecule)) {
                FRNA * RNA = dynamic_cast<FRNA *>(molecule);
                if (RNA->Type == "mRNA") {
                    IndexArray.push_back(Index);
                }
            }
            Index++;
        }
    }
    return IndexArray;
}

std::vector<std::string> FCompilerContext::GetNameListFromMoleculeList(std::string FClassName)
{ 
    std::vector<string> StrArray;

//    // update this code after studying template syntax
//    template <typename Derived, typename Base> 
//    static bool is_class_of(const Base*Node) {
//        Derived* DerivedNode = dynamic_cast<Derived *>(const_cast<Base *>(Node));
//        return DerivedNode != nullptr;
//    }

    if (FClassName == "Chromosome") {
        for (FMolecule * molecule : MoleculeList) {
            if (Utils::is_class_of<FChromosome, FMolecule>(molecule)) {
                StrArray.push_back(molecule->Name);
            }
        }
    } else if (FClassName == "Gene") {
        for (FMolecule * molecule : MoleculeList) {
            if (Utils::is_class_of<FGene, FMolecule>(molecule)) {
                StrArray.push_back(molecule->Name);
            }
        }
    } else if (FClassName == "RNA") {
        for (FMolecule * molecule : MoleculeList) {
            if (Utils::is_class_of<FRNA, FMolecule>(molecule)) {
                StrArray.push_back(molecule->Name);
            }
        }
    } else if (FClassName == "Protein") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FProtein, FMolecule>(molecule)) {
                StrArray.push_back(molecule->Name);
            }
        }
    } else if (FClassName == "Enzyme") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FEnzyme, FMolecule>(molecule)) {
                StrArray.push_back(molecule->Name);
            }
        }
    } else if (FClassName == "SmallMolecule") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FSmallMolecule, FMolecule>(molecule)) {
                StrArray.push_back(molecule->Name);
            }
        }
    } else if (FClassName == "Polymerase") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FPolymerase, FMolecule>(molecule)) {
                StrArray.push_back(molecule->Name);
            }
        }
    } else if (FClassName == "mRNA") {
        for (FMolecule * molecule :MoleculeList) {
            if (Utils::is_class_of<FRNA, FMolecule>(molecule)) {
                FRNA * RNA = dynamic_cast<FRNA *>(molecule);
                if (RNA->Type == "mRNA") {
                    StrArray.push_back(molecule->Name);
                }
            }
        }
    }
    return StrArray;
}

std::vector<int> FCompilerContext::GetIdx_EnzymeSubstrate_MoleculeList()
{ 
    std::vector<int> IndexArray;
    int Index;

    for (const FMolecule* molecule :MoleculeList) {
        if (Utils::is_class_of<FEnzyme, FMolecule>(molecule)) {
            auto enzyme = dynamic_cast<const FEnzyme *>(molecule);
            std::string EnzSub = enzyme->Substrate;
            Index = 0;
            for (auto& molecule : MoleculeList) {
                if (molecule->Name == EnzSub) {
                    IndexArray.push_back(Index);
                break;
                } 
                Index++;
            }
        }
    }
    return IndexArray;
}

// std::vector<int> FCompilerContext::GetIdx_PolymeraseSubstrate_MoleculeList()
// { 
//     std::vector<int> IndexArray;
//     int Index;
// 
//     for (const FMolecule* molecule :MoleculeList) {
//         if (Utils::is_class_of<FPolymerase, FMolecule>(molecule)) {
//             auto Polymerase = dynamic_cast<const FPolymerase *>(molecule);
//             std::string PolSub = Polymerase->Substrate;
//             Index = 0;
//             for (auto& molecule :MoleculeList) {
//                 if (molecule->Name == PolSub) {
//                     IndexArray.push_back(Index);
//                 break;
//                 } 
//                 Index++;
//             }
//         }
//     }
//     return IndexArray;
// }

std::vector<int> FCompilerContext::GetIdx_PolymeraseReactionSubstrate_ByPolymeraseName_MoleculeList(std::string InPolymeraseName)
{ 
    std::vector<int> IndexArray;
    int Index;

    std::vector<const FPolymeraseReaction *> PolymeraseReactionList = GetList_Polymerase_ReactionList();

    for (auto& PolymeraseReaction : PolymeraseReactionList){
        if (PolymeraseReaction->Polymerase == InPolymeraseName){
            for (auto& stoich : PolymeraseReaction->Stoichiometry) {
                Index = GetIdxByName_MoleculeList(stoich.first);
                IndexArray.push_back(Index);
            }
        }
    }
    return IndexArray;
}

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

//
//std::vector<float> FCompilerContext::GetInitialCountByStrList_MoleculeList(std::vector<std::string> StrList)
//{
//    std::vector<int> IndexArray;
//    int Index;
//
//    for (std::string Item : StrList){
//        Index = GetInitialCountByName_MoleculeList(Item);
//        IndexArray.push_back(Index);
//    }
//    return IndexArray;
//}

std::vector<int> FCompilerContext::GetIdxOfStrListFromStrList(std::vector<std::string> InputList, std::vector<std::string> RefList)
{
    std::vector<int> IndexArray;
    int Index;

    for (auto& InputItem : InputList) {
        Index = 0;
        for (auto& RefItem : RefList) {
            if (RefItem == InputItem) {
                break;
            } else {
                Index++;
            }
         }
        IndexArray.push_back(Index);
    }
    return IndexArray;
}

// std::vector<int> FCompilerContext::GetEnzSubstrateIdxFromAllSubstrates(std::vector<const FEnzyme *> EnzymeList)
// {
//     std::vector<std::string> EnzSubstrateList = GetSubstrateNames_EnzymeList(EnzymeList);
//     std::vector<std::string> AllSubstrateList = GetSubstrateNames_EnzymaticReactionList();
// 
//     return GetIdxListFromList(EnzSubstrateList, AllSubstrateList);
// }

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

