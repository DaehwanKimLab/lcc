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

float FCount::GetCount(float Time) {
    for (auto& range: Range) {
        if ((Time >= range.first.first) & (Time < range.first.second)) {
            return range.second;
        }
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

void FCompilerContext::PrintLists(std::ostream& os) 
{
    os << std::endl << "## Compiler Context Lists ##" << std::endl;
    if (!MoleculeList.empty()) {
        os << "  MoleculeList: " << std::endl << "  " << "  ";
        for (auto& item : MoleculeList){
            os << item->Name << ", ";
        }
        os << std::endl;
    }

    if (!PathwayList.empty()) { 
        os << "  PathwayList: " << std::endl << "  " << "  ";
        for (auto& item : PathwayList){
            os << item.Name << ", ";
        }
        os << std::endl;
    }
//    os << "  EnzymeList: " << std::endl << "  " << "  ";
//    for (auto& item : EnzymeList){
//        os << item.Name << ", ";
//    };
//    os << std::endl;
    if (!ReactionList.empty()) {
        os << "  ReactionList: " << std::endl << "  " << "  ";
        for (auto& item : ReactionList){
            os << item->Name << ", ";
        }
        os << std::endl;
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

std::vector<std::string> FCompilerContext::GetInhibitorNames_EnzymeList(std::vector<const FEnzyme *> EnzymeList) 
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymeList){
        StrList.push_back(item->Inhibitor);
    }
    return StrList;
}

std::vector<float> FCompilerContext::GetKis_EnzymeList(std::vector<const FEnzyme *> EnzymeList)
{
    std::vector<float> floatList;
    for (auto& item : EnzymeList){
        floatList.push_back(item->Ki);
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

std::vector<const FEnzyme *> FCompilerContext::GetSubList_EnzymeList(std::vector<const FEnzyme*> EnzymeList, std::string Type)
{
    std::vector<const FEnzyme *> SubList;
    
    for (auto& enzyme : EnzymeList){
        if (Type == "Enz_Standard") {
            if ((enzyme->k >= 0) & (enzyme->KM < 0) & (enzyme->Mode.empty())) {
                SubList.push_back(enzyme);
            }
        } else if (Type == "Enz_Standard_Inhibition_Allosteric") {
            if ((enzyme->k >= 0) & (enzyme->KM < 0) & (!enzyme->Inhibitor.empty()) & (enzyme->Activator.empty()) & (enzyme->Mode == "Allosteric")) {
                SubList.push_back(enzyme);
            }
        } else if (Type == "Enz_Standard_Activation_Allosteric") {
            if ((enzyme->k >= 0) & (enzyme->KM < 0) & (enzyme->Inhibitor.empty()) & (!enzyme->Activator.empty()) & (enzyme->Mode == "Allosteric")) {
                SubList.push_back(enzyme);
            }
        } else if (Type == "Enz_MichaelisMenten") {
            if ((enzyme->KM >= 0) & (enzyme->k < 0) & (enzyme->Mode.empty())) {
                SubList.push_back(enzyme);
            }
        } else if (Type == "Enz_MichaelisMenten_Inhibition_Competitive") {
            if ((enzyme->KM >= 0) & (enzyme->k < 0) & (!enzyme->Inhibitor.empty()) & (enzyme->Activator.empty()) & (enzyme->Mode == "Competitive")) {
                SubList.push_back(enzyme);
            }
        } else if (Type == "Enz_MichaelisMenten_Inhibition_Allosteric") {
            if ((enzyme->KM >= 0) & (enzyme->k < 0) & (!enzyme->Inhibitor.empty()) & (enzyme->Activator.empty()) & (enzyme->Mode == "Allosteric")) {
                SubList.push_back(enzyme);
            }
        } else if (Type == "Enz_MichaelisMenten_Activation_Allosteric") {
            if ((enzyme->KM >= 0) & (enzyme->k < 0) & (enzyme->Inhibitor.empty()) & (!enzyme->Activator.empty()) & (enzyme->Mode == "Allosteric")) {
                SubList.push_back(enzyme);
            }
        }
    }
    return SubList;
}

std::vector<const FStandardReaction *> FCompilerContext::GetList_Standard_ReactionList(std::string ReactionType)
{
    std::vector<const FStandardReaction *> SubList;

    std::vector<std::string> Inhibited = GetNames_ReactionList("Inhibited");
    std::vector<std::string> Activated = GetNames_ReactionList("Activated");

    for (const FReaction* Reaction: ReactionList) {
        if (Utils::is_class_of<FStandardReaction, FReaction>(Reaction)) {
            auto reaction = dynamic_cast<const FStandardReaction *>(Reaction);
            if (ReactionType == "Standard_All") {
                SubList.push_back(reaction);
            } else if (ReactionType == "Standard_Unregulated") {       
                if ((std::find(Inhibited.begin(), Inhibited.end(), reaction->Name) == Inhibited.end()) 
                  & (std::find(Activated.begin(), Activated.end(), reaction->Name) == Activated.end())) {
                    SubList.push_back(reaction);
                }
            } else if (ReactionType == "Standard_Inhibited") {
                if (std::find(Inhibited.begin(), Inhibited.end(), reaction->Name) != Inhibited.end()) {
                    SubList.push_back(reaction);
                }
            } else if (ReactionType == "Standard_Activated") {
                if (std::find(Activated.begin(), Activated.end(), reaction->Name) != Activated.end()) {
                    SubList.push_back(reaction);
                }
            }
        }
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
   
    for (auto& Reaction : ReactionList){
        if (Utils::is_class_of<FRegulatoryReaction, FReaction>(Reaction)) {
            auto reaction = dynamic_cast<const FRegulatoryReaction *>(Reaction);
            if (Type == "Inhibited") {
                if (reaction->Effect == "Inhibition") {
                    for (auto& stoich : reaction->Stoichiometry) {
                        if (stoich.second > 0) {
                            StrList.push_back(stoich.first);
                        }
                    }
                }
            } else if (Type == "Activated") {
                if (reaction->Effect == "Activation") {
                    for (auto& stoich : reaction->Stoichiometry) {
                        if (stoich.second > 0) {
                            StrList.push_back(stoich.first);
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

std::vector<std::vector<int>> FCompilerContext::GetStoichiometryMatrix(std::string Type)
{
    // Types: "Enz_MichaelisMenten", "Enz_Standard"

    std::vector<std::vector<int>> StoichMatrix;
    std::vector<int> Idx_Substrates = GetIdxForStoichiometryMatrix(Type);
    std::vector<const FEnzyme *> EnzymeList = GetList_Enzyme_MoleculeList();
    std::vector<const FSmallMolecule *> SMolList = GetList_SmallMolecule_MoleculeList();
    std::vector<const FEnzymaticReaction *> EnzymaticReactionList = GetList_Enzymatic_ReactionList();

    std::vector<const FEnzyme *> EnzymeList_Sub = GetSubList_EnzymeList(EnzymeList, Type);
 
    for (auto& enzyme : EnzymeList_Sub){
        for (auto& enzymaticReaction : EnzymaticReactionList){
            if (enzymaticReaction->Enzyme == enzyme->Name) {
                std::vector<int> CoeffArray(Idx_Substrates.size(), 0);
                for (auto& stoich : enzymaticReaction->Stoichiometry) {
				// std::cout << "enzymaticReaction-> Stoichiometry for loop. Now working on: " << stoich.first << endl;
                    std::string SubstrateName = stoich.first;
                    int Coeff = stoich.second;
                    int MolIdx = GetIdxByName_MoleculeList(SubstrateName);
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
				// std::cout << "pushback" << endl;
                StoichMatrix.push_back(CoeffArray);
            }
        }
    } 
    return StoichMatrix;
}

std::vector<int> FCompilerContext::GetIdxForStoichiometryMatrix(std::string Type)
{
    // Types: "Enz_MichaelisMenten", "Enz_Standard"
    std::vector<int> IdxList;

    std::vector<const FEnzyme *> EnzymeList = GetList_Enzyme_MoleculeList();
    std::vector<const FEnzymaticReaction *> EnzymaticReactionList = GetList_Enzymatic_ReactionList();
    // std::vector<const FSmallMolecule *> SMolList = GetList_SmallMolecule_MoleculeList();

    std::vector<const FEnzyme *> EnzymeList_Sub = GetSubList_EnzymeList(EnzymeList, Type);
 
    for (auto& enzyme : EnzymeList_Sub){
        for (auto& enzymaticReaction : EnzymaticReactionList) {
            if (enzymaticReaction->Enzyme == enzyme->Name) {
                for (auto& stoich : enzymaticReaction->Stoichiometry) {
                    int MolIdx = 0;
                    for (auto& mol : MoleculeList) {
                        if (mol->Name == stoich.first) {
                            if (std::find(IdxList.begin(), IdxList.end(), MolIdx) == IdxList.end()) {
//std::cout << "New mol idx added for : " << stoich.first << endl;
                                IdxList.push_back(MolIdx);
                                break;
                            } else {
//std::cout << "Redundant mol idx already exists for: " << stoich.first << endl;
                            }
                        }
                        MolIdx++;
                    //int MolIdx = GetIdxByName_MoleculeList(stoich.first);
                    //if (std::count(IdxList.begin(), IdxList.end(), MolIdx)) {
                    //    break;
                    //}
                    //IdxList.push_back(MolIdx);                        
                    }
                }
            }
        }
    }
    return IdxList;
}

std::vector<std::vector<int>> FCompilerContext::GetStoichiometryMatrix_PolymeraseReaction(std::vector<const FPolymeraseReaction *> PolymeraseReactionList)
{
    std::vector<std::vector<int>> StoichMatrix;
    std::vector<const FSmallMolecule *> SMolList = GetList_SmallMolecule_MoleculeList();

    for (auto& PolymeraseReaction : PolymeraseReactionList){
        std::vector<int> CoeffArray(SMolList.size(), 0);

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
    int Index = 0;
    for (auto& molecule : MoleculeList) {
        if (molecule->Name == InputName) {
            break;
        }
        Index++; 
    }
    return Index;
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

