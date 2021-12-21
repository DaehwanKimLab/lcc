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
        GeneTable.LoadFromTSV((InOption.DataPaths[0] + "/genes.tsv").c_str());
		ReactionTable.LoadFromTSV((InOption.DataPaths[0] + "/reactions.tsv").c_str());
        EnzymeTable.LoadFromTSV((InOption.DataPaths[0] + "/EnzymeDatabase.txt").c_str());
        PolymeraseTable.LoadFromTSV((InOption.DataPaths[0] + "/PolymeraseDatabase.txt").c_str());
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
    for (auto& Molecule : MoleculeList) {
        if (Molecule->Name == NewMolecule->Name){
            Addition = false;
            std::cout << "Redundant molecule " << NewMolecule->Name << " Found in MoleculeList" << std::endl;
            break;
        }
    }
    if (Addition) {
        MoleculeList.push_back(NewMolecule);
    }
}

void FCompilerContext::AddToReactionList(FReaction *NewReaction)
{
    bool Addition = true;
    // std::cout<< "Checking if " << NewMolecule->Name << "Exists in MoleculeList" << std::endl;
    for (auto& Reaction : ReactionList) {
        if (Reaction->Name == Reaction->Name){
            Addition = false;
            // std::cout << "Redundant molecule " << NewMolecule->Name << " Found in MoleculeList" << std::endl;
            break;
        }
    }
    if (Addition) {
        ReactionList.push_back(NewReaction);
    }
}

void FCompilerContext::PrintLists(std::ostream& os) 
{
    os << "## Compiler Context Lists ##" << std::endl;
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

const std::string FCompilerContext::QueryTable(const std::string& Name, const std::string& Property, FTable Table)
{
    for (auto& record : Table.Records) {
        if (record["Name"] == Name) {
            return record[Property];
        }          
    }
    std::cout << "Record not found in the database. Query: " << Name << std::endl;
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
    std::vector<float> FloatList;
    for (auto& item : EnzymeList){
        FloatList.push_back(item->kcat);
    }
    return FloatList;
}

std::vector<float> FCompilerContext::GetkMs_EnzymeList(std::vector<const FEnzyme *> EnzymeList) 
{
    std::vector<float> FloatList;
    for (auto& item : EnzymeList){
        FloatList.push_back(item->kM);
    }
    return FloatList;
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
            if ((std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) & (Stoich.second < 0)){
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
            if ((std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) & (Stoich.second > 0)){
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

std::vector<std::vector<int>> FCompilerContext::GetStoichiometryMatrix_EnzymaticReaction(std::vector<const FEnzymaticReaction *> EnzymaticReactionList)
{
    std::vector<std::vector<int>> StoichMatrix;
    std::vector<const FSmallMolecule *> SMolList = GetList_SmallMolecule_MoleculeList();

    for (auto& EnzymaticReaction : EnzymaticReactionList){
        std::vector<int> CoeffArray(SMolList.size(), 0);

        for (auto& Stoich : EnzymaticReaction->Stoichiometry){
            std::string SubstrateName = Stoich.first;
            int Coeff = Stoich.second;
	    int Index = 0;

            for (auto& Molecule : SMolList){
                // std::cout << "Searching from the List: " << Substrate << " | " << "Index: " << Index << endl;
                if (Molecule->Name == SubstrateName){
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
    for (auto& Molecule : MoleculeList) {
        if (Molecule->Name == InputName) {
            break;
        }
        Index++; 
    }
    return Index;
}

std::vector<const FEnzyme *> FCompilerContext::GetList_Enzyme_MoleculeList()
{
    std::vector<const FEnzyme *> EnzymeList;
    
    for (const FMolecule* Molecule : MoleculeList) {
        if (Utils::is_class_of<FEnzyme, FMolecule>(Molecule)) {
            auto Enzyme = dynamic_cast<const FEnzyme *>(Molecule);
            EnzymeList.push_back(Enzyme);
        }
    }
    return EnzymeList;
}

std::vector<const FSmallMolecule *> FCompilerContext::GetList_SmallMolecule_MoleculeList()
{
    std::vector<const FSmallMolecule *> SmallMoleculeList;
    
    for (const FMolecule* Molecule : MoleculeList) {
        if (Utils::is_class_of<FSmallMolecule, FMolecule>(Molecule)) {
            auto SmallMolecule = dynamic_cast<const FSmallMolecule *>(Molecule);
            SmallMoleculeList.push_back(SmallMolecule);
        }
    }
    return SmallMoleculeList;
}

std::vector<int> FCompilerContext::GetIdx_Enzyme_MoleculeList()
{ 
    std::vector<int> IndexArray;
    int Index;

    Index = 0;
    for (const FMolecule* Molecule : MoleculeList) {
        if (Utils::is_class_of<FEnzyme, FMolecule>(Molecule)) {
            IndexArray.push_back(Index);
        }
        Index++;
    }
    return IndexArray;
}

std::vector<int> FCompilerContext::GetIdx_EnzymeSubstrate_MoleculeList()
{ 
    std::vector<int> IndexArray;
    int Index;

    for (const FMolecule* Molecule : MoleculeList) {
        if (Utils::is_class_of<FEnzyme, FMolecule>(Molecule)) {
            auto Enzyme = dynamic_cast<const FEnzyme *>(Molecule);
            std::string EnzSub = Enzyme->Substrate;
            Index = 0;
            for (auto& Molecule : MoleculeList) {
                if (Molecule->Name == EnzSub) {
                    IndexArray.push_back(Index);
                break;
                } 
                Index++;
            }
        }
    }
    return IndexArray;
}

std::vector<int> FCompilerContext::GetIdx_SmallMolecule_MoleculeList()
{
    std::vector<int> IndexArray;
    int Index;
    
    for (const FMolecule* Molecule : MoleculeList) {
        if (Utils::is_class_of<FSmallMolecule, FMolecule>(Molecule)) {
            IndexArray.push_back(Index);
        }
        Index++;
    }
    return IndexArray;
}

std::vector<int> FCompilerContext::GetIdxListFromList(std::vector<std::string> InputList, std::vector<std::string> RefList)
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

