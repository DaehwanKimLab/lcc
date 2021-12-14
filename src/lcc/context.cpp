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

void FCompilerContext::PrintLists(std::ostream& os) 
{
    os << "Compiler Context Lists:" << std::endl;

    os << "  PathwayList: " << std::endl << "  " << "  ";
    for (auto& item : PathwayList){
        os << item.Name << ", ";
    };
    os << std::endl;

    os << "  EnzymeList: " << std::endl << "  " << "  ";
    for (auto& item : EnzymeList){
        os << item.Name << ", ";
    };
    os << std::endl;

    os << "  EnzymaticReactionList: " << std::endl << "  " << "  ";
    for (auto& item : EnzymaticReactionList){
        os << item.Name << ", ";
    };
    os << std::endl;
}

const std::string FCompilerContext::QueryEnzymeTable(const std::string& Name, const std::string& Property)
{
    for (auto& record : EnzymeTable.Records) {
        if (Name == record["EnzymeName"]) {
            return record[Property];
        }          
    }
    std::cout << "No such enzyme found in the database: " << Name << std::endl;
}

const std::string FCompilerContext::QueryReactionTable(const std::string& Name, const std::string& Property)
{
    for (auto& record : ReactionTable.Records) {
        if (Name == record["EnzymeName"]) {
            return record[Property];
        }          
    }
    std::cout << "No such reaction found in the database: " << Name << std::endl;
}

std::vector<std::string> FCompilerContext::GetNames_EnzymeList() 
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymeList){
        StrList.push_back(item.Name);
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetSubstrates_EnzymeList() 
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymeList){
        StrList.push_back(item.Substrate);
    }
    return StrList;
}

std::vector<float> FCompilerContext::Getkcats_EnzymeList()
{
    std::vector<float> FloatList;
    for (auto& item : EnzymeList){
        FloatList.push_back(item.kcat);
    }
    return FloatList;
}

std::vector<float> FCompilerContext::GetkMs_EnzymeList() 
{
    std::vector<float> FloatList;
    for (auto& item : EnzymeList){
        FloatList.push_back(item.kM);
    }
    return FloatList;
}

std::vector<std::string> FCompilerContext::GetNames_EnzymaticReactionList()
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymaticReactionList){
        StrList.push_back(item.Name);
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetSubstrateNames_EnzymaticReactionList()
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymaticReactionList){
        for (auto& Stoich : item.Stoichiometry){
            if (std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) {
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetReactantNames_EnzymaticReactionList()
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymaticReactionList){
        for (auto& Stoich : item.Stoichiometry){
            if ((std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) & (Stoich.second < 0)){
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetProductNames_EnzymaticReactionList()
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymaticReactionList){
        for (auto& Stoich : item.Stoichiometry){
            if ((std::find(StrList.begin(), StrList.end(), Stoich.first) == StrList.end()) & (Stoich.second > 0)){
                StrList.push_back(Stoich.first);
            }
        }
    }
    return StrList;
}

std::vector<std::string> FCompilerContext::GetEnzymeNames_EnzymaticReactionList()
{
    std::vector<std::string> StrList;
    for (auto& item : EnzymaticReactionList){
        StrList.push_back(item.Enzyme);
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

std::vector<std::vector<int>> FCompilerContext::GetStoichiometryMatrix()
{
    std::vector<std::vector<int>> StoichMatrix;
    std::vector<std::string> SubstrateList = GetSubstrateNames_EnzymaticReactionList();

    for (auto& EnzymaticReaction : EnzymaticReactionList){
        std::vector<int> CoeffArray(SubstrateList.size(), 0);

        for (auto& Stoich : EnzymaticReaction.Stoichiometry){
            std::string SubstrateName = Stoich.first;
            int Coeff = Stoich.second;
	    int Index = 0;

            std::cout << "NOW SEARCHING: " << SubstrateName << ", " << Coeff << " " << Index << endl;

            for (auto& Substrate : SubstrateList){
                std::cout << "Searching from the List: " << Substrate << " | " << "Index: " << Index << endl;
                if (Substrate == SubstrateName){
                    break;
                }
                Index++;
            }

            if (Index <= SubstrateList.size()) {
                std::cout << "Substrate indexing failed: " << SubstrateName << std::endl;
            } 

            CoeffArray[Index] = Coeff;
        }
        StoichMatrix.push_back(CoeffArray);         
    }
    return StoichMatrix;
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

std::vector<int> FCompilerContext::GetEnzSubstrateIdxFromAllSubstrates()
{
    std::vector<std::string> EnzSubstrateList = GetSubstrates_EnzymeList();
    std::vector<std::string> AllSubstrateList = GetSubstrateNames_EnzymaticReactionList();

    return GetIdxListFromList(EnzSubstrateList, AllSubstrateList);
}



