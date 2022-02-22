#include <iostream>
#include "state.h"

using namespace std;

void FState::SetVol(const int& InVol){
    Vol = InVol;
}

void FState::SetEnzCount(const std::string& Name, const int& Count){
    Enz2Count[Name] = Count;
}

void FState::SetSubCount(const std::string& Name, const int& Count){
    Sub2Count[Name] = Count;
}

int FState::GetEnzCount(const std::string& Name) {
    return Enz2Count[Name];
}

int FState::GetSubCount(const std::string& Name) {
    return Sub2Count[Name];
}

float FState::GetEnzConc(const std::string& Name) {
    return Enz2Count[Name] / Vol;
}

float FState::GetSubConc(const std::string& Name) {
    return Sub2Count[Name] / Vol;
}

void FState::PrintState() {
    std::cout << "EnzCounts:" << std::endl;
    for (auto& Enz : Enz2Count) {
       std::cout<< "  " << "PRINT HERE" << std::endl; 
    }
    
    std::cout << "SubCounts:" << std::endl;
    for (auto& Sub : Sub2Count) {
       std::cout<< "  " << "PRINT HERE" << std::endl; 
    }
}

std::vector<std::string> FState::ExportLegend() {
    std::vector<std::string> Legend;
    Legend.resize(Enz2Count.size() + Sub2Count.size() + 2); // 2 for step and vol

    int i = 0;
    
    Legend[0] = "SimStep"; i++;
    Legend[1] = "Vol"; i++;

    std::vector<float> EnzCounts;
    for (auto& KeyValue : Enz2Count){
        Legend[i] = KeyValue.first;
        i++;
    }

    std::vector<float> SubCounts;
    for (auto& KeyValue : Sub2Count){
        Legend[i] = KeyValue.first;
        i++;
    }

    return Legend;
}

std::vector<float> FState::ExportState() {
    std::vector<float> DataExport;
    DataExport.resize(Enz2Count.size() + Sub2Count.size() + 2); // 2 for step and vol

    int i = 0;
    
    DataExport[0] = SimStep; i++;
    DataExport[1] = Vol; i++;

    std::vector<float> EnzCounts;
    for (auto& KeyValue : Enz2Count){
        DataExport[i] = KeyValue.second;
        i++;
    }

    std::vector<float> SubCounts;
    for (auto& KeyValue : Sub2Count){
        DataExport[i] = KeyValue.second;
        i++;
    }

    SimStep++;

    return DataExport;
}

