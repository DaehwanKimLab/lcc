//#ifndef LCC_SIMULATION_H
#define LCC_SIMULATION_H

#include <string>
#include <vector>
#include <map>

// #include "option.h"
// #include "node.h"
#include "context.h"

float ConcentrationToCount(float Conc_Molecule, float Volume) {
    return Conc_Molecule * Volume;
};

float CountToConcentration(int Count_Molecule, float Volume) {
    return Count_Molecule / Volume;
};
 
float MichaelisMentenEqn(int Count_Enzyme, int Count_Substrate, float Volume ,float kcat, float kM) {
    float Conc_Enzyme = CountToConcentration(Count_Enzyme, Volume);
    float Conc_Substrate = CountToConcentration(Count_Substrate, Volume);

    float Rate = (kcat * Conc_Enzyme * Conc_Substrate) / (Conc_Substrate + kM);
    std::cout << Rate << std::endl;
    return Rate;
};

class FState {
public:

// private:
    std::map<std::string, int> Enz2Count;
    std::map<std::string, int> Sub2Count;
    float Vol;

    FState(){};

    void SetVol(const int& InVol){
        Vol = InVol;
    }

    void SetEnzCount(const std::string& Name, const int& Count){
        Enz2Count[Name] = Count;
    }

    void SetSubCount(const std::string& Name, const int& Count){
        Sub2Count[Name] = Count;
    }

    int GetEnzCount(const std::string& Name) {
        return Enz2Count[Name];
    }

    int GetSubCount(const std::string& Name) {
        return Enz2Count[Name];
    }

    float GetEnzConc(const std::string& Name) {
        return Enz2Count[Name] / Vol;
    }

    float GetSubConc(const std::string& Name) {
        return Enz2Count[Name] / Vol;
    }

    void PrintState() {
        std::cout << "EnzCounts:" << std::endl;
        for (auto& Enz : Enz2Count) {
           std::cout<< "  " << "PRINT HERE" << std::endl; 
        };
        
        std::cout << "SubCounts:" << std::endl;
        for (auto& Sub : Sub2Count) {
           std::cout<< "  " << "PRINT HERE" << std::endl; 
        };
    }
};

class FSimulation {
public:
    int N_SimSteps;
    FState State;

    FSimulation(){};

    int SetSimSteps(const int& SimSteps) {
        N_SimSteps = SimSteps;
    }

    // float 

    void Init(FState& State, const int& SimSteps) {
        SetSimSteps(SimSteps);

        // Set up all molecule counts for now
        State.SetEnzCount("GltA", 50);
        State.SetEnzCount("AcnA", 50);
        State.SetEnzCount("Icd", 50);
        State.SetEnzCount("SucA", 50);
        State.SetEnzCount("SucD", 50);
        State.SetEnzCount("Sdh", 50);
        State.SetEnzCount("FumA", 50);
        State.SetEnzCount("Mdh", 50);

        State.SetSubCount("oxaloacetate", 10000);
        State.SetSubCount("citrate", 10000);
        State.SetSubCount("isocitrate", 10000);
        State.SetSubCount("keto-glutarate", 10000);
        State.SetSubCount("succinyl-CoA", 10000);
        State.SetSubCount("succinate", 10000);
        State.SetSubCount("fumarate", 10000);
        State.SetSubCount("malate", 10000);

    }
    
    void Run(FState& State) {
//        CurrentSimStep = 0
//        for (CurrentSimStep < N_SimSteps) {
//            for (auto& Pathway : Context.PathwayList["TCA"]) {
//                for (auto& enzyme : Pathway->Enzyme){
//                    Context.EnzymeList[enzyme];
//                };
//            };
//
//            Rate = MichaelisMentenEqn() 
//
//            CurrentSimStep++;
//        };
//
//
    }

};







