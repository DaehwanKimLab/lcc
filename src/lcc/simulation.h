//#ifndef LCC_SIMULATION_H
#define LCC_SIMULATION_H

#include <string>
#include <vector>
#include <map>

// #include "option.h"
// #include "node.h"
#include "context.h"

float ConcToCount(float Conc_Molecule, float Volume) {
    return Conc_Molecule * Volume;
};

float CountToConc(int Count_Molecule, float Volume) {
    return Count_Molecule / Volume;
};
 
float MichaelisMentenEqn(float Conc_Enzyme, float Conc_Substrate, float kcat, float kM) {
    // float Conc_Enzyme = CountToConcentration(Count_Enzyme, Volume);
    // float Conc_Substrate = CountToConcentration(Count_Substrate, Volume);

    float Rate = (kcat * Conc_Enzyme * Conc_Substrate) / (Conc_Substrate + kM);
    // std::cout << "Rate: " << Rate << std::endl;
    return Rate;
};

class FState {
public:

// private:
    std::map<std::string, int> Enz2Count; // The count info may be merged but in the future sim.py
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
        return Sub2Count[Name];
    }

    float GetEnzConc(const std::string& Name) {
        return Enz2Count[Name] / Vol;
    }

    float GetSubConc(const std::string& Name) {
        return Sub2Count[Name] / Vol;
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

        State.SetVol(1);

        // Set up all molecule counts for now
        State.SetEnzCount("GltA", 1);
        State.SetEnzCount("AcnA", 7);
        State.SetEnzCount("Icd", 6);
        State.SetEnzCount("SucA", 1);
        State.SetEnzCount("SucD", 2);
        State.SetEnzCount("Sdh", 3);
        State.SetEnzCount("FumA", 3);
        State.SetEnzCount("Mdh", 5);

        State.SetSubCount("oxaloacetate", 100);
        State.SetSubCount("citrate", 100);
        State.SetSubCount("isocitrate", 100);
        State.SetSubCount("keto-glutarate", 100);
        State.SetSubCount("succinyl-CoA", 100);
        State.SetSubCount("succinate", 100);
        State.SetSubCount("fumarate", 100);
        State.SetSubCount("malate", 100);

        State.SetSubCount("acetyl-CoA", 100);
    }
    
    void Run(FState& State, FCompilerContext Context) {
        std::cout << "Simulation: " << std::endl;
        int CurrentSimStep = 0;

        while (CurrentSimStep < N_SimSteps) {
            std::cout << "Step: " << CurrentSimStep << std::endl;

            for (auto& Pathway : Context.PathwayList) {
                std::cout << "Pathway: " << Pathway.Name << std::endl;
                // Pathway contains a name (a string) and enzyme names (a vector of strings)

                for (auto& EnzymeName : Pathway.Sequence) {
                    std::cout << "  Enzyme: " << EnzymeName << std::endl; 
                    std::string Substrate;
                    float kcat;
                    float kM;

                    // loop enzyme list to find its substrate
                    for (auto& Enzyme: Context.EnzymeList) {
                        if (Enzyme.Name == EnzymeName) {
                            Substrate = Enzyme.Substrate;
                            kcat = Enzyme.kcat;
                            kM = Enzyme.kM;
                            std::cout << "  Substrate: " << Substrate << std::endl;
                            std::cout << "  kcat: " << kcat << std::endl;
                            std::cout << "  kM: " << kM << std::endl;
                        };
                    }; 

                    // std::cout << "  Substrate: " << Substrate << std::endl;
                    // std::cout << "  kcat: " << kcat << std::endl;
                    // std::cout << "  kM: " << kM << std::endl;

                    // get enzyme and substrate concentration
                    float EnzConc = State.GetEnzConc(EnzymeName);
                    float SubConc = State.GetSubConc(Substrate);

                    std::cout << "    EnzConc: " << EnzConc << std::endl;
                    std::cout << "    SubConc: " << SubConc << std::endl;
 
                    // calculate reaction rate
                    float Rate = MichaelisMentenEqn(EnzConc, SubConc, kcat, kM);
                    std::cout << "    Rate: " << Rate << std::endl;

                    // currently each pathway sub reaction occurs linearly. 
                    // set up delta count variable to change it to concurrent simulation in sim.py

                    // loop enzymatic reaction list to find the stoichiometry
                    for (auto& EnzymaticReaction: Context.EnzymaticReactionList) { 
                        if ((EnzymaticReaction.Enzyme == EnzymeName) & EnzymaticReaction.CheckIfReactant(Substrate)) {
   //                         std::map<std::string, int> Stoich = EnzymaticReaction.Stoichiometry;
//                            std::map<std::string, int>::iterator Stoich_Iter = Stoich.begin();
                            for (std::pair<std::string, int> stoich : EnzymaticReaction.Stoichiometry) {
                                std::string Molecule = stoich.first;
                                int Coeff = stoich.second;

                                // multiply reaction rate to stoichiometry to get delta concentration
                                float MolConc = State.GetSubConc(Molecule);
                                MolConc += Rate * Coeff;
                                int MolCount = ConcToCount(MolConc, State.Vol);
                                State.SetSubCount(Molecule, MolCount);
                            };
                        };
                    };
                };
            };
            CurrentSimStep++;
        };


    }

};







