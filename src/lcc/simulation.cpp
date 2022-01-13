#include <iostream>
#include <fstream>
#include "simulation.h"
#include "datamanager.h"

using namespace std;

extern FDataset Dataset;
extern FSimulation Simulation;
extern FDataManager DataManager;

float ConcToCount(float Conc_Molecule, float Volume) {
    return Conc_Molecule * Volume;
}

float CountToConc(int Count_Molecule, float Volume) {
    return Count_Molecule / Volume;
}
 
float MichaelisMentenEqn(float Conc_Enzyme, float Conc_Substrate, float kcat, float KM) {
    // float Conc_Enzyme = CountToConcentration(Count_Enzyme, Volume);
    // float Conc_Substrate = CountToConcentration(Count_Substrate, Volume);

    float Rate = (kcat * Conc_Enzyme * Conc_Substrate) / (Conc_Substrate + KM);
    // std::cout << "Rate: " << Rate << std::endl;
    return Rate;
}

void FDataset::PrintLegend(){
    std::cout << "DataStrip: " << std::endl;
    for (auto& legend : Legend) {
        std::cout << legend << ", ";
    }
    std::cout << std::endl;
}

void FDataset::PrintData(){
    std::cout << "DataStrip: " << std::endl;
    for (auto& data : Data) {
        std::cout << data << ", ";
    }
    std::cout << std::endl;
}

void FSimulation::SetSimSteps(const int& SimSteps) {
    N_SimSteps = SimSteps;
}

void FSimulation::Init(FState& State, FDataset& Dataset, FDataManager& DataManager, const int& SimSteps)
{
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

    State.SetSubCount("oxaloacetate", 500);
    State.SetSubCount("citrate", 500);
    State.SetSubCount("isocitrate", 500);
    State.SetSubCount("keto-glutarate", 500);
    State.SetSubCount("succinyl-CoA", 500);
    State.SetSubCount("succinate", 500);
    State.SetSubCount("fumarate", 500);
    State.SetSubCount("malate", 500);

    State.SetSubCount("acetyl-CoA", 500);
    State.SetSubCount("H2O", 500);
    State.SetSubCount("CoA", 500);
    State.SetSubCount("NAD+", 500);
    State.SetSubCount("NADH", 500);
    State.SetSubCount("CO2", 500);
    State.SetSubCount("ADP", 500);
    State.SetSubCount("ATP", 500);
    State.SetSubCount("Pi", 500);
    State.SetSubCount("FAD", 500);
    State.SetSubCount("FADH2", 500);
    State.SetSubCount("H+", 500);

    // Legend Export
    Dataset.Legend = State.ExportLegend();
    //Dataset.PrintLegend();
    DataManager.SetLegend(Dataset.Legend);

    // Data Export for Sim Step 0
    Dataset.Data = State.ExportState();
    //Dataset.PrintData();
    DataManager.Add(Dataset.Data);

}

void FSimulation::Run(FState& State, FCompilerContext& Context, FDataset& Dataset)
{
//    std::cout << "Simulation: " << std::endl;
//
//    // TODO: export this variable in simulation object instead of cell state Sim step
//    int CurrentSimStep = 0;
//
//    while (CurrentSimStep < N_SimSteps) {
//        std::cout << "Step: " << CurrentSimStep << std::endl;
//
//        for (auto& Pathway : Context.PathwayList) {
//            std::cout << "Pathway: " << Pathway.Name << std::endl;
//            // Pathway contains a name (a string) and enzyme names (a vector of strings)
//
//            for(auto& EnzymeName : Pathway.Sequence) {
//                std::cout << "  Enzyme: " << EnzymeName << "\t|"; //std::endl;
//                std::string Substrate;
//                float kcat;
//                float KM;
//
//                // loop enzyme list to find its substrate
//                for (auto& Enzyme: Context.EnzymeList) {
//                    if (Enzyme.Name == EnzymeName) {
//                        Substrate = Enzyme.Substrate;
//                        kcat = Enzyme.kcat;
//                        KM = Enzyme.KM;
//                        //std::cout << "  Substrate: " << Substrate << std::endl;
//                        //std::cout << "  kcat: " << kcat << std::endl;
//                        //std::cout << "  KM: " << KM << std::endl;
//                    }
//                }
//
//                // std::cout << "  Substrate: " << Substrate << std::endl;
//                // std::cout << "  kcat: " << kcat << std::endl;
//                // std::cout << "  KM: " << KM << std::endl;
//
//                // get enzyme and substrate concentration
//                float EnzConc = State.GetEnzConc(EnzymeName);
//                float SubConc = State.GetSubConc(Substrate);
//
//                std::cout << "    EnzConc: " << EnzConc << "\t|" ; //std::endl;
//                std::cout << "    SubConc: " << SubConc << "\t|" ; // std::endl;
//
//                // calculate reaction rate
//                float Rate = MichaelisMentenEqn(EnzConc, SubConc, kcat, KM);
//                std::cout << "    Rate: " << Rate << std::endl;
//
//                // currently each pathway sub reaction occurs linearly.
//                // set up delta count variable to change it to concurrent simulation in sim.py
//
//                // loop enzymatic reaction list to find the stoichiometry
//                for (auto& EnzymaticReaction: Context.EnzymaticReactionList) {
//                    if ((EnzymaticReaction.Enzyme == EnzymeName) & EnzymaticReaction.CheckIfReactant(Substrate)) {
////                         std::map<std::string, int> Stoich = EnzymaticReaction.Stoichiometry;
////                            std::map<std::string, int>::iterator Stoich_Iter = Stoich.begin();
//                        for (std::pair<std::string, int> stoich : EnzymaticReaction.Stoichiometry) {
//                            std::string Molecule = stoich.first;
//                            int Coeff = stoich.second;
//
//                            // multiply reaction rate to stoichiometry to get delta concentration
//                            float MolConc = State.GetSubConc(Molecule);
//                            MolConc += Rate * Coeff;
//                            int MolCount = ConcToCount(MolConc, State.Vol);
//                            State.SetSubCount(Molecule, MolCount);
//                        }
//                    }
//                } // Enzymatic Reaction for loop
//            } // Pathway.Sequence for loop
//        } // Context.Pathway loop
//
//        CurrentSimStep++;
//
//        // Data Export
//        Dataset.Data = State.ExportState();
//        //Dataset.PrintData();
//        DataManager.Add(Dataset.Data);
//
//    } // while loop
}
