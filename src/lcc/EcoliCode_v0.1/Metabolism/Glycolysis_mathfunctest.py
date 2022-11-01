import sys
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import math

PerturbationTag = "#"
PerturbationName = "ATP_Consumption"

class FNetwork():
    def __init__(self):
        # Molecule Concentration
        self.MolConc = dict()   # current concentration
        self.Dataset = dict()   # storing concentration
        self.Perturbation = dict()

        # Constants
        self.Constants = dict()

        # Simulation Parameters
        self.SimStep = 0
        self.TotalSimulationTime = 0
        self.SimulationTimeUnit = 0
        self.SteadyStateBrake = False
        self.SteadyStateCheckMolecule = None
        self.SteadyStateThresholdFactor = 0

    def SetSimulationParameters(self, TotalSimulationTime, SimulationTimeUnit, SteadyStateBrake, SteadyStateThresholdFactor):
        self.TotalSimulationTime = TotalSimulationTime
        self.SimulationTimeUnit = SimulationTimeUnit
        self.SteadyStateBrake = SteadyStateBrake
        self.SteadyStateThresholdFactor = SteadyStateThresholdFactor

    def PrintSimulationParameters(self):
        print("Total Simulation Time :", self.TotalSimulationTime, "s", " | Simulation Time Unit :", self.SimulationTimeUnit, "s", end="")

    def PrintSimulationSubject(self, Subject):
        print("\n\nSimulating " + Subject + "...")

    def InitializeDataset(self):
        for MolName, Conc in self.MolConc.items():
            self.Dataset[MolName] = [Conc]
        for PerturbationName, Conc in self.Perturbation.items():
            self.Dataset[PerturbationName] = [Conc]

    def InitializeSimStep(self):
        self.SimStep = 0

    def SetConstants(self, List_ConstantNames, List_ConstantValues):
        for ConstantName, ConstantValue in zip(List_ConstantNames, List_ConstantValues):
            self.Constants[ConstantName] = ConstantValue
            
    def AppendData(self):
        for MolName, Conc in self.MolConc.items():
            self.Dataset[MolName].append(Conc)
        for PerturbationName, Conc in self.Perturbation.items():
            self.Dataset[PerturbationName].append(Conc)

    def UpdateMolConc(self, NewMolConcentrations):
        # If concentrations getting below zero
        if self.CheckBelowZero(NewMolConcentrations):
            pass

        # If concentrations all above zero
        else:
            for MolName, Conc in NewMolConcentrations.items():
                self.MolConc[MolName] += NewMolConcentrations[MolName]

    def CheckBelowZero(self, NewMolConcentrations):
        bBelowZero = False
        for MolName, Conc in NewMolConcentrations.items():
            Test = self.MolConc[MolName] + NewMolConcentrations[MolName]
            if Test < 0:
                # print("Below zero: %s %s\t" % (MolName, str(Test)), end= " |\t")
                # self.PrintSimTimeAndMolConc()
                bBelowZero = True
                break
        return bBelowZero

    def PrintMolConc(self):
        for MolName, Conc in self.MolConc.items():
            print("[" + MolName + "] " + "{:.5f}".format(Conc), end=", \t")

    def PrintSimTime(self):
        print("\n", self.SimStep, "\t| ", "{:.2f}".format(self.SimStep * self.SimulationTimeUnit), end=" s | ")

    def PrintSimTimeAndMolConc(self):
        self.PrintSimTime()
        self.PrintMolConc()

    def TagPerturbationLabel(self):
        for PerturbationName, Conc in self.Perturbation.items():
            if self.Dataset[PerturbationName]:
                self.Dataset[PerturbationTag + PerturbationName] = self.Dataset[PerturbationName]
                del self.Dataset[PerturbationName]

    def Glycolysis(self, ATPProductionType, ATPConsumptionType):
        self.PrintSimulationSubject("Glycolysis" + "_" + ATPProductionType + "_" + ATPConsumptionType)
        self.SteadyStateCheckMolecule = "G6P"
        self.Initialize_Glycolysis()
        self.Run_Glycolysis(self.Convert2Int_ATPProductionType(ATPProductionType), self.Convert2Int_ATPConsumptionType(ATPConsumptionType))
        self.TagPerturbationLabel()
        return self.Dataset.copy()

    def Convert2Int_ATPProductionType(self, ATPProductionType):
        if ATPProductionType == "No_ATPProduction":
            return 0
        elif ATPProductionType == "ATP'=cg":
            return 1
        elif ATPProductionType == "ATP'=[ADP]/[ATP]*cg":
            return 2
        elif ATPProductionType == "dATP=[G6P]*cg":
            return 3
        else:
            print("[ERROR] Unsupported ATP Production Type for Glycolysis Model")
            sys.exit(1)

    def Convert2Int_ATPConsumptionType(self, ATPConsumptionType):
        if ATPConsumptionType == "No_ATPConsumption":
            return 0
        elif ATPConsumptionType == "Linear_ATPConsumption":
            return 1
        elif ATPConsumptionType == "Burst_ATPConsumption":
            return 2
        else:
            print("[ERROR] Unsupported ATP Consumption Type for Glycolysis Model")
            sys.exit(1)

    def Initialize_Glycolysis(self):
        # Initial Molecule Concentrations
        self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["ADP"]      = 5.6e-4
        self.MolConc["NAD"]      = 2.6e-3
        self.MolConc["Pi"]       = 6e-3   # undocumented
        self.MolConc["pyruvate"] = 3.9e-4
        self.MolConc["NADH"]     = 8.3e-3
        self.MolConc["ATP"]      = 9.6e-3

        # Perturbations
        self.Perturbation[PerturbationName] = 0

        # Constants
        self.Constants["cg"]     = 10

        self.InitializeDataset()
        self.InitializeSimStep()


    def Simulation_Glycolysis(self, ATPProductionType):
        d = dict()

        # Production Rate Type
        Input = 0
        if ATPProductionType == 1:
            Input = 0.01
        elif ATPProductionType == 2:
            Input = self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        elif ATPProductionType == 3:
            Input = self.Dataset["G6P"][-1] * 2

        d["ATP"] = (Input * self.Constants["cg"]) * self.SimulationTimeUnit
        d["ADP"] = - d["ATP"]
        d["Pi"] = -d["ATP"]
        d["NADH"] = d["ATP"]
        d["NAD"] = - d["NADH"]
        d["G6P"] = - d["ATP"] / 2
        d["pyruvate"] = - d["G6P"] * 2

        self.UpdateMolConc(d)

    def Simulation_ATPConsumption(self, ATPConsumptionType=0):
        d = dict()

        # Perturbation
        if ATPConsumptionType == 1:   # linear
            self.Perturbation[PerturbationName] = (self.TotalSimulationTime * 1.4) * self.SimulationTimeUnit
        elif ATPConsumptionType == 2:   # non-linear
            # self.Perturbation[PerturbationName] = (self.TotalSimulationTime * 0.1) * self.SimulationTimeUnit * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
            self.Perturbation[PerturbationName] = (self.TotalSimulationTime * 0.1) * self.SimulationTimeUnit * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4

        d["ATP"] = - self.Perturbation[PerturbationName]
        d["ADP"] = self.Perturbation[PerturbationName]
        d["Pi"] = self.Perturbation[PerturbationName]
        d["NADH"] = - self.Perturbation[PerturbationName]
        d["NAD"] = self.Perturbation[PerturbationName]

        self.UpdateMolConc(d)

    def Run_Glycolysis(self, ATPProductionType, ATPConsumptionType):
        self.PrintSimTimeAndMolConc()

        while self.SimStep * self.SimulationTimeUnit < self.TotalSimulationTime:
            # Debugging
            # self.PrintSimTimeAndMolConc()

            self.Simulation_Glycolysis(ATPProductionType)
            self.Simulation_ATPConsumption(ATPConsumptionType)
            self.AppendData()
            if self.SteadyStateBrake:
                if (abs(self.Dataset[self.SteadyStateCheckMolecule][-1] - self.Dataset[self.SteadyStateCheckMolecule][-2]) / self.Dataset[self.SteadyStateCheckMolecule][-1] < self.SteadyStateThresholdFactor):
                    break
            self.SimStep += 1

        self.PrintSimTimeAndMolConc()


class Plotter:
    def __init__(self):
        self.Filter_Inclusion = None
        self.Filter_Exclusion = None

    def SetFilter_Inclusion(self, List):
        self.Filter_Inclusion = List

    def SetFilter_Exclusion(self, List):
        self.Filter_Exclusion = List

    def FilterDatasets(self, Datasets):
        Datasets_Filtered = dict()
        for Key_Dataset, Dataset in Datasets.items():
            Dataset_Filtered = dict()
            for Key_Data, Data in Dataset.items():
                if self.Filter_Inclusion:
                    if (Key_Data in self.Filter_Inclusion) or (Key_Data[1:] in self.Filter_Inclusion):
                        Dataset_Filtered[Key_Data] = Data
                if self.Filter_Exclusion:
                    if (Key_Data not in self.Filter_Exclusion) or (Key_Data[1:] + Key_Data not in self.Filter_Exclusion):
                        Dataset_Filtered[Key_Data] = Data
            Datasets_Filtered[Key_Dataset] = Dataset_Filtered

        return Datasets_Filtered

    def PlotData(self, Datasets, SimulationTimeUnit):

        # Filter Datasets
        if self.Filter_Inclusion or self.Filter_Exclusion:
            Datasets = self.FilterDatasets(Datasets)

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.5, hspace=0.5)


        Time = None   # Universal X axis (time)
        Perturbation = False
        for Dataset in Datasets.values():
            for Data in Dataset.values():
                Time = [i * SimulationTimeUnit for i in range(len(Data))]
                break
            for Key in Dataset.keys():
                if Key[0] == PerturbationTag:
                    Perturbation = True
                    break

        # Plot data
        NPlotsInRows = len(Datasets)   # Default
        MaxNPlotsInRows = 3
        if len(Datasets) > 1:
            for Remainder in range(MaxNPlotsInRows):
                if len(Datasets) % (Remainder + 1) == 0:
                    NPlotsInRows = Remainder + 1

        for n, (Process, Dataset) in enumerate(Datasets.items()):
            ax1 = fig.add_subplot(math.ceil(len(Datasets) / NPlotsInRows), NPlotsInRows, n + 1)
            ax2 = None
            if Perturbation:
                ax2 = ax1.twinx()

            # Y axis (molecular concentrations)
            for MolName, Conc in Dataset.items():
                if MolName[0] != PerturbationTag:
                    ax1.plot(Time, Conc, label="[" + MolName + "]")

                else:
                    ax2.plot(Time, Conc, color='brown', label="[" + MolName[1:] + "]")
                    # ax2.plot(Time, Conc, label="[" + MolName[2:] + "]")

            ax1.set_title(Process)
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Concentration (mol L-1)')
            ax1.set_ylim([0, 0.015])
            ax1.legend(loc='upper left')
            if Perturbation:
                ax2.set_ylabel('Concentration (mol L-1)')
                ax2.set_ylim([0, 3e-7])
                ax2.legend(loc='upper right')

            # ax1.grid()

        plt.show()


def main():
    Datasets = dict()

    # Initialization
    Network = FNetwork()

    # Setting Simulation Parameters
    TotalSimulationTime = 1e-1   # s
    SimulationTimeUnit = 1e-6   # s
    SteadyStateBrake = False
    SteadyStateThresholdFactor = 0.0000001  # Steady state threshold
    Network.SetSimulationParameters(TotalSimulationTime, SimulationTimeUnit, SteadyStateBrake, SteadyStateThresholdFactor)
    Network.PrintSimulationParameters()

    # Glycolysis
    ModelName = "Glycolysis"
    ATPProductionTypes = ["ATP'=cg", "ATP'=[ADP]/[ATP]*cg"]
    # ATPProductionTypes = ["dATP=[G6P]*cg", "dATP=[G6P]/[ATP]*cg"]
    ATPConsumptionTypes = ["No_ATPConsumption", "Linear_ATPConsumption", "Burst_ATPConsumption"]
    # ATPConsumptionTypes = ["Linear_ATPConsumption", "Burst_ATPConsumption"]

    # Debugging

    for ATPProductionType in ATPProductionTypes:
        for ATPConsumptionType in ATPConsumptionTypes:
            Title = ModelName + "\n" + ATPProductionType + "\n" + ATPConsumptionType
            Datasets[Title] = Network.Glycolysis(ATPProductionType, ATPConsumptionType)

    # Plotting
    Plot = Plotter()

    # InclusionList = ["G6P", "ATP"]
    # InclusionList = [PerturbationName]
    # InclusionList = ["G6P", "ATP", "pyruvate", PerturbationName]
    # InclusionList = ["pyruvate", "ATP", PerturbationName]
    # InclusionList = ["G6P", "ATP", PerturbationName]
    InclusionList = ["ATP", PerturbationName]
    Plot.SetFilter_Inclusion(InclusionList)

    ExclusionList = ["Pi, NADH, NAD"]
    # Plot.SetFilter_Exclusion(ExclusionList)

    Plot.PlotData(Datasets, SimulationTimeUnit)


if __name__ == '__main__':
    main()
