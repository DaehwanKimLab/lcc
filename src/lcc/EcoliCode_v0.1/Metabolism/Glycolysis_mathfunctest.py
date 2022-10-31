import sys
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import math

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
        for MolName, Conc in NewMolConcentrations.items():
            self.MolConc[MolName] += NewMolConcentrations[MolName]

    def PrintMolConc(self):
        for MolName, Conc in self.MolConc.items():
            print("[" + MolName + "] " + "{:.5f}".format(Conc), end=", \t")

    def PrintSimTime(self):
        print("\n", self.SimStep, "\t| ", "{:.2f}".format(self.SimStep * self.SimulationTimeUnit), end=" s | ")

    def PrintSimTimeAndMolConc(self):
        self.PrintSimTime()
        self.PrintMolConc()

    def Glycolysis(self, ATPProductionType, ATPConsumptionType):
        self.PrintSimulationSubject("Glycolysis" + "_" + ATPProductionType + "_" + ATPConsumptionType)
        self.SteadyStateCheckMolecule = "G6P"
        self.Initialize_Glycolysis()
        self.Run_Glycolysis(self.Convert2Int_ATPProductionType(ATPProductionType), self.Convert2Int_ATPConsumptionType(ATPConsumptionType))
        return self.Dataset.copy()

    def Convert2Int_ATPProductionType(self, ATPProductionType):
        if ATPProductionType == "No_ATPProduction":
            return 0
        elif ATPProductionType == "dATP=[ADP]*cg":
            return 1
        elif ATPProductionType == "dATP=[ADP]/[ATP]*cg":
            return 2
        elif ATPProductionType == "dATP=[G6P]*cg":
            return 3
        elif ATPProductionType == "dATP=[G6P]/[ATP]*cg":
            return 4
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
        self.Perturbation["Process"] = 0

        # Constants
        self.Constants["cg"]     = 10

        self.InitializeDataset()
        self.InitializeSimStep()


    def Simulation_Glycolysis(self, ATPProductionType=2):
        d = dict()

        # Production Rate Type
        Input = 0
        if ATPProductionType == 1:
            Input = self.Dataset["ADP"][-1] * 10
        elif ATPProductionType == 2:
            Input = self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        elif ATPProductionType == 3:
            Input = self.Dataset["G6P"][-1] * 2
        elif ATPProductionType == 4:
            Input = self.Dataset["G6P"][-1] / self.Dataset["ATP"][-1] * 0.001

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
            self.Perturbation["Process"] = (self.TotalSimulationTime * 0.2) * self.SimulationTimeUnit
        elif ATPConsumptionType == 2:   # non-linear
            # self.Perturbation["Process"] = (self.TotalSimulationTime * 0.1) * self.SimulationTimeUnit * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
            self.Perturbation["Process"] = (self.TotalSimulationTime * 0.0125) * self.SimulationTimeUnit * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4

        d["ATP"] = - self.Perturbation["Process"]
        d["ADP"] = self.Perturbation["Process"]

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
                    if Key_Data in self.Filter_Inclusion:
                        Dataset_Filtered[Key_Data] = Data
                if self.Filter_Exclusion:
                    if Key_Data not in self.Filter_Exclusion:
                        Dataset_Filtered[Key_Data] = Data
            Datasets_Filtered[Key_Dataset] = Dataset_Filtered

        return Datasets_Filtered

    def PlotData(self, Datasets, SimulationTimeUnit):

        # Filter Datasets
        if self.Filter_Inclusion or self.Filter_Exclusion:
            Datasets = self.FilterDatasets(Datasets)

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.5)

        # Universal X axis (time)
        Time = None
        for Dataset in Datasets.values():
            for Data in Dataset.values():
                Time = [i * SimulationTimeUnit for i in range(len(Data))]
                break

        # Plot data
        Rows = 1
        if len(Datasets) > 1:
            Rows = 3
        for n, (Process, Dataset) in enumerate(Datasets.items()):
            ax = fig.add_subplot(math.ceil(len(Datasets) / Rows), Rows, n + 1)

            # Y axis (molecular concentrations)
            for MolName, Conc in Dataset.items():
                ax.plot(Time, Conc, label="[" + MolName + "]")

            ax.set_title(Process)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Concentration (mol L-1)')
            ax.legend(loc='upper center')
            ax.grid()

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
    ATPProductionTypes = ["dATP=[ADP]*cg", "dATP=[ADP]/[ATP]*cg"]
    # ATPProductionTypes = ["dATP=[G6P]*cg", "dATP=[G6P]/[ATP]*cg"]
    ATPConsumptionTypes = ["No_ATPConsumption", "Linear_ATPConsumption", "Burst_ATPConsumption"]
    for ATPProductionType in ATPProductionTypes:
        for ATPConsumptionType in ATPConsumptionTypes:
            Title = ModelName + "\n" + ATPProductionType + "\n" + ATPConsumptionType
            Datasets[Title] = Network.Glycolysis(ATPProductionType, ATPConsumptionType)

    # Plotting
    Plot = Plotter()

    # InclusionList = ["G6P", "ATP"]
    InclusionList = ["Process"]
    # Plot.SetFilter_Inclusion(InclusionList)

    ExclusionList = ["Pi, NADH, NAD"]
    Plot.SetFilter_Exclusion(ExclusionList)

    Plot.PlotData(Datasets, SimulationTimeUnit)


if __name__ == '__main__':
    main()
