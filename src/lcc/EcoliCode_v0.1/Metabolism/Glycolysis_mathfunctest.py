from argparse import ArgumentParser
import matplotlib.pyplot as plt
import math

class FNetwork():
    def __init__(self):
        # Molecule Concentration
        self.MolConc = dict()   # current concentration
        self.Dataset = dict()   # storing concentration

        # Constants
        self.Constants = dict()

        # Simulation Parameters
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
        print("Total Simulation Time :", self.TotalSimulationTime, "s", " | Simulation Time Unit :", self.SimulationTimeUnit, "s")

    def PrintSimulationSubject(self, Subject):
        print("Simulating " + Subject + "...")

    def InitializeDataset(self):
        for MolName, Conc in self.MolConc.items():
            self.Dataset[MolName] = [Conc]

    def SetConstants(self, List_ConstantNames, List_ConstantValues):
        for ConstantName, ConstantValue in zip(List_ConstantNames, List_ConstantValues):
            self.Constants[ConstantName] = ConstantValue
            
    def AppendData(self):
        for MolName, Conc in self.MolConc.items():
            self.Dataset[MolName].append(Conc)

    def UpdateMolConc(self, NewMolConcentrations):
        for MolName, Conc in self.MolConc.items():
            self.MolConc[MolName] += NewMolConcentrations[MolName]

    def PrintMolConc(self):
        for MolName, Conc in self.MolConc.items():
            print("[" + MolName + "] " + "{:.5f}".format(Conc), end=", \t")

    def PrintSimTime(self, SimStep):
        print("\n", SimStep, " | ", "{:.2f}".format(SimStep * self.SimulationTimeUnit), end=" s | ")

    # Glycolysis
    def Glycolysis(self):
        self.PrintSimulationSubject("Glycolysis")
        self.SteadyStateCheckMolecule = "G6P"
        self.Initialize_Glycolysis()
        self.Run_Glycolysis()
        return self.Dataset

    def Initialize_Glycolysis(self):
        # Initial Molecule Concentrations
        self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["ADP"]      = 5.6e-4
        self.MolConc["NAD"]      = 2.6e-3
        self.MolConc["Pi"]       = 6e-3   # undocumented
        self.MolConc["pyruvate"] = 3.9e-4
        self.MolConc["NADH"]     = 8.3e-3
        self.MolConc["ATP"]      = 9.6e-3

        # Constants
        self.Constants["cg"]     = 10
        self.InitializeDataset()

    def Simulation_Glycolysis(self):
        d = dict()

        d["ATP"] = (self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1] * self.Constants["cg"]) * self.SimulationTimeUnit
        d["ADP"] = - d["ATP"]
        d["Pi"] = -d["ATP"]
        d["NADH"] = d["ATP"]
        d["NAD"] = - d["NADH"]
        d["G6P"] = - d["ATP"] / 2
        d["pyruvate"] = - d["G6P"] * 2

        self.UpdateMolConc(d)

    def Run_Glycolysis(self):
        SimStep = 0
        while SimStep * self.SimulationTimeUnit < self.TotalSimulationTime:
            # Debug
            self.PrintSimTime(SimStep)
            self.PrintMolConc()

            self.Simulation_Glycolysis()
            self.AppendData()
            if self.SteadyStateBrake:
                if (abs(self.Dataset[self.SteadyStateCheckMolecule][-1] - self.Dataset[self.SteadyStateCheckMolecule][-2]) / self.Dataset[self.SteadyStateCheckMolecule][-1] < self.SteadyStateThresholdFactor):
                    break
            SimStep += 1
        

def PlotData(Datasets, SimulationTimeUnit):
    fig = plt.figure()
    fig.subplots_adjust(wspace=0.2, hspace=0.3)

    # Universal X axis (time)
    Time = None
    for Dataset in Datasets.values():
        for Data in Dataset.values():
            Time = [i * SimulationTimeUnit for i in range(len(Data))]
            break

    # Plot data
    Rows = 1
    if len(Datasets) > 1:
        Rows = 2
    for n, (Process, Dataset) in enumerate(Datasets.items()):
        ax = fig.add_subplot(math.ceil(len(Datasets) / Rows), Rows, n + 1)

        # Y axis (molecular concentrations)
        for MolName, Conc in Dataset.items():
            ax.plot(Time, Conc, label="[" + MolName + "]")

        ax.set_title(Process)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Concentration (mol L-1)')
        ax.legend(loc='center right')
        ax.grid()

    plt.show()


def main():
    Datasets = dict()

    # Initialization
    Network = FNetwork()

    # Setting Simulation Parameters
    TotalSimulationTime = 1e-2   # s
    SimulationTimeUnit = 1e-6   # s
    SteadyStateBrake = True
    SteadyStateThresholdFactor = 0.0000001  # Steady state threshold
    Network.SetSimulationParameters(TotalSimulationTime, SimulationTimeUnit, SteadyStateBrake, SteadyStateThresholdFactor)
    Network.PrintSimulationParameters()

    # Glycolysis
    Datasets["Glycolysis"] = Network.Glycolysis()
    Datasets["Glycolysis_ATPConsumption_Linear"] = Network.Glycolysis()

    # Plotting
    PlotData(Datasets, SimulationTimeUnit)


if __name__ == '__main__':
    main()
