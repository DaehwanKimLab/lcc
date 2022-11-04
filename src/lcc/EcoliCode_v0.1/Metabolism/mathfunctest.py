import random
import sys
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import math

PerturbationTag = "#"
PerturbationName = "ATP_Consumption"
G6P_Constant = 8.8e-3

class FSimulation():
    def __init__(self):
        # Molecule Concentration
        self.MolConc = dict()   # current concentration
        self.DeltaMolConc = dict()
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

        # Model Parameters
        self.Title = ""
        self.GlucoseAvailabilityType = 0
        self.GlucoseTransportType = 0
        self.ATPProductionType = 0
        self.ATPConsumptionType = 0
        self.AcetylCoAProductionType = 0

    def SetSimulationParameters(self, TotalSimulationTime, SimulationTimeUnit, SteadyStateBrake, SteadyStateThresholdFactor):
        self.TotalSimulationTime = TotalSimulationTime
        self.SimulationTimeUnit = SimulationTimeUnit
        self.SteadyStateBrake = SteadyStateBrake
        self.SteadyStateThresholdFactor = SteadyStateThresholdFactor

    def PrintSimulationParameters(self):
        print("Total Simulation Time :", self.TotalSimulationTime, "s", " | Simulation Time Unit :", self.SimulationTimeUnit, "s", end="")

    def PrintSimulationSubject(self):
        print("\n\nSimulating " + self.Title + "...")

    def InitializeDeltaMolConc(self):
        for MolName, Conc in self.MolConc.items():
            self.DeltaMolConc[MolName] = 0

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

    def UpdateMolConc(self):
        # If concentrations getting below zero
        if self.CheckBelowZero(self.DeltaMolConc):
            pass

        # If concentrations all above zero
        else:
            for MolName, Conc in self.DeltaMolConc.items():
                self.MolConc[MolName] += self.DeltaMolConc[MolName]
                self.DeltaMolConc[MolName] = 0

    def AddToDeltaMolConc(self, NewMolConcentrations):
        for MolName, Conc in NewMolConcentrations.items():
            self.DeltaMolConc[MolName] += NewMolConcentrations[MolName]

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

    def PrintDeltaMolConc(self):
        for MolName, Conc in self.DeltaMolConc.items():
            print("[" + MolName + "] " + "{:.5f}".format(Conc), end=", \t")

    def PrintSimTime(self):
        print("\n", self.SimStep, "\t| ", "{:.2f}".format(self.SimStep * self.SimulationTimeUnit), end=" s | ")

    def PrintSimTimeAndMolConc(self):
        self.PrintSimTime()
        self.PrintMolConc()
        self.PrintDeltaMolConc()

    def TagPerturbationLabel(self):
        for PerturbationName, Conc in self.Perturbation.items():
            if self.Dataset[PerturbationName]:
                self.Dataset[PerturbationTag + PerturbationName] = self.Dataset[PerturbationName]
                del self.Dataset[PerturbationName]

    def SetTitle(self, Title):
        self.Title = Title

    def SetGlucoseTransportType(self, Type):
        self.GlucoseTransportType = self.Convert2Int_GlucoseTransportType(Type)

    def SetGlucoseAvailabilityType(self, Type):
        self.GlucoseAvailabilityType = self.Convert2Int_GlucoseAvailabilityType(Type)

    def SetATPProductionType(self, Type):
        self.ATPProductionType = self.Convert2Int_ATPProductionType(Type)

    def SetATPConsumptionType(self, Type):
        self.ATPConsumptionType = self.Convert2Int_ATPConsumptionType(Type)

    def Sim(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCA=0, OxidativePhosphorylation=0, ATPConsumption=0):
        self.PrintSimulationSubject()
        self.SteadyStateCheckMolecule = "pyruvate"
        self.Initialize(Chemotaxis, Glycolysis, PyruvateOxidation, TCA, OxidativePhosphorylation)
        self.Run(Chemotaxis, Glycolysis, PyruvateOxidation, TCA, OxidativePhosphorylation, ATPConsumption)
        self.TagPerturbationLabel()
        return self.Dataset.copy()

    def Convert2Int_GlucoseAvailabilityType(self, GlucoseAvailabilityType):
        if GlucoseAvailabilityType == "No_GlucoseAvailability":
            return 0
        elif GlucoseAvailabilityType == "Uniform_GlucoseAvailability":
            return 1
        elif GlucoseAvailabilityType == "Linear_GlucoseAvailability":
            return 2
        elif GlucoseAvailabilityType == "NonLinear_GlucoseAvailability":
            return 3
        else:
            print("[ERROR] Unsupported Glucose Influx Type: %s" % GlucoseAvailabilityType)
            sys.exit(1)

    def Convert2Int_GlucoseTransportType(self, GlucoseTransportType):
        if GlucoseTransportType == "G6P'=0":
            return 0
        elif GlucoseTransportType == "G6P'=cgt":
            return 1
        elif GlucoseTransportType == "G6P'=1/[G6P]*cgt":
            return 2
        else:
            print("[ERROR] Unsupported glucose Transport Type: %s" % GlucoseTransportType)
            sys.exit(1)

    def Convert2Int_ATPProductionType(self, ATPProductionType):
        if ATPProductionType == "ATP'=0":
            return 0
        elif ATPProductionType == "ATP'=cg":
            return 1
        elif ATPProductionType == "ATP'=[ADP]/[ATP]*cg":
            return 2
        elif ATPProductionType == "dATP=[G6P]*cg":
            return 3
        else:
            print("[ERROR] Unsupported ATP Production Type for Glycolysis Model: %s" % ATPProductionType)
            sys.exit(1)

    def Convert2Int_ATPConsumptionType(self, ATPConsumptionType=0):
        if ATPConsumptionType == "No_ATPConsumption":
            return 0
        elif ATPConsumptionType == "Linear_ATPConsumption":
            return 1
        elif ATPConsumptionType == "Burst_ATPConsumption":
            return 2
        else:
            print("[ERROR] Unsupported ATP Consumption Type: %s" % ATPConsumptionType)
            sys.exit(1)

    def Initialize_Chemotaxis(self):
        # Initial Molecule Concentrations
        self.MolConc["G6P"] = 8.8e-3   # all Hexose-P combined
        self.MolConc["ATP"] = 9.6e-3
        self.MolConc["ADP"] = 5.6e-4

        # Perturbation
        self.Perturbation["glucose_Available"] = 0
        self.Perturbation["glucose_Influx"] = 0

        # Constants
        self.Constants["cc"] = 0.01
        self.Constants["cgt"] = 0.001

    def GetInputFactor_Chemotaxis(self):
        if self.GlucoseAvailabilityType == 0:
            return 0
        elif self.GlucoseAvailabilityType == 1:
            return 0.05
        elif self.GlucoseAvailabilityType == 2:
            return 0.000002 * self.SimStep
        elif self.GlucoseAvailabilityType == 3:
            return 0.05 * math.sin(self.SimStep/4000) + 0.5e-7 / self.SimulationTimeUnit

    def GetTransportFactor_Chemotaxis(self):
        if self.GlucoseTransportType == 0:
            return 0
        elif self.GlucoseTransportType == 1:
            return 0.1
        elif self.GlucoseTransportType == 2:
            return 0.01 / (self.Dataset["G6P"][-1] ** 2) * self.Constants["cgt"]

    def Simulation_Chemotaxis(self):
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_Chemotaxis()
        self.Perturbation["glucose_Available"] = InputFactor * self.SimulationTimeUnit

        TransportFactor = self.GetTransportFactor_Chemotaxis()
        TransportDemand = TransportFactor if TransportFactor <= 1 else 1
        self.Perturbation["glucose_Influx"] = self.Perturbation["glucose_Available"] * TransportDemand

        d["G6P"] = self.Perturbation["glucose_Influx"]
        d["ATP"] = - self.Constants["cc"] * self.SimulationTimeUnit
        d["ADP"] = - d["ATP"]   # use when G6P is constant

        self.AddToDeltaMolConc(d)

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
        if self.ATPConsumptionType:
            self.Perturbation["ATP_Consumption"] = 0

        # Constants
        self.Constants["cg"]     = 10

    def GetInputFactor_Glycolysis(self):
        if self.ATPProductionType == 1:
            return 0.01
        elif self.ATPProductionType == 2:
            return self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        elif self.ATPProductionType == 3:
            return self.Dataset["G6P"][-1] * 2

    def Simulation_Glycolysis(self):
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_Glycolysis()

        d["ATP"] = (InputFactor * self.Constants["cg"]) * self.SimulationTimeUnit
        d["ADP"] = - d["ATP"]
        d["Pi"] = - d["ATP"]
        d["NADH"] = d["ATP"]
        d["NAD"] = - d["NADH"]
        d["G6P"] = - d["ATP"] / 2 if d["ATP"] / 2 < self.Dataset["G6P"][-1] else self.Dataset["G6P"][-1]
        # d["pyruvate"] = - d["G6P"] * 2
        d["pyruvate"] = d["ATP"]   # use when G6P is constant

        self.AddToDeltaMolConc(d)

    def GetATPConsumptionRate(self):
        if self.ATPConsumptionType == 0:   # None
            return 0
        if self.ATPConsumptionType == 1:   # linear
            return (self.TotalSimulationTime * 1.4)
        elif self.ATPConsumptionType == 2:   # non-linear
            # self.Perturbation[PerturbationName] = (self.TotalSimulationTime * 0.1) * self.SimulationTimeUnit * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
            return (self.TotalSimulationTime * 0.1) * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4

    def Simulation_ATPConsumption(self):
        d = dict()

        self.Perturbation["ATP_Consumption"] = self.GetATPConsumptionRate() * self.SimulationTimeUnit

        d["ATP"] = - self.Perturbation["ATP_Consumption"]
        d["ADP"] = self.Perturbation["ATP_Consumption"]
        d["Pi"] = self.Perturbation["ATP_Consumption"]
        d["NADH"] = - self.Perturbation["ATP_Consumption"]
        d["NAD"] = self.Perturbation["ATP_Consumption"]

        self.AddToDeltaMolConc(d)

    def Initialize_PyruvateOxidation(self):
        # pyruvate + NAD + CoA --> acetyl-CoA + CO2 + NADH + H

        # Initial Molecule Concentrations
        # self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["pyruvate"] = 3.9e-4
        self.MolConc["acetyl-CoA"]  = 6.1e-4
        self.MolConc["NAD"]      = 2.6e-3
        self.MolConc["NADH"]     = 8.3e-3

        # Constants
        self.Constants["cp"]     = 0.9

    def GetInputFactor_PyruvateOxidation(self):
        if self.AcetylCoAProductionType == 1:
            return 0.01
        elif self.AcetylCoAProductionType == 2:
            return self.MolConc["pyruvate"][-1] / self.Dataset["ATP"][-1]
        elif self.AcetylCoAProductionType == 3:
            return self.MolConc["acetyl-CoA"][-1] * 2

    def Simulation_PyruvateOxidation(self):
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_PyruvateOxidation()

        d["acetyl-CoA"] = (InputFactor * self.Constants["cp"]) * self.SimulationTimeUnit
        d["pyruvate"] = - d["acetyl-CoA"]
        d["NAD"] = d["pyruvate"]
        d["NADH"] = - d["NAD"]

        self.AddToDeltaMolConc(d)

    def Initialize_TCA(self):
        # acetyl-CoA + 3 NAD + FAD + ADP + Pi + 2 H2O --> 2 CO2 + 3 NADH + FADH2 + ATP + 2 H + CoA

        # Initial Molecule Concentrations
        # self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["acetyl-CoA"]  = 6.1e-4
        self.MolConc["FAD"]     = 1.7e-4
        self.MolConc["NAD"]     = 2.6e-3
        self.MolConc["ADP"]     = 5.6e-4
        self.MolConc["Pi"]      = 6e-3   # undocumented
        self.MolConc["NADH"]    = 8.3e-3
        self.MolConc["FADH2"]   = 5e-4   # undocumented
        self.MolConc["ATP"]     = 9.6e-3
        self.MolConc["CoA"]     = 1.4e-3

        # Constants
        self.Constants["ct"]     = 10


    def Simulation_TCA(self):
        d = dict()

        # Production Rate Type
        ATPDemandFactor = self.GetATPDemandFactor_Glycolysis()




        d["acetyl-CoA"] = (self.Constants["cg"]) * self.SimulationTimeUnit
        d["ADP"] = - d["ATP"]
        d["Pi"] = -d["ATP"]
        d["NADH"] = d["ATP"]
        d["NAD"] = - d["NADH"]
        # d["G6P"] = - d["ATP"] / 2
        # d["pyruvate"] = - d["G6P"] * 2
        d["acetyl-CoA"] = d["ATP"]   # use when G6P is constant

        self.AddToDeltaMolConc(d)


    def Initialize_OxidativePhosphorylation(self):
        # acetyl-CoA + 3 NAD + FAD + ADP + Pi + 2 H2O --> 2 CO2 + 3 NADH + FADH2 + ATP + 2 H + CoA

        # Initial Molecule Concentrations
        # self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["acetyl-CoA"]  = 6.1e-4
        self.MolConc["FAD"]     = 1.7e-4
        self.MolConc["NAD"]     = 2.6e-3
        self.MolConc["ADP"]     = 5.6e-4
        self.MolConc["Pi"]      = 6e-3   # undocumented
        self.MolConc["NADH"]    = 8.3e-3
        self.MolConc["FADH2"]   = 5e-4   # undocumented
        self.MolConc["ATP"]     = 9.6e-3
        self.MolConc["CoA"]     = 1.4e-3

        # Constants
        self.Constants["ct"]     = 10

    def Simulation_OxidativePhosphorylation(self):
        d = dict()

        # Production Rate Type
        ATPDemandFactor = self.GetATPDemandFactor_Glycolysis()




        d["acetyl-CoA"] = (self.Constants["cg"]) * self.SimulationTimeUnit
        d["ADP"] = - d["ATP"]
        d["Pi"] = -d["ATP"]
        d["NADH"] = d["ATP"]
        d["NAD"] = - d["NADH"]
        # d["G6P"] = - d["ATP"] / 2
        # d["pyruvate"] = - d["G6P"] * 2
        d["acetyl-CoA"] = d["ATP"]   # use when G6P is constant

        self.AddToDeltaMolConc(d)

    def Initialize(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCA=0, OxidativePhosphorylation=0):
        if Chemotaxis:
            self.Initialize_Chemotaxis()
        if Glycolysis:
            self.Initialize_Glycolysis()
        if PyruvateOxidation:
            self.Initialize_PyruvateOxidation()
        if TCA:
            self.Initialize_TCA()
        if OxidativePhosphorylation:
            self.Initialize_OxidativePhosphorylation()
        self.InitializeDeltaMolConc()
        self.InitializeDataset()
        self.InitializeSimStep()


    def Run(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCA=0, OxidativePhosphorylation=0, ATPConsumption=0):
        self.PrintSimTimeAndMolConc()

        while self.SimStep * self.SimulationTimeUnit < self.TotalSimulationTime:
            # Debugging
            # self.PrintSimTimeAndMolConc()

            if Chemotaxis:
                self.Simulation_Chemotaxis()
            if Glycolysis:
                self.Simulation_Glycolysis()
            if PyruvateOxidation:
                self.Simulation_PyruvateOxidation()
            if TCA:
                self.Simulation_TCA()
            if OxidativePhosphorylation:
                self.Simulation_OxidativePhosphorylation()
            if ATPConsumption:
                self.Simulation_ATPConsumption()

            self.UpdateMolConc()
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
        Perturbation = 0
        PerturbationPlotColor = list()
        for Dataset in Datasets.values():
            for Data in Dataset.values():
                Time = [i * SimulationTimeUnit for i in range(len(Data))]
                break
            for Key in Dataset.keys():
                if Key[0] == PerturbationTag:
                    PerturbationPlotColor.append((random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 0.75)))
                    Perturbation += 1


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
            PerturbationIndex = 0
            for MolName, Conc in Dataset.items():
                if MolName[0] != PerturbationTag:
                    ax1.plot(Time, Conc, label="[" + MolName + "]")

                else:
                    ax2.plot(Time, Conc, color=PerturbationPlotColor[PerturbationIndex], label="[" + MolName[1:] + "]")
                    # ax2.plot(Time, Conc, label="[" + MolName[2:] + "]")
                    # print(PerturbationIndex)
                    PerturbationIndex += 1

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
    Simulation = FSimulation()

    # Setting Simulation Parameters
    TotalSimulationTime = 1e-1   # s
    SimulationTimeUnit = 1e-6   # s
    SteadyStateBrake = False
    SteadyStateThresholdFactor = 0.0000001  # Steady state threshold
    Simulation.SetSimulationParameters(TotalSimulationTime, SimulationTimeUnit, SteadyStateBrake, SteadyStateThresholdFactor)
    Simulation.PrintSimulationParameters()

    ModelSwitch = 11

    # Model Switch Options
    # 1:    Chemotaxis Unit Test
    # 2:    Glycolysis Unit Test
    # 11:   Chemotaxis + Glycolysis

    if ModelSwitch == 1:
        # Chemotaxis Only
        ModelName = "Chemotaxis"
        # GlucoseAvailabilityTypes = ["No_GlucoseAvailability", "Uniform_GlucoseAvailability", "Linear_GlucoseAvailability", "NonLinear_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Uniform_GlucoseAvailability", "Linear_GlucoseAvailability", "NonLinear_GlucoseAvailability"]
        # GlucoseAvailabilityTypes = ["Uniform_GlucoseAvailability", "Linear_GlucoseAvailability"]
        GlucoseTransportTypes = ["G6P'=0", "G6P'=cgt", "G6P'=1/[G6P]*cgt"]
        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                Title = ModelName + "\n" + GlucoseAvailabilityType + "\n" + GlucoseTransportType
                Simulation.SetTitle(Title)
                Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                Simulation.SetGlucoseTransportType(GlucoseTransportType)
                Datasets[Title] = Simulation.Sim(Chemotaxis=1)

    if ModelSwitch == 2:
        # Glycolysis Only
        ModelName = "Glycolysis"
        ATPProductionTypes = ["ATP'=cg", "ATP'=[ADP]/[ATP]*cg"]
        # ATPProductionTypes = ["dATP=[G6P]*cg", "dATP=[G6P]/[ATP]*cg"]
        ATPConsumptionTypes = ["No_ATPConsumption", "Linear_ATPConsumption", "Burst_ATPConsumption"]
        # ATPConsumptionTypes = ["Linear_ATPConsumption", "Burst_ATPConsumption"]

        for ATPProductionType in ATPProductionTypes:
            for ATPConsumptionType in ATPConsumptionTypes:
                Title = ModelName + "\n" + ATPProductionType + "\n" + ATPConsumptionType
                Simulation.SetTitle(Title)
                Simulation.SetATPProductionType(ATPProductionType)
                Simulation.SetATPConsumptionType(ATPConsumptionType)
                # Datasets[Title] = Simulation.Glycolysis()
                Datasets[Title] = Simulation.Sim(Glycolysis=1, ATPConsumption=1)


    if ModelSwitch == 11:
        # Chemotaxis + Glycolysis
        ModelName = "Chemotaxis + Glycolysis"
        GlucoseAvailabilityTypes = ["Uniform_GlucoseAvailability", "Linear_GlucoseAvailability", "NonLinear_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Uniform_GlucoseAvailability", "NonLinear_GlucoseAvailability"]
        # GlucoseAvailabilityTypes = ["No_GlucoseAvailability", "Uniform_GlucoseAvailability", "Linear_GlucoseAvailability", "NonLinear_GlucoseAvailability"]
        GlucoseTransportTypes = ["G6P'=1/[G6P]*cgt"]
        ATPProductionTypes = ["ATP'=[ADP]/[ATP]*cg"]
        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                for ATPProductionType in ATPProductionTypes:
                    Title = ModelName + "\n" + GlucoseAvailabilityType + "\n" + GlucoseTransportType + "\n" + ATPProductionType
                    Simulation.SetTitle(Title)
                    Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                    Simulation.SetGlucoseTransportType(GlucoseTransportType)
                    Simulation.SetATPProductionType(ATPProductionType)
                    Simulation.SetATPConsumptionType("Burst_ATPConsumption")
                    Datasets[Title] = Simulation.Sim(Chemotaxis=1, Glycolysis=1, ATPConsumption=1)


    # Plotting
    Plot = Plotter()

    InclusionList = ["ATP", PerturbationName]
    # InclusionList = ["ATP", "ADP", "pyruvate", PerturbationName]
    # Plot.SetFilter_Inclusion(InclusionList)

    ExclusionList = ["Pi, NADH, NAD", "pyruvate"]
    Plot.SetFilter_Exclusion(ExclusionList)

    Plot.PlotData(Datasets, SimulationTimeUnit)


if __name__ == '__main__':
    main()
