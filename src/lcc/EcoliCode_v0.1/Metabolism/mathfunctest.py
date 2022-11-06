import random
import sys
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import math

PerturbationTag = "#"
G6P_Constant = 8.8e-3

Key_GlucoseTransport = "GlucoseTransport"
Key_GlucoseAvailability = "GlucoseAvailability"
Key_G6PSink = "G6P_Sink"
Key_GlycolysisATP = "GlycolysisATP"
Key_ATPSink = "ATP_Sink"

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
        self.ModelType = dict()
        self.SupportedModelTypes2Int = dict()
        self.InitializeSupportedModelTypes()

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

    def InitializeModelType(self):
        self.ModelType[Key_GlucoseTransport] = 0
        self.ModelType[Key_GlucoseAvailability] = 0
        self.ModelType[Key_G6PSink] = 0
        self.ModelType[Key_GlycolysisATP] = 0
        self.ModelType[Key_ATPSink] = 0

    def InitializeSupportedModelTypes(self):
        # Glucose Availability
        self.SupportedModelTypes2Int["No_GlucoseAvailability"] = 0
        self.SupportedModelTypes2Int["Constant_GlucoseAvailability"] = 1
        self.SupportedModelTypes2Int["Increasing_GlucoseAvailability"] = 2
        self.SupportedModelTypes2Int["Fluctuating_GlucoseAvailability"] = 3

        # Glucose Transport
        self.SupportedModelTypes2Int["G6P'=0"] = 0
        self.SupportedModelTypes2Int["G6P'=cgt"] = 1
        self.SupportedModelTypes2Int["G6P'=0.0001/([G6P]*cgt)^2"] = 2

        # Arbitrary G6P Sink
        self.SupportedModelTypes2Int["No_G6PSink"] = 0
        self.SupportedModelTypes2Int["Linear_G6PSink"] = 1
        self.SupportedModelTypes2Int["Burst_G6PSink"] = 2

        # Glycolysis - ATP Production
        self.SupportedModelTypes2Int["ATP'=0"] = 0
        self.SupportedModelTypes2Int["ATP'=cg"] = 1
        self.SupportedModelTypes2Int["ATP'=[ADP]/[ATP]*cg"] = 2
        self.SupportedModelTypes2Int["dATP=[G6P]*cg"] = 3

        # Arbitrary ATP Sink
        self.SupportedModelTypes2Int["No_ATPSink"] = 0
        self.SupportedModelTypes2Int["Linear_ATPSink"] = 1
        self.SupportedModelTypes2Int["Burst_ATPSink"] = 2

    def GetModelType(self, Type):
        if Type in self.ModelType.keys():
            return self.ModelType[Type]
        else:
            print("[ERROR] Unrecognized Model: %s" % Type)
            sys.exit(1)

    def GetSupportedModelType2Int(self, Type):
        if Type in self.SupportedModelTypes2Int.keys():
            return self.SupportedModelTypes2Int[Type]
        else:
            print("[ERROR] Unsupported Model Type: %s" % Type)
            sys.exit(1)

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

    def SetGlucoseAvailabilityType(self, Type):
        self.ModelType[Key_GlucoseAvailability] = self.GetSupportedModelType2Int(Type)

    def SetGlucoseTransportType(self, Type):
        self.ModelType[Key_GlucoseTransport] = self.GetSupportedModelType2Int(Type)

    def SetG6PSinkType(self, Type):
        self.ModelType[Key_G6PSink] = self.GetSupportedModelType2Int(Type)

    def SetGlycolysisATPType(self, Type):
        self.ModelType[Key_GlycolysisATP] = self.GetSupportedModelType2Int(Type)

    def SetATPSinkType(self, Type):
        self.ModelType[Key_ATPSink] = self.GetSupportedModelType2Int(Type)

    def Sim(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCA=0, OxidativePhosphorylation=0, G6PSink=0, ATPSink=0):
        self.PrintSimulationSubject()
        self.SteadyStateCheckMolecule = "pyruvate"
        self.Initialize(Chemotaxis, Glycolysis, PyruvateOxidation, TCA, OxidativePhosphorylation)
        self.Run(Chemotaxis, Glycolysis, PyruvateOxidation, TCA, OxidativePhosphorylation, G6PSink, ATPSink)
        self.TagPerturbationLabel()
        return self.Dataset.copy()

    def Initialize_Chemotaxis(self):
        # Initial Molecule Concentrations
        self.MolConc["G6P"] = 8.8e-3   # all Hexose-P combined
        self.MolConc["ATP"] = 9.6e-3
        self.MolConc["ADP"] = 5.6e-4

        # Perturbations
        self.Perturbation["glucose"] = 0
        self.Perturbation["glucose_Influx"] = 0
        self.Perturbation["G6P_Sink"] = 0

        # Constants
        self.Constants["cc"] = 0.01
        self.Constants["cgt"] = 0.5

    def GetGlucoseAvailability_Chemotaxis(self):
        Model = self.GetModelType(Key_GlucoseAvailability)
        if Model == 1:
            return 0.15
        elif Model == 2:
            return 0.000003 * self.SimStep
        elif Model == 3:
            return 0.1 * math.sin(self.SimStep/4000) + 1e-7 / self.SimulationTimeUnit
        else:
            return 0

    def GetTransportFactor_Chemotaxis(self):
        Model = self.GetModelType(Key_GlucoseTransport)
        if Model == 1:
            return self.Constants["cgt"]
        elif Model == 2:
            return 0.0001 / (self.Dataset["G6P"][-1] ** 2) * self.Constants["cgt"]
        else:
            return 0

    def Simulation_Chemotaxis(self):
        d = dict()

        # Production Rate Type
        GlucoseAvailability = self.GetGlucoseAvailability_Chemotaxis()
        self.Perturbation["glucose"] = GlucoseAvailability * self.SimulationTimeUnit

        TransportFactor = self.GetTransportFactor_Chemotaxis()
        TransportDemand = TransportFactor if TransportFactor <= 1 else 1   # TODO: make the limit approaching smooth
        self.Perturbation["glucose_Influx"] = self.Perturbation["glucose"] * TransportDemand

        d["G6P"] = self.Perturbation["glucose_Influx"]
        d["ATP"] = - self.Constants["cc"] * self.SimulationTimeUnit
        d["ADP"] = - d["ATP"]   # use when G6P is constant

        self.AddToDeltaMolConc(d)

    def Initialize_Glycolysis(self):
        # G6P + 2 ADP + 2 NAD+ + 2 Pi --> 2 pyruvate + 2 NADH + 2 ATP

        # Initial Molecule Concentrations
        self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["ADP"]      = 5.6e-4
        self.MolConc["NAD"]      = 2.6e-3
        self.MolConc["Pi"]       = 6e-3   # undocumented
        self.MolConc["pyruvate"] = 3.9e-4
        self.MolConc["NADH"]     = 8.3e-3
        self.MolConc["ATP"]      = 9.6e-3

        # Perturbations
        self.Perturbation["ATP_Sink"] = 0

        # Constants
        self.Constants["cg"]     = 10

    def GetInputFactor_Glycolysis(self):
        Model = self.GetModelType(Key_GlycolysisATP)
        if Model == 1:
            return 0.01
        elif Model == 2:
            return self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        elif Model == 3:
            return self.Dataset["G6P"][-1] * 2
        else:
            return 0

    def Simulation_Glycolysis(self):
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_Glycolysis()
        dATP = (InputFactor * self.Constants["cg"]) * self.SimulationTimeUnit

        d["ATP"] = dATP if dATP / 2 < self.Dataset["G6P"][-1] else self.Dataset["G6P"][-1] * 2   # TODO: make the limit approaching smooth
        d["ADP"] = - d["ATP"]
        d["Pi"] = - d["ATP"]
        d["NADH"] = d["ATP"]
        d["NAD"] = - d["NADH"]
        d["G6P"] = - d["ATP"] / 2
        # d["pyruvate"] = - d["G6P"] * 2
        d["pyruvate"] = d["ATP"]   # use when G6P is constant

        self.AddToDeltaMolConc(d)

    def GetG6PSinkRate(self):
        Model = self.GetModelType(Key_G6PSink)
        if Model == 1:   # linear
            return (self.TotalSimulationTime * 1.4)
        elif Model == 2:   # non-linear
            # self.Perturbation[PerturbationName] = (self.TotalSimulationTime * 0.1) * self.SimulationTimeUnit * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
            return (self.TotalSimulationTime * 0.1) * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
        else:
            return 0

    def Simulation_G6PSink(self):
        d = dict()

        self.Perturbation["G6P_Sink"] = self.GetG6PSinkRate() * 0.3 * self.SimulationTimeUnit

        d["G6P"] = - self.Perturbation["G6P_Sink"]

        self.AddToDeltaMolConc(d)

    def GetATPSinkRate(self):
        Model = self.GetModelType(Key_ATPSink)
        if Model == 1:   # linear
            return (self.TotalSimulationTime * 1.4)
        elif Model == 2:   # non-linear
            # self.Perturbation[PerturbationName] = (self.TotalSimulationTime * 0.1) * self.SimulationTimeUnit * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
            return (self.TotalSimulationTime * 0.1) * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
        else:
            return 0

    def Simulation_ATPSink(self):
        d = dict()

        self.Perturbation["ATP_Sink"] = self.GetATPSinkRate() * self.SimulationTimeUnit

        d["ATP"] = - self.Perturbation["ATP_Sink"]
        d["ADP"] = self.Perturbation["ATP_Sink"]
        d["Pi"] = self.Perturbation["ATP_Sink"]
        d["NADH"] = - self.Perturbation["ATP_Sink"]
        d["NAD"] = self.Perturbation["ATP_Sink"]

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


    def Run(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCA=0, OxidativePhosphorylation=0, G6PSink=0, ATPSink=0):
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
            if G6PSink:
                self.Simulation_G6PSink()
            if ATPSink:
                self.Simulation_ATPSink()

            self.UpdateMolConc()
            self.AppendData()
            if self.SteadyStateBrake:
                if (abs(self.Dataset[self.SteadyStateCheckMolecule][-1] - self.Dataset[self.SteadyStateCheckMolecule][-2]) / self.Dataset[self.SteadyStateCheckMolecule][-1] < self.SteadyStateThresholdFactor):
                    break
            self.SimStep += 1

        self.PrintSimTimeAndMolConc()

class FModelRunner():
    def __init__(self, Simulation, Plotter):
        self.Simulation = Simulation
        self.Plot = Plotter
        self.Datasets = None

    def InitializeDatasets(self):
        self.Datasets = dict()

    def RunModel(self, ListOfModels):
        for Model in ListOfModels:
            if Model == 1:
                self.Model_Chemotaxis("No_G6PSink")
            if Model == 2:
                self.Model_Chemotaxis("Linear_G6PSink")
            if Model == 3:
                self.Model_Chemotaxis("Burst_G6PSink")
            if Model == 11:
                self.Model_Glycolysis()
            # if Model == 21:
            #     self.Model_PyruvateOxidation()
            if Model == 101:
                self.Model_Chemotaxis_Glycolysis()

    def Model_Chemotaxis(self, G6PSinkModel):
        self.InitializeDatasets()
        ModelName = "[Chemotaxis] glucose{out} --> G6P{in}"
        # GlucoseAvailabilityTypes = ["No_GlucoseAvailability", "Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Fluctuating_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability"]
        # GlucoseTransportTypes = ["G6P'=0", "G6P'=cgt", "G6P'=0.0001/([G6P]*cgt)^2"]
        GlucoseTransportTypes = []
        # G6PSinkTypes = ["No_G6PSink", "Linear_G6PSink", "Burst_G6PSink"]
        G6PSinkTypes = [G6PSinkModel]

        if G6PSinkModel == "No_G6PSink":
            # Chemotaxis without G6P Sink
            GlucoseTransportTypes = ["G6P'=0", "G6P'=cgt", "G6P'=0.0001/([G6P]*cgt)^2"]

        else:
            GlucoseTransportTypes = ["G6P'=cgt", "G6P'=0.0001/([G6P]*cgt)^2"]

        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                for G6PSinkType in G6PSinkTypes:
                    Subtitle = G6PSinkModel + "\n" + GlucoseAvailabilityType + "\n" + GlucoseTransportType
                    self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
                    self.Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                    self.Simulation.SetGlucoseTransportType(GlucoseTransportType)
                    self.Simulation.SetG6PSinkType(G6PSinkType)
                    self.Datasets[Subtitle] = self.Simulation.Sim(Chemotaxis=1, G6PSink=1)

        InclusionList = ["G6P", "glucose", "glucose_Influx", "G6P_Sink"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotData(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

    def Model_Glycolysis(self):
        self.InitializeDatasets()
        ModelName = "[Glycolysis] G6P + 2 ADP + 2 NAD+ + 2 Pi --> 2 pyruvate + 2 NADH + 2 ATP"
        # GlycolysisATPTypes = ["ATP'=0", "ATP'=cg", "ATP'=[ADP]/[ATP]*cg"]
        GlycolysisATPTypes = ["ATP'=cg", "ATP'=[ADP]/[ATP]*cg"]
        # GlycolysisATPTypes = ["dATP=[G6P]*cg", "dATP=[G6P]/[ATP]*cg"]
        ATPSinkTypes = ["No_ATPSink", "Linear_ATPSink", "Burst_ATPSink"]
        # ATPSinkTypes = ["Linear_ATPSink", "Burst_ATPSink"]

        for GlycolysisATPType in GlycolysisATPTypes:
            for ATPSinkType in ATPSinkTypes:
                Subtitle = GlycolysisATPType + "\n" + ATPSinkType
                self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
                self.Simulation.SetGlycolysisATPType(GlycolysisATPType)
                self.Simulation.SetATPSinkType(ATPSinkType)
                # Datasets[Title] = Simulation.Glycolysis()
                self.Datasets[Subtitle] = self.Simulation.Sim(Glycolysis=1, ATPSink=1)

        InclusionList = ["ATP", "G6P", "pyruvate", "ATP_Sink"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotData(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)


    def Model_Chemotaxis_Glycolysis(self):
        self.InitializeDatasets()
        ModelName = "Chemotaxis + Glycolysis"
        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Burst_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability"]
        GlucoseTransportTypes = ["G6P'=0.0001/([G6P]*cgt)^2"]
        # GlycolysisATPTypes = ["ATP'=cg", "ATP'=[ADP]/[ATP]*cg"]
        GlycolysisATPTypes = ["ATP'=[ADP]/[ATP]*cg"]
        ATPSinkTypes = ["Linear_ATPSink", "Burst_ATPSink"]

        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                for GlycolysisATPType in GlycolysisATPTypes:
                    for ATPSinkType in ATPSinkTypes:
                        Subtitle = GlucoseAvailabilityType + "\n" + GlucoseTransportType + "\n" + GlycolysisATPType + "\n" + ATPSinkType
                        self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
                        self.Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                        self.Simulation.SetGlucoseTransportType(GlucoseTransportType)
                        self.Simulation.SetGlycolysisATPType(GlycolysisATPType)
                        # self.Simulation.SetATPSinkType("Burst_ATPSink")
                        self.Simulation.SetATPSinkType(ATPSinkType)
                        self.Datasets[Subtitle] = self.Simulation.Sim(Chemotaxis=1, Glycolysis=1, ATPSink=1)

        PlottingPathway = "(Show all pathways)"
        InclusionList = ["ATP", "ATP_Sink", "glucose", "glucose_Influx"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotData(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + "\n" + PlottingPathway))

        PlottingPathway = "(Show Chemotaxis Only)"
        InclusionList = ["G6P", "glucose", "glucose_Influx"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotData(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + "\n" + PlottingPathway))

        PlottingPathway = "(Show Glycolysis Only)"
        InclusionList = ["ATP", "G6P", "pyruvate", "ATP_Sink"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotData(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + "\n" + PlottingPathway))


class FPlotter:
    def __init__(self):
        self.Filter_Inclusion = None
        self.Filter_Exclusion = None

    def SetFilters(self, InclusionList, ExclusionList):
        self.SetFilter_Inclusion(InclusionList)
        self.SetFilter_Exclusion(ExclusionList)

    def SetFilter_Inclusion(self, List):
        self.Filter_Inclusion = List

    def SetFilter_Exclusion(self, List):
        self.Filter_Exclusion = List

    def CheckToIncludeOrExclude(self, Key_Data):
        if self.Filter_Inclusion or self.Filter_Exclusion:
            if self.Filter_Inclusion:
                if (Key_Data in self.Filter_Inclusion) or (Key_Data[1:] in self.Filter_Inclusion):
                    return True
                else:
                    return False
            else:
                if (Key_Data in self.Filter_Exclusion) or (Key_Data[1:] in self.Filter_Exclusion):
                    return False
                else:
                    return True
        else:
            return True

    def FilterDatasets(self, Datasets):
        Datasets_Filtered = dict()
        for Key_Dataset, Dataset in Datasets.items():
            Dataset_Filtered = dict()
            for Key_Data, Data in Dataset.items():
                if self.CheckToIncludeOrExclude(Key_Data):
                    Dataset_Filtered[Key_Data] = Data
                    Datasets_Filtered[Key_Dataset] = Dataset_Filtered

        return Datasets_Filtered

    def PlotData(self, Datasets, SimulationTimeUnit, bSideLabel=True, SuperTitle=""):

        # Filter Datasets
        if self.Filter_Inclusion or self.Filter_Exclusion:
            Datasets = self.FilterDatasets(Datasets)

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.5, hspace=0.5)
        if SuperTitle:
            fig.suptitle(SuperTitle, fontsize=14)

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
                    line, = ax1.plot(Time, Conc, label="[" + MolName + "]")
                    if bSideLabel:
                        SelectedTimeFrameFromLeft = 0.1
                        ax1.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="left", va="bottom", color=line.get_color())
                        # ax1.text(Time[-1] * 1.01, Conc[-1], MolName, ha="left", va="bottom", color=line.get_color())
                        # ax1.text(Time[-1] * 1.1, Conc[-1], MolName + ": {}".format(Conc[-1]), va="center", color=line.get_color())

                else:
                    line, = ax2.plot(Time, Conc, color=PerturbationPlotColor[PerturbationIndex], label="[" + MolName[1:] + "]")
                    if bSideLabel:
                        SelectedTimeFrameFromLeft = 0.8
                        ax2.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName[1:], ha="center", va="bottom", color=line.get_color())
                    # ax2.plot(Time, Conc, label="[" + MolName[2:] + "]")
                    # print(PerturbationIndex)
                    PerturbationIndex += 1

            ax1.set_title(Process)
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Molecules: Concentration (mol L-1)')
            ax1.set_ylim([0, 0.015])
            if not bSideLabel:
                ax1.legend(loc='upper left')
            if Perturbation:
                ax2.set_ylabel('Event: Concentration (mol L-1)')
                ax2.set_ylim([0, 3e-7])
                if not bSideLabel:
                    ax2.legend(loc='upper right')

            # ax1.grid()

        plt.show()


def main():
    # Initialization
    Simulation = FSimulation()
    Plot = FPlotter()

    # Setting Simulation Parameters
    TotalSimulationTime = 1e-1   # s
    SimulationTimeUnit = 1e-6   # s
    SteadyStateBrake = False
    SteadyStateThresholdFactor = 0.0000001  # Steady state threshold

    Simulation.SetSimulationParameters(TotalSimulationTime, SimulationTimeUnit, SteadyStateBrake, SteadyStateThresholdFactor)
    Simulation.PrintSimulationParameters()

    # Model Runner
    ModelRunner = FModelRunner(Simulation, Plot)

    # Select Models
    # 1:    Chemotaxis unit test without G6PSink
    # 2:    Chemotaxis unit test with Linear G6PSink
    # 3:    Chemotaxis unit test with Burst G6PSink
    # 11:   Glycolysis unit test with ATPSink
    # 21:   Pyruvate Oxidation unit test
    # 101:  Chemotaxis + Glycolysis

    Models = [1]    # Chemotaxis unit test without G6PSink
    Models = [2]    # Chemotaxis unit test with Linear G6PSink
    Models = [3]    # Chemotaxis unit test with Burst G6PSink
    Models = [11]   # Glycolysis unit test with ATPSink
    Models = [21]   # Pyruvate Oxidation unit test
    Models = [101]  # Chemotaxis + Glycolysis
    # Models = [1, 2, 3, 11]
    ModelRunner.RunModel(Models)

if __name__ == '__main__':
    main()
