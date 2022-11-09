import random
import sys
import math

from src.lcc.plot2 import FPlotter

class FDiffusion():
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

    def SetSimulationParameters(self, TotalSimulationTime, SimulationTimeUnit, SteadyStateBrake=False, SteadyStateCheckMolecule=None, SteadyStateThresholdFactor=None):
        self.TotalSimulationTime = TotalSimulationTime
        self.SimulationTimeUnit = SimulationTimeUnit
        self.SteadyStateBrake = SteadyStateBrake
        self.SteadyStateCheckMolecule = SteadyStateCheckMolecule
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
        self.ModelType[Key_GlycolysisEndProduct] = 0

    def InitializeSupportedModelTypes(self):
        # Glucose Availability
        self.SupportedModelTypes2Int["No_GlucoseAvailability"] = 0
        self.SupportedModelTypes2Int["Constant_GlucoseAvailability"] = 1
        self.SupportedModelTypes2Int["Increasing_GlucoseAvailability"] = 2
        self.SupportedModelTypes2Int["Fluctuating_GlucoseAvailability"] = 3

        # Glucose Transport
        self.SupportedModelTypes2Int["G6P'=0"] = 0
        self.SupportedModelTypes2Int["G6P'=cgt"] = 1
        self.SupportedModelTypes2Int["G6P'=0.0001/(G6P*cgt)^2"] = 2

        # Arbitrary G6P Sink
        self.SupportedModelTypes2Int["No_G6PSink"] = 0
        self.SupportedModelTypes2Int["Linear_G6PSink"] = 1
        self.SupportedModelTypes2Int["Burst_G6PSink"] = 2

        # Glycolysis - ATP Production
        self.SupportedModelTypes2Int["ATP'=0"] = 0
        self.SupportedModelTypes2Int["ATP'=cg"] = 1
        self.SupportedModelTypes2Int["ATP'=ADP/ATP*cg"] = 2
        self.SupportedModelTypes2Int["ATP=G6P*cg"] = 3

        # Glycolysis: End product
        self.SupportedModelTypes2Int["pyruvate"] = 0
        self.SupportedModelTypes2Int["acetyl-CoA"] = 1

        # Pyruvate Oxidation - ATP Production
        self.SupportedModelTypes2Int["pyruvate'=0"] = 0
        self.SupportedModelTypes2Int["pyruvate'=cp"] = 1
        self.SupportedModelTypes2Int["pyruvate'=ADP/ATP*cp"] = 2   # Not testable without ATP and ADP
        self.SupportedModelTypes2Int["pyruvate'=pyruvate/acetyl-CoA*cp"] = 3

        # TCA Cycle: ATP Production
        # self.SupportedModelTypes2Int["ATP'=0"] = 0   # same as in Glycolysis ATP production
        self.SupportedModelTypes2Int["ATP'=ct"] = 1
        self.SupportedModelTypes2Int["ATP'=ADP/ATP*ct"] = 2
        # self.SupportedModelTypes2Int["ATP=G6P*cg"] = 3

        # TCA Cycle: ATP Production
        self.SupportedModelTypes2Int["Unfixed"] = 0
        self.SupportedModelTypes2Int["Fixed"] = 1

        # Oxidative Phosphorylation: ATP Production
        # self.SupportedModelTypes2Int["ATP'=0"] = 0   # same as in Glycolysis ATP production
        self.SupportedModelTypes2Int["ATP'=co"] = 1
        self.SupportedModelTypes2Int["ATP'=ADP/ATP*co"] = 2

        # Arbitrary ATP Sink
        self.SupportedModelTypes2Int["No_ATPSink"] = 0
        self.SupportedModelTypes2Int["Constant_ATPSink"] = 1
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

    def CorrectAndTagPerturbationLabel(self):
        for PerturbationName, Conc in self.Perturbation.items():
            self.Dataset[PerturbationName][0] = self.Dataset[PerturbationName][1]
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

    def SetPyruvateOxidationATPType(self, Type):
        self.ModelType[Key_PyruvateOxidationATP] = self.GetSupportedModelType2Int(Type)

    def SetAcetylCoALevel(self, Type):
        self.ModelType[Key_AcetylCoALevel] = Type

    def SetTCACycleATPType(self, Type):
        self.ModelType[Key_TCACycleATP] = self.GetSupportedModelType2Int(Type)

    def SetOxidativePhosphorylationATPType(self, Type):
        self.ModelType[Key_OxidativePhosphorylationATP] = self.GetSupportedModelType2Int(Type)

    def SetATPSinkType(self, Type):
        self.ModelType[Key_ATPSink] = self.GetSupportedModelType2Int(Type)

    def SetGlycolysisEndProduct(self, EndProduct):
        self.ModelType[Key_GlycolysisEndProduct] = self.GetSupportedModelType2Int(EndProduct)

    def Sim(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCACycle=0, OxidativePhosphorylation=0, G6PSink=0, ATPSink=0):
        self.PrintSimulationSubject()
        self.Initialize(Chemotaxis, Glycolysis, PyruvateOxidation, TCACycle, OxidativePhosphorylation)
        self.Run(Chemotaxis, Glycolysis, PyruvateOxidation, TCACycle, OxidativePhosphorylation, G6PSink, ATPSink)
        self.CorrectAndTagPerturbationLabel()
        return self.Dataset.copy()

    def Initialize_Chemotaxis(self):
        # Initial Molecule Concentrations
        self.MolConc["G6P"] = 8.8e-3   # all Hexose-P combined
        self.MolConc["ATP"] = 9.6e-3
        self.MolConc["ADP"] = 5.6e-4

        # Perturbations
        self.Perturbation["glucose{out}"] = 0
        self.Perturbation["glucose_Influx"] = 0
        self.Perturbation["G6P_Sink"] = 0

        # Constants
        self.Constants["cc"] = 0.01
        self.Constants["cgt"] = 0.5

    def GetGlucoseAvailability_Chemotaxis(self):
        Model = self.GetModelType(Key_GlucoseAvailability)
        if Model == 0:
            return 0
        elif Model == 1:
            return 0.15
        elif Model == 2:
            return 0.000003 * self.SimStep
        elif Model == 3:
            return 0.1 * math.sin(self.SimStep/4000) + 1e-7 / self.SimulationTimeUnit
        else:
            sys.exit(1)

    def GetTransportFactor_Chemotaxis(self):
        Model = self.GetModelType(Key_GlucoseTransport)
        if Model == 0:
            return 0
        elif Model == 1:
            return self.Constants["cgt"]
        elif Model == 2:
            return 0.0001 / (self.Dataset["G6P"][-1] ** 2) * self.Constants["cgt"]
        else:
            sys.exit(1)

    def Simulation_Chemotaxis(self):
        d = dict()

        # Production Rate Type
        GlucoseAvailability = self.GetGlucoseAvailability_Chemotaxis()
        self.Perturbation["glucose{out}"] = GlucoseAvailability * self.SimulationTimeUnit

        TransportFactor = self.GetTransportFactor_Chemotaxis()
        TransportDemand = TransportFactor if TransportFactor <= 1 else 1   # TODO: make the limit approaching smooth
        self.Perturbation["glucose_Influx"] = self.Perturbation["glucose{out}"] * TransportDemand

        d["G6P"] = self.Perturbation["glucose_Influx"]
        d["ATP"] = - self.Constants["cc"] * self.SimulationTimeUnit
        d["ADP"] = - d["ATP"]   # use when G6P is constant

        self.AddToDeltaMolConc(d)

    def Initialize_Glycolysis(self):
        # Without pyruvate oxidation
        # G6P + 2 ADP + 2 NAD+ + 2 Pi --> 2 pyruvate + 2 NADH + 2
        # With pyruvate oxidation
        # G6P + 2 ADP + 4 NAD+ + 2 Pi --> 2 acetyl-CoA + 4 NADH + 2 ATP

        # Initial Molecule Concentrations
        self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["ADP"]      = 5.6e-4
        self.MolConc["NAD"]      = 2.6e-3
        # self.MolConc["Pi"]       = 6e-3   # undocumented
        self.MolConc["NADH"]     = 8.3e-3
        self.MolConc["ATP"]      = 9.6e-3
        if self.ModelType[Key_GlycolysisEndProduct] == 0:
            self.MolConc["pyruvate"] = 3.9e-4
        elif self.ModelType[Key_GlycolysisEndProduct] == 1:
            self.MolConc["acetyl-CoA"] = 6.1e-4

        # Perturbations
        self.Perturbation["ATP_Sink"] = 0

        # Constants
        self.Constants["cg"]     = 10

    def GetInputFactor_Glycolysis(self):
        Model = self.GetModelType(Key_GlycolysisATP)
        if Model == 0:
            return 0
        elif Model == 1:
            return 0.01
        elif Model == 2:
            return self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        elif Model == 3:
            return self.Dataset["G6P"][-1] * 2
        else:
            sys.exit(1)

    def Simulation_Glycolysis(self):
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_Glycolysis()
        dATP = (InputFactor * self.Constants["cg"]) * self.SimulationTimeUnit

        d["ATP"] = dATP
        d["ADP"] = - d["ATP"]
        # d["Pi"] = - d["ATP"]
        d["G6P"] = - d["ATP"] / 2
        # d["pyruvate"] = - d["G6P"] * 2
        if self.ModelType[Key_GlycolysisEndProduct] == 0:
            d["pyruvate"] = d["ATP"]   # use when G6P is constant
            d["NADH"] = d["ATP"]
            d["NAD"] = - d["ATP"]

        elif self.ModelType[Key_GlycolysisEndProduct] == 1:
            d["acetyl-CoA"] = d["ATP"]   # use when G6P is constant
            d["NADH"] = d["ATP"] * 2
            d["NAD"] = - d["ATP"] * 2

        self.AddToDeltaMolConc(d)

    def GetG6PSinkRate(self):
        Model = self.GetModelType(Key_G6PSink)
        if Model == 0:
            return 0
        elif Model == 1:   # linear
            return (self.TotalSimulationTime * 1.4)
        elif Model == 2:   # non-linear
            # self.Perturbation[PerturbationName] = (self.TotalSimulationTime * 0.1) * self.SimulationTimeUnit * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
            return (self.TotalSimulationTime * 0.1) * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
        else:
            sys.exit(1)


    def Simulation_G6PSink(self):
        d = dict()

        self.Perturbation["G6P_Sink"] = self.GetG6PSinkRate() * 0.3 * self.SimulationTimeUnit

        d["G6P"] = - self.Perturbation["G6P_Sink"]

        self.AddToDeltaMolConc(d)

    def GetATPSinkRate(self):
        Model = self.GetModelType(Key_ATPSink)
        if Model == 0:
            return 0
        elif Model == 1:   # linear
            return (self.TotalSimulationTime * 1.4)
        elif Model == 2:   # non-linear
            # self.Perturbation[PerturbationName] = (self.TotalSimulationTime * 0.1) * self.SimulationTimeUnit * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
            return (self.TotalSimulationTime * 0.1) * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
        else:
            sys.exit(1)

    def Simulation_ATPSink(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCACycle=0, OxidativePhosphorylation=0):
        d = dict()

        dATP = self.GetATPSinkRate() * self.SimulationTimeUnit

        # Combinatorial case
        if Chemotaxis and Glycolysis and PyruvateOxidation and TCACycle:
            d["NADH"] = - dATP * 0.25
            d["NAD"] = dATP * 0.25
            d["FADH2"] = - dATP * 0.05
            d["FAD"] = dATP * 0.05
        if Chemotaxis and Glycolysis and PyruvateOxidation:
            pass
        elif Chemotaxis and Glycolysis:
            if self.ModelType[Key_GlycolysisEndProduct] == 0:
                d["NADH"] = - dATP
                d["NAD"] = dATP
            elif self.ModelType[Key_GlycolysisEndProduct] == 1:
                d["NADH"] = - dATP * 2
                d["NAD"] = dATP * 2
        # Single case
        elif Chemotaxis:
            pass
        elif Glycolysis:
            if self.ModelType[Key_GlycolysisEndProduct] == 0:
                d["NADH"] = - dATP
                d["NAD"] = dATP
            elif self.ModelType[Key_GlycolysisEndProduct] == 1:
                d["NADH"] = - dATP * 2
                d["NAD"] = dATP * 2
        elif TCACycle:
            d["NADH"] = - dATP * 0.5
            d["NAD"] = dATP * 0.5
            d["FADH2"] = - dATP * 0.1
            d["FAD"] = dATP * 0.1
        elif OxidativePhosphorylation:
            pass

        self.Perturbation["ATP_Sink"] = dATP
        d["ATP"] = - dATP
        d["ADP"] = dATP
        # d["Pi"] = dATP

        self.AddToDeltaMolConc(d)

    def Initialize_PyruvateOxidation(self):
        # With pyruvate oxidation
        # pyruvate + CoA-SH + NAD --> acetyl-CoA + NADH

        # Initial Molecule Concentrations
        self.MolConc["pyruvate"] = 3.9e-4
        self.MolConc["acetyl-CoA"] = 6.1e-4
        self.MolConc["NAD"]      = 2.6e-3
        self.MolConc["NADH"]     = 8.3e-3

        # Perturbations
        self.Perturbation["ATP_Sink"] = 0

        # Constants
        self.Constants["cp"]     = 1

    def GetInputFactor_PyruvateOxidation(self):
        Model = self.GetModelType(Key_PyruvateOxidationATP)
        if Model == 0:
            return 0
        elif Model == 1:
            return 0.01
        elif Model == 2:
            return self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        elif Model == 3:
            return self.Dataset["pyruvate"][-1] / self.Dataset["acetyl-CoA"][-1]
        else:
            sys.exit(1)

    def Simulation_PyruvateOxidation(self):
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_PyruvateOxidation()
        dpyruvate = (InputFactor * self.Constants["cp"]) * self.SimulationTimeUnit

        d["pyruvate"] = - dpyruvate
        d["acetyl-CoA"] = dpyruvate
        d["NADH"] = - dpyruvate
        d["NAD"] = dpyruvate

        self.AddToDeltaMolConc(d)

    def Initialize_TCACycle(self):
        # acetyl-CoA + 3 NAD + FAD + ADP + Pi + 2 H2O --> 2 CO2 + 3 NADH + FADH2 + ATP + 2 H + CoA

        # Initial Molecule Concentrations
        # self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["acetyl-CoA"]  = 6.1e-4
        self.MolConc["FAD"]     = 1.7e-4
        self.MolConc["NAD"]     = 2.6e-3
        self.MolConc["ADP"]     = 5.6e-4
        # self.MolConc["Pi"]      = 6e-3   # undocumented
        self.MolConc["NADH"]    = 8.3e-3
        self.MolConc["FADH2"]   = 5e-4   # undocumented
        self.MolConc["ATP"]     = 9.6e-3
        self.MolConc["CoA"]     = 1.4e-3

        # Perturbations
        self.Perturbation["ATP_Sink"] = 0

        # Constants
        self.Constants["ct"]     = 0.3

    def GetInputFactor_TCACycle(self):
        Model = self.GetModelType(Key_TCACycleATP)
        if Model == 0:
            return 0
        elif Model == 1:
            return 0.01
        elif Model == 2:
            return self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        elif Model == 3:
            return self.Dataset["G6P"][-1] * 2
        else:
            sys.exit(1)

    def Simulation_TCACycle(self):
        # acetyl-CoA + 3 NAD + FAD + ADP + Pi + 2 H2O --> 2 CO2 + 3 NADH + FADH2 + ATP + 2 H + CoA
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_TCACycle()
        dATP = (InputFactor * self.Constants["ct"]) * self.SimulationTimeUnit

        d["ATP"] = dATP
        d["ADP"] = - d["ATP"]
        # d["Pi"] = -d["ATP"]

        d["NADH"] = d["ATP"] * 3
        d["NAD"] = - d["NADH"]
        d["FADH2"] = d["ATP"]
        d["FAD"] = - d["FADH2"]

        if self.GetModelType(Key_AcetylCoALevel) == 0:
            d["acetyl-CoA"] = - d["ATP"]
            d["CoA"] = d["ATP"]

        self.AddToDeltaMolConc(d)

    def Initialize_OxidativePhosphorylation(self):
        # NADH + FADH2 + ADP + Pi --> NAD + FAD + 4 ATP

        # Initial Molecule Concentrations
        self.MolConc["FAD"]     = 1.7e-4
        self.MolConc["NAD"]     = 2.6e-3
        self.MolConc["ADP"]     = 5.6e-4
        # self.MolConc["Pi"]      = 6e-3   # undocumented
        self.MolConc["NADH"]    = 8.3e-3
        self.MolConc["FADH2"]   = 5e-4   # undocumented
        self.MolConc["ATP"]     = 9.6e-3

        # Perturbations
        self.Perturbation["ATP_Sink"] = 0

        # Constants
        self.Constants["co"]     = 1

    def GetInputFactor_OxidativePhosphorylation(self):
        Model = self.GetModelType(Key_OxidativePhosphorylationATP)
        if Model == 0:
            return 0
        elif Model == 1:
            return 0.01
        elif Model == 2:
            return self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        else:
            sys.exit(1)

    def Simulation_OxidativePhosphorylation(self):
        # NADH + FADH2 + 4 ADP + Pi --> NAD + FAD + 4 ATP
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_OxidativePhosphorylation()
        dATP = (InputFactor * self.Constants["co"]) * self.SimulationTimeUnit

        d["ATP"] = dATP
        d["ADP"] = - d["ATP"]

        d["NADH"] = d["ATP"] / 4
        d["NAD"] = - d["NADH"]
        d["FADH2"] = d["ATP"] /4
        d["FAD"] = - d["FADH2"]

        self.AddToDeltaMolConc(d)


    def Initialize(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCACycle=0, OxidativePhosphorylation=0):
        if Chemotaxis:
            self.Initialize_Chemotaxis()
        if Glycolysis:
            self.Initialize_Glycolysis()
        if PyruvateOxidation:
            self.Initialize_PyruvateOxidation()
        if TCACycle:
            self.Initialize_TCACycle()
        if OxidativePhosphorylation:
            self.Initialize_OxidativePhosphorylation()
        self.InitializeDeltaMolConc()
        self.InitializeDataset()
        self.InitializeSimStep()


    def Run(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCACycle=0, OxidativePhosphorylation=0, G6PSink=0, ATPSink=0):
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
            if TCACycle:
                self.Simulation_TCACycle()
            if OxidativePhosphorylation:
                self.Simulation_OxidativePhosphorylation()
            if G6PSink:
                self.Simulation_G6PSink()
            if ATPSink:
                self.Simulation_ATPSink(Chemotaxis, Glycolysis, PyruvateOxidation, TCACycle, OxidativePhosphorylation)

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
                self.Model_Chemotaxis(G6PSinkModel="No_G6PSink")
            elif Model == 2:
                self.Model_Chemotaxis(G6PSinkModel="Linear_G6PSink")
            elif Model == 3:
                self.Model_Chemotaxis(G6PSinkModel="Burst_G6PSink")
            elif Model == 11:
                self.Model_Glycolysis(EndProduct="pyruvate")
            elif Model == 12:
                self.Model_Glycolysis(EndProduct="acetyl-CoA")
            elif Model == 21:
                self.Model_PyruvateOxidation()
            elif Model == 31:
                self.Model_TCACycle(AcetylCoALevel="Unfixed")
            elif Model == 32:
                self.Model_TCACycle(AcetylCoALevel="Fixed")
            elif Model == 41:
                self.Model_OxidativePhosphorylation()
            elif Model == 101:
                self.Model_Chemotaxis_Glycolysis(EndProduct="pyruvate")
            elif Model == 102:
                self.Model_Chemotaxis_Glycolysis(EndProduct="acetyl-CoA")
            elif Model == 103:
                self.Model_Chemotaxis_Glycolysis_PyruvateOxidation()
            elif Model == 104:
                self.Model_Chemotaxis_Glycolysis_PyruvateOxidation_TCACycle()
            elif Model == 105:
                self.Model_Chemotaxis_Glycolysis_PyruvateOxidation_TCACycle_OxidativePhosphorylation()
            else:
                print("[ERROR] Unsupported Model: %s" % Model)
                sys.exit(1)


    def Model_Chemotaxis(self, G6PSinkModel="Linear_G6PSink"):
        self.InitializeDatasets()
        ModelName = "[Chemotaxis] glucose{out} --> G6P{in}"

        # GlucoseAvailabilityTypes = ["No_GlucoseAvailability", "Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Fluctuating_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability"]

        # GlucoseTransportTypes = ["G6P'=0", "G6P'=cgt", "G6P'=0.0001/(G6P*cgt)^2"]
        GlucoseTransportTypes = []

        # G6PSinkTypes = ["No_G6PSink", "Linear_G6PSink", "Burst_G6PSink"]
        G6PSinkTypes = [G6PSinkModel]

        if G6PSinkModel == "No_G6PSink":
            # Chemotaxis without G6P Sink
            GlucoseTransportTypes = ["G6P'=0", "G6P'=cgt", "G6P'=0.0001/(G6P*cgt)^2"]

        else:
            GlucoseTransportTypes = ["G6P'=cgt", "G6P'=0.0001/(G6P*cgt)^2"]

        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                for G6PSinkType in G6PSinkTypes:
                    Subtitle = G6PSinkModel + "\n" + GlucoseAvailabilityType + "\n" + GlucoseTransportType
                    self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
                    self.Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                    self.Simulation.SetGlucoseTransportType(GlucoseTransportType)
                    self.Simulation.SetG6PSinkType(G6PSinkType)
                    self.Datasets[Subtitle] = self.Simulation.Sim(Chemotaxis=1, G6PSink=1)

        InclusionList = ["G6P", "glucose{out}", "glucose_Influx", "G6P_Sink"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

    def Model_Glycolysis(self, EndProduct="pyruvate"):
        self.InitializeDatasets()
        ModelName = None
        if EndProduct == "pyruvate":
            ModelName = "[Glycolysis] G6P + 2 ADP + 2 NAD+ + 2 Pi --> 2 pyruvate + 2 NADH + 2 ATP"
        elif EndProduct =="acetyl-CoA":
            ModelName = "[Glycolysis] G6P + 2 ADP + 3 NAD+ + 2 Pi --> 2 acetyl-CoA + 4 NADH + 2 ATP"
        else:
            print("[ERROR] Unsupported End Product: %s" % EndProduct)
            sys.exit(1)

        # GlycolysisATPTypes = ["ATP'=0", "ATP'=cg", "ATP'=ADP/ATP*cg"]
        GlycolysisATPTypes = ["ATP'=cg", "ATP'=ADP/ATP*cg"]
        # GlycolysisATPTypes = ["ATP=G6P*cg", "ATP=G6P/[ATP]*cg"]

        ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]
        # ATPSinkTypes = ["Constant_ATPSink", "Burst_ATPSink"]

        for GlycolysisATPType in GlycolysisATPTypes:
            for ATPSinkType in ATPSinkTypes:
                Subtitle = GlycolysisATPType + "\n" + ATPSinkType
                self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
                self.Simulation.SetGlycolysisATPType(GlycolysisATPType)
                self.Simulation.SetATPSinkType(ATPSinkType)
                self.Simulation.SetGlycolysisEndProduct(EndProduct)
                # Datasets[Title] = Simulation.Glycolysis()
                self.Datasets[Subtitle] = self.Simulation.Sim(Glycolysis=1, ATPSink=1)

        InclusionList = ["ATP", "ATP_Sink"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

        InclusionList = ["ATP", "G6P", EndProduct, "ATP_Sink"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + "\n" + PlottingPathway))


    def Model_PyruvateOxidation(self):
        self.InitializeDatasets()
        ModelName = "[PyruvateOxidation] pyruvate + CoA-SH --> acetyl-CoA + NADH"

        # PyruvateOxidationATPTypes = ["pyruvate'=0", "pyruvate'=cp", "pyruvate'=ADP/ATP*cp", "pyruvate'=pyruvate/acetyl-CoA*cp"]
        PyruvateOxidationATPTypes = ["pyruvate'=0", "pyruvate'=cp", "pyruvate'=pyruvate/acetyl-CoA*cp"]
        for PyruvateOxidationATPType in PyruvateOxidationATPTypes:
            Subtitle = PyruvateOxidationATPType
            self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
            self.Simulation.SetPyruvateOxidationATPType(PyruvateOxidationATPType)
            self.Datasets[Subtitle] = self.Simulation.Sim(PyruvateOxidation=1)

        InclusionList = ["NADH"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + "\n" + PlottingPathway))

    def Model_TCACycle(self, AcetylCoALevel="Unfixed"):
        self.InitializeDatasets()
        ModelName = "[TCACycle] acetyl-CoA + 3 NAD + FAD + ADP + Pi --> 3 NADH + FADH2 + ATP + CoA"   # Removed H2O and CO2 and H
        # "[TCACycle] acetyl-CoA + 3 NAD + FAD + ADP + Pi + 2 H2O --> 2 CO2 + 3 NADH + FADH2 + ATP + 2 H + CoA"

        # TCACycleATPTypes = ["ATP'=0", "ATP'=ct", "ATP'=ADP/ATP*ct"]
        TCACycleATPTypes = ["ATP'=ct", "ATP'=ADP/ATP*ct"]

        # ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]
        ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]

        for TCACycleATPType in TCACycleATPTypes:
            for ATPSinkType in ATPSinkTypes:
                Subtitle = TCACycleATPType + "\n" + ATPSinkType
                self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
                self.Simulation.SetTCACycleATPType(TCACycleATPType)
                self.Simulation.SetATPSinkType(ATPSinkType)
                self.Simulation.SetAcetylCoALevel(AcetylCoALevel)
                self.Datasets[Subtitle] = self.Simulation.Sim(TCACycle=1, ATPSink=1)

        InclusionList = ["ATP", "ATP_Sink"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

        InclusionList = ["ATP", "acetyl-CoA", "CoA", "NADH", "FADH2", "ATP_Sink"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + "\n" + PlottingPathway))

    def Model_OxidativePhosphorylation(self):
        self.InitializeDatasets()
        ModelName = "[OxidativePhosphorylation] NADH + FADH2 + ADP + Pi --> NAD + FAD + 4 ATP"   # Not including H other biproducts

        # TCACycleATPTypes = ["ATP'=0", "ATP'=co", "ATP'=ADP/ATP*co"]
        OxidativePhosphorylationATPTypes = ["ATP'=co", "ATP'=ADP/ATP*co"]

        # ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]
        ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]

        for OxidativePhosphorylationATPType in OxidativePhosphorylationATPTypes:
            for ATPSinkType in ATPSinkTypes:
                Subtitle = OxidativePhosphorylationATPType + "\n" + ATPSinkType
                self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
                self.Simulation.SetOxidativePhosphorylationATPType(OxidativePhosphorylationATPType)
                self.Simulation.SetATPSinkType(ATPSinkType)
                self.Datasets[Subtitle] = self.Simulation.Sim(OxidativePhosphorylation=1, ATPSink=1)

        InclusionList = ["ATP", "ATP_Sink"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + "\n" + PlottingPathway))


    def Model_Chemotaxis_Glycolysis(self, EndProduct="pyruvate"):
        self.InitializeDatasets()
        ModelName = "Chemotaxis + Glycolysis"

        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Burst_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability"]

        GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        # GlycolysisATPTypes = ["ATP'=cg", "ATP'=ADP/ATP*cg"]
        GlycolysisATPTypes = ["ATP'=ADP/ATP*cg"]

        ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]

        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                for GlycolysisATPType in GlycolysisATPTypes:
                    for ATPSinkType in ATPSinkTypes:
                        Subtitle = GlucoseAvailabilityType + "\n" + GlucoseTransportType + "\n" + GlycolysisATPType + ATPSinkType + "\n"
                        self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
                        self.Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                        self.Simulation.SetGlucoseTransportType(GlucoseTransportType)
                        self.Simulation.SetGlycolysisATPType(GlycolysisATPType)
                        self.Simulation.SetATPSinkType(ATPSinkType)
                        self.Simulation.SetGlycolysisEndProduct(EndProduct)
                        self.Datasets[Subtitle] = self.Simulation.Sim(Chemotaxis=1, Glycolysis=1, ATPSink=1)

        # spacing = "\n"
        spacing = ", "

        PlottingPathway = "(Show all pathways)"
        InclusionList = ["ATP", "ATP_Sink", "glucose{out}", "glucose_Influx"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show Chemotaxis Only)"
        InclusionList = ["G6P", "glucose{out}", "glucose_Influx"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show Glycolysis Only)"
        InclusionList = ["ATP", "G6P", EndProduct, "ATP_Sink"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))


    def Model_Chemotaxis_Glycolysis_PyruvateOxidation(self, EndProduct="pyruvate"):
        self.InitializeDatasets()
        ModelName = "Chemotaxis + Glycolysis + Pyruvate Oxidation"

        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Burst_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability"]

        GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        # GlycolysisATPTypes = ["ATP'=cg", "ATP'=ADP/ATP*cg"]
        GlycolysisATPTypes = ["ATP'=ADP/ATP*cg"]

        # PyruvateOxidationATPTypes = ["pyruvate'=0", "pyruvate'=cp", "pyruvate'=pyruvate/acetyl-CoA*cp"]
        PyruvateOxidationATPTypes = ["pyruvate'=pyruvate/acetyl-CoA*cp"]

        ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]
        # ATPSinkTypes = ["No_ATPSink"]
        # ATPSinkTypes = ["Constant_ATPSink"]
        # ATPSinkTypes = ["Burst_ATPSink"]

        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                for GlycolysisATPType in GlycolysisATPTypes:
                    for PyruvateOxidationATPType in PyruvateOxidationATPTypes:
                        for ATPSinkType in ATPSinkTypes:
                            Subtitle = GlucoseAvailabilityType + ", " + GlucoseTransportType + "\n" + GlycolysisATPType + ", " + PyruvateOxidationATPType + "\n" + ATPSinkType
                            self.Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                            self.Simulation.SetGlucoseTransportType(GlucoseTransportType)
                            self.Simulation.SetGlycolysisATPType(GlycolysisATPType)
                            # self.Simulation.SetATPSinkType("Burst_ATPSink")
                            self.Simulation.SetPyruvateOxidationATPType(PyruvateOxidationATPType)
                            self.Simulation.SetATPSinkType(ATPSinkType)
                            self.Simulation.SetGlycolysisEndProduct(EndProduct)
                            self.Datasets[Subtitle] = self.Simulation.Sim(Chemotaxis=1, Glycolysis=1, PyruvateOxidation=1, ATPSink=1)

        # spacing = "\n"
        spacing = ", "

        PlottingPathway = "(Show all pathways)"
        InclusionList = ["ATP", "ATP_Sink", "glucose{out}", "glucose_Influx"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show Chemotaxis Only)"
        InclusionList = ["G6P", "glucose{out}", "glucose_Influx"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show Glycolysis Only)"
        InclusionList = ["ATP", "G6P", "pyruvate", "acetyl-CoA", "ATP_Sink"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))


    def Model_Chemotaxis_Glycolysis_PyruvateOxidation_TCACycle(self, EndProduct="pyruvate", AcetylCoALevel="Unfixed"):
        self.InitializeDatasets()
        ModelName = "Chemotaxis + Glycolysis + Pyruvate Oxidation + TCA Cycle"

        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Burst_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Increasing_GlucoseAvailability"]

        GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        # GlycolysisATPTypes = ["ATP'=cg", "ATP'=ADP/ATP*cg"]
        GlycolysisATPTypes = ["ATP'=ADP/ATP*cg"]

        # PyruvateOxidationATPTypes = ["pyruvate'=0", "pyruvate'=cp", "pyruvate'=pyruvate/acetyl-CoA*cp"]
        PyruvateOxidationATPTypes = ["pyruvate'=pyruvate/acetyl-CoA*cp"]

        # TCACycleATPTypes = ["ATP'=0", "ATP'=ct", "ATP'=ADP/ATP*ct"]
        TCACycleATPTypes = ["ATP'=ct", "ATP'=ADP/ATP*ct"]

        ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]
        # ATPSinkTypes = ["No_ATPSink"]
        # ATPSinkTypes = ["Constant_ATPSink"]
        # ATPSinkTypes = ["Burst_ATPSink"]

        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                for GlycolysisATPType in GlycolysisATPTypes:
                    for PyruvateOxidationATPType in PyruvateOxidationATPTypes:
                        for TCACycleATPType in TCACycleATPTypes:

                            for ATPSinkType in ATPSinkTypes:
                                Subtitle = GlucoseAvailabilityType + ", " + GlucoseTransportType + "\n" + GlycolysisATPType + ", " + PyruvateOxidationATPType + "\nTCA_ " + TCACycleATPType + ", " + ATPSinkType
                                self.Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                                self.Simulation.SetGlucoseTransportType(GlucoseTransportType)
                                self.Simulation.SetGlycolysisATPType(GlycolysisATPType)
                                self.Simulation.SetPyruvateOxidationATPType(PyruvateOxidationATPType)
                                self.Simulation.SetTCACycleATPType(TCACycleATPType)
                                self.Simulation.SetATPSinkType(ATPSinkType)
                                self.Simulation.SetGlycolysisEndProduct(EndProduct)
                                self.Simulation.SetAcetylCoALevel(AcetylCoALevel)
                                self.Datasets[Subtitle] = self.Simulation.Sim(Chemotaxis=1, Glycolysis=1, PyruvateOxidation=1, TCACycle=1, ATPSink=1)

        # spacing = "\n"
        spacing = ", "

        PlottingPathway = "(Show all pathways)"
        InclusionList = ["ATP", "ATP_Sink", "glucose{out}", "glucose_Influx", "NADH", "FADH2"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show Chemotaxis Only)"
        InclusionList = ["G6P", "glucose{out}", "glucose_Influx"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show Glycolysis + Pyruvate Oxidation Only)"
        InclusionList = ["ATP", "G6P", "pyruvate", "acetyl-CoA", "ATP_Sink"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show TCA Only)"
        InclusionList = ["ATP", "acetyl-CoA", "NADH", "FADH2", "ATP_Sink"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))


    def Model_Chemotaxis_Glycolysis_PyruvateOxidation_TCACycle_OxidativePhosphorylation(self, EndProduct="pyruvate", AcetylCoALevel="Unfixed"):
        self.InitializeDatasets()
        ModelName = "Chemotaxis + Glycolysis + Pyruvate Oxidation + TCA Cycle + Oxidative Phosphorylation"

        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Burst_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Increasing_GlucoseAvailability"]

        GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        # GlycolysisATPTypes = ["ATP'=cg", "ATP'=ADP/ATP*cg"]
        GlycolysisATPTypes = ["ATP'=ADP/ATP*cg"]

        # PyruvateOxidationATPTypes = ["pyruvate'=0", "pyruvate'=cp", "pyruvate'=pyruvate/acetyl-CoA*cp"]
        PyruvateOxidationATPTypes = ["pyruvate'=pyruvate/acetyl-CoA*cp"]

        # TCACycleATPTypes = ["ATP'=0", "ATP'=ct", "ATP'=ADP/ATP*ct"]
        TCACycleATPTypes = ["ATP'=ADP/ATP*ct"]

        OxidativePhosphorylationATPTypes = ["ATP'=co", "ATP'=ADP/ATP*co"]

        ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]
        # ATPSinkTypes = ["No_ATPSink"]
        # ATPSinkTypes = ["Constant_ATPSink"]
        # ATPSinkTypes = ["Burst_ATPSink"]

        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                for GlycolysisATPType in GlycolysisATPTypes:
                    for PyruvateOxidationATPType in PyruvateOxidationATPTypes:
                        for TCACycleATPType in TCACycleATPTypes:
                            for OxidativePhosphorylationATPType in OxidativePhosphorylationATPTypes:
                                for ATPSinkType in ATPSinkTypes:
                                    Subtitle = GlucoseAvailabilityType + ", " + GlucoseTransportType + "\n" + GlycolysisATPType + ", " + PyruvateOxidationATPType + "\nTCA_ " + TCACycleATPType + ", " + ATPSinkType + ", " + OxidativePhosphorylationATPType
                                    self.Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                                    self.Simulation.SetGlucoseTransportType(GlucoseTransportType)
                                    self.Simulation.SetGlycolysisATPType(GlycolysisATPType)
                                    self.Simulation.SetPyruvateOxidationATPType(PyruvateOxidationATPType)
                                    self.Simulation.SetTCACycleATPType(TCACycleATPType)
                                    self.Simulation.SetATPSinkType(ATPSinkType)
                                    self.Simulation.SetGlycolysisEndProduct(EndProduct)
                                    self.Simulation.SetAcetylCoALevel(AcetylCoALevel)
                                    self.Simulation.SetOxidativePhosphorylationATPType(OxidativePhosphorylationATPType)
                                    self.Datasets[Subtitle] = self.Simulation.Sim(Chemotaxis=1, Glycolysis=1, PyruvateOxidation=1, TCACycle=1, OxidativePhosphorylation=1, ATPSink=1)

        # spacing = "\n"
        spacing = ", "

        PlottingPathway = "(Show all pathways)"
        InclusionList = ["ATP", "ATP_Sink", "glucose{out}", "glucose_Influx", "NADH", "FADH2"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show Chemotaxis Only)"
        InclusionList = ["G6P", "glucose{out}", "glucose_Influx"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show Glycolysis + Pyruvate Oxidation Only)"
        InclusionList = ["ATP", "G6P", "pyruvate", "acetyl-CoA", "ATP_Sink"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

        PlottingPathway = "(Show TCA Only)"
        InclusionList = ["ATP", "acetyl-CoA", "NADH", "FADH2", "ATP_Sink"]
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

#
# class FPlotter:
#     def __init__(self):
#         self.Filter_Inclusion = None
#         self.Filter_Exclusion = None
#
#     def ResetFilters(self):
#         self.Filter_Inclusion = None
#         self.Filter_Exclusion = None
#
#     def SetFilters(self, InclusionList, ExclusionList):
#         self.ResetFilters()
#         self.SetFilter_Inclusion(InclusionList)
#         self.SetFilter_Exclusion(ExclusionList)
#
#     def SetFilter_Inclusion(self, List):
#         self.Filter_Inclusion = List
#
#     def SetFilter_Exclusion(self, List):
#         self.Filter_Exclusion = List
#
#     def CheckToIncludeOrExclude(self, Key_Data):
#         if self.Filter_Inclusion or self.Filter_Exclusion:
#             if self.Filter_Inclusion:
#                 if (Key_Data in self.Filter_Inclusion) or (Key_Data[1:] in self.Filter_Inclusion):
#                     return True
#                 else:
#                     return False
#             else:
#                 if (Key_Data in self.Filter_Exclusion) or (Key_Data[1:] in self.Filter_Exclusion):
#                     return False
#                 else:
#                     return True
#         else:
#             return True
#
#     def FilterDatasets(self, Datasets):
#         Datasets_Filtered = dict()
#         for Key_Dataset, Dataset in Datasets.items():
#             Dataset_Filtered = dict()
#             for Key_Data, Data in Dataset.items():
#                 if self.CheckToIncludeOrExclude(Key_Data):
#                     Dataset_Filtered[Key_Data] = Data
#                     Datasets_Filtered[Key_Dataset] = Dataset_Filtered
#
#         return Datasets_Filtered
#
#     def PlotDatasets(self, Datasets, SimulationTimeUnit, bSideLabel=True, SuperTitle=""):
#
#         # Filter Datasets
#         if self.Filter_Inclusion or self.Filter_Exclusion:
#             Datasets = self.FilterDatasets(Datasets)
#
#         fig = plt.figure()
#         fig.subplots_adjust(wspace=0.5, hspace=0.5)
#         if SuperTitle:
#             fig.suptitle(SuperTitle, fontsize=14)
#
#         Time = None   # Universal X axis (time)
#         Perturbation = 0
#         PerturbationPlotColor = list()
#         for Dataset in Datasets.values():
#             for Data in Dataset.values():
#                 Time = [i * SimulationTimeUnit for i in range(len(Data))]
#                 break
#             for Key in Dataset.keys():
#                 if Key[0] == PerturbationTag:
#                     PerturbationPlotColor.append((random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 0.75)))
#                     Perturbation += 1
#
#         # Plot data
#         NPlotsInRows = len(Datasets)   # Default
#         MaxNPlotsInRows = 3
#         if len(Datasets) > 1:
#             for Remainder in range(MaxNPlotsInRows):
#                 if len(Datasets) % (Remainder + 1) == 0:
#                     NPlotsInRows = Remainder + 1
#
#         for n, (Process, Dataset) in enumerate(Datasets.items()):
#             ax1 = fig.add_subplot(math.ceil(len(Datasets) / NPlotsInRows), NPlotsInRows, n + 1)
#             ax2 = None
#             if Perturbation:
#                 ax2 = ax1.twinx()
#
#             # Y axis (molecular concentrations)
#             PerturbationIndex = 0
#             for MolName, Conc in Dataset.items():
#                 if MolName[0] != PerturbationTag:
#                     line, = ax1.plot(Time, Conc, label="[" + MolName + "]")
#                     if bSideLabel:
#                         SelectedTimeFrameFromLeft = 0.1
#                         ax1.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="left", va="bottom", color=line.get_color())
#                         # ax1.text(Time[-1] * 1.01, Conc[-1], MolName, ha="left", va="bottom", color=line.get_color())
#                         # ax1.text(Time[-1] * 1.1, Conc[-1], MolName + ": {}".format(Conc[-1]), va="center", color=line.get_color())
#
#                 else:
#                     line, = ax2.plot(Time, Conc, color=PerturbationPlotColor[PerturbationIndex], label="[" + MolName[1:] + "]")
#                     if bSideLabel:
#                         SelectedTimeFrameFromLeft = 0.8
#                         ax2.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName[1:], ha="center", va="bottom", color=line.get_color())
#                     # ax2.plot(Time, Conc, label="[" + MolName[2:] + "]")
#                     # print(PerturbationIndex)
#                     PerturbationIndex += 1
#
#             ax1.set_title(Process)
#             ax1.set_xlabel('Time (s)')
#             ax1.set_ylabel('Molecules: Concentration (mol L-1)')
#             ax1.set_ylim([0, 0.015])
#             if not bSideLabel:
#                 ax1.legend(loc='upper left')
#             if Perturbation:
#                 ax2.set_ylabel('Event: Concentration (mol L-1)')
#                 ax2.set_ylim([0, 3e-7])
#                 if not bSideLabel:
#                     ax2.legend(loc='upper right')
#
#             # ax1.grid()
#
#         plt.show()


def main():
    # Initialization
    Simulation = FSimulation()
    Plot = FPlotter()

    # Setting Simulation Parameters
    TotalSimulationTime = 1e-1   # s
    SimulationTimeUnit = 1e-6   # s
    # SteadyStateBrake = False
    # SteadyStateCheckMolecule = "ATP"
    # SteadyStateThresholdFactor = 0.0000001  # Steady state threshold

    Simulation.SetSimulationParameters(TotalSimulationTime, SimulationTimeUnit)
    Simulation.PrintSimulationParameters()

    # Model Runner
    ModelRunner = FModelRunner(Simulation, Plot)

    # Select Models
    # Models = [1]    # Chemotaxis unit test without G6PSink
    # Models = [2]    # Chemotaxis unit test with Linear G6PSink
    # Models = [3]    # Chemotaxis unit test with Burst G6PSink
    # Models = [11]   # Glycolysis unit test without Pyruvate Oxidation
    # Models = [12]   # Glycolysis unit test with Pyruvate Oxidation
    # Models = [21]   # Pyruvate Oxidation unit test
    # Models = [31]   # TCACycle Cycle unit test without acetyl-CoA
    # Models = [32]   # TCACycle Cycle unit test with fixed acetyl-CoA
    # Models = [41]   # Oxidative Phosphorylation unit test
    # Models = [101]  # Chemotaxis + Glycolysis (end product: pyruvate)
    # Models = [102]  # Chemotaxis + Glycolysis (end product: acetyl-CoA)
    # Models = [103]  # Chemotaxis + Glycolysis + Pyruvate Oxidation
    # Models = [104]  # Chemotaxis + Glycolysis + Pyruvate Oxidation + TCA Cycle
    Models = [105]  # Chemotaxis + Glycolysis + Pyruvate Oxidation + TCA Cycle + Oxidative Phosphorylation
    # Models = [1, 2, 3, 11, 12, 21, 22, 31, 41]
    # Models = [101, 102, 103, 104, 105]
    ModelRunner.RunModel(Models)

if __name__ == '__main__':
    main()
