# BSD 3-Clause License
# © 2023, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import random
import sys
import math
import matplotlib.pyplot as plt
from src.lcc.plot2 import FPlotter

random.seed(1)

PerturbationTag = "#"

# Key names for simulation parameters
Key_GlucoseTransport = "GlucoseTransport"
Key_GlucoseAvailability = "GlucoseAvailability"
Key_G6PSink = "G6P_Sink"
Key_GlycolysisATP = "GlycolysisATP"
Key_TCACycleATP = "TCACycleATP"
Key_AcetylCoALevel = "AcetylCoALevel"
Key_PyruvateOxidationATP = "PyruvateOxidationATP"
Key_OxidativePhosphorylationATP = "OxidativePhosphorylationATP"
Key_ATPSink = "ATP_Sink"
Key_GlycolysisEndProduct = "GlycolysisEndProduct"
Key_OxidativePhosphorylationCapacity = "OxidativePhosphorylationCapacity"

# Physiological molecular concentrations
MolConc = dict()
MolConc["G6P"] = 8.8e-3  # all Hexose-P combined
MolConc["ATP"] = 9.6e-3
MolConc["ADP"] = 5.6e-4
MolConc["NAD"] = 2.6e-3
MolConc["NADH"] = 8.3e-3
MolConc["pyruvate"] = 3.9e-4
MolConc["acetyl-CoA"] = 6.1e-4
MolConc["FAD"] = 1.7e-4
MolConc["FADH2"] = 5e-4  # undocumented
MolConc["CoA"] = 1.4e-3
MolConc["Malate"] = 1.7e-3   # for TCA capacity

# Perturbation initialization
Perturbation = dict()
Perturbation["ATP_Sink"] = 0
Perturbation["glucose{out}"] = 0
Perturbation["glucose_Influx"] = 0
Perturbation["G6P_Sink"] = 0

# Constant Setup
TargetK = "cg"
Base = 2 # 2 is default
c_CellDivision = 1e-3
Capacity = dict()
Capacity["cg"] = c_CellDivision * 19
Capacity["cp"] = c_CellDivision * 3
Capacity["ct"] = c_CellDivision * 385
Capacity["co"] = c_CellDivision * 5000

# Debugging Switch
# DebugZero = True
DebugZero = False

# Stoichiometry Analysis:
# Glyc:     G6P + 2 ADP + 2 NAD -> 2 pyruvate + 2 NADH + 2 ATP
# PyrOx:    pyruvate + CoA + NAD -> acetyl-CoA + NADH
# TCA:      acetyl-CoA + 3 NAD + FAD + ADP -> 3 NADH + FADH2 + ATP + CoA
# OxPhos:   NADH + FADH2 + 4 ADP -> NAD + FAD + 4 ATP
# Overall:  G6P + 12 ADP + 8 NAD -> 12 ATP + 8 NADH

class FSimulation():
    def __init__(self):
        # Molecule Concentration
        self.MolConc = dict()   # current concentration
        self.DeltaMolConc = dict()
        self.Dataset = dict()   # storing concentration
        self.Perturbation = dict()
        self.Conjugation = dict()

        # Capacity
        self.Capacity = dict()

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
        self.InitializeSelectiveMolecules = list()
        self.InitializePerturbations = list()
        self.RestoreSelectiveMolecules = list()

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

    def InitializeConjugation(self):
        if 'ATP' in self.MolConc:
            self.Conjugation['ATP'] = ('ADP', self.MolConc['ATP'] + self.MolConc['ADP'])
        if 'NADH' in self.MolConc:
            self.Conjugation['NADH'] = ('NAD', self.MolConc['NADH'] + self.MolConc['NAD'])
        if 'FADH2' in self.MolConc:
            self.Conjugation['FADH2'] = ('FAD', self.MolConc['FADH2'] + self.MolConc['FAD'])

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
        self.SupportedModelTypes2Int["ATP'=(ADP/ATP)*cg"] = 2
        self.SupportedModelTypes2Int["ATP=G6P*cg"] = 3

        # Glycolysis: End product
        self.SupportedModelTypes2Int["pyruvate"] = 0
        self.SupportedModelTypes2Int["acetyl-CoA"] = 1

        # Pyruvate Oxidation - ATP Production
        self.SupportedModelTypes2Int["acetyl-CoA'=0"] = 0
        self.SupportedModelTypes2Int["acetyl-CoA'=cp"] = 1
        self.SupportedModelTypes2Int["acetyl-CoA'=(ADP/ATP)*cp"] = 2   # Not testable without ATP and ADP
        self.SupportedModelTypes2Int["acetyl-CoA'=(pyruvate/acetyl-CoA)*cp"] = 3
        self.SupportedModelTypes2Int["acetyl-CoA'=pyruvate*cp"] = 4

        # TCA Cycle: ATP Production
        # self.SupportedModelTypes2Int["ATP'=0"] = 0   # same as in Glycolysis ATP production
        self.SupportedModelTypes2Int["ATP'=ct"] = 1
        self.SupportedModelTypes2Int["ATP'=(ADP/ATP)*ct"] = 2
        self.SupportedModelTypes2Int["ATP'=(acetyl-CoA/CoA)*ct"] = 3
        self.SupportedModelTypes2Int["ATP'=acetyl-CoA*ct"] = 4

        # TCA Cycle: acetyl-CoA level control
        self.SupportedModelTypes2Int["Unfixed"] = 0
        self.SupportedModelTypes2Int["Fixed"] = 1

        # Oxidative Phosphorylation: ATP Production
        # self.SupportedModelTypes2Int["ATP'=0"] = 0   # same as in Glycolysis ATP production
        self.SupportedModelTypes2Int["ATP'=co"] = 1
        self.SupportedModelTypes2Int["ATP'=(ADP/ATP)*co"] = 2
        self.SupportedModelTypes2Int["ATP'=(1/ATP)*co"] = 3

        # Oxidative Phosphorylation: Capacity (co) determination
        self.SupportedModelTypes2Int["co=min(FADH2,NADH)"] = 0
        self.SupportedModelTypes2Int["co=min(FADH2,NADH)*0.1"] = 1
        self.SupportedModelTypes2Int["co=min(FADH2,NADH)*0.01"] = 2

        # Arbitrary ATP Sink
        self.SupportedModelTypes2Int["No_ATPSink"] = 0
        self.SupportedModelTypes2Int["Constant_ATPSink"] = 1
        self.SupportedModelTypes2Int["Burst_ATPSink"] = 2
        self.SupportedModelTypes2Int["CellDivision_ATPSink"] = 3

    def GetModelType(self, Type):
        if Type in self.ModelType.keys():
            return self.ModelType[Type]
        else:
            print("\n[Error] Unrecognized Model: %s" % Type)
            sys.exit(1)

    def GetSupportedModelType2Int(self, Type):
        if Type in self.SupportedModelTypes2Int.keys():
            return self.SupportedModelTypes2Int[Type]
        else:
            print("\n[Error] Unsupported Model Type: %s" % Type)
            sys.exit(1)

    def AdjustCapacityForCellDivisionSetting(self):
        self.Capacity["cg"] = Capacity["cg"]
        self.Capacity["cp"] = Capacity["cp"]
        # self.Capacity["ct"] = Capacity["ct"]
        # self.Capacity["co"] = Capacity["co"]

    def SetInitializeSelectiveMolecules(self, Molecules):
        for Molecule in Molecules:
            self.InitializeSelectiveMolecules.append(Molecule)

    def SetInitializePerturbations(self, Perturbations):
        for Perturbation in Perturbations:
            self.InitializePerturbations.append(Perturbation)

    def SetRestoreSelectiveMolecules(self, Molecules):
        for Molecule in Molecules:
            self.RestoreSelectiveMolecules.append(Molecule)

    def AppendData(self):
        for MolName, Conc in self.MolConc.items():
            self.Dataset[MolName].append(Conc)
        for PerturbationName, Conc in self.Perturbation.items():
            self.Dataset[PerturbationName].append(Conc)

    def UpdateMolConc(self):
        # If concentrations getting below zero
        if self.CheckBelowZero(self.DeltaMolConc):
            return True   # For KFinder

        # If concentrations all above zero
        else:
            for MolName, Conc in self.DeltaMolConc.items():
                self.MolConc[MolName] += self.DeltaMolConc[MolName]
                self.DeltaMolConc[MolName] = 0
            return False   # For KFinder

    def AddToDeltaMolConc(self, NewMolConcentrations):
        for MolName, Conc in NewMolConcentrations.items():
            self.DeltaMolConc[MolName] += NewMolConcentrations[MolName]

    def UpdateCapacity(self):
        if Key_OxidativePhosphorylationATP in self.ModelType:
            self.Capacity["co"] = min(self.MolConc["NADH"], self.MolConc["FADH2"]) * Capacity["co"]
        if Key_TCACycleATP in self.ModelType:
            self.Capacity["ct"] = self.MolConc["Malate"] * Capacity["ct"]

    def UpdateCapacityWithKFinderFactor(self, KeyFinderFactor):
        if KeyFinderFactor[0] == "co" or KeyFinderFactor[0] == "ct":
            self.Capacity[KeyFinderFactor[0]] *= KeyFinderFactor[1]
        else:
            self.Capacity[KeyFinderFactor[0]] = Capacity[KeyFinderFactor[0]] * KeyFinderFactor[1]
        # print(self.Capacity)

    def CheckBelowZero(self, NewMolConcentrations):
        bBelowZero = False
        for MolName, Conc in NewMolConcentrations.items():
            Test = self.MolConc[MolName] + NewMolConcentrations[MolName]
            if Test < 0:
                # print("Below zero: %s %s\t" % (MolName, str(Test)), end= " |\t")
                # self.PrintSimTimeAndMolConc()
                bBelowZero = True

                if DebugZero:
                    print("\nBelow Zero: ", MolName, "\tConcentration: ", self.MolConc[MolName], "\tdelta: ", NewMolConcentrations[MolName], end= "")
        return bBelowZero

    def CheckZeroSumConjugation(self):
        for Molecule, (Molecule_Conjugate, Conc_InitialSum) in self.Conjugation.items():
            Conc_Sum = self.MolConc[Molecule] + self.MolConc[Molecule_Conjugate]
            Ratio = abs(Conc_Sum / Conc_InitialSum)
            UpperBoundCheck = Ratio < 1 + 1e-6
            LowerBoundCheck = Ratio > 1 - 1e-6
            assert UpperBoundCheck and LowerBoundCheck , "Conjugation Pair: %s - %s" % (Molecule, Molecule_Conjugate) + " |\tInitial Sum: %s vs. \tCurrent Sum: %s" % (Conc_InitialSum, Conc_Sum)

    def PrintMolConc(self, Order=False):
        if Order:
            MolList = ["G6P", "ADP", "NAD", "pyruvate", "NADH", "ATP", "CoA", "acetyl-CoA", "FAD", "Malate", "FADH2"]   # matches lcc
            for MolName in MolList:
                print("[" + MolName + "] " + "{:.6f}".format(self.MolConc[MolName]), end=", \t")
        else:
            for MolName, Conc in self.MolConc.items():
                print("[" + MolName + "] " + "{:.6f}".format(Conc), end=", \t")

    def PrintDeltaMolConc(self):
        for MolName, Conc in self.DeltaMolConc.items():
            print("d[" + MolName + "] " + "{:.6f}".format(Conc), end=", \t")

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
        self.ModelType[Key_AcetylCoALevel] = self.GetSupportedModelType2Int(Type)

    def SetTCACycleATPType(self, Type):
        self.ModelType[Key_TCACycleATP] = self.GetSupportedModelType2Int(Type)

    def SetOxidativePhosphorylationATPType(self, Type):
        self.ModelType[Key_OxidativePhosphorylationATP] = self.GetSupportedModelType2Int(Type)

    def SetOxidativePhosphorylationCapacityType(self, Type):
        self.ModelType[Key_OxidativePhosphorylationCapacity] = self.GetSupportedModelType2Int(Type)

    def SetATPSinkType(self, Type):
        self.ModelType[Key_ATPSink] = self.GetSupportedModelType2Int(Type)

    def SetGlycolysisEndProduct(self, EndProduct):
        self.ModelType[Key_GlycolysisEndProduct] = self.GetSupportedModelType2Int(EndProduct)

    def Sim(self, KFinderFactor=("", 0), Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCACycle=0, OxidativePhosphorylation=0, G6PSink=0, ATPSink=0):
        self.PrintSimulationSubject()
        self.Initialize(Chemotaxis, Glycolysis, PyruvateOxidation, TCACycle, OxidativePhosphorylation)
        self.Run(KFinderFactor, Chemotaxis, Glycolysis, PyruvateOxidation, TCACycle, OxidativePhosphorylation, G6PSink, ATPSink)
        self.CorrectAndTagPerturbationLabel()
        return self.Dataset.copy()

    def Initialize_SelectiveMolecules(self):
        for Molecule in self.InitializeSelectiveMolecules:
            self.MolConc[Molecule] = MolConc[Molecule]

    def Initialize_Perturbations(self):
        for Molecule in self.InitializePerturbations:
            self.Perturbation[Molecule] = Perturbation[Molecule]

    def Restore_SelectiveMolecules(self):
        for Molecule in self.RestoreSelectiveMolecules:
            self.MolConc[Molecule] = MolConc[Molecule]

    def Initialize_Chemotaxis(self):
        # Initial Molecule Concentrations
        self.MolConc["G6P"] = 8.8e-3   # all Hexose-P combined
        self.MolConc["ATP"] = 9.6e-3
        self.MolConc["ADP"] = 5.6e-4

        # Perturbations
        self.Perturbation["glucose{out}"] = 0
        self.Perturbation["glucose_Influx"] = 0
        if Key_G6PSink in self.ModelType:
            self.Perturbation["G6P_Sink"] = 0

        # Capacity
        self.Capacity["cc"] = 0.01
        self.Capacity["cgt"] = 0.5

    def GetGlucoseAvailability_Chemotaxis(self):
        Model = self.GetModelType(Key_GlucoseAvailability)
        if Model == 0:
            return 0
        elif Model == 1:
            return 0.15
        elif Model == 2:
            # return 0.000003 * self.SimStep # Previous value
            return 0.000006 * self.SimStep
        elif Model == 3:
            return abs(0.15 * math.sin(self.SimStep/self.TotalSimulationTime))
        else:
            sys.exit(1)

    def GetTransportFactor_Chemotaxis(self):
        Model = self.GetModelType(Key_GlucoseTransport)
        if Model == 0:
            return 0
        elif Model == 1:
            return self.Capacity["cgt"]
        elif Model == 2:
            return 0.0001 / (self.Dataset["G6P"][-1] ** 2) * self.Capacity["cgt"]
        else:
            sys.exit(1)

    def Sim_Chemotaxis(self):
        d = dict()

        # Production Rate Type
        GlucoseAvailability = self.GetGlucoseAvailability_Chemotaxis()
        self.Perturbation["glucose{out}"] = GlucoseAvailability * self.SimulationTimeUnit

        TransportFactor = self.GetTransportFactor_Chemotaxis()
        TransportDemand = TransportFactor if TransportFactor <= 1 else 1   # TODO: make the limit approaching smooth
        self.Perturbation["glucose_Influx"] = self.Perturbation["glucose{out}"] * TransportDemand

        d["G6P"] = self.Perturbation["glucose_Influx"]
        d["ATP"] = - self.Capacity["cc"] * self.SimulationTimeUnit
        d["ADP"] = - d["ATP"]   # use when G6P is constant

        self.AddToDeltaMolConc(d)

    def Initialize_Glycolysis(self):
        # Without pyruvate oxidation
        # G6P + 2 ADP + 2 NAD -> 2 pyruvate + 2 NADH + 2 ATP
        # With pyruvate oxidation
        # G6P + 2 ADP + 4 NAD -> 2 acetyl-CoA + 4 NADH + 2 ATP

        # Initial Molecule Concentrations
        self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["ADP"]      = 5.6e-4
        self.MolConc["NAD"]      = 2.6e-3
        self.MolConc["NADH"]     = 8.3e-3
        self.MolConc["ATP"]      = 9.6e-3
        if self.ModelType[Key_GlycolysisEndProduct] == 0:
            self.MolConc["pyruvate"] = 3.9e-4
        elif self.ModelType[Key_GlycolysisEndProduct] == 1:
            self.MolConc["acetyl-CoA"] = 6.1e-4

        # Perturbations
        if Key_ATPSink in self.ModelType:
            self.Perturbation["ATP_Sink"] = 0
        self.Perturbation["ATP_Glycolysis"] = 0

        # Capacity
        self.Capacity["cg"] = 10

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

    def Sim_Glycolysis(self):
        # G6P + 2 ADP + 2 NAD -> 2 pyruvate + 2 NADH + 2 ATP
        # G6P + 2 ADP + 4 NAD -> 2 acetyl-CoA + 4 NADH + 2 ATP

        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_Glycolysis()
        dATP = (InputFactor * self.Capacity["cg"]) * self.SimulationTimeUnit

        d["ATP"] = dATP
        d["ADP"] = - d["ATP"]
        d["G6P"] = - d["ATP"] / 2
        if self.ModelType[Key_GlycolysisEndProduct] == 0:
            d["pyruvate"] = d["ATP"]
            d["NADH"] = d["ATP"]
            d["NAD"] = - d["ATP"]

        elif self.ModelType[Key_GlycolysisEndProduct] == 1:
            d["acetyl-CoA"] = d["ATP"]
            d["NADH"] = d["ATP"] * 2
            d["NAD"] = - d["ATP"] * 2

        self.Perturbation["ATP_Glycolysis"] = dATP

        self.AddToDeltaMolConc(d)

    def GetG6PSinkRate(self):
        Model = self.GetModelType(Key_G6PSink)
        if Model == 0:
            return 0
        elif Model == 1:   # linear
            return (self.TotalSimulationTime * 1.4)
        elif Model == 2:   # non-linear
            return (self.TotalSimulationTime * 0.1) * (self.TotalSimulationTime / 40 - (self.TotalSimulationTime / 2 - self.SimulationTimeUnit * self.SimStep) ** 2) * 1e4
        else:
            sys.exit(1)


    def Sim_G6PSink(self):
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
        elif Model == 3:   # cell division
            return 4.35E-03
        else:
            sys.exit(1)

    def Sim_ATPSink(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCACycle=0, OxidativePhosphorylation=0):
        d = dict()

        dATP = self.GetATPSinkRate() * self.SimulationTimeUnit

        # Combinatorial case
        if Key_ATPSink in self.ModelType and self.ModelType[Key_ATPSink] == 3:   # Cell Division
            d["NAD"] = dATP * 0.7725
            # d["NAD"] = dATP
            d["NADH"] = - d["NAD"]
            pass
        elif Chemotaxis and Glycolysis and PyruvateOxidation and TCACycle:
            d["NAD"] = dATP * 0.5
            d["NADH"] = - d["NAD"]
            d["FAD"] = dATP * 0.025
            d["FADH2"] = - d["FAD"]
        elif Chemotaxis and Glycolysis and PyruvateOxidation:
            d["NAD"] = dATP
            d["NADH"] = - d["NAD"]
        elif Chemotaxis and Glycolysis:
            if self.ModelType[Key_GlycolysisEndProduct] == 0:
                d["NAD"] = dATP
                d["NADH"] = - d["NAD"]
            elif self.ModelType[Key_GlycolysisEndProduct] == 1:
                d["NAD"] = dATP * 2
                d["NADH"] = - d["NAD"]
        # Single case
        elif Chemotaxis:
            pass
        elif Glycolysis:
            if self.ModelType[Key_GlycolysisEndProduct] == 0:
                d["NAD"] = dATP
                d["NADH"] = - d["NAD"]
            elif self.ModelType[Key_GlycolysisEndProduct] == 1:
                d["NAD"] = dATP * 2
                d["NADH"] = - d["NAD"]
        elif TCACycle:
            d["NAD"] = dATP * 0.5
            d["NADH"] = - d["NAD"]
            d["FAD"] = dATP * 0.1
            d["FADH2"] = - d["FAD"]
        elif OxidativePhosphorylation:
            pass

        self.Perturbation["ATP_Sink"] = dATP
        d["ATP"] = - dATP
        d["ADP"] = dATP

        self.AddToDeltaMolConc(d)

    def Initialize_PyruvateOxidation(self):
        # With pyruvate oxidation
        # pyruvate + CoA + NAD -> acetyl-CoA + NADH

        # Initial Molecule Concentrations
        self.MolConc["pyruvate"]    = 3.9e-4
        self.MolConc["acetyl-CoA"]  = 6.1e-4
        self.MolConc["NAD"]         = 2.6e-3
        self.MolConc["NADH"]        = 8.3e-3
        self.MolConc["CoA"]         = 1.4e-3

        # Capacity
        self.Capacity["cp"] = 0.1

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
        elif Model == 4:
            return self.Dataset["pyruvate"][-1]
        else:
            sys.exit(1)

    def Sim_PyruvateOxidation(self):
        # pyruvate + CoA + NAD -> acetyl-CoA + NADH
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_PyruvateOxidation()
        dacetylCoA = (InputFactor * self.Capacity["cp"]) * self.SimulationTimeUnit
        # print("\n", InputFactor, self.Capacity["cp"], dacetylCoA)

        d["acetyl-CoA"] = dacetylCoA
        d["pyruvate"] = - dacetylCoA
        d["CoA"] = - dacetylCoA
        d["NADH"] = dacetylCoA
        d["NAD"] = - dacetylCoA

        self.AddToDeltaMolConc(d)

    def Initialize_TCACycle(self):
        # acetyl-CoA + 3 NAD + FAD + ADP -> 3 NADH + FADH2 + ATP + CoA

        # Initial Molecule Concentrations
        # self.MolConc["G6P"]      = 8.8e-3   # all Hexose-P combined
        self.MolConc["acetyl-CoA"]  = 6.1e-4
        self.MolConc["FAD"]     = 1.7e-4
        self.MolConc["NAD"]     = 2.6e-3
        self.MolConc["ADP"]     = 5.6e-4
        self.MolConc["NADH"]    = 8.3e-3
        self.MolConc["FADH2"]   = 5e-4   # undocumented
        self.MolConc["ATP"]     = 9.6e-3
        self.MolConc["CoA"]     = 1.4e-3
        self.MolConc["Malate"]  = MolConc["Malate"]

        # Perturbations
        if Key_ATPSink in self.ModelType:
            self.Perturbation["ATP_Sink"] = 0
        self.Perturbation["ATP_TCA"] = 0

        # Capacity
        self.Capacity["ct"]     = 0.3


    def GetInputFactor_TCACycle(self):
        Model = self.GetModelType(Key_TCACycleATP)
        if Model == 0:
            return 0
        elif Model == 1:
            return 1
        elif Model == 2:
            return self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        elif Model == 3:
            return self.Dataset["acetyl-CoA"][-1] / self.Dataset["CoA"][-1]
        elif Model == 4:
            return self.Dataset["acetyl-CoA"][-1]
        else:
            sys.exit(1)

    def Sim_TCACycle(self):
        # acetyl-CoA + 3 NAD + FAD + ADP -> 3 NADH + FADH2 + ATP + CoA
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_TCACycle()
        dATP = (InputFactor * self.Capacity["ct"]) * self.SimulationTimeUnit
        # print("\n", InputFactor, self.Capacity["ct"], dATP)

        d["ATP"] = dATP
        d["ADP"] = - d["ATP"]

        d["NADH"] = d["ATP"] * 3
        d["NAD"] = - d["NADH"]
        d["FADH2"] = d["ATP"]
        d["FAD"] = - d["FADH2"]

        if self.GetModelType(Key_AcetylCoALevel) == 0:
            d["acetyl-CoA"] = - d["ATP"]
            d["CoA"] = d["ATP"]

        self.Perturbation["ATP_TCA"] = dATP

        self.AddToDeltaMolConc(d)

    def Initialize_OxidativePhosphorylation(self):
        # NADH + FADH2 + 4 ADP -> NAD + FAD + 4 ATP

        # Initial Molecule Concentrations
        self.MolConc["FAD"]     = 1.7e-4
        self.MolConc["NAD"]     = 2.6e-3
        self.MolConc["ADP"]     = 5.6e-4
        self.MolConc["NADH"]    = 8.3e-3
        self.MolConc["FADH2"]   = 5e-4   # undocumented
        self.MolConc["ATP"]     = 9.6e-3

        # Perturbations
        if Key_ATPSink in self.ModelType:
            self.Perturbation["ATP_Sink"] = 0
        self.Perturbation["ATP_OxPhos"] = 0

        # Capacity
        self.Capacity["co"]     = 1

    def GetInputFactor_OxidativePhosphorylation(self):
        Model = self.GetModelType(Key_OxidativePhosphorylationATP)
        if Model == 0:
            return 0
        elif Model == 1:
            return 1
        elif Model == 2:
            return self.Dataset["ADP"][-1] / self.Dataset["ATP"][-1]
        elif Model == 3:
            return 1 / self.Dataset["ATP"][-1]
        elif Model == 4:
            return 1
        else:
            sys.exit(1)

    def Sim_OxidativePhosphorylation(self):
        # NADH + FADH2 + 4 ADP -> NAD + FAD + 4 ATP
        d = dict()

        # Production Rate Type
        InputFactor = self.GetInputFactor_OxidativePhosphorylation()
        dATP = (InputFactor * self.Capacity["co"]) * self.SimulationTimeUnit
        # print("\n", InputFactor, self.Capacity["co"], dATP)

        d["ATP"] = dATP
        d["ADP"] = - dATP

        d["NADH"] = - dATP / 4
        d["NAD"] = dATP / 4
        d["FADH2"] = - dATP / 4
        d["FAD"] = dATP / 4

        self.Perturbation["ATP_OxPhos"] = dATP

        self.AddToDeltaMolConc(d)


    def Initialize(self, Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCACycle=0, OxidativePhosphorylation=0):
        if self.InitializeSelectiveMolecules:
            self.Initialize_SelectiveMolecules()
        if self.InitializePerturbations:
            self.Initialize_Perturbations()
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
        if Key_ATPSink in self.ModelType:
            if self.ModelType[Key_ATPSink] == 3:
                self.AdjustCapacityForCellDivisionSetting()

        self.InitializeDeltaMolConc()
        self.InitializeDataset()
        self.InitializeConjugation()
        self.InitializeSimStep()


    def Run(self, KFinderFactor=("", 0), Chemotaxis=0, Glycolysis=0, PyruvateOxidation=0, TCACycle=0,
            OxidativePhosphorylation=0, G6PSink=0, ATPSink=0):

        while self.SimStep * self.SimulationTimeUnit < self.TotalSimulationTime:
            # Debugging
            # if self.SimStep < 10:
            #     self.PrintSimTimeAndMolConc()

            if Key_OxidativePhosphorylationATP in self.ModelType or Key_TCACycleATP in self.ModelType:
                self.UpdateCapacity()
            if KFinderFactor[1]:
                self.UpdateCapacityWithKFinderFactor(KFinderFactor)

            if Chemotaxis:
                self.Sim_Chemotaxis()
            if Glycolysis:
                self.Sim_Glycolysis()
            if PyruvateOxidation:
                self.Sim_PyruvateOxidation()
            if TCACycle:
                self.Sim_TCACycle()
            if OxidativePhosphorylation:
                self.Sim_OxidativePhosphorylation()
            if G6PSink:
                self.Sim_G6PSink()
            if ATPSink:
                self.Sim_ATPSink(Chemotaxis, Glycolysis, PyruvateOxidation, TCACycle,
                                        OxidativePhosphorylation)

            # Debugging
            if self.SimStep < 10:
                self.PrintSimTimeAndMolConc()

            if self.RestoreSelectiveMolecules:
                self.Restore_SelectiveMolecules()

            self.CheckZeroSumConjugation()
            self.AppendData()

            bBelowZero = self.UpdateMolConc()
            if KFinderFactor[1] and bBelowZero:
                print("\n***** Molecular Concentration crashed with KFinderFactor: ", KFinderFactor)
                break

            if self.SteadyStateBrake:
                if (abs(self.Dataset[self.SteadyStateCheckMolecule][-1] -
                        self.Dataset[self.SteadyStateCheckMolecule][-2]) /
                        self.Dataset[self.SteadyStateCheckMolecule][-1] < self.SteadyStateThresholdFactor):
                    break
            self.SimStep += 1

        # Final Log
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
                self.Model_OxidativePhosphorylation(OP_CapacityTest=True)
            elif Model == 42:
                self.Model_OxidativePhosphorylation(OP_CapacityTest=False)
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
            elif Model == 201:
                self.Model_GlucoseMetabolism_CellDivision(ConstantByproducts=True)
            elif Model == 202:
                self.Model_GlucoseMetabolism_CellDivision(ConstantByproducts=False)
            elif Model == 203:
                self.Model_GlucoseMetabolism_CellDivision(ConstantByproducts=False, ConstantG6P=False)
            elif Model == 211:
                self.Model_GlucoseMetabolism_CellDivision_KFinder()
            else:
                print("\n[Error] Unsupported Model: %s" % Model)
                sys.exit(1)

    def RemoveEmptyDatasets(self):
        NewDatasets = dict()
        for Key, Dataset in self.Datasets.items():
            for Data in Dataset.values():
                if len(Data) > 2:
                    NewDatasets[Key] = Dataset
                    break
        self.Datasets = NewDatasets

    def Model_Chemotaxis(self, G6PSinkModel="Linear_G6PSink"):
        self.InitializeDatasets()
        ModelName = "[Chemotaxis] glucose{out} -> G6P{in}"

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
            ModelName = "[Glycolysis] G6P + 2 ADP + 2 NAD+ + 2 Pi -> 2 pyruvate + 2 NADH + 2 ATP"
        elif EndProduct =="acetyl-CoA":
            ModelName = "[Glycolysis] G6P + 2 ADP + 3 NAD+ + 2 Pi -> 2 acetyl-CoA + 4 NADH + 2 ATP"
        else:
            print("\n[Error] Unsupported End Product: %s" % EndProduct)
            sys.exit(1)

        # GlycolysisATPTypes = ["ATP'=0", "ATP'=cg", "ATP'=(ADP/ATP)*cg"]
        GlycolysisATPTypes = ["ATP'=cg", "ATP'=(ADP/ATP)*cg"]
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
        ModelName = "[PyruvateOxidation] pyruvate + CoA-SH -> acetyl-CoA + NADH"

        # PyruvateOxidationATPTypes = ["acetyl-CoA'=0", "acetyl-CoA'=cp", "acetyl-CoA'=(ADP/ATP)*cp", "acetyl-CoA'=(pyruvate/acetyl-CoA)*cp"]
        PyruvateOxidationATPTypes = ["acetyl-CoA'=0", "acetyl-CoA'=cp", "acetyl-CoA'=(pyruvate/acetyl-CoA)*cp", "acetyl-CoA'=pyruvate*cp"]
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
        ModelName = "[TCACycle] acetyl-CoA + 3 NAD + FAD + ADP -> 3 NADH + FADH2 + ATP + CoA"

        # TCACycleATPTypes = ["ATP'=0", "ATP'=ct", "ATP'=(ADP/ATP)*ct"]
        TCACycleATPTypes = ["ATP'=ct", "ATP'=(ADP/ATP)*ct"]

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

    def Model_OxidativePhosphorylation(self, OP_CapacityTest=True):
        self.InitializeDatasets()
        ModelName = "[OxidativePhosphorylation] NADH + FADH2 + 4 ADP -> NAD + FAD + 4 ATP"   # Not including H other biproducts

        if OP_CapacityTest:
            OxidativePhosphorylationATPType = "ATP'=co"
            OxidativePhosphorylationCapacityTypes = ["co=min(FADH2,NADH)", "co=min(FADH2,NADH)*0.1", "co=min(FADH2,NADH)*0.01"]
            ATPSinkTypes = ["No_ATPSink", "Constant_ATPSink", "Burst_ATPSink"]

            for OxidativePhosphorylationCapacityType in OxidativePhosphorylationCapacityTypes:
                self.Simulation.SetOxidativePhosphorylationCapacityType(OxidativePhosphorylationCapacityType)
                for ATPSinkType in ATPSinkTypes:
                    Subtitle = OxidativePhosphorylationATPType + "\n" + OxidativePhosphorylationCapacityType + "\n" + ATPSinkType
                    self.Simulation.SetTitle(ModelName + "\n" + Subtitle)
                    self.Simulation.SetOxidativePhosphorylationATPType(OxidativePhosphorylationATPType)
                    self.Simulation.SetATPSinkType(ATPSinkType)
                    self.Datasets[Subtitle] = self.Simulation.Sim(OxidativePhosphorylation=1, ATPSink=1)

        else:
            OxidativePhosphorylationATPTypes = ["ATP'=(ADP/ATP)*co", "ATP'=(1/ATP)*co"]
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

        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Fluctuating_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability"]

        GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        # GlycolysisATPTypes = ["ATP'=cg", "ATP'=(ADP/ATP)*cg"]
        GlycolysisATPTypes = ["ATP'=(ADP/ATP)*cg"]

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

        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Fluctuating_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability"]

        GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        # GlycolysisATPTypes = ["ATP'=cg", "ATP'=(ADP/ATP)*cg"]
        GlycolysisATPTypes = ["ATP'=(ADP/ATP)*cg"]

        # PyruvateOxidationATPTypes = ["acetyl-CoA'=0", "acetyl-CoA'=cp", "acetyl-CoA'=(pyruvate/acetyl-CoA)*cp"]
        PyruvateOxidationATPTypes = ["acetyl-CoA'=pyruvate*cp"]

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

        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Fluctuating_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Increasing_GlucoseAvailability"]

        GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        # GlycolysisATPTypes = ["ATP'=cg", "ATP'=(ADP/ATP)*cg"]
        GlycolysisATPTypes = ["ATP'=(ADP/ATP)*cg"]

        # PyruvateOxidationATPTypes = ["acetyl-CoA'=0", "acetyl-CoA'=cp", "acetyl-CoA'=(pyruvate/acetyl-CoA)*cp"]
        PyruvateOxidationATPTypes = ["acetyl-CoA'=pyruvate*cp"]

        # TCACycleATPTypes = ["ATP'=0", "ATP'=ct", "ATP'=(ADP/ATP)*ct"]
        TCACycleATPTypes = ["ATP'=ct", "ATP'=(ADP/ATP)*ct"]

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

        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Fluctuating_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Increasing_GlucoseAvailability"]

        GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        # GlycolysisATPTypes = ["ATP'=cg", "ATP'=(ADP/ATP)*cg"]
        GlycolysisATPTypes = ["ATP'=(ADP/ATP)*cg"]

        # PyruvateOxidationATPTypes = ["acetyl-CoA'=0", "acetyl-CoA'=cp", "acetyl-CoA'=(pyruvate/acetyl-CoA)*cp"]
        PyruvateOxidationATPTypes = ["acetyl-CoA'=pyruvate*cp"]

        # TCACycleATPTypes = ["ATP'=0", "ATP'=ct", "ATP'=(ADP/ATP)*ct"]
        TCACycleATPTypes = ["ATP'=(ADP/ATP)*ct"]

        OxidativePhosphorylationATPTypes = ["ATP'=co"]

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
                                    self.Simulation.SetOxidativePhosphorylationCapacityType("co=min(FADH2,NADH)")
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

    def Model_GlucoseMetabolism_CellDivision(self, ConstantByproducts=True, ConstantG6P=True):
        self.InitializeDatasets()
        ModelName = "Glucose Metabolism + CellDivision"

        GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Fluctuating_GlucoseAvailability"]
        GlucoseAvailabilityTypes = ["Increasing_GlucoseAvailability"]
        GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        GlycolysisATPTypes = ["ATP'=(ADP/ATP)*cg"]
        PyruvateOxidationATPTypes = ["acetyl-CoA'=(pyruvate/acetyl-CoA)*cp"]
        TCACycleATPTypes = ["ATP'=(ADP/ATP)*ct", "ATP'=(acetyl-CoA/CoA)*ct"]
        OxidativePhosphorylationATPTypes = ["ATP'=co"]
        OxidativePhosphorylationCapacityTypes = ["co=min(FADH2,NADH)", "co=min(FADH2,NADH)*0.1", "co=min(FADH2,NADH)*0.01"]
        ATPSinkTypes = ["CellDivision_ATPSink"]

        for GlucoseAvailabilityType in GlucoseAvailabilityTypes:
            for GlucoseTransportType in GlucoseTransportTypes:
                for GlycolysisATPType in GlycolysisATPTypes:
                    for PyruvateOxidationATPType in PyruvateOxidationATPTypes:
                        for TCACycleATPType in TCACycleATPTypes:
                            for OxidativePhosphorylationCapacityType in OxidativePhosphorylationCapacityTypes:
                                for ATPSinkType in ATPSinkTypes:
                                    # Subtitle = GlucoseAvailabilityType + ", " + GlucoseTransportType + "\n" + GlycolysisATPType + ", " + PyruvateOxidationATPType + "\nTCA_ " + TCACycleATPType + ", " + ATPSinkType + ", " + OxidativePhosphorylationATPType
                                    Subtitle = GlucoseAvailabilityType + ", " + GlucoseTransportType + "\n[Glyc]" + GlycolysisATPType + ", [PO]" + PyruvateOxidationATPType + "\n[TCA]" + TCACycleATPType + ", [OP]" + OxidativePhosphorylationCapacityType + "\n" + ATPSinkType

                                    self.Simulation.SetATPSinkType(ATPSinkType)
                                    self.Simulation.SetGlycolysisATPType(GlycolysisATPType)
                                    self.Simulation.SetPyruvateOxidationATPType(PyruvateOxidationATPType)
                                    self.Simulation.SetTCACycleATPType(TCACycleATPType)
                                    self.Simulation.SetATPSinkType(ATPSinkType)
                                    self.Simulation.SetGlycolysisEndProduct("pyruvate")
                                    self.Simulation.SetAcetylCoALevel("Unfixed")
                                    self.Simulation.SetOxidativePhosphorylationATPType("ATP'=co")
                                    self.Simulation.SetOxidativePhosphorylationCapacityType(OxidativePhosphorylationCapacityType)

                                    self.Simulation.SetInitializeSelectiveMolecules(["acetyl-CoA", "pyruvate","NAD", "NADH", "FAD", "FADH2", "G6P", "CoA"])
                                    if ConstantByproducts:
                                        self.Simulation.SetRestoreSelectiveMolecules(["acetyl-CoA", "pyruvate","NAD", "NADH", "FAD", "FADH2", "G6P", "CoA"])
                                        self.Datasets[Subtitle] = self.Simulation.Sim(Glycolysis=1, PyruvateOxidation=1, TCACycle=1, OxidativePhosphorylation=1, ATPSink=1)
                                    else:
                                        if ConstantG6P:
                                            self.Simulation.SetRestoreSelectiveMolecules(["G6P"])
                                            self.Datasets[Subtitle] = self.Simulation.Sim(Glycolysis=1, PyruvateOxidation=1, TCACycle=1, OxidativePhosphorylation=1, ATPSink=1)
                                        else:
                                            self.Simulation.SetGlucoseAvailabilityType(GlucoseAvailabilityType)
                                            self.Simulation.SetGlucoseTransportType(GlucoseTransportType)
                                            self.Datasets[Subtitle] = self.Simulation.Sim(Chemotaxis=1, Glycolysis=1, PyruvateOxidation=1, TCACycle=1, OxidativePhosphorylation=1, ATPSink=1)

        # spacing = "\n"
        spacing = ", "

        PlottingPathway = "ATP/ADP balance"
        InclusionList = ["ATP", "ATP_Sink", "ADP"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit,
                           SuperTitle=(ModelName + spacing + PlottingPathway))

    def Model_GlucoseMetabolism_CellDivision_KFinder(self):
        self.InitializeDatasets()
        ModelName = "Glucose Metabolism + CellDivision + KFinder (" + TargetK + ")"

        # GlucoseAvailabilityTypes = ["Constant_GlucoseAvailability", "Increasing_GlucoseAvailability", "Fluctuating_GlucoseAvailability"]
        # GlucoseAvailabilityTypes = ["Increasing_GlucoseAvailability"]
        # GlucoseTransportTypes = ["G6P'=0.0001/(G6P*cgt)^2"]

        GlycolysisATPTypes = ["ATP'=(ADP/ATP)*cg"]
        PyruvateOxidationATPTypes = ["acetyl-CoA'=(pyruvate/acetyl-CoA)*cp"]
        TCACycleATPTypes = ["ATP'=ct"]
        OxidativePhosphorylationATPTypes = ["ATP'=co"]
        OxidativePhosphorylationCapacityTypess = ["co=min(FADH2,NADH)"]
        ATPSinkTypes = ["CellDivision_ATPSink"]
        # ATPSinkTypes = ["No_ATPSink"]

        Subtitle = "[Glyc]" + GlycolysisATPTypes[0] + ", [PO]" + PyruvateOxidationATPTypes[0] + "\n[TCA]" + TCACycleATPTypes[0] + ", [OP]" + OxidativePhosphorylationATPTypes[0] + "\n" + ATPSinkTypes[0]

        self.Simulation.SetGlycolysisATPType(GlycolysisATPTypes[0])
        self.Simulation.SetPyruvateOxidationATPType(PyruvateOxidationATPTypes[0])
        self.Simulation.SetTCACycleATPType(TCACycleATPTypes[0])
        self.Simulation.SetATPSinkType(ATPSinkTypes[0])
        self.Simulation.SetGlycolysisEndProduct("pyruvate")
        self.Simulation.SetAcetylCoALevel("Unfixed")
        self.Simulation.SetOxidativePhosphorylationATPType(OxidativePhosphorylationATPTypes[0])
        self.Simulation.SetOxidativePhosphorylationCapacityType(OxidativePhosphorylationCapacityTypess[0])

        self.Simulation.SetInitializeSelectiveMolecules(
            ["acetyl-CoA", "pyruvate", "NAD", "NADH", "FAD", "FADH2", "G6P", "CoA"])

        # Key Settings here
        # self.Simulation.SetRestoreSelectiveMolecules(["G6P", "CoA", "pyruvate"])
        self.Simulation.SetRestoreSelectiveMolecules(["G6P", "CoA"])
        KFinderFactors = [1 / (Base ** x) for x in range(1)]
        for KFinderFactor in KFinderFactors:
            self.Datasets[Subtitle + ", KFactor=" + str(KFinderFactor)] = self.Simulation.Sim(KFinderFactor=(TargetK, KFinderFactor), Glycolysis=1, PyruvateOxidation=1, TCACycle=1, OxidativePhosphorylation=1, ATPSink=1)

        self.RemoveEmptyDatasets()

        # spacing = "\n"
        spacing = ", "

        PlottingPathway = "(ATP/ADP balance)"
        InclusionList = ["ATP", "ATP_Sink", "ADP", "ATP_Sink", "ATP_Glycolysis", "ATP_TCA", "ATP_OxPhos"]
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))
        #
        # PlottingPathway = "ATP Contribution"
        # InclusionList = ["ATP_Sink", "ATP_Glycolysis", "ATP_TCA", "ATP_OxPhos"]
        # ExclusionList = []
        # self.Plot.SetFilters(InclusionList, ExclusionList)
        # self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))

        # Debugging
        PlottingPathway = "(Show all molecules)"
        InclusionList = []
        ExclusionList = []
        self.Plot.SetFilters(InclusionList, ExclusionList)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=(ModelName + spacing + PlottingPathway))

def main():
    # Initialization
    Simulation = FSimulation()
    Plot = FPlotter()

    # Model Runner
    ModelRunner = FModelRunner(Simulation, Plot)

    # Unit test Models
    # Models = [1]    # Chemotaxis unit test without G6PSink
    # Models = [2]    # Chemotaxis unit test with Linear G6PSink
    # Models = [3]    # Chemotaxis unit test with Burst G6PSink
    Models = [11]   # Glycolysis unit test without Pyruvate Oxidation
    # Models = [12]   # Glycolysis unit test with Pyruvate Oxidation
    # Models = [21]   # Pyruvate Oxidation unit test
    # Models = [31]   # TCACycle Cycle unit test without acetyl-CoA
    # Models = [32]   # TCACycle Cycle unit test with fixed acetyl-CoA
    # Models = [41]   # Oxidative Phosphorylation with Capacity parameter testing
    # Models = [42]   # Oxidative Phosphorylation with Feedback models
    # Models = [101]  # Chemotaxis + Glycolysis (end product: pyruvate)
    # Models = [102]  # Chemotaxis + Glycolysis (end product: acetyl-CoA)
    # Models = [103]  # Chemotaxis + Glycolysis + Pyruvate Oxidation
    # Models = [104]  # Chemotaxis + Glycolysis + Pyruvate Oxidation + TCA Cycle
    # Models = [105]  # Chemotaxis + Glycolysis + Pyruvate Oxidation + TCA Cycle + Oxidative Phosphorylation
    # Models = [1, 2, 3, 11, 12, 21]
    # Models = [31, 32, 41, 42]   # Need revision due to new capacity implementation for cell division
    # Models = [101, 102, 103, 104, 105]   # not optimized but all pathways function.
    TotalSimulationTime = 1e-1   # s
    SimulationTimeUnit = 1e-6   # s

    # Cell Division Models
    # Models = [201]  # Glucose Metabolism + Cell Division (With Constant Byproducts)
    # Models = [202]  # Glucose Metabolism + Cell Division (With Constant G6P)
    # Models = [203]  # Glucose Metabolism + Cell Division (With Chemotaxis)
    # Models = [211]  # Glucose Metabolism + Cell Division (With KFinder)

    # Setting Simulation Parameters
    for Model in Models:
        if Model > 200:
            TotalSimulationTime = 60 * 40  # s
            SimulationTimeUnit = 0.01  # s
            break

    # SteadyStateBrake = False
    # SteadyStateCheckMolecule = "ATP"
    # SteadyStateThresholdFactor = 0.0000001  # Steady state threshold

    Simulation.SetSimulationParameters(TotalSimulationTime, SimulationTimeUnit)
    Simulation.PrintSimulationParameters()
    ModelRunner.RunModel(Models)

if __name__ == '__main__':
    main()