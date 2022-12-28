'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Goodwin Oscillator
Pulse Generator

'''

import random
import sys
import math
import matplotlib.pyplot as plt
from src.lcc.plot2 import FPlotter

'''Figure 7.17: The Goodwin oscillator. A. This simulation shows relaxation to sustained (limit cycle)
oscillations. B. A phase portrait showing convergence to a limit cycle in the three-dimensional phase space.

Parameter values are a = 360 (concentration · time−1), k = 1.368 (concentration), b = 1 (time−1), alpha = 1
(time−1), beta = 0.6 (time−1), gamma = 1 (time−1), delta = 0.8 (time−1), n = 12. Units are arbitrary.'''

random.seed(1)

# uM = 1e-6
# nM = 1e-9
# 
# Unit = nM

Mol2Count = 6.0221409e+23
Count2Mol = 1.0 / Mol2Count

# Debugging Switch
DebugPrint = True
# DebugPrint = False

# DebugZero = True
DebugZero = False

class FSimulation():
    def __init__(self):
        # Molecular Parameters
        self.MolConc = dict()   # current concentration
        self.DeltaMolConc = dict()
        self.Dataset = dict()   # storing concentration
        self.Const = dict()
        self.Conjugation = dict()

        # Simulation Parameters
        self.SimStep = 0
        self.TotalSimulationTime = 0
        self.SimulationTimeUnit = 0
        self.SteadyStateBrake = False
        self.SteadyStateCheckMolecule = None
        self.SteadyStateThresholdFactor = 0

        # Model Parameters
        self.Title = ""
        self.RestoreSelectiveMolecules = dict()

    def SetSimulationParameters(self, TotalSimulationTime, SimulationTimeUnit, SteadyStateBrake=False, SteadyStateCheckMolecule=None, SteadyStateThresholdFactor=None):
        self.TotalSimulationTime = TotalSimulationTime
        self.SimulationTimeUnit = SimulationTimeUnit
        self.SteadyStateBrake = SteadyStateBrake
        self.SteadyStateCheckMolecule = SteadyStateCheckMolecule
        self.SteadyStateThresholdFactor = SteadyStateThresholdFactor
        self.PrintSimulationParameters()

    def PrintSimulationParameters(self):
        print("Total Simulation Time :", self.TotalSimulationTime, "s", " | Simulation Time Unit :", self.SimulationTimeUnit, "s", end="")

    def PrintSimulationSubject(self):
        print("\n\nSimulating " + self.Title + "...")

    def InitializeMolConcConst(self, InitData):
        MolConc, Const = InitData
        for MolName, Conc in MolConc.items():
            self.MolConc[MolName] = Conc
        for ConstName, Value in Const.items():
            self.Const[ConstName] = Value

    def InitializeDeltaMolConc(self):
        for MolName, Conc in self.MolConc.items():
            self.DeltaMolConc[MolName] = 0

    def InitializeDataset(self):
        for MolName, Conc in self.MolConc.items():
            self.Dataset[MolName] = [Conc]
        # for PerturbationName, Conc in self.Perturbation.items():
        #     self.Dataset[PerturbationName] = [Conc]

    def InitializeConjugation(self):
        pass
    
    def InitializeSimStep(self):
        self.SimStep = 0

    def SetRestoreSelectiveMolecules(self, MoleculesToRestore):
        for Molecule, Quantity in MoleculesToRestore.items():
            self.RestoreSelectiveMolecules[Molecule] = Quantity

    def AppendData(self):
        for MolName, Conc in self.MolConc.items():
            self.Dataset[MolName].append(Conc)
        # for PerturbationName, Conc in self.Perturbation.items():
        #     self.Dataset[PerturbationName].append(Conc)

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

    def SetTitle(self, Title):
        self.Title = Title
        self.PrintSimulationSubject()

    def Restore_SelectiveMolecules(self):
        for Molecule, Quantity in self.RestoreSelectiveMolecules:
            self.MolConc[Molecule] = Quantity

    def Initialize(self, InitData):
        self.InitializeMolConcConst(InitData)
        self.InitializeDeltaMolConc()
        self.InitializeDataset()
        self.InitializeConjugation()
        self.InitializeSimStep()

    def Run(self, GoodwinOscillator=0, PulseGeneration=0):

        while self.SimStep * self.SimulationTimeUnit < self.TotalSimulationTime:

            if GoodwinOscillator:
                self.Sim_GoodwinOscillator()
            if PulseGeneration:
                self.Sim_PulseGeneration()

            # Debugging
            if DebugPrint:
                self.PrintSimTimeAndMolConc()

            if self.RestoreSelectiveMolecules:
                self.Restore_SelectiveMolecules()

            self.CheckZeroSumConjugation()
            self.AppendData()

            bBelowZero = self.UpdateMolConc()

            if self.SteadyStateBrake:
                if (abs(self.Dataset[self.SteadyStateCheckMolecule][-1] -
                        self.Dataset[self.SteadyStateCheckMolecule][-2]) /
                        self.Dataset[self.SteadyStateCheckMolecule][-1] < self.SteadyStateThresholdFactor):
                    break
            self.SimStep += 1

        # Final Log
        self.PrintSimTimeAndMolConc()
        return self.Dataset

    def Sim_GoodwinOscillator(self):
        d = dict()

        X = self.MolConc["X"]
        Y = self.MolConc["Y"]
        Z = self.MolConc["Z"]

        a = self.Const["a"]
        b = self.Const["b"]
        k = self.Const["k"]
        n = self.Const["n"]
        alpha = self.Const["alpha"]
        beta = self.Const["beta"]
        gamma = self.Const["gamma"]
        delta = self.Const["delta"]

        def Get_dX(a, b, k, n, X, Z):
            return a / ((k ** n) + (Z ** n)) - b * X

        def Get_dY(alpha, beta, X, Y):
            return (alpha * X) - (beta * Y)

        def Get_dZ(gamma, delta, Y, Z):
            return (gamma * Y) - (delta * Z)

        # Production Rate Type
        d["X"] = Get_dX(a, b, k, n, X, Z) * self.SimulationTimeUnit
        d["Y"] = Get_dY(alpha, beta, X, Y) * self.SimulationTimeUnit
        d["Z"] = Get_dZ(gamma, delta, Y, Z) * self.SimulationTimeUnit

        self.AddToDeltaMolConc(d)

    def Sim_PulseGeneration(self):
        d = dict()

        A = self.MolConc["A"]
        R = self.MolConc["R"]
        C = self.MolConc["C"]
        G = self.MolConc["G"]
        RT = self.MolConc["RT"]

        k1 = self.Const["k1"]
        k2 = self.Const["k2"]
        aC = self.Const["aC"]
        bC = self.Const["bC"]
        aG = self.Const["aG"]
        bG = self.Const["bG"]
        KR = self.Const["KR"]
        KC = self.Const["KC"]

        def Get_dR(k1, A, RT, R, k2):
            return k1 * A ** 2 * (RT - 2 * R) ** 2 - k2 * R

        def Get_dC(aC, R, KR, bC, C):
            return aC * (R / KR) / (1 + (R / KR)) - bC * C

        def Get_dG(aG, R, KR, C, KC):
            return aG * (R / KR) / (1 + (R / KR) + (C / KC) ** 2 + (R / KR) * (C / KC) ** 2) - bG * G

        d["R"] = Get_dR(k1, A, RT, R, k2) * self.SimulationTimeUnit
        d["C"] = Get_dC(aC, R, KR, bC, C) * self.SimulationTimeUnit
        d["G"] = Get_dG(aG, R, KR, C, KC) * self.SimulationTimeUnit

        self.AddToDeltaMolConc(d)


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
                self.Model_GoodwinOscillator()
            elif Model == 2:
                self.Model_PulseGeneration()
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

    def Initialize_GoodwinOscillator(self):
        # a=360, b=1, k=1.368, n=12, alpha=1, beta=0.6, gamma=1, delta=0.8, X=0, Y=0, Z=0

        MolConc = dict()
        MolConc["X"] = 0
        MolConc["Y"] = 0
        MolConc["Z"] = 0

        Const = dict()
        Const["a"] = 360
        Const["b"] = 1
        Const["k"] = 1.368
        Const["n"] = 12
        Const["alpha"] = 1
        Const["beta"] = 0.6
        Const["gamma"] = 1
        Const["delta"] = 0.8

        return MolConc, Const

    def Model_GoodwinOscillator(self):
        TotalSimulationTime = 50  # s
        SimulationTimeUnit = 1e-1  # s
        self.Simulation.SetSimulationParameters(TotalSimulationTime, SimulationTimeUnit)

        self.InitializeDatasets()
        
        ModelName = "Goodwin Oscillator"
        self.Simulation.SetTitle(ModelName)

        InitData = self.Initialize_GoodwinOscillator()
        self.Simulation.Initialize(InitData)
        
        self.Datasets[ModelName] = self.Simulation.Run(GoodwinOscillator=1)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

    def Initialize_PulseGeneration(self):
        # From the text:
        # k1 = 0.5 (min−1 · concentration−3), k2 = 0.02 (min−1 ·concentration−3), RT = 0.5 (concentration), aC = 0.5
        # (concentration ·min−1) KR = 0.02 (concentration), bC = 0.07 (min−1), aG = 80 (concentration ·min−1),
        # KC = 0.008 (concentration), bG = 0.07 (min−1). Concentration units are arbitrary.
        # page 269, caption of Figure 7.28, value of b_C should be 0.3 (not 0.07 as printed)

        MolConc = dict()

        MolConc["A"] = 10    # AHL
        MolConc["R"] = 0  # AHL-LuxR complex
        MolConc["C"] = 0    # cI
        MolConc["G"] = 0 # GFP
        MolConc["RT"] = 0.5 # active LuxR-AHL complexes

        Const = dict()
        Const["k1"] = 0.5
        Const["k2"] = 0.02
        Const["aC"] = 0.5
        Const["bC"] = 0.3
        Const["aG"] = 80
        Const["bG"] = 0.07
        Const["KR"] = 0.02
        Const["KC"] = 0.008

        return MolConc, Const

    def Model_PulseGeneration(self):
        TotalSimulationTime = 50  # s
        SimulationTimeUnit = 1e-2  # s
        self.Simulation.SetSimulationParameters(TotalSimulationTime, SimulationTimeUnit)

        self.InitializeDatasets()

        ModelName = "Pulse Generation"
        self.Simulation.SetTitle(ModelName)

        InitData = self.Initialize_PulseGeneration()
        self.Simulation.Initialize(InitData)

        self.Datasets[ModelName] = self.Simulation.Run(PulseGeneration=1)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

def main():
    # Initialization
    Simulation = FSimulation()
    Plot = FPlotter()
    ModelRunner = FModelRunner(Simulation, Plot)

    # Unit test Models
    # Models = [1]    # Goodwin Oscillator
    Models = [2]    # Pulse Generation

    ModelRunner.RunModel(Models)

if __name__ == '__main__':
    main()