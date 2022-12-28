'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Goodwin Oscillator
Network  7.16,   p.213
Model   7.22,    p.212
Figure  7.17,    p.214

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

    def Run(self, GoodwinOscillator=0):

        while self.SimStep * self.SimulationTimeUnit < self.TotalSimulationTime:

            if GoodwinOscillator:
                self.Sim_GoodwinOscillator()

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
        
        ModelName = "GoodwinOscillator"
        self.Simulation.SetTitle(ModelName)

        InitData = self.Initialize_GoodwinOscillator()
        self.Simulation.Initialize(InitData)
        
        self.Datasets[ModelName] = self.Simulation.Run(GoodwinOscillator=1)
        self.Plot.PlotDatasets(self.Datasets, self.Simulation.SimulationTimeUnit, SuperTitle=ModelName)

def main():
    # Initialization
    Simulation = FSimulation()
    Plot = FPlotter()
    ModelRunner = FModelRunner(Simulation, Plot)

    # Unit test Models
    Models = [1]    # Chemotaxis unit test without G6PSink
    ModelRunner.RunModel(Models)

if __name__ == '__main__':
    main()