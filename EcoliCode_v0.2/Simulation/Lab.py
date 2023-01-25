# BSD 3-Clause License
# © 2023, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import sys
from datetime import datetime
import random, copy
import numpy as np
import math
import matplotlib.pyplot as plt

from Simulator import FSimulator
from Ecoli import FEcoliSimulator
import Plotter

class FColonySimulator(FSimulator):
    def __init__(self, EcoliSim):
        self.EcoliSims = [[1, EcoliSim]]
        self.DataPoints = []

    def SimulateDelta(self, DeltaTime = 1.0):
        NewEcoliSims = []
        for i, (NumEcoli, Sim) in enumerate(self.EcoliSims):
            Sim.SimulateDelta(DeltaTime)
            NumCellDivisions = Sim.GetNumCellDivisions()
            if NumCellDivisions > 0:
                NumEcoli *= math.pow(2, NumCellDivisions)
            self.EcoliSims[i][0] = NumEcoli
            if len(self.EcoliSims) < 20:
                if NumEcoli > 1:
                    NumEcoli1 = int(NumEcoli / 2)
                    NumEcoli2 = NumEcoli - NumEcoli1
                    self.EcoliSims[i][0] = NumEcoli1
                    Sim2 = copy.deepcopy(Sim)
                    Sim2.SetProgress(random.uniform(-0.3, 0.3))
                    NewEcoliSims.append([NumEcoli2, Sim2])
        self.EcoliSims += NewEcoliSims
                    
        self.AddToDataset()

    def Info(self):
        None

    def AddToDataset(self):
        TotalNumEcoli = 0
        for NumEcoli, _ in self.EcoliSims:
            TotalNumEcoli += NumEcoli
        self.DataPoints.append(TotalNumEcoli)
        
    def GetDataset(self):
        return self.DataPoints

    def GetNumEcoli(self):
        return self.DataPoints[-1] if len(self.DataPoints) > 0 else 0


class FPopulationSimulator(FSimulator):
    DEFAULT_COLONY = [
        FEcoliSimulator.DNAREPLICATIONRATE,
        FEcoliSimulator.PROTEINSYNTHESISRATE,
        FEcoliSimulator.CYTOKINESISRATE
    ]

    DEFAULT_POPULATION = [DEFAULT_COLONY]

    TEST1_POPULATION = [
        DEFAULT_COLONY,
        [
            FEcoliSimulator.DNAREPLICATIONRATE,
            FEcoliSimulator.PROTEINSYNTHESISRATE * 0.5,
            FEcoliSimulator.CYTOKINESISRATE
        ],
        [
            FEcoliSimulator.DNAREPLICATIONRATE,
            FEcoliSimulator.PROTEINSYNTHESISRATE,
            FEcoliSimulator.CYTOKINESISRATE * 0.5
        ]
    ]

    def __init__(self, Population):
        self.Colonies = [
            FColonySimulator(FEcoliSimulator(Colony[0], Colony[1], Colony[2])) for Colony in Population
        ]
        self.Dataset = {}

    def Simulate(self, TotalTime = 24 * 60.0 * 60.0, DeltaTime = 1.0):
        Iter = 0
        while Iter < TotalTime / DeltaTime:
            for ColonySim in self.Colonies:
                ColonySim.SimulateDelta(DeltaTime)

            self.AddToDataset()

            Iter += 1
            if Iter % 100 == 0:
                print("\n")
                print("-- Iteration {} --".format(Iter))
                self.Info()

    def Info(self):
        for i, Colony in enumerate(self.Colonies):
            print("Colony {:>2}: {:>10}".format(i, Colony.GetNumEcoli()))

    def AddToDataset(self):
        for i, Colony in enumerate(self.Colonies):
            if i not in self.Dataset:
                self.Dataset[i] = []
            self.Dataset[i].append(Colony.GetNumEcoli())
        
    def GetDataset(self):
        PlotDataset = {}
        TotalData = None
        for i, Data in self.Dataset.items():
            Name = "Colony {:>2}".format(i)
            PlotDataset[Name] = Data
            if TotalData:
                TotalData = [a + b for a, b in zip(TotalData, Data)]
            else:
                TotalData = Data
        PlotDataset["Population"] = TotalData        

        return {"E. coli Growth": PlotDataset}


class FPetridishSimulator(FSimulator):
    def __init__(self, Population):
        self.PopSim = FPopulationSimulator(Population)

    def Simulate(self, TotalTime = 10, DeltaTime = 0.01):
        self.PopSim.Simulate(TotalTime, DeltaTime)

    def GetDataset(self):
        return self.PopSim.GetDataset()


class FExperimentSimulator(FSimulator):
    def __init__(self):
        # self.Control = FPetridishSimulator(FPopulationSimulator.DEFAULT_POPULATION)
        self.Control = FPetridishSimulator(FPopulationSimulator.TEST1_POPULATION)
        self.NumTest = 30
        self.Tests = [] 
        for i in range(self.NumTest):
            self.Tests.append(FPetridishSimulator([[
                FEcoliSimulator.DNAREPLICATIONRATE,
                FEcoliSimulator.PROTEINSYNTHESISRATE * (self.NumTest - max(0, i - self.NumTest / 1.2)) / self.NumTest,
                FEcoliSimulator.CYTOKINESISRATE
            ]]))

    def Simulate(self, TotalTime = 10, DeltaTime = 0.01):
        self.Control.Simulate(TotalTime, DeltaTime)
        for Test in self.Tests:
            Test.Simulate(TotalTime, DeltaTime)
                              
    def GetDataset(self):
        return self.Control.GetDataset()

    def GetGrowthDataset(self):
        XSet, YSet = [], []
        ControlDataset = list(self.Control.GetDataset().values())[0]
        ControlDataset = ControlDataset["Population"]
        MaxGrowth = ControlDataset[-1]
        for i, Test in enumerate(self.Tests):
            TestDataset = list(Test.GetDataset().values())[0]
            TestDataset = TestDataset["Population"]
            TestGrowth = TestDataset[-1]
            XSet.append(i / self.NumTest)
            YSet.append(TestGrowth / MaxGrowth)

        return {"E. coli Growth": [XSet, YSet]}


if __name__ == '__main__':
    Sim = FExperimentSimulator()
    Sim.Plot = True

    TotalTime = 6 * 60.0 * 60.0
    DeltaTime = 1.0

    Sim.Initialize()
    Sim.Simulate(TotalTime=TotalTime, DeltaTime=DeltaTime)

    if Sim.Plot:
        Datasets = Sim.GetDataset()
        Plot = Plotter.FGrowthPlotter()
        Plot.PlotDatasets(Datasets, DeltaTime=DeltaTime, Unitless=True)

        GrowthDatasets = Sim.GetGrowthDataset()
        GrowthPlot = Plotter.FGeneRepressionPlotter()
        GrowthPlot.PlotDatasets(GrowthDatasets)



