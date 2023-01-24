# BSD 3-Clause License
# © 2022, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import sys
from datetime import datetime
import random, copy
import numpy as np
import math
import matplotlib.pyplot as plt


class FSimulator():
    def __init__(self):
        None

    def Initialize(self):
        None

    def Simulate(self, TotalTime = 10, DeltaTime = 0.01):
        None

    def GetDataset(self):
        return None


class FEcoliSimulator(FSimulator):
    DNAREPLICATIONRATE = 1.0 / (16 * 60)
    PROTEINSYNTHESISRATE = 1.0 / (16 * 60)
    CYTOKINESISRATE = 1.0 / ( 6 * 60)

    def __init__(self, 
                 DNAReplicationRate = DNAREPLICATIONRATE,
                 ProteinSynthesisRate = PROTEINSYNTHESISRATE,
                 CytoKinesisRate = CYTOKINESISRATE):
        self.DNAReplicationRate = DNAReplicationRate
        self.DNAReplicationProgress = 0.0
        self.ProteinSynthesisRate = ProteinSynthesisRate
        self.ProteinSynthesisProgress = 0.0
        self.CytoKinesisRate = CytoKinesisRate
        self.CytoKinesisProgress = 0.0

        self.NumCellDivisions = 0

    def SimulateDelta(self, DeltaTime = 1.0):
        if self.DNAReplicationProgress >= 1.0 and \
           self.ProteinSynthesisProgress >= 1.0:
            self.CytoKinesisProgress += self.CytoKinesisRate * DeltaTime
        else:
            if self.DNAReplicationProgress < 1.0:
                self.DNAReplicationProgress += self.DNAReplicationRate * DeltaTime
            if self.ProteinSynthesisProgress < 1.0:
                self.ProteinSynthesisProgress += self.ProteinSynthesisRate * DeltaTime

        self.NumCellDivisions = 0
        if self.CytoKinesisProgress >= 1.0:
            self.DNAReplicationProgress = 0.0
            self.ProteinSynthesisProgress = 0.0
            self.CytoKinesisProgress = 0.0
            self.NumCellDivisions = 1 

    def SetProgress(self, Progress):
        self.DNAReplicationProgress = Progress
        self.ProteinSynthesisProgress = Progress
        self.CytoKinesisProgress = 0.0

    def GetNumCellDivisions(self):
        return self.NumCellDivisions

    
if __name__ == '__main__':
    Sim = FEcoliSimulator()
    TotalTime = 6 * 60.0 * 60.0
    DeltaTime = 1.0

    Sim.Initialize()
    Sim.Simulate(TotalTime=TotalTime, DeltaTime=DeltaTime)

    """
    Sim.Plot = True
    if Sim.Plot:
        Datasets = Sim.GetDataset()
        Plot = FPlotter()
        Plot.PlotDatasets(Datasets, DeltaTime=DeltaTime)

        GrowthDatasets = Sim.GetGrowthDataset()
        GrowthPlot = FGrowthPlotter()
        GrowthPlot.PlotDatasets(GrowthDatasets)
    """



