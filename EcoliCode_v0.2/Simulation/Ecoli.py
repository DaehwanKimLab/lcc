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
import Metabolism

class FEcoliSimulator(FSimulator):
    DNAREPLICATIONRATE = 1.0 / (16 * 60)
    PROTEINSYNTHESISRATE = 1.0 / (16 * 60)
    CYTOKINESISRATE = 1.0 / ( 6 * 60)

    def __init__(self, 
                 DNAReplicationRate = DNAREPLICATIONRATE,
                 ProteinSynthesisRate = PROTEINSYNTHESISRATE,
                 CytoKinesisRate = CYTOKINESISRATE):
        super().__init__()

        self.DNAReplicationRate = DNAReplicationRate
        self.DNAReplicationProgress = 0.0
        self.ProteinSynthesisRate = ProteinSynthesisRate
        self.ProteinSynthesisProgress = 0.0
        self.CytoKinesisRate = CytoKinesisRate
        self.CytoKinesisProgress = 0.0

        self.NumCellDivisions = 0

        self.Sim = Metabolism.ReactionSimulator()
        Metabolism.EcoliInfo.Info()
        ATPConsumption_Sec = Metabolism.EcoliInfo.ECM_CellDivision_Sec

        # self.Sim.AddReaction(Metabolism.NADPlusSynthesis())
        # self.Sim.AddReaction(Metabolism.NADPPlusSynthesis())
        # self.Sim.AddReaction(Metabolism.CoASynthesis())
        self.Sim.AddReaction(Metabolism.Glycolysis())
        self.Sim.AddReaction(Metabolism.PyruvateOxidation())
        self.Sim.AddReaction(Metabolism.TCACycle())
        self.Sim.AddReaction(Metabolism.NADH_OxidativePhosphorylation())
        self.Sim.AddReaction(Metabolism.FADH2_OxidativePhosphorylation())
        self.Sim.AddReaction(Metabolism.ATPControl(-ATPConsumption_Sec))

        # Set permanent molecules
        PermanentMolecules = [
            # "G6P",
            # "pyruvate",
            # "CoA-SH",
            # "NADH",
            # "NAD+",
            # "FADH2",
            # "FAD",
            # "CoA-SH",
            # "oxaloacetate",
        ]
        self.Sim.SetPermanentMolecules(PermanentMolecules)

        # Debugging options
        # Sim.Debug_Reaction = True
        self.Sim.Debug_Info = 100
        self.Sim.Plot = True
        
        # Set initial molecule concentrations
        Molecules = {}
        # Molecules["ATP"] = 1.0 * 1e-3
        # Molecules["ADP"] = 8.6 * 1e-3
        # Molecules["G6P"] = 50 * 1e-3
        
        # Execute simulation
        self.Sim.Initialize(Molecules)

        # Simulation parameters
        self.TotalTime = self.Sim.Molecules["G6P"] * 32 / max(1e-3, ATPConsumption_Sec) + 200
        self.DeltaTime = 0.01

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

        # self.Sim.SimulateDelta(DeltaTime)

    def SetProgress(self, Progress):
        self.DNAReplicationProgress = Progress
        self.ProteinSynthesisProgress = Progress
        self.CytoKinesisProgress = 0.0

    def GetNumCellDivisions(self):
        return self.NumCellDivisions

    def Info(self):
        print("DNA Replication Progress:   {:<.3f}".format(self.DNAReplicationProgress))
        print("Protien Synthesis Progress: {:<.3f}".format(self.ProteinSynthesisProgress))
        print("CytoKinesis Progress:       {:<.3f}".format(self.CytoKinesisProgress))

    
if __name__ == '__main__':
    Sim = FEcoliSimulator()
    
    TotalTime = 22 * 60.0
    DeltaTime = 0.01

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



