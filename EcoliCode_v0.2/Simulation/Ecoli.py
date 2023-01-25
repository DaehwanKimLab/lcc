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
from Plotter import FPlotter

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

        self.Sim.AddReaction(Metabolism.DNAReplication())

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

        # Set initial molecule concentrations
        Molecules = {
            # "ATP": 1.0 * 1e-3,
            # "ADP": 8.6 * 1e-3,
            # "G6P": 50  * 1e-3,
        }
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

        # DK - temporary workaround
        if DeltaTime <= 0.02:
            self.SimulateDelta_Detail(DeltaTime)


    def SimulateDelta_Detail(self, DeltaTime):
        self.Sim.SimulateDelta(DeltaTime)

    def SetProgress(self, Progress):
        self.DNAReplicationProgress = Progress
        self.ProteinSynthesisProgress = Progress
        self.CytoKinesisProgress = 0.0

    def GetNumCellDivisions(self):
        return self.NumCellDivisions

    def GetDataset(self):
        return self.Sim.GetDataset()

    def KnownMolConc(self):
        return self.Sim.KnownMolConc

    def Info(self):
        print("DNA Replication Progress:   {:<.3f}".format(self.DNAReplicationProgress))
        print("Protien Synthesis Progress: {:<.3f}".format(self.ProteinSynthesisProgress))
        print("CytoKinesis Progress:       {:<.3f}".format(self.CytoKinesisProgress))

        self.Sim.Info()

        for Reaction in self.Sim.Reactions:
            Progress = Reaction.GetProgress()
            if Progress > 0.0:
                print("{}: {:<.3f}".format(Reaction.ReactionName, Progress))

    
if __name__ == '__main__':
    Sim = FEcoliSimulator()
    
    TotalTime = Sim.TotalTime
    DeltaTime = Sim.DeltaTime

    Sim.Initialize()
    Sim.Simulate(TotalTime=TotalTime, DeltaTime=DeltaTime)

    Sim.Plot = True
    if Sim.Plot:
        Datasets = Sim.GetDataset()
        KnownMolConc = Sim.KnownMolConc()
        Plot = FPlotter()

        Plot.SetKnownMolConc(copy.deepcopy(KnownMolConc))
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, All=True, Export='pdf')

        Plot.SetKnownMolConc(copy.deepcopy(KnownMolConc))
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, Multiscale=True, Export='pdf')

        Plot.SetKnownMolConc(copy.deepcopy(KnownMolConc))
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, Individual=True, MolRange=True, Export='pdf')
