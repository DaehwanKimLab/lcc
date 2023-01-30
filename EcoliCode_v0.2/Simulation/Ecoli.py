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
    DNAREPLICATIONRATE = Metabolism.EcoliInfo.DNAReplicationRate
    PROTEINSYNTHESISRATE = Metabolism.EcoliInfo.ProteinSynthesisRate
    CYTOKINESISRATE = 1.0 / (6 * 60)

    def __init__(self, 
                 DNAReplicationRate = DNAREPLICATIONRATE,
                 ProteinSynthesisRate = PROTEINSYNTHESISRATE,
                 CytoKinesisRate = CYTOKINESISRATE,
                 PermanentMolecules = [],
                 InitialMolecules = {}):
        super().__init__()

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
        # self.Sim.AddReaction(Metabolism.ATPControl(-ATPConsumption_Sec))

        self.Sim.AddReaction(Metabolism.dNTPSynthesis(DNAReplicationRate))
        self.DNAReplication = Metabolism.DNAReplication(DNAReplicationRate)
        self.Sim.AddReaction(self.DNAReplication)

        self.Sim.AddReaction(Metabolism.AASynthesis(ProteinSynthesisRate))
        # DK - debugging purposes
        self.ProteinSynthesis = Metabolism.ProteinSynthesis(ProteinSynthesisRate)
        self.Sim.AddReaction(self.ProteinSynthesis)

        self.Sim.SetPermanentMolecules(PermanentMolecules)
        self.Sim.Initialize(InitialMolecules)

        # Simulation parameters
        self.TotalTime = self.Sim.Molecules["G6P"] * 32 / max(1e-3, ATPConsumption_Sec) + 200
        self.DeltaTime = 0.01

    def SimulateDelta(self, DeltaTime = 0.01):
        self.Sim.SimulateDelta(DeltaTime)

    def GetProgress(self):
        return self.DNAReplication.GetProgress()

    def SetProgress(self, Progress):
        self.DNAReplication.SetProgress(Progress)
        self.ProteinSynthesis.SetProgress(Progress)

    def GetDataset(self):
        return self.Sim.GetDataset()

    def KnownMolConc(self):
        return self.Sim.KnownMolConc

    def PrintReactions(self):
        self.Sim.PrintReactions()

    def Info(self):
        print("DNA Replication Progress:   {:<.3f}".format(self.DNAReplication.GetProgress()))
        print("Protein Synthesis Progress:   {:<.3f}".format(self.ProteinSynthesis.GetProgress()))
        print()

        self.Sim.Info()

    
if __name__ == '__main__':
    PermanentMolecules = [
        # "G6P",
    ]

    InitialMolecules = {
        # "ATP": 1.0 * 1e-3,
    }

    Sim = FEcoliSimulator(PermanentMolecules = PermanentMolecules,
                          InitialMolecules = InitialMolecules)

    Sim.PrintReactions()

    TotalTime = Sim.TotalTime
    DeltaTime = Sim.DeltaTime
    Sim.Simulate(TotalTime = TotalTime,
                 DeltaTime = DeltaTime)

    Sim.Plot = True
    if Sim.Plot:
        Datasets = Sim.GetDataset()
        KnownMolConc = Sim.KnownMolConc()
        Plot = FPlotter()

        Plot.SetKnownMolConc(copy.deepcopy(KnownMolConc))
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Unitless=False, Multiscale=True, Export='')
        Plot.SetKnownMolConc(copy.deepcopy(KnownMolConc))
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Unitless=False, All=False, Individual=True, MolRange=True)
