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
                 Perturbation = {},
                 UserSetInitialMolecules = {}):
        super().__init__()

        self.Sim = Metabolism.ReactionSimulator()
        Metabolism.EcoliInfo.Info()
        ATPConsumption_Sec = Metabolism.EcoliInfo.ECM_CellDivision_Sec

        # self.Debug_Info = 1

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
        self.ProteinSynthesis = Metabolism.ProteinSynthesis(ProteinSynthesisRate)
        self.Sim.AddReaction(self.ProteinSynthesis)

        self.Sim.SetPermanentMolecules(PermanentMolecules)
        self.Sim.SetPerturbation(Perturbation)
        self.Sim.Initialize(UserSetInitialMolecules)

        self.TotalTime = self.Sim.Molecules["G6P"] * 32 / max(1e-3, ATPConsumption_Sec) + 200
        self.DeltaTime = 0.01

    def SimulateDelta(self, DeltaTime = 0.01):
        self.Sim.Iter = self.Iter
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
    KnownMolConc = Metabolism.EcoliInfo.OpenKnownMolConc()

    Sim = FEcoliSimulator(
        PermanentMolecules = [
            "G6P",
        ],
        Perturbation = {   # Set Perturbation (time: {mol, conc})
            # 50  : {
            #     "PfkA": KnownMolConc["PfkA"][0] * 0.02,
            #     "AceE": KnownMolConc["AceE"][0] * 0.02,
            # },
            # 150 : {
            #     "PfkA": KnownMolConc["PfkA"][0] * 1,
            #     "AceE": KnownMolConc["AceE"][0] * 0.2,
            # },
        }
    )

    Sim.PrintReactions()

    TotalTime = Sim.TotalTime
    DeltaTime = Sim.DeltaTime
    Sim.Simulate(TotalTime = TotalTime,
                 DeltaTime = DeltaTime)

    Sim.Plot = True
    if Sim.Plot:
        Datasets = Sim.GetDataset()
        Plot = FPlotter()

        Plot.SetKnownMolConc(Metabolism.EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Unitless=False, Multiscale=True, Export='')

        Plot.SetKnownMolConc(Metabolism.EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Unitless=False, All=False, Individual=True, MolRange=True)
