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
                 UserSetInitialMolecules = {},
                 Perturbation = {}):
        super().__init__()

        self.Sim = Metabolism.ReactionSimulator()
        Metabolism.EcoliInfo.Info()
        ATPConsumption_Sec = Metabolism.EcoliInfo.ECM_CellDivision_Sec

        # self.Debug_Info = 1

        # Select Pathways
        Pathways = [
            "Central Carbon Metabolism",
            # "ATP Consumption by Cell Division",
            "Folic Acid",
            "Purine Metabolism",
            "DNA Replication",
            "Protein Synthesis",
            # "Volume Expansion",
        ]

        if "Central Carbon Metabolism" in Pathways:
            self.Sim.AddReaction(Metabolism.Glycolysis())
            self.Sim.AddReaction(Metabolism.PyruvateOxidation())
            self.Sim.AddReaction(Metabolism.TCACycle())
            self.Sim.AddReaction(Metabolism.NADH_OxidativePhosphorylation())
            self.Sim.AddReaction(Metabolism.FADH2_OxidativePhosphorylation())

        if "ATP Consumption by Cell Division" in Pathways:
            self.Sim.AddReaction(Metabolism.ATPControl(-ATPConsumption_Sec))

        # Cofactors
        # self.Sim.AddReaction(Metabolism.NADPlusSynthesis())
        # self.Sim.AddReaction(Metabolism.NADPPlusSynthesis())
        # self.Sim.AddReaction(Metabolism.CoASynthesis())

        ExcludedMolecules = list()

        # Otto et al., 2022 - Repression of folA
        Debug_folA = False
        folA_ExpressionFactor = 1.0

        if Debug_folA:
            ExcludedMolecules.append("THF")
            ExcludedMolecules.append("5-methyl-THF")
            folA_ExpressionFactor = 0.05
            # folA_ExpressionFactor = 0.0
            UserSetInitialMolecules["dTTP"] = 0.01e-3

        # Otto et al., 2022 - Repression of purN and purL
        Debug_purN_purL = False
        purN_ExpressionFactor = 1.0
        purL_ExpressionFactor = 1.0

        if Debug_purN_purL:
            ExcludedMolecules.append("GAR")
            ExcludedMolecules.append("FGAR")
            ExcludedMolecules.append("FGAM")
            purN_ExpressionFactor = 0.2
            purL_ExpressionFactor = 0.2

        if ("Folic Acid" in Pathways) and ("Purine Metabolism" not in Pathways):

            self.Sim.AddReaction(Metabolism.DHFSynthesis())
            self.Sim.AddReaction(Metabolism.THFSynthesisByFolA(ExpressionFactor = folA_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.FiveMethylTHFSynthesis())
            self.Sim.AddReaction(Metabolism.dTTPSynthesis(DNAReplicationRate))
            self.Sim.AddReaction(Metabolism.dUTPSynthesis(DNAReplicationRate))
            self.DNAReplication = Metabolism.DNAReplication(DNAReplicationRate, BuildingBlocks=["dTTP"])

        elif ("Folic Acid" not in Pathways) and ("Purine Metabolism" in Pathways):

            self.Sim.AddReaction(Metabolism.PRPPSynthesis(Rate = 3e-5))
            self.Sim.AddReaction(Metabolism.GARSynthesis(Rate = 3e-5))
            self.Sim.AddReaction(Metabolism.FGARSynthesisByPurN(Rate = 3e-5, ExpressionFactor = purN_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.FGAMSynthesisByPurL(Rate = 3e-5, ExpressionFactor = purL_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.PurineSynthesis(DNAReplicationRate * 0.5))  # dATP, dGTP
            self.DNAReplication = Metabolism.DNAReplication(DNAReplicationRate, BuildingBlocks=["dATP", "dGTP"])

        elif ("Folic Acid" in Pathways) and ("Purine Metabolism" in Pathways):    # overwrite with new rates and building blocks

            self.Sim.AddReaction(Metabolism.DHFSynthesis())
            self.Sim.AddReaction(Metabolism.THFSynthesisByFolA(ExpressionFactor=folA_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.FiveMethylTHFSynthesis())
            self.Sim.AddReaction(Metabolism.dUTPSynthesis(DNAReplicationRate * 0.25))
            self.Sim.AddReaction(Metabolism.dTTPSynthesis(DNAReplicationRate * 0.25))
            self.Sim.AddReaction(Metabolism.dCTPSynthesis(DNAReplicationRate * 0.25))
            self.Sim.AddReaction(Metabolism.PurineSynthesis(DNAReplicationRate * 0.25))  # dATP, dGTP
            self.Sim.AddReaction(Metabolism.PRPPSynthesis())
            self.Sim.AddReaction(Metabolism.GARSynthesis())
            self.Sim.AddReaction(Metabolism.FGARSynthesisByPurN(ExpressionFactor = purN_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.FGAMSynthesisByPurL(ExpressionFactor = purL_ExpressionFactor))
            self.DNAReplication = Metabolism.DNAReplication(DNAReplicationRate, BuildingBlocks=["dATP", "dCTP", "dGTP", "dTTP"])

        else:
            self.Sim.AddReaction(Metabolism.dNTPSynthesis(DNAReplicationRate))  # dATP
            self.DNAReplication = Metabolism.DNAReplication(DNAReplicationRate, BuildingBlocks=["dATP"])

        # Processes
        self.VolumeExpansion = Metabolism.VolumeExpansion(
            ExcludedMolecules = ExcludedMolecules,
            Rate = Metabolism.EcoliInfo.VolumeExpansionRate
        )
        self.Sim.AddReaction(self.VolumeExpansion)

        if "DNA Replication" in Pathways:
            self.Sim.AddReaction(self.DNAReplication)

        if "Protein Synthesis" in Pathways:
            self.Sim.AddReaction(Metabolism.AASynthesis(ProteinSynthesisRate))
            self.ProteinSynthesis = Metabolism.ProteinSynthesis(ProteinSynthesisRate, BuildingBlocks=["glutamine"])
            self.Sim.AddReaction(self.ProteinSynthesis)

        self.Sim.SetPermanentMolecules(PermanentMolecules)
        self.Sim.SetPerturbation(Perturbation)
        self.Sim.Initialize(UserSetInitialMolecules)

        self.TotalTime = self.Sim.Molecules["G6P"] * 32 / max(1e-3, ATPConsumption_Sec) + 200
        self.DeltaTime = 0.01

    def SimulateDelta(self, DeltaTime = 0.01):
        self.Sim.Iter = self.Iter
        self.Sim.SimulateDelta(DeltaTime)

        # temporary code - let's sync DNA replication and Volume expansion
        self.VolumeExpansion.SetProgress(self.DNAReplication.GetProgress())

    def GetProgress(self):
        return min(self.DNAReplication.GetProgress(),
                   self.ProteinSynthesis.GetProgress(),
                   self.VolumeExpansion.GetProgress())

    def SetProgress(self, Progress):
        self.DNAReplication.SetProgress(Progress)
        self.ProteinSynthesis.SetProgress(Progress)
        self.VolmeExpansion.SetProgress(Progress)

    def GetDataset(self):
        return self.Sim.GetDataset()

    def KnownMolConc(self):
        return self.Sim.KnownMolConc

    def PrintReactions(self):
        self.Sim.PrintReactions()

    def Info(self):
        print("DNA Replication Progress:   {:<.3f}".format(self.DNAReplication.GetProgress()))
        print("Protein Synthesis Progress: {:<.3f}".format(self.ProteinSynthesis.GetProgress()))
        print("Volume Expansion Progress:  {:<.3f}".format(self.VolumeExpansion.GetProgress()))
        print()

        self.Sim.Info()

    def Summary(self):
        self.Sim.Summary()
    
if __name__ == '__main__':
    KnownMolConc = Metabolism.EcoliInfo.OpenKnownMolConc()

    Sim = FEcoliSimulator(
        PermanentMolecules = [
            "G6P",
        ],
        UserSetInitialMolecules = {
        },
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
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Unitless=False, Multiscale=True, Export='', Log=10)

        # """
        Plot.SetKnownMolConc(Metabolism.EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Unitless=False, All=False, Individual=True, MolRange=True)
        # """
