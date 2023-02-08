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
                 Perturbations = {},
                 PermanentMolecules = [],
                 UserSetInitialMolecules = {},
                 ):
        super().__init__()

        self.Sim = Metabolism.ReactionSimulator()
        Metabolism.EcoliInfo.Info()
        ATPConsumption_Sec = Metabolism.EcoliInfo.ECM_CellDivision_Sec
        DNAReplicationRateInCount = DNAReplicationRate
        DNAReplicationRateInConc = DNAReplicationRate * Metabolism.EcoliInfo.C2M

        # self.Debug_Info = 1

        # Select Pathways
        Pathways = [
            "Central Carbon Metabolism",
            # "ATP Consumption by Cell Division",
            "Folic Acid",
            "Purine Metabolism",
            "DNA Replication",
            "Protein Synthesis",
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
        folA_ExpressionFactor = 1.0
        thyA_ExpressionFactor = 1.0
        glyA_ExpressionFactor = 1.0
        if "folA" in Perturbations or "thyA" in Perturbations or "glyA" in Perturbations:
            ExcludedMolecules.append("THF")
            ExcludedMolecules.append("5-methyl-THF")
            folA_ExpressionFactor = Perturbations["folA"]
            thyA_ExpressionFactor = Perturbations["thyA"]
            glyA_ExpressionFactor = Perturbations["glyA"]
            UserSetInitialMolecules["dTTP"] = 0.01e-3

        # Otto et al., 2022 - Repression of purN and purL
        purN_ExpressionFactor = 1.0
        purL_ExpressionFactor = 1.0

        if "purN" in Perturbations:
            ExcludedMolecules.append("GAR")
            ExcludedMolecules.append("FGAR")
            ExcludedMolecules.append("FGAM")
            purN_ExpressionFactor = Perturbations["purN"]
            purL_ExpressionFactor = Perturbations["purL"]

        if ("Folic Acid" in Pathways) and ("Purine Metabolism" not in Pathways):

            self.Sim.AddReaction(Metabolism.DHFSynthesis())
            self.Sim.AddReaction(Metabolism.THFSynthesisByFolA(Rate = DNAReplicationRateInConc, ExpressionFactor = folA_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.FiveMethylTHFSynthesisByGlyA(DNAReplicationRateInConc))
            self.Sim.AddReaction(Metabolism.dTTPSynthesisByThyA(DNAReplicationRateInCount))
            self.Sim.AddReaction(Metabolism.dUTPSynthesis(DNAReplicationRateInCount))
            self.DNAReplication = Metabolism.DNAReplication(DNAReplicationRate, BuildingBlocks=["dTTP"])

        elif ("Folic Acid" not in Pathways) and ("Purine Metabolism" in Pathways):

            assert False
            self.Sim.AddReaction(Metabolism.PRPPSynthesis(Rate = 3e-5))
            self.Sim.AddReaction(Metabolism.GARSynthesis(Rate = 3e-5))
            self.Sim.AddReaction(Metabolism.FGARSynthesisByPurN(Rate = 3e-5, ExpressionFactor = purN_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.FGAMSynthesisByPurL(Rate = 3e-5, ExpressionFactor = purL_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.PurineSynthesis(DNAReplicationRate * 0.5))  # dATP, dGTP
            self.DNAReplication = Metabolism.DNAReplication(DNAReplicationRate, BuildingBlocks=["dATP", "dGTP"])

        elif ("Folic Acid" in Pathways) and ("Purine Metabolism" in Pathways):    # overwrite with new rates and building blocks

            self.Sim.AddReaction(Metabolism.DHFSynthesis())
            self.Sim.AddReaction(Metabolism.THFSynthesisByFolA(Rate = DNAReplicationRateInConc, ExpressionFactor = folA_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.FiveTenMethyleneTHFSynthesisByGlyA(Rate = DNAReplicationRateInConc, ExpressionFactor = glyA_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.FiveMethylTHFSynthesis(DNAReplicationRateInConc))
            self.Sim.AddReaction(Metabolism.dUTPSynthesis(DNAReplicationRateInCount))
            self.Sim.AddReaction(Metabolism.dTTPSynthesisByThyA(Rate = DNAReplicationRateInCount * 0.25, ExpressionFactor = thyA_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.dCTPSynthesis(DNAReplicationRateInCount * 0.25))
            self.Sim.AddReaction(Metabolism.PurineSynthesis(DNAReplicationRateInCount * 0.25))  # dATP, dGTP
            self.Sim.AddReaction(Metabolism.PRPPSynthesis())
            self.Sim.AddReaction(Metabolism.GARSynthesis())
            self.Sim.AddReaction(Metabolism.FGARSynthesisByPurN(ExpressionFactor = purN_ExpressionFactor))
            self.Sim.AddReaction(Metabolism.FGAMSynthesisByPurL(ExpressionFactor = purL_ExpressionFactor))
            self.DNAReplication = Metabolism.DNAReplication(DNAReplicationRate, BuildingBlocks=["dATP", "dCTP", "dGTP", "dTTP"])

        else:
            self.Sim.AddReaction(Metabolism.dNTPSynthesis())  # dATP
            self.DNAReplication = Metabolism.DNAReplication(BuildingBlocks=["dATP"])

        # Processes
        if "DNA Replication" in Pathways:
            self.Sim.AddReaction(self.DNAReplication)

        if "Protein Synthesis" in Pathways:
            self.Sim.AddReaction(Metabolism.AASynthesis(ProteinSynthesisRate))
            self.ProteinSynthesis = Metabolism.ProteinSynthesis(ProteinSynthesisRate, BuildingBlocks=["glutamine"])
            self.Sim.AddReaction(self.ProteinSynthesis)

        self.VolumeExpansion = Metabolism.VolumeExpansion(
            ExcludedMolecules = ExcludedMolecules,
            Rate = Metabolism.EcoliInfo.VolumeExpansionRate
        )
        self.Sim.AddReaction(self.VolumeExpansion)

        self.Sim.SetPermanentMolecules(PermanentMolecules)
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
        self.VolumeExpansion.SetProgress(Progress)

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
        Perturbations = {
            # "folA": 0.1,
            # "purN": 0.02,
        },
        PermanentMolecules = [
            "G6P",
        ],
        UserSetInitialMolecules = {
        },
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

        """
        Plot.SetKnownMolConc(Metabolism.EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Unitless=False, All=False, Individual=True, MolRange=True)
        """
