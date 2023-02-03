# BSD 3-Clause License
# © 2023, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

"""
Reference

Metabolite concentrations, fluxes and free energies imply efficient enzyme usage
Junyoung O Park, Sara A Rubin, Yi-Fan Xu, Daniel Amador-Noguez, Jing Fan, Tomer Shlomi & Joshua D Rabinowitz
Nature Chemical Biology volume 12, pages482–489 (2016)
https://www.nature.com/articles/nchembio.2077#Sec21
"""

import sys
from datetime import datetime
import numpy as np
import math
import csv
import copy
import os

from Simulator import FSimulator
from Plotter import FPlotter

AvogadroNum = 6.022e23
GlobalKineticScale = 1

def Conc2Str(Conc):
    AbsConc = abs(Conc)
    if AbsConc >= 1e-1:
        Str = "{:.3f}  M".format(Conc)
    elif AbsConc >= 1e-4:
        Str = "{:.3f} mM".format(Conc * 1e3)
    elif AbsConc >= 1e-7:
        Str = "{:.3f} uM".format(Conc * 1e6)
    elif AbsConc >= 1e-10:
        Str = "{:.3f} nM".format(Conc * 1e9)
    else:
        Str = "{:.3f} pM".format(Conc * 1e12)
    return Str


class EcoliInfo:
    DuplicationTime_LogPhase = 20 * 60

    GenomeSize = 4.5e6
    ProteomeSize = 3e6 * 300

    # KEY PARAMETERS #
    MiniEcoli = True
    if MiniEcoli:
        New_DuplicationTime_LogPhase = 10
        ScaleFactor = New_DuplicationTime_LogPhase / DuplicationTime_LogPhase
        GenomeSize *= ScaleFactor
        ProteomeSize *= ScaleFactor
        DuplicationTime_LogPhase = New_DuplicationTime_LogPhase

    DNAReplicationRate = GenomeSize / DuplicationTime_LogPhase
    ATPConsumptionPerdNTPExtension = 3

    ProteinSynthesisRate = ProteomeSize / DuplicationTime_LogPhase
    ATPConsumptionPerAAExtension = 10

    # EC stands for Energy Consumption in ATP molecules (count)
    ECC_DNAReplication = GenomeSize * 2 * 3   # 4.5Mbp genome (double strand), 3 ATPs per 1 nucleotide extension
    ECC_ProteinSynthesis = ProteomeSize * 10  # 3M proteins (each 300aa), 10 ATPs per 1 amino acid extension
    ECC_Cytokinesis = 0
    ECC_CellDivision = ECC_DNAReplication + ECC_ProteinSynthesis + ECC_Cytokinesis

    # Convert count to M assuming E. coli volume is 1um^3
    Volume = 1e-15
    C2M = 1 / AvogadroNum / Volume
    M2C = 1 / C2M

    # ECM stands for Energy Consumption in ATP (M = mol/L)
    ECM_CellDivision = ECC_CellDivision * C2M
    ECM_CellDivision_Sec = ECM_CellDivision / DuplicationTime_LogPhase

    # ECGM stands for Energy Consumption in Glucose (M)
    ECGM_CellDivision = ECM_CellDivision / 32

    # EGM stands for Energy Generation in ATP (M)
    EGM_Glycolysis_Sec = ECM_CellDivision_Sec * 0.9
    EGM_PyruvateOxidation_Sec = EGM_Glycolysis_Sec
    EGM_TCACycle_Sec = ECM_CellDivision_Sec * 0.5
    EGM_OxidativePhosphorylation_NADH_Sec = ECM_CellDivision_Sec * 2.5   # Needs to run fast to maintain NAD+ to support Glycolysis
    EGM_OxidativePhosphorylation_FADH2_Sec = ECM_CellDivision_Sec * 0.5

    def Info():
        print("Molecule Consumption - Glucose:                                {:>10}".format(
            Conc2Str(EcoliInfo.ECGM_CellDivision)))
        print("Energy Consumption   - Cell Division per sec:                  {:>10}".format(
            Conc2Str(EcoliInfo.ECM_CellDivision_Sec)))
        print("Energy Generation    - Glycolysis per sec:                     {:>10}".format(
            Conc2Str(EcoliInfo.EGM_Glycolysis_Sec)))
        print("Energy Generation    - TCA & Oxidative Phosporylation per sec: {:>10}".format(
            Conc2Str(EcoliInfo.EGM_TCACycle_Sec + EcoliInfo.EGM_OxidativePhosphorylation_NADH_Sec)))
        print("")

    def GetDNAReplicationTime():
        return EcoliInfo.GenomeSize / float(EcoliInfo.DNAReplicationRate)

    def GetProteinSynthesisTime():
        return EcoliInfo.ProteomeSize / float(EcoliInfo.ProteinSynthesisRate)

    def OpenKnownMolConc():
        def LoadTSVDatabase(db_fname):
            db = None
            with open(db_fname) as fp:
                csv_reader = csv.reader(fp, delimiter='\t')
                list_of_rows = list(csv_reader)
                db = list_of_rows[1:]
            return db

        def ParseMetaboliteCon(db_KnownMolConc):

            def MetaboliteSynonyms(Name):
                if Name == "glucose-6-phosphate":
                    return "G6P"
                elif Name == "glyceraldehyde-3-phosphate":
                    return "GAP"
                elif Name == "coenzyme-A":
                    return "CoA-SH"
                # elif " (assumed 1 / 2 ile+leu)" in Name:
                #     return Name.replace(" (assumed 1 / 2 ile+leu)", "")
                elif Name == "ribose-5-phosphate":
                    return "R5P"
                elif Name == "ribulose-5-phosphate":
                    return "Ru5P"
                elif Name == "erythrose-4-phosphate":
                    return "E4P"
                else:
                    return Name

            KnownMolConc = dict()

            for i, Value in enumerate(db_KnownMolConc):
                Metabolite, ConcInEcoli, LowerBound, UpperBound = Value
                if ConcInEcoli == '-':
                    continue
                Metabolite = Metabolite.replace('[c]', '').replace('[m]', '')
                Metabolite = MetaboliteSynonyms(Metabolite)
                KnownMolConc[Metabolite] = [float(ConcInEcoli), float(LowerBound), float(UpperBound)]

            return KnownMolConc

        db_KnownMolConc = LoadTSVDatabase(os.path.join(os.path.dirname(__file__), "MetaboliteConcentrations.tsv"))
        KnownMolConc = ParseMetaboliteCon(db_KnownMolConc)

        # FADH2 missing data
        KnownMolConc["FADH2"] = [(KnownMolConc["NADH"][0] / KnownMolConc["NAD+"][0]) * KnownMolConc["FAD"][0],
                                 (KnownMolConc["NADH"][0] / KnownMolConc["NAD+"][0]) * KnownMolConc["FAD"][1],
                                 (KnownMolConc["NADH"][0] / KnownMolConc["NAD+"][0]) * KnownMolConc["FAD"][2]]

        return KnownMolConc


class Reaction:
    def __init__(self):
        self.ReactionName = ""
        self.Input = dict()
        self.Output = dict()
        self.Stoich = dict()
        self.Capacity = 0
        self.CapacityConstant = 0
        self.Regulators = dict()     # Capacity modulator, enzyme name: critical gene expression threshold for full function

    def Specification(self, Molecules, InitCond):
        return {}

    def Callback(self, dMolecules):
        None

    # Progress is between 0.0 and 1.0
    def GetProgress(self):
        return 0.0

    def SetProgress(self, Progress):
        None

    def UpdateCapacity(self, Molecules, RateLimitingCofactors):
        CapacityFactors = 1 # set to 1M by default
        for i in range(len(RateLimitingCofactors)):
             CapacityFactors = min(CapacityFactors, Molecules[RateLimitingCofactors[i]])
        self.Capacity = CapacityFactors * self.CapacityConstant

    def ApplyRegulation(self, Molecules):
        RegulatedCapacity = 1
        for regulator, threshold in self.Regulators.items():
            if Molecules[regulator] < threshold:
                RegulatedCapacity *= (Molecules[regulator] / threshold)
        return RegulatedCapacity * self.Capacity

    def GetChemicalEquationStr(self):

        def AddMolCoeff(Eq, MolCoeffDict):
            for n, (mol, coeff) in enumerate(MolCoeffDict.items()):
                if n != 0:
                    Eq += " + "
                if coeff > 1:
                    Eq += str(coeff) + " "
                Eq += mol

            return Eq

        Eq = ""
        Eq = AddMolCoeff(Eq, self.Input)
        Eq += " -> "
        Eq = AddMolCoeff(Eq, self.Output)
        return Eq


class NADPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"nicotinamide": 1, "PRPP": 1, "ATP": 2}
        self.Output = {"NAD+": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["NAD+"] / (InitCond["NAD+"] + Molecules["NAD+"]) * self.CapacityConstant
        return "NAD+", VO


class NADPPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"NAD+": 1, "ATP": 1}
        self.Output = {"NADP+": 1, "ADP": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["NADP+"] / (InitCond["NADP+"] + Molecules["NADP+"]) * self.CapacityConstant
        return "NADP+", VO


class CoASynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"pantothenate": 1, "cysteine": 1, "ATP": 4}
        # self.Input = {"pantothenate": 1, "cysteine": 1, "ATP": 3, "CTP": 1, "CO2": 1}
        self.Output = {"CoA-SH": 1, "ADP": 2, "AMP": 1}
        # self.Output = {"CoA-SH": 1, "ADP": 2, "CMP": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["CoA-SH"] / (InitCond["CoA-SH"] + Molecules["CoA-SH"]) * self.CapacityConstant
        return "CoA-SH", VO


class ACPSynthesis(Reaction):
    '''
    Enzymes     Reactions
    Acs:        acetate + ATP + CoA = acetyl-CoA + AMP + diphosphate (not the focus in this pathway)
    AcpS, AcpT: apo-[ACP] + CoA = 3',5'-ADP + holo-[ACP] (the focus in this pathway)
    AcpH:       H2O + holo-[ACP] = 4'-phosphopantetheine + apo-[ACP] + H+ (not the focus in this pathway)
    '''
    def __init__(self):
        super().__init__()
        self.Input = {"apo-[ACP]": 1, "CoA-SH": 1}
        self.Output = {"3',5'-ADP;": 1, "holo-[ACP]": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["holo-[ACP]"] / (InitCond["holo-[ACP]"] + Molecules["holo-[ACP]"]) * self.CapacityConstant
        return "holo-[ACP]", VO


class IPPSynthesis(Reaction):
    # reaction IPPsynthesis(pyruvate + GAP + NADPH + CTP + ATP -> IPP + NADP+ + ADP + CMP)
    def __init__(self):
        super().__init__()
        self.Input = {"pyruvate": 1, "GAP": 1, "NADPH": 1, "ATP": 1, "CTP": 1}
        self.Output = {"IPP": 1, "NADP+": 1, "ADP": 1, "CMP": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["pyruvate"] / (InitCond["pyruvate"] + Molecules["pyruvate"]) * self.CapacityConstant
        return "IPP", VO


class UbiquinoneSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"chorismate": 1, "IPP": 8}
        self.Output = {"UQ": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["chorismate"] / (InitCond["chorismate"] + Molecules["chorismate"]) * self.CapacityConstant
        return "UQ", VO


# Vitamins
class ThiamineSynthesis(Reaction):
    # PRPP + 2 glutamine + glycine + 2 SAM + 5 ATP + tyrosine + cysteine + pyruvate + GAP
    # -> TDP + AMP + 4 ADP + formate + 4-cresol + 2 methionine + 2 5’-deoxyadenosine + alanine + 2 glutamate
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Thiamine Synthesis'
        # self.Input = {"PRPP": 1, "SAM": 2, "ATP": 5, "glutamine": 2, "glycine": 1, "tyrosine": 1, "cysteine": 1, "pyruvate": 1, "GAP": 1}
        # self.Input = {"PRPP": 1, "ATP": 5, "pyruvate": 1, "GAP": 1}
        self.Input = {"PRPP": 1, "pyruvate": 2}   # 1 NADH and 2 ATP lost by substituting GAP with pyruvate.
        # self.Output = {"TDP": 1, "AMP": 1, "ADP": 4, "formate": 1, "4-cresol": 1, "methionine": 2, "5’-deoxyadenosine": 2, "alanine": 1, "glutamate": 2}
        # self.Output = {"TDP": 1, "AMP": 1, "ADP": 4, "formate": 1}
        self.Output = {"TDP": 1}   # adjusted ADP according to the ATP loss by GAP to pyruvate substitution
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["PRPP"] / (InitCond["PRPP"] + Molecules["PRPP"]) * self.CapacityConstant
        return "TDP", VO


class RiboflavinSynthesis(Reaction):
    # Ru5P + GTP + NADP+ -> riboflavin + 2 formate + NH3 + NADPH
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Riboflavin Synthesis'
        self.Input = {"Ru5P": 1, "GTP": 1, "NADP+": 1}
        # self.Output = {"riboflavin": 1, "formate": 2, "NH3": 1, "NADPH": 1}
        self.Output = {"riboflavin": 1, "NADPH": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["Ru5P"] / (InitCond["Ru5P"] + Molecules["Ru5P"]) * self.CapacityConstant
        return "riboflavin", VO


class PantotheateSynthesis(Reaction):
    # Ru5P + GTP + NADP+ -> riboflavin + 2 formate + NH3 + NADPH
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Pantotheate Synthesis'
        self.Input = {"pyruvate": 2, "NADPH": 2, "aspartate": 1, "ATP": 1}
        # self.Output = {"riboflavin": 1, "formate": 2, "NH3": 1, "NADPH": 1}
        self.Output = {"pantotheate": 1, "NADP+": 2, "AMP": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["pyruvate"] / (InitCond["pyruvate"] + Molecules["pyruvate"]) * self.CapacityConstant
        return "pantotheate", VO


class PyridoxalPhosphateSynthesis(Reaction):
    # E4P + 3 NAD+ + glutamate + pyruvate + GAP -> PLP + 3 NADH + 2-oxoglutarate)
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Pyridoxal Phosphate Synthesis'
        self.Input = {"E4P": 1, "NAD+": 3, "glutamate": 1, "pyruvate": 1, "GAP": 1}
        self.Output = {"PLP": 1, "NADH": 3, "a-ketoglutarate": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["E4P"] / (InitCond["E4P"] + Molecules["E4P"]) * self.CapacityConstant
        return "PLP", VO


class BiotinSynthesis(Reaction):
    # pimeloyl-[ACP] + alanine + AMTOB + ATP -> biotin + [ACP] + AMP + 5’-deoxyadenosine
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Biotin Synthesis'
        self.Input = {"pimeloyl-[ACP]": 1, "alanine": 1, "AMTOB": 1, "ATP": 1}
        self.Output = {"biotin": 1, "apo-[ACP]": 1, "AMP": 1, "5’-deoxyadenosine": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["pimeloyl-[ACP]"] / (InitCond["pimeloyl-[ACP]"] + Molecules["pimeloyl-[ACP]"]) * self.CapacityConstant
        return "biotin", VO


class FolateSynthesis(Reaction):
    # pABA + GTP + 2 ATP + glutamate + 2 NADPH -> tetrahydrofolate + 2 NADP+ + AMP
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Folate Synthesis'
        self.Input = {"pABA": 1, "GTP": 1, "ATP": 2, "glutamate": 1, "NADPH": 2}
        self.Output = {"tetrahydrofolate": 1, "NADP+": 2, "AMP": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["pABA"] / (InitCond["pABA"] + Molecules["pABA"]) * self.CapacityConstant
        return "tetrahydrofolate", VO


# Central Carbon Metabolism
class Glycolysis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Glycolysis'
        self.Input = {"G6P": 1, "ADP": 2, "NAD+": 2}
        self.Output = {"pyruvate": 2, "NADH": 2, "ATP": 2}
        self.CapacityConstant = EcoliInfo.EGM_Glycolysis_Sec

    def Specification(self, Molecules, InitCond):
        MinNADPlusConc = InitCond["NAD+"] / 10.0
        if Molecules["NAD+"] <= MinNADPlusConc:
            return "ATP", 0.0

        MinATPConc = 0.1e-3
        VO = Molecules["ADP"] / max(Molecules["ATP"], MinATPConc) * self.CapacityConstant
        return "ATP", VO


class OxidativePentosePhosphatePathway(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'OxidativePPP'
        self.Input = {"G6P": 1, "NADP+": 2}
        self.Output = {"Ru5P": 1, "NADPH": 2}
        self.CapacityConstant = EcoliInfo.EGM_PyruvateOxidation_Sec

    def Specification(self, Molecules, InitCond):
        VO = Molecules["G6P"] / (InitCond["G6P"] + Molecules["G6P"]) * self.CapacityConstant
        return "Ru5P", VO


class PyruvateOxidation(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Pyruvate Oxidation'
        self.Input = {"pyruvate": 1, "NAD+": 1, "CoA-SH": 1}
        self.Output = {"acetyl-CoA": 1, "NADH": 1}
        self.CapacityConstant = EcoliInfo.EGM_PyruvateOxidation_Sec

    def Specification(self, Molecules, InitCond):
        VO = Molecules["pyruvate"] / (InitCond["pyruvate"] + Molecules["pyruvate"]) * self.CapacityConstant
        return "acetyl-CoA", VO


class TCACycle(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'TCA cycle'
        self.Input = {"acetyl-CoA": 1, "NAD+": 3, "FAD": 1, "ADP": 1}
        self.Output = {"NADH": 3, "FADH2": 1, "ATP": 1, "CoA-SH": 1}
        self.CapacityConstant = EcoliInfo.EGM_TCACycle_Sec

    def Specification(self, Molecules, InitCond):
        VO = Molecules["acetyl-CoA"] / (InitCond["acetyl-CoA"] + Molecules["acetyl-CoA"]) * self.CapacityConstant
        return "ATP", VO


class NADH_OxidativePhosphorylation(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'NADH OxidativePhosphorylation'
        self.Input = {"NADH": 1, "ADP": 2.5}
        self.Output = {"NAD+": 1, "ATP": 2.5}
        self.CapacityConstant = EcoliInfo.EGM_OxidativePhosphorylation_NADH_Sec

    def Specification(self, Molecules, InitCond):
        VO = Molecules["NADH"] / (InitCond["NADH"] + Molecules["NADH"]) * self.CapacityConstant
        return "ATP", VO


class FADH2_OxidativePhosphorylation(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'FADH2 OxidativePhosphorylation'
        self.Input = {"FADH2": 1, "ADP": 1.5}
        self.Output = {"FAD": 1, "ATP": 1.5}
        self.CapacityConstant = EcoliInfo.EGM_OxidativePhosphorylation_FADH2_Sec

    def Specification(self, Molecules, InitCond):
        VO = Molecules["FADH2"] / (InitCond["FADH2"] + Molecules["FADH2"]) * self.CapacityConstant
        return "ATP", VO


class dNTPSynthesis(Reaction):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRate):
        super().__init__()
        self.ReactionName = 'dNTP Synthesis'
        self.Input = {}
        self.Output = {"dATP": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "dATP", self.Rate * EcoliInfo.C2M


class AASynthesis(Reaction):
    def __init__(self, Rate = EcoliInfo.ProteinSynthesisRate):
        super().__init__()
        self.ReactionName = 'AA Synthesis'
        self.Input = {}
        self.Output = {"glutamine": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "glutamine", self.Rate * EcoliInfo.C2M


class ATPControl(Reaction):
    def __init__(self, ControlRate=-4.35E-03):   # Cell Division ATP consumption: c = 4.35E-03
        super().__init__()
        self.ReactionName = 'ATP Control'
        self.Input = {"ATP": 1}
        self.Output = {"ADP": 1}
        self.Capacity = ControlRate

    def Specification(self, Molecules, InitCond):
        return "ATP", self.Capacity


class NADHControl(Reaction):
    def __init__(self, ControlRate=-1E-05):
        super().__init__()
        self.ReactionName = 'NADH Control'
        self.Input = {"NADH": 1}
        self.Output = {"NAD+": 1}
        self.Capacity = ControlRate

    def Specification(self, Molecules, InitCond):
        return "NADH", self.Capacity


class FADH2Control(Reaction):
    def __init__(self, ControlRate=-1E-05):
        super().__init__()
        self.ReactionName = 'FADH2 Control'
        self.Input = {"FADH2": 1}
        self.Output = {"FAD": 1}
        self.Capacity = ControlRate

    def Specification(self, Molecules, InitCond):
        return "FADH2", self.Capacity


# For Debugging
class MolControl(Reaction):
    def __init__(self, MolName, ControlRate=-1e-5):
        super().__init__()
        self.ReactionName = MolName + ' Consumption'
        self.Input = {MolName: 1}
        # self.Output = {"ATP": 1}

        self.ConsumedMol = MolName
        self.Capacity = ControlRate

    def Specification(self, Molecules, InitCond):
        return self.ConsumedMol, self.Capacity


class Process(Reaction):
    def __init__(self):
        super().__init__()
        self.BuildingBlocks = dict() # gets incorporated into input
        self.EnergyConsumption = dict()   # ATP Consumption, gets incorporated into stoichiometry
        self.Rate = 0.0
        self.Progress = 0.0
        self.MaxProgress = 0.0

    def GetProgress(self):
        return self.Progress / self.MaxProgress

    def SetProgress(self, Progress):
        self.Progress = self.MaxProgress * Progress


class DNAReplication(Process):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRate):
        super().__init__()
        self.ReactionName = "DNA Replication"
        # self.BuildingBlocks = {"dATP": 1, "dCTP": 1, "dGTP": 1, "dTTP": 1}
        self.BuildingBlocks = {"dATP": 1}
        self.EnergyConsumption = EcoliInfo.ATPConsumptionPerdNTPExtension
        self.Regulators = {}
        self.Rate = Rate
        self.Progress = 0.0
        self.MaxProgress = EcoliInfo.GenomeSize

        self.Input = {"ATP": self.EnergyConsumption, "dATP": 1}
        self.Output = {"ADP": self.EnergyConsumption}

    def Specification(self, Molecules, InitCond):
        return "dATP", -self.Rate * EcoliInfo.C2M

    def Callback(self, dMolecules):
        assert "dATP" in dMolecules
        dElongation = -dMolecules["dATP"] * EcoliInfo.M2C
        assert dElongation >= 0
        self.Progress += dElongation


class ProteinSynthesis(Process):
    def __init__(self, Rate = EcoliInfo.ProteinSynthesisRate):
        super().__init__()
        self.ReactionName = "Protein Synthesis"
        self.BuildingBlocks = {"glutamine": 1}
        self.EnergyConsumption = EcoliInfo.ATPConsumptionPerAAExtension
        self.Regulators = {}
        self.Rate = Rate
        self.Progress = 0.0
        self.MaxProgress = EcoliInfo.ProteomeSize

        self.Input = {"ATP": self.EnergyConsumption, "glutamine": 1}
        self.Output = {"ADP": self.EnergyConsumption}

    def Specification(self, Molecules, InitCond):
        return "glutamine", -self.Rate * EcoliInfo.C2M

    def Callback(self, dMolecules):
        assert "glutamine" in dMolecules
        dElongation = -dMolecules["glutamine"] * EcoliInfo.M2C
        assert dElongation >= 0
        self.Progress += dElongation


class ReactionSimulator(FSimulator):
    def __init__(self):
        super().__init__()
        self.Reactions = []
        self.Dataset = {}

        self.KnownMolConc = EcoliInfo.OpenKnownMolConc()

        self.Molecules = dict()
        self.InitialConditions = dict()
        self.PermanentMolecules = list()

        self.dMolecules = dict()

        self.Plot = False
        self.ExportToPDF = False
        self.Debug_Info = 1
        self.Debug_Reaction = False

    def AddPermanentMolecules(self, PermanentMolecules):
        for Molecule in PermanentMolecules:
            self.PermanentMolecules.append(Molecule)

    def SetPermanentMolecules(self, PermanentMolecules):
        self.PermanentMolecules = PermanentMolecules

    def Initialize(self, Molecules={}):
        self.InitializeStoich()
        self.InitializeMolecules(Molecules)
        self.InitializeDataset()
        self.ApplyGlobalKineticCapacity()

    def InitializeStoich(self):
        for Reaction in self.Reactions:
            for Mol, Coeff in Reaction.Input.items():
                Reaction.Stoich[Mol] = -Coeff
            for Mol, Coeff in Reaction.Output.items():
                Reaction.Stoich[Mol] = Coeff

    def GetKnownMolConc(self, Molecule):
        if Molecule in self.KnownMolConc:
            return self.KnownMolConc[Molecule][0]
        else:
            print(Molecule, " is not found in the KnownMolConc database")
            return 0.0

    def GetReactionNames(self, Eq=False, Alignment=False):
        ReactionNames = ""
        for Reaction in self.Reactions:
            ReactionLabel = ""
            ReactionName = "[" + Reaction.ReactionName + "]"
            if Eq:
                ReactionEq = Reaction.GetChemicalEquationStr()
                if Alignment:
                    ReactionLabel += "{:<35} {}".format(ReactionName, ReactionEq)
                else:
                    ReactionLabel += ReactionName + ReactionEq
            else:
                ReactionLabel += ReactionName
            ReactionNames += ReactionLabel + "\n"
        return str(ReactionNames)

    def PrintReactions(self):
        print('-- Reactions --')
        print(self.GetReactionNames(Eq=True, Alignment=True))

    def InitializeMolecules(self, Molecules={}):
        AllMolecules = dict()
        for Reaction in self.Reactions:
            for Molecule, Coeff in Reaction.Stoich.items():
                if Molecule not in AllMolecules:
                    AllMolecules[Molecule] = self.GetKnownMolConc(Molecule)

        self.InitialConditions = dict(sorted(AllMolecules.items()))
        self.Molecules = self.InitialConditions.copy()
        for Molecule, Conc in Molecules.items():
            self.Molecules[Molecule] = Conc

        self.Molecules_Min = {}
        for Mol in ["ADP", "pyruvate", "acetyl-CoA", "NADH", "FADH2"]:
            if Mol in self.InitialConditions:
                self.Molecules_Min[Mol] = self.InitialConditions[Mol]

    def AddReaction(self, Reaction):
        self.Reactions.append(Reaction)

    def GetMolCount(self, Conc):
        return Conc * EcoliInfo.Volume * AvogadroNum

    def GetConc(self, Count):
        return Count / EcoliInfo.Volume / AvogadroNum

    def AdjustRefdCon(self, Reaction, RefMol, RefdConc):
        # Compare dConc of reference molecule to input concentrations and adjust reference dConc
        if len(Reaction.Input) == 0:
            return RefdConc        

        RefCoeff = Reaction.Stoich[RefMol]
        UnitdConc = RefdConc / RefCoeff

        Out = 10  # 10 Molar to begin comparison
        for Mol, Coeff in Reaction.Input.items():
            MinConc = self.Molecules_Min[Mol] if Mol in self.Molecules_Min else 0.0
            assert self.Molecules[Mol] >= MinConc * 0.999999
            AdjustedConc = (self.Molecules[Mol] - MinConc) / Coeff
            Out = min(Out, min(AdjustedConc, UnitdConc))

            # The following is possible due to multiplication/division in floating numbers
            if Out * Coeff > self.Molecules[Mol] - MinConc:
                Out *= 0.999999
            assert Out * Coeff <= self.Molecules[Mol] - MinConc

        return Out * RefCoeff

    def DeterminedConc(self, Reaction, RefMol, RefdConc):
        UnitdConc = RefdConc / Reaction.Stoich[RefMol]
        dConc = dict()
        for Mol, Coeff in Reaction.Stoich.items():
            dConc[Mol] = UnitdConc * Coeff

        return dConc

    def RunReaction(self, Reaction, DeltaTime):
        RefMol, RefdConc = Reaction.Specification(self.Molecules, self.InitialConditions)
        RefdConc *= DeltaTime
        RefdConc = self.AdjustRefdCon(Reaction, RefMol, RefdConc)
        return self.DeterminedConc(Reaction, RefMol, RefdConc)

    def UpdateMolecules(self, dMolecules):
        for dMolecule, dConc in dMolecules.items():
            assert dMolecule in self.Molecules
            assert dConc + self.Molecules[dMolecule] >= 0
            self.Molecules[dMolecule] += dConc
            self.Molecules[dMolecule] = max(0, self.Molecules[dMolecule])

    def RestorePermanentMolecules(self):
        for Molecule in self.PermanentMolecules:
            if Molecule in self.Molecules:
                self.Molecules[Molecule] = self.InitialConditions[Molecule]

    def CheckZeroConcentrations(self):
        for Molecule, Conc in self.Molecules.items():
            if Conc < self.GetConc(1):
                self.Molecules[Molecule] = 0   # Make sure that division by zero would not happen in reaction specifications
                # self.Molecules[Molecule] = self.GetConc(0.1)
                # assert self.Molecules[Molecule] > 0

    def InitializeDataset(self):
        for Molecule, Conc in self.Molecules.items():
            self.Dataset[Molecule] = [Conc]

    def AddToDataset(self):
        for Molecule, Conc in self.Molecules.items():
            self.Dataset[Molecule].append(Conc)

    def ApplyGlobalKineticCapacity(self):
        for Reaction in self.Reactions:
            Reaction.Capacity *= GlobalKineticScale
            Reaction.CapacityConstant *= GlobalKineticScale

    def Info(self):
        HeadStr = "{:<16}: {:>10} {:>10}  ".format("", "Initial", "Current")
        for ReactionName, dMolecules in self.dMolecules.items():
            HeadStr += " {:<10.10}".format(ReactionName)
        print(HeadStr)
        for Molecule, Conc in self.Molecules.items():
            InitConc = self.InitialConditions[Molecule]
            InitConcStr = Conc2Str(InitConc)
            ConcStr = Conc2Str(Conc)
            LineStr = "{:<16}: {:>10} {:>10}  ".format(Molecule, InitConcStr, ConcStr)
            for ReactionName, dMolecules in self.dMolecules.items():
                if Molecule in dMolecules:
                    Conc2 = dMolecules[Molecule]
                    ConcStr2 = Conc2Str(Conc2)
                else:
                    ConcStr2 = ""
                LineStr += " {:>10}".format(ConcStr2)
            print(LineStr)

    def Summary(self):
        print("{:<16} {:>12} {:>12} {:>13}  {:<12}".format("", "Initial", "Final", "dConc", "OutOfRange"))
        for Molecule, Conc in self.Molecules.items():
            dConc = Conc - self.InitialConditions[Molecule]

            InitConc = Conc2Str(self.InitialConditions[Molecule])
            FinalConc = Conc2Str(Conc)
            FinaldConc = Conc2Str(dConc) if Molecule not in self.PermanentMolecules else "-"
            if dConc < 0:
                FinaldConc = "\033[91m {:>12}\033[00m".format(FinaldConc)
            elif dConc > 0:
                FinaldConc = "\033[92m {:>12}\033[00m".format(FinaldConc)

            OutOfRange = ''
            if Conc < self.KnownMolConc[Molecule][1]:
                OutOfRange = "\033[91m {:<12}\033[00m".format('DEPLETED')
            elif Conc > self.KnownMolConc[Molecule][2]:
                OutOfRange = "\033[92m {:<12}\033[00m".format('ACCUMULATED')

            print("{:16} {:>12} {:>12} {:>12} {:<12}".format(Molecule, InitConc, FinalConc, FinaldConc, OutOfRange))

        if self.Plot:
            print("\nPlot: ON")
        else:
            print("\nPlot: OFF")

    def PrintMolConc(self, MolConcDict, End='\t'):
        for Molecule, Conc in MolConcDict.items():
            ConStr = Conc2Str(Conc)
            print("\t{}: {}".format(Molecule, ConStr), end=End)
        print()

    def SimulateDelta(self, DeltaTime):
        self.dMolecules = {}
        for Reaction in self.Reactions:
            dMolecules = self.RunReaction(Reaction, DeltaTime)
            if self.Debug_Reaction:
                print("[{}] {}".format(Reaction.ReactionName, Reaction.GetChemicalEquationStr()))
                self.PrintMolConc(dMolecules, End=', ')
            self.dMolecules[Reaction.ReactionName] = dMolecules

        while True:
            dMolecules_Sum = {}
            for ReactionName, dMolecules in self.dMolecules.items():
                for Mol, dConc in dMolecules.items():
                    if Mol not in dMolecules_Sum:
                        dMolecules_Sum[Mol] = [0, 0, []]   # Plus_dConc, Minus_dConc, ReactionNameList
                    if dConc >= 0:
                        dMolecules_Sum[Mol][0] += dConc
                    else:
                        dMolecules_Sum[Mol][1] += dConc
                        dMolecules_Sum[Mol][2].append(ReactionName)


            Updated = False
            for Mol, Value in dMolecules_Sum.items():
                Conc = self.Molecules[Mol]
                MinConc = self.Molecules_Min[Mol] if Mol in self.Molecules_Min else 0.0
                Plus_dConc, Minus_dConc, ReactionNameList = Value
                dConc = Plus_dConc + Minus_dConc

                if Conc + dConc >= MinConc * 0.999999:
                    continue

                Updated = True
                UpdateFactor = (Minus_dConc + (MinConc - Conc - dConc)) / Minus_dConc * 0.999999
                assert UpdateFactor < 1.0

                """
                print(Mol)
                print("Conc:", Conc2Str(Conc))
                print("MinConc:", Conc2Str(MinConc))
                print("dConc:", Conc2Str(dConc))
                print("Plus_dConc:", Conc2Str(Plus_dConc))
                print("Minus_dConc:", Conc2Str(Minus_dConc))
                print("UpdateFactor:", UpdateFactor)
                """

                for ReactionName in ReactionNameList:
                    dMolecules = self.dMolecules[ReactionName]
                    for Mol, Conc in dMolecules.items():
                        dMolecules[Mol] *= UpdateFactor

            if not Updated:
                break

        Sum_dMolecules = {}
        for Reaction in self.Reactions:
            dMolecules = self.dMolecules[Reaction.ReactionName]
            Reaction.Callback(dMolecules)
            for Mol, dConc in dMolecules.items():
                if Mol in Sum_dMolecules:
                    Sum_dMolecules[Mol] += dConc
                else:
                    Sum_dMolecules[Mol] = dConc

        self.UpdateMolecules(Sum_dMolecules)
        self.CheckZeroConcentrations()

        if len(self.PermanentMolecules):
            self.RestorePermanentMolecules()

        self.AddToDataset()        

    def GetDataset(self):
        return {self.GetReactionNames(Eq=False): self.Dataset}


    def UnitTest(self, ReactionName):
        # a-ketoGlutarate synthesis unit test
        if ReactionName == "Glycolysis":
            Sim.AddReaction(Glycolysis())
            Sim.AddReaction(ATPControl(-1e-5))

        elif ReactionName == "PyruvateOxidation":
            Sim.AddReaction(PyruvateOxidation())
            Sim.AddReaction(MolControl("acetyl-CoA", -3e-6))

        else:
            print("\nNo unit test has been selected\n")
            sys.exit(1)


def GetUnitTestReactions():
    '''
    Unit tests available:
    '''
    return ["Glycolysis",
            "PyruvateOxidation",
            "TCACycle_a-KetoGlutarateSynthesis",
            "TCACycle_SuccinylCoASynthesis",
            "TCACycle_FumarateSynthesis",
            "TCACycle_OxaloacetateSynthesis",
            "PyruvateCarboxylation",
            "OxidativePhosphorylation",
            ]


if __name__ == '__main__':
    Sim = ReactionSimulator()

    RunUnitTest = False
    # RunUnitTest = True

    EcoliInfo.Info()
    ATPConsumption_Sec = EcoliInfo.ECM_CellDivision_Sec
    # ATPConsumption_Sec = 0

    if RunUnitTest:
        UnitTestReactions = "Glycolysis"
        Sim.UnitTest(UnitTestReactions)

    else:
        # Add Reactions
        # Sim.AddReaction(NADPlusSynthesis())
        # Sim.AddReaction(NADPPlusSynthesis())
        # Sim.AddReaction(CoASynthesis())
        Sim.AddReaction(Glycolysis())
        Sim.AddReaction(PyruvateOxidation())
        Sim.AddReaction(TCACycle())
        # Sim.AddReaction(TCACycle_alphaKetoGlutarateSynthesis())
        # Sim.AddReaction(TCACycle_SuccinylCoASynthesis())
        # Sim.AddReaction(TCACycle_FumarateSynthesis())
        # Sim.AddReaction(TCACycle_OxaloacetateSynthesis())
        # Sim.AddReaction(PyruvateCarboxylation())
        Sim.AddReaction(NADH_OxidativePhosphorylation())
        Sim.AddReaction(FADH2_OxidativePhosphorylation())
        Sim.AddReaction(ATPControl(-ATPConsumption_Sec))

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
        Sim.SetPermanentMolecules(PermanentMolecules)

    # Debugging options
    # Sim.Debug_Reaction = True
    Sim.Debug_Info = 100
    Sim.Plot = True
    # Sim.ExportToPDF = True

    # Set initial molecule concentrations
    Molecules = {}
    # Molecules["ATP"] = 1.0 * 1e-3
    # Molecules["ADP"] = 8.6 * 1e-3
    # Molecules["G6P"] = 50 * 1e-3

    # Execute simulation
    Sim.Initialize(Molecules)

    # Simulation parameters
    TotalTime = Sim.Molecules["G6P"] * 32 / max(1e-3, ATPConsumption_Sec) + 200
    DeltaTime = 0.01

    Sim.Simulate(TotalTime=TotalTime, DeltaTime=DeltaTime)

    # Plot
    Datasets = Sim.GetDataset()
    Plot = FPlotter()

    if Sim.Plot:
        Plot.SetKnownMolConc(EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Unitless=False, Multiscale=True)
        # Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, All=True, Multiscale=True, Include_nM=True)
        Plot.SetKnownMolConc(EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='right', Unitless=False, All=False, Individual=True, MolRange=True)

    if Sim.ExportToPDF:
        Plot.SetKnownMolConc(EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Export='pdf')
        Plot.SetKnownMolConc(EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='right', All=False, Individual=True, MolRange=True, Export='pdf')

