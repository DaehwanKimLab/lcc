# BSD 3-Clause License
# © 2023, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

"""
Reference

Metabolite concentrations, fluxes and free energies imply efficient enzyme usage
Junyoung O Park, Sara A Rubin, Yi-Fan Xu, Daniel Amador-Noguez, Jing Fan, Tomer Shlomi & Joshua D Rabinowitz
Nature Chemical Biology volume 12, pages482–489 (2016)
https://www.nature.com/articles/nchembio.2077#Sec21

Measuring and modeling energy and power consumption in living microbial cells with a synthetic ATP reporter
Yijie Deng, Douglas Raymond Beahm, Steven Ionov & Rahul Sarpeshkar
BMC Biology volume 19, Article number: 101 (2021)
https://doi.org/10.1186/s12915-021-01023-2
->  Bacteria turn over their cellular ATP pool a few times per second during the exponential phase and
    slow this rate by ~ 2–5-fold in lag and stationary phases.

The rate of turnover of the adenosine triphosphate pool of Escherichia coli growing aerobically in simple defined media
W. H. Holms, I. D. Hamilton & A. G. Robertson
Archiv für Mikrobiologie volume 83, pages95–109 (1972)
https://doi.org/10.1007/BF00425016
->  During logarithmic growth the rate of turnover of the ATP pool falls within the range 250–450 times min-1.

Fundamental limits on the rate of bacterial growth and their influence on proteomic composition
Nathan M. Belliveau 1 9, Griffin Chure 2 9 10, Christina L. Hueschen 3, Hernan G. Garcia 4, Jane Kondev 5, Daniel S. Fisher 6, Julie A. Theriot 1, Rob Phillips
Cell Systems Volume 12, Issue 9, Pages 924-944.e2, (2021)
https://doi.org/10.1016/j.cels.2021.06.002
->  In our estimate of ATP production above we found that a cell demands about 5 × 109.
    ATP per cell cycle or ≈ 1 × 10e6 ATP/s. With a cell volume of roughly 1 fL (BNID: 100004),
    this corresponds to about 2 × 10e10 ATP per fL of cell volume, in line with previous estimates (Stouthamer, 1973)
    and within 3–4 fold of more extensive calculations (Feist et al., 2007; Szenk et al., 2017).
    # DL: 1 × 10e6 ATP/s seems to be a typo. It should be just 1 × 10e6 ATP.
->  Maximum growth rate is determined by the rate of ribosomal synthesis

A theoretical study on the amount of ATP required for synthesis of microbial cell material
Antonie Leeuwenhoek, 39 (1973), pp. 545-565, 10.1007/BF02578899

A.M. Feist, C.S. Henry, J.L. Reed, M. Krummenacker, A.R. Joyce, P.D. Karp, L.J. Broadbelt, V. Hatzimanikatis, B.Ø. Palsson
A genome-scale metabolic reconstruction for Escherichia coli K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information
Mol. Syst. Biol., 3 (2007), p. 121, 10.1038/msb4100155

M. Szenk, K.A. Dill, A.M.R. de Graff
Why do fast-growing bacteria enter overflow metabolism? Testing the membrane real estate hypothesis
Cell Syst., 5 (2017), pp. 95-104, 10.1016/j.cels.2017.06.005

"""
import json
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
    # Convert count to M assuming E. coli volume is 1um^3
    Volume = 1e-15
    C2M = 1 / AvogadroNum / Volume
    M2C = 1 / C2M

    DuplicationTime_LogPhase = 20 * 60

    GenomeSize = 4.5e6
    ProteomeSize = 3e6 * 300

    # KEY PARAMETERS #
    MiniEcoli = False
    if MiniEcoli:
        New_DuplicationTime_LogPhase = 10
        ScaleFactor = New_DuplicationTime_LogPhase / DuplicationTime_LogPhase
        GenomeSize *= ScaleFactor
        ProteomeSize *= ScaleFactor
        DuplicationTime_LogPhase = New_DuplicationTime_LogPhase

    DNAReplicationRate = GenomeSize / DuplicationTime_LogPhase
    DNAReplicationRateInConc = DNAReplicationRate * C2M
    ATPConsumptionPerdNTPExtension = 3

    ProteinSynthesisRate = ProteomeSize / DuplicationTime_LogPhase
    ATPConsumptionPerAAExtension = 10

    VolumeExpansionRate = DNAReplicationRate
    VolumeExpansionSize = GenomeSize

    # ECC stands for Energy Consumption in ATP molecules (count)
    ECC_DNAReplication = GenomeSize * 2 * 3   # 4.5Mbp genome (double strand), 3 ATPs per 1 nucleotide extension
    ECC_ProteinSynthesis = ProteomeSize * 10  # 3M proteins (each 300aa), 10 ATPs per 1 amino acid extension
    ECC_Cytokinesis = 0
    ECC_VolumeExpansion = GenomeSize
    ECC_CellDivision = ECC_DNAReplication + ECC_ProteinSynthesis + ECC_Cytokinesis + ECC_VolumeExpansion

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

    # NCC stands for Nucleotide Consumption in molecules (count)
    NCC_CellDivision = GenomeSize * 2

    # NCM stands for Nucleotide Consumption in concentration (M = mol/L)
    NCM_CellDivision = NCC_CellDivision * C2M
    NCM_CellDivision_Sec = NCM_CellDivision / DuplicationTime_LogPhase

    # ACC stands for AA Consumption in molecules (count)
    ACC_CellDivision = ProteomeSize

    # ACM stands for AA Consumption in concentration (M = mol/L)
    ACM_CellDivision = ACC_CellDivision * C2M
    ACM_CellDivision_Sec = ACM_CellDivision / DuplicationTime_LogPhase

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

            KnownMetConc = dict()

            for i, Value in enumerate(db_KnownMolConc):
                Metabolite, ConcInEcoli, LowerBound, UpperBound = Value
                if ConcInEcoli == '-':
                    continue
                Metabolite = Metabolite.replace('[c]', '').replace('[m]', '')
                Metabolite = MetaboliteSynonyms(Metabolite)
                KnownMetConc[Metabolite] = [float(ConcInEcoli), float(LowerBound), float(UpperBound), "Park et al., (2016) Nat Chem Biol"]

            return KnownMetConc

        def ParseProteinCon(db_KnownMolConc):

            KnownProtConc = dict()

            for i, Value in enumerate(db_KnownMolConc):
                # Value has a list of gene names and molecule counts in different media.
                GeneName = Value[0]
                Count = Value[1] # in glucose medium; use [2] for LB medium
                if Count == 'NA':
                    continue
                ProteinName = GeneName[0].upper() + GeneName[1:]
                Conc = float(Count) * EcoliInfo.C2M
                KnownProtConc[ProteinName] = [Conc, Conc * 0.5, Conc * 1.5, "Schmidt et al., (2016) Nat Biotech; in glucose media"]

            return KnownProtConc

        db_KnownMolConc = LoadTSVDatabase(os.path.join(os.path.dirname(__file__), "MetaboliteConcentrations.tsv"))
        db_AdditionalKnownMolConc = LoadTSVDatabase(os.path.join(os.path.dirname(__file__), "AdditionalMetaboliteConcentrations.tsv"))
        db_KnownMolConc += db_AdditionalKnownMolConc
        KnownMetConc = ParseMetaboliteCon(db_KnownMolConc)
        db_KnownProtConc = LoadTSVDatabase(os.path.join(os.path.dirname(__file__), "ProteinConcentrations.tsv"))
        KnownProtConc = ParseProteinCon(db_KnownProtConc)

        KnownMolConc = {**KnownMetConc, **KnownProtConc}
        assert len(KnownMolConc) == len(KnownMetConc) + len(KnownProtConc), "Duplicate molecules exist in 'KnownMetConc' and 'KnownProtConc'"

        # Chemical info (no concentration)
        db_KnownMolConc_ecmdb = json.load(open(os.path.join(os.path.dirname(__file__), "ecmdb.json"), encoding='utf8'))


        return KnownMolConc

    def PrintKnownMolConc():
        KnownMolConc = EcoliInfo.OpenKnownMolConc()

        for mol, concs in KnownMolConc.items():
            Conc = concs[0]
            print(mol + " = " + Conc2Str(Conc) + ";")

class Reaction:
    MaxConc = 10
    MinConc = 0

    def __init__(self):
        self.ReactionName = ""
        self.Input = dict()
        self.Output = dict()
        self.Stoich = dict()
        self.CapacityConstant = 0
        self.Rate = 0

    def Specification(self, Molecules, InitCond):
        return {}

    def Callback(self, dMolecules):
        None

    def GetMaxConc(self, MolName, Molecules, InitCond):
        return Reaction.MaxConc

    def GetMinConc(self, MolName, Molecules, InitCond):
        return Reaction.MinConc

    # Progress is between 0.0 and 1.0
    def GetProgress(self):
        return 0.0

    def SetProgress(self, Progress):
        None

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

    def GetAllParameters(self):
        CapacityConstantOut = 0
        RateOut = 0

        if self.CapacityConstant > 0:
            CapacityConstantOut = self.CapacityConstant
        if self.Rate > 0:
            RateOut = self.Rate

        Parameters = " CapConst : " + str(CapacityConstantOut)
        Parameters += "\t Rate : " + str(RateOut)

        return Parameters

class NADPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'NAD+ Synthesis'
        self.Input = {"nicotinamide": 1, "PRPP": 1, "ATP": 2}
        self.Output = {"NAD+": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["nicotinamide"] / (InitCond["nicotinamide"] + Molecules["nicotinamide"]) * self.CapacityConstant
        return "NAD+", VO


class NADPPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'NADP+ Synthesis'
        self.Input = {"NAD+": 1, "ATP": 1}
        self.Output = {"NADP+": 1, "ADP": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["NAD+"] / (InitCond["NAD+"] + Molecules["NAD+"]) * self.CapacityConstant
        return "NADP+", VO


class CoASynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'CoA Synthesis'
        self.Input = {"pantothenate": 1, "cysteine": 1, "ATP": 4}
        # self.Input = {"pantothenate": 1, "cysteine": 1, "ATP": 3, "CTP": 1, "CO2": 1}
        self.Output = {"CoA-SH": 1, "ADP": 2, "AMP": 1}
        # self.Output = {"CoA-SH": 1, "ADP": 2, "CMP": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["pantothenate"] / (InitCond["pantothenate"] + Molecules["pantothenate"]) * self.CapacityConstant
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
        self.ReactionName = 'ACP Synthesis'
        self.Input = {"AcpP": 1, "CoA-SH": 1}
        self.Output = {"3',5'-ADP": 1, "acyl-ACP": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["AcpP"] / (InitCond["AcpP"] + Molecules["AcpP"]) * self.CapacityConstant
        return "acyl-AcpP", VO


class IPPSynthesis(Reaction):
    # reaction IPPsynthesis(pyruvate + GAP + NADPH + CTP + ATP -> IPP + NADP+ + ADP + CMP)
    def __init__(self):
        super().__init__()
        self.ReactionName = 'IPP Synthesis'
        self.Input = {"pyruvate": 1, "GAP": 1, "NADPH": 1, "ATP": 1, "CTP": 1}
        self.Output = {"IPP": 1, "NADP+": 1, "ADP": 1, "CMP": 1}

    def Specification(self, Molecules, InitCond):
        VO = Molecules["pyruvate"] / (InitCond["pyruvate"] + Molecules["pyruvate"]) * self.CapacityConstant
        return "IPP", VO


class UbiquinoneSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Ubiquinone Synthesis'
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
        self.Input = {"PRPP": 1, "ATP": 5, "pyruvate": 1, "GAP": 1}
        # self.Output = {"TDP": 1, "AMP": 1, "ADP": 4, "formate": 1, "4-cresol": 1, "methionine": 2, "5’-deoxyadenosine": 2, "alanine": 1, "glutamate": 2}
        self.Output = {"TDP": 1, "AMP": 1, "ADP": 4, "formate": 1}
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
        self.Input = {"pimeloyl-AcpP": 1, "alanine": 1, "AMTOB": 1, "ATP": 1}
        self.Output = {"biotin": 1, "AcpP": 1, "AMP": 1, "5’-deoxyadenosine": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["pimeloyl-AcpP"] / (InitCond["pimeloyl-AcpP"] + Molecules["pimeloyl-AcpP"]) * self.CapacityConstant
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


# PRPP: Phosphoribosyl pyrophosphate
class PRPPSynthesis(Reaction):   # Oxidative PPP
    def __init__(self, Rate = 1.5e-5):
        super().__init__()
        self.ReactionName = 'PRPP Synthesis'
        self.Input = {"G6P": 1, "ATP": 1}
        # self.Output = {"PRPP": 1, "AMP": 1}
        self.Output = {"PRPP": 1, "ADP": 1}   # temporary
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "PRPP", min(self.Rate, Molecules["G6P"])

    def GetMaxConc(self, MolName, Molecules, InitCond):
        if MolName == "PRPP":
            return InitCond[MolName]
        return Reaction.MaxConc


# GAR: Glycinamide ribonucleotide
class GARSynthesis(Reaction):
    def __init__(self, Rate = 1.5e-5):
        super().__init__()
        self.ReactionName = "GAR Synthesis"
        # self.Input = {"PRPP": 1, "glutamate": 1}
        self.Input = {"PRPP": 1}
        # self.Output = {"GAR": 1, "glutamine": 1}
        self.Output = {"GAR": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "GAR", min(self.Rate, Molecules["PRPP"])

    def GetMaxConc(self, MolName, Molecules, InitCond):
        if MolName == "GAR":
            return InitCond[MolName]
        return Reaction.MaxConc


# FGAR: Formylglycinamide ribonucleotide
class FGARSynthesisByPurN(Reaction):
    def __init__(self, Rate = 1.5e-5, ExpressionFactor = 1.0):
        super().__init__()
        self.ReactionName = "FGAR Synthesis"
        self.Input = {"GAR": 1, "10-formyl-THF": 1}
        self.Output = {"FGAR": 1, "THF": 1}
        self.Rate = Rate * ExpressionFactor

    def Specification(self, Molecules, InitCond):
        return "FGAR", min(self.Rate, Molecules["GAR"])

    def GetMaxConc(self, MolName, Molecules, InitCond):
        if MolName == "FGAR":
            return InitCond[MolName]
        return Reaction.MaxConc


# FGAM: Formylglycinamidine ribonucleotide
class FGAMSynthesisByPurL(Reaction):
    def __init__(self, Rate = 1.5e-5, ExpressionFactor = 1.0):
        super().__init__()
        self.ReactionName = "FGAM Synthesis"
        # self.Input = {"FGAR": 1, "glutamine":1, "ATP": 1}
        self.Input = {"FGAR": 1, "ATP": 1}
        # self.Output = {"FGAM": 1, "glutamate": 1, "ADP": 1}
        self.Output = {"FGAM": 1, "ADP": 1}
        self.Rate = Rate * ExpressionFactor

    def Specification(self, Molecules, InitCond):
        VO = min(self.Rate, Molecules["FGAR"])
        return "FGAM", min(self.Rate, Molecules["FGAR"])

    def GetMaxConc(self, MolName, Molecules, InitCond):
        if MolName == "FGAM":
            return InitCond[MolName]
        return Reaction.MaxConc


class PurineSynthesis(Reaction):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRate * 0.25):
        super().__init__()
        self.ReactionName = "PurineSynthesis"
        self.Input = {"FGAM": 2}
        self.Output = {"dATP": 1, "dGTP": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "dATP", min(self.Rate * EcoliInfo.C2M, Molecules["FGAM"])


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


# Nucleotide Synthesis
class dNTPSynthesis(Reaction):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRate):
        super().__init__()
        self.ReactionName = 'dNTP Synthesis'
        self.Input = {}
        self.Output = {"dATP": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "dATP", self.Rate * EcoliInfo.C2M


class dCTPSynthesis(Reaction):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRate * 0.25):
        super().__init__()
        self.ReactionName = 'dCTP Synthesis'
        self.Input = {"PRPP": 1}
        self.Output = {"dCTP": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "dCTP", min(self.Rate * EcoliInfo.C2M, Molecules["PRPP"])

    
class dUTPSynthesis(Reaction):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRate * 0.25):
        super().__init__()
        self.ReactionName = 'dUTP Synthesis'
        self.Input = {}
        self.Output = {"dUTP": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "dUTP", self.Rate * EcoliInfo.C2M

    def GetMaxConc(self, MolName, Molecules, InitCond):
        if MolName == "dUTP":
            return InitCond[MolName]
        return Reaction.MaxConc


# AA Synthesis
class AASynthesis(Reaction):
    def __init__(self, Rate = EcoliInfo.ProteinSynthesisRate):
        super().__init__()
        self.ReactionName = 'AA Synthesis'
        self.Input = {}
        self.Output = {"glutamine": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "glutamine", self.Rate * EcoliInfo.C2M


class AASynthesis_Asp(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Glu Synthesis'
        self.Input = {"G6P": 1}
        self.Output = {"glutamate": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["G6P"] / (InitCond["G6P"] + Molecules["G6P"]) * self.CapacityConstant
        return "glutamate", VO


class AASynthesis_GluDerivatives(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Gln, Pro, Arg Synthesis'
        self.Input = {"glutamate": 3}
        self.Output = {"glutamine": 1, "proline": 1, "arginine": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["glutamate"] / (InitCond["glutamate"] + Molecules["glutamate"]) * self.CapacityConstant
        return "glutamine", VO


class AASynthesis_Asp(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Asp Synthesis'
        self.Input = {"G6P": 1}
        self.Output = {"aspartate": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["G6P"] / (InitCond["G6P"] + Molecules["G6P"]) * self.CapacityConstant
        return "aspartate", VO


class AASynthesis_AspDerivatives(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Asn, Met, Thr, Lys Synthesis'
        self.Input = {"aspartate": 4}
        self.Output = {"asparagine": 1, "methionine": 1, "threonine": 1, "lysine": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["aspartate"] / (InitCond["aspartate"] + Molecules["aspartate"]) * self.CapacityConstant
        return "asparagine", VO


class AASynthesis_Aromatic(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Phe, Tyr, Trp Synthesis'
        self.Input = {"G6P": 3}
        self.Output = {"phenylalanine": 1, "tyrosine": 1, "tryptophan": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["G6P"] / (InitCond["G6P"] + Molecules["G6P"]) * self.CapacityConstant
        return "phenylalanine", VO


class AASynthesis_GluDerivatives(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Gln, Pro, Arg Synthesis'
        self.Input = {"glutamate": 3}
        self.Output = {"glutamine": 1, "proline": 1, "arginine": 1}
        self.CapacityConstant = 0

    def Specification(self, Molecules, InitCond):
        VO = Molecules["glutamate"] / (InitCond["glutamate"] + Molecules["glutamate"]) * self.CapacityConstant
        return "glutamine", VO


# DHF: Dihydrofolate
class DHFSynthesis(Reaction):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRateInConc):
        super().__init__()
        self.ReactionName = "DHF Synthesis"
        self.Input = {}
        self.Output = {"DHF": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "DHF", self.Rate

    def GetMaxConc(self, MolName, Molecules, InitCond):
        if MolName == "DHF":
            return InitCond[MolName]
        return Reaction.MaxConc


# THF: Tetrahydrofolate
class THFSynthesisByFolA(Reaction):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRateInConc, ExpressionFactor = 1.0):
        super().__init__()
        self.ReactionName = "THF Synthesis"
        self.Input = {"DHF": 1}
        self.Output = {"THF": 1}
        self.Rate = Rate * ExpressionFactor

    def Specification(self, Molecules, InitCond):
        return "THF", min(self.Rate, Molecules["DHF"])

    def GetMaxConc(self, MolName, Molecules, InitCond):
        if MolName == "THF":
            return (Molecules[MolName] + Molecules["DHF"]) / 2
        return Reaction.MaxConc


# 5,10-methylene-THF
class FiveTenMethyleneTHFSynthesisByGlyA(Reaction):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRate * 2, ExpressionFactor = 1.0):
        super().__init__()
        self.ReactionName = "5,10-methylene-THF Synthesis"
        # DL - self.Input = {"THF": 1, "serine": 1}
        # DL - self.Output = {"5,10-methylene-THF": 1, "glycine": 1}
        self.Input = {"THF": 1}
        self.Output = {"5,10-methylene-THF": 1}
        self.Rate = Rate * ExpressionFactor

    def Specification(self, Molecules, InitCond):
        return "5,10-methylene-THF", min(self.Rate, Molecules["THF"])

    def GetMaxConc(self, MolName, Molecules, InitCond):
        if MolName == "5,10-methylene-THF":
            return (Molecules[MolName] + Molecules["THF"]) / 2
        return Reaction.MaxConc


# 10-formyl-THF
class TenFormylTHFSynthesis(Reaction):
    def __init__(self, Rate = 1.5e-5):
        super().__init__()
        self.ReactionName = "10-formyl-THF Synthesis"
        self.Input = {"5,10-methylene-THF": 1}
        self.Output = {"10-formyl-THF": 1}
        self.Rate = Rate

    def Specification(self, Molecules, InitCond):
        return "10-formyl-THF", min(self.Rate, Molecules["5,10-methylene-THF"])

    def GetMaxConc(self, MolName, Molecules, InitCond):
        if MolName == "10-formyl-THF":
            return (Molecules[MolName] + Molecules["5,10-methylene-THF"]) / 2
        return Reaction.MaxConc


class dTTPSynthesisByThyA(Reaction):
    def __init__(self, Rate = EcoliInfo.DNAReplicationRate * 0.25, ExpressionFactor = 1.0):
        super().__init__()
        self.ReactionName = "dTTP Synthesis"
        # DK - self.Input = {"dUMP": 1, "5-methyl-THF": 1}
        # DK - self.Output = {"dTMP": 1}
        self.Input = {"dUTP": 1, "5,10-methylene-THF": 1}
        self.Output = {"dTTP": 1, "DHF": 1}
        self.Rate = Rate * ExpressionFactor

    def Specification(self, Molecules, InitCond):
        return "dTTP", min(self.Rate * EcoliInfo.C2M, Molecules["5,10-methylene-THF"])


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
    def __init__(self, Rate = EcoliInfo.DNAReplicationRate, BuildingBlocks = ["dATP", "dCTP", "dGTP", "dTTP"]):
        super().__init__()
        self.ReactionName = "DNA Replication"
        self.EnergyConsumption = EcoliInfo.ATPConsumptionPerdNTPExtension
        self.Rate = Rate
        self.Progress = 0.0
        self.MaxProgress = EcoliInfo.GenomeSize

        self.BuildingBlocks = BuildingBlocks
        assert len(self.BuildingBlocks) > 0, 'No dNTP consumed in DNA replication'
        for buildingblock in self.BuildingBlocks:
            self.Input[buildingblock] = 1

        self.Input["ATP"] = self.EnergyConsumption * len(self.BuildingBlocks)
        self.Output["ADP"] = self.EnergyConsumption * len(self.BuildingBlocks)

    def Specification(self, Molecules, InitCond):
        return self.BuildingBlocks[0], -self.Rate * EcoliInfo.C2M / len(self.BuildingBlocks)

    def Callback(self, dMolecules):
        assert self.BuildingBlocks[0] in dMolecules
        dElongation = -dMolecules[self.BuildingBlocks[0]] * EcoliInfo.M2C * len(self.BuildingBlocks)
        assert dElongation >= 0
        self.Progress += dElongation


class ProteinSynthesis(Process):
    def __init__(self, Rate = EcoliInfo.ProteinSynthesisRate, BuildingBlocks = []):
        super().__init__()
        self.ReactionName = "Protein Synthesis"
        self.EnergyConsumption = EcoliInfo.ATPConsumptionPerAAExtension
        self.Rate = Rate
        self.Progress = 0.0
        self.MaxProgress = EcoliInfo.ProteomeSize

        self.BuildingBlocks = BuildingBlocks
        assert len(BuildingBlocks) > 0, 'No AA consumed in DNA replication'
        for buildingblock in BuildingBlocks:
            self.Input[buildingblock] = 1

        self.Input["ATP"] = self.EnergyConsumption * len(self.BuildingBlocks)
        self.Output["ADP"] = self.EnergyConsumption * len(self.BuildingBlocks)

    def Specification(self, Molecules, InitCond):
        return "glutamine", -self.Rate * EcoliInfo.C2M / len(self.BuildingBlocks)

    def Callback(self, dMolecules):
        assert "glutamine" in dMolecules
        dElongation = -dMolecules["glutamine"] * EcoliInfo.M2C * len(self.BuildingBlocks)
        assert dElongation >= 0
        self.Progress += dElongation


class VolumeExpansion(Process):
    def __init__(self,
                 ExcludedMolecules = [],
                 Rate = EcoliInfo.VolumeExpansionRate,
                 MaxProgress = EcoliInfo.VolumeExpansionSize):
        super().__init__()
        self.ReactionName = "Cell Volume Expansion"
        self.Rate = Rate
        self.Progress = 0.0
        self.MaxProgress = MaxProgress
        self.T = self.MaxProgress / self.Rate
        self.ln2_over_T = math.log(2, math.e) / self.T
        self.ExcludedMolecules = ExcludedMolecules

        self.Input = {"ATP": 1}
        self.Output = {"ADP": 1}
        for Mol in ExcludedMolecules:
            self.Input[Mol] = 1

    def Specification(self, Molecules, InitCond):
        dATPConc = -self.Rate * EcoliInfo.C2M
        for Mol in self.ExcludedMolecules:
            InitConc = InitCond[Mol]
            t = self.Progress / self.MaxProgress * self.T
            dMolConc = -InitConc * self.ln2_over_T * math.pow(math.e, -self.ln2_over_T * t)
            self.Stoich[Mol] = -dMolConc / dATPConc

        return "ATP", dATPConc


class ReactionSimulator(FSimulator):
    def __init__(self):
        super().__init__()
        self.Reactions = []
        self.Dataset = {}

        self.KnownMolConc = EcoliInfo.OpenKnownMolConc()

        self.Molecules = dict()
        self.InitialConditions = dict()   # Subset of self.KnownMolConc
        self.PermanentMolecules = list()

        self.dMolecules = dict()
        self.SortedReactionNames = list()

        self.Plot = False
        self.ExportToPDF = False
        self.Debug_Info = 1
        self.Debug_Reaction = False

    def SetPermanentMolecules(self, PermanentMolecules):
        self.PermanentMolecules = PermanentMolecules

    def Initialize(self, UserSetInitialMolecules={}):
        self.InitializeStoich()
        self.InitializeMolecules(UserSetInitialMolecules)
        self.InitializeDataset()
        self.ApplyGlobalKineticCapacityScale()

        self.SortedReactionNames = self.GetSortedReactionNames()

    def SortReactions(self):
        SortedReactions = []
        ReactionNames = []
        RxnIdx = {}

        for i, reaction in enumerate(self.Reactions):
            RxnIdx[reaction.ReactionName] = i
            ReactionNames.append(reaction.ReactionName)

        ReactionNames.sort()
        for rxnname in ReactionNames:
            rxnidx = RxnIdx[rxnname]
            SortedReactions.append(self.Reactions[rxnidx])

        self.Reactions = SortedReactions

    def GetSortedReactionNames(self):
        ReactionNames = []
        for i, reaction in enumerate(self.Reactions):
            ReactionNames.append(reaction.ReactionName)
        ReactionNames.sort()
        return ReactionNames

    def InitializeStoich(self):
        for Reaction in self.Reactions:
            for Mol, Coeff in Reaction.Input.items():
                Reaction.Stoich[Mol] = -Coeff
            for Mol, Coeff in Reaction.Output.items():
                Reaction.Stoich[Mol] = Coeff

    def GetKnownMolConc(self, Molecule):

        def EstimateMolConc(Molecule):
            if Molecule == "FADH2":
                return (self.KnownMolConc["NADH"][0] / self.KnownMolConc["NAD+"][0]) * self.KnownMolConc["FAD"][0]
            elif Molecule == "MolName":
                return
            else:
                return 0

        if Molecule in self.KnownMolConc:
            return self.KnownMolConc[Molecule][0]
        else:
            EstConc = EstimateMolConc(Molecule)
            if EstConc > 0:
                print("[Known Concentration Not Found] {:>20} concentration is set to {} (Estimated)".format(Molecule, Conc2Str(EstConc)))
                return EstConc
            else:
                DefaultConc = 0.00001
                print("[Known Concentration Not Found] {:>20} concentration is set to {} (Default)".format(Molecule, Conc2Str(DefaultConc)))
                return DefaultConc   # Default 0.01 mM

    def CheckIfInKnownMolConc(self, Molecule):
        return True if Molecule in self.KnownMolConc else False

    def GetReactionNames(self, Eq=False, Alignment=False):
        ReactionNames = ""
        for Reaction in self.Reactions:
            ReactionLabel = ""
            ReactionName = "[" + Reaction.ReactionName + "]"
            if Eq:
                ReactionEq = ""
                # ReactionEq = Reaction.GetChemicalEquationStr()
                ReactionParameters = Reaction.GetAllParameters()
                if Alignment:
                    ReactionLabel += "{:<35} {}".format(ReactionName, ReactionEq) + ReactionParameters
                else:
                    ReactionLabel += ReactionName + ReactionEq + ReactionParameters
            else:
                ReactionLabel += ReactionName
            ReactionNames += ReactionLabel + "\n"
        return str(ReactionNames)

    def PrintReactions(self):
        print('-- Reactions --')
        print(self.GetReactionNames(Eq=True, Alignment=True))

    def InitializeMolecules(self, UserSetInitialMolecules={}):
        AllMolecules = dict()
        for Reaction in self.Reactions:
            for Molecule, Coeff in Reaction.Stoich.items():
                if Molecule not in AllMolecules:
                    AllMolecules[Molecule] = self.GetKnownMolConc(Molecule)

        self.InitialConditions = dict(sorted(AllMolecules.items()))
        self.Molecules = self.InitialConditions.copy()
        for Molecule, Conc in UserSetInitialMolecules.items():
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
        RefConc = self.Molecules[RefMol]

        # Compare dConc of reference molecule to input concentrations and adjust reference dConc
        RefCoeff = Reaction.Stoich[RefMol]
        UnitdConc = RefdConc / RefCoeff
        Out = UnitdConc

        # DL: Debug
        Debug_AdjustRefdCon = False
        Debug = ""
        Debug_Min = ""
        Debug_Max = ""
        MolsInStoich = [name for name in Reaction.Input.keys()]
        if Debug_AdjustRefdCon:
            Debug = "{:<10} | {:<10} {:>10}".format(Reaction.ReactionName[:10], RefMol[:10], Conc2Str(UnitdConc))
            MolsInStoich.sort()

        for i in range(len(Reaction.Input)):
            Mol = MolsInStoich[i]
            Coeff = Reaction.Input[Mol]

            MinConc = self.Molecules_Min[Mol] if Mol in self.Molecules_Min else 0.0
            assert self.Molecules[Mol] >= MinConc * 0.999999
            AdjustedConc = (self.Molecules[Mol] - MinConc) / Coeff
            Out = min(Out, AdjustedConc)

            if Debug_AdjustRefdCon:
                Debug_Min += "\n{:>23}  min {:>6} {:>10} {:>10} {:>10} {:>10}".format(Mol[:10], Coeff, Conc2Str(MinConc), Conc2Str(self.Molecules[Mol]), Conc2Str(AdjustedConc), Conc2Str(Out))

            # The following is possible due to multiplication/division in floating numbers
            if Out * Coeff > self.Molecules[Mol] - MinConc:
                Out *= 0.999999

                if Debug_AdjustRefdCon:
                    Debug_Min += " float adjusted"

            assert Out * Coeff <= self.Molecules[Mol] - MinConc

        MolsInStoich = [name for name in Reaction.Output.keys()]
        if Debug_AdjustRefdCon:
            Debug += " -> {:>10}".format(Conc2Str(Out))
            MolsInStoich.sort()

        for i in range(len(Reaction.Output)):
            Mol = MolsInStoich[i]
            Coeff = Reaction.Output[Mol]

            MaxConc = Reaction.GetMaxConc(Mol, self.Molecules, self.InitialConditions)
            AdjustedConc = max(0, (MaxConc - self.Molecules[Mol]) / Coeff)
            Out = min(Out, AdjustedConc)

            if Debug_AdjustRefdCon:
                Debug_Max += "\n{:>23}  max {:>6} {:>10} {:>10} {:>10} {:>10}".format(Mol[:10], Coeff, Conc2Str(MaxConc), Conc2Str(self.Molecules[Mol]), Conc2Str(AdjustedConc), Conc2Str(Out))

            # The following is possible due to multiplication/division in floating numbers
            if Out * Coeff > max(0, MaxConc - self.Molecules[Mol]):
                Out *= 0.999999

                if Debug_AdjustRefdCon:
                    Debug_Max += " float adjusted"

            assert Out * Coeff <= max(0, MaxConc - self.Molecules[Mol])

        if Debug_AdjustRefdCon:
            Debug += " -> {:>10}".format(Conc2Str(Out))
            print(Debug, Debug_Min, Debug_Max)

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

    def ApplyGlobalKineticCapacityScale(self):
        for Reaction in self.Reactions:
            Reaction.CapacityConstant *= GlobalKineticScale

    def Info(self):
        HeadStr = "{:<18}: {:>10} {:>10} |".format("", "Initial", "Current")
        for i in range(len(self.dMolecules)):
            HeadStr += " {:<10.10}".format(self.SortedReactionNames[i])
        print(HeadStr)
        for Molecule, Conc in self.Molecules.items():
            InitConc = self.InitialConditions[Molecule]
            InitConcStr = Conc2Str(InitConc)
            ConcStr = Conc2Str(Conc)
            LineStr = "{:<18}: {:>10} {:>10} |".format(Molecule, InitConcStr, ConcStr)
            for i in range(len(self.dMolecules)):
                dMolecules = self.dMolecules[self.SortedReactionNames[i]]
            # for ReactionName, dMolecules in self.dMolecules.items():
                if Molecule in dMolecules:
                    Conc2 = dMolecules[Molecule]
                    ConcStr2 = Conc2Str(Conc2)
                else:
                    ConcStr2 = ""
                LineStr += " {:>10}".format(ConcStr2)
            print(LineStr)

    def Summary(self):
        print("{:<20} {:>12} {:>12} {:>12}  {:<12}".format("", "Initial", "Final", "dConc", "OutOfRange"))
        for Molecule, Conc in self.Molecules.items():
            dConc = Conc - self.InitialConditions[Molecule]

            InitConc = Conc2Str(self.InitialConditions[Molecule])
            FinalConc = Conc2Str(Conc)
            FinaldConc = Conc2Str(dConc) if Molecule not in self.PermanentMolecules else "-"
            if dConc < 0:
                FinaldConc = "\033[91m{:>12}\033[00m".format(FinaldConc)
            elif dConc > 0:
                FinaldConc = "\033[92m{:>12}\033[00m".format(FinaldConc)

            OutOfRange = ''
            if self.CheckIfInKnownMolConc(Molecule):
                if Conc < self.KnownMolConc[Molecule][1]:
                    OutOfRange = "\033[91m{:<12}\033[00m".format('DEPLETED')
                elif Conc > self.KnownMolConc[Molecule][2]:
                    OutOfRange = "\033[92m{:<12}\033[00m".format('ACCUMULATED')
            else:
                OutOfRange = "{:<12}".format('-')

            print("{:20} {:>12} {:>12} {:>12} {:<12}".format(Molecule, InitConc, FinalConc, FinaldConc, OutOfRange))

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
        if ReactionName == "Glycolysis":
            self.AddReaction(Glycolysis())
            self.AddReaction(ATPControl(-1e-5))

        elif ReactionName == "PyruvateOxidation":
            self.AddReaction(PyruvateOxidation())
            self.AddReaction(MolControl("acetyl-CoA", -3e-6))

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
    KnownMolConc = EcoliInfo.OpenKnownMolConc()
    Sim = ReactionSimulator()
    
    RunUnitTest = False
    # RunUnitTest = True

    ATPConsumption_Sec = EcoliInfo.ECM_CellDivision_Sec
    # ATPConsumption_Sec = 0

    if RunUnitTest:
        UnitTestReactions = "Glycolysis"
        Sim.UnitTest(UnitTestReactions)

    else:
        # Energy Consumption
        # Sim.AddReaction(ATPControl(-ATPConsumption_Sec))

        # Central Carbon Metabolism
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

        # Cofactors
        # Sim.AddReaction(NADPlusSynthesis())
        # Sim.AddReaction(NADPPlusSynthesis())
        # Sim.AddReaction(CoASynthesis())
        # Sim.AddReaction(ACPSynthesis())
        # Sim.AddReaction(IPPSynthesis())
        # Sim.AddReaction(UbiquinoneSynthesis())
        # Sim.AddReaction(ThiamineSynthesis())
        # Sim.AddReaction(RiboflavinSynthesis())
        # Sim.AddReaction(PantotheateSynthesis())
        # Sim.AddReaction(PyridoxalPhosphateSynthesis())
        # Sim.AddReaction(BiotinSynthesis())
        # Sim.AddReaction(FolateSynthesis())

        # Nucleotide Synthesis
        # Sim.AddReaction(OxidativePPP())
        Sim.AddReaction(PRPPSynthesis())

        Sim.AddReaction(dNTPSynthesis())
        Sim.AddReaction(DNAReplication(BuildingBlocks=["dATP"]))

        Sim.AddReaction(AASynthesis())
        Sim.AddReaction(ProteinSynthesis(BuildingBlocks=["glutamine"]))

        # Amino

        # Set permanent molecules
        PermanentMolecules = [
            "G6P",
        ]
        Sim.SetPermanentMolecules(PermanentMolecules)

    # Debugging options
    # Sim.Debug_Reaction = True
    # Sim.Debug_Info = 10000
    Sim.Plot = True
    # Sim.ExportToPDF = True

    # Set initial molecule concentrations
    UserSetInitialMolecules = {}
    # UserSetMolecules["ATP"] = 1.0 * 1e-3
    # UserSetMolecules["ADP"] = 8.6 * 1e-3
    # UserSetMolecules["G6P"] = 50 * 1e-3

    Sim.PrintReactions()

    # Execute simulation
    Sim.Initialize(UserSetInitialMolecules)

    # Simulation parameters
    TotalTime = Sim.Molecules["G6P"] * 32 / max(1e-3, ATPConsumption_Sec) + 200
    DeltaTime = 0.01

    Sim.Simulate(TotalTime=TotalTime, DeltaTime=DeltaTime)

    # Plot
    Datasets = Sim.GetDataset()
    Plot = FPlotter()

    if Sim.Plot:
        Plot.SetKnownMolConc(EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Unitless=False, Multiscale=True, Log='e')
        # Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, All=True, Multiscale=True, Include_nM=True)
        Plot.SetKnownMolConc(EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='right', Unitless=False, All=False, Individual=True, MolRange=True)

    if Sim.ExportToPDF:
        Plot.SetKnownMolConc(EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='both', Export='pdf')
        Plot.SetKnownMolConc(EcoliInfo.OpenKnownMolConc())
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, bSideLabel='right', All=False, Individual=True, MolRange=True, Export='pdf')

