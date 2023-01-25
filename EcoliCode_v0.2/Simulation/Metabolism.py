# BSD 3-Clause License
# © 2022, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import sys
from datetime import datetime
import numpy as np
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import csv
import copy

from Simulator import FSimulator

NA = 6e23
CytoVol = 1e-15
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


class FEcoliInfo():
    def __init__(self):
        # EC stands for Energy Consumption in ATP molecules (count)
        self.ECC_DNAReplication = 4.5 * 1e6 * 2  # 4.5Mbp genome (double strand)
        self.ECC_ProteinSynthesis = 3 * 1e6 * 300  # 3M proteins (each 300aa)
        self.ECC_Cytokinesis = 10 * 1e6
        self.ECC_CellDivision = self.ECC_DNAReplication + self.ECC_ProteinSynthesis + self.ECC_Cytokinesis

        # Convert count to M assuming E. coli volume is 1um^3
        self.Volume = 1e-15
        self.C2M = 1 / AvogadroNum / self.Volume

        # ECM stands for Energy Consumption in ATP (M = mol/L)
        self.ECM_CellDivision = self.ECC_CellDivision * self.C2M
        self.ECM_CellDivision_Sec = self.ECM_CellDivision / (20 * 60)  # 20 minutes assumed for one cycle of cell division

        # ECGM stands for Energy Consumption in Glucose (M)
        self.ECGM_CellDivision = self.ECM_CellDivision / 32

        # EGM stands for Energy Generation in ATP (M)
        self.EGM_Glycolysis_Sec = self.ECM_CellDivision_Sec * 10  # Glycolysis is assumed to be fast (10 times faster than necessary for cell division)
        self.EGM_TCA_Sec = self.ECM_CellDivision_Sec * 1.2
        self.EGM_OxidativePhosphorylation_Sec = self.ECM_CellDivision_Sec * 1.2

    def Info(self):
        print("Molecule Consumption - Glucose:                                {:>10}".format(
            Conc2Str(self.ECGM_CellDivision)))
        print("Energy Consumption   - Cell Division per sec:                  {:>10}".format(
            Conc2Str(self.ECM_CellDivision_Sec)))
        print("Energy Generation    - Glycolysis per sec:                     {:>10}".format(
            Conc2Str(self.EGM_Glycolysis_Sec)))
        print("Energy Generation    - TCA & Oxidative Phosporylation per sec: {:>10}".format(
            Conc2Str(self.EGM_TCA_Sec)))
        print("")


class FPlotter:
    def __init__(self):
        self.Filter_Inclusion = None
        self.Filter_Exclusion = None
        self.KnownMolConc = dict()
        self.ConsistentColorDict = dict()

    def SetKnownMolConc(self, KnownMolConc):
        self.KnownMolConc = KnownMolConc
        # https://stackoverflow.com/questions/12957582/plot-yerr-xerr-as-shaded-region-rather-than-error-bars

    def ResetFilters(self):
        self.Filter_Inclusion = None
        self.Filter_Exclusion = None

    def SetFilters(self, InclusionList, ExclusionList):
        self.ResetFilters()
        self.SetFilter_Inclusion(InclusionList)
        self.SetFilter_Exclusion(ExclusionList)

    def SetFilter_Inclusion(self, List):
        self.Filter_Inclusion = List

    def SetFilter_Exclusion(self, List):
        self.Filter_Exclusion = List

    def CheckToIncludeOrExclude(self, Key_Data):
        if self.Filter_Inclusion or self.Filter_Exclusion:
            if self.Filter_Inclusion:
                if (Key_Data in self.Filter_Inclusion) or (Key_Data[1:] in self.Filter_Inclusion):
                    return True
                else:
                    return False
            else:
                if (Key_Data in self.Filter_Exclusion) or (Key_Data[1:] in self.Filter_Exclusion):
                    return False
                else:
                    return True
        else:
            return True

    def FilterDatasets(self, Datasets):
        Datasets_Filtered = dict()
        for Key_Dataset, Dataset in Datasets.items():
            Dataset_Filtered = dict()
            for Key_Data, Data in Dataset.items():
                if self.CheckToIncludeOrExclude(Key_Data):
                    Dataset_Filtered[Key_Data] = Data
                    Datasets_Filtered[Key_Dataset] = Dataset_Filtered

        return Datasets_Filtered

    def PlotDatasets(self, Datasets, DeltaTime=1.0, bSideLabel=True, SuperTitle="", All=False, Multiscale=False, Individual=False, MolRange=False, Export=''):
        # Filter Datasets
        if self.Filter_Inclusion or self.Filter_Exclusion:
            Datasets = self.FilterDatasets(Datasets)

        def ExtractTime(Datasets, SimulationTimeUnit):
            for n, (Process, Dataset) in enumerate(Datasets.items()):
                for Data in Dataset.values():
                    Time = [i * SimulationTimeUnit for i in range(len(Data))]
                    return Time

        def ApplyGlobalUnit(Dataset):
            Unit = 'a.u.'
            array_unit = {
                'nM': 1e-9,
                'uM': 1e-6,
                'mM': 1e-3
            }

            DatasetArray = np.empty([0, 0])
            DatasetKeyIndex = dict()
            for i, (Key, Data) in enumerate(Dataset.items()):
                DatasetKeyIndex[Key] = i
                if DatasetArray.shape[0] == 0:
                    DatasetArray = np.array(Data, ndmin=2)
                else:
                    DatasetArray = np.vstack([DatasetArray, np.array(Data, ndmin=2)])

            for UnitText, UnitValue in array_unit.items():
                if np.any(DatasetArray > UnitValue):
                    # print('Unit has been set to', UnitText)
                    Unit = UnitText

            # Convert data to appropriate unit
            if Unit in array_unit:
                DatasetArray = DatasetArray / array_unit[Unit]

                # print('Final unit:', Unit)
                for Key, i in DatasetKeyIndex.items():
                    Dataset[Key] = DatasetArray[i].tolist()
                    AdjustUnitOnKnownConc(Key, array_unit[Unit])

            return Unit, Dataset

        def GetSubPlottingInfo(Datasets, MaxNPlotsInRows=4):
            ScaleFactor = 0
            if All:
                ScaleFactor += 1
            if Multiscale:
                ScaleFactor += 3
            if Individual:
                for Dataset in Datasets.values():
                    ScaleFactor += len(Dataset)

            NPlotsInRows = ScaleFactor   # Default
            if ScaleFactor > MaxNPlotsInRows * 2:
                NPlotsInRows = MaxNPlotsInRows
            else:
                for Remainder in range(MaxNPlotsInRows):
                    if ScaleFactor % (Remainder + 1) == 0:
                        NPlotsInRows = Remainder + 1

            NPlotsInColumn = math.ceil(ScaleFactor / NPlotsInRows)

            return ScaleFactor, NPlotsInRows, NPlotsInColumn

        def AdjustUnitOnKnownConc(MolName, Factor):
            self.KnownMolConc[MolName][0] /= Factor
            self.KnownMolConc[MolName][1] /= Factor
            self.KnownMolConc[MolName][2] /= Factor

        def GetYMax(Datasets):
            YMax = 0
            for Process, Dataset in Datasets.items():
                Dataset_Copy = Dataset.copy()
                UnitTxt, Dataset_Copy = ApplyGlobalUnit(Dataset_Copy)
                # Y axis (molecular concentrations)
                for MolName, Conc in Dataset_Copy.items():
                    YMaxData = max(Conc) * 1.1
                    if YMaxData > YMax:
                        YMax = YMaxData

        def PlotDataset(Dataset, Process, UnitTxt):
            # Y axis (molecular concentrations)
            YMax = 0
            for MolName, Conc in Dataset.items():

                line, = ax1.plot(Time, Conc, label="[" + MolName + "]")

                # Save assigned color for each molecule
                if MolName not in self.ConsistentColorDict:
                    self.ConsistentColorDict[MolName] = line.get_color()
                else:
                    line._color = self.ConsistentColorDict[MolName]

                # Display Molecule Range
                if MolRange:
                    ax1.fill_between(Time, self.KnownMolConc[MolName][1], self.KnownMolConc[MolName][2], alpha=0.2, facecolor=self.ConsistentColorDict[MolName])

                # Mol labeling on the curve
                if bSideLabel:
                    SelectedTimeFrameFromLeft = 0.99
                    ax1.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="right", va="bottom", color=self.ConsistentColorDict[MolName])
                    SelectedTimeFrameFromLeft = 0
                    ax1.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="left", va="bottom", color=self.ConsistentColorDict[MolName])

                YMaxData = max(Conc) * 1.1
                if YMaxData > YMax:
                    if Individual:
                        YMax = YMaxData
                    else:
                        YMax = min(1000 * 1.1, YMaxData)

            ax1.set_title(Process)
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Conc (' + UnitTxt + ')')
            ax1.set_ylim(ymin=0)
            ax1.set_ylim(ymax=YMax)

            if not bSideLabel:
                ax1.legend(loc='upper left')
            # ax1.grid()

        def PrintToPDF():
            DateTime = datetime.now().strftime('%Y%m%d-%H%M')
            PlotType = ''
            if All:
                PlotType += '_All'
            if Multiscale:
                PlotType += '_Multiscale'
            if Individual:
                PlotType += '_Individual'
            PDFFileName = 'MolConc_' + DateTime + PlotType +'.pdf'
            pdf = PdfPages(PDFFileName)
            pdf.savefig()
            pdf.close()

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.5, hspace=0.5)
        if SuperTitle:
            fig.suptitle(SuperTitle, fontsize=14)

        # # global ymax
        # if not YMax:
        #     YMax = GetYMax(Datasets)

        # Extract Global Time from the Dataset
        Time = ExtractTime(Datasets, DeltaTime)

        # Determine subplotting map
        ScaleFactor, NPlotsInRows, NPlotsInColumn = GetSubPlottingInfo(Datasets)

        # Plot data
        SubplotID = 1
        if All:
            for n, (Process, Dataset) in enumerate(Datasets.items()):
                UnitTxt, Dataset = ApplyGlobalUnit(Dataset)
                ax1 = fig.add_subplot(NPlotsInColumn, NPlotsInRows, SubplotID)
                PlotDataset(Dataset, Process, UnitTxt)
                SubplotID += 1

        # Multiscale split datasets
        if Multiscale:
            def SplitScales(Dataset):
                Dataset_mM = dict()
                Dataset_uM = dict()
                Dataset_nM = dict()
                for mol, conc in Dataset.items():
                    if conc[0] < 1e-6:
                        Dataset_nM[mol] = (np.array(conc) / 1e-9).tolist()
                        AdjustUnitOnKnownConc(mol, 1e-9)
                    elif conc[0] < 1e-3:
                        Dataset_uM[mol] = (np.array(conc) / 1e-6).tolist()
                        AdjustUnitOnKnownConc(mol, 1e-6)
                    else:
                        Dataset_mM[mol] = (np.array(conc) / 1e-3).tolist()
                        AdjustUnitOnKnownConc(mol, 1e-3)
                return (Dataset_mM, Dataset_uM, Dataset_nM)

            for n, (Process, Dataset) in enumerate(Datasets.items()):
                Dataset_Unit = SplitScales(Dataset)
                Units = ['mM', 'uM', 'nM']

                for Unit, dataset in zip(Units, Dataset_Unit):
                    ax1 = fig.add_subplot(NPlotsInColumn, NPlotsInRows, SubplotID)
                    PlotDataset(dataset, Unit + ' range', Unit)
                    SubplotID += 1

        # Individual Datasets
        if Individual:
            def ConvertToIndividualDataset(Dataset):
                AdjustedDatasets = dict()
                Units = dict()
                for mol, conc in Dataset.items():
                    if conc[0] < 1e-6:
                        AdjustedDatasets[mol] = {mol: (np.array(conc) / 1e-9).tolist()}
                        Units[mol] = "nM"
                        AdjustUnitOnKnownConc(mol, 1e-9)
                    elif conc[0] < 1e-3:
                        AdjustedDatasets[mol] = {mol: (np.array(conc) / 1e-6).tolist()}
                        Units[mol] = "uM"
                        AdjustUnitOnKnownConc(mol, 1e-6)
                    else:
                        AdjustedDatasets[mol] = {mol: (np.array(conc) / 1e-3).tolist()}
                        Units[mol] = "mM"
                        AdjustUnitOnKnownConc(mol, 1e-3)
                return Units, AdjustedDatasets

            for n, (Process, Dataset) in enumerate(Datasets.items()):
                Unit, AdjustedDatasets = ConvertToIndividualDataset(Dataset)
                for MolName, Dataset in AdjustedDatasets.items():
                    ax1 = fig.add_subplot(NPlotsInColumn, NPlotsInRows, SubplotID)
                    PlotDataset(Dataset, MolName, Unit[MolName])
                    SubplotID += 1

        if Export == 'pdf':
            PrintToPDF()

        plt.show()

    # def PrintToPDF(self):
    #     # Create the PdfPages object to which we will save the pages:
    #     # The with statement makes sure that the PdfPages object is closed properly at
    #     # the end of the block, even if an Exception occurs.
    #     with PdfPages('multipage_pdf.pdf') as pdf:
    #         plt.figure(figsize=(3, 3))
    #         plt.plot(range(7), [3, 1, 4, 1, 5, 9, 2], 'r-o')
    #         plt.title('Page One')
    #         pdf.savefig()  # saves the current figure into a pdf page
    #         plt.close()
    #
    #         # if LaTeX is not installed or error caught, change to `False`
    #         plt.rcParams['text.usetex'] = True
    #         plt.figure(figsize=(8, 6))
    #         x = np.arange(0, 5, 0.1)
    #         plt.plot(x, np.sin(x), 'b-')
    #         plt.title('Page Two')
    #         pdf.attach_note("plot of sin(x)")  # attach metadata (as pdf note) to page
    #         pdf.savefig()
    #         plt.close()
    #
    #         plt.rcParams['text.usetex'] = False
    #         fig = plt.figure(figsize=(4, 5))
    #         plt.plot(x, x ** 2, 'ko')
    #         plt.title('Page Three')
    #         pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    #         plt.close()
    #
    #         # We can also set the file's metadata via the PdfPages object:
    #         d = pdf.infodict()
    #         d['Title'] = 'Multipage PDF Example'
    #         d['Author'] = 'Jouni K. Sepp\xe4nen'
    #         d['Subject'] = 'How to create a multipage pdf file and set its metadata'
    #         d['Keywords'] = 'PdfPages multipage keywords author title subject'
    #         d['CreationDate'] = datetime.datetime(2009, 11, 13)
    #         d['ModDate'] = datetime.datetime.today()

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

    # Specification Models
    def MassAction_ReactantHomeostasis(self, Reactant, Molecules, InitCond):
        '''
        mass action driven by the rate-limiting reactant, maintaining the initial condition of the reactant
        '''
        return max(0, Molecules[Reactant] - InitCond[Reactant]) * Molecules[Reactant] * self.Capacity

    def ProductInhibition_ProductHomeostasis(self, Product, Molecules, InitCond):
        '''
        product inhibition maintaining the initial condition of the product
        '''
        return max(0, InitCond[Product] - Molecules[Product]) / Molecules[Product] * self.Capacity

    def ProductInhibition_ProductHomeostasis_MassAction(self, Reactant, Product, Molecules, InitCond):
        '''
        product inhibition, maintaining the initial condition of the product
        + mass action driven by the rate-limiting reactant
        '''
        return max(0, InitCond[Product] - Molecules[Product]) * Molecules[Reactant] / Molecules[Product] * self.Capacity

    def Homeostasis(self, Reactant, Product, Molecules, InitCond):
        '''
        product inhibition, maintaining the initial condition of the product
        + mass action driven by the rate-limiting reactant, maintaining the initial condition of the reactant
        '''
        return Molecules[Reactant] / Molecules[Product] * self.Capacity
        # return (max(0, Molecules[Reactant] - InitCond[Reactant]) + max(0, InitCond[Product] - Molecules[Product])) * Molecules[Reactant] / Molecules[Product] * self.Capacity

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

    def DisplayChemicalEquation(self):

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

# reaction NAD+Synthesis(nicotinamide + PRPP + 2 ATP -> NAD+)
class NADPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"nicotinamide": 1, "PRPP": 1, "ATP": 2}
        self.Output = {"NAD+": 1}

    def Specification(self, Molecules, InitCond):
        return {}

# reaction NADP+Synthesis(NAD+ + ATP -> NADP+ + ADP)
class NADPPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"NAD+": 1, "ATP": 1}
        self.Output = {"NADP+": 1, "ADP": 1}

    def Specification(self, Molecules, InitCond):
        return {}

# reaction CoASynthesis(pantothenate + cysteine + 3 ATP + CTP + CO2 -> CoA + ADP + CMP + 3 PPi)
class CoASynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"pantothenate": 1, "cysteine": 1, "ATP": 3, "CTP": 1, "CO2": 1}
        self.Output = {"CoA": 1, "ADP": 1, "CMP": 1, "PPi": 3}

    def Specification(self, Molecules, InitCond):
        return {}

class Glycolysis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Glycolysis'
        self.Input = {"G6P": 1, "ADP": 2, "NAD+": 2}
        self.Output = {"pyruvate": 2, "NADH": 2, "ATP": 2}
        self.CapacityConstant = 12.7 * 1e-3

    def Specification(self, Molecules, InitCond):
        # dATP = (Molecules["ADP"] / Molecules["ATP"]) * c
        # dATP = (abs(Molecules["ADP"] - InitCond["ADP"]) / Molecules["ATP"]) * c   # For homeostasis of ADP level
        # dATP = (Molecules["ADP"] * abs(Molecules["ATP"] - InitCond["ATP"]) / Molecules["ATP"]) * c   # For homeostasis of ATP level
        # dATP = (max(0, InitCond["ATP"] - Molecules["ATP"]) / Molecules["ATP"]) * Molecules["G6P"] * self.Capacity # Product Inhibition
        # return "ATP", self.Homeostasis("G6P", "ATP", Molecules, InitCond)
        # return "ATP", (max(0, Molecules["ADP"] - 0.5 * 1e-3) / Molecules["ATP"]) * min(1e-3, Molecules["G6P"]) * self.Capacity
        VO = Molecules["ADP"] / Molecules["ATP"] * self.CapacityConstant
        VMax = max(0, Molecules["ADP"] - InitCond["ADP"])
        self.Capacity = min(VO, VMax)
        return "ATP", self.Capacity

class PyruvateOxidation(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Pyruvate Oxidation'
        self.Input = {"pyruvate": 1, "NAD+": 1, "CoA-SH": 1}
        self.Output = {"acetyl-CoA": 1, "NADH": 1}
        self.CapacityConstant = 0.2 * 1e-3

    def Specification(self, Molecules, InitCond):
        # return "acetyl-CoA", self.Homeostasis("pyruvate", "acetyl-CoA", Molecules, InitCond)
        VO = Molecules["pyruvate"] / (InitCond["pyruvate"] + Molecules["pyruvate"]) * self.CapacityConstant
        VMax = max(0, Molecules["pyruvate"] - InitCond["pyruvate"])
        self.Capacity = min(VO, VMax)
        return "acetyl-CoA", self.Capacity


class TCACycle(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'TCA cycle'
        self.Input = {"acetyl-CoA": 1, "NAD+": 3, "FAD": 1, "ADP": 1}
        self.Output = {"NADH": 3, "FADH2": 1, "ATP": 1, "CoA-SH": 1}
        self.CapacityConstant = 0.2 * 1e-3

    def Specification(self, Molecules, InitCond):
        # self.Capacity = min(Molecules["oxaloacetate"], Molecules["acetyl-CoA"]) * 1e2
        # return "a-ketoglutarate", self.ProductInhibition_ProductHomeostasis("a-ketoglutarate", Molecules, InitCond)
        VO = Molecules["acetyl-CoA"] / (InitCond["acetyl-CoA"] + Molecules["acetyl-CoA"]) * self.CapacityConstant
        VMax = max(0, Molecules["acetyl-CoA"] - InitCond["acetyl-CoA"])
        self.Capacity = min(VO, VMax)
        return "ATP", self.Capacity


class TCACycle_alphaKetoGlutarateSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'TCA_a-ketoglutarate Synthesis'
        self.Input = {"acetyl-CoA": 1, "oxaloacetate": 1, "CoA-SH": 1, "NAD+": 1}
        self.Output = {"CoA-SH": 1, "a-ketoglutarate": 1, "NADH": 1}
        self.CapacityConstant = 1e2

    def Specification(self, Molecules, InitCond):
        self.UpdateCapacity(Molecules, ["oxaloacetate", "acetyl-CoA"])
        return "a-ketoglutarate", self.ProductInhibition_ProductHomeostasis("a-ketoglutarate", Molecules, InitCond)


class TCACycle_SuccinylCoASynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'TCA_Succinyl-CoA Synthesis'
        self.Input = {"a-ketoglutarate": 1, "NAD+": 1, "CoA-SH": 1}
        self.Output = {"succinyl-CoA": 1, "NADH": 1}
        self.Capacity = 1

    def Specification(self, Molecules, InitCond):
        return "succinyl-CoA", self.Homeostasis("a-ketoglutarate", "succinyl-CoA", Molecules, InitCond)


class TCACycle_FumarateSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'TCA_Fumarate Synthesis'
        self.Input = {"succinyl-CoA": 1, "ADP": 1, "FAD": 1}
        self.Output = {"fumarate": 1, "CoA-SH": 1, "ATP": 1, "FADH2": 1}  # TODO: Update from ATP to GTP later
        self.Capacity = 1e1

    def Specification(self, Molecules, InitCond):
        return "fumarate", self.Homeostasis("succinyl-CoA", "fumarate", Molecules, InitCond)


class TCACycle_OxaloacetateSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'TCA_Oxaloacetate Synthesis'
        self.Input = {"fumarate": 1, "NAD+": 1}
        self.Output = {"oxaloacetate": 1, "NADH": 1}
        self.Capacity = 1e-1

    def Specification(self, Molecules, InitCond):
        return "oxaloacetate", self.Homeostasis("fumarate", "oxaloacetate", Molecules, InitCond)


class PyruvateCarboxylation(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'PyruvateCarboxylation'
        self.Input = {"pyruvate": 1, "ATP": 1}
        self.Output = {"oxaloacetate": 1, "ADP": 1}
        self.Capacity = 1e-3

    def Specification(self, Molecules, InitCond):
        return "oxaloacetate", self.Homeostasis("pyruvate", "oxaloacetate", Molecules, InitCond)


class NADH_OxidativePhosphorylation(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'NADH OxidativePhosphorylation'
        self.Input = {"NADH": 1, "ADP": 2.5}
        self.Output = {"NAD+": 1, "ATP": 2.5}
        self.CapacityConstant = 1.2e-3

    def Specification(self, Molecules, InitCond):
        VO = Molecules["NADH"] / (InitCond["NADH"] + Molecules["NADH"]) * self.CapacityConstant
        VMax = max(0, Molecules["NADH"] - InitCond["NADH"])
        self.Capacity = min(VO, VMax)
        return "ATP", self.Capacity


class FADH2_OxidativePhosphorylation(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'FADH2 OxidativePhosphorylation'
        self.Input = {"FADH2": 1, "ADP": 1.5}
        self.Output = {"FAD": 1, "ATP": 1.5}
        self.CapacityConstant = 0.15e-3

    def Specification(self, Molecules, InitCond):
        VO = Molecules["FADH2"] / (InitCond["FADH2"] + Molecules["FADH2"]) * self.CapacityConstant
        VMax = max(0, Molecules["FADH2"] - InitCond["FADH2"])
        self.Capacity = min(VO, VMax)
        return "ATP", self.Capacity


class ATPControl(Reaction):
    def __init__(self, ControlRate=-4.35E-03):
        # Cell Division ATP consumption: c = 4.35E-03

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
        self.MaxRate = 0.0
        self.Progress = 0.0
        self.MaxProgress = 0.0


class DNAReplication(Process):
    def __init__(self):
        super().__init__()
        # self.BuildingBlocks = {"dATP": 1, "dCTP": 1, "dGTP": 1, "dTTP": 1}
        self.BuildingBlocks = {"dATP": 1}
        self.EnergyConsumption = 3 # 3 ATPs per 1 nucleotide extension
        self.Regulators = {}
        self.Rate = 0.0
        self.MaxRate = 4000.0   # total of 4000 bp synthesis during the exponential phase
        self.Progress = 0.0
        self.MaxProgress = 4.5e6


class ProteinSynthesis(Process):
    def __init__(self):
        super().__init__()
        # self.BuildingBlocks = {"dATP": 1, "dCTP": 1, "dGTP": 1, "dTTP": 1}
        self.BuildingBlocks = {"Glutamate": 1}
        self.EnergyConsumption = 3 # 3 ATPs per 1 nucleotide extension
        self.Regulators = {}
        self.Rate = 0.0
        self.MaxRate = 2000.0
        self.Progress = 0.0
        self.MaxProgress = 4.5e6


class ReactionSimulator(FSimulator):
    def __init__(self):
        self.Reactions = []
        self.Dataset = {}
        self.Iter = 0

        self.KnownMolConc = self.OpenKnownMolConc()

        self.Molecules = dict()
        self.InitialConditions = dict()
        self.PermanentMolecules = list()

        self.dMolecules = dict()

        self.Plot = False
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
        if self.Plot:
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

    def GetReactionNames(self, Eq=False):
        ReactionNames = ""
        for Reaction in self.Reactions:
            ReactionLabel = "[" + Reaction.ReactionName + "]"
            if Eq:
                ReactionLabel += " " + Reaction.DisplayChemicalEquation()
            ReactionNames += ReactionLabel + "\n"
        return str(ReactionNames)

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

    def AddReaction(self, Reaction):
        self.Reactions.append(Reaction)

    def GetMolCount(self, Conc):
        return Conc * EcoliInfo.Volume * AvogadroNum

    def GetConc(self, Count):
        return Count / EcoliInfo.Volume / AvogadroNum

    def AdjustRefdCon(self, Reaction, RefMol, RefdConc):
        ''' Compare dConc of reference molecule to input concentrations and adjust reference dConc '''
        RefCoeff = Reaction.Stoich[RefMol]
        UnitdConc = RefdConc / RefCoeff

        Out = 10  # 10 Molar to begin comparison
        for Mol, Coeff in Reaction.Input.items():
            AdjustedConc = self.Molecules[Mol] / Coeff

            Out = min(Out, min(AdjustedConc, UnitdConc))

            # The following is possible due to multiplication/division in floating numbers
            if Out * Coeff > self.Molecules[Mol]:
                Out *= 0.999999
            assert Out * Coeff <= self.Molecules[Mol]

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
        # DEBUG
        # print(self.Molecules, dMolecules)
        for dMolecule, dConc in dMolecules.items():
            assert dMolecule in self.Molecules
            assert dConc + self.Molecules[dMolecule] >= 0, \
                'Iter {}\t | {} \t| Conc:{}, \t dConc:{}'.format(self.Iter,
                                                                 dMolecule,
                                                                 Conc2Str(
                                                                     self.Molecules[
                                                                         dMolecule]),
                                                                 Conc2Str(
                                                                     dConc))
            self.Molecules[dMolecule] += dConc
            self.Molecules[dMolecule] = max(0, self.Molecules[dMolecule])

    def RestorePermanentMolecules(self):
        for Molecule in self.PermanentMolecules:
            if Molecule in self.Molecules:
                self.Molecules[Molecule] = self.InitialConditions[Molecule]

    def CheckZeroConcentrations(self):
        for Molecule, Conc in self.Molecules.items():
            if Conc < self.GetConc(1):  # if
                self.Molecules[Molecule] = self.GetConc(0.1)
                assert self.Molecules[Molecule] > 0

    def InitializeDataset(self):
        for Molecule, Conc in self.Molecules.items():
            self.Dataset[Molecule] = [Conc]

    def AddToDataset(self):
        for Molecule, Conc in self.Molecules.items():
            self.Dataset[Molecule].append(Conc)

    def ApplyGlobalKineticCapacity(self):
        for Reaction in self.Reactions:
            Reaction.Capacity *= GlobalKineticScale

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
        print("{:<16} {:>10} {:>10} {:>10} {:<10} {:<10}".format("", "Initial", "Final", "dConc", "Depleted", "Accumulated"))
        for Molecule, Conc in self.Molecules.items():
            dConc = Conc - self.InitialConditions[Molecule]

            InitConc = Conc2Str(self.InitialConditions[Molecule])
            FinalConc = Conc2Str(Conc)
            FinaldConc = Conc2Str(dConc) if Molecule not in self.PermanentMolecules else "-"
            Depletion = 'DEPLETED' if Conc < self.KnownMolConc[Molecule][1] else ''
            Accumulation = 'ACCUMULATED' if Conc > self.KnownMolConc[Molecule][2] else ''

            print("{:16} {:>10} {:>10} {:>10} {:<10} {:<10}".format(Molecule, InitConc, FinalConc, FinaldConc, Depletion, Accumulation))

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
        self.dMolecules = dict()
        for Reaction in self.Reactions:
            dMolecules = self.RunReaction(Reaction, DeltaTime)
            if self.Debug_Reaction:
                print("[{}] {}".format(Reaction.ReactionName, Reaction.DisplayChemicalEquation()))
                self.PrintMolConc(dMolecules, End=', ')
            self.UpdateMolecules(dMolecules)
            self.CheckZeroConcentrations()
            self.dMolecules[Reaction.ReactionName] = dMolecules

        if len(self.PermanentMolecules):
            self.RestorePermanentMolecules()

        if self.Plot:
            self.AddToDataset()        

    def GetDataset(self):
        return {self.GetReactionNames(Eq=True): self.Dataset}

    def LoadTSVDatabase(self, db_fname):
        db = None
        with open(db_fname) as fp:
            csv_reader = csv.reader(fp, delimiter='\t')
            list_of_rows = list(csv_reader)
            db = list_of_rows[1:]
        return db

    def OpenKnownMolConc(self):

        def ParseMetaboliteCon(db_KnownMolConc):

            def MetaboliteSynonyms(Name):
                if Name == "glucose-6-phosphate":
                    return "G6P"
                elif Name == "glyceraldehyde-3-phosphate":
                    return "G3P"
                elif Name == "coenzyme-A":
                    return "CoA-SH"
                # elif " (assumed 1 / 2 ile+leu)" in Name:
                #     return Name.replace(" (assumed 1 / 2 ile+leu)", "")
                elif Name == "ribose-5-phosphate":
                    return "R5P"
                elif Name == "ribulose-5-phosphate":
                    return "Ru5P"
                else:
                    return Name

            KnownMolConc = dict()

            for i, Value in enumerate(db_KnownMolConc):
                Metabolite, ConcInEcoli, LowerBound, UpperBound = Value
                # data = dict()
                # data['ConcInEcoli'] = ConcInEcoli
                # data['LowerBound'] = LowerBound
                # data['UpperBound'] = UpperBound
                # KnownMolConc[Metabolite] = data
                if ConcInEcoli == '-':
                    continue
                Metabolite = Metabolite.replace('[c]', '').replace('[m]', '')
                Metabolite = MetaboliteSynonyms(Metabolite)
                KnownMolConc[Metabolite] = [float(ConcInEcoli), float(LowerBound), float(UpperBound)]

            return KnownMolConc

        db_KnownMolConc = self.LoadTSVDatabase("MetaboliteConcentrations.tsv")
        KnownMolConc = ParseMetaboliteCon(db_KnownMolConc)

        # FADH2 missing data
        KnownMolConc["FADH2"] = [(KnownMolConc["NADH"][0] / KnownMolConc["NAD+"][0]) * KnownMolConc["FAD"][0],
                                 (KnownMolConc["NADH"][0] / KnownMolConc["NAD+"][0]) * KnownMolConc["FAD"][1],
                                 (KnownMolConc["NADH"][0] / KnownMolConc["NAD+"][0]) * KnownMolConc["FAD"][2]]

        return KnownMolConc

    def UnitTest(self, ReactionName):
        # a-ketoGlutarate synthesis unit test
        if ReactionName == "Glycolysis":
            Sim.AddReaction(Glycolysis())
            Sim.AddReaction(ATPControl(-1e-5))

        elif ReactionName == "PyruvateOxidation":
            Sim.AddReaction(PyruvateOxidation())
            Sim.AddReaction(MolControl("acetyl-CoA", -3e-6))

        elif ReactionName == "TCACycle_a-KetoGlutarateSynthesis":
            Sim.AddReaction(TCACycle_alphaKetoGlutarateSynthesis())
            Sim.AddReaction(MolControl("a-ketoglutarate", -3e-7))
            Sim.AddPermanentMolecules(["oxaloacetate"])

        elif ReactionName == "TCACycle_SuccinylCoASynthesis":
            Sim.AddReaction(TCACycle_SuccinylCoASynthesis())
            Sim.AddReaction(MolControl("succinyl-CoA", -2e-6))

        elif ReactionName == "TCACycle_FumarateSynthesis":
            Sim.AddReaction(TCACycle_FumarateSynthesis())
            Sim.AddReaction(MolControl("fumarate", -1e-5))
            Sim.AddPermanentMolecules(["succinyl-CoA", "ADP", "ATP", "FAD", "FADH2"])

        elif ReactionName == "TCACycle_OxaloacetateSynthesis":
            Sim.AddReaction(TCACycle_OxaloacetateSynthesis())
            Sim.AddReaction(MolControl("oxaloacetate", -5e-7))

        elif ReactionName == "PyruvateCarboxylation":
            Sim.AddReaction(PyruvateCarboxylation())
            Sim.AddReaction(MolControl("oxaloacetate", -5e-8))
            # Sim.AddPermanentMolecules(["succinyl-CoA", "ADP", "ATP", "FAD", "FADH2"])

        # elif ReactionName == "OxidativePhosphorylation":
        #     Sim.AddReaction(OxidativePhosphorylation())
        #     Sim.AddReaction(ATPControl(-4.35E-04))
        #     Sim.AddReaction(NADHControl(1.0874e-4))
        #     Sim.AddReaction(FADH2Control(1.0875e-4))

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

    EcoliInfo = FEcoliInfo()
    EcoliInfo.Info()
    ATPConsumption_Sec = EcoliInfo.ECM_CellDivision_Sec
    # ATPConsumption_Sec = 0

    if RunUnitTest:
        # UnitTestReactions = "OxidativePhosphorylation"
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
    if Sim.Plot:
        Datasets = Sim.GetDataset()
        Plot = FPlotter()

        Plot.SetKnownMolConc(copy.deepcopy(Sim.KnownMolConc))
        # Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, All=True)
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, All=True, Export='pdf')
        #
        Plot.SetKnownMolConc(copy.deepcopy(Sim.KnownMolConc))
        # Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, Multiscale=True)
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, Multiscale=True, Export='pdf')

        Plot.SetKnownMolConc(copy.deepcopy(Sim.KnownMolConc))
        # Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, Individual=True, MolRange=True)
        Plot.PlotDatasets(copy.deepcopy(Datasets), DeltaTime=DeltaTime, Individual=True, MolRange=True, Export='pdf')

"""
Reference

Metabolite concentrations, fluxes and free energies imply efficient enzyme usage
Junyoung O Park, Sara A Rubin, Yi-Fan Xu, Daniel Amador-Noguez, Jing Fan, Tomer Shlomi & Joshua D Rabinowitz
Nature Chemical Biology volume 12, pages482–489 (2016)
https://www.nature.com/articles/nchembio.2077#Sec21
"""
