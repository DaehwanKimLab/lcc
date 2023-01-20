# BSD 3-Clause License
# © 2022, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import sys
from datetime import datetime
import numpy as np
import math
import matplotlib.pyplot as plt
import csv

NA = 6e23
CytoVol = 1e-15

GlobalKineticScale = 1

class FPlotter:
    def __init__(self):
        self.Filter_Inclusion = None
        self.Filter_Exclusion = None

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

    def PlotDatasets(self, Datasets, DeltaTime=1.0, YMin=0, YMax=0, bSideLabel=True, SuperTitle="", Multiscale=False):
        # Filter Datasets
        if self.Filter_Inclusion or self.Filter_Exclusion:
            Datasets = self.FilterDatasets(Datasets)

        def ExtractTime(Dataset, SimulationTimeUnit):
            Time = None
            for Data in Dataset.values():
                Time = [i * SimulationTimeUnit for i in range(len(Data))]
                break
            return Time, Dataset

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

            return Unit, Dataset

        def GetSubPlottingInfo(Datasets, MaxNPlotsInRows=3, ScaleFactor=1):
            NPlotsInRows = len(Datasets) * ScaleFactor # Default
            if len(Datasets) > 1:
                for Remainder in range(MaxNPlotsInRows):
                    if len(Datasets) % (Remainder + 1) == 0:
                        NPlotsInRows = Remainder + 1
            return NPlotsInRows

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

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.5, hspace=0.5)
        if SuperTitle:
            fig.suptitle(SuperTitle, fontsize=14)

        ConsistentColorDict = dict()

        # # global ymax
        # if not YMax:
        #     YMax = GetYMax(Datasets)

        def PlotDataset(Dataset, Process, UnitTxt):
            # Y axis (molecular concentrations)
            YMax = 0
            for MolName, Conc in Dataset.items():

                line, = ax1.plot(Time, Conc, label="[" + MolName + "]")

                # Save assigned color for each molecule
                if MolName not in ConsistentColorDict:
                    ConsistentColorDict[MolName] = line.get_color()
                else:
                    line._color = ConsistentColorDict[MolName]

                # Mol labeling on the curve
                if bSideLabel:
                    SelectedTimeFrameFromLeft = 0.99
                    ax1.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="right", va="bottom", color=ConsistentColorDict[MolName])
                    SelectedTimeFrameFromLeft = 0
                    ax1.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="left", va="bottom", color=ConsistentColorDict[MolName])

                YMaxData = max(Conc) * 1.1
                if YMaxData > YMax:
                    YMax = min(1000 * 1.1, YMaxData)

            ax1.set_title(Process)
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Molecules: Concentration (' + UnitTxt + ')')
            ax1.set_ylim(ymin=0)
            ax1.set_ylim(ymax=YMax)

            if not bSideLabel:
                ax1.legend(loc='upper left')
            # ax1.grid()

        # Plot data
        if Multiscale:
            ScaleFactor = 4
            NPlotsInRows = GetSubPlottingInfo(Datasets, ScaleFactor=ScaleFactor)
            NPlotsInColumn = math.ceil(len(Datasets) / NPlotsInRows)

            def SplitScales(Dataset):
                Dataset_mM = dict()
                Dataset_uM = dict()
                Dataset_nM = dict()
                for mol, conc in Dataset.items():
                    if conc[0] < 1e-6:
                        Dataset_nM[mol] = (np.array(conc) / 1e-9).tolist()
                    elif conc[0] < 1e-3:
                        Dataset_uM[mol] = (np.array(conc) / 1e-6).tolist()
                    else:
                        Dataset_mM[mol] = (np.array(conc) / 1e-3).tolist()
                return (Dataset_mM, Dataset_uM, Dataset_nM)

            for n, (Process, Dataset) in enumerate(Datasets.items()):
                Time, Dataset = ExtractTime(Dataset, DeltaTime)
                Dataset_Unit = SplitScales(Dataset.copy())
                Units = ['mM', 'uM', 'nM']

                # All Data
                i = 1
                UnitTxt, Dataset = ApplyGlobalUnit(Dataset)
                ax1 = fig.add_subplot(NPlotsInColumn, NPlotsInRows, ScaleFactor * n + i)
                PlotDataset(Dataset, Process, UnitTxt)

                # Split scale ranges
                for Unit, dataset in zip(Units, Dataset_Unit):
                    i += 1
                    ax1 = fig.add_subplot(NPlotsInColumn, NPlotsInRows, ScaleFactor * n + i)
                    PlotDataset(dataset, Unit + ' range', Unit)


        else:
            NPlotsInRows = GetSubPlottingInfo(Datasets)

            for n, (Process, Dataset) in enumerate(Datasets.items()):
                Time, Dataset = ExtractTime(Dataset, DeltaTime)
                UnitTxt, Dataset = ApplyGlobalUnit(Dataset)
                ax1 = fig.add_subplot(math.ceil(len(Datasets) / NPlotsInRows), NPlotsInRows, n + 1)

                PlotDataset(Dataset, Process, UnitTxt)

        plt.show()

class Reaction:
    def __init__(self):
        self.ReactionName = ""
        self.Input = dict()
        self.Output = dict()
        self.Stoich = dict()
        self.Capacity = 0

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
        return (max(0, Molecules[Reactant] - InitCond[Reactant]) + max(0, InitCond[Product] - Molecules[Product])) * Molecules[Reactant] / Molecules[Product] * self.Capacity

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
        self.Output = {"pyruvate": 1, "NADH": 2, "ATP": 2}
        self.Capacity = 0.5

    def Specification(self, Molecules, InitCond):
        # dATP = (Molecules["ADP"] / Molecules["ATP"]) * c
        # dATP = (abs(Molecules["ADP"] - InitCond["ADP"]) / Molecules["ATP"]) * c   # For homeostasis of ADP level
        # dATP = (Molecules["ADP"] * abs(Molecules["ATP"] - InitCond["ATP"]) / Molecules["ATP"]) * c   # For homeostasis of ATP level
        # dATP = (max(0, InitCond["ATP"] - Molecules["ATP"]) / Molecules["ATP"]) * Molecules["G6P"] * self.Capacity # Product Inhibition
        return "ATP", self.Homeostasis("G6P", "ATP", Molecules, InitCond)

class PyruvateOxidation(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Pyruvate Oxidation'
        self.Input = {"pyruvate": 1, "NAD+": 1, "CoA-SH": 1}
        self.Output = {"acetyl-CoA": 1, "NADH": 1}
        self.Capacity = 1

    def Specification(self, Molecules, InitCond):
        return "acetyl-CoA", self.Homeostasis("pyruvate", "acetyl-CoA", Molecules, InitCond)


class TCACycle_alphaKetoGlutarateSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'TCA_a-ketoglutarate Synthesis'
        self.Input = {"acetyl-CoA": 1, "oxaloacetate": 1, "CoA-SH": 1, "NAD+": 1}
        self.Output = {"CoA-SH": 1, "a-ketoglutarate": 1, "NADH": 1}
        self.Capacity = 0

    def Specification(self, Molecules, InitCond):
        self.Capacity = min(Molecules["oxaloacetate"], Molecules["acetyl-CoA"]) * 1e2
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
        self.Output = {"fumarate": 1, "CoA-SH": 1, "ATP": 1, "FADH2": 1} # TODO: Update from ATP to GTP later
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

class OxidativePhosphorylation(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'OxidativePhosphorylation'
        self.Input = {"NADH": 1, "FADH2": 1, "ADP": 4}
        self.Output = {"NAD+": 1, "FAD": 1, "ATP": 4}
        self.Capacity = 0

    def Specification(self, Molecules, InitCond):
        self.Capacity = min(Molecules["NADH"], Molecules["FADH2"]) * 1e3
        return "ATP", self.ProductInhibition_ProductHomeostasis("ATP", Molecules, InitCond)

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


class Simulator():
    def __init__(self):
        self.Reactions = []
        self.Dataset = {}
        self.Iter = 0

        self.KnownMolConc = self.OpenKnownMolConc()

        self.Molecules = dict()
        self.InitialConditions = dict()
        self.PermanentMolecules = list()

        self.Plot = False
        self.Debug_Info = 1
        self.Debug_Reaction = False

    def AddPermanentMolecules(self, PermanentMolecules):
        for Molecule in PermanentMolecules:
            self.PermanentMolecules.append(Molecule)

    def SetPermanentMolecules(self, PermanentMolecules):
        self.PermanentMolecules = PermanentMolecules

    def Initialize(self):
        self.InitializeStoich()
        self.InitializeMolecules()
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
            return self.KnownMolConc[Molecule]
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

    def InitializeMolecules(self):
        AllMolecules = dict()
        for Reaction in self.Reactions:
            for Molecule, Coeff in Reaction.Stoich.items():
                if Molecule not in AllMolecules:
                    AllMolecules[Molecule] = self.GetKnownMolConc(Molecule)

        self.Molecules = dict(sorted(AllMolecules.items()))
        self.InitialConditions = self.Molecules.copy()

    def AddReaction(self, Reaction):
        self.Reactions.append(Reaction)

    def GetMolCount(self, Conc):
        return Conc * CytoVol * NA

    def GetConc(self, Count):
        return Count / CytoVol / NA

    def AdjustRefdCon(self, Reaction, RefMol, RefdConc):
        ''' Compare dConc of reference molecule to input concentrations and adjust reference dConc '''
        RefCoeff = Reaction.Stoich[RefMol]
        UnitdConc = RefdConc / RefCoeff

        Out = 10   # 1 Molar to begin comparison
        for Mol, Coeff in Reaction.Input.items():
            AdjustedConc = self.Molecules[Mol] / Coeff

            Out = min(Out, min(AdjustedConc, UnitdConc))

            # DEBUG
            # if self.Molecules[Mol] < 1e-8:
            #     print(Mol, self.Molecules[Mol], AdjustedConc, UnitdConc, Out)

        return Out * RefCoeff

    def DeterminedConc(self, Reaction, RefMol, RefdConc):
        UnitdConc = RefdConc / Reaction.Stoich[RefMol]
        dConc = dict()
        for Mol, Coeff in Reaction.Stoich.items():
            dConc[Mol] = UnitdConc * Coeff

            # DEBUG
            # print(Mol, dConc[Mol])

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
            assert dConc + self.Molecules[dMolecule] >= 0.0, 'Iter {}\t | {} \t| Conc:{}, \t dConc:{}'.format(self.Iter, dMolecule, self.Conc2Str(self.Molecules[dMolecule]), self.Conc2Str(dConc))
            self.Molecules[dMolecule] += dConc

    def RestorePermanentMolecules(self):
        for Molecule in self.PermanentMolecules:
            if Molecule in self.Molecules:
                self.Molecules[Molecule] = self.InitialConditions[Molecule]

    def CheckZeroConcentrations(self):
        for Molecule, Conc in self.Molecules.items():
            if Conc < self.GetConc(1): # if
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

    def Conc2Str(self, Conc):
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

    def Info(self):
        for Molecule, Conc in self.Molecules.items():
            ConStr = self.Conc2Str(Conc)
            print("{:<16}: {:>10}".format(Molecule, ConStr))

    def Summary(self):
        print("{:<16} {:>10} {:>10} {:>10} {:<10}".format("", "Initial", "Final", "dConc", "(Permanent)"))
        for Molecule, Conc in self.Molecules.items():
            dConc = Conc - self.InitialConditions[Molecule]

            InitConc = self.Conc2Str(self.InitialConditions[Molecule])
            FinalConc = self.Conc2Str(Conc)
            FinaldConc = self.Conc2Str(dConc)

            print("{:16} {:>10} {:>10} {:>10} {}".format(Molecule, InitConc, FinalConc, FinaldConc, "*" if Molecule in self.PermanentMolecules else ""))
        if self.Plot:
            print("\nPlot: ON")
        else:
            print("\nPlot: OFF")

    def PrintMolConc(self, MolConcDict, End='\t'):
        for Molecule, Conc in MolConcDict.items():
            ConStr = self.Conc2Str(Conc)
            print("\t{}: {}".format(Molecule, ConStr), end=End)
        print()

    def Simulate(self, TotalTime = 10, DeltaTime = 0.01,):
        print("-- Initial Conditions --")
        self.Info()
        print()

        while self.Iter < TotalTime / DeltaTime:
            # # Debug
            # if self.Iter == 1942:
            #     print("")

            for Reaction in self.Reactions:
                dMolecules = self.RunReaction(Reaction, DeltaTime)
                if self.Debug_Reaction:
                    print("[{}] {}".format(Reaction.ReactionName, Reaction.DisplayChemicalEquation()))
                    self.PrintMolConc(dMolecules, End=', ')
                self.UpdateMolecules(dMolecules)
                self.CheckZeroConcentrations()

            if len(self.PermanentMolecules):
                self.RestorePermanentMolecules()

            if self.Plot:
                self.AddToDataset()

            self.Iter += 1

            if self.Iter % self.Debug_Info == 0:
                print()
                print("-- Iteration {} --".format(self.Iter))
                self.Info()
                print()

        print("\n")
        print("-- Summary --")
        self.Summary()

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
                KnownMolConc[Metabolite] = float(ConcInEcoli)
            return KnownMolConc

        db_KnownMolConc = self.LoadTSVDatabase("MetaboliteConcentrations.tsv")
        KnownMolConc = ParseMetaboliteCon(db_KnownMolConc)

        KnownMolConc["FADH2"] = (KnownMolConc["NADH"] / KnownMolConc["NAD+"]) * KnownMolConc["FAD"]

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

        elif ReactionName == "OxidativePhosphorylation":
            Sim.AddReaction(OxidativePhosphorylation())
            Sim.AddReaction(ATPControl(-4.35E-04))
            Sim.AddReaction(NADHControl(1.0874e-4))
            Sim.AddReaction(FADH2Control(1.0875e-4))

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
    Sim = Simulator()

    RunUnitTest = False
    # RunUnitTest = True

    if RunUnitTest:
        UnitTestReactions = "OxidativePhosphorylation"
        Sim.UnitTest(UnitTestReactions)

    else:
        # Add Reactions
        # Sim.AddReaction(NADPlusSynthesis())
        # Sim.AddReaction(NADPPlusSynthesis())
        # Sim.AddReaction(CoASynthesis())
        Sim.AddReaction(Glycolysis())
        Sim.AddReaction(PyruvateOxidation())
        Sim.AddReaction(TCACycle_alphaKetoGlutarateSynthesis())
        Sim.AddReaction(TCACycle_SuccinylCoASynthesis())
        Sim.AddReaction(TCACycle_FumarateSynthesis())
        Sim.AddReaction(TCACycle_OxaloacetateSynthesis())
        # Sim.AddReaction(PyruvateCarboxylation())
        Sim.AddReaction(OxidativePhosphorylation())
        Sim.AddReaction(ATPControl())
        # Sim.AddReaction(ATPControl())
        # Sim.AddReaction(ATPControl())

        # Set permanent molecules
        PermanentMolecules = [
            "G6P",
            # "pyruvate",
            "CoA-SH",
            "NADH",
            "NAD+",
            "FADH2",
            "FAD",
            # "CoA-SH",
            # "oxaloacetate",
        ]
        Sim.SetPermanentMolecules(PermanentMolecules)

    # Debugging options
    # Sim.Debug_Reaction = True
    Sim.Debug_Info = 100
    Sim.Plot = True

    # Simulation parameters
    TotalTime = 100
    DeltaTime = 0.01

    # Execute simulation
    Sim.Initialize()
    Sim.Simulate(TotalTime=TotalTime, DeltaTime=DeltaTime)

    # Plot
    if Sim.Plot:
        Datasets = Sim.GetDataset()
        Plot = FPlotter()
        Plot.PlotDatasets(Datasets, DeltaTime=DeltaTime, Multiscale=True)

"""
Reference

Metabolite concentrations, fluxes and free energies imply efficient enzyme usage
Junyoung O Park, Sara A Rubin, Yi-Fan Xu, Daniel Amador-Noguez, Jing Fan, Tomer Shlomi & Joshua D Rabinowitz
Nature Chemical Biology volume 12, pages482–489 (2016)
https://www.nature.com/articles/nchembio.2077#Sec21
"""
