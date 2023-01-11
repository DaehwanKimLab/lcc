# BSD 3-Clause License
# © 2022, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import sys
from datetime import datetime
import numpy as np
import math
import matplotlib.pyplot as plt


class Reaction:
    def __init__(self):
        self.ReactionName = ""
        self.Input = dict()
        self.Output = dict()
        self.Stoich = dict()

    def Specification(self, Molecules):
        return {}

# reaction NAD+Synthesis(nicotinamide + PRPP + 2 ATP -> NAD+)
class NADPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"nicotinamide": 1, "PRPP": 1, "ATP": 2}
        self.Output = {"NAD+": 1}

    def Specification(self, Molecules):
        return {}

# reaction NADP+Synthesis(NAD+ + ATP -> NADP+ + ADP)
class NADPPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"NAD+": 1, "ATP": 1}
        self.Output = {"NADP+": 1, "ADP": 1}

    def Specification(self, Molecules):
        return {}

# reaction CoASynthesis(pantothenate + cysteine + 3 ATP + CTP + CO2 -> CoA + ADP + CMP + 3 PPi)
class CoASynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"pantothenate": 1, "cysteine": 1, "ATP": 3, "CTP": 1, "CO2": 1}
        self.Output = {"CoA": 1, "ADP": 1, "CMP": 1, "PPi": 3}

    def Specification(self, Molecules):
        return {}

"""
reaction Glycolysis(G6P + 2 ADP + 2 NAD+ -> 2 pyruvate + 2 NADH + 2 ATP)
specification
{
    float cg = 0.019M;
    ATP'  = (ADP / ATP) * cg;
    ADP'  = - ATP';
    NADH' = ATP';
    NAD' = - NADH';
    G6P'  = - ATP' / 2;
    pyruvate' = G6P' * 2;
}

G6P[:] = 8.8mM;
ATP = 9.6mM;
ADP = 0.56mM;
NADH = 83uM;
NAD = 2.6mM;
pyruvate = 0.39mM;
"""
class Glycolysis(Reaction):
    def __init__(self):
        super().__init__()
        self.ReactionName = 'Glycolysis'
        self.Input = {"G6P": 1, "ADP": 2, "NAD+": 2}
        self.Output = {"pyruvate": 1, "NADH": 2, "ATP": 2}

    def Specification(self, Molecules):
        c = 0.019
        ATP = Molecules["ATP"]
        ADP = Molecules["ADP"]
        dATP = (ADP / ATP) * c
        # dATP = (abs(ADP - 0.00056) / ATP) * c   # For homeostasis of ADP level
        # dATP = (ADP * abs(ATP - 0.0096) / ATP) * c   # For homeostasis of ATP level
        # dATP = ((ADP / ATP) - (0.00056 / 0.0096)) * ATP * c   # For homeostasis of ADP to ATP ratio
        return "ATP", dATP

class ATPConsumption(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"ATP": 1}
        self.Output = {"ADP": 1}

    def Specification(self, Molecules):
        # Cell Division ATP consumption: 4.35E-03
        c = 4.35E-03
        # c = 4.35E-05
        dATP = -c
        return "ATP", dATP

# For Debugging
class ADPConsumption(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"ADP": 1}
        self.Output = {"ATP": 1}

    def Specification(self, Molecules):
        c = 4.35E-03
        dADP = -c
        return "ADP", dADP


class Simulator():
    def __init__(self):
        self.Reactions = []
        self.Dataset = {}

        self.PermanentMolecules = {
            "G6P"      : 8.8  * 1e-3,
            }

        self.Molecules = dict()

        self.KnownMolConc = {
            "G6P": 8.8 * 1e-3,
            "ATP": 9.6 * 1e-3,
            "ADP": 0.56 * 1e-3,
            "NADH": 83 * 1e-6,
            "NAD+": 2.6 * 1e-3,
            "pyruvate": 0.39 * 1e-3,
        }

        # Metabolite concentrations, fluxes and free energies imply efficient enzyme usage
        # Junyoung O Park, Sara A Rubin, Yi-Fan Xu, Daniel Amador-Noguez, Jing Fan, Tomer Shlomi & Joshua D Rabinowitz
        # Nature Chemical Biology volume 12, pages482–489 (2016)
        # https://www.nature.com/articles/nchembio.2077#Sec21

        self.InitialConditions = {}

        self.Plot = False

    def Initialize(self):
        self.InitializeStoich()
        self.InitializeMolecules()
        if self.Plot:
            self.InitializeDataset()

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
            return 0.0

    def GetReactionNames(self):
        ReactionNames = []
        for Reaction in self.Reactions:
            ReactionNames.append(Reaction.ReactionName)
        return str(ReactionNames)

    def InitializeMolecules(self):
        for Reaction in self.Reactions:
            for Molecule, Coeff in Reaction.Stoich.items():
                if Molecule not in self.Molecules:
                    self.Molecules[Molecule] = self.GetKnownMolConc(Molecule)
        self.InitialConditions = self.Molecules.copy()

    def AddReaction(self, Reaction):
        self.Reactions.append(Reaction)

    def AdjustRefdCon(self, Reaction, RefMol, RefdConc):
        ''' Compare dConc of reference molecule to input concentrations and adjust reference dConc '''
        RefCoeff = Reaction.Stoich[RefMol]
        UnitdConc = RefdConc / RefCoeff

        Out = 0
        for Mol, Coeff in Reaction.Input.items():
            AdjustedConc = self.Molecules[Mol] / Coeff
            Out = min(AdjustedConc, UnitdConc)

        return Out * RefCoeff

    def DeterminedConc(self, Reaction, RefMol, RefdConc):
        UnitdConc = RefdConc / Reaction.Stoich[RefMol]
        dConc = dict()
        for Mol, Coeff in Reaction.Stoich.items():
            dConc[Mol] = UnitdConc * Coeff
        return dConc

    def RunReaction(self, Reaction):
        RefMol, RefdConc = Reaction.Specification(self.Molecules)
        RefdConc = self.AdjustRefdCon(Reaction, RefMol, RefdConc)
        return self.DeterminedConc(Reaction, RefMol, RefdConc)

    def UpdateMolecules(self, dMolecules, DeltaTime):
        for dMolecule, dConc in dMolecules.items():
            assert dMolecule in self.Molecules
            dConc = dConc * DeltaTime            
            assert dConc + self.Molecules[dMolecule] >= 0.0, '{} \t| Conc:{}, \t dConc:{}'.format(dMolecule, self.Conc2Str(self.Molecules[dMolecule]), self.Conc2Str(dConc))
            self.Molecules[dMolecule] += dConc

    def Adjust(self):
        for Molecule, Conc in self.PermanentMolecules.items():
            assert Molecule in self.Molecules
            self.Molecules[Molecule] = Conc

    def InitializeDataset(self):
        for Molecule, Conc in self.Molecules.items():
            self.Dataset[Molecule] = [Conc]

    def AddToDataset(self):
        for Molecule, Conc in self.Molecules.items():
            self.Dataset[Molecule].append(Conc)

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

    def Simulate(self, TotalTime = 10, DeltaTime = 0.01,):
        print("-- Initial Conditions --")
        self.Info()

        Iter = 0
        while Iter < TotalTime / DeltaTime:
            # # Debug
            # if Iter == 220:
            #     print("\n")
            #     print("Debugging:", Iter)

            for Reaction in self.Reactions:
                dMolecules = self.RunReaction(Reaction)
                self.UpdateMolecules(dMolecules, DeltaTime)

            self.Adjust()

            if self.Plot:
                self.AddToDataset()

            Iter += 1
            print("\n")
            print("-- Iteration {} --".format(Iter))

            self.Info()

        print("\n")
        print("-- Summary --")
        self.Summary()

        if Sim.Plot:
            Datasets = {self.GetReactionNames(): self.Dataset}

            Plot = FPlotter()
            Plot.PlotDatasets(Datasets, DeltaTime=DeltaTime)

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

    def PlotDatasets(self, Datasets, DeltaTime=1.0, YMin=0, YMax=0, bSideLabel=True, SuperTitle=""):
        # Filter Datasets
        if self.Filter_Inclusion or self.Filter_Exclusion:
            Datasets = self.FilterDatasets(Datasets)

        def ExtractTime(Dataset, SimulationTimeUnit):
            Time = None
            for Data in Dataset.values():
                Time = [i * SimulationTimeUnit for i in range(len(Data))]
                break
            return Time, Dataset

        def ApplyUnit(Dataset):
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

        def GetSubPlottingInfo(Datasets, MaxNPlotsInRows=3):
            NPlotsInRows = len(Datasets)  # Default
            if len(Datasets) > 1:
                for Remainder in range(MaxNPlotsInRows):
                    if len(Datasets) % (Remainder + 1) == 0:
                        NPlotsInRows = Remainder + 1
            return NPlotsInRows

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.5, hspace=0.5)
        if SuperTitle:
            fig.suptitle(SuperTitle, fontsize=14)
        NPlotsInRows = GetSubPlottingInfo(Datasets)

        # global ymax
        YMax = YMax
        for Process, Dataset in Datasets.items():
            Dataset_Copy = Dataset.copy()
            UnitTxt, Dataset_Copy = ApplyUnit(Dataset_Copy)
            # Y axis (molecular concentrations)
            for MolName, Conc in Dataset_Copy.items():
                YMaxData = max(Conc) * 1.1
                if YMaxData > YMax:
                    YMax = YMaxData

        # Plot data
        for n, (Process, Dataset) in enumerate(Datasets.items()):
            Time, Dataset = ExtractTime(Dataset, DeltaTime)
            UnitTxt, Dataset = ApplyUnit(Dataset)
            ax1 = fig.add_subplot(math.ceil(len(Datasets) / NPlotsInRows), NPlotsInRows, n + 1)

            # Y axis (molecular concentrations)
            PerturbationIndex = 0
            for MolName, Conc in Dataset.items():

                # YMax Debug
                YMaxData = max(Conc) * 1.1
                print('Plot', MolName, YMaxData)

                line, = ax1.plot(Time, Conc, label="[" + MolName + "]")
                if bSideLabel:
                    SelectedTimeFrameFromLeft = 0.9
                    ax1.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="left", va="bottom", color=line.get_color())
                    # ax1.text(Time[-1] * 1.01, Conc[-1], MolName, ha="left", va="bottom", color=line.get_color())
                    # ax1.text(Time[-1] * 1.1, Conc[-1], MolName + ": {}".format(Conc[-1]), va="center", color=line.get_color())

            ax1.set_title(Process)
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Molecules: Concentration (' + UnitTxt + ')')
            ax1.set_ylim(ymax=YMax)
            if YMin:
                ax1.set_ylim(ymin=YMin)
            else:
                ax1.set_ylim(ymin=0)

            if not bSideLabel:
                ax1.legend(loc='upper left')

            # ax1.grid()
        plt.show()


if __name__ == '__main__':
    Sim = Simulator()

    # Add reactions
    # Sim.AddReaction(NADPlusSynthesis())
    # Sim.AddReaction(NADPPlusSynthesis())
    # Sim.AddReaction(CoASynthesis())
    Sim.AddReaction(Glycolysis())
    # Sim.AddReaction(ATPConsumption())
    # Sim.AddReaction(ADPConsumption())

    Sim.Plot = True
    TotalTime = 10
    DeltaTime = 0.01

    Sim.Initialize()
    Sim.Simulate(TotalTime=TotalTime, DeltaTime=DeltaTime)
