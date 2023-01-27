# BSD 3-Clause License
# © 2023, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import sys
from datetime import datetime
import copy, csv, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class FPlotter:
    def __init__(self):
        self.Filter_Inclusion = None
        self.Filter_Exclusion = None
        self.KnownMolConc = dict()
        self.ConsistentColorDict = dict()

        self.SubplotID = 1
        self.XLabel = None
        self.YLabel = None

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

    def ResetSubplotID(self):
        self.SubplotID = 1

    def PlotDatasets(self, Datasets, DeltaTime=1.0, bSideLabel='right', SuperTitle="", MaxNPlotsInRows=3, Unitless=True, All=True, Multiscale=False, Include_nM=False, Individual=False, MolRange=False, Export=''):
        # Filter Datasets
        if self.Filter_Inclusion or self.Filter_Exclusion:
            Datasets = self.FilterDatasets(Datasets)

        def ExtractTime(Datasets, SimulationTimeUnit):
            for Process, Dataset in Datasets.items():
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

        def GetSubPlottingInfo(Datasets, MaxNPlotsInRows):
            ScaleFactor = 0
            if All:
                ScaleFactor += 1
            if Multiscale:
                ScaleFactor += 3 if Include_nM else 2
            if Individual:
                for Dataset in Datasets.values():
                    ScaleFactor += len(Dataset)

            NPlotsInRows = ScaleFactor   # Default
            if ScaleFactor > MaxNPlotsInRows:
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

        def PlotDataset(Dataset, Process, UnitTxt, Plt):
            # Y axis (molecular concentrations)
            YMax = 0
            for MolName, Conc in Dataset.items():

                line, = Plt.plot(Time, Conc, label="[" + MolName + "]")

                # Save assigned color for each molecule
                if MolName not in self.ConsistentColorDict:
                    self.ConsistentColorDict[MolName] = line.get_color()
                else:
                    line._color = self.ConsistentColorDict[MolName]

                # Display Molecule Range
                if MolRange:
                    Plt.fill_between(Time, self.KnownMolConc[MolName][1], self.KnownMolConc[MolName][2], alpha=0.2, facecolor=self.ConsistentColorDict[MolName])

                # Mol labeling on the curve
                if bSideLabel == 'right' or bSideLabel == 'both':
                    SelectedTimeFrameFromLeft = 0.99
                    Plt.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="right", va="bottom", color=self.ConsistentColorDict[MolName])
                if bSideLabel == 'left' or bSideLabel == 'both':
                    SelectedTimeFrameFromLeft = 0
                    Plt.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="left", va="bottom", color=self.ConsistentColorDict[MolName])

                YMaxData = max(Conc) * 1.1
                if YMaxData > YMax:
                    if Multiscale:
                        YMax = min(1000 * 1.1, YMaxData)
                    else:
                        YMax = YMaxData

            Plt.set_title(Process)
            if self.XLabel:
                Plt.set_xlabel(self.XLabel)
            else:
                Plt.set_xlabel('Time (s)')
            if self.YLabel:
                Plt.set_ylabel(self.YLabel)
            else:
                Plt.set_ylabel('Conc (' + UnitTxt + ')')
            Plt.set_ylim(ymin=0)
            Plt.set_ylim(ymax=YMax)
            # Plt.grid()

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
        ScaleFactor, NPlotsInRows, NPlotsInColumn = GetSubPlottingInfo(Datasets, MaxNPlotsInRows)

        # Plot data
        def PlotAll():
            for Process, Dataset in copy.deepcopy(Datasets).items():
                UnitTxt = ''
                if not Unitless:
                    UnitTxt, Dataset = ApplyGlobalUnit(Dataset)
                Plt = fig.add_subplot(NPlotsInColumn, NPlotsInRows, self.SubplotID)
                PlotDataset(Dataset, Process, UnitTxt, Plt)
                self.SubplotID += 1

        def PlotMultiscale():
            assert not Unitless, "To plot multiscale, set unitless option False"

            def SplitScales(Dataset):
                Dataset_mM = dict()
                Dataset_uM = dict()
                Dataset_nM = dict()
                for mol, conc in Dataset.items():
                    if conc[0] < 1e-6:
                        if Include_nM:
                            Dataset_nM[mol] = (np.array(conc) / 1e-9).tolist()
                            AdjustUnitOnKnownConc(mol, 1e-9)
                    elif conc[0] < 1e-3:
                        Dataset_uM[mol] = (np.array(conc) / 1e-6).tolist()
                        AdjustUnitOnKnownConc(mol, 1e-6)
                    else:
                        Dataset_mM[mol] = (np.array(conc) / 1e-3).tolist()
                        AdjustUnitOnKnownConc(mol, 1e-3)
                return (Dataset_mM, Dataset_uM, Dataset_nM)

            for Process, Dataset in copy.deepcopy(Datasets).items():
                Dataset_Unit = SplitScales(Dataset)
                Units = ['mM', 'uM', 'nM']

                for Unit, dataset in zip(Units, Dataset_Unit):
                    if not Include_nM and Unit == 'nM':
                        continue
                    Plt = fig.add_subplot(NPlotsInColumn, NPlotsInRows, self.SubplotID)
                    PlotDataset(dataset, Unit + ' range', Unit, Plt)
                    self.SubplotID += 1

        def PlotIndividual():

            def ConvertDatasetToIndividualDatasets(Dataset):
                AdjustedDatasets = dict()
                Units = dict()
                for mol, conc in Dataset.items():
                    if Unitless:
                        AdjustedDatasets[mol] = {mol: conc}
                        Units[mol] = ''
                    elif conc[0] < 1e-6:
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

            for Process, Dataset in copy.deepcopy(Datasets).items():
                Unit, AdjustedDatasets = ConvertDatasetToIndividualDatasets(Dataset)
                for MolName, Dataset in AdjustedDatasets.items():
                    Plt = fig.add_subplot(NPlotsInColumn, NPlotsInRows, self.SubplotID)
                    PlotDataset(Dataset, MolName, Unit[MolName], Plt)
                    self.SubplotID += 1

        self.ResetSubplotID()

        if All:
            PlotAll()

        if Multiscale:
            PlotMultiscale()

        if Individual:
            PlotIndividual()

        if Export == 'pdf':
            PrintToPDF()

        plt.show()


class FGrowthPlotter(FPlotter):
    def __init__(self):
        super().__init__()
        self.XLabel = 'Time (s)'
        self.YLabel = 'Number of E. coli'


class FGeneRepressionPlotter():
    def PlotDatasets(self, Datasets):
        for Title, Dataset in Datasets.items():
            XSet, YSet = Dataset
            plt.title(Title)
            plt.ylabel("Relative Growth")
            plt.xlabel("Repression Efficiency")
            plt.plot(XSet, YSet, "b")
            plt.plot(XSet, YSet, "go")
            plt.show()

