import os, sys
from argparse import ArgumentParser
import csv
import matplotlib.pyplot as plt
import math
import numpy as np
import random

random.seed(1)
SaveFilename = None
NA = 6.0221409e+23
PerturbationTag = "#"

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


    def PlotDatasets(self, Datasets, SimulationTimeUnit=1, bSideLabel=True, SuperTitle=""):
        # Filter Datasets
        if self.Filter_Inclusion or self.Filter_Exclusion:
            Datasets = self.FilterDatasets(Datasets)

        def ExtractTime(Dataset, SimulationTimeUnit):
            Time = None
            if 'Time' in Dataset:
                Time = Dataset['Time']
                del Dataset['Time']
            else:
                for Data in Dataset.values():
                    Time = [i * SimulationTimeUnit for i in range(len(Data))]
                    break
            return Time, Dataset

        def GetPerturbation(Dataset):
            Perturbation = 0
            PerturbationPlotColor = list()
            for Key in Dataset.keys():
                if Key[0] == PerturbationTag:
                    PerturbationPlotColor.append((random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 0.75)))
                    Perturbation += 1
            return Perturbation, PerturbationPlotColor

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
                if DatasetArray.shape[0] == 0:
                    DatasetArray = np.array(Data, ndmin=2)
                    continue
                DatasetKeyIndex[Key] = i
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

        # Plot data
        for n, (Process, Dataset) in enumerate(Datasets.items()):
            Time, Dataset = ExtractTime(Dataset, SimulationTimeUnit)
            Perturbation, PerturbationPlotColor = GetPerturbation(Dataset)
            UnitTxt, Dataset = ApplyUnit(Dataset)
            ax1 = fig.add_subplot(math.ceil(len(Datasets) / NPlotsInRows), NPlotsInRows, n + 1)
            ax2 = None
            if Perturbation:
                ax2 = ax1.twinx()

            # Y axis (molecular concentrations)
            PerturbationIndex = 0
            for MolName, Conc in Dataset.items():
                if MolName[0] != PerturbationTag:

                    # Debug
                    # print(Process, "\t", MolName)

                    line, = ax1.plot(Time, Conc, label="[" + MolName + "]")
                    if bSideLabel:
                        SelectedTimeFrameFromLeft = 0.1
                        if Perturbation == 0:
                            SelectedTimeFrameFromLeft = 0.9
                        ax1.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="left", va="bottom", color=line.get_color())
                        # ax1.text(Time[-1] * 1.01, Conc[-1], MolName, ha="left", va="bottom", color=line.get_color())
                        # ax1.text(Time[-1] * 1.1, Conc[-1], MolName + ": {}".format(Conc[-1]), va="center", color=line.get_color())

                else:
                    line, = ax2.plot(Time, Conc, color=PerturbationPlotColor[PerturbationIndex], label="[" + MolName[1:] + "]")
                    if bSideLabel:
                        SelectedTimeFrameFromLeft = 0.8
                        ax2.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName[1:], ha="center", va="bottom", color=line.get_color())
                    # ax2.plot(Time, Conc, label="[" + MolName[2:] + "]")
                    # print(PerturbationIndex)
                    PerturbationIndex += 1

            ax1.set_title(Process)
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Molecules: Concentration (' + UnitTxt + ')')
            # ax1.set_ylim([0, 0.015])
            if not bSideLabel:
                ax1.legend(loc='upper left')
            if Perturbation:
                ax2.set_ylabel('Events: Concentration (' + UnitTxt + ')')
                # ax2.set_ylim([0, 3e-7])
                if not bSideLabel:
                    ax2.legend(loc='upper right')

            # ax1.grid()

        if SaveFilename:
            plt.savefig(SaveFilename)
        else:
            plt.show()

    def LoadRawData(self, Data_Dir):
        Datasets = dict()

        def Parse_tsv(FilePath, FileName):
            Fullpath = FilePath + '/' + FileName
            # print(fname)

            def RemoveEmptyRows(ListOfRows):
                RowsWithInfo = list()
                for Row in ListOfRows:
                    if Row != list():
                        RowsWithInfo.append(Row)
                return RowsWithInfo

            with open(Fullpath) as fp:
                Reader = csv.reader(fp, delimiter='\t')
                Dataset = RemoveEmptyRows(list(Reader))

                Names = Dataset[0]
                Data = np.array(Dataset[1:]).transpose()

                Dataset = dict()
                for i, Name in enumerate(Names):
                    Vol = 1   # Temporary
                    if Data[i] == list():   # empty data
                        continue
                    elif Name in ['SimStep', 'SimTime', 'Time', 'Time (s)']:
                        Dataset['Time'] = [float(x) for x in Data[i]]
                        continue
                    elif Name == 'Vol':
                        continue
                        # Dataset[Name] = [float(x) for x in Data[i]]   # TODO: Volume is temporarily blocked off
                    elif Name == 'Pseudo':  # Pseudo Molecule from old lcc
                        continue
                    else:   # TODO: Modify according to the final lcc data format and take options?
                        # For Absolute count
                        # Dataset[Name] = [float(x) for x in Data[i]]

                        # From absolute count to concentrations
                        Dataset[Name] = [float(x) / NA / Vol for x in Data[i]]

                Datasets[FileName] = Dataset

        for FileName in os.listdir(Data_Dir):
            if FileName.endswith('.tsv'):
                Parse_tsv(Data_Dir, FileName)

        return Datasets

def main(SaveToFile = None):
    global SaveFilename
    if SaveToFile is not None:
        SaveFilename = SaveToFile

    Plotter = FPlotter()

    Data_Dir = '.'
    Datasets = Plotter.LoadRawData(Data_Dir)
    # for FileName, Dataset in Datasets.items():
    #     PlotData(Dataset)
    Plotter.PlotDatasets(Datasets)

if __name__ == '__main__':
    parser = ArgumentParser(description = 'LCC plot')
    parser.add_argument('--save-fig',
            dest='save_fname',
            type=str,
            help='Save figure to file')

    args = parser.parse_args()

    main(args.save_fname)

