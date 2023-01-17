# BSD 3-Clause License
# © 2022, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import sys
from datetime import datetime
import random, copy
import numpy as np
import math
import matplotlib.pyplot as plt


class Simulator():
    def __init__(self):
        None

    def Initialize(self):
        None

    def Simulate(self, TotalTime = 10, DeltaTime = 0.01):
        None

    def GetDataset(self):
        return None


class EcoliSimulator(Simulator):
    DNAREPLICATIONRATE = 1.0 / (16 * 60)
    PROTEINSYNTHESISRATE = 1.0 / (16 * 60)
    CYTOKINESISRATE = 1.0 / ( 6 * 60)

    def __init__(self, 
                 DNAReplicationRate = DNAREPLICATIONRATE,
                 ProteinSynthesisRate = PROTEINSYNTHESISRATE,
                 CytoKinesisRate = CYTOKINESISRATE):
        self.DNAReplicationRate = DNAReplicationRate
        self.DNAReplicationProgress = 0.0
        self.ProteinSynthesisRate = ProteinSynthesisRate
        self.ProteinSynthesisProgress = 0.0
        self.CytoKinesisRate = CytoKinesisRate
        self.CytoKinesisProgress = 0.0

        self.NumCellDivisions = 0

    def SimulateDelta(self, DeltaTime = 1.0):
        if self.DNAReplicationProgress >= 1.0 and \
           self.ProteinSynthesisProgress >= 1.0:
            self.CytoKinesisProgress += self.CytoKinesisRate * DeltaTime
        else:
            if self.DNAReplicationProgress < 1.0:
                self.DNAReplicationProgress += self.DNAReplicationRate * DeltaTime
            if self.ProteinSynthesisProgress < 1.0:
                self.ProteinSynthesisProgress += self.ProteinSynthesisRate * DeltaTime

        self.NumCellDivisions = 0
        if self.CytoKinesisProgress >= 1.0:
            self.DNAReplicationProgress = 0.0
            self.ProteinSynthesisProgress = 0.0
            self.CytoKinesisProgress = 0.0
            self.NumCellDivisions = 1 

    def SetProgress(self, Progress):
        self.DNAReplicationProgress = Progress
        self.ProteinSynthesisProgress = Progress
        self.CytoKinesisProgress = 0.0

    def GetNumCellDivisions(self):
        return self.NumCellDivisions

    
class ColonySimulator(Simulator):
    def __init__(self, EcoliSim):
        self.EcoliSims = [[1, EcoliSim]]
        self.DataPoints = []

    def SimulateDelta(self, DeltaTime = 1.0):
        for i, (NumEcoli, Sim) in enumerate(self.EcoliSims):
            Sim.SimulateDelta(DeltaTime)
            NumCellDivisions = Sim.GetNumCellDivisions()
            if NumCellDivisions > 0:
                NumEcoli *= math.pow(2, NumCellDivisions)
            self.EcoliSims[i][0] = NumEcoli
            if len(self.EcoliSims) < 20:
                if NumEcoli > 1:
                    NumEcoli1 = int(NumEcoli / 2)
                    NumEcoli2 = NumEcoli - NumEcoli1
                    self.EcoliSims[i][0] = NumEcoli1
                    Sim2 = copy.deepcopy(Sim)
                    Sim2.SetProgress(random.uniform(0, 0.5))
                    self.EcoliSims.append([NumEcoli2, Sim2])
                    

        self.AddToDataset()

    def Info(self):
        None

    def AddToDataset(self):
        TotalNumEcoli = 0
        for NumEcoli, _ in self.EcoliSims:
            TotalNumEcoli += NumEcoli
        self.DataPoints.append(TotalNumEcoli)
        
    def GetDataset(self):
        return self.DataPoints

    def GetNumEcoli(self):
        return self.DataPoints[-1] if len(self.DataPoints) > 0 else 0


class PopulationSimulator(Simulator):
    def __init__(self):
        self.Colonies = [
            ColonySimulator(EcoliSimulator(EcoliSimulator.DNAREPLICATIONRATE,
                                           EcoliSimulator.PROTEINSYNTHESISRATE,
                                           EcoliSimulator.CYTOKINESISRATE)),
            ColonySimulator(EcoliSimulator(EcoliSimulator.DNAREPLICATIONRATE,
                                           EcoliSimulator.PROTEINSYNTHESISRATE * 0.5,
                                           EcoliSimulator.CYTOKINESISRATE)),
            ColonySimulator(EcoliSimulator(EcoliSimulator.DNAREPLICATIONRATE,
                                           EcoliSimulator.PROTEINSYNTHESISRATE,
                                           EcoliSimulator.CYTOKINESISRATE * 0.5))
        ]
        self.Dataset = {}

    def Simulate(self, TotalTime = 24 * 60.0 * 60.0, DeltaTime = 1.0):
        Iter = 0
        while Iter < TotalTime / DeltaTime:
            for ColonySim in self.Colonies:
                ColonySim.SimulateDelta(DeltaTime)

            self.AddToDataset()

            Iter += 1
            if Iter % 100 == 0:
                print("\n")
                print("-- Iteration {} --".format(Iter))
                self.Info()

    def Info(self):
        for i, Colony in enumerate(self.Colonies):
            print("Colony {:>2}: {:>10}".format(i, Colony.GetNumEcoli()))

    def AddToDataset(self):
        for i, Colony in enumerate(self.Colonies):
            if i not in self.Dataset:
                self.Dataset[i] = []
            self.Dataset[i].append(Colony.GetNumEcoli())
        
    def GetDataset(self):
        PlotDataset = {}
        TotalData = None
        for i, Data in self.Dataset.items():
            Name = "Colony {:>2}".format(i)
            PlotDataset[Name] = Data
            if TotalData:
                TotalData = [a + b for a, b in zip(TotalData, Data)]
            else:
                TotalData = Data
        PlotDataset["Population"] = TotalData        

        return {"E. coli Growth": PlotDataset}


class PetridishSimulator(Simulator):
    def __init__(self):
        self.PopSim = PopulationSimulator()

    def Simulate(self, TotalTime = 10, DeltaTime = 0.01):
        self.PopSim.Simulate(TotalTime, DeltaTime)

    def GetDataset(self):
        return self.PopSim.GetDataset()

    
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
                'mM': 1e-3,
                'count' : 1,
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

                line, = ax1.plot(Time, Conc, label="[" + MolName + "]")
                if bSideLabel:
                    SelectedTimeFrameFromLeft = 0.9
                    ax1.text(Time[-1] * SelectedTimeFrameFromLeft, Conc[int(len(Time) * SelectedTimeFrameFromLeft)] * 1.02, MolName, ha="left", va="bottom", color=line.get_color())

            ax1.set_title(Process)
            ax1.set_xlabel('Time (s)')
            # ax1.set_ylabel('Molecules: Concentration (' + UnitTxt + ')')
            ax1.set_ylabel('Number of E. coli')
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
    Sim = PetridishSimulator()
    Sim.Plot = True

    TotalTime = 6 * 60.0 * 60.0
    DeltaTime = 1.0

    Sim.Initialize()
    Sim.Simulate(TotalTime=TotalTime, DeltaTime=DeltaTime)

    if Sim.Plot:
        Datasets = Sim.GetDataset()
        Plot = FPlotter()
        Plot.PlotDatasets(Datasets, DeltaTime=DeltaTime)
