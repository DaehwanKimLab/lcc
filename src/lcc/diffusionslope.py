import random
import sys
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np

# from src.lcc.plot2 import FPlotter

class FDiffusion():
    def __init__(self):
        # Simulation Parameters
        self.SimStep = 0
        self.TotalSimulationTime = 0
        self.TimeResolution = 0
        self.DistanceResolution = 0

        # Model Parameters
        self.Title = ""
        self.Type = ""
        self.m = 0   # slope
        self.b = 0   # y-intercept
        self.TotalDistance = 0
        self.TotalQuantity = 0
        self.DiffusionRate = 0

        # Data
        self.Dataset = dict()
        self.IntegralDataset = dict()

    def SetSimulationParameters(self, TotalSimulationTime, TimeResolution):
        self.TotalSimulationTime = TotalSimulationTime
        self.TimeResolution = TimeResolution

    def PrintSimulationParameters(self):
        print("\n-------------------- Simulation Parameters --------------------")
        print("Total Simulation Time :", self.TotalSimulationTime, "s", " | Time Resolution :", self.TimeResolution, "s")
        print("\n-------------------- Run Simulation --------------------")

    def PrintSimulationSubject(self):
        print("\n\nSimulating " + self.Title + "...")

    def InitializeSimStep(self):
        self.SimStep = 0

    def InitializeModelParameters(self):
        self.Title = ""
        self.Type = ""
        self.m = 0
        self.b = 0
        self.TotalDistance = 0
        self.DistanceResolution = 0
        self.TotalQuantity = 0
        self.DiffusionRate = 0

    def PrintDataset(self):
        print(self.Dataset)

    def PrintSimTime(self):
        print("\n", self.SimStep, "\t| ", "{:.2f}".format(self.SimStep * self.TimeResolution), end=" s | ")

    def PrintSimTimeAndDataset(self):
        self.PrintSimTime()
        self.PrintDataset()

    def SetTitleAndType(self, Title, Type):
        self.Title = Title
        self.Type = Type

    def SetDiffusionParameters(self, TotalDistance, DistanceResolution, TotalQuantity, DiffusionRate):
        self.TotalDistance = TotalDistance
        self.DistanceResolution = DistanceResolution
        self.TotalQuantity = TotalQuantity
        self.DiffusionRate = DiffusionRate

    def SetFunctionParameters(self, m, b):
        self.m = m
        self.b = b

    def GetFunctionParameters_Old(self, Time):
        CompletionFactor = Time / self.TotalSimulationTime
        m = (CompletionFactor / 2) - 0.5
        b = self.TotalQuantity / self.TotalDistance * (1 - 0.5 * CompletionFactor)
        return (m, b)

    def GetFunctionParameters(self, Time, DiffusionRate=2, ClosedSystem=True, Debug=False):
        m = - self.TotalQuantity
        # CompletionFactor = Time / self.TotalSimulationTime # 0.1
        # FinalIntercept = self.TotalQuantity / self.TotalDistance
        if ClosedSystem:
            m += DiffusionRate * Time * 20

            if Debug:
                print(m, - m * 2 * self.TotalQuantity)

            if m >= 0:
                if Debug:
                    print("reverted slope")
                return (0, self.TotalQuantity / self.TotalDistance)
            b = math.sqrt(- m * 2 * self.TotalQuantity)

            def TestIntegral(m, b):
                # Integral Test and correction
                Data = self.Simulation_LinearFunction(m, b)
                Integral = IntegralApproximation(Data, 0, self.TotalDistance)
                Integrity = Integral / self.TotalQuantity

                if Debug:
                    print("Integrity: ", Integrity)

                if Integrity < 0.999:
                    b += b * 0.01
                    TestIntegral(m, b)
                return b

            b = TestIntegral(m, b)
            print("Time: ", Time, " s\tm: ", m, "\tb: ", b)
            return (m, b)

        else:
            print("Open system is not supported")
            sys.exit(1)

    def GetDistanceRange(self):
        return np.arange(0, self.TotalDistance, self.DistanceResolution)

    def Simulation_LinearFunction(self, m, b):

        def LinearFunction(x, m, b):
            return m * x + b

        DistanceRange = self.GetDistanceRange()
        Data = np.zeros(DistanceRange.shape[0])
        for i in range(DistanceRange.shape[0]):
            Output = LinearFunction(DistanceRange[i], m, b)
            if Output < 0:
                break
            else:
                Data[i] = Output

        return Data.tolist()

    def Initialize(self):
        self.InitializeModelParameters()
        self.InitializeSimStep()

    def Run(self, message=False):
        self.PrintSimulationParameters()
        while self.SimStep * self.TimeResolution <= self.TotalSimulationTime:
            Time = self.SimStep * self.TimeResolution
            Slope, Intercept = self.GetFunctionParameters(Time, self.DiffusionRate)
            self.SetFunctionParameters(Slope, Intercept)

            # Generate Data
            Data = None
            if self.Type == 'Linear':
                Data = self.Simulation_LinearFunction(self.m, self.b)
            else:
                print("no type was specified")
            self.Dataset[str(Time) + 's'] = Data

            # Debugging
            if message:
                self.PrintSimTimeAndDataset()

            self.SimStep += 1

    def IntegralCheck(self, message=False):

        for Title, Data in self.Dataset.items():
            Integral = IntegralApproximation(Data, 0, self.TotalDistance)
            self.IntegralDataset[Title] = Integral

        if message:
            self.PrintIntegralDataset()

    def PrintIntegralDataset(self):
        print("\n-------------------- Integral Validation --------------------")
        for Title, Data in self.IntegralDataset.items():
            print(Title, "Integral:", Data)

    def GetDataset(self):
        return self.Dataset.copy()

def IntegralApproximation(Data, Start, End):
    return (End - Start) * np.mean(Data)

# Modified from plot2.py
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

    def PlotDatasets(self, Datasets, XRange, SuperTitle=""):
        print("\n-------------------- Plotting... --------------------")

        # Filter Datasets
        if self.Filter_Inclusion or self.Filter_Exclusion:
            Datasets = self.FilterDatasets(Datasets)

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
            ax1 = fig.add_subplot(math.ceil(len(Datasets) / NPlotsInRows), NPlotsInRows, n + 1)

            for TimeStamp, Conc in Dataset.items():
                # Debug
                # print(Process, "\t", MolName)

                line, = ax1.plot(XRange, Conc, label="[" + TimeStamp + "]")
                Integral = IntegralApproximation(Conc, 0, XRange[-1])
                SelectedTimeFrameFromLeft = 0.01
                Label = TimeStamp + ", Integral: " + str(Integral)
                ax1.text(XRange[-1] * SelectedTimeFrameFromLeft, Conc[int(len(XRange) * SelectedTimeFrameFromLeft)] * 1.02, Label, ha="left", va="bottom", color=line.get_color())

            ax1.set_title(Process)
            ax1.set_xlabel('Distance (a.u.)')
            # ax1.set_ylim([0, 0.015])
            ax1.set_ylabel('Concentration (a.u.)')
            # ax1.legend(loc='upper left')
            ax1.grid()

        if SaveFilename:
            plt.savefig(SaveFilename)
        else:
            plt.show()

def main():
    Diffusion = FDiffusion()
    Plot = FPlotter()

    # Setting Simulation Parameters
    TotalSimulationTime = 30   # s
    TimeResolution = 0.5
    Diffusion.SetSimulationParameters(TotalSimulationTime, TimeResolution)

    # Model Parameters
    Type = 'Linear'
    TotalDistance = 10
    DistanceResolution = 0.0001
    TotalQuantity = 100
    DiffusionRate = 0.2
    Diffusion.SetDiffusionParameters(TotalDistance, DistanceResolution, TotalQuantity, DiffusionRate)
    Diffusion.SetTitleAndType(Type, Type)

    # Run model
    Diffusion.Run(message=False)

    # Integral check
    Diffusion.IntegralCheck(message=True)

    # Plot datasets
    Datasets = dict()
    Datasets[Type] = Diffusion.GetDataset()
    DistanceRange = Diffusion.GetDistanceRange()
    Plot.PlotDatasets(Datasets, DistanceRange)

if __name__ == '__main__':
    main()
