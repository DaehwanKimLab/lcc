'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Bacterial Chemotaxis
Network  6.13,   p.163
Model   6.12,    p.162
Figure  6.14,    p.164

'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

"""
Figure 6.14: Behavior of the chemotaxis signal transduction pathway. Starting at a low level of chemoattractant
ligand ([L] = 20), the system responds to a doubling of ligand at time t = 10 with an immediate
drop in CheA activity (corresponding to a reduction in tumbling) followed by a return to almost the original
nominal activity level. A second doubling of ligand concentration at time t = 30 produces a similar effect.
Parameter values are, (in time−1): k1 = 200, k2 = 1, k3 = 1, k4 = 1, k5 = 0.05, k−1 = 1, k−2 = 1, k−3 = 1,
k−4 = 1, k−5 = 0.005; (in concentration): kM1 = 1, kM2 = 1. Units are arbitrary.
"""

"""
https://www.math.uwaterloo.ca/~bingalls/MMSB/Code/matlab/chemotaxis.m
cheR=1;
L=20;
%assign initial condition species vector S=[Am AmL A AL B BP]
%these values were determined by running a previous simulation to steady
%state with L=20.
S0=[0.0360    1.5593    0.0595    0.3504    0.7356    0.2644];
"""

uM = 1e-6
nM = 1e-9

Unit = nM

Mol2Count = 6.0221409e+23
Count2Mol = 1.0 / Mol2Count


def dAm_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B,
                            BP):
    return krev1 * R - (k1 * BP) * Am / (kM1 + Am) - k3 * Am * L + krev3 * AmL
    # return min(krev1 * R, A) - (k1 * BP) * Am / (kM1 + Am) - k3 * Am * L + krev3 * AmL


def dAmL_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B,
                             BP):
    return krev2 * R - (k2 * BP) * AmL / (kM2 + AmL) + k3 * Am * L - krev3 * AmL
    # return min(krev2 * R, AL) - (k2 * BP) * AmL / (kM2 + AmL) + k3 * Am * L - krev3 * AmL


def dA_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B,
                           BP):
    return -krev1 * R + (k1 * BP) * Am / (kM1 + Am) - k4 * A * L + krev4 * AL
    # return -min(krev1 * R, A) + (k1 * BP) * Am / (kM1 + Am) - k4 * A * L + krev4 * AL


def dAL_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B,
                            BP):
    return -krev2 * R + (k2 * BP) * AmL / (kM2 + AmL) + k4 * A * L - krev4 * AL
    # return -min(krev2 * R, AL) + (k2 * BP) * AmL / (kM2 + AmL) + k4 * A * L - krev4 * AL


def dB_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B,
                           BP):
    return -k5 * Am * B + krev5 * BP


def dBP_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B,
                            BP):
    return k5 * Am * B - krev5 * BP

class FModel():
    def __init__(self,
                 k1=200, k2=1.0, k3=1 / Unit, k4=1 / Unit, k5=0.05 / Unit, kM1=1 * Unit, kM2=1 * Unit,  # for Ingalls model
                 L=20 * Unit, R=5 * Unit, A=500 * Unit, B=0.1 * Unit,  # for Ingalls model
                 krev1=1, krev2=1, krev3=1, krev4=1, krev5=0.005,
                 Am=0, AL=0, AmL=0, BP=0):

        R = 1 * Unit

        Am=0.0360 * Unit
        AmL=1.5593 * Unit
        A=0.0595 * Unit
        AL=0.3504 * Unit
        B=0.7356 * Unit
        BP=0.2644 * Unit

        # Initial Concentrations
        self.L = L
        self.R = R
        self.Am = Am
        self.AmL = AmL
        self.A = A
        self.AL = AL
        self.B = B
        self.BP = BP

        # Constants
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4
        self.k5 = k5
        self.krev1 = krev1
        self.krev2 = krev2
        self.krev3 = krev3
        self.krev4 = krev4
        self.krev5 = krev5
        self.kM1 = kM1
        self.kM2 = kM2

        # Dataset
        self.Data_L = list()
        self.Data_R = list()
        self.Data_Am = list()
        self.Data_AmL = list()
        self.Data_A = list()
        self.Data_AL = list()
        self.Data_B = list()
        self.Data_BP = list()

        # Set initial values
        self.InitializeSimStepZero()

        self.TimeResolution = 100

        # Debugging purposes
        self.Time_LigandInduction = list()
        self.Amount_LigandInduction = list()
        self.Homeostasis = list()
        self.HomeostasisCheckInterval = 100

        # DoublePlot
        self.DoublePlotSwitch = False
        self.ShowTrueUnits = False

    def InitializeSimStepZero(self):
        self.AppendData(self.L, self.R, self.Am, self.AmL, self.A, self.AL, self.B, self.BP)

    def SetLigandInitialConcentration(self, Input):
        self.L = Input

    def AppendData(self, L, R, Am, AmL, A, AL, B, BP):
        self.Data_L.append(L)
        self.Data_R.append(R)
        self.Data_Am.append(Am)
        self.Data_AmL.append(AmL)
        self.Data_A.append(A)
        self.Data_AL.append(AL)
        self.Data_B.append(B)
        self.Data_BP.append(BP)

    def AddToHomeostasis(self, Value):
        if Value not in self.Homeostasis:
            self.Homeostasis.append(Value)

    def PrintHomeostasis(self):
        print("# Homeostasis points: (Time, [Am]")
        print(self.Homeostasis)

    def Run(self, SimSteps=1000000, TimeResolution=100):
        Flat = 0.000001  # Steady state threshold

        self.TimeResolution = TimeResolution

        i = 0

        # Debugging purposes
        self.Time_LigandInduction = []
        if self.Time_LigandInduction:
            ManualLigandInduction = True
        else:
            ManualLigandInduction = False

        N_LigandInduction = 2
        LigandFoldChange = 2
        # TimeToChangeLigand = [range(100, 300)]
        # LigandFoldChange = 1.001

        # Testing different homeostasis levels by ligand level change
        # self.SetLigandInitialConcentration(0.1 * nM)
        # self.SetLigandInitialConcentration(1 * nM)
        # self.SetLigandInitialConcentration(10 * nM)

        while i < SimSteps:
            # for homeostasis recording
            PrevAm = self.Data_Am[-1]

            # Calculate new

            Time = i / self.TimeResolution
            if Time in self.Time_LigandInduction:
                if ManualLigandInduction:
                    self.L = self.L * LigandFoldChange
                else:
                    index = self.Time_LigandInduction.index(Time)
                    # DK - debugging purposes
                    # self.L = self.L * 1.001
                    self.L = self.Amount_LigandInduction[index]

            L = self.L
            self.Simulate(L)

            # DK - debugging purposes
            Am = self.Data_Am[-1]
            print("Time: {}, Steps: {}, Am: {}".format(Time, i, Am))

            # for homeostasis recording
            if not ManualLigandInduction:
                if Time % self.HomeostasisCheckInterval == 0:
                    if PrevAm > 0 and abs(Am - PrevAm) / PrevAm < 1e-9:
                        if len(self.Homeostasis) < N_LigandInduction:
                            Time_LigandInduction = Time + 50
                            Amount_LigandInduction = self.L * LigandFoldChange
                            self.Time_LigandInduction.append(Time + 50)
                            self.Amount_LigandInduction.append(self.L * 2)
                        else:
                            SimSteps = i
                        self.AddToHomeostasis((Time, Am))

            i += 1

        self.PrintHomeostasis()

    def Simulate(self, L, SimUnitTime=None):
        R = self.R
        Am = self.Data_Am[-1]
        AmL = self.Data_AmL[-1]
        A = self.Data_A[-1]
        AL = self.Data_AL[-1]
        B = self.Data_B[-1]
        BP = self.Data_BP[-1]

        if SimUnitTime == None:
            SimUnitTime = 1.0 / self.TimeResolution

        dAm = dAm_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3,
                                      self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B,
                                      BP) * SimUnitTime
        dAmL = dAmL_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3,
                                        self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B,
                                        BP) * SimUnitTime
        dA = dA_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3,
                                    self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B,
                                    BP) * SimUnitTime
        dAL = dAL_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3,
                                      self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B,
                                      BP) * SimUnitTime
        dB = dB_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3,
                                    self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B,
                                    BP) * SimUnitTime
        dBP = dBP_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3,
                                      self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B,
                                      BP) * SimUnitTime

        self.AppendData(L, R, Am + dAm, AmL + dAmL, A + dA, AL + dAL, B + dB, BP + dBP)

        return Am

    def PlotData(self):
        L = self.Data_L
        A = self.Data_A
        Am = self.Data_Am
        AL = self.Data_AL
        AmL = self.Data_AmL
        B = self.Data_B
        BP = self.Data_BP

        UnitTxt = '(a.u.)'

        if self.ShowTrueUnits:
            Am = [x / Unit for x in Am]
            L = [x / Unit for x in L]
            if Unit == uM:
                UnitTxt = '(uM)'
            elif Unit == nM:
                UnitTxt = '(nM)'

        T = list(np.array(range(len(L))) / self.TimeResolution)

        # Plotting options
        if self.DoublePlotSwitch:
            self.DoublePlot(T, Am, L, UnitTxt)

        # Debugging purposes
        else:
            fig = plt.figure()
            fig.subplots_adjust(wspace=0.2, hspace=0.3)

            PlotBegin = 0

            # Dynamics
            ax1 = fig.add_subplot(2, 2, 1)
            ax1.plot(T[PlotBegin:], Am[PlotBegin:], 'b-', label="[Am]")  # RL
            ax1.set_title('Bacterial Chemotaxis')
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Concentration ' + UnitTxt)
            ax1.set_yticks([0, 50, 0.01, 0.04])
            ax1.legend(loc='upper right')
            ax1.grid()

            ax2 = fig.add_subplot(2, 2, 2)
            ax2.plot(T[PlotBegin:], Am[PlotBegin:], 'b-', label="[Am]")  # RL
            ax2.plot(T[PlotBegin:], A[PlotBegin:], 'r-', label="[A]")  # RL
            ax2.set_ylabel('Concentration ' + UnitTxt)
            ax2.legend(loc='upper right')
            ax2.grid()

            ax3 = fig.add_subplot(2, 2, 3)
            ax3.plot(T[PlotBegin:], AmL[PlotBegin:], 'b-', label="[AmL]")  # RL
            ax3.plot(T[PlotBegin:], AL[PlotBegin:], 'r-', label="[AL]")  # RL
            ax3.set_xlabel('Time (a.u.)')
            ax3.set_ylabel('Concentration ' + UnitTxt)
            ax3.legend(loc='upper right')
            ax3.grid()

            ax4 = fig.add_subplot(2, 2, 4)
            ax4.plot(T[PlotBegin:], BP[PlotBegin:], 'r-', label="[BP]")  # RL
            ax4.set_xlabel('Time (a.u.)')
            ax4.set_ylabel('Concentration ' + UnitTxt)
            ax4.legend(loc='upper right')
            ax4.grid()

            plt.show()

    def DoublePlot(self, T, Am, L, UnitTxt):
        fig = plt.figure()

        PlotBegin = int(self.Homeostasis[0][0] - self.HomeostasisCheckInterval) * self.TimeResolution

        # Dynamics
        ax1 = plt.axes()
        ax1.set_ylim(min(Am[PlotBegin:]) * 0.2, max(Am[PlotBegin:]) * 1.3)
        ax1.plot(T[PlotBegin:], Am[PlotBegin:], 'r-', label="[Am]")  # RL
        ax1.set_title('Bacterial Chemotaxis')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('AM Concentration ' + UnitTxt)
        ax1.legend(loc='upper right')
        ax1.grid()

        ax1T = ax1.twinx()
        ax1T.set_ylim(min(L[PlotBegin:]) * 0.5, max(L[PlotBegin:]) * 4)
        ax1T.plot(T[PlotBegin:], L[PlotBegin:], 'b-', label="[Glucose]")  # RL
        ax1T.set_ylabel('Glucose Concentration ' + UnitTxt)
        ax1T.legend(loc='lower right')

        plt.show()

def main():
    Model = FModel()
    Model.Run()
    Model.PlotData()

if __name__ == '__main__':
    main()