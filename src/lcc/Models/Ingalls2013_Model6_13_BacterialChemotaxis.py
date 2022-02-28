'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Two Component Signaling
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

nM = 1e-9
Mol2Count = 6.0221409e+23
Count2Mol = 1.0 / Mol2Count

def dAm_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return min(krev1 * R, A) - (k1 * BP) * Am / (kM1 + Am) - k3 * Am * L + krev3 * AmL

def dAmL_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return min(krev2 * R, AL) - (k2 * BP) * AmL / (kM2 + AmL) + k3 * Am * L - krev3 * AmL
    
def dA_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return -min(krev1 * R, A) + (k1 * BP) * Am / (kM1 + Am) - k4 * A * L + krev4 * AL
    
def dAL_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return -min(krev2 * R, AL) + (k2 * BP) * AmL / (kM2 + AmL) + k4 * A * L - krev4 * AL

def dB_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return -k5 * Am * B + krev5 * BP

def dBP_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return k5 * Am * B - krev5 * BP

class FModel():
    def __init__(self,
                 k1=200, k2=1, k3=1/nM, k4=1/nM, k5=0.05/nM, kM1=1*nM, kM2=1*nM,
                 L=20*nM, R=5*nM, A=500*nM, B=0.1*nM,
                 krev1=1, krev2=1, krev3=1, krev4=1, krev5=0.005,
                 Am=0, AL=0, AmL=0, BP=0):

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

    def InitializeSimStepZero(self):
        self.AppendData(self.L, self.R, self.Am, self.AmL, self.A, self.AL, self.B, self.BP)

    def AppendData(self, L, R, Am, AmL, A, AL, B, BP):
        self.Data_L.append(L)
        self.Data_R.append(R)
        self.Data_Am.append(Am)
        self.Data_AmL.append(AmL)
        self.Data_A.append(A)
        self.Data_AL.append(AL)
        self.Data_B.append(B)
        self.Data_BP.append(BP)

    def Run(self, SimSteps=100000, TimeResolution=100):
        Flat = 0.000001 # Steady state threshold

        self.TimeResolution = TimeResolution

        i = 0
        while i < SimSteps:
            # Calculate new

            Time = i / self.TimeResolution
            if Time in [100, 300]:
                self.L = self.L * 2

            L = self.L
            self.Simulate(L)

            # DK - debugging purposes
            Am = self.Data_Am[-1]
            print("Time: {}, Steps: {}, Am: {}".format(Time, i, Am))

            i += 1

    def Simulate(self, L, SimUnitTime = None):
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
        X = self.Data_L
        A = self.Data_A
        Am = self.Data_Am
        AL = self.Data_AL
        AmL = self.Data_AmL
        B = self.Data_B
        BP = self.Data_BP

        T = list(np.array(range(len(X))) / self.TimeResolution)

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)

        # PlotBegin = 5000
        PlotBegin = 0

        # Dynamics
        ax1 = fig.add_subplot(2, 2, 1)
        ax1.plot(T[PlotBegin:], Am[PlotBegin:], 'b-', label="[Am]") # RL
        ax1.set_title('Bacterial Chemotaxis')
        ax1.set_xlabel('Time (a.u.)')
        ax1.set_ylabel('Concentration (a.u.)')
        ax1.legend(loc='upper right')
        ax1.grid()

        ax2 = fig.add_subplot(2, 2, 2)
        ax2.plot(T[PlotBegin:], Am[PlotBegin:], 'b-', label="[Am]") # RL
        ax2.plot(T[PlotBegin:], A[PlotBegin:], 'r-', label="[A]") # RL
        ax2.set_title('Bacterial Chemotaxis')
        ax2.set_xlabel('Time (a.u.)')
        ax2.set_ylabel('Concentration (a.u.)')
        ax2.legend(loc='upper right')
        ax2.grid()

        ax3 = fig.add_subplot(2, 2, 3)
        ax3.plot(T[PlotBegin:], AmL[PlotBegin:], 'b-', label="[AmL]") # RL
        ax3.plot(T[PlotBegin:], AL[PlotBegin:], 'r-', label="[AL]") # RL
        ax3.set_title('Bacterial Chemotaxis')
        ax3.set_xlabel('Time (a.u.)')
        ax3.set_ylabel('Concentration (a.u.)')
        ax3.legend(loc='upper right')
        ax3.grid()

        ax4 = fig.add_subplot(2, 2, 4)
        ax4.plot(T[PlotBegin:], BP[PlotBegin:], 'r-', label="[BP]") # RL
        ax4.set_title('Bacterial Chemotaxis')
        ax4.set_xlabel('Time (a.u.)')
        ax4.set_ylabel('Concentration (a.u.)')
        ax4.legend(loc='upper right')
        ax4.grid()

        plt.show()

def main():
    Model = FModel()
    Model.Run()
    Model.PlotData()

if __name__ == '__main__':

    main()
