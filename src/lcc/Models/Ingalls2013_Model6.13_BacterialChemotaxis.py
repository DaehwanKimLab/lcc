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

def dAm_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return krev1 * R * A - (k1 * BP * Am) / (kM1 + Am) - k3 * Am * L + krev3 * AmL

def dAmL_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return krev2 * R * A - (k2 * BP * AmL) / (kM2 + AmL) + k3 * Am * L - krev3 * AmL

def dA_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return -krev1 * R * A + (k1 * BP * Am) / (kM1 + Am) - k4 * A * L + krev4 * AL

def dAL_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return -krev2 * R * A + (k2 * BP * AmL) / (kM2 + AmL) + k4 * A * L - krev4 * AL

def dB_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return -k5 * Am * B + krev5 * BP

def dBP_NumericalSimulation(k1, k2, k3, k4, k5, krev1, krev2, krev3, krev4, krev5, kM1, kM2, L, R, Am, AmL, A, AL, B, BP):
    return k5 * Am * B - krev5 * BP


class FModel():
    def __init__(self, k1=200, k2=1, k3=1, k4=1, k5=0.05, krev1=1, krev2=1, krev3=1, krev4=1, krev5=0.005, kM1=1, kM2=1, L=20, R=5, Am=0, AmL=0, A=1, AL=0, B=0.1, BP=0):
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

        self.TimeResolution = 1

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

    def Run(self, SimSteps=10000, TimeResolution=100):
        Flat = 0.000001 # Steady state threshold

        self.TimeResolution = TimeResolution

        i = 0
        while i < SimSteps:
            # Calculate new

            Time = i / self.TimeResolution
            if ((Time == 10) or (Time == 30)):
                self.L = self.L * 2

            L = self.L
            R = self.R
            Am = self.Data_Am[-1]
            AmL= self.Data_AmL[-1]
            A = self.Data_A[-1]
            AL = self.Data_AL[-1]
            B = self.Data_B[-1]
            BP = self.Data_BP[-1]

            dAm = dAm_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3, self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B, BP) / TimeResolution
            dAmL = dAmL_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3, self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B, BP) / TimeResolution
            dA = dA_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3, self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B, BP) / TimeResolution
            dAL = dAL_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3, self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B, BP) / TimeResolution
            dB = dB_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3, self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B, BP) / TimeResolution
            dBP = dBP_NumericalSimulation(self.k1, self.k2, self.k3, self.k4, self.k5, self.krev1, self.krev2, self.krev3, self.krev4, self.krev5, self.kM1, self.kM2, L, R, Am, AmL, A, AL, B, BP) / TimeResolution

            self.AppendData(L, R, self.Data_Am[-1] + dAm, self.Data_AmL[-1] + dAmL, self.Data_A[-1] + dA, self.Data_AL[-1] + dAL, self.Data_B[-1] + dB, self.Data_BP[-1] + dBP)

            i += 1

    def PlotData(self):
        X = self.Data_L
        Y = self.Data_Am

        T = list(np.array(range(len(X))) / self.TimeResolution)

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)

        # Dynamics
        ax1 = fig.add_subplot(1, 1, 1)

        # ax1.plot(T, X, 'r-', label="[Ligand]") # L
        ax1.plot(T, Y, 'b-', label="[Active CheA]") # RL
        ax1.set_title('Bacterial Chemotaxis')
        ax1.set_xlabel('Time (a.u.)')
        ax1.set_ylabel('Concentration (a.u.)')
        ax1.legend(loc='upper right')
        ax1.grid()

        plt.show()

def main():
    Model = FModel()
    Model.Run()
    Model.PlotData()

if __name__ == '__main__':

    main()