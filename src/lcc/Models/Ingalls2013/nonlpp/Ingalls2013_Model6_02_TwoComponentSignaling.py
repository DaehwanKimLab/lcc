'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Two Component Signaling
Network  6.2,   p.151
Model   6.2,    p.151
Figure  6.3,    p.152

'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dR_NumericalSimulation(k1, krev1, k2, k3, R, L, RL, P, pP):
    return -k1 * R * L + krev1 * RL

def dRL_NumericalSimulation(k1, krev1, k2, k3, R, L, RL, P, pP):
    return k1 * R * L - krev1 * RL

def dP_NumericalSimulation(k1, krev1, k2, k3, R, L, RL, P, pP):
    return -k2 * P * RL + k3 * pP

def dpP_NumericalSimulation(k1, krev1, k2, k3, R, L, RL, P, pP):
    return k2 * P * RL - k3 * pP

class FModel():
    def __init__(self, Ink1=5, Inkrev1=1, Ink2=6, Ink3=3, InRL=0, InTotal_R=2, InpP=0, InTotal_P=8):
        # Initial Concentrations
        self.RL = InRL
        self.pP = InpP

        self.Total_R = InTotal_R
        self.Total_P = InTotal_P

        # Constants
        self.k1 = Ink1
        self.krev1 = Inkrev1
        self.k2 = Ink2
        self.k3 = Ink3

        # Dataset
        self.Data_L = list()
        self.Data_RL = list()
        self.Data_pP = list()

        self.Data_Time = list()

        # Set initial values
        self.InitializeSimStepZero()

    def InitializeSimStepZero(self):
        self.AppendData(0, self.RL, self.pP)

    def AppendData(self, L, RL, pP):
        self.Data_L.append(L)
        self.Data_RL.append(RL)
        self.Data_pP.append(pP)

    def Run(self, SimSteps=1000, TimeResolution=100):
        Flat = 0.000001 # Steady state threshold

        i = 0
        while i < SimSteps:
            # Calculate new
            L = None
            Time = i / TimeResolution
            if ((Time >= 1) & (Time <3)):
                L = 3
            else:
                L = 0

            RL = self.Data_RL[-1]
            R = self.Total_R - RL
            pP = self.Data_pP[-1]
            P = self.Total_P - pP

            dRL = dRL_NumericalSimulation(self.k1, self.krev1, self.k2, self.k3, R, L, RL, P, pP) / TimeResolution
            dpP = dpP_NumericalSimulation(self.k1, self.krev1, self.k2, self.k3, R, L, RL, P, pP) / TimeResolution

            self.AppendData(L, self.Data_RL[-1] + dRL, self.Data_pP[-1] + dpP)
            self.Data_Time.append(i/TimeResolution)

            i += 1

    def PlotData(self):
        X = self.Data_L
        Y = self.Data_RL
        Z = self.Data_pP

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)

        # Dynamics
        ax1 = fig.add_subplot(1, 1, 1)

        ax1.plot(X, 'r-', label="[Ligand]") # L
        ax1.plot(Y, 'b-', label="[Ligand-bound Receptor]") # RL
        ax1.plot(Z, 'g-', label="[Phosphorylated Target Protein]")
        ax1.set_title('Two Component Signaling')
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