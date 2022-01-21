'''
Ingalls 2013 Mathematical Modelling in Systems Biology

G-protein Signaling
Network  6.4,   p.154
Model    6.4,   p.154
Figure   6.5,   p.154

'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dRL_NumericalSimulation(kRL, kRLm, kGa, kGd0, kG1, R, L, RL, G, Ga, Gbg, Gd):
    return kRL * R * L - kRLm * RL

def dG_NumericalSimulation(kRL, kRLm, kGa, kGd0, kG1, R, L, RL, G, Ga, Gbg, Gd):
    return -kGa * G * RL + kG1 * Gd * Gbg

def dGa_NumericalSimulation(kRL, kRLm, kGa, kGd0, kG1, R, L, RL, G, Ga, Gbg, Gd):
    return kGa * G * RL - kGd0 * Ga

def dGbg_NumericalSimulation(kRL, kRLm, kGa, kGd0, kG1, R, L, RL, G, Ga, Gbg, Gd):
    return kGa * G * RL - kG1 * Gd * Gbg

def dGd_NumericalSimulation(kRL, kRLm, kGa, kGd0, kG1, R, L, RL, G, Ga, Gbg, Gd):
    return kGd0 * Ga - kG1 * Gd * Gbg


"""
kRL = 2×10−3nM−1 s−1, kRLm =10−2 s−1, kGa = 10−5 (molecules per cell)−1 s−1, kGd0 = 0.004 s−1, kG1 = 1 (molecules per cell)−1 s−1. 
The total G-protein population is 10000 molecules per cell, while the total receptor population is 4000 molecules per cell.
"""

class FModel():
    def __init__(self, InkRL=2e6, InkRLm=1e-2, InkGa=1e-4, InkGd0=0.004e-9, InkG1=1e-9, InTotal_R=4000, InTotal_G=10000):
        self.Vol = 1

        # Initial Concentrations
        self.Total_R = InTotal_R
        self.Total_G = InTotal_G

        # Constants
        self.kRL = InkRL
        self.kRLm = InkRLm
        self.kGa = InkGa
        self.kGd0 = InkGd0
        self.kG1 = InkG1

        # Dataset
        self.Data_L = list()
        self.Data_RL = list()
        self.Data_G = list()
        self.Data_Ga = list()
        self.Data_Gbg = list()
        self.Data_Gd = list()

        # Set initial values
        self.InitializeSimStepZero()

    def InitializeSimStepZero(self):
        self.AppendData(0, 0, self.Total_G, 0, 0, 0)

    def AppendData(self, L, RL, G, Ga, Gbg, Gd):
        self.Data_L.append(L)
        self.Data_RL.append(RL)
        self.Data_G.append(G)
        self.Data_Ga.append(Ga)
        self.Data_Gbg.append(Gbg)
        self.Data_Gd.append(Gd)

    def Run(self, SimSteps=1500, TimeResolution=1):
        Flat = 0.000001 # Steady state threshold

        i = 0
        while i < SimSteps:
            # Calculate new
            L = None
            Time = i / TimeResolution
            if ((Time >= 100) & (Time <700)):
                L = 1e-9   # Concentration
            else:
                L = 0

            # Get Concentrations
            RL = self.Data_RL[-1] / self.Vol
            R = (self.Total_R - self.Data_RL[-1]) / self.Vol
            G = self.Data_G[-1] / self.Vol
            Ga = self.Data_Ga[-1] / self.Vol
            Gbg = self.Data_Gbg[-1] / self.Vol
            Gd = self.Data_Gd[-1] / self.Vol

            if Time > 100:
                print(RL, R, G, Ga, Gbg, Gd)

            # Calculate delta Count (not concentration)
            dRL = int(dRL_NumericalSimulation(self.kRL, self.kRLm, self.kGa, self.kGd0, self.kG1, R, L, RL, G, Ga, Gbg, Gd) * self.Vol / TimeResolution)
            dG = int(dG_NumericalSimulation(self.kRL, self.kRLm, self.kGa, self.kGd0, self.kG1, R, L, RL, G, Ga, Gbg, Gd) * self.Vol / TimeResolution)
            dGa = int(dGa_NumericalSimulation(self.kRL, self.kRLm, self.kGa, self.kGd0, self.kG1, R, L, RL, G, Ga, Gbg, Gd) * self.Vol / TimeResolution)
            dGbg = int(dGbg_NumericalSimulation(self.kRL, self.kRLm, self.kGa, self.kGd0, self.kG1, R, L, RL, G, Ga, Gbg, Gd) * self.Vol / TimeResolution)
            dGd = int(dGd_NumericalSimulation(self.kRL, self.kRLm, self.kGa, self.kGd0, self.kG1, R, L, RL, G, Ga, Gbg, Gd) * self.Vol / TimeResolution)

            self.AppendData(L, self.Data_RL[-1] + dRL, self.Data_G[-1] + dG, self.Data_Ga[-1] + dGa, self.Data_Gbg[-1] + dGbg, self.Data_Gd[-1] + dGd)

            i += 1

    def PlotData(self):
        # X = self.Data_L
        Y = self.Data_RL
        Z = self.Data_Ga

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)

        # Dynamics
        ax1 = fig.add_subplot(1, 1, 1)

        # ax1.plot(X, 'r-', label="[Ligand]") # L
        ax1.plot(Y, 'b-', label="[Ligand-bound Receptor]") # RL
        ax1.plot(Z, 'g-', label="[G_alpha]")
        ax1.set_title('G-protein Signaling')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Count')
        ax1.legend(loc='upper left')
        ax1.grid()

        plt.show()

def main():
    Model = FModel()
    Model.Run()
    Model.PlotData()

if __name__ == '__main__':

    main()