'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Product Inhibition
Network  4.1,   p.78
Model   4.1,    p.77
Figure  4.2,    p.78

'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt

def dS1_NumericalSimulation(k1, k3, k5, K, n, S1, S2):
    return (k1 / (1 + (S2 / K) ** n)) - (k3 * S1) - (k5 * S1)

def dS2_NumericalSimulation(k2, k4, k5, S1, S2):
    return k2 + k5 * S1 - k4 * S2

class FModel():
    def __init__(self, InS1=0, InS2=0, Ink1=20, Ink2=5, Ink3=5, Ink4=5, Ink5=2, InK=1, Inn=4):
        # Initial Concentrations
        self.S1 = InS1
        self.S2 = InS2

        # Constants
        self.k1 = Ink1
        self.k2 = Ink2
        self.k3 = Ink3
        self.k4 = Ink4
        self.k5 = Ink5
        self.K = InK
        self.n = Inn

        # Dataset
        self.Data_S1 = list()
        self.Data_S2 = list()

        # Set initial values
        self.InitializeSimStepZero()

    def InitializeSimStepZero(self):
        self.AppendData(self.S1, self.S2)

    def AppendData(self, S1, S2):
        self.Data_S1.append(S1)
        self.Data_S2.append(S2)

    def Run(self, SimSteps=500, TimeResolution=100):
        Flat = 0.00000001 # Steady state threshold

        i = 0
        while i < SimSteps:
            # Calculate new
            S1 = self.Data_S1[-1]
            S2 = self.Data_S2[-1]

            dS1 = dS1_NumericalSimulation(self.k1, self.k3, self.k5, self.K, self.n, S1, S2) / TimeResolution
            dS2 = dS2_NumericalSimulation(self.k2, self.k4, self.k5, S1, S2) / TimeResolution

            self.AppendData(self.Data_S1[-1] + dS1, self.Data_S2[-1] + dS2)

            if (abs(self.Data_S1[-1] - self.Data_S1[-2]) < Flat):
                break
            i += 1

    def PlotData(self):
        X = self.Data_S1
        Y = self.Data_S2

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)

        # Dynamics
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)

        ax1.plot(X, 'r-', label="[S1]")
        ax1.plot(Y, 'b-', label="[S2]")
        ax1.set_title('Product Inhibition: Dynamics')
        ax1.set_xlabel('Time(a.u.)')
        ax1.set_ylabel('Concentration(a.u.)')
        ax1.legend(loc='upper right')
        ax1.grid()

        # Phase plane
        ax2.plot(X, Y, color="blue")
        ax2.set_title('Product Inhibition: Phase plane')
        ax2.set_xlabel("[S1]")
        ax2.set_ylabel("[S2]")
        ax2.grid()

        plt.show()

def main():
    Model = FModel()
    Model.Run()
    Model.PlotData()

if __name__ == '__main__':

    main()