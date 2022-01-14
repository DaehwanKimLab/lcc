'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Limit Cycle Oscillations
Network  4.14,   p.93
Model   4.10,    p.93
Figure  4.15,    p.94

'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt

def dS1_NumericalSimulation(k0, k1, K, n, S1, S2):
    return k0 - (k1 * (1 + (S2 / K) ** n) * S1)

def dS2_NumericalSimulation(k1, k2, K, n, S1, S2):
    return (k1 * (1 + (S2 / K) ** n) * S1) - (k2 * S2)

class FModel():
    def __init__(self, InS1=1.5, InS2=1, Ink0=8, Ink1=1, Ink2=5, InK=1, Inn=2):   # Inn=2 is too strong
    # def __init__(self, InS1=1.5, InS2=1, Ink0=8, Ink1=1, Ink2=11, InK=1, Inn=1):   # Inn=2 is too strong
        # Initial Concentrations
        self.S1 = InS1
        self.S2 = InS2

        # Constants
        self.k0 = Ink0
        self.k1 = Ink1
        self.k2 = Ink2
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

    def Run(self, SimSteps=5000, TimeResolution=100):
        Flat = 0.000001 # Steady state threshold

        i = 0
        while i < SimSteps:
            # Calculate new
            S1 = self.Data_S1[-1]
            S2 = self.Data_S2[-1]

            dS1 = dS1_NumericalSimulation(self.k0, self.k1, self.K, self.n, S1, S2) / TimeResolution
            dS2 = dS2_NumericalSimulation(self.k1, self.k2, self.K, self.n, S1, S2) / TimeResolution

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
        ax1.set_title('Dynamics')
        ax1.set_xlabel('Time(a.u.)')
        ax1.set_ylabel('Concentration(a.u.)')
        ax1.legend(loc='upper left')
        ax1.grid()

        # Phase plane
        ax2.plot(X, Y, color="blue")
        ax2.set_title('Phase plane')
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