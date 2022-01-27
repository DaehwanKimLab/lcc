'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Activation-Inactivation Cycle
Network  6.6,   p.156
Model   6.1,    p.156
Figure  6.7,    p.157
'''

from argparse import ArgumentParser
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

def GetWstar(v1, v2, K1, K2):   # Eq 7 from Goldbeter and Koshland 1981, PNAS: https://www.pnas.org/content/pnas/78/11/6840.full.pdf
    return ((v1 / v2 - 1) - K2 * (K1 / K2 + v1 / v2) + ((v1 / v2 - 1 - K2 * (K1 / K2 + v1 / v2)) ** 2 + 4 * K2 * (v1 / v2 - 1) * (v1 / v2)) ** (1 / 2)) / (2 * (v1 / v2 - 1))

class FNetwork():
    def __init__(self, k1=1, k2=1, E2=1):
        # Constants
        self.k1 = k1   # range of integers tested
        self.k2 = k2   # range of integers tested
        self.E2 = E2
        self.v2 = self.k2 * self.E2

        # Dataset
        self.Data_E1 = None
        self.Data_Wstar = None
        self.Dataset = list()

    def ReinitializeData(self):
        self.Data_E1 = list()
        self.Data_Wstar = list()

    def AppendData(self, X, Y):
        self.Data_E1.append(X)
        self.Data_Wstar.append(Y)

    def Model(self, K, XRange):

        for i in range(XRange):
            E1 = (3 / XRange) * i   # E1
            v1 = self.k1 * E1
            Wstar = GetWstar(v1, self.v2, K, K)

            self.AppendData(E1, Wstar)

    def Run(self, XRange=5000):

        Karray = [0.01, 0.1, 1]
        for K in Karray:

            self.ReinitializeData()
            self.Model(K, XRange)
            self.Dataset.append([self.Data_E1, self.Data_Wstar])

    def PlotData(self):
        # Dynamics
        Karray = [0.01, 0.1, 1]

        for K, [X, Y] in zip(Karray, self.Dataset):
            plt.plot(X, Y, label="K1=K2=%s" % str(K))
            # plt.plot(np.log(X), Y, label="K1=K2=%s" % str(K))
            plt.title('Activation-inactivation cycle')
            plt.xlabel('Concentration of activating enzyme (a.u.)')
            plt.ylabel('Fraction of active protein')
            plt.legend(loc='upper left')
            plt.grid()

        plt.show()

def main():
    Network = FNetwork()
    Network.Run()
    Network.PlotData()

if __name__ == '__main__':

    main()
