'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Cooperativity
Network  3.16,   p.49
Figure 3.10,    p.62

'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import numpy as np

def HillFunction(X, K, n):
    return X**n / (K**n + X**n)

class FNetwork():
    def __init__(self, Inn=4):
        # Constants
        self.n = Inn   # range of integers tested

        # Dataset
        self.Data_X = None
        self.Data_Y = None
        self.Dataset = list()

    def ReinitializeData(self):
        self.Data_X = list()
        self.Data_Y = list()

    def AppendData(self, X, Y):
        self.Data_X.append(X)
        self.Data_Y.append(Y)

    def Model(self, n, XRange):
        for i in range(XRange):
            X = i
            Y = HillFunction(X, (n ** 2) * 5, n)
            self.AppendData(X, Y)

    def Run(self, XRange=200):

        for n in range(self.n):

            self.ReinitializeData()
            self.Model(n+1, XRange)
            self.Dataset.append([self.Data_X, self.Data_Y])

    def PlotData(self):
        # Dynamics
        for n, [X, Y] in enumerate(self.Dataset):


            plt.plot(X, Y, label="n = %s" % str(n+1))
            # plt.plot(np.log(X), Y, label="n = %s" % str(n+1))
            plt.title('Hill functions')
            plt.xlabel('X: Ligand concentration (a.u.)')
            plt.ylabel('Y: Fraction of Bound Protein')
            plt.legend(loc='upper left')
            plt.grid()

        plt.show()

def main():
    Network = FNetwork()
    Network.Run()
    Network.PlotData()

if __name__ == '__main__':

    main()