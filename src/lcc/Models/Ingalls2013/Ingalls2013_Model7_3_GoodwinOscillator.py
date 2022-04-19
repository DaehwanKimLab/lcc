'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Goodwin Oscillator
Network  7.16,   p.213
Model   7.22,    p.212
Figure  7.17,    p.214

'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

'''Figure 7.17: The Goodwin oscillator. A. This simulation shows relaxation to sustained (limit cycle)
oscillations. B. A phase portrait showing convergence to a limit cycle in the three-dimensional phase space.

Parameter values are a = 360 (concentration · time−1), k = 1.368 (concentration), b = 1 (time−1), alpha = 1
(time−1), beta = 0.6 (time−1), gamma = 1 (time−1), delta = 0.8 (time−1), n = 12. Units are arbitrary.'''

uM = 1e-6
nM = 1e-9

Unit = nM

Mol2Count = 6.0221409e+23
Count2Mol = 1.0 / Mol2Count



def dX_NumericalSimulation(a, b, k, n, alpha, beta, gamma, delta, X, Y, Z):
    return a / ((k ** n) + (Z ** n)) - b * X

def dY_NumericalSimulation(a, b, k, n, alpha, beta, gamma, delta, X, Y, Z):
    return (alpha * X) - (beta * Y)

def dZ_NumericalSimulation(a, b, k, n, alpha, beta, gamma, delta, X, Y, Z):
    return (gamma * Y) - (delta * Z)


class FModel():
    def __init__(self, a=360, b=1, k=1.368, n=12, alpha=1, beta=0.6, gamma=1, delta=0.8, X=0, Y=0, Z=0):

        # Initial Concentrations
        self.X = X
        self.Y = Y
        self.Z = Z

        # Constants
        self.a = a
        self.b = b
        self.k = k
        self.n = n
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta

        # Dataset
        self.Data_X = list()
        self.Data_Y = list()
        self.Data_Z = list()

        # Set initial values
        self.InitializeSimStepZero()

        self.TimeResolution = 100

        self.ShowTrueUnits = True

    def InitializeSimStepZero(self):
        self.AppendData(self.X, self.Y, self.Z)

    def AppendData(self, X, Y, Z):
        self.Data_X.append(X)
        self.Data_Y.append(Y)
        self.Data_Z.append(Z)

    def Run(self, SimSteps=3500, TimeResolution=100):
        Flat = 0.000001  # Steady state threshold

        self.TimeResolution = TimeResolution

        i = 0

        while i < SimSteps:
            Time = i / self.TimeResolution
            self.Simulate()

            # Debugging
            self.Debugging(i)

            i += 1

    def Debugging(self, SimStep):
        print("SimStep: {} |\tX: {} |\tY: {} |\tZ: {}".format(SimStep, self.Data_X[-1], self.Data_Y[-1], self.Data_Z[-1]))

    def Simulate(self, SimUnitTime=None):
        X = self.Data_X[-1]
        Y = self.Data_Y[-1]
        Z = self.Data_Z[-1]

        if SimUnitTime == None:
            SimUnitTime = 1.0 / self.TimeResolution

        dX = dX_NumericalSimulation(self.a, self.b, self.k, self.n, self.alpha, self.beta, self.gamma, self.delta, X, Y, Z) * SimUnitTime
        dY = dY_NumericalSimulation(self.a, self.b, self.k, self.n, self.alpha, self.beta, self.gamma, self.delta, X, Y, Z) * SimUnitTime
        dZ = dZ_NumericalSimulation(self.a, self.b, self.k, self.n, self.alpha, self.beta, self.gamma, self.delta, X, Y, Z) * SimUnitTime

        self.AppendData(X + dX, Y + dY, Z + dZ)

    def PlotData(self):
        X = self.Data_X
        Y = self.Data_Y
        Z = self.Data_Z

        UnitTxt = '(a.u.)'
        #
        # if self.ShowTrueUnits:
        #     Am = [x / Unit for x in Am]
        #     L = [x / Unit for x in L]
        #     if Unit == uM:
        #         UnitTxt = '(uM)'
        #     elif Unit == nM:
        #         UnitTxt = '(nM)'

        T = list(np.array(range(len(X))) / self.TimeResolution)

        # Plotting options
        self.Plot(T, X, Y, Z, UnitTxt)

        # # Debugging purposes
        # else:
        #     fig = plt.figure()
        #     fig.subplots_adjust(wspace=0.2, hspace=0.3)
        #
        #     PlotBegin = 5000
        #
        #     # Dynamics
        #     ax1 = fig.add_subplot(2, 2, 1)
        #     ax1.plot(T[PlotBegin:], Am[PlotBegin:], 'b-', label="[Am]")  # RL
        #     ax1.set_title('Bacterial Chemotaxis')
        #     ax1.set_xlabel('Time (s)')
        #     ax1.set_ylabel('Concentration ' + UnitTxt)
        #     ax1.legend(loc='upper right')
        #     ax1.grid()
        #
        #     ax2 = fig.add_subplot(2, 2, 2)
        #     ax2.plot(T[PlotBegin:], Am[PlotBegin:], 'b-', label="[Am]")  # RL
        #     ax2.plot(T[PlotBegin:], A[PlotBegin:], 'r-', label="[A]")  # RL
        #     ax2.set_ylabel('Concentration ' + UnitTxt)
        #     ax2.legend(loc='upper right')
        #     ax2.grid()
        #
        #     ax3 = fig.add_subplot(2, 2, 3)
        #     ax3.plot(T[PlotBegin:], AmL[PlotBegin:], 'b-', label="[AmL]")  # RL
        #     ax3.plot(T[PlotBegin:], AL[PlotBegin:], 'r-', label="[AL]")  # RL
        #     ax3.set_xlabel('Time (a.u.)')
        #     ax3.set_ylabel('Concentration ' + UnitTxt)
        #     ax3.legend(loc='upper right')
        #     ax3.grid()
        #
        #     ax4 = fig.add_subplot(2, 2, 4)
        #     ax4.plot(T[PlotBegin:], BP[PlotBegin:], 'r-', label="[BP]")  # RL
        #     ax4.set_xlabel('Time (a.u.)')
        #     ax4.set_ylabel('Concentration ' + UnitTxt)
        #     ax4.legend(loc='upper right')
        #     ax4.grid()
        #
        #     plt.show()

    def Plot(self, T, X, Y, Z, UnitTxt):
        fig = plt.figure()

        PlotBegin = 0

        # Dynamics
        ax1 = plt.axes()
        ax1.plot(T[PlotBegin:], X[PlotBegin:], 'r-', label="mRNA X")
        ax1.plot(T[PlotBegin:], Y[PlotBegin:], 'g-', label="enzyme Y")
        ax1.plot(T[PlotBegin:], Z[PlotBegin:], 'b-', label="metabolite Z")
        ax1.set_title('Goodwin Oscillator')
        ax1.set_xlabel('Time (a.u.)')
        ax1.set_ylabel('Concentration ' + UnitTxt)
        ax1.legend(loc='upper right')
        ax1.grid()

        plt.show()

def main():
    Model = FModel()
    Model.Run()
    Model.PlotData()

if __name__ == '__main__':
    main()