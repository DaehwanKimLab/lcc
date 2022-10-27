'''
Ingalls 2013 Mathematical Modelling in Systems Biology

Limit Cycle Oscillations
Network  4.14,   p.93
Model   4.10,    p.93
Figure  4.16,    p.95

'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dS1_NumericalSimulation(k0, k1, K, n, S1, S2):
    return k0 - (k1 * (1 + (S2 / K) ** n) * S1)

def dS2_NumericalSimulation(k1, k2, K, n, S1, S2):
    return (k1 * (1 + (S2 / K) ** n) * S1) - (k2 * S2)

class FModel():
    def __init__(self, InS1=1.5, InS2=1, Ink0=8, Ink1=1, Ink2=5, InK=1, Inn=2.5):
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

        self.Data_Time = list()

        # Set initial values
        self.InitializeSimStepZero()

    def InitializeSimStepZero(self):
        self.AppendData(self.S1, self.S2)

    def AppendData(self, S1, S2):
        self.Data_S1.append(S1)
        self.Data_S2.append(S2)

    def Run(self, SimSteps=1000, TimeResolution=100):
        Flat = 0.000001 # Steady state threshold

        i = 0
        while i < SimSteps:
            # Calculate new
            S1 = self.Data_S1[-1]
            S2 = self.Data_S2[-1]

            dS1 = dS1_NumericalSimulation(self.k0, self.k1, self.K, self.n, S1, S2) / TimeResolution
            dS2 = dS2_NumericalSimulation(self.k1, self.k2, self.K, self.n, S1, S2) / TimeResolution

            self.AppendData(self.Data_S1[-1] + dS1, self.Data_S2[-1] + dS2)
            self.Data_Time.append(i/TimeResolution)

            if (abs(self.Data_S1[-1] - self.Data_S1[-2]) < Flat):
                break
            i += 1

        return self.Data_S1, self.Data_S2

    def PlotData(self):
        X = self.Data_S1
        Y = self.Data_S2

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)

        # Dynamics
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)

        ax1.plot(X, 'b-', label="[S1]")
        ax1.plot(Y, 'r-', label="[S2]")
        ax1.set_title('Limit Cycle Oscillations (Sustained)')
        ax1.set_xlabel('Time (a.u.)')
        ax1.set_ylabel('Concentration (a.u.)')
        ax1.legend(loc='upper right')
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


# Animation option
Animate = False
# Animate = True

PlotType = 1  # Dynamics
# PlotType = 2 # PhasePlane

if Animate:
    Model = FModel()
    X, Y = Model.Run()

    assert len(X) == len(Y)
    skip = 3

    X = X[::skip]
    Y = Y[::skip]

    i_data = []
    x_data = []
    y_data = []

    fig, ax = plt.subplots()
    line1, = ax.plot(0, 0, label="[Substrate]", color='g')
    line2, = ax.plot(0, 0, label="[Product]", color='r')
    line3, = ax.plot(0, 0, color='b')

    def animation_frame(i):
        i_data.append(i)
        x_data.append(X[i])
        y_data.append(Y[i])

        if PlotType == 1:
            # dynamics
            ax.set_xlim(0, i + 10)
            ax.set_ylim(0, max(X) * 1.1)
            ax.set_ylabel('Concentration (a.u.)')
            ax.set_xlabel('Time (a.u.)')
            ax.set_title('Persistent Oscillation: Dynamics')
            ax.legend(loc='lower right')
            line1.set_data(i_data, x_data)
            line2.set_data(i_data, y_data)
            return [line1, line2]

        if PlotType == 2:
            # phase plane
            ax.set_xlim(0, max(X) * 1.1)
            ax.set_ylim(0.5, max(Y) * 1.1)
            line3.set_data(x_data, y_data)
            ax.set_ylabel('[Substrate] (a.u.)')
            ax.set_xlabel('[Product] (a.u.)')
            ax.set_title('Persistent Oscillation: Phase Plane')
            return line3,

    animation = FuncAnimation(fig, animation_frame, frames=1000, interval=1)
    plt.show()