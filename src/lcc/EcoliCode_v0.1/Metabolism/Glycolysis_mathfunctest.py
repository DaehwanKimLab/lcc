from argparse import ArgumentParser
import matplotlib.pyplot as plt


class FNetwork():
    def __init__(self,
                 G6P=20,
                 ADP=5,
                 NAD=10,
                 Pi=15,
                 pyruvate=0,
                 NADH=1,
                 ATP=10,
                 cg=10):

        # Initial Concentrations
        self.G6P = G6P
        self.ADP = ADP
        self.NAD = NAD
        self.Pi = Pi
        self.pyruvate = pyruvate
        self.NADH = NADH
        self.ATP = ATP

        # Constants
        self.cg = cg

        # Dataset
        self.Data_G6P = list()
        self.Data_ADP = list()
        self.Data_NAD = list()
        self.Data_Pi = list()
        self.Data_pyruvate = list()
        self.Data_NADH = list()
        self.Data_ATP = list()

        # Set initial values
        self.InitializeSimStepZero()

    def InitializeSimStepZero(self):
        Data = self.G6P, \
               self.ADP, \
               self.NAD, \
               self.Pi, \
               self.pyruvate, \
               self.NADH, \
               self.ATP
        self.AppendData(Data)

    def AppendData(self, Data):
        self.Data_G6P.append(Data[0])
        self.Data_ADP.append(Data[1])
        self.Data_NAD.append(Data[2])
        self.Data_Pi.append(Data[3])
        self.Data_pyruvate.append(Data[4])
        self.Data_NADH.append(Data[5])
        self.Data_ATP.append(Data[6])

    def Simulation(self, t):
        G6P = self.Data_G6P[-1]
        ADP = self.Data_ADP[-1]
        NAD = self.Data_NAD[-1]
        Pi = self.Data_Pi[-1]
        pyruvate = self.Data_pyruvate[-1]
        NADH = self.Data_NADH[-1]
        ATP = self.Data_ATP[-1]

        dATP = ADP / ATP * self.cg * t
        dADP = - dATP
        dPi = -dATP
        dNADH = dATP
        dNAD = -dNADH
        dG6P = - dATP / 2
        dpyruvate = - dG6P * 2

        return G6P + dG6P, \
               ADP + dADP, \
               NAD + dNAD, \
               Pi + dPi, \
               pyruvate + dpyruvate, \
               NADH + dNADH, \
               ATP + dATP

    def Run(self, SimSteps, TimeResolution):
        i = 0
        Flat = 0.000001  # Steady state threshold

        while i < SimSteps:
            NewData = self.Simulation(1 / TimeResolution)
            self.AppendData(NewData)
            if (abs(self.Data_G6P[-1] - self.Data_G6P[-2]) < Flat):
                break
            i += 1

    def PlotData(self):
        G6P = self.Data_G6P
        ADP = self.Data_ADP
        NAD = self.Data_NAD
        Pi = self.Data_Pi
        pyruvate = self.Data_pyruvate
        NADH = self.Data_NADH
        ATP = self.Data_ATP

        plt.plot(G6P, 'red', label="[G6P]")
        plt.plot(ADP, 'green', label="[ADP]")
        plt.plot(NAD, 'brown', label="[NAD]")
        plt.plot(Pi, 'yellow', label="[Pi]")
        plt.plot(pyruvate, 'magenta', label="[pyruvate]")
        plt.plot(NADH, 'orange', label="[NADH]")
        plt.plot(ATP, 'blue', label="[ATP]")
        plt.title('Glycolysis')
        plt.xlabel('Time (a.u.)')
        plt.ylabel('Concentration (a.u.)')
        plt.legend(loc='center right')
        plt.grid()

        plt.show()


def main():
    Network = FNetwork()
    Network.Run(1000, 100)
    Network.PlotData()


if __name__ == '__main__':
    main()
