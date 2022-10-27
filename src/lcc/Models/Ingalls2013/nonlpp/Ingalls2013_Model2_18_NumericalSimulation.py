'''
Ingalls 2012 Mathematical Modelling in Systems Biology

Numerical Simulation
Network  2.18,   p.35
Model   2.19,   p.35

'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import UnitTestHelper

def MMequation(Vmax, S, KM):
    return (Vmax * S) / (KM + S)

class FModel():
    def __init__(self, InA=0, InB=10, Ink1=9, Inkrev1=12, Ink2=2):
        # Initial Concentrations
        self.A = InA
        self.B = InB
        self.C = 0

        # Constants
        self.k1 = Ink1
        self.krev1 = Inkrev1
        self.k2 = Ink2

        # Dataset
        self.Data_A = list()
        self.Data_B = list()
        self.Data_C = list()
        self.Data_Time = list()

        # Set initial values
        self.InitializeSimStepZero()

    def InitializeSimStepZero(self):
        self.AppendData()

    def UpdateConcentration(self, dA, dB, dC):
        self.A += dA
        self.B += dB
        self.C += dC

    def AppendData(self):
        self.Data_A.append(self.A)
        self.Data_B.append(self.B)
        self.Data_C.append(self.C)

    def Derivative_A(self):
        return -self.k1 * self.A + self.krev1 * self.B

    def Derivative_B(self):
        return self.k1 * self.A - self.krev1 * self.B - self.k2 * self.B

    def Derivative_C(self):
        return self.k2 * self.B

    def Run(self, SimSteps=500, TimeResolution=100):
        i = 0
        while i < SimSteps:
            # Calculate new
            dA = self.Derivative_A() / TimeResolution
            dB = self.Derivative_B() / TimeResolution
            dC = self.Derivative_C() / TimeResolution

            self.UpdateConcentration(dA, dB, dC)
            self.AppendData()
            self.Data_Time.append(i/TimeResolution)

            i += 1

    def PlotData(self):
        X = self.Data_A
        Y = self.Data_B
        Z = self.Data_C

        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)

        # Dynamics
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)

        ax1.plot(X, 'r-', label="[A]")
        ax1.plot(Y, 'b-', label="[B]")
        ax1.plot(Z, 'g-', label="[C]")
        ax1.set_title('Dynamics')
        ax1.set_xlabel('Time(a.u.)')
        ax1.set_ylabel('Concentration(a.u.)')
        ax1.legend(loc='upper left')
        ax1.grid()

        # Phase plane
        ax2.plot(X, Y, color="blue")
        ax2.set_title('Phase plane')
        ax2.set_xlabel("[A]")
        ax2.set_ylabel("[B]")
        ax2.grid()

        plt.show()

def main(args):
    Model = FModel()
    Model.Run()

    if args.print_data:
        UnitTestHelper.PrintData(Model)
    else:
        Model.PlotData()

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--print-data', dest='print_data', action='store_true', default=False,
                        help="Print simulation data to screen")
    args = parser.parse_args()
    main(args)
