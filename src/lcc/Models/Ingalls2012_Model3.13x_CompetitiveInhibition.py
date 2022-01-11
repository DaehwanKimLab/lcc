'''
Ingalls 2012 Mathematical Modelling in Systems Biology

3.2.1 Competitive Inhibition
Network  x  p.55
Equation    x   p.56
Figure 3.6, p.57


'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt

def MMEquation_CompetitiveInhibition(Vmax, S, kM, I, ki):
    return (Vmax * S) / (kM * (1 + I/ki) + S)

class FNetwork():
    def __init__(self, InE=1, InI=1, Ink1=5, Inkrev1=1, Ink2=8, Ink3=2, Inkrev3=1):
        # Initial Concentrations
        self.E = InE
        self.I = InI
        
        # Constants
        self.k1 = Ink1
        self.krev1 = Inkrev1
        self.k2 = Ink2
        self.k3 = Ink3
        self.krev3 = Inkrev3

        self.kM = (self.krev1 + self.k2) / self.k1
        self.ki = self.krev3/self.k3

        # Dataset
        self.Data_S = list()
        self.Data_Rate = list()
        self.Dataset = list()

    def InitializeSimStepZero(self):
        self.Data_S.append(0)
        self.Data_Rate.append(0)

    # def UpdateConcentration(self, dS, dE, dC, dP):
    #     self.S += dS
    #     self.E += dE
    #     self.C += dC
    #     self.P += dP
    def Model(self, S, IFoldChange, TimeResolution):
        I = self.I * IFoldChange
        Vmax = self.k2 * self.E

        Rate = MMEquation_CompetitiveInhibition(Vmax, S, self.kM, I, self.ki) / TimeResolution
        self.AppendData(S, Rate)
        
    def AppendData(self, S, Rate):
        self.Data_S.append(S)
        self.Data_Rate.append(Rate)

    def ReinitializeData(self):
        self.Data_S = list()
        self.Data_Rate = list()
        self.InitializeSimStepZero()

    def Run(self, SRange=100, TimeResolution=1, IFoldChange=1, N_Datasets=1):
        for n in range(N_Datasets):
            self.ReinitializeData()
            for S in range(SRange):
                if S == 0:
                    continue
                self.Model(S, n * IFoldChange, TimeResolution)
            self.Dataset.append([self.Data_S, self.Data_Rate])

    def PlotData(self, IFoldChange=1):
        # Dynamics
        for n, [S, Rate] in enumerate(self.Dataset):
            plt.plot(S, Rate, label="[I] = %s" % str(self.I * IFoldChange * n))
            # plt.set_title('Hill functions')
            # plt.set_xlabel('X: Ligand concentration (a.u.)')
            # plt.set_ylabel('Y: Fraction of Bound Protein')
            plt.legend(loc='upper left')
            plt.grid()

        plt.show()

def main():
    Network = FNetwork()
    IFoldChange = 5
    N_Datasets = 4
    Network.Run(IFoldChange=IFoldChange, N_Datasets=N_Datasets)
    Network.PlotData(IFoldChange=IFoldChange)
    
if __name__ == '__main__':

    main()