'''
Ingalls 2013 Mathematical Modelling in Systems Biology

MichaelisMenten
Network  3.2,   p.49
Figure 3.3, p.50


'''
from argparse import ArgumentParser
import matplotlib.pyplot as plt

def MMequation(Vmax, S, kM):
    return (Vmax * S) / (kM + S)

def Derivative_S(S, E, C, P, k1, krev1, k2):
    return -(k1 * S * E) + (krev1 * C)

def Derivative_E(S, E, C, P, k1, krev1, k2):
    return (krev1 * C) - (k1 * S * E) + (k2 * C)

def Derivative_C(S, E, C, P, k1, krev1, k2):
    return -(krev1 * C) + (k1 * S * E) - (k2 * C)

def Derivative_P(S, E, C, P, k1, krev1, k2):
    return k2 * C

class FNetwork():
    def __init__(self, InS=5, InE=1, InC=0, InP=0, Ink1=30, Inkrev1=1, Ink2=10):
        # Initial Concentrations
        self.S = InS
        self.E = InE
        self.C = InC
        self.P = InP
        
        # Constants
        self.k1 = Ink1
        self.krev1 = Inkrev1
        self.k2 = Ink2
        
        self.kM = (self.krev1 + self.k2) / self.k1

        # Dataset
        self.Data_S = list()
        self.Data_E = list()
        self.Data_C = list()
        self.Data_P = list()

        self.Data_S_MM = list()
        self.Data_P_MM = list()

        # Set initial values
        self.InitializeSimStepZero()

    def InitializeSimStepZero(self):
        self.Data_S.append(self.S)
        self.Data_E.append(self.E)
        self.Data_C.append(self.C)
        self.Data_P.append(self.P)        
        
        self.Data_S_MM.append(self.S)
        self.Data_P_MM.append(self.P)

    # def UpdateConcentration(self, dS, dE, dC, dP):
    #     self.S += dS
    #     self.E += dE
    #     self.C += dC
    #     self.P += dP

    def AppendData_Full(self, S, E, C, P):
        self.Data_S.append(S)
        self.Data_E.append(E)
        self.Data_C.append(C)
        self.Data_P.append(P)
     
    def Model_Full(self, TimeResolution):
        S = self.Data_S[-1]
        E = self.Data_E[-1]
        C = self.Data_C[-1]
        P = self.Data_P[-1]
        
        dS = Derivative_S(S, E, C, P, self.k1, self.krev1, self.k2) / TimeResolution
        dE = Derivative_E(S, E, C, P, self.k1, self.krev1, self.k2) / TimeResolution
        dC = Derivative_C(S, E, C, P, self.k1, self.krev1, self.k2) / TimeResolution
        dP = Derivative_P(S, E, C, P, self.k1, self.krev1, self.k2) / TimeResolution
        
        self.AppendData_Full(S + dS, E + dE, C + dC, P + dP)
        
    def AppendData_Reduced(self, S, P):
        self.Data_S_MM.append(S)
        self.Data_P_MM.append(P)
        
    def Model_Reduced(self, TimeResolution):
        S = self.Data_S_MM[-1]
        P = self.Data_P_MM[-1]
        Vmax = self.k2 * self.E

        Rate_MM = MMequation(Vmax, S, self.kM) / TimeResolution
        
        self.AppendData_Reduced(S - Rate_MM, P + Rate_MM)

    def Run(self, SimSteps=50000, TimeResolution=500):
        i = 0
        Flat = 0.00000001 # Steady state threshold

        while i < SimSteps:
            # Calculate new
            self.Model_Full(TimeResolution)
            self.Model_Reduced(TimeResolution)
            
            if (abs(self.Data_S[-1] - self.Data_S[-2]) < Flat):
                break
            i += 1

    def PlotData(self):
        S = self.Data_S
        E = self.Data_E
        C = self.Data_C
        P = self.Data_P

        S_MM = self.Data_S_MM
        P_MM = self.Data_P_MM
    
        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)

        # Dynamics
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)

        ax1.plot(S, 'r-', label="[S]")
        ax1.plot(E, 'b-', label="[E]")
        ax1.plot(C, 'g-', label="[C]")
        ax1.plot(P, 'm-', label="[P]")
        ax1.set_title('Dynamics')
        ax1.set_xlabel('Time(a.u.)')
        ax1.set_ylabel('Concentration(a.u.)')
        ax1.legend(loc='upper left')
        ax1.grid()
        
        ax2.plot(S, 'r-', label="[S]_full")
        ax2.plot(P, 'b-', label="[P]_full")
        ax2.plot(S_MM, 'g-', label="[C]_reduced")
        ax2.plot(P_MM, 'm-', label="[P]_reduced")
        ax2.set_title('Dynamics')
        ax2.set_xlabel('Time(a.u.)')
        ax1.set_ylabel('Concentration(a.u.)')
        ax2.legend(loc='upper left')
        ax2.grid()
    
        plt.show()

def main():
    Network = FNetwork()
    Network.Run()
    Network.PlotData()
    
if __name__ == '__main__':

    main()