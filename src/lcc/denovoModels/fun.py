import numpy as np

class EryingPolaniEquation():
    def __init__(self, deltaG = 30500, transmissionCoefficient = 1, temp = 310.15):

        # Constants:
        self.kBolzmann =  1.380649e-23 # J/K
        self.R = 8.31446261815324 # gas constant: J*K^-1*mol^-1
        self.hPlanck = 6.62607015e-34 # J/Hz
        
        # Inputs
        self.transmissionCoefficient = transmissionCoefficient 
        self.deltaG = deltaG #In J/mol
        self.T = temp # Default 37C in Kelvin
        
    def SolveForReactionRate(self):
        # Erying-polani equation
        reactionRate = (self.transmissionCoefficient * self.kBolzmann * self.T)/self.hPlanck * np.exp(-self.deltaG/(self.R * self.T)) 
        print(f'Reaction rate k: {reactionRate}')

def main():
    Eq = EryingPolaniEquation()
    Eq.SolveForReactionRate()

if __name__ == '__main__':
    main()