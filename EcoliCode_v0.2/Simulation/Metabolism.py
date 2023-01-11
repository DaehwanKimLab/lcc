# BSD 3-Clause License
# © 2022, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import sys
from datetime import datetime
import numpy as np
import math


class Reaction:
    def __init__(self):
        self.Input = dict()
        self.Output = dict()
        self.Stoich = dict()

    def Specification(self, Molecules):
        return {}

# reaction NAD+Synthesis(nicotinamide + PRPP + 2 ATP -> NAD+)
class NADPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"nicotinamide": 1, "PRPP": 1, "ATP": 2}
        self.Output = {"NAD+": 1}

    def Specification(self, Molecules):
        return {}

# reaction NADP+Synthesis(NAD+ + ATP -> NADP+ + ADP)
class NADPPlusSynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"NAD+": 1, "ATP": 1}
        self.Output = {"NADP+": 1, "ADP": 1}

    def Specification(self, Molecules):
        return {}

# reaction CoASynthesis(pantothenate + cysteine + 3 ATP + CTP + CO2 -> CoA + ADP + CMP + 3 PPi)
class CoASynthesis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"pantothenate": 1, "cysteine": 1, "ATP": 3, "CTP": 1, "CO2": 1}
        self.Output = {"CoA": 1, "ADP": 1, "CMP": 1, "PPi": 3}

    def Specification(self, Molecules):
        return {}

"""
reaction Glycolysis(G6P + 2 ADP + 2 NAD+ -> 2 pyruvate + 2 NADH + 2 ATP)
specification
{
    float cg = 0.019M;
    ATP'  = (ADP / ATP) * cg;
    ADP'  = - ATP';
    NADH' = ATP';
    NAD' = - NADH';
    G6P'  = - ATP' / 2;
    pyruvate' = G6P' * 2;
}

G6P[:] = 8.8mM;
ATP = 9.6mM;
ADP = 0.56mM;
NADH = 83uM;
NAD = 2.6mM;
pyruvate = 0.39mM;
"""
class Glycolysis(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"G6P": 1, "ADP": 2, "NAD+": 2}
        self.Output = {"pyruvate": 1, "NADH": 2, "ATP": 2}

    def Specification(self, Molecules):
        c = 0.019
        ATP = Molecules["ATP"]
        ADP = Molecules["ADP"]
        dATP = (ADP / ATP) * c
        # dATP = (abs(ADP - 0.00056) / ATP) * c   # For homeostasis of ADP level
        # dATP = (ADP * abs(ATP - 0.0096) / ATP) * c   # For homeostasis of ATP level
        # dATP = ((ADP / ATP) - (0.00056 / 0.0096)) * ATP * c   # For homeostasis of ADP to ATP ratio
        return "ATP", dATP

class ATPConsumption(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"ATP": 1}
        self.Output = {"ADP": 1}

    def Specification(self, Molecules):
        # Cell Division ATP consumption: 4.35E-03
        c = 4.35E-03
        # c = 4.35E-05
        dATP = -c
        return "ATP", dATP

# For Debugging
class ADPConsumption(Reaction):
    def __init__(self):
        super().__init__()
        self.Input = {"ADP": 1}
        self.Output = {"ATP": 1}

    def Specification(self, Molecules):
        c = 4.35E-03
        dADP = -c
        return "ADP", dADP


class Simulator():
    def __init__(self):
        self.Reactions = []

        self.PermanentMolecules = {
            "G6P"      : 8.8  * 1e-3,
            }

        self.Molecules = {
            "G6P"      : 8.8  * 1e-3,
            "ATP"      : 9.6  * 1e-3,
            "ADP"      : 0.56 * 1e-3,
            "NADH"     : 83   * 1e-6,
            "NAD+"     : 2.6  * 1e-3,
            "pyruvate" : 0.39 * 1e-3,
        }

        # Metabolite concentrations, fluxes and free energies imply efficient enzyme usage
        # Junyoung O Park, Sara A Rubin, Yi-Fan Xu, Daniel Amador-Noguez, Jing Fan, Tomer Shlomi & Joshua D Rabinowitz
        # Nature Chemical Biology volume 12, pages482–489 (2016)
        # https://www.nature.com/articles/nchembio.2077#Sec21

        self.InitialConditions = self.Molecules.copy()

    def Initialize(self):
        self.InitializeStoich()
        self.InitializeMolecules()

    def InitializeStoich(self):
        for Reaction in self.Reactions:
            for Mol, Coeff in Reaction.Input.items():
                Reaction.Stoich[Mol] = -Coeff
            for Mol, Coeff in Reaction.Output.items():
                Reaction.Stoich[Mol] = Coeff

    def InitializeMolecules(self):
        AdditionalMolecules = set()
        for Reaction in self.Reactions:
            for Molecule, Num in Reaction.Stoich.items():
                AdditionalMolecules.add(Molecule)

        for Molecule in AdditionalMolecules:
            if Molecule not in self.Molecules:
                self.Molecules[Molecule] = 0.0

    def AddReaction(self, Reaction):
        self.Reactions.append(Reaction)

    def AdjustRefdCon(self, Reaction, RefMol, RefdConc):
        ''' Compare dConc of reference molecule to input concentrations and adjust reference dConc '''
        RefCoeff = Reaction.Stoich[RefMol]
        UnitdConc = RefdConc / RefCoeff

        Out = 0
        for Mol, Coeff in Reaction.Input.items():
            AdjustedConc = self.Molecules[Mol] / Coeff
            Out = min(AdjustedConc, UnitdConc)

        return Out * RefCoeff

    def DeterminedConc(self, Reaction, RefMol, RefdConc):
        UnitdConc = RefdConc / Reaction.Stoich[RefMol]
        dConc = dict()
        for Mol, Coeff in Reaction.Stoich.items():
            dConc[Mol] = UnitdConc * Coeff
        return dConc

    def RunReaction(self, Reaction):
        RefMol, RefdConc = Reaction.Specification(self.Molecules)
        RefdConc = self.AdjustRefdCon(Reaction, RefMol, RefdConc)
        return self.DeterminedConc(Reaction, RefMol, RefdConc)

    def UpdateMolecules(self, dMolecules, DeltaTime):
        for dMolecule, dConc in dMolecules.items():
            assert dMolecule in self.Molecules
            dConc = dConc * DeltaTime            
            assert dConc + self.Molecules[dMolecule] >= 0.0, '{} \t| Conc:{}, \t dConc:{}'.format(dMolecule, self.Conc2Str(self.Molecules[dMolecule]), self.Conc2Str(dConc))
            self.Molecules[dMolecule] += dConc

    def Adjust(self):
        for Molecule, Conc in self.PermanentMolecules.items():
            assert Molecule in self.Molecules
            self.Molecules[Molecule] = Conc

    def Conc2Str(self, Conc):
        AbsConc = abs(Conc)
        if AbsConc >= 1e-1:
            Str = "{:.3f}  M".format(Conc)
        elif AbsConc >= 1e-4:
            Str = "{:.3f} mM".format(Conc * 1e3)
        elif AbsConc >= 1e-7:
            Str = "{:.3f} uM".format(Conc * 1e6)
        elif AbsConc >= 1e-10:
            Str = "{:.3f} nM".format(Conc * 1e9)
        else:
            Str = "{:.3f} pM".format(Conc * 1e12)
        return Str

    def Info(self):
        for Molecule, Conc in self.Molecules.items():
            ConStr = self.Conc2Str(Conc)
            print("{:<16}: {:>10}".format(Molecule, ConStr))

    def Summary(self):
        print("{:<16} {:>10} {:>10} {:>10} {:<10}".format("", "Initial", "Final", "dConc", "(Permanent)"))
        for Molecule, Conc in self.Molecules.items():
            dConc = Conc - self.InitialConditions[Molecule]

            InitConc = self.Conc2Str(self.InitialConditions[Molecule])
            FinalConc = self.Conc2Str(Conc)
            FinaldConc = self.Conc2Str(dConc)

            print("{:16} {:>10} {:>10} {:>10} {}".format(Molecule, InitConc, FinalConc, FinaldConc, "*" if Molecule in self.PermanentMolecules else ""))

    def Simulate(self, DeltaTime = 0.01, TotalTime = 10):
        print("-- Initial Conditions --")
        self.Info()

        Iter = 0
        while Iter < TotalTime:
            # # Debug
            # if Iter == 220:
            #     print("\n")
            #     print("Debugging:", Iter)

            for Reaction in self.Reactions:
                dMolecules = self.RunReaction(Reaction)
                self.UpdateMolecules(dMolecules, DeltaTime)

            self.Adjust()

            Iter += 1
            print("\n")
            print("-- Iteration {} --".format(Iter))

            self.Info()

        print("\n")
        print("-- Summary --")
        self.Summary()


if __name__ == '__main__':
    Sim = Simulator()

    # Add reactions

    # Sim.AddReaction(NADPlusSynthesis())
    # Sim.AddReaction(NADPPlusSynthesis())
    # Sim.AddReaction(CoASynthesis())
    Sim.AddReaction(Glycolysis())
    # Sim.AddReaction(ATPConsumption())
    # Sim.AddReaction(ADPConsumption())

    Sim.Initialize()
    Sim.Simulate(TotalTime=1000)
