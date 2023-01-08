# BSD 3-Clause License
# © 2022, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee and Daehwan Kim

import sys
from datetime import datetime
import numpy as np
import math


class Reaction:
    def __init__(self):
        None

    def Specification(self, Molecules):
        return {}


# reaction NAD+Synthesis(nicotinamide + PRPP + 2 ATP -> NAD+)
class NADPlusSynthesis(Reaction):
    def __init__(self):
        self.Input = [["nicotinamide", 1], ["PRPP", 1], ["ATP", 2]]
        self.Output = [["NAD+", 1]]

    def Specification(self, Molecules):
        return {}

# reaction NADP+Synthesis(NAD+ + ATP -> NADP+ + ADP)
class NADPPlusSynthesis(Reaction):
    def __init__(self):
        self.Input = [["NAD+", 1], ["ATP", 1]]
        self.Output = [["NADP+", 1], ["ADP", 1]]

    def Specification(self, Molecules):
        return {}

# reaction CoASynthesis(pantothenate + cysteine + 3 ATP + CTP + CO2 -> CoA + ADP + CMP + 3 PPi)
class CoASynthesis(Reaction):
    def __init__(self):
        self.Input = [["pantothenate", 1], ["cysteine", 1], ["ATP", 3], ["CTP", 1], ["CO2", 1]]
        self.Output = [["CoA", 1], ["ADP", 1], ["CMP", 1], ["PPi", 3]]

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
NADH = 8.3mM;
NAD = 2.6mM;
pyruvate = 0.39mM;
"""
class Glycolysis(Reaction):
    def __init__(self):
        self.Input = [["G6P", 1], ["ADP", 2], ["NAD+", 2]]
        self.Output = [["pyruvate", 1], ["NADH", 2], ["ATP", 2]]

    def Specification(self, Molecules):
        cg = 0.019
        ATP = Molecules["ATP"]
        ADP = Molecules["ADP"]
        dATP = (ADP / ATP) * cg
        dADP = -dATP
        dNADH = dATP
        dNADPlus = -dNADH
        dG6P = -dATP / 2
        dpyruvate = dG6P * 2
        return {"ATP" : dATP, "ADP" : dADP, "NADH" : dNADH, "NAD+" : dNADPlus,
                "G6P" : dG6P, "pyruvate" : dpyruvate}


class Simulator():
    def __init__(self):
        self.Reactions = []
        self.Reactions.append(NADPlusSynthesis())
        self.Reactions.append(NADPPlusSynthesis())
        self.Reactions.append(CoASynthesis())
        self.Reactions.append(Glycolysis())

        self.PermanentMolecules = {
            "G6P"      : 8.8  * 1e-3,
            }

        self.Molecules = {
            "G6P"      : 8.8  * 1e-3,
            "ATP"      : 9.6  * 1e-3,
            "ADP"      : 0.56 * 1e-3,
            "NADH"     : 8.3  * 1e-3,
            "NAD+"     : 2.6  * 1e-3,
            "pyruvate" : 0.39 * 1e-3
        }

        self.InitializeMolecules()

    def InitializeMolecules(self):
        AdditionalMolecules = set()
        for Reaction in self.Reactions:
            for Molecule, Num in Reaction.Input + Reaction.Output:
                AdditionalMolecules.add(Molecule)

        for Molecule in AdditionalMolecules:
            if Molecule not in self.Molecules:
                self.Molecules[Molecule] = 0.0

    def UpdateMolecules(self, dMolecules, DeltaTime):
        for dMolecule, dConc in dMolecules.items():
            assert dMolecule in self.Molecules
            dConc = dConc * DeltaTime            
            assert dConc + self.Molecules[dMolecule] >= 0.0
            self.Molecules[dMolecule] += dConc

    def Adjust(self):
        for Molecule, Conc in self.PermanentMolecules.items():
            assert Molecule in self.Molecules
            self.Molecules[Molecule] = Conc

    def Info(self):
        def Conc2Str(Conc):
            if Conc >= 1e-1:
                Str = "{:.2f}".format(Conc)
            elif Conc >= 1e-4:
                Str = "{:.2f}mM".format(Conc * 1e3)
            elif Conc >= 1e-7:
                Str = "{:.2f}uM".format(Conc * 1e6)
            else:
                Str = "{:.2f}nM".format(Conc * 1e9)
            return Str

        for Molecule, Conc in self.Molecules.items():
            ConStr = Conc2Str(Conc)
            print("{:<16}: {}".format(Molecule, ConStr))
            

    def Simulate(self, DeltaTime = 0.01):
        self.Info()

        Iter = 0
        while Iter < 5:
            for Reaction in self.Reactions:
                dMolecules = Reaction.Specification(self.Molecules)
                self.UpdateMolecules(dMolecules, DeltaTime)

            self.Adjust()

            Iter += 1
            print("\n\n")
            print("Iteration {}".format(Iter))

            self.Info()

            


if __name__ == '__main__':
    Sim = Simulator()
    Sim.Simulate()
