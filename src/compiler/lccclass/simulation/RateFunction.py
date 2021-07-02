'''
MasterProcess class

This class contains all of the elementary reaction information, including the following items:
    Inputs:
    - Participating molecule's index
    - Stoichiometry of molecules in the reaction
    - Rate of the reaction
    Outputs:
    - Master matrices of the inputs above

    "S" matrix
    Row: Reaction
    Column: Molecule
    Matrix Value: Stoichiometry

    "r" matrix
    Row: Reaction
    Matrix Value: Rate

'''
import abc

def Write_RateFunction(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("class FRateFunction():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.Rate", 0)
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def DetermineRate(self):"):
            Writer.Pass_____()
            Writer.BlankLine()
