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

def Write_BiochemicalReactionRateFunction(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("class FBiochemicalReactionRateFunction():"):
        with Writer.Function_("__init__"):
            # Writer.Variable_("self.Stoich", 0)
            # Writer.Variable_("self.Rate", 0)
            # Writer.Variable_("self.Count", 0)
            # Writer.Variable_("self.DeltaCount", 0)
            Writer.Pass_____()
            Writer.BlankLine()

