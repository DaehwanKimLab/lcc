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


def Write_MasterProcess_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def Simulation_Init():"):
        Writer.Variable_("self.SimulationStep", 100) #
        Writer.Variable_("self.CellCycle", 1)