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


def Write_RateGaugeModelOnly(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FRateGaugeModelOnly(FReactionExecution):"):
        with Writer.Function_("__init__"):
            Writer.Statement("super().__init__()")
            Writer.BlankLine()
