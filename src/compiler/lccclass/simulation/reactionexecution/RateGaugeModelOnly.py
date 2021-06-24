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


def Write_RateGaugeModelOnly(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("class FRateGaugeModelOnly(FReactionExecution):"):
        with Writer.Statement("def __init__(self):"):
            Writer.Statement("super().__init__()")
            Writer.BlankLine()

        with Writer.Statement("def MultiplyMatrix(self):"):
            # Multiplies all Reaction Stoichiometry and Rate matrices.
            Writer.Statement("self.DeltaCount = tf.linalg.matmul(self.Stoich, self.Rate)")
            Writer.DebugPVar("self.DeltaCount")
            Writer.BlankLine()

