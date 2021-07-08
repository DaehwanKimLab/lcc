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

def Write_ReactionExecution(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FReactionExecution():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.Stoich", 0)
            Writer.Variable_("self.Rate", 0)
            Writer.Variable_("self.Count", 0)
            Writer.Variable_("self.DeltaCount", 0)

            Writer.BlankLine()

        with Writer.Statement("def LoadStoichMatrix(self, Stoich):"):
            Writer.Overwrite("self.Stoich", "Stoich")
            Writer.BlankLine()

        with Writer.Statement("def LoadRateMatrix(self, Rate):"):
            Writer.Overwrite("self.Rate", "Rate")
            Writer.BlankLine()

        with Writer.Statement("def LoadCountMatrix(self, Count):"):
            Writer.Overwrite("self.Count", "Count")
            Writer.BlankLine()

        # RunReactions method may be overwritten
        with Writer.Statement("def MultiplyRXNMatrices(self):"):
            # Multiplies all Reaction Stoichiometry and Rate matrices.
            Writer.Statement("self.DeltaCount = tf.linalg.matmul(self.Stoich, self.Rate)")
            Writer.DebugPVar("self.DeltaCount")
            Writer.BlankLine()

        with Writer.Statement("def RoundDeltaCountMatrix(self):"):
            # Multiplies all Reaction Stoichiometry and Rate matrices.
            Writer.OperRound("self.DeltaCount", "self.DeltaCount")
            Writer.Cast_____("self.DeltaCount", "self.DeltaCount", 'int32')
            Writer.BlankLine()

        with Writer.Statement("def AddCountMatrices(self):"):
            Writer.OperMXAdd("FinalCount", "self.Count", "self.DeltaCount")
            Writer.DebugPVar("FinalCount")
            Writer.ReturnVar("FinalCount")
            Writer.BlankLine()

        with Writer.Statement("def RunReactions(self):"):
            Writer.Statement("self.MultiplyRXNMatrices()")
            Writer.Statement("self.RoundDeltaCountMatrix()")