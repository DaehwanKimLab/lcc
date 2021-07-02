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

def Write_BiologicalEventRateFunction(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("class FBiologicalEventRateFunction():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.Rate_Min", 0)
            Writer.Variable_("self.Rate_Max", 0)
            Writer.BlankLine()

            Writer.Statement("super().__init__()")
            Writer.BlankLine()

        with Writer.Statement("def LoadRateMinMax(self, Rate_Min, Rate_Max):"):
            Writer.Reshape__("self.Rate_Min", "Rate_Min", [])
            Writer.Reshape__("self.Rate_Max", "Rate_Max", [])
            Writer.BlankLine()

        with Writer.Statement("def DetermineRate(self):"):
            Writer.Statement("self.Rate = tf.random.uniform(shape=[1], minval=self.Rate_Min, maxval=self.Rate_Max, dtype=tf.int32)")
            Writer.ReturnVar("self.Rate")
            Writer.BlankLine()
