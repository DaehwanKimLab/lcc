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

def Write_PolymerizationRateFunction(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("class FPolymerizationRateFunction():"):
        with Writer.Statement("def __init__(self):"):

            Writer.Variable_("self.Mean", 0)
            # For normal distribution
            Writer.Variable_("self.SD", 0) # Standard deviation
            # For uniform distribution
            # Writer.Variable_("self.Range", 0)  # Mean +,- Range/2
            Writer.BlankLine()

            Writer.Statement("super().__init__()")
            Writer.BlankLine()

        with Writer.Statement("def LoadRateMean(self, Mean):"):
            Writer.Reshape__("self.Mean", "Mean", [])
            Writer.Cast_____("self.Mean", "self.Mean", 'float32')
            Writer.BlankLine()

        with Writer.Statement("def LoadRateSD(self, SD):"):
            Writer.Reshape__("self.SD", "SD", [])
            Writer.Cast_____("self.SD", "self.SD", 'float32')
            Writer.BlankLine()

        # with Writer.Statement("def LoadRateMinMax(self, Range):"):
        #     Writer.Reshape__("self.Range", "Range", [])
        #     Writer.BlankLine()

        # with Writer.Statement("def DetermineRateUniform(self):"):
        #     Writer.Statement("Min = Mean - Range/2")
        #     Writer.Statement("Max = Mean + Range/2")
        #     Writer.Statement("self.Rate = tf.random.uniform(shape=[1], minval=self.Min, maxval=self.Max, dtype=tf.int32)")
        #     Writer.ReturnVar("self.Rate")
        #     Writer.BlankLine()

        with Writer.Statement("def DetermineRate(self):"):
            Writer.Statement("self.Rate = tf.random.normal(shape=[1], mean=self.Mean, stddev=self.SD)")
            Writer.OperRound("self.Rate", "self.Rate")
            Writer.ReturnVar("self.Rate")
            Writer.BlankLine()
