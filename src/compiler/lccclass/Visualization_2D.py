
# Comp is a short hand for CompilerData
def Write_Visulization_2D(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FVisualization_2D():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Statement("# Define temporary variables for visualization purposes")
            Writer.Variable_("self.dNTPs", 0)
            Writer.Variable_("self.NTPs", 0)
            Writer.Variable_("self.AAs", 0)
            Writer.Variable_("self.Produced_RNAs", 0)
            Writer.Variable_("self.Degraded_RNAs", 0)
            Writer.Variable_("self.Produced_Proteins", 0)
            Writer.Variable_("self.Degraded_Proteins", 0)

            Writer.BlankLine()