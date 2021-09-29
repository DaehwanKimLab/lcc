
# Comp is a short hand for CompilerData
def Write_Visulization_2D(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FVisualization_2D():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Statement("# Define temporary variables for visualization purposes")
            Writer.BlankLine()

        with Writer.Statement("def VisualizeData(self, X, Y, XLabel, YLabel, Legends, Title)")



