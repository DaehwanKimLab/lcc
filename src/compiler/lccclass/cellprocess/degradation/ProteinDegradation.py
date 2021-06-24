# protein damage
# ssrA-mediated protein degradation by ClpXP and ClpAP proteases

def Write_DegPRT(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FDegPRT(FCellProcess):"):
        with Writer.Statement("def __init__(self):"):
            Writer.BlankLine()
            Writer.Statement("super().__init__()")

            Writer.BlankLine()

        # Abstract Methods for CellProcess
        Writer.AbsMethod()
        with Writer.Statement("def InitProcess(self):"):
            # Call AddElementaryProcess
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def LoopProcess(self):"):
            # Call UpdateReactionRate
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def AddElementaryProcess(self):"):
            # Call GetReactionMolIndex, GetReactionStoich, GetReactionRate methods
            # Call AddToMasterReactionStoichs, AddToMasterReactionRates
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def GetReactionMolIndex(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def GetReactionStoich(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def GetReactionRate(self):"):
            # AdjustRate
            Writer.Pass_____()
            Writer.BlankLine()