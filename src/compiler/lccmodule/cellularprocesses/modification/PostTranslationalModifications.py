# posttranslational modifications
# protein localization

def Write_ModPRT(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FModPRT():"):
        with Writer.Statement("def __init__(self):"):
            Writer.BlankLine()
            Writer.Statement("super().__init__()")

            Writer.BlankLine()

        # Abstract Methods for CellProcess
        Writer.AbsMethod()
        with Writer.Statement("def InitProcess(self, Cel, Cst, Env):"):
            # Call AddElementaryProcess
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def LoopProcess(self, Cel, Cst, Env, Sim):"):
            # Call UpdateReactionRate
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def AddElementaryProcess(self, Cel, Cst, Env):"):
            # Call GetReactionMolIndex, GetReactionStoich, GetReactionRate methods
            # Call AddToMasterReactionStoichs, AddToMasterReactionRates
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def GetReactionMolIndex(self, Cel):"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def GetReactionStoich(self, Cel):"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def GetReactionRate(self, Cel):"):
            # AdjustRate
            Writer.Pass_____()
            Writer.BlankLine()