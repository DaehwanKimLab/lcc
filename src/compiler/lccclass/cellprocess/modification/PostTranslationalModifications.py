# posttranslational modifications
# protein localization

def Write_ModPRT(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FModPRT():"):
        with Writer.Function_("__init__"):
            Writer.BlankLine()
            Writer.Statement("super().__init__()")

            Writer.BlankLine()

        # Abstract Methods for CellProcess
        Writer.AbsMethod()
        with Writer.Function_("InitProcess"):
            # Call AddElementaryProcess
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Function_("LoopProcess"):
            # Call UpdateReactionRate
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Function_("AddElementaryProcess"):
            # Call GetReactionMolIndex, GetReactionStoich, GetReactionRate methods
            # Call AddToMasterReactionStoichs, AddToMasterReactionRates
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Function_("GetReactionMolIndex"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Function_("GetReactionStoich"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Function_("GetReactionRate"):
            # AdjustRate
            Writer.Pass_____()
            Writer.BlankLine()