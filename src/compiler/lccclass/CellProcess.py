# Interface for all lcc modules, a.k.a. cellular processes


# Comp is a short hand for CompilerData
def Write_CellProcess(Writer):
    with Writer.Statement("class FCellProcess():"):
        with Writer.Statement("def __init__(self, Cel, Cst, Env, Exe):"):
            Writer.LinkClObj('Cel')
            Writer.LinkClObj('Cst')
            Writer.LinkClObj('Env')

            Writer.LinkClObj('Exe')
            Writer.BlankLine()

            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def ExecuteProcess(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def AddToDeltaCounts(self, MolIdxs, MolCounts):"):
            Writer.OperScAdd("self.Cel.DeltaCounts", "MolIdxs", "MolCounts")
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.Pass_____()
            Writer.BlankLine()


        # Writer.AbsMethod()
        # with Writer.Statement("def AddToStoichiometryMatrix(self):"):
        #     Writer.Pass_____()
        #     Writer.BlankLine()
        #
        # Writer.AbsMethod()
        # with Writer.Statement("def CalculateRate(self):"):
        #     Writer.Pass_____()
        #     Writer.BlankLine()
        #
        # Writer.AbsMethod()
        # with Writer.Statement("def AddToRateMatrix(self):"):
        #     Writer.Pass_____()
        #     Writer.BlankLine()
        #
        # Writer.AbsMethod()
        # with Writer.Statement("def PrintMolCounts(self):"):
        #     Writer.Pass_____()
        #     Writer.BlankLine()
        #
        # Writer.TF_Graph_()
        # with Writer.Statement("def SetUpStoichiometryMatrix(self):"):
        #     Writer.Statement("self.AddToStoichiometryMatrix()")
        #     Writer.BlankLine()
        #
        # Writer.TF_Graph_()
        # with Writer.Statement("def UpdateRates(self):"):
        #     Writer.Statement("self.CalculateRate()")
        #     Writer.Statement("self.AddToRateMatrix()")
        #     Writer.BlankLine()
