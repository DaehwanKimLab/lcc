# Metabolism

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):
    # ProGen.GenerateCellProcess(Writer, ProcessID)

    Idx_NADH = Comp.Master.ID2Idx_Master['NADH[c]']
    Idx_NADPH = Comp.Master.ID2Idx_Master['NADPH[c]']

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.Count_MetabolitesInitial", 0)
            Writer.BlankLine()

        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.Cel.Idx_NADH", Idx_NADH)
            Writer.Variable_("self.Cel.Idx_NADPH", Idx_NADPH)
            Writer.BlankLine()

            Writer.Statement("self.GetMetaboliteCounts()")
            Writer.BlankLine()

        with Writer.Statement("def GetMetaboliteCounts(self):"):
            Writer.Gather___("self.Cel.Count_MetabolitesInitial", "self.Cel.Count_Master", "self.Cel.Idx_Master_Metabolites")
            Writer.BlankLine()

        with Writer.Statement("def ReplenishMetabolites(self):"):
            Writer.ScatNdUpd("self.Cel.Counts", "self.Cel.Counts", "self.Cel.Idx_Master_Metabolites", "self.Cel.Count_MetabolitesInitial")
            Writer.BlankLine()

        with Writer.Statement("def Message_ReplenishMetabolites(self):"):
            Writer.PrintStrg("\t- Metabolites have been replenished.")
            Writer.BlankLine()

        with Writer.Statement("def ExecuteProcess(self):"):
            Writer.Pass_____()
            Writer.BlankLine()



# def SetUpReactions(ProGen):
#     Reactions = list()
#
#     Reaction = None
#
#     Reaction_SetUp = Reaction
#     Reactions.append(Reaction_SetUp)
#
#     return Reactions
#

