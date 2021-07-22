# Metabolism

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):
    # ProGen.GenerateCellProcess(Writer, ProcessID)

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Statement("def Init_ProcessSpecificVariables()"):
            Writer.Variable_("self.MetaboliteIdxs", 0)
            Writer.Variable_("self.MetaboliteCountsInitial", 0)
            Writer.Statement("self.GetMetaboliteIdxsAndCounts()")
            Writer.BlankLine()

        with Writer.Statement("def GetMetaboliteIdxsAndCounts(self):"):
            MetaboliteIdxs = list()
            MetaboliteCountsInitial = list()
            for MolIdx in range(Comp.Master.NUniq_Master):
                if Comp.Master.Type_Master[MolIdx] == 'Metabolite':
                    MetaboliteIdxs.append(MolIdx)
                    MetaboliteCountsInitial.append(Comp.Master.Count_Master[MolIdx])

            Writer.Variable_("self.MetaboliteIdxs", MetaboliteIdxs)
            Writer.Reshape__("self.MetaboliteIdxs", "self.MetaboliteIdxs", [1, -1])
            Writer.Statement("self.Exe.LoadMetaboliteIdxs(self.MetaboliteIdxs)")

            Writer.Variable_("self.MetaboliteCountsInitial", MetaboliteCountsInitial)
            Writer.Reshape__("self.MetaboliteCountsInitial", "self.MetaboliteCountsInitial", [-1])
            Writer.Statement("self.Exe.LoadMetaboliteCountsInitial(self.MetaboliteCountsInitial)")
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

