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


def Write_RGM_MetaboliteReplenish(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FRGM_MetaboliteReplenish(FReactionExecution):"):
        with Writer.Statement("def __init__(self):"):

            Writer.Variable_("self.MetaboliteIdxs", 0)
            Writer.Variable_("self.MetaboliteCountsInitial", 0)

            Writer.Statement("super().__init__()")
            Writer.BlankLine()

        with Writer.Statement("def LoadMetaboliteIdxs(self, Idxs):"):
            Writer.Overwrite("self.MetaboliteIdxs", "Idxs")
            Writer.BlankLine()

        with Writer.Statement("def LoadMetaboliteCountsInitial(self, Counts):"):
            Writer.Overwrite("self.MetaboliteCountsInitial", "Counts")
            Writer.BlankLine()

        with Writer.Statement("def ReplenishMetabolites(self, FinalCount):"):
            Writer.ScatNdUpd("FinalCount", "FinalCount", "self.MetaboliteIdxs", "self.MetaboliteCountsInitial")
            Writer.Reshape__("FinalCount_Replenished", "FinalCount", [-1, 1])
            Writer.ReturnVar("FinalCount_Replenished")
            Writer.BlankLine()

        with Writer.Statement("def AddCountMatrices(self):"):
            Writer.Add______("FinalCount", "self.Count", "self.DeltaCount")
            Writer.ReturnVar("FinalCount")
            Writer.BlankLine()