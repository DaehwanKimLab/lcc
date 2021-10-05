import numpy as np

'''

Ce

'''

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):

    # Chromosome indices
    Idx_Ch_Original = [0]   # will stay as 1

    Idx_DNAStrand_Proxy_Left = [0]
    Idx_DNAStrand_Proxy_Right = [2]


    # Tasks to do (temporary hard coding solution)

    # Conditioning
    # Chromosomes Len == Max

    # Chromosome counts
    # New Original = Old Round 1 (unchanged count)
    # New Round 1 = Reset to 0
    # New Round 2 = Reset to 0

    # Chromosome Len to 0

    # Gene Counts / 2

    # Promoter Counts (not implemented yet)

    # RNA, Protein, Complex, Metabolite Counts / 2 (floordiv)

    # RNA, Protein Len
    # generate a mosaic random binary table


    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.CellDivisionSwitch", 0)
            Writer.BlankLine()

        with Writer.Statement("def ExecuteCellDivision(self):"):
            Writer.Statement("self.CellDivisionSwitch = self.EvaluateCellStateForCellDivision()")
            Writer.Statement("self.Cel.Len_ChromosomesReplicating = self.GenomeSegregation(self.Cel.Len_ChromosomesReplicating)")

            CountMatrices = [
                "self.Cel.Counts",
                ]
            for CountMatrix in CountMatrices:
                Writer.Statement("%s = self.SplitCountMatrix(%s)" % (CountMatrix, CountMatrix))

            LengthMatrices = [
                "self.Cel.Len_RNAsNascent",
                "self.Cel.Len_ProteinsNascent"
                ]

            for LenMatrix in LengthMatrices:
                Writer.Statement("%s = self.SplitLengthMatrices(%s)" % (LenMatrix, LenMatrix))
            Writer.BlankLine()


        with Writer.Statement("def GenomeSegregation(self, Len_ChromosomesReplicating):"):
            Writer.Overwrite("InheretedChromosomesReplicating", "Len_ChromosomesReplicating[1,:]")
            Writer.Reshape__("InheretedChromosomesReplicating", "InheretedChromosomesReplicating", [1, 4])
            Writer.VarFill__("NegOnes", [2, 4], -1)
            Writer.Concat___("Len_ChromosomesReplicating_New", "InheretedChromosomesReplicating", "NegOnes")

            Writer.Comment__("Update replicating chromosome length matrix")
            Writer.Multiply_("ToSubtractToMakeLenZero", "Len_ChromosomesReplicating", "self.CellDivisionSwitch")
            Writer.Subtract_("Len_ChromosomesReplicating_Updated", "Len_ChromosomesReplicating", "ToSubtractToMakeLenZero")

            Writer.Multiply_("ToReplaceToUpdateLen", "Len_ChromosomesReplicating_New", "self.CellDivisionSwitch")
            Writer.Add______("Len_ChromosomesReplicating_Updated", "Len_ChromosomesReplicating_Updated", "ToReplaceToUpdateLen")
            Writer.ReturnVar("Len_ChromosomesReplicating_Updated")
            Writer.BlankLine()

        with Writer.Statement("def EvaluateCellStateForCellDivision(self):"):
            Writer.Statement("DNAReplicationState = self.Cel.Len_ChromosomesReplicating[0, :]")
            Writer.Equal____("Bool_CellReplicationState", "DNAReplicationState", "self.Cel.Len_ChromosomesReplicatingMax")
            Writer.ReduceAll("Bool_CellDivisionSwitch", "Bool_CellReplicationState")
            Writer.BoolToBin("Bin_CellDivisionSwitch", "Bool_CellDivisionSwitch")
            Writer.ReturnVar("Bin_CellDivisionSwitch")
            Writer.BlankLine()

        with Writer.Statement("def SplitMatrixInValues(self, Matrix):"):
            Writer.FloorDiv_("Matrix_Halved", "Matrix", 2)
            Writer.ReturnVar("Matrix_Halved")
            Writer.BlankLine()

        with Writer.Statement("def SplitCountMatrix(self, CountMatrix):"):
            Writer.Statement("Matrix_Halved = self.SplitMatrixInValues(CountMatrix)")
            Writer.Multiply_("Matrix_Halved_Switch", "Matrix_Halved", "self.CellDivisionSwitch")
            Writer.Subtract_("CountMatrix_Updated", "CountMatrix", "Matrix_Halved_Switch")
            Writer.ReturnVar("CountMatrix_Updated")
            Writer.BlankLine()

        with Writer.Statement("def SplitMatrixInElements(self, Matrix):"):
            Writer.Shape____("Shape_Matrix", "Matrix")
            Writer.RndNumUni("RandomElements", Shape="Shape_Matrix", MinVal="0", MaxVal="2", Type='int32')
            Writer.Multiply_("Matrix_RandomElements", "Matrix", "RandomElements")
            Writer.Comment__("Correct 0 to -1")
            Writer.Statement("Matrix_RandomElements_Corrected = self.ResetZerosToNegOnes(Matrix_RandomElements)")
            Writer.ReturnVar("Matrix_RandomElements_Corrected")
            Writer.BlankLine()

        with Writer.Statement("def SplitLengthMatrices(self, LenMatrix):"):
            Writer.Subtract_("CellDivisionCounterSwitch", 1, "self.CellDivisionSwitch")
            Writer.Statement("Matrix_RandElements = self.SplitMatrixInElements(LenMatrix)")
            Writer.Statement("LenMatrix_Updated = tf.math.add(tf.math.multiply(CellDivisionCounterSwitch, LenMatrix), tf.math.multiply(self.CellDivisionSwitch, Matrix_RandElements))")
            Writer.ReturnVar("LenMatrix_Updated")
            Writer.BlankLine()

        with Writer.Statement("def TestCellDivision(self):"):
            Writer.Comment__("Double all of the cell state counts.")
            Writer.Statement("self.Cel.Counts = tf.math.multiply(self.Cel.Counts, 2)")
            Writer.ScatNdUpd("self.Cel.Counts", "self.Cel.Counts", "[[0]]", "[1]")

            Writer.Comment__("Set DNA replication state close to the completion.")
            Writer.Statement("self.Cel.Len_ChromosomesReplicating = tf.constant([[2319000, 2319000, 2319000, 2319000], [ -1, -1, -1, -1], [-1, -1, -1, -1]])")

            # Writer.Comment__("Forcing CellDivision without any molecular manipulation")
            # Writer.Statement("self.CellDivisionSwitch = tf.constant(1)")
            Writer.BlankLine()

        with Writer.Statement("def PrintMessage(self):"):
            Writer.PrintStrg("===== Cell Division ===== ")
            Writer.PrintStVa("Cell Division",
                             "self.CellDivisionSwitch")
            Writer.BlankLine()


