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
            Writer.ScatNdAdd("self.Cel.DeltaCounts", "MolIdxs", "MolCounts")
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.PrintStrg("Summary message not implemented yet")
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def ViewProcessDebuggingMessages(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def DebuggingVar(self, VarBefore, VarAfter):"):
            Writer.Statement("print('[Debug] %s' % VarAfter)")
            Writer.Subtract_("Delta", "VarBefore", "VarAfter")
            Writer.PrintVaVa("Delta")
            Writer.BlankLine()

        with Writer.Statement("def GetBinToAddNewCountToLenMatrix(self, LengthMatrix, CountMatrixToAdd):"):
            Writer.Comment__("Generate boolean and binary matrices for element availability in nascent proteins length matrix.")
            Writer.Less_____("Bool_LenAvailable", "LengthMatrix", 0)
            Writer.BoolToBin("Bin_LenAvailable", "Bool_LenAvailable")
            Writer.BlankLine()

            Writer.Comment__("Generate a cumulative sum matrix for labeling the available elements with count.")
            Writer.CumSum___("LenCumsum", "Bin_LenAvailable", 1)
            Writer.BlankLine()

            Writer.Comment__("Generate a boolean matrix of the count label matrix > 0, < count.")
            Writer.Greater__("Bool_LenCumsumGreaterThanZero", "LenCumsum", 0)
            Writer.Reshape__("CountMatrixToAdd", "CountMatrixToAdd", [-1, 1])
            Writer.LessEq___("Bool_LenCumsumLessThanOrEqualToCountMatrixToAdd", "LenCumsum", "CountMatrixToAdd")
            Writer.LogicAnd_("Bool_LenCumsum", "Bool_LenCumsumGreaterThanZero", "Bool_LenCumsumLessThanOrEqualToCountMatrixToAdd")
            Writer.BlankLine()

            Writer.Comment__("Generate a binary matrix for elements satisfying both availability and count conditions.")
            Writer.LogicAnd_("Bool_LenSelected", "Bool_LenAvailable", "Bool_LenCumsum")
            Writer.BoolToBin("Bin_LenSelected", "Bool_LenSelected")
            Writer.ReturnVar("Bin_LenSelected")
            Writer.BlankLine()

        with Writer.Statement("def ConvertCountToBinary(self, Count, N_Columns):"):
            Writer.Comment__("Generate binary pair vector")
            Writer.Statement("Length_Count = tf.shape(Count)")
            Writer.Multiply_("Length_BinaryPairArray", "Length_Count", "2")
            Writer.VarTile__("BinaryPairArray", "[1, 0]", "Length_BinaryPairArray")
            Writer.BlankLine()

            Writer.Comment__("Generate binary count vector")
            Writer.Reshape__("Count_2D", "Count", "[1, -1]")
            Writer.Subtract_("Count_Inv", "N_Columns", "Count_2D")
            Writer.Concat___("BinaryPairCount", "Count", "Count_Inv", "0")
            Writer.Transpose("BinaryPairCount")
            Writer.Reshape__("BinaryPairCount", "BinaryPairCount", "-1")
            Writer.BlankLine()

            Writer.Comment__("Generate Generate binary matrix")
            Writer.VarRepeat("BinaryRepeatArray", "BinaryPairArray", "BinaryPairCount")
            Writer.Reshape__("BinaryMatrix", "BinaryRepeatArray", "[-1, N_Columns]")
            Writer.Statement("return BinaryMatrix")
            Writer.BlankLine()

        with Writer.Statement("def PickRandomIndexFromPool_Uniform(self, N_MoleculesToDistribute, Indices):"):
            Writer.Comment__("Select random values")
            Writer.RndNumUni("RandomValues", Shape="N_MoleculesToDistribute", MaxVal="len(Indices)", Type='int32')
            Writer.Comment__("Select indices corresponding to the random values")
            Writer.Gather___("Idx_Random", "Indices", "RandomValues")
            Writer.ReturnVar("Idx_Random")
            Writer.BlankLine()

        with Writer.Statement("def PickRandomIndexFromPool_Weighted(self, N_MoleculesToDistribute, Indices, Weights):"):
            Writer.Reshape__("Weights", "Weights", -1)
            Writer.Reshape__("Indices", "Indices", -1)

            # Generate a scalar value for the shape of Weights
            Writer.Statement("N_RandomValues = tf.shape(Indices)[0]")

            # Cumsum on Counts
            Writer.CumSum___("WeightCumSum", "Weights", 0)

            # Generate random number ranging between first and last elements of cumsum int values
            Writer.RndNumUni("RandomValues", Shape="N_MoleculesToDistribute", MinVal="WeightCumSum[0]", MaxVal="WeightCumSum[-1]", Type='int32')

            # Repeat random number array then reshape to generate matrix shaped with [Count shape, -1]
            Writer.VarRepeat("RandomValues_RepeatArray", "RandomValues", "N_RandomValues")
            Writer.Reshape__("RandomValues_Matrix", "RandomValues_RepeatArray", "[-1, N_RandomValues]")

            # Boolean operation on (repeated random number array < Cumsum (2D matrix)) for which count range it falls
            Writer.Less_____("Bool_RandomValuesLessThanWeightCumSum", "RandomValues_Matrix", "WeightCumSum")

            # Boolean operation on (Counts > 0) for idxs with counts
            Writer.Greater__("Bool_RandomValuesGreaterThanZero", "RandomValues_Matrix", 0)

            # Logical_and on (repeated random number array < Cumsum, Counts > 0) then convert to binary
            Writer.LogicAnd_("Bool_RandomValuesSelected", "Bool_RandomValuesLessThanWeightCumSum", "Bool_RandomValuesGreaterThanZero")
            Writer.BoolToBin("Bin_RandomValuesSelected", "Bool_RandomValuesSelected")

            # Argmax on the binary to get the first indices of 1's in each random number generated.
            Writer.ArgMax___("Idx_Random", "Bin_RandomValuesSelected", Axis=1)

            Writer.ReturnVar("Idx_Random")
            Writer.BlankLine()

        with Writer.Statement("def DetermineAmountOfElongation(self, Len_Matrix, Rate_Elongation, MaxLen_Matrix):"):
            Writer.Comment__("Elongate Length Matrix")
            # Get a binary matrix to elongate
            Writer.ConvToBin("Bin_LengthElongating", "Len_Matrix", ">=", 0)

            # Apply elongation rate to the elongating elements (overelongated not counted yet)
            Writer.Multiply_("Rate_Elongation_Matrix", "Bin_LengthElongating",
                             "Rate_Elongation")
            Writer.Add______("Len_Matrix_Elongated", "Len_Matrix",
                             "Rate_Elongation_Matrix")
            Writer.BlankLine()

            Writer.Comment__("Correct Over-elongated Elements")
            # Compare to len MAX to determine overelongated (> MAX)
            Writer.Reshape__("MaxLen_2D", "MaxLen_Matrix", [-1, 1])
            Writer.ConvToBin("Bin_GreaterThanZero", "Len_Matrix_Elongated", ">", "MaxLen_2D")
            Writer.Subtract_("Len_Matrix_Elongated_MaxSubtracted", "Len_Matrix_Elongated", "MaxLen_2D")
            Writer.Multiply_("Len_Matrix_OverElongated", "Len_Matrix_Elongated_MaxSubtracted", "Bin_GreaterThanZero")

            # Update nascent Protein length by removing over elongated length
            Writer.Subtract_("Len_Matrix_Elongated_Corrected", "Len_Matrix_Elongated",
                             "Len_Matrix_OverElongated")
            # now all completed elongation has its max value
            Writer.BlankLine()

            Writer.ReturnVar("Len_Matrix_Elongated_Corrected")
            Writer.BlankLine()

        with Writer.Statement("def ResetLengthOfCompletedElongation(self, Len_Matrix, MaxLen_Matrix):"):
            # Reset all NMAX lengths to -1 by subtracting max length from completed
            Writer.ConvToBin("Bin_Len_MatrixMax", "Len_Matrix", "==", "MaxLen_Matrix")
            Writer.Multiply_("MaxLen_Matrix_MaxElementsOnly", "MaxLen_Matrix", "Bin_Len_MatrixMax")
            Writer.Subtract_("Len_Matrix_MaxToZero", "Len_Matrix", "MaxLen_Matrix_MaxElementsOnly")
            Writer.BlankLine()

            # Check len matrix integrity
            Writer.AsrtGrEq_("Len_Matrix_MaxToZero", -1)
            Writer.BlankLine()

            Writer.ReturnVar("Len_Matrix_MaxToZero")
            Writer.BlankLine()

        with Writer.Statement("def ResetZerosToNegOnes(self, Len_Matrix):"):
            # Replace 0's with -1 to indicate no RNAP binding
            Writer.ConvToBin("Ones", "Len_Matrix", "==", 0)
            Writer.Subtract_("Len_Matrix_ZeroToNegOne", "Len_Matrix", "Ones")
            Writer.BlankLine()

            # Check len matrix integrity
            Writer.AsrtNoEq_("Len_Matrix_ZeroToNegOne", 0)
            Writer.AsrtGrEq_("Len_Matrix_ZeroToNegOne", -1)
            Writer.BlankLine()

            Writer.ReturnVar("Len_Matrix_ZeroToNegOne")
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
