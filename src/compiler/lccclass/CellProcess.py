# Interface for all lcc modules, a.k.a. cellular processes


# Comp is a short hand for CompilerData
def Write_CellProcess(Writer):
    with Writer.Statement("class FCellProcess():"):
        with Writer.Function_("__init__", "Cel", "Cst", "Env", "Exe"):
            Writer.LinkClObj('Cel')
            Writer.LinkClObj('Cst')
            Writer.LinkClObj('Env')

            Writer.LinkClObj('Exe')
            Writer.BlankLine()

            Writer.InitOnes_("self.One", 1, 'int32')
            Writer.InitOnes_("self.One_Float", 1, 'float32')
            Writer.InitZeros("self.Zero", 1, 'int32')
            Writer.InitZeros("self.Zero_Float", 1, 'float32')
            Writer.Variable_("self.WeightResolution", 2 ** 10, 'tf.float32')   # For random number generator
            Writer.BlankLine()

        # Abstract methods
        Writer.AbsMethod()
        with Writer.Function_("Init_ProcessSpecificVariables"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Function_("SetUp_ProcessSpecificVariables"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Function_("ExecuteProcess"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Function_("ViewProcessSummary"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Function_("ViewProcessDebuggingMessages"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Function_("SaveData"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Function_("Pass"):
            Writer.Pass_____()
            Writer.BlankLine()


        # Get/Add Counts/DeltaCounts
        with Writer.Function_("GetCounts", "MolIdxs"):
            Writer.ReturnVar("self.Cel.GetCounts(MolIdxs)")
            Writer.BlankLine()

        with Writer.Function_("GetCounts_Float", "MolIdxs"):
            Writer.Statement("Count = self.Cel.GetCounts(MolIdxs)")
            Writer.Cast_____("Count_Float", "Count", 'float32')
            Writer.ReturnVar("Count_Float")
            Writer.BlankLine()

        with Writer.Function_("GetCounts_Rate", "MolIdxs", "Rate"):
            Writer.Statement("Count = self.Cel.GetCounts(MolIdxs)")
            Writer.Statement("Count_RateApplied = self.DetermineAmountFromRate(Count, Rate)")
            Writer.ReturnVar("Count_RateApplied")
            Writer.BlankLine()

        with Writer.Function_("GetDeltaCounts", "MolIdxs"):
            Writer.ReturnVar("self.Cel.GetDeltaCounts(MolIdxs)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts", "MolIdxs", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(MolIdxs, MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_RndIdx", "MolIdxs"):
            Writer.OnesLike_("MolCounts", "MolIdxs", 'int32')
            Writer.Statement("self.Cel.AddToDeltaCounts(MolIdxs, MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_dNTPs", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(self.Cel.Idx_dNTPs, MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_NTPs", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(self.Cel.Idx_NTPs, MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_AAs", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(self.Cel.Idx_AAs, MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_DehydrationSynthesis", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(self.Cel.Idx_H2O, MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_Hydrolysis", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(self.Cel.Idx_H2O, -MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_PiRelease", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(self.Cel.Idx_Pi, -MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_Phosphorolysis", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(self.Cel.Idx_Pi, -MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_ATPConsumption", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(self.Cel.Idx_ATP, -MolCounts)")
            Writer.BlankLine()

        with Writer.Function_("AddToDeltaCounts_PPiRelease", "MolCounts"):
            Writer.Statement("self.Cel.AddToDeltaCounts(self.Cel.Idx_PPi, MolCounts)")
            Writer.BlankLine()

        # Conversion
        with Writer.Function_("ConvertConcToCount", "Conc"):
            Writer.Multiply_("Mole", "Conc", "self.Cel.Vol")
            Writer.Multiply_("Count", "Mole", "self.Cst.NA")
            Writer.Cast_____("Count", "Count", 'int32')
            Writer.ReturnVar("Count")
            Writer.BlankLine()

        with Writer.Function_("GetConc", "MolIdxs"):
            Writer.Statement("Count = self.Cel.GetCounts(MolIdxs)")
            Writer.Cast_____("Count", "Count", 'float32')
            Writer.Divide___("Mole", "Count", "self.Cst.NA")
            Writer.Divide___("Conc", "Mole", "self.Cel.Vol")
            Writer.ReturnVar("Conc")
            Writer.BlankLine()

        # with Writer.Statement("def AddToDeltaCounts_Hydrolysis(self, MolCounts):"):
        #     Writer.ScatNdAdd("self.Cel.DeltaCounts", "self.Cel.DeltaCounts", "self.Cel.Idx_H2O", "MolCounts")
        #     Writer.BlankLine()

        with Writer.Function_("DebuggingVar", "VarBefore", "VarAfter"):
            Writer.Statement("print('[Debug] %s' % VarAfter)")
            Writer.Subtract_("Delta", "VarBefore", "VarAfter")
            Writer.PrintVaVa("Delta")
            Writer.BlankLine()

        # Determine Amounts
        with Writer.Function_("DetermineAmountOfBuildingBlocks", "Len_Polymers", "Freq_Monomers"):
            Writer.Statement("NUniq_Monomers = tf.shape(Freq_Monomers, out_type='int32')[0]")
            Writer.Statement("NUniq_Polymers = tf.shape(Len_Polymers, out_type='int32')[0]")
            Writer.BlankLine()

            Writer.Comment__("GetRawMonomerConsumption")
            # Prepare NT elongation length for each mRNA matrix for multiplication
            Writer.ReshType_("Len_Polymers_Float", "Len_Polymers", '[NUniq_Polymers, -1]', 'float32')
            # Monomer Frequency per mRNA * elongating Polypeptide (in NT count) per mRNA
            Writer.MatrixMul("MonomerConsumption_Raw", "Freq_Monomers",
                             "Len_Polymers_Float")
            Writer.BlankLine()

            Writer.Comment__("GetRoundedMonomerConsumption")
            Writer.RoundInt_("MonomerConsumption_Rounded", "MonomerConsumption_Raw")
            Writer.BlankLine()

            Writer.Comment__("CalculateDiscrepancy")
            # Determine the difference between rounded vs corrected elongation
            Writer.ReduceSum("LenSum_AfterRounding", "MonomerConsumption_Rounded")
            Writer.Cast_____("LenSum_AfterRounding", "LenSum_AfterRounding", 'int32')
            Writer.ReduceSum("LenSum_BeforeRounding", "Len_Polymers")
            Writer.Subtract_("Discrepancy", "LenSum_BeforeRounding", "LenSum_AfterRounding")
            Writer.BlankLine()

            Writer.Comment__("AdjustMonomerConsumption")
            # Get equal amount of Monomers if discrepancy is greater than or equal to 4 or less than or equal to -4.
            Writer.FloorDiv_("N_MonomerSets", "Discrepancy", "NUniq_Monomers")
            Writer.VarRepeat("N_MonomerSets", "N_MonomerSets", "NUniq_Monomers")
            Writer.Add______("MonomerConsumption_Rounded_MissingSet", "MonomerConsumption_Rounded", "N_MonomerSets")
            Writer.BlankLine()

            # Get random Monomer for the remainder (replace with weighted random Monomer based on Monomer frequency in all mRNAs)
            Writer.Remainder("N_MonomerRemainder", "Discrepancy", "NUniq_Monomers")
            Writer.Reshape__("N_MonomerRemainder", "N_MonomerRemainder", -1)
            Writer.InitZeros("MonomerConsumption_MissingRemainder", "NUniq_Monomers", 'int32')
            Writer.RndNumUni("Idx_Remainder", "N_MonomerRemainder", "0", "NUniq_Monomers")
            Writer.InitOnes_("OnesForRemainder", "N_MonomerRemainder", 'int32')
            Writer.ScatNdAdd("MonomerConsumption_MissingRemainder", "MonomerConsumption_MissingRemainder", "Idx_Remainder", "OnesForRemainder")
            Writer.BlankLine()

            # Calculate adjusted Monomer Consumption
            Writer.Add______("Count_Monomers", "MonomerConsumption_Rounded_MissingSet", "MonomerConsumption_MissingRemainder")
            Writer.ReturnVar("Count_Monomers")
            Writer.BlankLine()

        with Writer.Function_("GetBinToAddNewCountToLenMatrix", "LengthMatrix", "CountMatrixToAdd"):
            Writer.Comment__("Generate boolean and binary matrices for element availability in nascent proteins length matrix.")
            Writer.Less_____("Bool_LenAvailable", "LengthMatrix", 0)
            Writer.BoolToBin("Bin_LenAvailable", "Bool_LenAvailable")
            Writer.BlankLine()

        # Determine X number of 0 containing elements
        with Writer.Function_("GetBinToAddNewCountToLenMatrix", "LengthMatrix", "CountMatrixToAdd"):
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

        # Expand 1-dimentional values to 2-dimentional 1's
        with Writer.Function_("ConvertCountToBinary", "Count", "N_Columns"):
            Writer.Comment__("Generate binary pair vector")
            Writer.Statement("Length_Count = tf.shape(Count)")
            Writer.Multiply_("Length_BinaryPairArray", "Length_Count", "2")
            Writer.VarTile__("BinaryPairArray", "[1, 0]", "Length_BinaryPairArray")
            Writer.BlankLine()

            Writer.Comment__("Generate binary count vector")
            Writer.Reshape__("Count_2D", "Count", "[1, -1]")
            Writer.Subtract_("Count_Inv", "N_Columns", "Count_2D")
            Writer.Concat___("BinaryPairCount", "Count", "Count_Inv", "0")
            Writer.Transpose("BinaryPairCount", "BinaryPairCount")
            Writer.Reshape__("BinaryPairCount", "BinaryPairCount", "-1")
            Writer.BlankLine()

            Writer.Comment__("Generate Generate binary matrix")
            Writer.VarRepeat("BinaryRepeatArray", "BinaryPairArray", "BinaryPairCount")
            Writer.Reshape__("BinaryMatrix", "BinaryRepeatArray", "[-1, N_Columns]")
            Writer.ReturnVar("BinaryMatrix")
            Writer.BlankLine()

        with Writer.Function_("IdxFromLocalToReference", "Idx_Local", "Idx_Reference"):
            Writer.Gather___("Idx_Converted", "Idx_Reference", "Idx_Local")
            Writer.ReturnVar("Idx_Converted")
            Writer.BlankLine()

        with Writer.Function_("DetermineAmountOfElongation", "Len_Matrix", "Rate_Elongation", "MaxLen_Matrix"):
            Writer.Comment__("Elongate Length Matrix")
            Writer.Comment__("Get a binary matrix to elongate")
            Writer.ConvToBin("Bin_LengthElongating", "Len_Matrix", ">=", 0)

            Writer.Comment__("Apply elongation rate to the elongating elements (overelongated not counted yet)")
            Writer.Multiply_("Rate_Elongation_Matrix", "Bin_LengthElongating",
                             "Rate_Elongation")
            Writer.Add______("Len_Matrix_Elongated", "Len_Matrix",
                             "Rate_Elongation_Matrix")
            Writer.BlankLine()

            Writer.Comment__("Correct Over-elongated Elements")
            Writer.Comment__("Compare to len MAX to determine overelongated (> MAX)")
            # Writer.Reshape__("MaxLen_2D", "MaxLen_Matrix", [-1, 1])
            Writer.ConvToBin("Bin_GreaterThanZero", "Len_Matrix_Elongated", ">", "MaxLen_Matrix")
            Writer.Subtract_("Len_Matrix_Elongated_MaxSubtracted", "Len_Matrix_Elongated", "MaxLen_Matrix")
            Writer.Multiply_("Len_Matrix_OverElongated", "Len_Matrix_Elongated_MaxSubtracted", "Bin_GreaterThanZero")

            Writer.Comment__("Update nascent Protein length by removing over elongated length")
            Writer.Subtract_("Len_Matrix_Elongated_Corrected", "Len_Matrix_Elongated",
                             "Len_Matrix_OverElongated")

            Writer.Comment__("now all completed elongation has its max value")
            Writer.BlankLine()

            Writer.ReturnVar("Len_Matrix_Elongated_Corrected")
            Writer.BlankLine()

        with Writer.Function_("ResetLengthOfCompletedElongation", "Len_Matrix", "MaxLen_Matrix"):
            Writer.Comment__("Reset all NMAX lengths to -1 by subtracting max length from completed")
            Writer.ConvToBin("Bin_Len_MatrixMax", "Len_Matrix", "==", "MaxLen_Matrix")
            Writer.Multiply_("MaxLen_Matrix_MaxElementsOnly", "MaxLen_Matrix", "Bin_Len_MatrixMax")
            Writer.Subtract_("Len_Matrix_MaxToZero", "Len_Matrix", "MaxLen_Matrix_MaxElementsOnly")
            Writer.BlankLine()

            Writer.Comment__("Check len matrix integrity")
            Writer.AsrtGrEq_("Len_Matrix_MaxToZero", -1)
            Writer.BlankLine()

            Writer.ReturnVar("Len_Matrix_MaxToZero")
            Writer.BlankLine()
        # with Writer.Statement("def SqueezeDistributionRangeNegAndPosOne(self, Input):"):
        #     Writer.ReduceMax("Max", "Input")
        #     Writer.ReduceMin("Min", "Input")
        #     Writer.Add______("MaxMinAdded", "Max", "Min")
        #     Writer.Cast_____("MaxMinAdded", "MaxMinAdded", 'float32')
        #     Writer.Divide___("Median", "MaxMinAdded", 2)
        #     Writer.Subtract_("Input_ZeroCentered", "Input", "Median")
        #     Writer.ReduceMax("Max_ZeroCentered", "Input_ZeroCentered")
        #     Writer.Divide___("Input_ZeroCenteredNormalizedToOne", "Input_ZeroCentered", "Max_ZeroCentered")
        #     Writer.ReturnVar("Input_ZeroCenteredNormalizedToOne")
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def StretchDistributionRangeToNegAndPosValue(self, Input, Value):"):
        #     Writer.Statement("Input_NormalizedToOne = self.SqueezeDistributionRangeNegAndPosOne(Input)")
        #     Writer.Multiply_("Input_Stretched", "Input_NormalizedToOne", "Value")
        #     Writer.ReturnVar("Input_Stretched")
        #     Writer.BlankLine()


        # Matrix operation tools
        with Writer.Function_("DetermineAmountFromRate", "Count", "Rate"):
            Writer.Cast_____("Count_Float", "Count", 'float32')
            Writer.Multiply_("Count_RateApplied", "Count_Float", "Rate")
            Writer.RoundInt_("Count_RateApplied_Rounded", "Count_RateApplied")
            Writer.ReturnVar("Count_RateApplied_Rounded")
            Writer.BlankLine()

        with Writer.Function_("ConvertToRatio", "Count_Matrix"):
            Writer.ReduceSum("ReducedSum_Count", "Count_Matrix")
            Writer.Divide___("Freq_Count", "Count_Matrix", "ReducedSum_Count")
            Writer.Transpose("Freq_Count", "Freq_Count")
            Writer.ReturnVar("Freq_Count")
            Writer.BlankLine()

        with Writer.Function_("PickRandomIndexFromPool_Uniform", "N_MoleculesToDistribute", "Indices"):
            Writer.Comment__("Select random values")
            Writer.RndNumUni("RandomValues", Shape="N_MoleculesToDistribute", MaxVal="len(Indices)", Type='int32')
            Writer.Comment__("Select indices corresponding to the random values")
            Writer.Gather___("Idx_Random", "Indices", "RandomValues")
            Writer.ReturnVar("Idx_Random")
            Writer.BlankLine()

        with Writer.Function_("PickRandomIndexFromPool_Weighted_Local", "N_MoleculesToDistribute", "Indices", "Weights"):
            Writer.Reshape__("Weights", "Weights", -1)
            Writer.Reshape__("Indices", "Indices", -1)

            Writer.Comment__("Generate a scalar value for the shape of Weights")
            Writer.Statement("N_RandomValues = tf.shape(Indices)[0]")

            Writer.Comment__("Cumsum on Counts")
            Writer.CumSum___("WeightCumSum", "Weights", 0)

            Writer.Comment__("Generate random number ranging between first and last elements of cumsum int values")
            Writer.RndNumUni("RandomValues", Shape="N_MoleculesToDistribute", MinVal="WeightCumSum[0]", MaxVal="WeightCumSum[-1]", Type='int32')

            Writer.Comment__("Repeat random number array then reshape to generate matrix shaped with [Count shape, -1]")
            Writer.VarRepeat("RandomValues_RepeatArray", "RandomValues", "N_RandomValues")
            Writer.Reshape__("RandomValues_Matrix", "RandomValues_RepeatArray", "[-1, N_RandomValues]")

            Writer.ConvToBin("Bin_RandomValuesLessThanWeightCumSum", "RandomValues_Matrix", "<", "WeightCumSum")

            Writer.Comment__("Argmax on the binary to get the first indices of 1's in each random number generated.")
            Writer.ArgMax___("Idx_Random", "Bin_RandomValuesLessThanWeightCumSum", Axis=1)

            Writer.ReturnVar("Idx_Random")
            Writer.BlankLine()

        with Writer.Function_("PickRandomIndexFromPool_Weighted_Global", "N_MoleculesToDistribute", "Indices", "Weights"):
            Writer.Statement("Idx_Random_Local = self.PickRandomIndexFromPool_Weighted_Local(N_MoleculesToDistribute, Indices, Weights)")
            Writer.Statement("Idx_Random = self.IdxFromLocalToReference(Idx_Random_Local, Indices)")
            # Writer.Gather___("Idx_Random", "Indices", "Idx_Random_Local")
            Writer.ReturnVar("Idx_Random")
            Writer.BlankLine()


        with Writer.Function_("ResetZerosToNegOnes", "Len_Matrix"):
            Writer.Comment__("Replace 0's with -1 to indicate no processor binding")
            Writer.Statement("Len_Matrix_ZeroToNegOne = self.ResetZerosToSpecificValue(Len_Matrix, -1)")

            Writer.Comment__("Check len matrix integrity")
            Writer.AsrtNoEq_("Len_Matrix_ZeroToNegOne", 0)
            Writer.AsrtGrEq_("Len_Matrix_ZeroToNegOne", -1)
            Writer.BlankLine()

            Writer.ReturnVar("Len_Matrix_ZeroToNegOne")
            Writer.BlankLine()

        with Writer.Function_("ResetZerosToSpecificValue", "Matrix", "Value"):
            Writer.ConvToBin("Ones_ForZeros", "Matrix", "==", 0)
            Writer.Multiply_("Values_ForZeros", "Ones_ForZeros", "Value")
            Writer.Add______("Matrix_ZeroToValue", "Matrix", "Values_ForZeros")
            Writer.ReturnVar("Matrix_ZeroToValue")
            Writer.BlankLine()

        with Writer.Function_("ResetZerosToSpecificValue_Float", "Matrix", "Value"):
            Writer.ConvToBin("Ones_ForZeros", "Matrix", "==", 0)
            Writer.Multiply_("Values_ForZeros", "Ones_ForZeros", "Value")
            Writer.Cast_____("Matrix", "Matrix", 'float32')
            Writer.Cast_____("Values_ForZeros", "Values_ForZeros", 'float32')
            Writer.Add______("Matrix_ZeroToValue", "Matrix", "Values_ForZeros")
            Writer.ReturnVar("Matrix_ZeroToValue")
            Writer.BlankLine()

        with Writer.Function_("RemoveZeroElements", "Matrix"):
            Writer.NotEqual_("Bool_NonZeros", "Matrix", 0)
            Writer.BoolMask_("Matrix_ZerosRemoved", "Matrix", "Bool_NonZeros")
            Writer.ReturnVar("Matrix_ZerosRemoved")
            Writer.BlankLine()

        with Writer.Function_("GenerateCountMatrixForSelectedIndices", "Indices", "NUniq_Indices"):
            Writer.InitZeros("Zeros", "NUniq_Indices", 'int32')
            Writer.OnesLike_("Ones_Indices", "Indices", 'int32')
            Writer.ScatNdAdd("Count_Indices", "Zeros", "Indices", "Ones_Indices")
            Writer.ReturnVar("Count_Indices")
            Writer.BlankLine()

        with Writer.Function_("GetCountAfterApplyingDeltaCount", "Idx_Master"):
            Writer.Statement("Count = self.GetCounts(Idx_Master)")
            Writer.Statement("DeltaCount = self.GetDeltaCounts(Idx_Master)")
            Writer.Add______("Count_AfterApplyingDeltaCount", "Count", "DeltaCount")
            Writer.ReturnVar("Count_AfterApplyingDeltaCount")
            Writer.BlankLine()

        with Writer.Function_("IdentifyCountBelowZeroAfterRemoval", "Idx_Master", "Count_ToRemove"):
            Writer.Statement("Count_Current = self.GetCountAfterApplyingDeltaCount(Idx_Master)")
            Writer.Subtract_("Count_AfterRemoval", "Count_Current", "Count_ToRemove")
            Writer.ConvToBin("Bin_BelowZeroAfterRemoval", "Count_AfterRemoval", "<", 0)
            Writer.Multiply_("Count_BelowZeroAfterRemoval", "Count_AfterRemoval", "Bin_BelowZeroAfterRemoval")
            Writer.ReturnVar("-Count_BelowZeroAfterRemoval")
            Writer.BlankLine()

        with Writer.Function_("CorrectCountGettingBelowZeroAfterRemoval", "Idx_Master", "Count_ToRemove"):
            Writer.Statement("Count_ToRescue = self.IdentifyCountBelowZeroAfterRemoval(Idx_Master, Count_ToRemove)")
            Writer.Subtract_("Count_ToRemove_Corrected", "Count_ToRemove", "Count_ToRescue")
            Writer.ReturnVar("Count_ToRemove_Corrected")
            Writer.BlankLine()

        # CorrectCountGettingBelowZeroAfterRemoval Version 1
        # with Writer.Statement("def CorrectCountGettingBelowZeroAfterRemoval(self, Idx_Master, Idx_Removal, Count_Removal):"):
        #     Writer.Statement("Count_ToRescue = self.IdentifyCountBelowZeroAfterRemoval(Idx_Master, Count_Removal)")
        #     Writer.Subtract_("Count_Removal", "Count_Removal", "Count_ToRescue")
        #     Writer.Greater__("Bool_Count_Removal", "Count_Removal", 0)
        #     Writer.BoolMask_("Count_Removal_Adjusted", "Count_Removal", "Bool_Count_Removal")
        #     Writer.BoolMask_("Idx_RemovalAdjusted", "Idx_Removal", "Bool_Count_Removal")
        #     Writer.ReturnVar("Idx_RemovalAdjusted", "Count_Removal_Adjusted")
        #     Writer.BlankLine()

        with Writer.Function_("SqueezeDistributionRangeZeroAndOne", "Input"):
            Writer.AsrtGrEq_("Input", 0, "tf.where(Input < 0)")
            Writer.Cast_____("Input_Float", "Input", 'float32')
            Writer.ReduceMax("Max", "Input_Float")
            Writer.DivNoNan_("Input_Squeezed", "Input_Float", "Max")
            # Writer.Reshape__("Input_Squeezed", "Input_Squeezed", [1, -1])
            Writer.ReturnVar("Input_Squeezed")
            Writer.BlankLine()

        with Writer.Function_("DistributeCountsInRatio", "Counts_ToDistribute", "Count_TotalToAdjustTo"):
            Writer.Comment__("Prepare variables for float calculations")
            Writer.ReduceSum("Counts_ToDistribute_Total", "Counts_ToDistribute")

            Writer.Cast_____("Counts_ToDistribute_Total_Float", "Counts_ToDistribute_Total", 'float32')
            Writer.Cast_____("Count_TotalToAdjustTo_Float", "Count_TotalToAdjustTo", 'float32')
            Writer.Cast_____("Counts_ToDistribute_Float", "Counts_ToDistribute", 'float32')
            Writer.BlankLine()

            Writer.Comment__("Apply Ratio to counts")
            Writer.Divide___("RateToAdjust", "Count_TotalToAdjustTo_Float", "Counts_ToDistribute_Total_Float")
            Writer.Multiply_("Counts_Distributed_Float", "Counts_ToDistribute_Float", "RateToAdjust")
            Writer.RoundInt_("Counts_Distributed", "Counts_Distributed_Float")
            Writer.BlankLine()

            Writer.Comment__("Reality Check")
            Writer.ReduceSum("Counts_Distributed_Total", "Counts_Distributed")
            Writer.AsrtGrEq_("Count_TotalToAdjustTo", "Counts_Distributed_Total")
            Writer.BlankLine()

            Writer.ReturnVar("Counts_Distributed")
            Writer.BlankLine()

