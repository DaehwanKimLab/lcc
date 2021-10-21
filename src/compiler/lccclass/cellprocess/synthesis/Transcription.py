# transcription
# use sigma factor
import os
import datetime

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):

    # Molecule indices for Molecular IDs
    Idx_RNAP = Comp.Master.ID2Idx_Master[Comp.Complex.Name2ID_Complexes['RNA polymerase, core enzyme']]

    # Set up Idx to Idx system in cell.py

    # Temporary parameters
    Rate_RNAElongation = 60
    Rate_RNAP_Active = 0.22  # 22% of Free RNAP is active

    NMax_RNAPsPerGene = 10

    # Temporary references
    NUniq_Genes = Comp.Gene.NUniq_Genes
    NUniq_RNAs = Comp.RNA.NUniq_RNAs
    assert NUniq_Genes == NUniq_RNAs

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Function_("Init_ProcessSpecificVariables"):
            Writer.Variable_("self.Idx_RNAs", 0)
            Writer.Variable_("self.Idx_RndRNAsNascent", 0)
            Writer.Variable_("self.Idx_RNAP_Active", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Rate_RNAP_Active", 0)
            Writer.Variable_("self.Rate_RNAElongation", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_ElongatedTotal", 0)
            Writer.Variable_("self.Count_ElongationCompletionTotal", 0)
            Writer.Variable_("self.Count_NTPConsumptionTotal", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Len_RNAsNascentMax", 0)
            # For temporary export
            Writer.Variable_("self.Count_NTPConsumption", 0)
            Writer.BlankLine()
            
        # Override the abstract method
        with Writer.Function_("SetUp_ProcessSpecificVariables"):
            Writer.VarRange_("self.Idx_RNAs", 0, NUniq_RNAs)
            Writer.BlankLine()

            Writer.Variable_("self.Rate_RNAP_Active", Rate_RNAP_Active)
            Writer.Variable_("self.Rate_RNAElongation", Rate_RNAElongation)
            # TODO: Transcriptional efficiency may be imported from the literature (see flat data)

            Writer.VarFill__("self.Cel.Len_RNAsNascent", [NUniq_RNAs, NMax_RNAPsPerGene], -1)
            Writer.Reshape__("self.Len_RNAsNascentMax", "self.Cel.Len_RNAs", [-1, 1])
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Function_("ExecuteProcess"):
            Writer.Statement("self.Initiation()")
            Writer.Statement("self.Elongation()")
            Writer.Statement("self.Termination()")
            Writer.BlankLine()

        with Writer.Function_("DetermineRNAPToBindPromoter"):
            Writer.Comment__("Total count of active RNAP to bind")
            Writer.Statement("Count_RNAP_CoreEnzyme = self.GetCounts(self.Cel.Idx_RNAP_CoreEnzyme)")
            Writer.Statement("Count_RNAP_HoloEnzymes = self.GetCounts(self.Cel.Idx_RNAP_HoloEnzymes)")
            Writer.ReduceSum("Count_RNAP_HoloEnzymes_Total", "Count_RNAP_HoloEnzymes")
            Writer.Add______("Count_RNAP", "Count_RNAP_CoreEnzyme", "Count_RNAP_HoloEnzymes_Total")
            Writer.Cast_____("Count_RNAP_Float", "Count_RNAP", 'float32')
            Writer.Multiply_("Count_RNAP_ToBindPromoter_Total", "Count_RNAP_Float", "self.Rate_RNAP_Active")
            Writer.RoundInt_("Count_RNAP_ToBindPromoter_Total", "Count_RNAP_ToBindPromoter_Total")

            Writer.Comment__("Subtract DNA-bound RNAP to bind")
            Writer.ConvToBin("Bin_Len_RNAsNascent", "self.Cel.Len_RNAsNascent", ">=", 0)
            Writer.ReduceSum("Count_RNAP_BoundToDNA", "Bin_Len_RNAsNascent")
            Writer.Subtract_("Count_RNAP_ToBindPromoter", "Count_RNAP_ToBindPromoter_Total", "Count_RNAP_BoundToDNA")
            Writer.ReturnVar("Count_RNAP_ToBindPromoter")
            Writer.BlankLine()

        with Writer.Function_("SelectRNAsToTranscribe", "Count_RNAP_ToBindPromoter", "Freq_TranscriptionPerGene"):
            Writer.Comment__("Apply Gene Counts in transcription frequency")
            Writer.Statement("Count_Genes = self.GetCounts_Float(self.Cel.Idx_Master_Genes)")
            Writer.Reshape__("Count_Genes", "Count_Genes", [-1, 1])
            Writer.Multiply_("Freq_TranscriptionPerGene_WithCopyNumber", "Freq_TranscriptionPerGene", "Count_Genes")
            Writer.BlankLine()

            Writer.Comment__("Apply RNAP available to bind promoters")
            Writer.ReshType_("Count_RNAP_ToBindPromoter_Float", "Count_RNAP_ToBindPromoter", [-1, 1], 'float32')
            Writer.MatrixMul("Weight_TranscriptionPerGene", "Freq_TranscriptionPerGene_WithCopyNumber", "Count_RNAP_ToBindPromoter_Float")
            Writer.BlankLine()

            Writer.Comment__("Transform Transcription Weight")
            Writer.Statement("Weight_TranscriptionPerGene_Squeezed = self.SqueezeDistributionRangeZeroAndOne(Weight_TranscriptionPerGene)")
            Writer.Statement("Weight_TranscriptionPerGene_SqueezedWithZerosToOnes = self.ResetZerosToSpecificValue_Float(Weight_TranscriptionPerGene_Squeezed, 1)")
            Writer.Reshape__("Weight_TranscriptionPerGene_SqueezedWithZerosToOnes", "Weight_TranscriptionPerGene_SqueezedWithZerosToOnes", [-1, 1])
            Writer.ReduceMul("Weight_TranscriptionsPossible", "Weight_TranscriptionPerGene_SqueezedWithZerosToOnes", 1)
            Writer.Multiply_("Weight_TranscriptionsPossible_Stretched", "Weight_TranscriptionsPossible", "self.WeightResolution")
            Writer.RoundInt_("Weight_TranscriptionsPossible_StretchedInt", "Weight_TranscriptionsPossible_Stretched")
            Writer.BlankLine()

            Writer.Comment__("Give a minimum weight to the reactions rounded to 0 by adding 1 to all")
            Writer.Add______("Weight_TranscriptionsPossible_Adjusted", "Weight_TranscriptionsPossible_StretchedInt", 1)
            Writer.BlankLine()

            Writer.ReduceSum("N_RNAP_ToBindPromoter", "Count_RNAP_ToBindPromoter")
            Writer.Reshape__("N_RNAP_ToBindPromoter", "N_RNAP_ToBindPromoter", -1)
            Writer.Statement("Idx_RndRNAsNascent = self.PickRandomIndexFromPool_Weighted_Local(N_RNAP_ToBindPromoter, self.Idx_RNAs, Weight_TranscriptionsPossible_Adjusted)")
            Writer.ReturnVar("Idx_RndRNAsNascent")
            Writer.BlankLine()

        with Writer.Function_("AddNewCountsToNascentRNALengthMatrix", "Len_Matrix", "Idx_NewSelected"):
            # To turn -1 to 0 where nascent RNAs are to be made in the nascent RNA length table.
            Writer.InitZeros("Count_RNAsNascentNew", NUniq_RNAs, 'int32')
            Writer.OnesLike_("Ones", "Idx_NewSelected", 'int32')
            Writer.ScatNdAdd("Count_RNAsNascentNew", "Count_RNAsNascentNew", "Idx_NewSelected", "Ones")

            Writer.Statement("Bin_AddNewCount = self.GetBinToAddNewCountToLenMatrix(Len_Matrix, Count_RNAsNascentNew)")
            Writer.Add______("Len_Matrix_NewCountAdded", "Len_Matrix", "Bin_AddNewCount")
            Writer.ReturnVar("Len_Matrix_NewCountAdded")
            Writer.BlankLine()

        with Writer.Function_("UpdateNascentRNALengths", "Len_RNAsNascent_Update"):
            Writer.Overwrite("self.Cel.Len_RNAsNascent", "Len_RNAsNascent_Update")
            Writer.BlankLine()

        with Writer.Function_("ReleaseInitiationFactors", "Count_RNAP_Active_ToBindPromoter"):
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_RNAP_HoloEnzymes, -Count_RNAP_Active_ToBindPromoter)")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_SigmaFactors, Count_RNAP_Active_ToBindPromoter)")
            Writer.BlankLine()

        with Writer.Function_("Initiation"):
            Writer.Comment__("Distribute active RNAPs to genes")
            if 'SigmaFactor' in ProGen.Dict_CellProcesses:
                Writer.Statement("Count_RNAP_ToBindPromoters = self.GetCounts(self.Cel.Idx_RNAP_HoloEnzymes)")
            else:
                Writer.Statement("Count_RNAP_ToBindPromoters = self.DetermineRNAPToBindPromoter()")

            Writer.Statement("self.Idx_RndRNAsNascent = self.SelectRNAsToTranscribe(Count_RNAP_ToBindPromoters, self.Cel.Freq_TranscriptionPerGene)")
            Writer.BlankLine()

            Writer.Comment__("Update the nascent RNA length matrix")
            Writer.Statement("Len_RNANascent_NewToZero = self.AddNewCountsToNascentRNALengthMatrix(self.Cel.Len_RNAsNascent, self.Idx_RndRNAsNascent)")
            Writer.Statement("self.UpdateNascentRNALengths(Len_RNANascent_NewToZero)")
            Writer.BlankLine()

            if 'SigmaFactor' in ProGen.Dict_CellProcesses:
                Writer.Comment__("Update initiation factor counts")
                Writer.Statement("self.ReleaseInitiationFactors(Count_RNAP_ToBindPromoters)")
                Writer.BlankLine()
            else:
                pass

        with Writer.Function_("DetermineNTPConsumption", "Len_ToElongate"):
            Writer.Statement("NTPConsumption = self.DetermineAmountOfBuildingBlocks(Len_ToElongate, self.Cel.Freq_NTsInRNAs)")

            # Return the adjusted NTP Consumption
            Writer.ReduceSum("TotalNTPConsumptionFinal", "NTPConsumption")
            Writer.ReduceSum("TotalNTPConsumptionInitial", "Len_ToElongate")
            Writer.AsrtEq___("TotalNTPConsumptionFinal", "TotalNTPConsumptionInitial")
            Writer.ReturnVar("NTPConsumption")
            Writer.BlankLine()

        with Writer.Function_("ConsumeNTPs", "Len_ToElongate"):
            Writer.Statement("NTPConsumption = self.DetermineNTPConsumption(Len_ToElongate)")
            Writer.Statement("self.Count_NTPConsumption = NTPConsumption")
            Writer.ReduceSum("self.Count_NTPConsumptionTotal", "NTPConsumption")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_NTPs, -NTPConsumption)")
            Writer.BlankLine()

        with Writer.Function_("ReleasePPi", "Count_TotalLengthOfElongation"):
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_PPi, Count_TotalLengthOfElongation)")
            Writer.BlankLine()

        with Writer.Function_("UpdateByproducts", "Len_ToElongate"):
            Writer.ReduceSum("Count_TotalLengthOfElongationPerRNA", "Len_ToElongate", 1)
            Writer.Statement("self.ConsumeNTPs(Count_TotalLengthOfElongationPerRNA)")
            Writer.ReduceSum("Count_TotalLengthOfElongationTotal", "Len_ToElongate")
            Writer.Statement("self.ReleasePPi(Count_TotalLengthOfElongationTotal)")
            Writer.BlankLine()

        with Writer.Function_("Elongation"):
            Writer.Statement("Len_RNAsNascent_Elongated = self.DetermineAmountOfElongation(self.Cel.Len_RNAsNascent, self.Rate_RNAElongation, self.Len_RNAsNascentMax)")
            Writer.Statement("Len_ToElongate = Len_RNAsNascent_Elongated - self.Cel.Len_RNAsNascent")
            Writer.ConvToBin("Bin_Elongated", "Len_ToElongate", ">", 0)
            Writer.ReduceSum("self.Count_ElongatedTotal", "Bin_Elongated")
            Writer.Statement("self.UpdateNascentRNALengths(Len_RNAsNascent_Elongated)")
            Writer.Statement("self.UpdateByproducts(Len_ToElongate)")
            Writer.BlankLine()

        # with Writer.Statement("def GetIndicesOfRNAsCompletedElongation(self):"):
        #     # Get indices of RNAs completed elongation. Currently not used
        #     Writer.GetIdx___("Idx_RNAElongationCompleted", "self.Cel.Count_RNAElongationCompletedPerRNA", ">", 0)
        #     Writer.Gather___("self.Cel.Idx_RNAElongationCompleted", "self.Cel.Idx_RNAsNascent",
        #                      "Idx_RNAElongationCompleted")
        #     Writer.BlankLine()

        with Writer.Function_("IncrementCountOfRNAs", "Count_ElongationCompletionPerRNA"):
            # Increment count of RNAs completed elongation
            Writer.Statement(
                "self.AddToDeltaCounts(self.Cel.Idx_Master_RNAs, Count_ElongationCompletionPerRNA)")
            Writer.BlankLine()

        with Writer.Function_("Termination"):
            # Reset nascent RNA length matrix
            Writer.Statement("Len_RNAsNascent_ElongationCompletedToZero = self.ResetLengthOfCompletedElongation(self.Cel.Len_RNAsNascent, self.Len_RNAsNascentMax)")
            Writer.Statement("Len_RNAsNascent_ElongationCompletedToNegOne = self.ResetZerosToNegOnes(Len_RNAsNascent_ElongationCompletedToZero)")
            Writer.Statement("self.UpdateNascentRNALengths(Len_RNAsNascent_ElongationCompletedToNegOne)")
            Writer.BlankLine()

            # Update RNA and RNAP core counts
            Writer.ConvToBin("Bin_ElongationCompleted", "Len_RNAsNascent_ElongationCompletedToZero", "==", 0)
            Writer.ReduceSum("Count_ElongationCompletionPerRNA", "Bin_ElongationCompleted", 1)
            Writer.ReduceSum("self.Count_ElongationCompletionTotal", "Count_ElongationCompletionPerRNA")
            Writer.Statement("self.IncrementCountOfRNAs(Count_ElongationCompletionPerRNA)")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_RNAP_CoreEnzyme, self.Count_ElongationCompletionTotal)")
            Writer.BlankLine()

        with Writer.Function_("ViewProcessSummary"):
            Writer.PrintStrg("===== Transcription ===== ")
            Writer.PrintStVa("# of RNAPs",
                             "self.GetCounts(self.Idx_RNAP_Active)")
            Writer.PrintStVa("# of RNAPs Newly Bound To Gene",
                             "tf.shape(self.Idx_RndRNAsNascent)[0]")
            Writer.PrintStVa("# of Nascent RNAs Elongated (= # of Active RNAPs)",
                             "self.Count_ElongatedTotal")
            Writer.PrintStVa("# of New RNAs Generated",
                             "self.Count_ElongationCompletionTotal")
            Writer.PrintStVa("# of NTP consumption",
                             "self.Count_NTPConsumptionTotal")
            Writer.BlankLine()
