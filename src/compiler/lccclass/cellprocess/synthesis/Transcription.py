# transcription
# use sigma factor


def Write_CellProcess(Writer, Comp, ProGen, ProcessID):
    # ProGen.GenerateCellProcess(Writer, ProcessID)

    # Molecule indices for Molecular IDs
    Idx_RNAP = Comp.Master.ID2Idx_Master[Comp.Complex.Name2ID_Complexes['RNA polymerase, core enzyme']]

    Idx_NTPs = ProGen.BuildingBlockIdxs('NTP')
    Idx_PPi = Comp.Master.ID2Idx_Master['PPI[c]']

    # Set up Idx to Idx system in cell.py

    # Temporary parameters
    Rate_RNAElongation = 60
    Rate_RNAPActive = 0.22  # 22% of Free RNAP is active
    Rate_RNAPActiveCanBind = 0.5  # 50% of Free RNAP binds to the promoter
    NMax_RNAPsPerGene = 10

    # Temporary references
    NUniq_Genes = Comp.Gene.NUniq_Genes
    NUniq_RNAs = Comp.RNA.NUniq_RNAs
    assert NUniq_Genes == NUniq_RNAs

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        # Override the abstract method
        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.Cel.Idx_RNAP", Idx_RNAP)
            Writer.Variable_("self.Cel.Idx_NTPs", Idx_NTPs)
            Writer.Variable_("self.Cel.Idx_PPi", Idx_PPi)
            Writer.BlankLine()

            Writer.Variable_("self.Cel.Rate_RNAPActive", Rate_RNAPActive)  # Not implemented yet
            Writer.Variable_("self.Cel.Rate_RNAPActiveCanBind", Rate_RNAPActiveCanBind)
            Writer.Variable_("self.Cel.Rate_RNAElongation", Rate_RNAElongation)
            Writer.Variable_("self.Cel.Rate_RNAElongation_Matrix", Rate_RNAElongation, [NUniq_Genes, NMax_RNAPsPerGene])
            Writer.BlankLine()

            Writer.InitZeros("self.Cel.Len_RNAsNascentInitial", [NUniq_RNAs, NMax_RNAPsPerGene], 'int32')
            Writer.VarRepeat("self.Cel.Len_RNAsNascentMax", "self.Cel.Len_RNAs", NMax_RNAPsPerGene)
            Writer.Reshape__("self.Cel.Len_RNAsNascentMax", "self.Cel.Len_RNAsNascentMax", [NUniq_RNAs, NMax_RNAPsPerGene])
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def ExecuteProcess(self):"):
            Writer.Statement("self.TranscriptionInitiation()")
            Writer.Statement("self.TranscriptElongation()")
            Writer.Statement("self.TranscriptionTermination()")
            Writer.BlankLine()

        with Writer.Statement("def SigmaFactorBinding(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def RetrieveTotalRNAPCount(self):"):
            Writer.Statement("self.Cel.Count_RNAP = tf.constant(self.Cel.Counts[%s])  # %s: Index for RNAP" % (Idx_RNAP, Idx_RNAP))
            Writer.BlankLine()

        with Writer.Statement("def CalculateActiveRNAP(self):"):
            Writer.Cast_____("Count_RNAP_Float", "self.Cel.Count_RNAP", 'float32')
            Writer.OperElMul("self.Cel.Count_RNAPActive", "Count_RNAP_Float", "self.Cel.Rate_RNAPActive")  # For summary, active RNAP
            Writer.RoundInt_("self.Cel.Count_RNAPActive", "self.Cel.Count_RNAPActive")
            Writer.BlankLine()

        with Writer.Statement("def CalculateActiveRNAPCanBind(self):"):
            # Apply the rate of active RNAP binding (hypothetical)
            Writer.Cast_____("Count_RNAPActive_Float", "self.Cel.Count_RNAPActive", 'float32')
            Writer.OperElMul("self.Cel.Count_RNAPActiveCanBind", "Count_RNAPActive_Float", "self.Cel.Rate_RNAPActiveCanBind")  # For summary, active RNAP portion that can bind to promoter
            Writer.RoundInt_("self.Cel.Count_RNAPActiveCanBind", "self.Cel.Count_RNAPActiveCanBind")
            Writer.BlankLine()

        with Writer.Statement("def DetermineRNAPWillBind(self):"):
            # RNAP to bind in this step
            Writer.OperElSub("self.Cel.Count_RNAPWillBind", "self.Cel.Count_RNAPActiveCanBind", "self.Cel.Count_RNAPBound")
            Writer.OperElGr_("Bool_RNAPWillBind", "self.Cel.Count_RNAPActiveCanBind", "self.Cel.Count_RNAPBound")
            Writer.BoolToBin("Bin_RNAPWillBind", "Bool_RNAPWillBind")
            Writer.OperElMul("self.Cel.Count_RNAPWillBind", "self.Cel.Count_RNAPWillBind", "Bin_RNAPWillBind")   # For Summary, active RNAP to bind to promoter this step
            Writer.BlankLine()

        with Writer.Statement("def DetermineRNAPToBind(self):"):
            Writer.Statement("self.RetrieveTotalRNAPCount()")
            Writer.Statement("self.CalculateActiveRNAP()")
            Writer.Statement("self.CalculateActiveRNAPCanBind()")
            Writer.Statement("self.DetermineRNAPWillBind()")
            Writer.BlankLine()

        with Writer.Statement("def SelectRNAsToTranscribe(self):"):
            # Replace with weighted random later
            Writer.RndIdxUni("self.Cel.Idx_RndRNAsNascent", "self.Cel.Count_RNAPWillBind", "self.Cel.Idx_Master_RNAsNascent")

            # Weighted random distribution later
            # Writer.OperGathr("self.Cel.Count_GeneCopies", "self.Cel.Counts", "self.Cel.Idx_Master_Genes")
            # Writer.RndIdxWgh("self.Cel.Idx_RndRNAsNascent", "self.Cel.Count_RNAP", "self.Cel.Count_GeneCopies")
            Writer.BlankLine()

        with Writer.Statement("def IncrementCountOfNascentRNAs(self):"):
            Writer.InitOnes_("OnesForRndRNAs", "self.Cel.Count_RNAPWillBind", 'int32')
            # Distribute RNAP to selected RNA Nascent
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_RndRNAsNascent, OnesForRndRNAs)")
            Writer.BlankLine()

        with Writer.Statement("def IncrementBoundRNAPs(self):"):
            Writer.OperElAdd("self.Cel.Count_RNAPBound", "self.Cel.Count_RNAPBound", "self.Cel.Count_RNAPWillBind")
            Writer.BlankLine()

        with Writer.Statement("def DetermineRNAPUnbound(self):"):
            Writer.OperElSub("self.Cel.Count_RNAPUnbound", "self.Cel.Count_RNAP", "self.Cel.Count_RNAPBound")
            Writer.BlankLine()

        with Writer.Statement("def DistributeRNAPToGenes(self):"):
            Writer.Statement("self.DetermineRNAPToBind()")
            Writer.Statement("self.SelectRNAsToTranscribe()")
            Writer.Statement("self.IncrementCountOfNascentRNAs()")
            Writer.Statement("self.IncrementBoundRNAPs()")
            Writer.Statement("self.DetermineRNAPUnbound()")
            Writer.BlankLine()

        with Writer.Statement("def TranscriptionInitiation(self):"):
            Writer.Statement("self.SigmaFactorBinding()  # Not implemented")
            Writer.Statement("self.DistributeRNAPToGenes()")
            Writer.BlankLine()

        with Writer.Statement("def GetNascentRNACounts(self):"):
            Writer.OperGathr("self.Cel.Count_RNAsNascent_Matrix", "self.Cel.Counts", "self.Cel.Idx_Master_RNAsNascent")
            Writer.OperRdSum("self.Cel.Count_RNAsNascentInitialTotal", "self.Cel.Count_RNAsNascent_Matrix")
            # Temporary code
            # Writer.InitOnes_("self.Cel.Count_RNAsNascent_Matrix", NUniq_RNAs, 'int32')
            Writer.BlankLine()

        with Writer.Statement("def DetermineNascentRNAIndices(self):"):
            # Get indices of Nascent RNA in length table by using tf.range? or repeat (Ones x repeat by Count then turn into boolean)
            # TO DO: use a single API for the routine to fill 1s from the left then zeros

            Writer.Variable_("RNAPArrayForAllTranscripts", 0)
            with Writer.Statement("for i, Count in enumerate(self.Cel.Count_RNAsNascent_Matrix):"):
                Writer.InitOnes_("RNAP_Present", "Count", 'int32')
                Writer.InitZeros("RNAP_Abscent", "%s - Count" % NMax_RNAPsPerGene, 'int32')
                Writer.OperCncat("RNAPArrayForCurrentTranscript", "RNAP_Present", "RNAP_Abscent")
                with Writer.Statement("if i == 0:"):
                    Writer.Statement("RNAPArrayForAllTranscripts = RNAPArrayForCurrentTranscript")
                    Writer.Statement("continue")
                Writer.OperCncat("RNAPArrayForAllTranscripts", "RNAPArrayForAllTranscripts", "RNAPArrayForCurrentTranscript")
            Writer.Reshape__("RNAPArrayForAllTranscripts", "RNAPArrayForAllTranscripts", [-1, NMax_RNAPsPerGene])

            Writer.OperElGr_("self.Cel.Bool_RNAsNascentElongating", "RNAPArrayForAllTranscripts", 0)
            Writer.Cast_____("self.Cel.Bin_RNAsNascentElongating", "self.Cel.Bool_RNAsNascentElongating", 'int32')
            Writer.OperRdSum("self.Cel.Count_RNAsNascentElongatingTotal", "self.Cel.Bin_RNAsNascentElongating")  # For summary, the number of all elongating RNAs with RNAP
            Writer.BlankLine()

        with Writer.Statement("def GetNascentRNAsIndicesInLengthMatrix(self):"):
            Writer.Statement("self.GetNascentRNACounts()")
            Writer.Statement("self.DetermineNascentRNAIndices()")
            Writer.BlankLine()

        with Writer.Statement("def ElongateNascentRNAsInLengthMatrix(self):"):
            # Rate Matrix * bin(boolean) + length table (overelongated not counted yet)
            Writer.OperElMul("self.Cel.Rate_RNAElongation_Matrix", "self.Cel.Bin_RNAsNascentElongating", "self.Cel.Rate_RNAElongation")
            Writer.OperElAdd("self.Cel.Len_RNAsNascentElongated", "self.Cel.Len_RNAsNascentInitial", "self.Cel.Rate_RNAElongation_Matrix")
            Writer.BlankLine()

        with Writer.Statement("def DetermineOverElongatedNascentRNAs(self):"):
            # compare with len MAX to get boolean matrix for overelongated (>= MAX)
            Writer.OperElGr_("self.Cel.Bool_RNAsNascentOverElongated", "self.Cel.Len_RNAsNascentElongated", "self.Cel.Len_RNAsNascentMax")
            Writer.BoolToBin("self.Cel.Bin_RNAsNascentOverElongated", "self.Cel.Bool_RNAsNascentOverElongated")
            Writer.OperRdSum("self.Cel.Count_RNAsNascentOverElongatedTotal", "self.Cel.Bin_RNAsNascentOverElongated")  # For summary, the number of completed/overelongated RNAs
            Writer.BlankLine()

        with Writer.Statement("def DetermineAmountOfOverElongatedNTPs(self):"):
            Writer.OperElMul("MaxLengthForOverElongated", "self.Cel.Bin_RNAsNascentOverElongated", "self.Cel.Len_RNAsNascentMax")
            Writer.OperElMul("LengthForOverElongatedOnly", "self.Cel.Bin_RNAsNascentOverElongated", "self.Cel.Len_RNAsNascentElongated")
            # Get the over-elongated length of nascent RNAs
            Writer.OperElSub("self.Cel.Len_RNAsNascentOverElongated", "LengthForOverElongatedOnly", "MaxLengthForOverElongated")
            Writer.OperRdSum("self.Cel.Count_NTPsOverElongated", "self.Cel.Len_RNAsNascentOverElongated")  # For summary, over consumped NTPs to be corrected
            Writer.BlankLine()

        with Writer.Statement("def CorrectOverElongationInLengthMatrix(self):"):
            # Update nascent RNA length by removing over elongated length
            Writer.OperElSub("self.Cel.Len_RNAsNascentAdjusted", "self.Cel.Len_RNAsNascentElongated", "self.Cel.Len_RNAsNascentOverElongated")
            # now all completed elongation has its max value

            Writer.BlankLine()

        with Writer.Statement("def CorrectElongatedLengthForRNAs(self):"):
            # Initial elongation rate
            Writer.OperElMul("Rate_RNAElongation_Matrix_Initial", "self.Cel.Rate_RNAElongation_Matrix", "self.Cel.Bin_RNAsNascentElongating")
            # Adjusted elongation rate (initial - elongated)
            Writer.OperElSub("self.Cel.Rate_RNAElongation_Matrix_Corrected", "Rate_RNAElongation_Matrix_Initial", "self.Cel.Len_RNAsNascentOverElongated")
            # Calculate the length of elongation per gene
            Writer.OperRdSum("self.Cel.Count_RNAElongationLengthPerRNA", "self.Cel.Rate_RNAElongation_Matrix_Corrected", 1)  # For summary, total elongation length per Gene
            # Calculate the total length of elongation
            Writer.OperRdSum("self.Cel.Count_RNAElongationLengthTotal", "self.Cel.Rate_RNAElongation_Matrix_Corrected")  # For summary, total elongation length in this sim step
            Writer.BlankLine()

        with Writer.Statement("def ElongateTranscripts(self):"):
            # Nascent RNAs will be elongated.
            Writer.Statement("self.GetNascentRNAsIndicesInLengthMatrix()")
            Writer.Statement("self.ElongateNascentRNAsInLengthMatrix()")
            Writer.Statement("self.DetermineOverElongatedNascentRNAs()")
            Writer.Statement("self.DetermineAmountOfOverElongatedNTPs()")
            Writer.Statement("self.CorrectOverElongationInLengthMatrix()")
            Writer.Statement("self.CorrectElongatedLengthForRNAs()")
            Writer.BlankLine()

        with Writer.Statement("def GetRawNTPConsumption(self):"):
            # Prepare NT elongation length for each gene matrix for multiplication
            Writer.Reshape__("self.Cel.Count_RNAElongationLengthPerRNA", "self.Cel.Count_RNAElongationLengthPerRNA", [-1, 1])
            Writer.Cast_____("self.Cel.Count_RNAElongationLengthPerRNA", "self.Cel.Count_RNAElongationLengthPerRNA", 'float32')
            # NTP Frequency per gene * elongating transcript (in NT count) per gene
            Writer.OperMXMul("NTPConsumption_Raw", "self.Cel.Freq_NTsInRNAs", "self.Cel.Count_RNAElongationLengthPerRNA")
            Writer.ReturnVar("NTPConsumption_Raw")
            Writer.BlankLine()

        with Writer.Statement("def GetRoundedNTPConsumption(self, NTPConsumption_Raw):"):
            Writer.RoundInt_("NTPConsumption_Rounded", "NTPConsumption_Raw")
            Writer.ReturnVar("NTPConsumption_Rounded")
            Writer.BlankLine()

        with Writer.Statement("def CalculateDiscrepancy(self, NTPConsumption_Rounded):"):
            # Determine the difference between rounded vs corrected elongation
            Writer.OperRdSum("Sum_Rounded", "NTPConsumption_Rounded")
            Writer.Cast_____("Sum_Rounded", "Sum_Rounded", 'int32')
            Writer.OperRdSum("Sum_CorrectedElongationRate", "self.Cel.Rate_RNAElongation_Matrix_Corrected")
            Writer.OperElSub("DeltaSum", "Sum_CorrectedElongationRate", "Sum_Rounded")
            Writer.ReturnVar("DeltaSum")
            Writer.BlankLine()

        with Writer.Statement("def AdjustNTPConsumption(self, NTPConsumption_Rounded, Discrepancy):"):
            # Get equal amount of NTPs if discrepancy is greater than or equal to 4 or less than or eqqal to -4.
            Writer.OperElQuo("N_NTPSets", "Discrepancy", "self.Cel.NUniq_NTPs")
            Writer.VarRepeat("N_NTPSets", "N_NTPSets", "self.Cel.NUniq_NTPs")
            Writer.OperElAdd("NTPConsumption_MissingSet", "NTPConsumption_Rounded", "N_NTPSets")
            Writer.BlankLine()

            # Get random NTP for the remainder (replace with weighted random NTP based on ACGU in all genes)
            Writer.OperElRem("N_NTPRemainder", "Discrepancy", "self.Cel.NUniq_NTPs")
            Writer.Reshape__("N_NTPRemainder", "N_NTPRemainder", -1)
            Writer.InitZeros("NTPConsumption_MissingRemainder", "self.Cel.NUniq_NTPs", 'int32')
            Writer.RndNumUni("Idx_Remainder", "N_NTPRemainder", "0", "self.Cel.NUniq_NTPs")
            Writer.InitOnes_("OnesForRemainder", "N_NTPRemainder", 'int32')
            Writer.OperScAdd("NTPConsumption_MissingRemainder", "Idx_Remainder", "OnesForRemainder")
            Writer.BlankLine()

            # Return the adjusted NTP Consumption
            Writer.OperElAdd("NTPConsumption_Adjusted", "NTPConsumption_MissingSet", "NTPConsumption_MissingRemainder")
            Writer.OperRdSum("TotalNTPConsumption", "NTPConsumption_Adjusted")
            Writer.AsrtElEq_("TotalNTPConsumption", "self.Cel.Count_RNAElongationLengthTotal")
            Writer.ReturnVar("NTPConsumption_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def DetermineNTPConsumption(self):"):
            Writer.Statement("NTPConsumption_Raw = self.GetRawNTPConsumption()")
            Writer.Statement("NTPConsumption_Rounded = self.GetRoundedNTPConsumption(NTPConsumption_Raw)")
            Writer.Statement("Discrepancy = self.CalculateDiscrepancy(NTPConsumption_Rounded)")
            Writer.Statement("NTPConsumption_Adjusted = self.AdjustNTPConsumption(NTPConsumption_Rounded, Discrepancy)")
            Writer.ReturnVar("NTPConsumption_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def ConsumeNTPs(self):"):
            Writer.Statement("self.Cel.Count_RNAElongationNTPConsumption = self.DetermineNTPConsumption()")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_NTPs, -self.Cel.Count_RNAElongationNTPConsumption)")
            Writer.BlankLine()

        with Writer.Statement("def ReleasePPi(self):"):
            Writer.Overwrite("self.Cel.Count_RNAElongationPPiProduction", "self.Cel.Count_RNAElongationLengthTotal")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_PPi, self.Cel.Count_RNAElongationPPiProduction)")
            Writer.BlankLine()

        with Writer.Statement("def TranscriptElongation(self):"):
            Writer.Statement("self.ElongateTranscripts()")
            Writer.Statement("self.ConsumeNTPs()")
            Writer.Statement("self.ReleasePPi()")
            Writer.BlankLine()

        with Writer.Statement("def IdentifyCompletedRNAElongation(self):"):
            # identify all NMAX values by comparison (boolean)
            Writer.OperElEq_("self.Cel.Bool_RNAsNascentElongationCompleted", "self.Cel.Len_RNAsNascentAdjusted", "self.Cel.Len_RNAsNascentMax")
            Writer.BoolToBin("self.Cel.Bin_RNAsNascentElongationCompleted", "self.Cel.Bool_RNAsNascentElongationCompleted")
            Writer.BlankLine()

        with Writer.Statement("def CountCompletedRNAElongation(self):"):
            # Count all completed per RNA and total
            Writer.OperRdSum("self.Cel.Count_RNAElongationCompletedPerRNA", "self.Cel.Bin_RNAsNascentElongationCompleted", 1)
            Writer.OperRdSum("self.Cel.Count_RNAElongationCompletedTotal", "self.Cel.Bin_RNAsNascentElongationCompleted")   # For summary, total completed elongation
            Writer.BlankLine()

        with Writer.Statement("def ResetLengthOfCompletedNascentRNAs(self):"):
            # Reset all NMAX lengths to zero by subtracting max length from completed
            Writer.OperElMul("self.Cel.Len_RNAsNascentCompleted", "self.Cel.Len_RNAsNascentMax", "self.Cel.Bin_RNAsNascentElongationCompleted")
            Writer.OperElSub("self.Cel.Len_RNAsNascentFinal", "self.Cel.Len_RNAsNascentAdjusted", "self.Cel.Len_RNAsNascentCompleted")

            # Check reset status
            Writer.OperElMul("CheckCompletionReset", "self.Cel.Len_RNAsNascentFinal", "self.Cel.Bin_RNAsNascentElongationCompleted")
            Writer.AsrtElEq_("CheckCompletionReset", 0)
            Writer.BlankLine()

        with Writer.Statement("def GetIndicesOfRNAsCompletedElongation(self):"):
            # Get indices of RNAs completed elongation. Currently not used
            Writer.GetIdxGr_("Idx_RNAElongationCompleted", "self.Cel.Count_RNAElongationCompletedPerRNA", 0)
            Writer.OperGathr("self.Cel.Idx_RNAElongationCompleted", "self.Cel.Idx_RNAsNascent", "Idx_RNAElongationCompleted")
            Writer.BlankLine()

        with Writer.Statement("def ReleaseRNAP(self):"):
            Writer.Overwrite("self.Cel.Count_RNAPReleased", "self.Cel.Count_RNAElongationCompletedTotal")
            Writer.OperElSub("self.Cel.Count_RNAPBound", "self.Cel.Count_RNAPBound", "self.Cel.Count_RNAPReleased")
            Writer.OperElAdd("self.Cel.Count_RNAPUnbound", "self.Cel.Count_RNAPUnbound", "self.Cel.Count_RNAPReleased")
            Writer.OperElAdd("SumOfBoundUnbound", "self.Cel.Count_RNAPBound", "self.Cel.Count_RNAPUnbound")
            Writer.AsrtElEq_("self.Cel.Count_RNAP", "SumOfBoundUnbound")
            Writer.BlankLine()

        with Writer.Statement("def DeductCountOfNascentRNAs(self):"):
            # Deduct count of nascent RNAs completed elongation
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_Master_RNAsNascent, -self.Cel.Count_RNAElongationCompletedPerRNA)")
            Writer.BlankLine()

        with Writer.Statement("def IncrementCountOfRNAs(self):"):
            # Increment count of RNAs completed elongation
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_Master_RNAs, self.Cel.Count_RNAElongationCompletedPerRNA)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateLengthOfNascentRNAsMatrix(self):"):
            # Overwrite self.Cel.Len_RNAsNascentInitial with self.Cel.Len_RNAsNascentFinal
            Writer.Overwrite("self.Cel.Len_RNAsNascentInitial", "self.Cel.Len_RNAsNascentFinal")
            Writer.BlankLine()

        with Writer.Statement("def TranscriptionTermination(self):"):
            Writer.Statement("self.IdentifyCompletedRNAElongation()")
            Writer.Statement("self.CountCompletedRNAElongation()")
            Writer.Statement("self.ResetLengthOfCompletedNascentRNAs()")
            # Writer.Statement("self.GetIndicesOfRNAsCompletedElongation()")
            Writer.Statement("self.ReleaseRNAP()")
            Writer.Statement("self.DeductCountOfNascentRNAs()")
            Writer.Statement("self.IncrementCountOfRNAs()")
            Writer.Statement("self.UpdateLengthOfNascentRNAsMatrix()")
            Writer.BlankLine()
            
        with Writer.Statement("def Debugging_RNA_Indices(self, Str):"):
            Writer.Statement("print(Str)")
            Writer.BlankLine()

            Writer.Comment__("RNA length Indices")
            Writer.PrintStrg("RNA length Indices")

            Writer.PrintVar_("self.Cel.Idx_RndRNAsNascent")

            Writer.PrtIdxGr_("self.Cel.Len_RNAsNascentInitial", 0)
            Writer.PrtIdxGr_("self.Cel.Len_RNAsNascentElongated", 0)
            Writer.PrtIdxGr_("self.Cel.Len_RNAsNascentOverElongated", 0)

            Writer.PrtIdxGr_("self.Cel.Rate_RNAElongation_Matrix_Corrected", 0)

            Writer.PrtIdxGr_("self.Cel.Bin_RNAsNascentElongating", 0)
            Writer.PrtIdxGr_("self.Cel.Bin_RNAsNascentOverElongated", 0)
            Writer.PrtIdxGr_("self.Cel.Bin_RNAsNascentElongationCompleted", 0)

            Writer.BlankLine()

        with Writer.Statement("def Debugging_Count_RNAsNascent(self, Str):"):
            Writer.Statement("print(Str)")
            Writer.BlankLine()
            
            Writer.Comment__("Nascent RNAs")
            Writer.PrintStrg("Nascent RNAs")
            
            Writer.OperRdSum("Reduced_Sum_Of_self_Cel_Count_RNAsNascent", "self.Cel.Count_RNAsNascent")
            Writer.PrintVaVa("Reduced_Sum_Of_self_Cel_Count_RNAsNascent")
            Writer.BlankLine()
            
            Writer.OperGathr("Counts_RNAsNascent", "self.Cel.Counts", "self.Cel.Idx_Master_RNAsNascent")
            Writer.OperRdSum("Reduced_Sum_Of_Counts_RNAsNascent", "Counts_RNAsNascent")
            Writer.PrintVaVa("Reduced_Sum_Of_Counts_RNAsNascent")
            Writer.BlankLine()
            
            Writer.OperGathr("DeltaCounts_RNAsNascent", "self.Cel.DeltaCounts", "self.Cel.Idx_Master_RNAsNascent")
            Writer.OperRdSum("Reduced_Sum_Of_DeltaCounts_RNAsNascent", "DeltaCounts_RNAsNascent")
            Writer.PrintVaVa("Reduced_Sum_Of_DeltaCounts_RNAsNascent")
            Writer.BlankLine()

            Writer.Comment__("RNAs (full)")
            Writer.PrintStrg("RNAs (full)")

            Writer.OperRdSum("Reduced_Sum_Of_self_Cel_Count_RNAs", "self.Cel.Count_RNAs")
            Writer.PrintVaVa("Reduced_Sum_Of_self_Cel_Count_RNAs")
            Writer.BlankLine()

            Writer.OperGathr("Counts_RNAs", "self.Cel.Counts", "self.Cel.Idx_Master_RNAs")
            Writer.OperRdSum("Reduced_Sum_Of_Counts_RNAs", "Counts_RNAs")
            Writer.PrintVaVa("Reduced_Sum_Of_Counts_RNAs")
            Writer.BlankLine()

            Writer.OperGathr("DeltaCounts_RNAs", "self.Cel.DeltaCounts", "self.Cel.Idx_Master_RNAs")
            Writer.OperRdSum("Reduced_Sum_Of_DeltaCounts_RNAs", "DeltaCounts_RNAs")
            Writer.PrintVaVa("Reduced_Sum_Of_DeltaCounts_RNAs")
            Writer.BlankLine()

        with Writer.Statement("def ViewProcessSummary(self):"):

            Writer.PrintStrg("===== Transcription Initiation ===== ")
            # Number of Total RNAP
            Writer.PrintStVa("# of Total RNAPs",
                             "self.Cel.Count_RNAP")
            # Number of Active RNAP
            Writer.PrintStVa("# of active RNAPs",
                             "self.Cel.Count_RNAPActive")
            # Number of RNAP binding
            Writer.PrintStVa("# of RNAPs available that can bind a promoter",
                             "self.Cel.Count_RNAPActiveCanBind")
            # Number of new RNAP binding in this step (= new nascent RNAs to elongate next step)
            Writer.PrintStVa("# of RNAPs that newly bind a promoter in this step",
                             "self.Cel.Count_RNAPWillBind")
            Writer.BlankLine()

            Writer.PrintStrg("===== Transcript Elongation ===== ")
            # Number of Nascent RNAs elongated
            Writer.PrintStVa("# of all nascent RNAs elongating",
                             "self.Cel.Count_RNAsNascentElongatingTotal")
            # Total elongation length of RNAs
            Writer.PrintStVa("Total elongation length of RNAs (nt)",
                             "self.Cel.Count_RNAElongationLengthTotal")
            # Total NTP consumption and PPi production
            Writer.PrintStVa("Total NTP consumption [ATP, CTP, GTP, UTP]",
                             "self.Cel.Count_RNAElongationNTPConsumption")
            Writer.PrintStVa("Total PPi production",
                             "self.Cel.Count_RNAElongationPPiProduction")
            Writer.BlankLine()

            Writer.PrintStrg("===== Transcription Termination ===== ")
            # Number of RNA elongation Completed
            Writer.PrintStVa("# of RNA Elongation Completed",
                             "self.Cel.Count_RNAElongationCompletedTotal")
            # Number of RNAPs released
            Writer.PrintStVa("# of RNAPs released",
                             "self.Cel.Count_RNAPReleased")
            # Number of RNAP bound to DNA
            Writer.PrintStVa("# of total RNAPs bound to DNA",
                             "self.Cel.Count_RNAPBound")
            # Number of RNAP freely floating
            Writer.PrintStVa("# of total RNAPs unbound floating",
                             "self.Cel.Count_RNAPUnbound")
            Writer.BlankLine()
