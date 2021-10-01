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
    Rate_RNAPActive = 0.22  # 22% of Free RNAP is active
    Rate_RNAPActiveCanBind = 0.5  # 50% of Free RNAP binds to the promoter

    NMax_RNAPsPerGene = 10

    # Temporary references
    NUniq_Genes = Comp.Gene.NUniq_Genes
    NUniq_RNAs = Comp.RNA.NUniq_RNAs
    assert NUniq_Genes == NUniq_RNAs


    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.Idx_RNAs", 0)
            Writer.Variable_("self.Idx_RndRNAsNascent", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Rate_RNAPActive", 0)
            Writer.Variable_("self.Rate_RNAPActiveCanBind", 0)
            Writer.Variable_("self.Rate_RNAElongation", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_RNAP", 0)
            Writer.Variable_("self.Count_ElongatedTotal", 0)
            Writer.Variable_("self.Count_ElongationCompletionTotal", 0)
            Writer.Variable_("self.Count_NTPConsumptionTotal", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Len_RNAsNascentMax", 0)
            # For temporary export
            Writer.Variable_("self.Count_NTPConsumption", 0)
            Writer.BlankLine()
            
        # Override the abstract method
        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):
            Writer.VarRange_("self.Idx_RNAs", 0, NUniq_RNAs)
            
            Writer.Variable_("self.Cel.Idx_RNAP", Idx_RNAP)
            Writer.BlankLine()

            Writer.Variable_("self.Rate_RNAPActive", Rate_RNAPActive)
            Writer.Variable_("self.Rate_RNAPActiveCanBind", Rate_RNAPActiveCanBind)
            Writer.Variable_("self.Rate_RNAElongation", Rate_RNAElongation)
            # TODO: Transcriptional efficiency may be imported from the literature (see flat data)

            Writer.VarFill__("self.Cel.Len_RNAsNascent", [NUniq_RNAs, NMax_RNAPsPerGene], -1)
            Writer.Reshape__("self.Len_RNAsNascentMax", "self.Cel.Len_RNAs", [-1, 1])
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def ExecuteProcess(self):"):
            Writer.Statement("self.Initiation()")
            Writer.Statement("self.Elongation()")
            Writer.Statement("self.Termination()")
            Writer.BlankLine()

        with Writer.Statement("def SigmaFactorBinding(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def DetermineRNAPToBind(self):"):
            Writer.Comment__("RetrieveTotalRNAPCount")
            Writer.Statement(
                "Count_RNAP = tf.constant(self.Cel.Counts[%s])  # %s: Index for RNAP" % (
                Idx_RNAP, Idx_RNAP))
            Writer.Overwrite("self.Count_RNAP", "Count_RNAP")
            Writer.BlankLine()

            Writer.Comment__("CalculateActiveRNAP")
            Writer.Cast_____("Count_RNAP_Float", "Count_RNAP", 'float32')
            Writer.Multiply_("Count_RNAPActive", "Count_RNAP_Float",
                             "self.Rate_RNAPActive")  # For summary, active RNAP
            Writer.RoundInt_("Count_RNAPActive", "Count_RNAPActive")
            Writer.BlankLine()

            Writer.Comment__("CalculateActiveRNAPCanBind")
            Writer.Cast_____("Count_RNAPActive_Float", "Count_RNAPActive", 'float32')
            Writer.Multiply_("Count_RNAPActiveCanBind", "Count_RNAPActive_Float",
                             "self.Rate_RNAPActiveCanBind")  # For summary, active RNAP portion that can bind to promoter
            Writer.RoundInt_("Count_RNAPActiveCanBind", "Count_RNAPActiveCanBind")
            Writer.BlankLine()

            Writer.Comment__("DetermineRNAPWillBind")
            Writer.ConvToBin("Bin_RNAPBound", "self.Cel.Len_RNAsNascent", ">=", 0)
            Writer.ReduceSum("Count_RNAPBound", "Bin_RNAPBound")
            Writer.Subtract_("Count_RNAPWillBind", "Count_RNAPActiveCanBind",
                             "Count_RNAPBound")
            Writer.ConvToBin("Bin_RNAPWillBind", "Count_RNAPActiveCanBind", ">",
                             "Count_RNAPBound")
            Writer.Multiply_("Count_RNAPWillBind", "Count_RNAPWillBind",
                             "Bin_RNAPWillBind")  # For Summary, active RNAP to bind to promoter this step
            Writer.ReturnVar("Count_RNAPWillBind")
            Writer.BlankLine()

        with Writer.Statement("def SelectRNAsToTranscribe(self, Count_RNAPWillBind):"):
            # Weighted random distribution later
            Writer.Gather___("Count_Genes", "self.Cel.Counts", "self.Cel.Idx_Master_Genes")
            Writer.Statement("Idx_RndRNAsNascent = self.PickRandomIndexFromPool_Weighted_Local(Count_RNAPWillBind, self.Idx_RNAs, Count_Genes)")
            Writer.ReturnVar("Idx_RndRNAsNascent")
            Writer.BlankLine()

        with Writer.Statement("def AddNewCountsToNascentRNALengthMatrix(self, Len_Matrix, Idx_NewSelected):"):
            # To turn -1 to 0 where nascent RNAs are to be made in the nascent RNA length table.
            Writer.InitZeros("Count_RNAsNascentNew", NUniq_RNAs, 'int32')
            Writer.OnesLike_("Ones", "Idx_NewSelected", 'int32')
            Writer.ScatNdAdd("Count_RNAsNascentNew", "Count_RNAsNascentNew", "Idx_NewSelected", "Ones")

            Writer.Statement("Bin_AddNewCount = self.GetBinToAddNewCountToLenMatrix(Len_Matrix, Count_RNAsNascentNew)")
            Writer.Add______("Len_Matrix_NewCountAdded", "Len_Matrix", "Bin_AddNewCount")
            Writer.ReturnVar("Len_Matrix_NewCountAdded")
            Writer.BlankLine()

        with Writer.Statement("def DistributeRNAPsToGenes(self):"):
            Writer.Statement("Count_RNAPWillBind = self.DetermineRNAPToBind()")
            Writer.Statement("Idx_RndRNAsNascent = self.SelectRNAsToTranscribe(Count_RNAPWillBind)")
            Writer.ReturnVar("Idx_RndRNAsNascent")
            Writer.BlankLine()

        with Writer.Statement("def UpdateNascentRNALengths(self, Len_RNAsNascent_Update):"):
            Writer.Overwrite("self.Cel.Len_RNAsNascent", "Len_RNAsNascent_Update")
            Writer.BlankLine()

        with Writer.Statement("def Initiation(self):"):
            Writer.Statement("self.SigmaFactorBinding()  # Not implemented")
            Writer.Statement("self.Idx_RndRNAsNascent = self.DistributeRNAPsToGenes()")

            # Update the nascent RNA length matrix
            Writer.Statement("Len_RNANascent_NewToZero = self.AddNewCountsToNascentRNALengthMatrix(self.Cel.Len_RNAsNascent, self.Idx_RndRNAsNascent)")
            Writer.Statement("self.UpdateNascentRNALengths(Len_RNANascent_NewToZero)")
            Writer.BlankLine()

        with Writer.Statement("def DetermineNTPConsumption(self, Len_ToElongate):"):
            Writer.Statement("NTPConsumption = self.DetermineAmountOfBuildingBlocks(Len_ToElongate, self.Cel.Freq_NTsInRNAs)")

            # Return the adjusted NTP Consumption
            Writer.ReduceSum("TotalNTPConsumptionFinal", "NTPConsumption")
            Writer.ReduceSum("TotalNTPConsumptionInitial", "Len_ToElongate")
            Writer.AsrtEq___("TotalNTPConsumptionFinal", "TotalNTPConsumptionInitial")
            Writer.ReturnVar("NTPConsumption")
            Writer.BlankLine()

        with Writer.Statement("def ConsumeNTPs(self, Len_ToElongate):"):
            Writer.Statement("NTPConsumption = self.DetermineNTPConsumption(Len_ToElongate)")
            Writer.Statement("self.Count_NTPConsumption = NTPConsumption")
            Writer.ReduceSum("self.Count_NTPConsumptionTotal", "NTPConsumption")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_NTPs, -NTPConsumption)")
            Writer.BlankLine()

        with Writer.Statement("def ReleasePPi(self, Count_TotalLengthOfElongation):"):
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_PPi, Count_TotalLengthOfElongation)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateByproducts(self, Len_ToElongate):"):
            Writer.ReduceSum("Count_TotalLengthOfElongationPerRNA", "Len_ToElongate", 1)
            Writer.Statement("self.ConsumeNTPs(Count_TotalLengthOfElongationPerRNA)")
            Writer.ReduceSum("Count_TotalLengthOfElongationTotal", "Len_ToElongate")
            Writer.Statement("self.ReleasePPi(Count_TotalLengthOfElongationTotal)")
            Writer.BlankLine()

        with Writer.Statement("def Elongation(self):"):
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

        with Writer.Statement("def IncrementCountOfRNAs(self, Count_ElongationCompletionPerRNA):"):
            # Increment count of RNAs completed elongation
            Writer.Statement(
                "self.AddToDeltaCounts(self.Cel.Idx_Master_RNAs, Count_ElongationCompletionPerRNA)")
            Writer.BlankLine()

        with Writer.Statement("def Termination(self):"):
            # Reset nascent RNA length matrix
            Writer.Statement("Len_RNAsNascent_ElongationCompletedToZero = self.ResetLengthOfCompletedElongation(self.Cel.Len_RNAsNascent, self.Len_RNAsNascentMax)")
            Writer.Statement("Len_RNAsNascent_ElongationCompletedToNegOne = self.ResetZerosToNegOnes(Len_RNAsNascent_ElongationCompletedToZero)")
            Writer.Statement("self.UpdateNascentRNALengths(Len_RNAsNascent_ElongationCompletedToNegOne)")
            Writer.BlankLine()

            # Update RNA counts
            Writer.ConvToBin("Bin_ElongationCompleted", "Len_RNAsNascent_ElongationCompletedToZero", "==", 0)
            Writer.ReduceSum("Count_ElongationCompletionPerRNA", "Bin_ElongationCompleted", 1)
            Writer.ReduceSum("self.Count_ElongationCompletionTotal", "Count_ElongationCompletionPerRNA")
            Writer.Statement("self.IncrementCountOfRNAs(Count_ElongationCompletionPerRNA)")
            Writer.BlankLine()

        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.PrintStrg("===== Transcription ===== ")
            Writer.PrintStVa("# of RNAPs",
                             "self.Count_RNAP")
            Writer.PrintStVa("# of RNAPs Newly Bound To Gene",
                             "tf.shape(self.Idx_RndRNAsNascent)[0]")
            Writer.PrintStVa("# of Nascent RNAs Elongated (= # of Active RNAPs)",
                             "self.Count_ElongatedTotal")
            Writer.PrintStVa("# of New RNAs Generated",
                             "self.Count_ElongationCompletionTotal")
            Writer.PrintStVa("# of NTP consumption",
                             "self.Count_NTPConsumptionTotal")
            Writer.BlankLine()

