

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):

    # Molecule indices for Molecular IDs
    Idx_Ribosome30S = Comp.Master.ID2Idx_Master['CPLX0-3953']
    Idx_Ribosome50S = Comp.Master.ID2Idx_Master['CPLX0-3956']

    Idx_AAs = ProGen.BuildingBlockIdxs('AA')
    Idx_PPi = Comp.Master.ID2Idx_Master['PPI[c]']

    Idx_SelenoCysteineInAAs = 19
    Idx_AAsLocalAssignmentNoSelenoCysteine = list()
    for i in range(len(Idx_AAs)):
        if i == Idx_SelenoCysteineInAAs:
            continue
        Idx_AAsLocalAssignmentNoSelenoCysteine.append(i)

    # Set up Idx to Idx system in cell.py

    # Temporary parameters
    Rate_ProteinElongation = 20
    Rate_RibosomeActive = 0.22  # 22% of Free Ribosome is active CHECK
    Rate_RibosomeActiveCanBind = 0.5  # 50% of Free Ribosome binds to the promoter CHECK

    NMax_RibosomesPermRNA = 10

    # Temporary references
    NUniq_mRNAs = Comp.RNA.NUniq_mRNAs
    NUniq_Proteins = Comp.Protein.NUniq_Proteins
    assert NUniq_mRNAs == NUniq_Proteins

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.Idx_Proteins", 0)
            Writer.Variable_("self.Idx_RndProteinsNascent", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Bin_ProteinsNascentElongating", 0)
            Writer.Variable_("self.Bin_ProteinsNascentOverElongated", 0)
            Writer.Variable_("self.Bin_ProteinsNascentElongationCompleted", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_ProteinsNascentNew", 0)
            Writer.Variable_("self.Rate_ProteinElongation_Matrix", 0)
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):
            Writer.VarRange_("self.Idx_Proteins", 0, NUniq_Proteins)

            Writer.Variable_("self.Cel.Idx_Ribosome30S", Idx_Ribosome30S)
            Writer.Variable_("self.Cel.Idx_Ribosome50S", Idx_Ribosome50S)

            Writer.Variable_("self.Cel.Idx_AAs", Idx_AAs)
            Writer.Variable_("self.Cel.Idx_PPi", Idx_PPi)
            Writer.Variable_("self.Cel.Idx_SelenoCysteineInAAs", Idx_SelenoCysteineInAAs)
            Writer.Variable_("self.Cel.Idx_AAsLocalAssignmentNoSelenoCysteine", Idx_AAsLocalAssignmentNoSelenoCysteine)
            Writer.BlankLine()

            Writer.Variable_("self.Cel.Rate_RibosomeActive", Rate_RibosomeActive)  # Not implemented yet
            Writer.Variable_("self.Cel.Rate_RibosomeActiveCanBind", Rate_RibosomeActiveCanBind)
            Writer.Variable_("self.Cel.Rate_ProteinElongation", Rate_ProteinElongation)

            # TODO: translational efficiency may be imported from the literature (see flat data)
            Writer.Variable_("self.Cel.Rate_ProteinElongation_Matrix", Rate_ProteinElongation, Shape=[NUniq_mRNAs, NMax_RibosomesPermRNA])
            Writer.BlankLine()

            Writer.InitZeros("self.Cel.Len_ProteinsNascentInitial", [NUniq_Proteins, NMax_RibosomesPermRNA], 'int32')

            Writer.VarFill__("self.Cel.Len_ProteinsNascent", [NUniq_Proteins, NMax_RibosomesPermRNA], -1)
            Writer.Overwrite("self.Cel.Len_ProteinsNascentMax", "self.Cel.Len_Proteins")
            Writer.BlankLine()

            Writer.InitZeros("self.Count_ProteinsNascentNew", NUniq_Proteins, 'int32')
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

        with Writer.Statement("def RetrieveTotalRibosomeCount(self):"):
            Writer.Statement(
                "self.Cel.Count_Ribosome30S = tf.constant(self.Cel.Counts[%s])  # %s: Index for Ribosome" % (Idx_Ribosome30S, Idx_Ribosome30S))
            Writer.Statement(
                "self.Cel.Count_Ribosome50S = tf.constant(self.Cel.Counts[%s])  # %s: Index for Ribosome" % (Idx_Ribosome50S, Idx_Ribosome50S))
            # Ribosome70S count is taken from the lower value between 30S and 50S
            Writer.Minimum__("self.Cel.Count_Ribosome70S", "self.Cel.Count_Ribosome30S", "self.Cel.Count_Ribosome50S")
            # Assume Ribosome count = Ribosome70S count
            Writer.Overwrite("self.Cel.Count_Ribosome", "self.Cel.Count_Ribosome70S")
            Writer.BlankLine()

        with Writer.Statement("def CalculateActiveRibosome(self):"):
            Writer.Cast_____("Count_Ribosome_Float", "self.Cel.Count_Ribosome", 'float32')
            Writer.Multiply_("self.Cel.Count_RibosomeActive", "Count_Ribosome_Float",
                             "self.Cel.Rate_RibosomeActive")  # For summary, active Ribosome
            Writer.RoundInt_("self.Cel.Count_RibosomeActive", "self.Cel.Count_RibosomeActive")
            Writer.BlankLine()

        with Writer.Statement("def CalculateActiveRibosomeCanBind(self):"):
            # Apply the rate of active Ribosome binding (hypothetical)
            Writer.Cast_____("Count_RibosomeActive_Float", "self.Cel.Count_RibosomeActive", 'float32')
            Writer.Multiply_("self.Cel.Count_RibosomeActiveCanBind", "Count_RibosomeActive_Float",
                             "self.Cel.Rate_RibosomeActiveCanBind")  # For summary, active Ribosome portion that can bind to promoter
            Writer.RoundInt_("self.Cel.Count_RibosomeActiveCanBind", "self.Cel.Count_RibosomeActiveCanBind")
            Writer.BlankLine()

        with Writer.Statement("def DetermineRibosomeWillBind(self):"):
            # Ribosome to bind in this step
            Writer.Subtract_("self.Cel.Count_RibosomeWillBind", "self.Cel.Count_RibosomeActiveCanBind",
                             "self.Cel.Count_RibosomeBound")
            Writer.ConvToBin("Bin_RibosomeWillBind", "self.Cel.Count_RibosomeActiveCanBind", ">", "self.Cel.Count_RibosomeBound")
            Writer.Multiply_("self.Cel.Count_RibosomeWillBind", "self.Cel.Count_RibosomeWillBind",
                             "Bin_RibosomeWillBind")  # For Summary, active Ribosome to bind to promoter this step
            Writer.BlankLine()

        with Writer.Statement("def DetermineRibosomeToBind(self):"):
            Writer.Statement("self.RetrieveTotalRibosomeCount()")
            Writer.Statement("self.CalculateActiveRibosome()")
            Writer.Statement("self.CalculateActiveRibosomeCanBind()")
            Writer.Statement("self.DetermineRibosomeWillBind()")
            Writer.BlankLine()

        with Writer.Statement("def SelectProteinsToTranscribe(self):"):
            # # Weighted random distribution later
            Writer.Gather___("Count_mRNAs", "self.Cel.Counts", "self.Cel.Idx_Master_mRNAs")
            Writer.RndIdxWg_("self.Idx_RndProteinsNascent", "self.Cel.Count_RibosomeWillBind", "self.Idx_Proteins", "Count_mRNAs")
            Writer.BlankLine()

        with Writer.Statement("def GetCountOfNewNascentProteins(self):"):
            Writer.InitZeros("self.Count_ProteinsNascentNew", NUniq_Proteins, 'int32')
            Writer.InitOnes_("OnesForRndProteins", "self.Cel.Count_RibosomeWillBind", 'int32')
            Writer.ScatNdUpd("self.Count_ProteinsNascentNew", "self.Idx_RndProteinsNascent", "OnesForRndProteins")
            Writer.BlankLine()

        with Writer.Statement("def UpdateBoundAndUnboundRibosomes(self):"):
            Writer.Add______("self.Cel.Count_RibosomeBound", "self.Cel.Count_RibosomeBound", "self.Cel.Count_RibosomeWillBind")
            Writer.Subtract_("self.Cel.Count_RibosomeUnbound", "self.Cel.Count_Ribosome", "self.Cel.Count_RibosomeBound")
            Writer.BlankLine()

        with Writer.Statement("def DistributeRibosomeTomRNAs(self):"):
            Writer.Statement("self.DetermineRibosomeToBind()")
            Writer.Statement("self.SelectProteinsToTranscribe()")
            Writer.Statement("self.GetCountOfNewNascentProteins()")
            Writer.Statement("self.UpdateBoundAndUnboundRibosomes()")
            Writer.BlankLine()

        with Writer.Statement("def Initiation(self):"):
            Writer.Statement("self.SigmaFactorBinding()  # Not implemented")
            Writer.Statement("self.DistributeRibosomeTomRNAs()")
            Writer.BlankLine()

        with Writer.Statement("def GetNascentProteinCounts(self):"):
            Writer.Gather___("self.Cel.Count_ProteinsNascent_Matrix", "self.Cel.Counts", "self.Cel.Idx_Master_ProteinsNascent")
            Writer.ReduceSum("self.Cel.Count_ProteinsNascentInitialTotal", "self.Cel.Count_ProteinsNascent_Matrix")
            # Temporary code
            # Writer.InitOnes_("self.Cel.Count_ProteinsNascent_Matrix", NUniq_Proteins, 'int32')
            Writer.BlankLine()

        with Writer.Statement("def DetermineNascentProteinIndices(self):"):
            # Generate boolean and binary matrices for element availability in nascent proteins length matrix.
            Writer.Less_____("Bool_ProteinsNascentAvailable", "self.Cel.Len_ProteinsNascent", 0)
            Writer.BoolToBin("Bin_ProteinsNascentAvailable", "Bool_ProteinsNascentAvailable")

            # Generate a cumulative sum matrix for labeling the available elements with count.
            Writer.CumSum___("Cumsum_ProteinsNascent", "Bin_ProteinsNascentAvailable", 1)

            # Generate a boolean matrix of the count label matrix > 0, < count.
            Writer.Greater__("Bool_CumsumGreaterThanZero", "Cumsum_ProteinsNascent", 0)
            Writer.Reshape__("self.Count_ProteinsNascentNew", "self.Count_ProteinsNascentNew", [-1, 1])
            Writer.LessEq___("Bool_CumsumLessThanOrEqualToCount", "Cumsum_ProteinsNascent", "self.Count_ProteinsNascentNew")
            Writer.LogicAnd_("Bool_Cumsum_ProteinsNascent", "Bool_CumsumGreaterThanZero", "Bool_CumsumLessThanOrEqualToCount")

            # Generate a binary matrix for elements satisfying both availability and count conditions.
            Writer.LogicAnd_("Bool_ProteinsNascentNew", "Bool_ProteinsNascentAvailable", "Bool_Cumsum_ProteinsNascent")
            Writer.BoolToBin("Bin_ProteinsNascentNew", "Bool_ProteinsNascentNew")

            # Turn -1 to 0 where nascent proteins are to be made in the nascent protein length table.
            Writer.Add______("self.Cel.Len_ProteinsNascent", "self.Cel.Len_ProteinsNascent", "Bin_ProteinsNascentNew")
            Writer.BlankLine()

        with Writer.Statement("def ElongateNascentProteinsInLengthMatrix(self):"):
            # Get a binary matrix for nascent proteins to elongate
            Writer.ConvToBin("self.Bin_ProteinsNascentElongating", "self.Cel.Len_ProteinsNascent", ">=", 0)

            # Rate Matrix * bin(boolean) + length table (overelongated not counted yet)
            Writer.Multiply_("Rate_ProteinElongation_Matrix", "self.Bin_ProteinsNascentElongating",
                             "self.Cel.Rate_ProteinElongation")
            Writer.Add______("self.Cel.Len_ProteinsNascentElongated", "self.Cel.Len_ProteinsNascent",
                             "Rate_ProteinElongation_Matrix")
            Writer.BlankLine()

        with Writer.Statement("def DetermineOverElongatedNascentProteins(self):"):
            # compare with len MAX to get boolean matrix for overelongated (>= MAX)
            Writer.Reshape__("self.Cel.Len_ProteinsNascentMax", "self.Cel.Len_ProteinsNascentMax", [-1, 1])
            Writer.ConvToBin("self.Bin_ProteinsNascentOverElongated", "self.Cel.Len_ProteinsNascentElongated", ">",
                             "self.Cel.Len_ProteinsNascentMax")
            Writer.ReduceSum("self.Cel.Count_ProteinsNascentOverElongatedTotal",
                             "self.Bin_ProteinsNascentOverElongated")  # For summary, the number of completed/overelongated Proteins
            Writer.BlankLine()

        with Writer.Statement("def DetermineAmountOfOverElongatedAAs(self):"):
            Writer.Multiply_("MaxLengthForOverElongated", "self.Bin_ProteinsNascentOverElongated",
                             "self.Cel.Len_ProteinsNascentMax")
            Writer.Multiply_("LengthForOverElongatedOnly", "self.Bin_ProteinsNascentOverElongated",
                             "self.Cel.Len_ProteinsNascentElongated")
            # Get the over-elongated length of nascent Proteins
            Writer.Subtract_("self.Cel.Len_ProteinsNascentOverElongated", "LengthForOverElongatedOnly",
                             "MaxLengthForOverElongated")
            Writer.ReduceSum("self.Cel.Count_AAsOverElongated",
                             "self.Cel.Len_ProteinsNascentOverElongated")  # For summary, over consumped AAs to be corrected
            Writer.BlankLine()

        with Writer.Statement("def CorrectOverElongationInLengthMatrix(self):"):
            # Update nascent Protein length by removing over elongated length
            Writer.Subtract_("self.Cel.Len_ProteinsNascentAdjusted", "self.Cel.Len_ProteinsNascentElongated",
                             "self.Cel.Len_ProteinsNascentOverElongated")
            # now all completed elongation has its max value

            Writer.BlankLine()

        with Writer.Statement("def CorrectElongatedLengthForProteins(self):"):
            # Initial elongation rate
            Writer.Multiply_("Rate_ProteinElongation_Matrix_Initial", "self.Cel.Rate_ProteinElongation_Matrix",
                             "self.Bin_ProteinsNascentElongating")
            # Adjusted elongation rate (initial - elongated)
            Writer.Subtract_("self.Rate_ProteinElongation_Matrix", "Rate_ProteinElongation_Matrix_Initial",
                             "self.Cel.Len_ProteinsNascentOverElongated")
            # Calculate the length of elongation per mRNA
            Writer.ReduceSum("self.Cel.Count_ProteinElongationLengthPerProtein", "self.Rate_ProteinElongation_Matrix",
                             1)  # For summary, total elongation length per mRNA
            # Calculate the total length of elongation
            Writer.ReduceSum("self.Cel.Count_ProteinElongationLengthTotal",
                             "self.Rate_ProteinElongation_Matrix")  # For summary, total elongation length in this sim step
            Writer.BlankLine()

        with Writer.Statement("def ElongatePolypeptides(self):"):
            # Nascent Proteins will be elongated.
            Writer.Statement("self.GetNascentProteinCounts()")
            Writer.Statement("self.DetermineNascentProteinIndices()")
            Writer.Statement("self.ElongateNascentProteinsInLengthMatrix()")
            Writer.Statement("self.DetermineOverElongatedNascentProteins()")
            Writer.Statement("self.DetermineAmountOfOverElongatedAAs()")
            Writer.Statement("self.CorrectOverElongationInLengthMatrix()")
            Writer.Statement("self.CorrectElongatedLengthForProteins()")
            Writer.BlankLine()

        with Writer.Statement("def GetRawAAConsumption(self):"):
            # Prepare NT elongation length for each mRNA matrix for multiplication
            Writer.Reshape__("self.Cel.Count_ProteinElongationLengthPerProtein", "self.Cel.Count_ProteinElongationLengthPerProtein",
                             [-1, 1])
            Writer.Cast_____("self.Cel.Count_ProteinElongationLengthPerProtein", "self.Cel.Count_ProteinElongationLengthPerProtein",
                             'float32')
            # AA Frequency per mRNA * elongating Polypeptide (in NT count) per mRNA
            Writer.MatrixMul("AAConsumption_Raw", "self.Cel.Freq_AAsInProteins",
                             "self.Cel.Count_ProteinElongationLengthPerProtein")
            Writer.ReturnVar("AAConsumption_Raw")
            Writer.BlankLine()

        with Writer.Statement("def GetRoundedAAConsumption(self, AAConsumption_Raw):"):
            Writer.RoundInt_("AAConsumption_Rounded", "AAConsumption_Raw")
            Writer.ReturnVar("AAConsumption_Rounded")
            Writer.BlankLine()

        with Writer.Statement("def CalculateDiscrepancy(self, AAConsumption_Rounded):"):
            # Determine the difference between rounded vs corrected elongation
            Writer.ReduceSum("Sum_Rounded", "AAConsumption_Rounded")
            Writer.Cast_____("Sum_Rounded", "Sum_Rounded", 'int32')
            Writer.ReduceSum("Sum_CorrectedElongationRate", "self.Rate_ProteinElongation_Matrix")
            Writer.Subtract_("DeltaSum", "Sum_CorrectedElongationRate", "Sum_Rounded")
            Writer.ReturnVar("DeltaSum")
            Writer.BlankLine()

        with Writer.Statement("def AdjustAAConsumption(self, AAConsumption_Rounded, Discrepancy):"):
            # Get equal amount of AAs if discrepancy is greater than or equal to 4 or less than or equal to -4.
            Writer.FloorDiv_("N_AASets", "Discrepancy", "self.Cel.NUniq_AAs")
            Writer.VarRepeat("N_AASets", "N_AASets", "self.Cel.NUniq_AAs")
            Writer.Add______("AAConsumption_MissingSet", "AAConsumption_Rounded", "N_AASets")
            Writer.BlankLine()

            # Get random AA for the remainder (replace with weighted random AA based on aa frequency in all mRNAs)
            Writer.Remainder("N_AARemainder", "Discrepancy", "self.Cel.NUniq_AAs")
            Writer.Reshape__("N_AARemainder", "N_AARemainder", -1)
            Writer.InitZeros("AAConsumption_MissingRemainder", "self.Cel.NUniq_AAs", 'int32')
            Writer.RndNumUni("Idx_Remainder", "N_AARemainder", "0", "self.Cel.NUniq_AAs")
            Writer.InitOnes_("OnesForRemainder", "N_AARemainder", 'int32')
            Writer.ScatNdAdd("AAConsumption_MissingRemainder", "Idx_Remainder", "OnesForRemainder")
            Writer.BlankLine()

            # Calculate adjusted AA Consumption
            Writer.Add______("AAConsumption_Adjusted", "AAConsumption_MissingSet", "AAConsumption_MissingRemainder")

            # Correction for Seleno-cysteine when -1
            Writer.Gather___("Count_SelenoCysteine", "AAConsumption_Adjusted", "self.Cel.Idx_SelenoCysteineInAAs")
            Writer.ConvToBin("Bin_SelenoCysteine", "Count_SelenoCysteine", "<", 0)
            Writer.Multiply_("CountToAdjust", "Bin_SelenoCysteine", "self.Cel.One")
            Writer.ScatNdAdd("AAConsumption_Adjusted", "self.Cel.Idx_SelenoCysteineInAAs", "CountToAdjust")


            Writer.RndIdxUni("Idx_RndAA", [1],
                             "self.Cel.Idx_AAsLocalAssignmentNoSelenoCysteine")
            Writer.NegValue_("CountToAdjust_Neg", "CountToAdjust")
            Writer.ScatNdAdd("AAConsumption_Adjusted", "Idx_RndAA", "CountToAdjust_Neg")

            # Return the adjusted AA Consumption
            Writer.ReduceSum("TotalAAConsumption", "AAConsumption_Adjusted")
            Writer.AsrtElEq_("TotalAAConsumption", "self.Cel.Count_ProteinElongationLengthTotal")
            Writer.ReturnVar("AAConsumption_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def DetermineAAConsumption(self):"):
            Writer.Statement("AAConsumption_Raw = self.GetRawAAConsumption()")
            Writer.Statement("AAConsumption_Rounded = self.GetRoundedAAConsumption(AAConsumption_Raw)")
            Writer.Statement("Discrepancy = self.CalculateDiscrepancy(AAConsumption_Rounded)")
            Writer.Statement("AAConsumption_Adjusted = self.AdjustAAConsumption(AAConsumption_Rounded, Discrepancy)")
            Writer.ReturnVar("AAConsumption_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def ConsumeAAs(self):"):
            Writer.Statement("self.Cel.Count_ProteinElongationAAConsumption = self.DetermineAAConsumption()")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_AAs, -self.Cel.Count_ProteinElongationAAConsumption)")
            Writer.BlankLine()

        with Writer.Statement("def ReleasePPi(self):"):
            Writer.Overwrite("self.Cel.Count_ProteinElongationPPiProduction", "self.Cel.Count_ProteinElongationLengthTotal")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_PPi, self.Cel.Count_ProteinElongationPPiProduction)")
            Writer.BlankLine()

        with Writer.Statement("def Elongation(self):"):
            Writer.Statement("self.ElongatePolypeptides()")
            Writer.Statement("self.ConsumeAAs()")
            Writer.Statement("self.ReleasePPi()")
            Writer.BlankLine()

        with Writer.Statement("def IdentifyCompletedProteinElongation(self):"):
            # identify all NMAX values by comparison (boolean)
            Writer.ConvToBin("self.Bin_ProteinsNascentElongationCompleted","self.Cel.Len_ProteinsNascentAdjusted", "==",
                             "self.Cel.Len_ProteinsNascentMax")
            Writer.BlankLine()

        with Writer.Statement("def CountCompletedProteinElongation(self):"):
            # Count all completed per Protein and total
            Writer.ReduceSum("self.Cel.Count_ProteinElongationCompletedPerProtein",
                             "self.Bin_ProteinsNascentElongationCompleted", 1)
            Writer.ReduceSum("self.Cel.Count_ProteinElongationCompletedTotal",
                             "self.Bin_ProteinsNascentElongationCompleted")  # For summary, total completed elongation
            Writer.BlankLine()

        with Writer.Statement("def ResetLengthOfCompletedNascentProteins(self):"):
            # Reset all NMAX lengths to zero by subtracting max length from completed
            Writer.Multiply_("self.Cel.Len_ProteinsNascentCompleted", "self.Cel.Len_ProteinsNascentMax",
                             "self.Bin_ProteinsNascentElongationCompleted")
            Writer.Subtract_("self.Cel.Len_ProteinsNascentFinal", "self.Cel.Len_ProteinsNascentAdjusted",
                             "self.Cel.Len_ProteinsNascentCompleted")

            # Check reset status
            Writer.Multiply_("CheckCompletionReset", "self.Cel.Len_ProteinsNascentFinal",
                             "self.Bin_ProteinsNascentElongationCompleted")
            Writer.AsrtElEq_("CheckCompletionReset", 0)
            Writer.BlankLine()

        with Writer.Statement("def GetIndicesOfProteinsCompletedElongation(self):"):
            # Get indices of Proteins completed elongation. Currently not used
            Writer.GetIdx___("Idx_ProteinElongationCompleted", "self.Cel.Count_ProteinElongationCompletedPerProtein", ">", 0)
            Writer.Gather___("self.Cel.Idx_ProteinElongationCompleted", "self.Cel.Idx_ProteinsNascent",
                             "Idx_ProteinElongationCompleted")
            Writer.BlankLine()

        with Writer.Statement("def ReleaseRibosome(self):"):
            Writer.Overwrite("self.Cel.Count_RibosomeReleased", "self.Cel.Count_ProteinElongationCompletedTotal")
            Writer.Subtract_("self.Cel.Count_RibosomeBound", "self.Cel.Count_RibosomeBound", "self.Cel.Count_RibosomeReleased")
            Writer.Add______("self.Cel.Count_RibosomeUnbound", "self.Cel.Count_RibosomeUnbound", "self.Cel.Count_RibosomeReleased")
            Writer.Add______("SumOfBoundUnbound", "self.Cel.Count_RibosomeBound", "self.Cel.Count_RibosomeUnbound")
            Writer.AsrtElEq_("self.Cel.Count_Ribosome", "SumOfBoundUnbound")
            Writer.BlankLine()

        with Writer.Statement("def IncrementCountOfProteins(self):"):
            # Increment count of Proteins completed elongation
            Writer.Statement(
                "self.AddToDeltaCounts(self.Cel.Idx_Master_Proteins, self.Cel.Count_ProteinElongationCompletedPerProtein)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateLengthOfNascentProteinsMatrix(self):"):
            # Overwrite self.Cel.Len_ProteinsNascentInitial with self.Cel.Len_ProteinsNascentFinal
            Writer.Overwrite("self.Cel.Len_ProteinsNascent", "self.Cel.Len_ProteinsNascentFinal")
            Writer.BlankLine()

        with Writer.Statement("def Termination(self):"):
            Writer.Statement("self.IdentifyCompletedProteinElongation()")
            Writer.Statement("self.CountCompletedProteinElongation()")
            Writer.Statement("self.ResetLengthOfCompletedNascentProteins()")
            # Writer.Statement("self.GetIndicesOfProteinsCompletedElongation()")
            Writer.Statement("self.ReleaseRibosome()")
            Writer.Statement("self.IncrementCountOfProteins()")
            Writer.Statement("self.UpdateLengthOfNascentProteinsMatrix()")
            Writer.BlankLine()

        with Writer.Statement("def Debugging_Protein_Indices(self, Str):"):
            Writer.Statement("print(Str)")
            Writer.BlankLine()

            Writer.Comment__("Protein length Indices")
            Writer.PrintStrg("Protein length Indices")

            Writer.PrintVar_("self.Idx_RndProteinsNascent")

            Writer.PrtIdxGr_("self.Cel.Len_ProteinsNascentInitial", 0)
            Writer.PrtIdxGr_("self.Cel.Len_ProteinsNascentElongated", 0)
            Writer.PrtIdxGr_("self.Cel.Len_ProteinsNascentOverElongated", 0)

            Writer.PrtIdxGr_("self.Rate_ProteinElongation_Matrix", 0)

            Writer.PrtIdxGr_("self.Bin_ProteinsNascentElongating", 0)
            Writer.PrtIdxGr_("self.Bin_ProteinsNascentOverElongated", 0)
            Writer.PrtIdxGr_("self.Bin_ProteinsNascentElongationCompleted", 0)

            Writer.BlankLine()

        with Writer.Statement("def Debugging_Count_ProteinsNascent(self, Str):"):
            Writer.Statement("print(Str)")
            Writer.BlankLine()

            Writer.Comment__("Nascent Proteins")
            Writer.PrintStrg("Nascent Proteins")

            Writer.ReduceSum("Reduced_Sum_Of_self_Cel_Count_ProteinsNascent", "self.Cel.Count_ProteinsNascent")
            Writer.PrintVaVa("Reduced_Sum_Of_self_Cel_Count_ProteinsNascent")
            Writer.BlankLine()

            Writer.Gather___("Counts_ProteinsNascent", "self.Cel.Counts", "self.Cel.Idx_Master_ProteinsNascent")
            Writer.ReduceSum("Reduced_Sum_Of_Counts_ProteinsNascent", "Counts_ProteinsNascent")
            Writer.PrintVaVa("Reduced_Sum_Of_Counts_ProteinsNascent")
            Writer.BlankLine()

            Writer.Gather___("DeltaCounts_ProteinsNascent", "self.Cel.DeltaCounts", "self.Cel.Idx_Master_ProteinsNascent")
            Writer.ReduceSum("Reduced_Sum_Of_DeltaCounts_ProteinsNascent", "DeltaCounts_ProteinsNascent")
            Writer.PrintVaVa("Reduced_Sum_Of_DeltaCounts_ProteinsNascent")
            Writer.BlankLine()

            Writer.Comment__("Proteins (full)")
            Writer.PrintStrg("Proteins (full)")

            Writer.ReduceSum("Reduced_Sum_Of_self_Cel_Count_Proteins", "self.Cel.Count_Proteins")
            Writer.PrintVaVa("Reduced_Sum_Of_self_Cel_Count_Proteins")
            Writer.BlankLine()

            Writer.Gather___("Counts_Proteins", "self.Cel.Counts", "self.Cel.Idx_Master_Proteins")
            Writer.ReduceSum("Reduced_Sum_Of_Counts_Proteins", "Counts_Proteins")
            Writer.PrintVaVa("Reduced_Sum_Of_Counts_Proteins")
            Writer.BlankLine()

            Writer.Gather___("DeltaCounts_Proteins", "self.Cel.DeltaCounts", "self.Cel.Idx_Master_Proteins")
            Writer.ReduceSum("Reduced_Sum_Of_DeltaCounts_Proteins", "DeltaCounts_Proteins")
            Writer.PrintVaVa("Reduced_Sum_Of_DeltaCounts_Proteins")
            Writer.BlankLine()

        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.PrintStrg("===== Translation Initiation ===== ")
            # Number of Total Ribosome
            Writer.PrintStVa("# of Total Ribosomes",
                             "self.Cel.Count_Ribosome")
            # Number of Active Ribosome
            Writer.PrintStVa("# of Active Ribosomes",
                             "self.Cel.Count_RibosomeActive")
            # Number of Ribosome binding
            Writer.PrintStVa("# of Ribosomes available that can bind a promoter",
                             "self.Cel.Count_RibosomeActiveCanBind")
            # Number of new Ribosome binding in this step (= new nascent Proteins to elongate next step)
            Writer.PrintStVa("# of Ribosomes that newly bind a promoter in this step",
                             "self.Cel.Count_RibosomeWillBind")
            Writer.BlankLine()

            Writer.PrintStrg("===== Polypeptide Elongation ===== ")
            # Number of Nascent Proteins elongated
            Writer.PrintStVa("# of All Nascent Proteins Elongating",
                             "self.Cel.Count_ProteinsNascentElongatingTotal")
            # Total elongation length of Proteins
            Writer.PrintStVa("Total Elongation Length of Proteins (aa)",
                             "self.Cel.Count_ProteinElongationLengthTotal")
            # Total AA consumption and PPi production
            Writer.PrintStVa("Total AA Consumption [A,R,N,D,C,E,Q,G,H,I,L,K,M,F,P,S,T,W,Y,O,V]",
                             "self.Cel.Count_ProteinElongationAAConsumption")
            Writer.PrintStVa("Total PPi Production",
                             "self.Cel.Count_ProteinElongationPPiProduction")
            Writer.BlankLine()

            Writer.PrintStrg("===== Translation Termination ===== ")
            # Number of Protein elongation Completed
            Writer.PrintStVa("# of Protein Elongation Completed",
                             "self.Cel.Count_ProteinElongationCompletedTotal")
            # Number of Ribosomes released
            Writer.PrintStVa("# of Ribosomes Released",
                             "self.Cel.Count_RibosomeReleased")
            # Number of Ribosome bound to DNA
            Writer.PrintStVa("# of Total Ribosomes Bound to mRNA",
                             "self.Cel.Count_RibosomeBound")
            # Number of Ribosome freely floating
            Writer.PrintStVa("# of Total Ribosomes Unbound Floating",
                             "self.Cel.Count_RibosomeUnbound")
            Writer.BlankLine()
