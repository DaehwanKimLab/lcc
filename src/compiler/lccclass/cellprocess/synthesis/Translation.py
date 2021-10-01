

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):

    # Molecule indices for Molecular IDs
    Idx_Ribosome30S = Comp.Master.ID2Idx_Master['CPLX0-3953']
    Idx_Ribosome50S = Comp.Master.ID2Idx_Master['CPLX0-3956']

    Idx_AAs = ProGen.BuildingBlockIdxs('AAs')

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
            Writer.Variable_("self.Rate_RibosomeActive", 0)
            Writer.Variable_("self.Rate_RibosomeActiveCanBind", 0)
            Writer.Variable_("self.Rate_ProteinElongation", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_Ribosome", 0)
            Writer.Variable_("self.Count_ElongatedTotal", 0)
            Writer.Variable_("self.Count_ElongationCompletionTotal", 0)
            Writer.Variable_("self.Count_AAConsumptionTotal", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Len_ProteinsNascentMax", 0)
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):
            Writer.VarRange_("self.Idx_Proteins", 0, NUniq_Proteins)

            Writer.Variable_("self.Cel.Idx_Ribosome30S", Idx_Ribosome30S)
            Writer.Variable_("self.Cel.Idx_Ribosome50S", Idx_Ribosome50S)

            Writer.Variable_("self.Cel.Idx_SelenoCysteineInAAs", Idx_SelenoCysteineInAAs)
            Writer.Variable_("self.Cel.Idx_AAsLocalAssignmentNoSelenoCysteine", Idx_AAsLocalAssignmentNoSelenoCysteine)
            Writer.BlankLine()

            Writer.Variable_("self.Rate_RibosomeActive", Rate_RibosomeActive)
            Writer.Variable_("self.Rate_RibosomeActiveCanBind", Rate_RibosomeActiveCanBind)
            Writer.Variable_("self.Rate_ProteinElongation", Rate_ProteinElongation)
            # TODO: translational efficiency may be imported from the literature (see flat data)

            Writer.VarFill__("self.Cel.Len_ProteinsNascent", [NUniq_Proteins, NMax_RibosomesPermRNA], -1)
            Writer.Reshape__("self.Len_ProteinsNascentMax", "self.Cel.Len_Proteins", [-1, 1])
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def ExecuteProcess(self):"):
            Writer.Statement("self.Initiation()")
            Writer.Statement("self.Elongation()")
            Writer.Statement("self.Termination()")
            Writer.BlankLine()

        with Writer.Statement("def DetermineRibosomeToBind(self):"):

            Writer.Comment__("RetrieveTotalRibosomeCount")
            Writer.Statement(
                "Count_Ribosome30S = tf.constant(self.Cel.Counts[%s])  # %s: Index for Ribosome" % (
                Idx_Ribosome30S, Idx_Ribosome30S))
            Writer.Statement(
                "Count_Ribosome50S = tf.constant(self.Cel.Counts[%s])  # %s: Index for Ribosome" % (
                Idx_Ribosome50S, Idx_Ribosome50S))

            Writer.Comment__("Ribosome70S count is taken from the lower value between 30S and 50S")
            Writer.Minimum__("Count_Ribosome70S", "Count_Ribosome30S",
                             "Count_Ribosome50S")

            Writer.Comment__("Assume Ribosome count = Ribosome70S count")
            Writer.Overwrite("Count_Ribosome", "Count_Ribosome70S")
            Writer.Overwrite("self.Count_Ribosome", "Count_Ribosome70S")
            Writer.BlankLine()

            Writer.Comment__("CalculateActiveRibosome")
            Writer.Cast_____("Count_Ribosome_Float", "Count_Ribosome", 'float32')
            Writer.Multiply_("Count_RibosomeActive", "Count_Ribosome_Float",
                             "self.Rate_RibosomeActive")  # For summary, active Ribosome
            Writer.RoundInt_("Count_RibosomeActive", "Count_RibosomeActive")
            Writer.BlankLine()

            Writer.Comment__("CalculateActiveRibosomeCanBind")
            Writer.Cast_____("Count_RibosomeActive_Float", "Count_RibosomeActive", 'float32')
            Writer.Multiply_("Count_RibosomeActiveCanBind", "Count_RibosomeActive_Float",
                             "self.Rate_RibosomeActiveCanBind")  # For summary, active Ribosome portion that can bind to promoter
            Writer.RoundInt_("Count_RibosomeActiveCanBind", "Count_RibosomeActiveCanBind")

            Writer.Comment__("DetermineRibosomeWillBind")
            Writer.ConvToBin("Bin_RibosomeBound", "self.Cel.Len_ProteinsNascent", ">=", 0)
            Writer.ReduceSum("Count_RibosomeBound", "Bin_RibosomeBound")
            Writer.Subtract_("Count_RibosomeWillBind", "Count_RibosomeActiveCanBind",
                             "Count_RibosomeBound")
            Writer.ConvToBin("Bin_RibosomeWillBind", "Count_RibosomeActiveCanBind", ">",
                             "Count_RibosomeBound")
            Writer.Multiply_("Count_RibosomeWillBind", "Count_RibosomeWillBind",
                             "Bin_RibosomeWillBind")  # For Summary, active Ribosome to bind to promoter this step
            Writer.ReturnVar("Count_RibosomeWillBind")
            Writer.BlankLine()

        with Writer.Statement("def SelectProteinsToTranslate(self, Count_RibosomeWillBind):"):
            # Weighted random distribution later
            Writer.Gather___("Count_mRNAs", "self.Cel.Counts", "self.Cel.Idx_Master_mRNAs")
            Writer.Statement("Idx_RndProteinsNascent = self.PickRandomIndexFromPool_Weighted_Local(Count_RibosomeWillBind, self.Idx_Proteins, Count_mRNAs)")
            Writer.ReturnVar("Idx_RndProteinsNascent")
            Writer.BlankLine()

        with Writer.Statement("def AddNewCountsToNascentProteinLengthMatrix(self, Len_Matrix, Idx_NewSelected):"):
            # To turn -1 to 0 where nascent proteins are to be made in the nascent protein length table.
            Writer.InitZeros("Count_ProteinsNascentNew", NUniq_Proteins, 'int32')
            Writer.OnesLike_("Ones", "Idx_NewSelected", 'int32')
            Writer.ScatNdAdd("Count_ProteinsNascentNew", "Count_ProteinsNascentNew", "Idx_NewSelected", "Ones")

            Writer.Statement("Bin_AddNewCount = self.GetBinToAddNewCountToLenMatrix(Len_Matrix, Count_ProteinsNascentNew)")
            Writer.Add______("Len_Matrix_NewCountAdded", "Len_Matrix", "Bin_AddNewCount")
            Writer.ReturnVar("Len_Matrix_NewCountAdded")
            Writer.BlankLine()

        with Writer.Statement("def DistributeRibosomesTomRNAs(self):"):
            Writer.Statement("Count_RibosomeWillBind = self.DetermineRibosomeToBind()")
            Writer.Statement("Idx_RndProteinsNascent = self.SelectProteinsToTranslate(Count_RibosomeWillBind)")
            Writer.ReturnVar("Idx_RndProteinsNascent")

            # TODO: RNAP-bound mRNAs are committed to translation, resistant to RNA degradation. Make a separate matrix?
            Writer.BlankLine()

        with Writer.Statement("def UpdateNascentProteinLengths(self, Len_ProteinsNascent_Update):"):
            Writer.Overwrite("self.Cel.Len_ProteinsNascent", "Len_ProteinsNascent_Update")
            Writer.BlankLine()

        with Writer.Statement("def Initiation(self):"):
            Writer.Statement("self.Idx_RndProteinsNascent = self.DistributeRibosomesTomRNAs()")

            # Update the nascent protein length matrix
            Writer.Statement("Len_ProteinNascent_NewToZero = self.AddNewCountsToNascentProteinLengthMatrix(self.Cel.Len_ProteinsNascent, self.Idx_RndProteinsNascent)")
            Writer.Statement("self.UpdateNascentProteinLengths(Len_ProteinNascent_NewToZero)")
            Writer.BlankLine()

        with Writer.Statement("def DetermineAAConsumption(self, Len_ToElongate):"):
            Writer.Statement("AAConsumption = self.DetermineAmountOfBuildingBlocks(Len_ToElongate, self.Cel.Freq_AAsInProteins)")

            # Correction for Seleno-cysteine when -1
            Writer.Gather___("Count_SelenoCysteine", "AAConsumption", "self.Cel.Idx_SelenoCysteineInAAs")
            Writer.ConvToBin("Bin_SelenoCysteine", "Count_SelenoCysteine", "<", 0)
            Writer.Multiply_("CountToAdjust", "Bin_SelenoCysteine", "self.Cel.One")
            Writer.ScatNdAdd("AAConsumption_Adjusted", "AAConsumption", "self.Cel.Idx_SelenoCysteineInAAs", "CountToAdjust")

            Writer.RndIdxUni("Idx_RndAA", [1],
                             "self.Cel.Idx_AAsLocalAssignmentNoSelenoCysteine")
            Writer.Negative_("CountToAdjust_Neg", "CountToAdjust")
            Writer.ScatNdAdd("AAConsumption_Adjusted", "AAConsumption_Adjusted", "Idx_RndAA", "CountToAdjust_Neg")
            Writer.BlankLine()

            # Return the adjusted AA Consumption
            Writer.ReduceSum("TotalAAConsumptionFinal", "AAConsumption_Adjusted")
            Writer.ReduceSum("TotalAAConsumptionInitial", "Len_ToElongate")
            Writer.AsrtEq___("TotalAAConsumptionFinal", "TotalAAConsumptionInitial")
            Writer.ReturnVar("AAConsumption_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def ConsumeAAs(self, Len_ToElongate):"):
            Writer.Statement("AAConsumption = self.DetermineAAConsumption(Len_ToElongate)")
            Writer.ReduceSum("self.Count_AAConsumptionTotal", "AAConsumption")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_AAs, -AAConsumption)")
            Writer.BlankLine()

        with Writer.Statement("def ReleasePPi(self, Count_TotalLengthOfElongation):"):
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_PPi, Count_TotalLengthOfElongation)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateByproducts(self, Len_ToElongate):"):
            Writer.ReduceSum("Count_TotalLengthOfElongationPerProtein", "Len_ToElongate", 1)
            Writer.Statement("self.ConsumeAAs(Count_TotalLengthOfElongationPerProtein)")
            Writer.ReduceSum("Count_TotalLengthOfElongationTotal", "Len_ToElongate")
            Writer.Statement("self.ReleasePPi(Count_TotalLengthOfElongationTotal)")
            Writer.BlankLine()

        with Writer.Statement("def Elongation(self):"):
            Writer.Statement("Len_ProteinsNascent_Elongated = self.DetermineAmountOfElongation(self.Cel.Len_ProteinsNascent, self.Rate_ProteinElongation, self.Len_ProteinsNascentMax)")
            Writer.Statement("Len_ToElongate = Len_ProteinsNascent_Elongated - self.Cel.Len_ProteinsNascent")
            Writer.ConvToBin("Bin_Elongated", "Len_ToElongate", ">", 0)
            Writer.ReduceSum("self.Count_ElongatedTotal", "Bin_Elongated")
            Writer.Statement("self.UpdateNascentProteinLengths(Len_ProteinsNascent_Elongated)")
            Writer.Statement("self.UpdateByproducts(Len_ToElongate)")
            Writer.BlankLine()

        # with Writer.Statement("def GetIndicesOfProteinsCompletedElongation(self):"):
        #     # Get indices of Proteins completed elongation. Currently not used
        #     Writer.GetIdx___("Idx_ProteinElongationCompleted", "self.Cel.Count_ProteinElongationCompletedPerProtein", ">", 0)
        #     Writer.Gather___("self.Cel.Idx_ProteinElongationCompleted", "self.Cel.Idx_ProteinsNascent",
        #                      "Idx_ProteinElongationCompleted")
        #     Writer.BlankLine()

        with Writer.Statement("def IncrementCountOfProteins(self, Count_ElongationCompletionPerProtein):"):
            # Increment count of Proteins completed elongation
            Writer.Statement(
                "self.AddToDeltaCounts(self.Cel.Idx_Master_Proteins, Count_ElongationCompletionPerProtein)")
            Writer.BlankLine()

        with Writer.Statement("def Termination(self):"):
            # Reset nascent protein length matrix
            Writer.Statement("Len_ProteinsNascent_ElongationCompletedToZero = self.ResetLengthOfCompletedElongation(self.Cel.Len_ProteinsNascent, self.Len_ProteinsNascentMax)")
            Writer.Statement("Len_ProteinsNascent_ElongationCompletedToNegOne = self.ResetZerosToNegOnes(Len_ProteinsNascent_ElongationCompletedToZero)")
            Writer.Statement("self.UpdateNascentProteinLengths(Len_ProteinsNascent_ElongationCompletedToNegOne)")
            Writer.BlankLine()

            # Update protein counts
            Writer.ConvToBin("Bin_ElongationCompleted", "Len_ProteinsNascent_ElongationCompletedToZero", "==", 0)
            Writer.ReduceSum("Count_ElongationCompletionPerProtein", "Bin_ElongationCompleted", 1)
            Writer.ReduceSum("self.Count_ElongationCompletionTotal", "Count_ElongationCompletionPerProtein")
            Writer.Statement("self.IncrementCountOfProteins(Count_ElongationCompletionPerProtein)")
            Writer.BlankLine()

        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.PrintStrg("===== Translation ===== ")
            Writer.PrintStVa("# of Ribosomes",
                             "self.Count_Ribosome")
            Writer.PrintStVa("# of Ribosomes Newly Bound To mRNA",
                             "tf.shape(self.Idx_RndProteinsNascent)[0]")
            Writer.PrintStVa("# of Nascent Proteins Elongated (= # of Active Ribosomes)",
                             "self.Count_ElongatedTotal")
            Writer.PrintStVa("# of New Proteins Generated",
                             "self.Count_ElongationCompletionTotal")
            Writer.PrintStVa("# of AA consumption",
                             "self.Count_AAConsumptionTotal")
            Writer.BlankLine()

            # Writer.PrintStrg("===== Polypeptide Elongation ===== ")
            # # Number of Nascent Proteins elongated
            # Writer.PrintStVa("# of All Nascent Proteins Elongating",
            #                  "self.Cel.Count_ProteinsNascentElongatingTotal")
            # # Total elongation length of Proteins
            # Writer.PrintStVa("Total Elongation Length of Proteins (aa)",
            #                  "self.Cel.Count_ProteinElongationLengthTotal")
            # # Total AA consumption and PPi production
            # Writer.PrintStVa("Total AA Consumption [A,R,N,D,C,E,Q,G,H,I,L,K,M,F,P,S,T,W,Y,O,V]",
            #                  "self.Cel.Count_ProteinElongationAAConsumption")
            # Writer.PrintStVa("Total PPi Production",
            #                  "self.Cel.Count_ProteinElongationPPiProduction")
            # Writer.BlankLine()
            #
            # Writer.PrintStrg("===== Translation Termination ===== ")
            # # Number of Protein elongation Completed
            # Writer.PrintStVa("# of Protein Elongation Completed",
            #                  "self.Cel.Count_ProteinElongationCompletedTotal")
            # # Number of Ribosomes released
            # Writer.PrintStVa("# of Ribosomes Released",
            #                  "self.Cel.Count_RibosomeReleased")
            # # Number of Ribosome bound to DNA
            # Writer.PrintStVa("# of Total Ribosomes Bound to mRNA",
            #                  "self.Cel.Count_RibosomeBound")
            # # Number of Ribosome freely floating
            # Writer.PrintStVa("# of Total Ribosomes Unbound Floating",
            #                  "self.Cel.Count_RibosomeUnbound")
            # Writer.BlankLine()

        # with Writer.Statement("def GetNascentProteinCounts(self):"):
        #     # CURRENTLY NOT USED - MAYBE FOR SUMMARY?
        #     Writer.ReduceSum("Count_ProteinsNascentPerProtein", "self.Cel.Len_ProteinsNascent", 1)
        #     Writer.ReduceSum("Count_ProteinsNascentTotal", "self.Cel.Len_ProteinsNascent", 0)
        #     Writer.BlankLine()