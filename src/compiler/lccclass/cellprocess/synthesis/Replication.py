# dna replication: DNA Polymerase III holozyme,  DNA Polymerase II (back-up)
'''
https://en.wikipedia.org/wiki/DNA_polymerase_III_holoenzyme
The replisome is composed of the following:

2 DNA Pol III enzymes, each comprising α, ε and θ subunits. (It has been proven that there is a third copy of Pol III at the replisome.[1])
the α subunit (encoded by the dnaE gene) has the polymerase activity.
the ε subunit (dnaQ) has 3'→5' exonuclease activity.
the θ subunit (holE) stimulates the ε subunit's proofreading.
2 β units (dnaN) which act as sliding DNA clamps, they keep the polymerase bound to the DNA.
2 τ units (dnaX) which act to dimerize two of the core enzymes (α, ε, and θ subunits).
1 γ unit (also dnaX) which acts as a clamp loader for the lagging strand Okazaki fragments, helping the two β subunits to form a unit and bind to DNA. The γ unit is made up of 5 γ subunits which include 3 γ subunits, 1 δ subunit (holA), and 1 δ' subunit (holB). The δ is involved in copying of the lagging strand.
Χ (holC) and Ψ (holD) which form a 1:1 complex and bind to γ or τ. X can also mediate the switch from RNA primer to DNA.[2]

1000 nucleotides per second
'''

'''
When Chr bp + 1,
+ 2 PPi
- 1 (use ACGT (frequency table for stoichiometry))

rate: + 1000 bp / 1 second
'''

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):

    # Initiators
    Idx_DndNTP = Comp.Master.ID2Idx_Master[Comp.Protein.Name2ID_Proteins['DndNTP']]
    Idx_DnaB = Comp.Master.ID2Idx_Master[Comp.Protein.Name2ID_Proteins['DnaB']]
    Idx_DnaC = Comp.Master.ID2Idx_Master[Comp.Protein.Name2ID_Proteins['DnaC']]

    # Molecule indices for Molecular IDs
    Idx_dNTPs = ProGen.BuildingBlockIdxs('dNTPs')
    Idx_PPi = Comp.Master.ID2Idx_Master['PPI[c]']

    # Chromosome indices
    Idx_Ch_Original = [0]
    Idx_Ch_Replicating = [1, 2, 3]

    # # Enzyme for leading strand
    # Idx_DNAPIII = Comp.Master.ID2Idx_Master[Comp.Complex.Name2ID_Complexes['DNA polymerase III, core enzyme']]
    #
    # # Enzyme for lagging strand
    # Idx_DNAPI = Comp.Master.ID2Idx_Master[Comp.Protein.Name2ID_Proteins["DNA polymerase I, 5' --> 3' polymerase, 5' --> 3'  and 3' --> 5' exonuclease"]]
    # Idx_DNAPrimase = Comp.Master.ID2Idx_Master[Comp.Protein.Name2ID_Proteins['DNA primase']]
    #
    # # helicase, ligase, topoisomerase, SSBs (single-strand binding proteins)
    #

    NUniq_Chromosomes = Comp.Chromosome.NUniq_Chromosomes
    NMax_Chromosomes = Comp.Chromosome.NMax_Chromosomes

    N_ReplicatingStrands = 4   # [Left_Leading, Left_Lagging, Right_Leading, Right_Lagging]

    Len_ChromosomeReplicatingMax = list()
    Len_ChromosomeReplicatingMaxA = round(Comp.Chromosome.Len_ChromosomesInGenome / 2)
    Len_ChromosomeReplicatingMaxB = Comp.Chromosome.Len_ChromosomesInGenome - Len_ChromosomeReplicatingMaxA
    for i in range(N_ReplicatingStrands):
        if i < N_ReplicatingStrands / 2:
            Len_ChromosomeReplicatingMax.append(Len_ChromosomeReplicatingMaxA)
        else:
            Len_ChromosomeReplicatingMax.append(Len_ChromosomeReplicatingMaxB)

    Rate_DNAReplication = 1000  # nt per second, accounting for both directions

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.Bool_ChromosomesReplicatingElongating", 0)
            Writer.Variable_("self.Bool_ChromosomesReplicatingOverElongated", 0)
            Writer.Variable_("self.Bool_ChromosomesReplicatingElongationCompleted", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Bin_ChromosomesReplicatingElongating", 0)
            Writer.Variable_("self.Bin_ChromosomesReplicatingOverElongated", 0)
            Writer.Variable_("self.Bin_ChromosomesReplicatingElongationCompleted", 0)
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.Cel.Idx_Ch_Original", Idx_Ch_Original)
            Writer.Variable_("self.Cel.Idx_Ch_Replicating", Idx_Ch_Replicating)

            Writer.Variable_("self.Cel.Idx_dNTPs", Idx_dNTPs)
            Writer.Variable_("self.Cel.Idx_PPi", Idx_PPi)

            Writer.Variable_("self.Cel.Rate_DNAReplication", Rate_DNAReplication)  # Not implemented yet
            Writer.Variable_("self.Cel.Rate_DNAReplication_Matrix", Rate_DNAReplication, [NUniq_Chromosomes, N_ReplicatingStrands])
            Writer.BlankLine()

            Writer.InitZeros("self.Cel.Len_ChromosomesReplicatingInitial", [NUniq_Chromosomes, N_ReplicatingStrands], 'int32')
            Writer.VarRepeat("self.Cel.Len_ChromosomesReplicatingMax", "self.Cel.Len_Chromosomes", N_ReplicatingStrands)
            Writer.Reshape__("self.Cel.Len_ChromosomesReplicatingMax", "self.Cel.Len_ChromosomesReplicatingMax",
                             [NUniq_Chromosomes, N_ReplicatingStrands])
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def ExecuteProcess(self):"):
            Writer.Statement("self.DNAReplicationInitiation()   # Not implemented")
            Writer.Statement("self.DNAElongation()")
            Writer.Statement("self.DNAReplicationTermination()")
            Writer.BlankLine()

        with Writer.Statement("def DNAReplicationInitiation(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def GetReplicatingChromosomeCounts(self):"):
            Writer.OperGathr("self.Cel.Count_ChromosomesReplicating_Matrix", "self.Cel.Counts", "self.Cel.Idx_Master_ChromosomesReplicating")
            Writer.OperRdSum("self.Cel.Count_ChromosomesReplicatingInitialTotal", "self.Cel.Count_ChromosomesReplicating_Matrix")
            Writer.BlankLine()

        with Writer.Statement("def DetermineReplicatingChromosomeIndices(self):"):
            # Get indices of Replicating Chromosome in length table by using tf.range? or repeat (Ones x repeat by Count then turn into boolean)
            # TO DO: use a single API for the routine to fill 1s from the left then zeros

            Writer.Variable_("ReplisomeArrayForAllDNAStrands", 0)
            with Writer.Statement("for i, Count in enumerate(self.Cel.Count_ChromosomesReplicating_Matrix):"):
                Writer.InitOnes_("Replisome_Present", "Count", 'int32')
                Writer.InitZeros("Replisome_Abscent", "%s - Count" % N_ReplicatingStrands, 'int32')
                Writer.OperCncat("ReplisomeArrayForCurrentDNAStrand", "Replisome_Present", "Replisome_Abscent")
                with Writer.Statement("if i == 0:"):
                    Writer.Statement("ReplisomeArrayForAllDNAStrands = ReplisomeArrayForCurrentDNAStrand")
                    Writer.Statement("continue")
                Writer.OperCncat("ReplisomeArrayForAllDNAStrands", "ReplisomeArrayForAllDNAStrands",
                                 "ReplisomeArrayForCurrentDNAStrand")
            Writer.Reshape__("ReplisomeArrayForAllDNAStrands", "ReplisomeArrayForAllDNAStrands", [-1, N_ReplicatingStrands])

            Writer.OperElGr_("self.Bool_ChromosomesReplicatingElongating", "ReplisomeArrayForAllDNAStrands", 0)
            Writer.Cast_____("self.Bin_ChromosomesReplicatingElongating", "self.Bool_ChromosomesReplicatingElongating", 'int32')
            Writer.OperRdSum("self.Cel.Count_ChromosomesReplicatingElongatingTotal",
                             "self.Bin_ChromosomesReplicatingElongating")  # For summary, the number of all elongating Chromosomes with Replisome
            Writer.BlankLine()

        with Writer.Statement("def GetReplicatingChromosomesIndicesInLengthMatrix(self):"):
            Writer.Statement("self.GetReplicatingChromosomeCounts()")
            Writer.Statement("self.DetermineReplicatingChromosomeIndices()")
            Writer.BlankLine()


        with Writer.Statement("def ElongateReplicatingChromosomesInLengthMatrix(self):"):
            # Rate Matrix * bin(boolean) + length table (overelongated not counted yet)
            Writer.OperElMul("self.Cel.Rate_ChromosomeElongation_Matrix", "self.Bin_ChromosomesReplicatingElongating",
                             "self.Cel.Rate_ChromosomeElongation")
            Writer.OperElAdd("self.Cel.Len_ChromosomesReplicatingElongated", "self.Cel.Len_ChromosomesReplicatingInitial",
                             "self.Cel.Rate_ChromosomeElongation_Matrix")
            Writer.BlankLine()

        with Writer.Statement("def DetermineOverElongatedReplicatingChromosomes(self):"):
            # compare with len MAX to get boolean matrix for overelongated (>= MAX)
            Writer.OperElGr_("self.Bool_ChromosomesReplicatingOverElongated", "self.Cel.Len_ChromosomesReplicatingElongated",
                             "self.Cel.Len_ChromosomesReplicatingMax")
            Writer.BoolToBin("self.Bin_ChromosomesReplicatingOverElongated", "self.Bool_ChromosomesReplicatingOverElongated")
            Writer.OperRdSum("self.Cel.Count_ChromosomesReplicatingOverElongatedTotal",
                             "self.Bin_ChromosomesReplicatingOverElongated")  # For summary, the number of completed/overelongated Chromosomes
            Writer.BlankLine()

        with Writer.Statement("def DetermineAmountOfOverElongateddNTPs(self):"):
            Writer.OperElMul("MaxLengthForOverElongated", "self.Bin_ChromosomesReplicatingOverElongated",
                             "self.Cel.Len_ChromosomesReplicatingMax")
            Writer.OperElMul("LengthForOverElongatedOnly", "self.Bin_ChromosomesReplicatingOverElongated",
                             "self.Cel.Len_ChromosomesReplicatingElongated")
            # Get the over-elongated length of Replicating Chromosomes
            Writer.OperElSub("self.Cel.Len_ChromosomesReplicatingOverElongated", "LengthForOverElongatedOnly",
                             "MaxLengthForOverElongated")
            Writer.OperRdSum("self.Cel.Count_dNTPsOverElongated",
                             "self.Cel.Len_ChromosomesReplicatingOverElongated")  # For summary, over consumped dNTPs to be corrected
            Writer.BlankLine()

        with Writer.Statement("def CorrectOverElongationInLengthMatrix(self):"):
            # Update Replicating Chromosome length by removing over elongated length
            Writer.OperElSub("self.Cel.Len_ChromosomesReplicatingAdjusted", "self.Cel.Len_ChromosomesReplicatingElongated",
                             "self.Cel.Len_ChromosomesReplicatingOverElongated")
            # now all completed elongation has its max value

            Writer.BlankLine()

        with Writer.Statement("def CorrectElongatedLengthForChromosomes(self):"):
            # Initial elongation rate
            Writer.OperElMul("Rate_ChromosomeElongation_Matrix_Initial", "self.Cel.Rate_ChromosomeElongation_Matrix",
                             "self.Bin_ChromosomesReplicatingElongating")
            # Adjusted elongation rate (initial - elongated)
            Writer.OperElSub("self.Cel.Rate_ChromosomeElongation_Matrix_Corrected", "Rate_ChromosomeElongation_Matrix_Initial",
                             "self.Cel.Len_ChromosomesReplicatingOverElongated")
            # Calculate the length of elongation per mRNA
            Writer.OperRdSum("self.Cel.Count_ChromosomeElongationLengthPerChromosome", "self.Cel.Rate_ChromosomeElongation_Matrix_Corrected",
                             1)  # For summary, total elongation length per mRNA
            # Calculate the total length of elongation
            Writer.OperRdSum("self.Cel.Count_ChromosomeElongationLengthTotal",
                             "self.Cel.Rate_ChromosomeElongation_Matrix_Corrected")  # For summary, total elongation length in this sim step
            Writer.BlankLine()

        with Writer.Statement("def ElongateDNAStrands(self):"):
            # Replicating Chromosomes will be elongated.
            Writer.Statement("self.GetReplicatingChromosomesIndicesInLengthMatrix()")
            Writer.Statement("self.ElongateReplicatingChromosomesInLengthMatrix()")
            Writer.Statement("self.DetermineOverElongatedReplicatingChromosomes()")
            Writer.Statement("self.DetermineAmountOfOverElongateddNTPs()")
            Writer.Statement("self.CorrectOverElongationInLengthMatrix()")
            Writer.Statement("self.CorrectElongatedLengthForChromosomes()")
            Writer.BlankLine()

        with Writer.Statement("def GetRawdNTPConsumption(self):"):
            # Prepare NT elongation length for each mRNA matrix for multiplication
            Writer.Reshape__("self.Cel.Count_ChromosomeElongationLengthPerChromosome", "self.Cel.Count_ChromosomeElongationLengthPerChromosome",
                             [-1, 1])
            Writer.Cast_____("self.Cel.Count_ChromosomeElongationLengthPerChromosome", "self.Cel.Count_ChromosomeElongationLengthPerChromosome",
                             'float32')
            # dNTP Frequency per mRNA * elongating DNAStrand (in NT count) per mRNA
            Writer.OperMXMul("dNTPConsumption_Raw", "self.Cel.Freq_dNTPsInChromosomes",
                             "self.Cel.Count_ChromosomeElongationLengthPerChromosome")
            Writer.ReturnVar("dNTPConsumption_Raw")
            Writer.BlankLine()

        with Writer.Statement("def GetRoundeddNTPConsumption(self, dNTPConsumption_Raw):"):
            Writer.RoundInt_("dNTPConsumption_Rounded", "dNTPConsumption_Raw")
            Writer.ReturnVar("dNTPConsumption_Rounded")
            Writer.BlankLine()

        with Writer.Statement("def CalculateDiscrepancy(self, dNTPConsumption_Rounded):"):
            # Determine the difference between rounded vs corrected elongation
            Writer.OperRdSum("Sum_Rounded", "dNTPConsumption_Rounded")
            Writer.Cast_____("Sum_Rounded", "Sum_Rounded", 'int32')
            Writer.OperRdSum("Sum_CorrectedElongationRate", "self.Cel.Rate_ChromosomeElongation_Matrix_Corrected")
            Writer.OperElSub("DeltaSum", "Sum_CorrectedElongationRate", "Sum_Rounded")
            Writer.ReturnVar("DeltaSum")
            Writer.BlankLine()

        with Writer.Statement("def AdjustdNTPConsumption(self, dNTPConsumption_Rounded, Discrepancy):"):
            # Get equal amount of dNTPs if discrepancy is greater than or equal to 4 or less than or equal to -4.
            Writer.OperElQuo("N_dNTPsets", "Discrepancy", "self.Cel.NUniq_dNTPs")
            Writer.VarRepeat("N_dNTPsets", "N_dNTPsets", "self.Cel.NUniq_dNTPs")
            Writer.OperElAdd("dNTPConsumption_MissingSet", "dNTPConsumption_Rounded", "N_dNTPsets")
            Writer.BlankLine()

            # Get random dNTP for the remainder (replace with weighted random dNTP based on dNTP frequency in all mRNAs)
            Writer.OperElRem("N_dNTPRemainder", "Discrepancy", "self.Cel.NUniq_dNTPs")
            Writer.Reshape__("N_dNTPRemainder", "N_dNTPRemainder", -1)
            Writer.InitZeros("dNTPConsumption_MissingRemainder", "self.Cel.NUniq_dNTPs", 'int32')
            Writer.RndNumUni("Idx_Remainder", "N_dNTPRemainder", "0", "self.Cel.NUniq_dNTPs")
            Writer.InitOnes_("OnesForRemainder", "N_dNTPRemainder", 'int32')
            Writer.OperScAdd("dNTPConsumption_MissingRemainder", "Idx_Remainder", "OnesForRemainder")
            Writer.BlankLine()

            # Calculate adjusted dNTP Consumption
            Writer.OperElAdd("dNTPConsumption_Adjusted", "dNTPConsumption_MissingSet", "dNTPConsumption_MissingRemainder")

            # Return the adjusted dNTP Consumption
            Writer.OperRdSum("TotaldNTPConsumption", "dNTPConsumption_Adjusted")
            Writer.AsrtElEq_("TotaldNTPConsumption", "self.Cel.Count_ChromosomeElongationLengthTotal")
            Writer.ReturnVar("dNTPConsumption_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def DeterminedNTPConsumption(self):"):
            Writer.Statement("dNTPConsumption_Raw = self.GetRawdNTPConsumption()")
            Writer.Statement("dNTPConsumption_Rounded = self.GetRoundeddNTPConsumption(dNTPConsumption_Raw)")
            Writer.Statement("Discrepancy = self.CalculateDiscrepancy(dNTPConsumption_Rounded)")
            Writer.Statement("dNTPConsumption_Adjusted = self.AdjustdNTPConsumption(dNTPConsumption_Rounded, Discrepancy)")
            Writer.ReturnVar("dNTPConsumption_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def ConsumedNTPs(self):"):
            Writer.Statement("self.Cel.Count_ChromosomeElongationdNTPConsumption = self.DeterminedNTPConsumption()")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_dNTPs, -self.Cel.Count_ChromosomeElongationdNTPConsumption)")
            Writer.BlankLine()

        with Writer.Statement("def ReleasePPi(self):"):
            Writer.Overwrite("self.Cel.Count_ChromosomeElongationPPiProduction", "self.Cel.Count_ChromosomeElongationLengthTotal")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_PPi, self.Cel.Count_ChromosomeElongationPPiProduction)")
            Writer.BlankLine()

        with Writer.Statement("def DNAStrandElongation(self):"):
            Writer.Statement("self.ElongateDNAStrands()")
            Writer.Statement("self.ConsumedNTPs()")
            Writer.Statement("self.ReleasePPi()")
            Writer.BlankLine()

        with Writer.Statement("def IdentifyCompletedChromosomeElongation(self):"):
            # identify all NMAX values by comparison (boolean)
            Writer.OperElEq_("self.Bool_ChromosomesReplicatingElongationCompleted", "self.Cel.Len_ChromosomesReplicatingAdjusted",
                             "self.Cel.Len_ChromosomesReplicatingMax")
            Writer.BoolToBin("self.Bin_ChromosomesReplicatingElongationCompleted",
                             "self.Bool_ChromosomesReplicatingElongationCompleted")
            Writer.BlankLine()

        with Writer.Statement("def CountCompletedChromosomeElongation(self):"):
            # Count all completed per Chromosome and total
            Writer.OperRdSum("self.Cel.Count_ChromosomeElongationCompletedPerChromosome",
                             "self.Bin_ChromosomesReplicatingElongationCompleted", 1)
            Writer.OperRdSum("self.Cel.Count_ChromosomeElongationCompletedTotal",
                             "self.Bin_ChromosomesReplicatingElongationCompleted")  # For summary, total completed elongation
            Writer.BlankLine()

        with Writer.Statement("def ResetLengthOfCompletedReplicatingChromosomes(self):"):
            # Reset all NMAX lengths to zero by subtracting max length from completed
            Writer.OperElMul("self.Cel.Len_ChromosomesReplicatingCompleted", "self.Cel.Len_ChromosomesReplicatingMax",
                             "self.Bin_ChromosomesReplicatingElongationCompleted")
            Writer.OperElSub("self.Cel.Len_ChromosomesReplicatingFinal", "self.Cel.Len_ChromosomesReplicatingAdjusted",
                             "self.Cel.Len_ChromosomesReplicatingCompleted")

            # Check reset status
            Writer.OperElMul("CheckCompletionReset", "self.Cel.Len_ChromosomesReplicatingFinal",
                             "self.Bin_ChromosomesReplicatingElongationCompleted")
            Writer.AsrtElEq_("CheckCompletionReset", 0)
            Writer.BlankLine()

        with Writer.Statement("def GetIndicesOfChromosomesCompletedElongation(self):"):
            # Get indices of Chromosomes completed elongation. Currently not used
            Writer.GetIdxGr_("Idx_ChromosomeElongationCompleted", "self.Cel.Count_ChromosomeElongationCompletedPerChromosome", 0)
            Writer.OperGathr("self.Cel.Idx_ChromosomeElongationCompleted", "self.Cel.Idx_ChromosomesReplicating",
                             "Idx_ChromosomeElongationCompleted")
            Writer.BlankLine()

        with Writer.Statement("def ReleaseReplisome(self):"):
            Writer.Overwrite("self.Cel.Count_ReplisomeReleased", "self.Cel.Count_ChromosomeElongationCompletedTotal")
            Writer.OperElSub("self.Cel.Count_ReplisomeBound", "self.Cel.Count_ReplisomeBound", "self.Cel.Count_ReplisomeReleased")
            Writer.OperElAdd("self.Cel.Count_ReplisomeUnbound", "self.Cel.Count_ReplisomeUnbound", "self.Cel.Count_ReplisomeReleased")
            Writer.OperElAdd("SumOfBoundUnbound", "self.Cel.Count_ReplisomeBound", "self.Cel.Count_ReplisomeUnbound")
            Writer.AsrtElEq_("self.Cel.Count_Replisome", "SumOfBoundUnbound")
            Writer.BlankLine()

        with Writer.Statement("def DeductCountOfReplicatingChromosomes(self):"):
            # Deduct count of Replicating Chromosomes completed elongation
            Writer.Statement(
                "self.AddToDeltaCounts(self.Cel.Idx_Master_ChromosomesReplicating, -self.Cel.Count_ChromosomeElongationCompletedPerChromosome)")
            Writer.BlankLine()

        with Writer.Statement("def IncrementCountOfChromosomes(self):"):
            # Increment count of Chromosomes completed elongation
            Writer.Statement(
                "self.AddToDeltaCounts(self.Cel.Idx_Master_Chromosomes, self.Cel.Count_ChromosomeElongationCompletedPerChromosome)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateLengthOfReplicatingChromosomesMatrix(self):"):
            # Overwrite self.Cel.Len_ChromosomesReplicatingInitial with self.Cel.Len_ChromosomesReplicatingFinal
            Writer.Overwrite("self.Cel.Len_ChromosomesReplicatingInitial", "self.Cel.Len_ChromosomesReplicatingFinal")
            Writer.BlankLine()

        with Writer.Statement("def DNAReplicationTermination(self):"):
            Writer.Statement("self.IdentifyCompletedChromosomeElongation()")
            Writer.Statement("self.CountCompletedChromosomeElongation()")
            Writer.Statement("self.ResetLengthOfCompletedReplicatingChromosomes()")
            # Writer.Statement("self.GetIndicesOfChromosomesCompletedElongation()")
            Writer.Statement("self.ReleaseReplisome()")
            Writer.Statement("self.DeductCountOfReplicatingChromosomes()")
            Writer.Statement("self.IncrementCountOfChromosomes()")
            Writer.Statement("self.UpdateLengthOfReplicatingChromosomesMatrix()")
            Writer.BlankLine()

        with Writer.Statement("def Debugging_Chromosome_Indices(self, Str):"):
            Writer.Statement("print(Str)")
            Writer.BlankLine()

            Writer.Comment__("Chromosome length Indices")
            Writer.PrintStrg("Chromosome length Indices")

            Writer.PrintVar_("self.Cel.Idx_RndChromosomesReplicating")

            Writer.PrtIdxGr_("self.Cel.Len_ChromosomesReplicatingInitial", 0)
            Writer.PrtIdxGr_("self.Cel.Len_ChromosomesReplicatingElongated", 0)
            Writer.PrtIdxGr_("self.Cel.Len_ChromosomesReplicatingOverElongated", 0)

            Writer.PrtIdxGr_("self.Cel.Rate_ChromosomeElongation_Matrix_Corrected", 0)

            Writer.PrtIdxGr_("self.Bin_ChromosomesReplicatingElongating", 0)
            Writer.PrtIdxGr_("self.Bin_ChromosomesReplicatingOverElongated", 0)
            Writer.PrtIdxGr_("self.Bin_ChromosomesReplicatingElongationCompleted", 0)

            Writer.BlankLine()

        with Writer.Statement("def Debugging_Count_ChromosomesReplicating(self, Str):"):
            Writer.Statement("print(Str)")
            Writer.BlankLine()

            Writer.Comment__("Replicating Chromosomes")
            Writer.PrintStrg("Replicating Chromosomes")

            Writer.OperRdSum("Reduced_Sum_Of_self_Cel_Count_ChromosomesReplicating", "self.Cel.Count_ChromosomesReplicating")
            Writer.PrintVaVa("Reduced_Sum_Of_self_Cel_Count_ChromosomesReplicating")
            Writer.BlankLine()

            Writer.OperGathr("Counts_ChromosomesReplicating", "self.Cel.Counts", "self.Cel.Idx_Master_ChromosomesReplicating")
            Writer.OperRdSum("Reduced_Sum_Of_Counts_ChromosomesReplicating", "Counts_ChromosomesReplicating")
            Writer.PrintVaVa("Reduced_Sum_Of_Counts_ChromosomesReplicating")
            Writer.BlankLine()

            Writer.OperGathr("DeltaCounts_ChromosomesReplicating", "self.Cel.DeltaCounts", "self.Cel.Idx_Master_ChromosomesReplicating")
            Writer.OperRdSum("Reduced_Sum_Of_DeltaCounts_ChromosomesReplicating", "DeltaCounts_ChromosomesReplicating")
            Writer.PrintVaVa("Reduced_Sum_Of_DeltaCounts_ChromosomesReplicating")
            Writer.BlankLine()

            Writer.Comment__("Chromosomes (full)")
            Writer.PrintStrg("Chromosomes (full)")

            Writer.OperRdSum("Reduced_Sum_Of_self_Cel_Count_Chromosomes", "self.Cel.Count_Chromosomes")
            Writer.PrintVaVa("Reduced_Sum_Of_self_Cel_Count_Chromosomes")
            Writer.BlankLine()

            Writer.OperGathr("Counts_Chromosomes", "self.Cel.Counts", "self.Cel.Idx_Master_Chromosomes")
            Writer.OperRdSum("Reduced_Sum_Of_Counts_Chromosomes", "Counts_Chromosomes")
            Writer.PrintVaVa("Reduced_Sum_Of_Counts_Chromosomes")
            Writer.BlankLine()

            Writer.OperGathr("DeltaCounts_Chromosomes", "self.Cel.DeltaCounts", "self.Cel.Idx_Master_Chromosomes")
            Writer.OperRdSum("Reduced_Sum_Of_DeltaCounts_Chromosomes", "DeltaCounts_Chromosomes")
            Writer.PrintVaVa("Reduced_Sum_Of_DeltaCounts_Chromosomes")
            Writer.BlankLine()

        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.PrintStrg("===== DNAReplication Initiation ===== ")
            # Number of Total Replisome
            Writer.PrintStVa("# of Total Replisomes",
                             "self.Cel.Count_Replisome")
            # Number of Active Replisome
            Writer.PrintStVa("# of active Replisomes",
                             "self.Cel.Count_ReplisomeActive")
            # Number of Replisome binding
            Writer.PrintStVa("# of Replisomes available that can bind a promoter",
                             "self.Cel.Count_ReplisomeActiveCanBind")
            # Number of new Replisome binding in this step (= new Replicating Chromosomes to elongate next step)
            Writer.PrintStVa("# of Replisomes that newly bind a promoter in this step",
                             "self.Cel.Count_ReplisomeWillBind")
            Writer.BlankLine()

            Writer.PrintStrg("===== DNAStrand Elongation ===== ")
            # Number of Replicating Chromosomes elongated
            Writer.PrintStVa("# of all Replicating Chromosomes elongating",
                             "self.Cel.Count_ChromosomesReplicatingElongatingTotal")
            # Total elongation length of Chromosomes
            Writer.PrintStVa("Total elongation length of Chromosomes (nt)",
                             "self.Cel.Count_ChromosomeElongationLengthTotal")
            # Total dNTP consumption and PPi production
            Writer.PrintStVa("Total dNTP consumption [ATP, CTP, GTP, UTP]",
                             "self.Cel.Count_ChromosomeElongationdNTPConsumption")
            Writer.PrintStVa("Total PPi production",
                             "self.Cel.Count_ChromosomeElongationPPiProduction")
            Writer.BlankLine()

            Writer.PrintStrg("===== DNAReplication Termination ===== ")
            # Number of Chromosome elongation Completed
            Writer.PrintStVa("# of Chromosome Elongation Completed",
                             "self.Cel.Count_ChromosomeElongationCompletedTotal")
            # Number of Replisomes released
            Writer.PrintStVa("# of Replisomes released",
                             "self.Cel.Count_ReplisomeReleased")
            # Number of Replisome bound to DNA
            Writer.PrintStVa("# of total Replisomes bound to DNA",
                             "self.Cel.Count_ReplisomeBound")
            # Number of Replisome freely floating
            Writer.PrintStVa("# of total Replisomes unbound floating",
                             "self.Cel.Count_ReplisomeUnbound")
            Writer.BlankLine()


# def SetUpReactions(ProGen):
#     Reactions = list()
#
#     Reaction = dict()
#     # Reaction No.1
#
#     RXNType = 'Polymerization'
#     RXNEquation = '2 dNTP -> 1 ChromosomeSize + 2 PPi'
#     RXNRate = '1000 +- 100 events per second'
#     RXNTrigger = 'DNA polymerase III, core enzyme >= 4'
#
#     # The final product would be the following:
#
#     # Type = 'Polymerization' or 'Biochemical Reaction'
#
#     # Type 'Polymerization'
#     # Stoich_MolIDs = ['MolID #1', 'MolID #2', etc],
#     #                   where MolIDs of reactants and products must exist in 'Cel.ID_Master'
#     # Stoich_Coeffs = [Coeff for Mol #1, Coeff for Mol #2, etc],
#     #                   where negative and positive integers are used for reactants and products, respectively
#     # Rate_Min = integer
#     # Rate_Max = integer
#     # Rate_Distribution = 'Normal', 'Uniform',
#     #
#     # Trigger_MolIDs = ['MolID #1', 'MolID #2', etc],
#     #                   where conditional presence of MolIDs (or environmental condition to be implemented)
#     # Trigger_Thresholds = ['string of integer count for MolID #1', etc],
#     #                   where there many be many triggers to satisfy
#
#     # DNAPI = ProGen.Comp.Protein.Name2ID_Proteins["DNA polymerase I, 5' --> 3' polymerase, 5' --> 3'  and 3' --> 5' exonuclease"]
#     # DNAPrimase = ProGen.Comp.Protein.Name2ID_Proteins['DNA primase']
#     DNAPIII = ProGen.Comp.Complex.Name2ID_Complexes['DNA polymerase III, core enzyme']
#
#     # Temporary parameters
#     ChromosomeNumber = 1
#     ChromosomeRepNumber = 1
#
#     # This is for leading
#
#     for Attribute in ProGen.Comp.Chromosome.ReplicatingChromosomeAttributes:
#         Type = 'Polymerization'
#
#         ChromosomeState = 'Partial_' + Attribute
#         ChromosomeID = 'Ch%d_Rep%d_%s' % (ChromosomeNumber, ChromosomeRepNumber, ChromosomeState)
#
#         Stoich_MolIDs = ['dNTP', ChromosomeID, 'PPI[c]']
#         Stoich_Coeffs = [-2, 1, 2]
#
#         Rate_Mean = 1000  # basepairs per second, accounting for both directions
#         Rate_SD = 100
#         Rate_UnitTime = 'Second'
#
#         # Trigger_MolIDs = [DNAPIII, DNAPrimase, DNAPI]  # 'DNA polymerase III, core enzyme'
#         # Trigger_Thresholds = ['4', '4', '4']
#         # Trigger_Conditions = ['>=', '>=', '>=']  # Greater than or equal to
#
#         # Generate a reaction dictionary with above inputs
#         Reaction['Type'] = Type
#         Reaction['Stoichiometry'] = [Stoich_MolIDs, Stoich_Coeffs]
#         Reaction['Rate'] = [Rate_Mean, Rate_SD, Rate_UnitTime]
#         # Reaction['Trigger'] = [Trigger_MolIDs, Trigger_Thresholds, Trigger_Conditions]
#
#         Reaction_SetUp = ProGen.SetUpReaction(Reaction)
#         # MolIdxs have been added to Stoichiometry and Trigger variables
#         Reactions.append(Reaction_SetUp)
#
#     return Reactions



