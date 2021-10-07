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

'''
Review Article
Published: 04 June 2019
The bacterial cell cycle, chromosome inheritance and cell growth
Rodrigo Reyes-Lamothe & David J. Sherratt 
Nature Reviews Microbiology volume 17, pages467–478 (2019)
https://www.nature.com/articles/s41579-019-0212-7


initiation of E. coli DNA replication occurs at a fixed cell volume per origin (oriC), 
independent of birth size and growth rate in individual cells20,21,22, in line with Donachie’s original paper.14

14. Donachie, W. D. Relationship between cell size and time of initiation of DNA replication. Nature 219, 1077–1079 (1968).
20. Si, F. et al. Invariance of initiation mass and predictability of cell size in Escherichia coli. Curr. Biol. 27, 1278–1287 (2017).Return to ref 20 in article
21. Zheng, H. et al. Interrogating the Escherichia coli cell cycle by cell dimension perturbations. Proc. Natl Acad. Sci. USA 113, 15000–15005 (2016).
22. Wallden, M., Fange, D., Lundius, E. G., Baltekin, O. & Elf, J. The synchronization of replication and division cycles in individual E. coli cells. Cell 166, 729–739 (2016).

First, initiation occurs at a fixed volume per origin, marking the start of the cell cycle. 

Second, the relation between volume increase and initiating capacity during growth makes models for initiation incompatible with a mechanism 
whereby regulation occurs only through changes in the copy number of proteins, including the initiator.
'''

import numpy as np

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):

    # Initiators
    Idx_DnaA = Comp.Master.ID2Idx_Master[Comp.Protein.ID2ID_Gene2Protein[Comp.Gene.Sym2ID_Genes['dnaA']]]
    Idx_DnaB = Comp.Master.ID2Idx_Master[Comp.Protein.ID2ID_Gene2Protein[Comp.Gene.Sym2ID_Genes['dnaB']]]
    Idx_DnaC = Comp.Master.ID2Idx_Master[Comp.Protein.ID2ID_Gene2Protein[Comp.Gene.Sym2ID_Genes['dnaC']]]

    # Molecule indices for Molecular IDs

    # Chromosome indices
    Idx_Ch_Original = [0]
    Idx_Ch_Replicating = [1, 2, 3]

    Idx_DNAStrand_Proxy_Left = [0]
    Idx_DNAStrand_Proxy_Right = [2]

    # # Enzyme for leading strand
    # Idx_DNAPIII = Comp.Master.ID2Idx_Master[Comp.Complex.Name2ID_Complexes['DNA polymerase III, core enzyme']]
    #
    # # Enzyme for lagging strand
    # Idx_DNAPI = Comp.Master.ID2Idx_Master[Comp.Protein.Name2ID_Proteins["DNA polymerase I, 5' --> 3' polymerase, 5' --> 3'  and 3' --> 5' exonuclease"]]
    # Idx_DNAPrimase = Comp.Master.ID2Idx_Master[Comp.Protein.Name2ID_Proteins['DNA primase']]
    #
    # # helicase, ligase, topoisomerase, SSBs (single-strand binding proteins)
    #

    Bin_ChromosomesReplicating_Slow = [1, 0, 0]
    Bin_ChromosomesReplicating_Fast = [1, 1, 1]

    Conc_DnaAATP_Slow = 20 * 10**-9   # (10~30nM)
    Conc_DnaAATP_Fast = 80 * 10**-9   # (60~100nM)

    NUniq_Chromosomes = Comp.Chromosome.NUniq_ChromosomesInGenome
    NMax_Chromosomes = Comp.Chromosome.NMax_Chromosomes
    N_ChromosomesReplicating = NMax_Chromosomes - NUniq_Chromosomes   # All - original set

    N_ReplicatingStrands = 4   # [Left_Leading, Left_Lagging, Right_Leading, Right_Lagging]

    # TODO: Expand this part to handle more than a single chromosome
    Len_ChromosomesReplicatingMax = list()

    Len_Chromosome = int(Comp.Chromosome.Len_ChromosomesInGenome)
    Len_ChromosomesReplicatingMaxL = round(Len_Chromosome / 2)
    Len_ChromosomesReplicatingMaxR = Len_Chromosome - Len_ChromosomesReplicatingMaxL
    for i in range(N_ReplicatingStrands):
        if i < N_ReplicatingStrands / 2:
            Len_ChromosomesReplicatingMax.append(Len_ChromosomesReplicatingMaxL)
        else:
            Len_ChromosomesReplicatingMax.append(Len_ChromosomesReplicatingMaxR)

    Rate_DNAReplication = 1000  # nt per second, accounting for both directions

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Function_("Init_ProcessSpecificVariables"):
            Writer.Variable_("self.Idx_DNAStrand_Leftward", 0)
            Writer.Variable_("self.Idx_DNAStrand_Rightward", 0)
            Writer.Variable_("self.Idx_GenesReplicated", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Conc_DnaAATP_Slow", 0)
            Writer.Variable_("self.Conc_DnaAATP_Fast", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Bin_ChromosomesReplicating_Slow", 0)
            Writer.Variable_("self.Bin_ChromosomesReplicating_Fast", 0)
            Writer.BlankLine()

            Writer.Variable_("self.N_GenesReplicated", 0)
            Writer.Variable_("self.Rate_DNAReplication", 0)
            Writer.Variable_("self.Len_ChromosomesReplicatingInitial", 0)
            Writer.Variable_("self.Count_dNTPConsumption", 0)
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Function_("SetUp_ProcessSpecificVariables"):
            Writer.Variable_("self.Cel.Idx_DnaA", Idx_DnaA)
            Writer.Variable_("self.Cel.Idx_DnaB", Idx_DnaB)
            Writer.Variable_("self.Cel.Idx_DnaC", Idx_DnaC)

            Writer.Variable_("self.Cel.Idx_Ch_Original", Idx_Ch_Original)
            Writer.Variable_("self.Cel.Idx_Ch_Replicating", Idx_Ch_Replicating)

            Writer.Variable_("self.Idx_DNAStrand_Proxy_Left", Idx_DNAStrand_Proxy_Left)
            Writer.Variable_("self.Idx_DNAStrand_Proxy_Right", Idx_DNAStrand_Proxy_Right)

            Writer.Variable_("self.Rate_DNAReplication", Rate_DNAReplication)
            Writer.BlankLine()

            Writer.Variable_("self.Conc_DnaAATP_Slow", Conc_DnaAATP_Slow)
            Writer.Variable_("self.Conc_DnaAATP_Fast", Conc_DnaAATP_Fast)
            Writer.BlankLine()

            # Presets
            Writer.Variable_("self.Bin_ChromosomesReplicating_Slow", Bin_ChromosomesReplicating_Slow)
            Writer.Variable_("self.Bin_ChromosomesReplicating_Fast", Bin_ChromosomesReplicating_Fast)
            Writer.BlankLine()

            Writer.VarFill__("self.Cel.Len_ChromosomesReplicating", [N_ChromosomesReplicating, N_ReplicatingStrands], -1)
            Writer.Cast_____("self.Cel.Len_ChromosomesReplicating", "self.Cel.Len_ChromosomesReplicating", 'int32')
            Writer.Variable_("self.Cel.Len_ChromosomesReplicatingMax", Len_ChromosomesReplicatingMax)
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Function_("ExecuteProcess"):
            Writer.Statement("self.Initiation()")
            Writer.Statement("self.Elongation()")
            Writer.Statement("self.Termination()")
            Writer.BlankLine()

        with Writer.Function_("DetermineChromosomeReplicatingState"):
            # # Get DnaA counts
            # Writer.Statement("Count_DnaA = self.GetCounts(self.Cel.Idx_DnaA)")
            # Writer.Statement("Count_DnaB = self.GetCounts(self.Cel.Idx_DnaB)")
            # Writer.Statement("Count_DnaC = self.GetCounts(self.Cel.Idx_DnaC)")
            #
            # # Calculate Threshold values for DnaA count for slow vs fast replicating state
            # Writer.Statement("Count_DnaAATP_Slow = self.ConvertConcToCount(self.Conc_DnaAATP_Slow)")
            # Writer.Statement("Count_DnaAATP_Fast = self.ConvertConcToCount(self.Conc_DnaAATP_Fast)")
            #
            # Writer.ConvToBin("Bin_ChromosomesReplicating_Slow", "Count_DnaA", ">", "Count_DnaAATP_Slow")
            # Writer.ConvToBin("Bin_ChromosomesReplicating_Fast", "Count_DnaA", ">", "Count_DnaAATP_Fast")
            #
            # Writer.Concat___("Bin_ChromosomesReplicating", "Bin_ChromosomesReplicating_Slow", "Bin_ChromosomesReplicating_Fast")
            # Writer.Concat___("Bin_ChromosomesReplicating", "Bin_ChromosomesReplicating", "Bin_ChromosomesReplicating_Fast")
            # Writer.BlankLine()

            # TODO: Currently set to slow growth only
            Writer.Overwrite("Bin_ChromosomesReplicating", "self.Bin_ChromosomesReplicating_Slow")
            Writer.ReturnVar("Bin_ChromosomesReplicating")
            Writer.BlankLine()


        with Writer.Function_("Initiation"):
            # TODO: Replisome regulated replication initiation mechanism needs to be implemented
            Writer.Statement("Bin_ChromosomesReplicating = self.DetermineChromosomeReplicatingState()")
            Writer.Transpose("Bin_ChromosomesReplicating", "Bin_ChromosomesReplicating")

            # Check and prepare Len_ChromosomesReplicating to match the replicating state
            Writer.VarFill__("TurnOnReplication", [N_ChromosomesReplicating, N_ReplicatingStrands], -1)
            Writer.Add______("TurnOnReplication", "TurnOnReplication", "Bin_ChromosomesReplicating")

            Writer.ConvToBin("AddToReplication", "self.Cel.Len_ChromosomesReplicating", "<", "TurnOnReplication")
            Writer.Add______("self.Cel.Len_ChromosomesReplicating", "self.Cel.Len_ChromosomesReplicating", "AddToReplication")
            Writer.BlankLine()

        # with Writer.Statement("def GetReplicatingChromosomeCounts(self):"):
        #     Writer.Gather___("self.Cel.Count_ChromosomesReplicating_Matrix", "self.Cel.Counts", "self.Cel.Idx_Ch_Replicating")
        #     Writer.ReduceSum("self.Cel.Count_ChromosomesReplicatingInitialTotal", "self.Cel.Count_ChromosomesReplicating_Matrix")   # for summary
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def GetReplicatingChromosomesIndicesInLengthMatrix(self):"):
        #     Writer.Statement("self.GetReplicatingChromosomeCounts()")
        #     Writer.Statement("self.DetermineReplicatingChromosomeIndices()")
        #     Writer.BlankLine()


        with Writer.Function_("UpdateChromosomesReplicatingLengths", "Len_ChromosomesReplicating_Update"):
            Writer.Overwrite("self.Cel.Len_ChromosomesReplicating", "Len_ChromosomesReplicating_Update")
            Writer.BlankLine()

        with Writer.Function_("DeterminedNTPConsumption", "Len_ToElongate"):
            Writer.Statement("dNTPConsumption = self.DetermineAmountOfBuildingBlocks(Len_ToElongate, self.Cel.Freq_NTsInChromosomesReplicating)")

            # Return the adjusted dNTP Consumption
            Writer.ReduceSum("TotaldNTPConsumptionFinal", "dNTPConsumption")
            Writer.ReduceSum("TotaldNTPConsumptionInitial", "Len_ToElongate")
            Writer.AsrtEq___("TotaldNTPConsumptionFinal", "TotaldNTPConsumptionInitial")
            Writer.ReturnVar("dNTPConsumption")
            Writer.BlankLine()

        with Writer.Function_("ConsumedNTPs", "Len_ToElongate"):
            Writer.Statement("dNTPConsumption = self.DeterminedNTPConsumption(Len_ToElongate)")
            Writer.Statement("self.Count_dNTPConsumption = dNTPConsumption")
            Writer.ReduceSum("self.Count_dNTPConsumptionTotal", "dNTPConsumption")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_dNTPs, -dNTPConsumption)")
            Writer.BlankLine()

        with Writer.Function_("ReleasePPi", "Count_TotalLengthOfElongation"):
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_PPi, Count_TotalLengthOfElongation)")
            Writer.BlankLine()

        with Writer.Function_("UpdateByproducts", "Len_ToElongate"):
            Writer.ReduceSum("Count_TotalLengthOfElongationPerChromosome", "Len_ToElongate", 1)
            Writer.Statement("self.ConsumedNTPs(Count_TotalLengthOfElongationPerChromosome)")
            Writer.ReduceSum("Count_TotalLengthOfElongationTotal", "Len_ToElongate")
            Writer.Statement("self.ReleasePPi(Count_TotalLengthOfElongationTotal)")
            Writer.BlankLine()

        with Writer.Function_("IncrementGeneCount", "Len_ChromosomesReplicatingInitial", "Len_ChromosomesReplicatingFinal"):
            # Determine Start and End positions.
            # TODO: Expand this part to handle more than a single chromosome
            Writer.Gather___("Len_Initial_Left", "Len_ChromosomesReplicatingInitial", "self.Idx_DNAStrand_Proxy_Left")
            Writer.Gather___("Len_Initial_Right", "Len_ChromosomesReplicatingInitial", "self.Idx_DNAStrand_Proxy_Right")
            Writer.BlankLine()
            Writer.Gather___("Len_Final_Left", "Len_ChromosomesReplicatingFinal", "self.Idx_DNAStrand_Proxy_Left")
            Writer.Gather___("Len_Final_Right", "Len_ChromosomesReplicatingFinal", "self.Idx_DNAStrand_Proxy_Right")
            Writer.BlankLine()

            # Check if reindexed gene coordinates is between Start and End positions.
            Writer.GreaterEq("Bool_Initial_Left", "self.Cel.Coord_Genes_Reindexed_Leftward", "Len_Initial_Left")
            Writer.Less_____("Bool_Final_Left", "self.Cel.Coord_Genes_Reindexed_Leftward", "Len_Final_Left")
            Writer.LogicAnd_("Bool_Left", "Bool_Initial_Left", "Bool_Final_Left")
            Writer.BlankLine()
            Writer.GreaterEq("Bool_Initial_Right", "self.Cel.Coord_Genes_Reindexed_Rightward", "Len_Initial_Right")
            Writer.Less_____("Bool_Final_Right", "self.Cel.Coord_Genes_Reindexed_Rightward", "Len_Final_Right")
            Writer.LogicAnd_("Bool_Right", "Bool_Initial_Right", "Bool_Final_Right")
            Writer.BlankLine()
            Writer.LogicOr__("Bool_GenesReplicated", "Bool_Left", "Bool_Right")
            Writer.Reshape__("Bool_GenesReplicated", "Bool_GenesReplicated", -1)
            Writer.GenIdx___("Idx_GenesReplicated", "Bool_GenesReplicated")
            Writer.BlankLine()

            # Convert local gene indexing to master
            Writer.Gather___("self.Idx_GenesReplicated", "self.Cel.Idx_Master_Genes", "Idx_GenesReplicated")

            # Count the number of genes replicated in the current step
            Writer.NonZeros_("self.N_GenesReplicated", "Bool_GenesReplicated")

            # Get an array of ones to match the number of genes replicated
            Writer.InitOnes_("Ones", "self.N_GenesReplicated", "int32")

            # Increment Gene count in the delta
            Writer.Statement("self.AddToDeltaCounts(self.Idx_GenesReplicated, Ones)")
            Writer.BlankLine()

        with Writer.Function_("Elongation"):
            Writer.Overwrite("self.Len_ChromosomesReplicatingInitial", "self.Cel.Len_ChromosomesReplicating")
            Writer.Statement("Len_ChromosomesReplicating_Elongated = self.DetermineAmountOfElongation(self.Cel.Len_ChromosomesReplicating, self.Rate_DNAReplication, self.Cel.Len_ChromosomesReplicatingMax)")
            Writer.Statement("Len_ToElongate = Len_ChromosomesReplicating_Elongated - self.Cel.Len_ChromosomesReplicating")
            Writer.Statement("self.UpdateChromosomesReplicatingLengths(Len_ChromosomesReplicating_Elongated)")
            Writer.Statement("self.UpdateByproducts(Len_ToElongate)")
            Writer.Statement("self.IncrementGeneCount(self.Len_ChromosomesReplicatingInitial, self.Cel.Len_ChromosomesReplicating)")
            Writer.BlankLine()

        with Writer.Function_("DetermineFullyReplicatedChromosomes"):
            # identify all NMAX values by comparison to get a binary table for elongation completion
            Writer.ConvToBin("Bin_ChromosomesReplicatingElongationCompleted", "self.Cel.Len_ChromosomesReplicating", "==", "self.Cel.Len_ChromosomesReplicatingMax")
            # Count all completed per Chromosome and total
            Writer.ReduceSum("Count_DNAStrandElongationCompletedTotal",
                             "Bin_ChromosomesReplicatingElongationCompleted")  # For summary, total completed elongation
            Writer.ReturnVar("Count_DNAStrandElongationCompletedTotal")
            Writer.BlankLine()

        # with Writer.Statement("def ResetLengthOfCompletedReplicatingChromosomes(self):"):
        #     # Reset all NMAX lengths to zero by subtracting max length from completed
        #     Writer.Multiply_("self.Cel.Len_ChromosomesReplicatingCompleted", "self.Cel.Len_ChromosomesReplicatingMax",
        #                      "self.Bin_ChromosomesReplicatingElongationCompleted")
        #     Writer.Subtract_("self.Cel.Len_ChromosomesReplicatingFinal", "self.Cel.Len_ChromosomesReplicatingAdjusted",
        #                      "self.Cel.Len_ChromosomesReplicatingCompleted")
        #
        #     # Check reset status
        #     Writer.Multiply_("CheckCompletionReset", "self.Cel.Len_ChromosomesReplicatingFinal",
        #                      "self.Bin_ChromosomesReplicatingElongationCompleted")
        #     Writer.AsrtEq___("CheckCompletionReset", 0)
        #     Writer.BlankLine()

        # with Writer.Statement("def GetIndicesOfChromosomesCompletedElongation(self):"):
        #     # Get indices of Chromosomes completed elongation. Currently not used
        #     Writer.GetIdxGr_("Idx_ChromosomeElongationCompleted", "self.Cel.Count_DNAStrandElongationCompletedPerChromosome", 0)
        #     Writer.Gather___("self.Cel.Idx_ChromosomeElongationCompleted", "self.Cel.Idx_ChromosomesReplicating",
        #                      "Idx_ChromosomeElongationCompleted")
        #     Writer.BlankLine()

        with Writer.Function_("ReleaseReplisome", "Count_DNAStrandsFullReplicatedTotal"):
            # Writer.Overwrite("self.Cel.Count_ReplisomeReleased", Count_DNAStrandsFullReplicatedTotal)
            # Writer.Subtract_("self.Cel.Count_ReplisomeBound", "self.Cel.Count_ReplisomeBound", "self.Cel.Count_ReplisomeReleased")
            # Writer.Add______("self.Cel.Count_ReplisomeUnbound", "self.Cel.Count_ReplisomeUnbound", "self.Cel.Count_ReplisomeReleased")
            # Writer.Add______("SumOfBoundUnbound", "self.Cel.Count_ReplisomeBound", "self.Cel.Count_ReplisomeUnbound")
            # Writer.AsrtEq___("self.Cel.Count_Replisome", "SumOfBoundUnbound")
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Function_("Termination"):
            Writer.ConvToBin("Bin_DNAStrandsFullReplicated", "self.Cel.Len_ChromosomesReplicating", "==", "self.Cel.Len_ChromosomesReplicatingMax")
            Writer.ReduceSum("Count_DNAStrandsFullReplicatedPerChromosome", "Bin_DNAStrandsFullReplicated", 1)
            Writer.ConvToBin("Bin_ChromosomesFullReplicated", "Count_DNAStrandsFullReplicatedPerChromosome", "==", N_ReplicatingStrands)
            Writer.ReduceSum("Count_ChromosomesFullReplicatedTotal", "Bin_ChromosomesFullReplicated")
            Writer.BlankLine()

            Writer.ReduceSum("Count_DNAStrandsFullReplicatedTotal", "Bin_DNAStrandsFullReplicated")
            Writer.Statement("self.ReleaseReplisome(Count_DNAStrandsFullReplicatedTotal)   # To be implemented")
            # The length of the completed replicating chromosomes is handled by the cell division
            Writer.BlankLine()

        with Writer.Function_("DeterminePercentReplicationCompletion"):
            Writer.ConvToBin("Bin_ChromosomesReplicating", "self.Cel.Len_ChromosomesReplicating", ">=", 0)
            Writer.Multiply_("Len_ChromosomesReplicating_ToEvaluate", "self.Cel.Len_ChromosomesReplicating", "Bin_ChromosomesReplicating")
            Writer.ReduceSum("LengthTotal_NTs", "Len_ChromosomesReplicating_ToEvaluate", 1)
            Writer.Cast_____("LengthTotal_NTs", "LengthTotal_NTs", 'float32')
            Writer.Divide___("LengthTotal_BPs", "LengthTotal_NTs", 2)
            Writer.Divide___("FractionCompletion", "LengthTotal_BPs", Comp.Chromosome.Len_ChromosomesInGenome[0])
            Writer.Multiply_("PercentCompletion", "FractionCompletion", 100)
            Writer.ReturnVar("PercentCompletion")
            Writer.BlankLine()

        with Writer.Function_("ViewProcessSummary"):
            # Writer.PrintStrg("===== DNAReplication Initiation ===== ")
            # # Number of Total Replisome
            # Writer.PrintStVa("# of Total Replisomes",
            #                  "self.Cel.Count_Replisome")
            # # Number of Active Replisome
            # Writer.PrintStVa("# of active Replisomes",
            #                  "self.Cel.Count_ReplisomeActive")
            # # Number of Replisome binding
            # Writer.PrintStVa("# of Replisomes available that can bind a promoter",
            #                  "self.Cel.Count_ReplisomeActiveCanBind")
            # # Number of new Replisome binding in this step (= new Replicating Chromosomes to elongate next step)
            # Writer.PrintStVa("# of Replisomes that newly bind a promoter in this step",
            #                  "self.Cel.Count_ReplisomeWillBind")
            # Writer.BlankLine()

            Writer.PrintStrg("===== DNA Strand Elongation ===== ")
            # # Number of Replicating Chromosomes elongated
            # Writer.PrintStVa("# of Replicating Chromosomes",
            #                  "self.Cel.Count_ChromosomesReplicatingInitialTotal")
            # Total elongation length of Chromosomes
            Writer.Subtract_("DeltaLen_ChromosomesReplicating", "self.Cel.Len_ChromosomesReplicating[:,:2]", "self.Len_ChromosomesReplicatingInitial[:,:2]")
            Writer.ReduceSum("DeltaLen_ChromosomesReplicatingTotal", "DeltaLen_ChromosomesReplicating")
            Writer.PrintStVa("Total Elongation Length of Replicating Chromosomes (bp)",
                             "DeltaLen_ChromosomesReplicatingTotal")
            # # Resulting length of elongating DNA strands of the Chromosomes
            # Writer.PrintStVa("Resulting Length of Elongating DNA Strands of the Chromosomes (nt)",
            #                  "self.Cel.Len_ChromosomesReplicatingFinal[0]")
            # % Replication completion
            Writer.PrintStVa("% Replication completion",
                             "self.DeterminePercentReplicationCompletion()")
            # New gene copies
            Writer.PrintStVa("# of new genes copied",
                             "self.N_GenesReplicated")
            # # New gene copy indices
            # Writer.PrintStVa("Indices of new genes copied",
            #                  "self.Idx_GenesReplicated")

            # Total dNTP consumption and PPi production
            Writer.PrintStVa("Total dNTP consumption [ATP, CTP, GTP, UTP]",
                             "self.Count_dNTPConsumption")
            # Writer.PrintStVa("Total PPi production",
            #                  "self.Count_DNAStrandElongationPPiProduction")
            Writer.BlankLine()

            # Writer.PrintStrg("===== DNA Replication Termination ===== ")
            # Number of Chromosome elongation Completed
            Writer.Equal____("Bool_ChromosomesReplicatingElongationCompleted",
                             "self.Cel.Len_ChromosomesReplicating", "self.Cel.Len_ChromosomesReplicatingMax")
            Writer.ReduceAll("Bool_ChromosomesReplicatingCompletedTotal", "Bool_ChromosomesReplicatingElongationCompleted", 1)
            Writer.NonZeros_("Count_ChromosomesReplicatingCompletedTotal", "Bool_ChromosomesReplicatingCompletedTotal")
            Writer.PrintStVa("# of Chromosome Elongation Completed",
                             "Count_ChromosomesReplicatingCompletedTotal")
            # # Number of Replisomes released
            # Writer.PrintStVa("# of Replisomes released",
            #                  "self.Cel.Count_ReplisomeReleased")
            # # Number of Replisome bound to DNA
            # Writer.PrintStVa("# of total Replisomes bound to DNA",
            #                  "self.Cel.Count_ReplisomeBound")
            # # Number of Replisome freely floating
            # Writer.PrintStVa("# of total Replisomes unbound floating",
            #                  "self.Cel.Count_ReplisomeUnbound")
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



