import numpy as np
import os, sys
from os import listdir


# Abstract Methods for CellProcess
def Write_CellState(Writer, Comp, ProGen):
    Writer.BlankLine()
    with Writer.Statement("class FCellState():"):
        with Writer.Statement("def __init__(self):"):

            Writer.Variable_("self.Species", 0) # Index for Exact or closest species
            Writer.Variable_("self.ID", 0) # Not implemented yet
            Writer.Variable_("self.Vol", 0) # Not implemented yet
            Writer.BlankLine()

            # Writer.Variable_("self.MX_Stoichiometries", 0)  # Stoichiometry matrix for all elementary reactions
            # Writer.Variable_("self.MX_Rates", 0)  # Rate matrix for all elementary reactions

            Writer.Variable_("self.Counts", 0)  # Counts matrix for all molecules
            Writer.Variable_("self.DeltaCounts", 0)
            Writer.Variable_("self.MWs", 0)  # MW matrix for molecules
            Writer.BlankLine()

            Writer.Statement("self.ImportCompilerData()")
            Writer.Statement("self.TransposeFreqMatrices()")
            Writer.BlankLine()

            # Initialize variables for (organization purpose)
            Writer.Statement("self.InitializeVariablesForCellProcesses()")
            Writer.BlankLine()

            Writer.Statement("self.CheckCountsPos()")
            Writer.BlankLine()

        with Writer.Statement("def ImportCompilerData(self):"):
            # Load CompilerData.
            SavePath = os.path.realpath(Comp.SavePath)
            SaveFiles = listdir(SavePath)
            FileTypes = ['npy']

            SavedDataType4Int = ['Coord', 'Count', 'Dir', 'Idx', 'Len', 'NUniq', 'Rev']
            SavedDataType4Float = ['Freq', 'MW']
            SavedDataTypeDict = dict()
            for SavedDataType in SavedDataType4Int:
                SavedDataTypeDict[SavedDataType] = 'int32'
            for SavedDataType in SavedDataType4Float:
                SavedDataTypeDict[SavedDataType] = 'float32'
            for SaveFile in SaveFiles:
                # Load only if File type is .npy and only contains numerical values
                SavedFileType = SaveFile.split('.')[-1]
                SavedDataType = SaveFile.split('_')[0]
                if (SavedFileType in FileTypes) & (SavedDataType in SavedDataTypeDict):
                    DataType = SavedDataTypeDict[SavedDataType]
                    SaveFilePath = os.path.join(SavePath, SaveFile)
                    VariableName = SaveFile.split('.')[0]
                    Writer.LoadSaved(SaveFilePath, VariableName, DataType)
            Writer.BlankLine()

            Writer.Statement("self.Counts = self.Count_Master")
            Writer.Statement("self.N_Counts = len(self.Counts)")
            Writer.InitZeros("self.DeltaCounts", "self.N_Counts")
            Writer.Statement("self.MWs = self.MW_Master")

            # Temporary code
            # E coli cell volume: 0.7 um3 (average), which is 7e-16 liters
            Writer.Variable_("self.Vol", 7e-16)  # TO BE REPLACED AND MOVED INTO SIMULATION
            Writer.BlankLine()

        with Writer.Statement("def InitializeDeltaCounts(self):"):
            Writer.InitZeros("self.DeltaCounts", "self.N_Counts")
            Writer.BlankLine()

        with Writer.Statement("def ClearDeltaCounts(self):"):
            Writer.Statement("self.InitializeDeltaCounts()")
            Writer.BlankLine()

        with Writer.Statement("def ShowDeltaCounts(self):"):
            Writer.OperElNEq("Bool_CheckPointNonZero", "self.DeltaCounts", "tf.constant(0)")
            Writer.PrintStrg("[DeltaCounts] Non-zeros:")
            Writer.PrintStVa("\tIndex", "tf.reshape(tf.where(Bool_CheckPointNonZero), -1)")
            Writer.PrintStVa("\tValue", "tf.reshape(tf.gather(self.DeltaCounts, tf.where(Bool_CheckPointNonZero)), -1)")
            Writer.BlankLine()

        with Writer.Statement("def CheckDeltaCountsNeg(self):"):
            Writer.OperElLe_("Bool_CheckPointLessThanZero", "self.DeltaCounts", "tf.constant(0)")
            Writer.PrintStrg("[DeltaCounts] Negatives:")
            Writer.PrintStVa("\tIndex", "tf.reshape(tf.where(Bool_CheckPointLessThanZero), -1)")
            Writer.PrintStVa("\tValue", "tf.reshape(tf.gather(self.DeltaCounts, tf.where(Bool_CheckPointLessThanZero)), -1)")
            Writer.BlankLine()

        with Writer.Statement("def CheckCountsPos(self):"):
            if Writer.Switch4SoftCheckCounts:
                Writer.Comment__("Soft Checkpoint")
                Writer.OperElLe_("Bool_SoftCheckPoint", "self.Counts", "tf.constant(0)")
                Writer.PrintStrg("[Counts] Negatives:")
                Writer.PrintStVa("\tIndex", "tf.reshape(tf.where(Bool_SoftCheckPoint), -1)")
                Writer.PrintStVa("\tValue", "tf.reshape(tf.gather(self.Counts, tf.where(Bool_SoftCheckPoint)), -1)")
                Writer.BlankLine()

            if Writer.Switch4HardCheckCounts:
                Writer.Comment__("Hard Checkpoint")
                Writer.AsrtNoNeg("self.Counts")
                Writer.BlankLine()

            else:
                Writer.Pass_____()
                Writer.BlankLine()

        with Writer.Statement("def TransposeFreqMatrices(self):"):
            Writer.Transpose("self.Freq_NTsInChromosomes")
            Writer.Transpose("self.Freq_NTsInChromosomesReplicating")
            Writer.Transpose("self.Freq_NTsInChromosomesInGenome")
            Writer.Transpose("self.Freq_NTsInRNAs")
            Writer.Transpose("self.Freq_AAsInProteins")
            Writer.BlankLine()

        with Writer.Statement("def InitializeMatrices(self):"):
            Writer.Reshape__("self.Counts", "self.Count_Master", [-1, 1])
            Writer.Reshape__("self.MWs", "self.MW_Master", [-1, 1])
            Writer.BlankLine()

        with Writer.Statement("def InitializeVariablesForCellProcesses(self):"):
            # Popular
            Writer.Comment__("Popular")
            Writer.Variable_("self.One", 1)
            Writer.Variable_("self.Zero", 0)

            Writer.Variable_("self.Idx_PPi", 0)

            Writer.BlankLine()

            # Replication
            Writer.Comment__("Replication")
            Writer.Variable_("self.Idx_Ch_Original", 0)
            Writer.Variable_("self.Idx_Ch_Replicating", 0)
            Writer.Variable_("self.Idx_dNTPs", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Rate_DNAReplication", 0)
            Writer.Variable_("self.Rate_DNAReplication_Matrix", 0)
            Writer.BlankLine()
            # Writer.Variable_("self.Count_Replisome", 0)
            # Writer.Variable_("self.Count_ReplisomeActive", 0)
            # Writer.Variable_("self.Count_ReplisomeActiveCanBind", 0)
            # Writer.Variable_("self.Count_ReplisomeActiveToBindThisStep", 0)
            # Writer.Variable_("self.Count_ReplisomeBound", 0)
            # Writer.Variable_("self.Count_ReplisomeUnbound", 0)
            # Writer.Variable_("self.Count_ReplisomeReleased", 0)
            Writer.Variable_("self.Count_ChromosomesReplicating_Matrix", 0)
            Writer.Variable_("self.Count_ChromosomesReplicatingInitialTotal", 0)
            Writer.Variable_("self.Count_ChromosomesReplicatingElongatingTotal", 0)
            Writer.Variable_("self.Count_ChromosomesReplicatingOverElongatedTotal", 0)
            Writer.Variable_("self.Count_ChromosomesReplicatingFinalTotal", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_DNAStrandElongationNTLengthFinalPerChromosome", 0)
            Writer.Variable_("self.Count_DNAStrandElongationNTLengthFinalTotal", 0)
            Writer.Variable_("self.Count_DNAStrandElongationdNTPConsumption", 0)
            Writer.Variable_("self.Count_DNAStrandElongationPPiProduction", 0)
            Writer.Variable_("self.Count_DNAStrandElongationCompletedPerChromosome", 0)
            Writer.Variable_("self.Count_DNAStrandElongationCompletedTotal", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_ChromosomeElongationLengthBPPerChromosome", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_dNTPsOverElongated", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Len_ChromosomesOriginal", 0)
            Writer.Variable_("self.Len_ChromosomesReplicatingInitial", 0)
            Writer.Variable_("self.Len_ChromosomesReplicatingMax", 0)
            Writer.Variable_("self.Len_ChromosomesReplicatingElongated", 0)
            Writer.Variable_("self.Len_ChromosomesReplicatingOverElongated", 0)
            Writer.Variable_("self.Len_ChromosomesReplicatingAdjusted", 0)
            Writer.Variable_("self.Len_ChromosomesReplicatingCompleted", 0)
            Writer.Variable_("self.Len_ChromosomesReplicatingFinal", 0)

            # Transcription
            Writer.Comment__("Transcription")
            Writer.BlankLine()
            Writer.Variable_("self.Idx_RNAP", 0)
            Writer.Variable_("self.Idx_NTPs", 0)
            Writer.Variable_("self.Idx_RndRNAsNascent", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Rate_RNAPActive", 0)  # Not implemented yet
            Writer.Variable_("self.Rate_RNAPActiveCanBind", 0)
            Writer.Variable_("self.Rate_RNAElongation", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_RNAP", 0)
            Writer.Variable_("self.Count_RNAPActive", 0)
            Writer.Variable_("self.Count_RNAPActiveCanBind", 0)
            Writer.Variable_("self.Count_RNAPActiveToBindThisStep", 0)
            Writer.Variable_("self.Count_RNAPBound", 0)
            Writer.Variable_("self.Count_RNAPUnbound", 0)
            Writer.Variable_("self.Count_RNAPReleased", 0)
            Writer.Variable_("self.Count_RNAsNascent_Matrix", 0)
            Writer.Variable_("self.Count_RNAsNascentInitialTotal", 0)
            Writer.Variable_("self.Count_RNAsNascentElongatingTotal", 0)
            Writer.Variable_("self.Count_RNAsNascentOverElongatedTotal", 0)
            Writer.Variable_("self.Count_RNAsNascentFinalTotal", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_RNAElongationLengthFinalPerRNA", 0)
            Writer.Variable_("self.Count_RNAElongationLengthFinalTotal", 0)
            Writer.Variable_("self.Count_RNAElongationNTPConsumption", 0)
            Writer.Variable_("self.Count_RNAElongationPPiProduction", 0)
            Writer.Variable_("self.Count_RNAElongationCompletedPerRNA", 0)
            Writer.Variable_("self.Count_RNAElongationCompletedTotal", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_NTPsOverElongated", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Len_RNAsNascentInitial", 0)
            Writer.Variable_("self.Len_RNAsNascentMax", 0)
            Writer.Variable_("self.Len_RNAsNascentElongated", 0)
            Writer.Variable_("self.Len_RNAsNascentOverElongated", 0)
            Writer.Variable_("self.Len_RNAsNascentAdjusted", 0)
            Writer.Variable_("self.Len_RNAsNascentCompleted", 0)
            Writer.Variable_("self.Len_RNAsNascentFinal", 0)
            Writer.BlankLine()

            # Translation
            Writer.Comment__("Translation")
            Writer.BlankLine()
            Writer.Variable_("self.Idx_Ribosome30S", 0)
            Writer.Variable_("self.Idx_Ribosome50S", 0)
            Writer.Variable_("self.Idx_Ribosome70S", 0)
            Writer.Variable_("self.Idx_AAs", 0)
            Writer.Variable_("self.Idx_SelenoCysteineInAAs", 0)
            Writer.Variable_("self.Idx_AAsLocalAssignmentNoSelenoCysteine", 0)
            Writer.Variable_("self.Idx_RndProteinsNascent", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Rate_RibosomeActive", 0)  # Not implemented yet
            Writer.Variable_("self.Rate_RibosomeActiveCanBind", 0)
            Writer.Variable_("self.Rate_ProteinElongation", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_mRNACopies", 0)
            Writer.Variable_("self.Count_Ribosome30S", 0)
            Writer.Variable_("self.Count_Ribosome50S", 0)
            Writer.Variable_("self.Count_Ribosome70S", 0)
            Writer.Variable_("self.Count_Ribosome", 0)
            Writer.Variable_("self.Count_RibosomeActive", 0)
            Writer.Variable_("self.Count_RibosomeActiveCanBind", 0)
            Writer.Variable_("self.Count_RibosomeActiveToBindThisStep", 0)
            Writer.Variable_("self.Count_RibosomeBound", 0)
            Writer.Variable_("self.Count_RibosomeUnbound", 0)
            Writer.Variable_("self.Count_RibosomeReleased", 0)
            Writer.Variable_("self.Count_ProteinsNascent_Matrix", 0)
            Writer.Variable_("self.Count_ProteinsNascentInitialTotal", 0)
            Writer.Variable_("self.Count_ProteinsNascentElongatingTotal", 0)
            Writer.Variable_("self.Count_ProteinsNascentOverElongatedTotal", 0)
            Writer.Variable_("self.Count_ProteinsNascentFinalTotal", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_ProteinElongationLengthFinalPerProtein", 0)
            Writer.Variable_("self.Count_ProteinElongationLengthFinalTotal", 0)
            Writer.Variable_("self.Count_ProteinElongationAAConsumption", 0)
            Writer.Variable_("self.Count_ProteinElongationPPiProduction", 0)
            Writer.Variable_("self.Count_ProteinElongationCompletedPerProtein", 0)
            Writer.Variable_("self.Count_ProteinElongationCompletedTotal", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_AAsOverElongated", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Len_ProteinsNascentInitial", 0)
            Writer.Variable_("self.Len_ProteinsNascentMax", 0)
            Writer.Variable_("self.Len_ProteinsNascentElongated", 0)
            Writer.Variable_("self.Len_ProteinsNascentOverElongated", 0)
            Writer.Variable_("self.Len_ProteinsNascentAdjusted", 0)
            Writer.Variable_("self.Len_ProteinsNascentCompleted", 0)
            Writer.Variable_("self.Len_ProteinsNascentFinal", 0)
            Writer.BlankLine()


            # Metabolism
            Writer.Comment__("Metabolism")
            Writer.BlankLine()
            Writer.Variable_("self.Count_MetabolitesInitial", 0)


        #     Writer.Statement("self.InitializeStoichiometryMatrix()")
        #     Writer.Statement("self.InitializeRateMatrix()")
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def FinalizeReactionMatrices(self):"):
        #     Writer.Statement("self.FinalizeStoichiometryMatrix()")
        #     Writer.Statement("self.FinalizeRateMatrix()")
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def FinalizeStoichiometryMatrix(self):"):
        #     # Transposition of RXN stoichiometry matrix is necessary for proper matrix multiplication
        #     Writer.Transpose("self.MX_Stoichiometries")
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def FinalizeRateMatrix(self):"):
        #     # Reshaping of RXN rate matrix to be 2-dimentional is necessary for proper matrix multiplication
        #     Writer.Reshape__("self.MX_Rates", "self.MX_Rates", [-1, 1])
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def InitializeStoichiometryMatrix(self):"):
        #     # Reset self.Stoichs to a 2D array of zeros for concatenation
        #     Writer.InitZeros("self.MX_Stoichiometries", [1, Comp.Master.NUniq_Master])
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def InitializeRateMatrix(self):"):
        #     # Reset self.Rates to a 2D array of zero for concatenation
        #     Writer.InitZeros("self.MX_Rates", [1, 1])  # Rate matrix for all elementary reactions
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def ClearRateMatrix(self):"):
        #     # Reset self.Rates to the initialized state.
        #     Writer.Statement("self.InitializeRateMatrix()")
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def AddToStoichiometryMatrix(self, Stoichiometry):"):
        #     Writer.OperCncat("self.MX_Stoichiometries", "self.MX_Stoichiometries", "Stoichiometry")
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def AddToRateMatrix(self, Rate):"):
        #     Writer.OperCncat("self.MX_Rates", "self.MX_Rates", "Rate")
        #     Writer.BlankLine()