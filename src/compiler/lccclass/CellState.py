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
            Writer.Variable_("self.One", 1)
            Writer.BlankLine()

            # Replication
            Writer.BlankLine()

            # Transcription
            Writer.Variable_("self.Idx_RNAP", 0)
            Writer.Variable_("self.Idx_NTPs", 0)
            Writer.Variable_("self.Idx_PPi", 0)
            Writer.Variable_("self.Idx_RndRNAsNascent", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Rate_RNAPActive", 0)  # Not implemented yet
            Writer.Variable_("self.Rate_RNAPActiveCanBind", 0)
            Writer.Variable_("self.Rate_RNAElongation", 0)
            Writer.Variable_("self.Rate_RNAElongation_Scalar", 0)
            Writer.Variable_("self.Rate_RNAElongation_Vector", 0)
            Writer.Variable_("self.Rate_RNAElongation_Vector_Corrected", 0)
            Writer.Variable_("self.Rate_RNAElongation_Matrix", 0)
            Writer.Variable_("self.Rate_RNAElongation_Matrix_Corrected", 0)
            Writer.Variable_("self.Rate_RNAOverElongation", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_GeneCopies", 0)
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
            Writer.Variable_("self.Bool_RNAsNascentElongating", 0)
            Writer.Variable_("self.Bool_RNAsNascentOverElongated", 0)
            Writer.Variable_("self.Bool_RNAsNascentElongationCompleted", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Bin_RNAsNascentElongating", 0)
            Writer.Variable_("self.Bin_RNAsNascentOverElongated", 0)
            Writer.Variable_("self.Bin_RNAsNascentElongationCompleted", 0)
            Writer.BlankLine()

            # Translation
            Writer.Variable_("self.Idx_Ribosome50S", 0)
            Writer.Variable_("self.Idx_Ribosome30S", 0)
            Writer.Variable_("self.Idx_AAs", 0)
            Writer.Variable_("self.Idx_H2O", 0)
            Writer.Variable_("self.Idx_RndProteinsNascent", 0)
            Writer.Variable_("self.Idx_ProteinElongationCompleted", 0)
            Writer.Variable_("self.Rate_RibosomeActive", 0)  # Not implemented yet
            Writer.Variable_("self.Rate_RibosomeFreeBinding", 0)
            Writer.Variable_("self.Rate_ProteinElongation", 0)
            Writer.Variable_("self.Rate_ProteinElongation_Scalar", 0)
            Writer.Variable_("self.Rate_ProteinElongation_Vector", 0)
            Writer.Variable_("self.Rate_ProteinElongation_Vector_Corrected", 0)
            Writer.Variable_("self.Rate_ProteinElongation_Matrix", 0)
            Writer.Variable_("self.Rate_ProteinElongation_Matrix_Corrected", 0)
            Writer.Variable_("self.Rate_ProteinOverElongation", 0)
            Writer.Variable_("self.Count_GeneCopies", 0)
            Writer.Variable_("self.Count_RibosomeFree", 0)
            Writer.Variable_("self.Count_RibosomeFreeBinding", 0)
            Writer.Variable_("self.Count_RibosomeReleased", 0)
            Writer.Variable_("self.Count_ProteinsNascent", 0)
            Writer.Variable_("self.Count_ProteinsNascentElongating", 0)
            Writer.Variable_("self.Count_ProteinsNascentOverElongated", 0)
            Writer.Variable_("self.Count_ProteinElongationLengthPerRNA", 0)
            Writer.Variable_("self.Count_ProteinElongationLengthTotal", 0)
            Writer.Variable_("self.Count_ProteinElongationAAConsumption", 0)
            Writer.Variable_("self.Count_ProteinElongationH2OProduction", 0)
            Writer.Variable_("self.Count_ProteinElongationCompletedPerProtein", 0)
            Writer.Variable_("self.Count_ProteinElongationCompletedTotal", 0)
            Writer.Variable_("self.Len_ProteinsNascentElongated", 0)
            Writer.Variable_("self.Len_ProteinsNascentOverElongated", 0)
            Writer.Variable_("self.Len_ProteinsNascent", 0)
            Writer.Variable_("self.Len_ProteinsNascentMax", 0)
            Writer.Variable_("self.Bool_ProteinsNascentElongating", 0)
            Writer.Variable_("self.Bool_ProteinsNascentOverElongated", 0)
            Writer.Variable_("self.Bool_ProteinsNascentElongationCompleted", 0)
            Writer.Variable_("self.Bin_ProteinsNascentElongating", 0)
            Writer.Variable_("self.Bin_ProteinsNascentOverElongated", 0)
            Writer.Variable_("self.Bin_ProteinsNascentElongationCompleted", 0)
            Writer.BlankLine()
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