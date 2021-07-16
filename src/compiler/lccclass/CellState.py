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
            Writer.InitZeros("self.DeltaCounts", "len(self.Counts)")
            Writer.Statement("self.MWs = self.MW_Master")

            # Temporary code
            # E coli cell volume: 0.7 um3 (average), which is 7e-16 liters
            Writer.Variable_("self.Vol", 7e-16)  # TO BE REPLACED AND MOVED INTO SIMULATION
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
            Writer.Variable_("self.Idx_NTPs", 0)
            Writer.Variable_("self.Idx_PPi", 0)
            Writer.Variable_("self.Idx_RndRNAsNascent", 0)
            Writer.Variable_("self.Rate_RNAElongation", 0)
            Writer.Variable_("self.Rate_RNAElongation_Scalar", 0)
            Writer.Variable_("self.Rate_RNAElongation_Vector", 0)
            Writer.Variable_("self.Rate_RNAElongation_Vector_Corrected", 0)
            Writer.Variable_("self.Rate_RNAElongation_Matrix", 0)
            Writer.Variable_("self.Rate_RNAElongation_Matrix_Corrected", 0)
            Writer.Variable_("self.Rate_RNAOverElongation", 0)
            Writer.Variable_("self.Count_GeneCopies", 0)
            Writer.Variable_("self.Count_RNAPFree", 0)
            Writer.Variable_("self.Count_RNAsNascent", 0)
            Writer.Variable_("self.Count_RNAsNascentElongating", 0)
            Writer.Variable_("self.Count_RNAsNascentOverElongated", 0)
            Writer.Variable_("self.Count_RNAElongationLengthPerGene", 0)
            Writer.Variable_("self.Count_RNAElongationLengthTotal", 0)
            Writer.Variable_("self.Count_ElongationCompletedPerRNA", 0)
            Writer.Variable_("self.Len_RNAsNascentElongated", 0)
            Writer.Variable_("self.Len_RNAsNascentOverElongated", 0)
            Writer.Variable_("self.Len_RNAsNascent", 0)
            Writer.Variable_("self.Len_RNAsNascentMax", 0)
            Writer.Variable_("self.Bool_RNAsNascentElongating", 0)
            Writer.Variable_("self.Bool_RNAsNascentOverElongated", 0)
            Writer.Variable_("self.Bin_RNAsNascentElongating", 0)
            Writer.Variable_("self.Bin_RNAsNascentOverElongated", 0)
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