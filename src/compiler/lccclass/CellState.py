import numpy as np
import os, sys
from os import listdir


# Abstract Methods for CellProcess
def Write_CellState(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FCellState():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.Species", 0) # Index for Exact or closest species
            Writer.Variable_("self.ID", 0) # Not implemented yet
            Writer.Variable_("self.Vol", 0) # Not implemented yet


            Writer.Variable_("self.Stoichs", 0)  # Stoichiometry matrix for all elementary reactions
            Writer.Variable_("self.Rates", 0)  # Rate matrix for all elementary reactions
            Writer.Variable_("self.Counts", 0)  # Counts matrix for all molecules
            Writer.Variable_("self.MWs", 0)  # MW matrix for molecules
            Writer.Variable_("self.DeltaCounts", 0)

            Writer.Statement("self.Initialize()")
            Writer.BlankLine()

        with Writer.Statement("def Initialize(self):"):
            # Load CompilerData.
            SavePath = os.path.realpath(Comp.SavePath)
            SaveFiles = listdir(SavePath)
            FileTypes = ['npy']

            SavedDataType4Int = ['Coord', 'Count', 'Dir', 'Idx', 'Len', 'NUniq', 'Rev']
            SavedDataType4Float = ['Freq', 'MW']
            SavedDataTypeDict = {}
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

        with Writer.Statement("def SetUpMatrices(self):"):
            Writer.InitZeros("self.Stoichs", [1, Comp.Master.NUniq_Master])
            Writer.InitZeros("self.Rates", [1, 1])  # Rate matrix for all elementary reactions
            Writer.Overwrite("self.Counts", "self.Count_Master")
            Writer.Reshape__("self.Counts", "self.Counts", [-1, 1])
            Writer.Overwrite("self.MWs", "self.MW_Master")
            Writer.Reshape__("self.MWs", "self.MWs", [-1, 1])
            Writer.BlankLine()

        with Writer.Statement("def FinalizeMatrices(self):"):
            # Transposition of RXN stoichiometry matrix is necessary for proper matrix multiplication
            Writer.Transpose("self.Stoichs")
            # Reshaping of RXN rate matrix to be 2-dimentional is necessary for proper matrix multiplication
            Writer.Statement("self.FinalizeRateMatrix()")
            Writer.BlankLine()

        with Writer.Statement("def FinalizeRateMatrix(self):"):
            # Reshaping of RXN rate matrix to be 2-dimentional is necessary for proper matrix multiplication
            Writer.Reshape__("self.Rates", "self.Rates", [-1, 1])
            Writer.BlankLine()


