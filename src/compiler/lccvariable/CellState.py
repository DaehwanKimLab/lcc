import numpy as np
import os, sys
from os import listdir

def Write_CellState_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("class FCellState():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.Species", 0) # Index for Exact or closest species
            Writer.Variable_("self.CellID", 0) # Not implemented yet
            Writer.Variable_("self.CellVol", 0) # Not implemented yet
            Writer.Statement("self.Initialize()")
            Writer.BlankLine()

        with Writer.Statement("def Initialize(self):"):
            # Load CompilerData.
            SavePath = os.path.realpath(CompilerData.SavePath)
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
                    Writer.LoadSaved(SavePath, SaveFile, DataType)
