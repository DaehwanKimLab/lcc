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

            Writer.Variable_("self.Len_ChsReplicating", 0)
            Writer.Variable_("self.Len_RNAsNascent", 0)
            Writer.Variable_("self.Len_ProteinsNascent", 0)

            Writer.Variable_("self.Count_NTPs_RNAsCleaved", 0)
            Writer.Variable_("self.Count_AAS_ProteinsCleaved", 0)

            Writer.Variable_("self.MWs", 0)  # MW matrix for molecules
            Writer.BlankLine()

            Writer.Variable_("self.Coeff_Complexation", 0)
            Writer.Variable_("self.Coeff_Equilibrium", 0)
            Writer.Variable_("self.Coeff_Reaction", 0)

            Writer.Statement("self.LoadStaticCompilerData()")
            Writer.Statement("self.InitializeMatrices()")
            Writer.Statement("self.TransposeFreqNCountMatrices()")
            Writer.BlankLine()

            # Initialize variables for (organization purpose)
            Writer.Statement("self.InitializeVariablesForCellProcesses()")
            Writer.BlankLine()

            Writer.Statement("self.CheckCountsPos()")
            Writer.BlankLine()

        with Writer.Statement("def LoadStaticCompilerData(self):"):
            # Load CompilerData.
            SavePath = os.path.realpath(Comp.SavePath)
            SaveFiles = listdir(SavePath)
            FileTypes = ['npy']

            SavedDataType4Int = ['Coeff', 'Coord', 'Count', 'Dir', 'Idx', 'Len', 'NUniq', 'Rev']
            SavedDataType4Float = ['Const', 'Freq', 'MW']
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

        with Writer.Statement("def InitializeDeltaCounts(self):"):
            Writer.InitZeros("self.DeltaCounts", "self.N_Counts")
            Writer.BlankLine()

        with Writer.Statement("def ClearDeltaCounts(self):"):
            Writer.Statement("self.InitializeDeltaCounts()")
            Writer.BlankLine()

        with Writer.Statement("def ShowDeltaCounts(self):"):
            Writer.NotEqual_("Bool_CheckPointNonZero", "self.DeltaCounts", "tf.constant(0)")
            Writer.PrintStrg("[DeltaCounts] Non-zeros:")
            Writer.PrintStVa("\tIndex", "tf.reshape(tf.where(Bool_CheckPointNonZero), -1)")
            Writer.PrintStVa("\tValue", "tf.reshape(tf.gather(self.DeltaCounts, tf.where(Bool_CheckPointNonZero)), -1)")
            Writer.BlankLine()

        with Writer.Statement("def CheckDeltaCountsNeg(self):"):
            Writer.Less_____("Bool_CheckPointLessThanZero", "self.DeltaCounts", "tf.constant(0)")
            Writer.PrintStrg("[DeltaCounts] Negatives:")
            Writer.PrintStVa("\tIndex", "tf.reshape(tf.where(Bool_CheckPointLessThanZero), -1)")
            Writer.PrintStVa("\tValue", "tf.reshape(tf.gather(self.DeltaCounts, tf.where(Bool_CheckPointLessThanZero)), -1)")
            Writer.BlankLine()

        with Writer.Statement("def CheckCountsPos(self):"):
            if Writer.Switch4SoftCheckCounts:
                Writer.Comment__("Soft Checkpoint")
                Writer.Less_____("Bool_SoftCheckPoint", "self.Counts", "tf.constant(0)")
                Writer.PrintStrg("[Counts] Negatives:")
                Writer.PrintStVa("\tIndex", "tf.reshape(tf.where(Bool_SoftCheckPoint), -1)")
                Writer.PrintStVa("\tValue", "tf.reshape(tf.gather(self.Counts, tf.where(Bool_SoftCheckPoint)), -1)")
                Writer.BlankLine()

            if Writer.Switch4HardCheckCounts:
                Writer.Comment__("Hard Checkpoint")

                Writer.ConvToBin("Bin_NegPos", "self.Counts", "<", "0")
                Writer.GenIdxCnd("Idx_NegPos", "Bin_NegPos", "==", "1")
                Writer.ShapeAxis("N_NegPos", "Idx_NegPos", 0)
                Writer.AsrtNoNeg("self.Counts", "'# of Negative Values in Cel.Counts: %s, Indices: %s' % (N_NegPos, Idx_NegPos)")
                Writer.BlankLine()

            else:
                Writer.Pass_____()
                Writer.BlankLine()

        with Writer.Statement("def InitializeMatrices(self):"):
            Writer.Statement("self.Counts = self.Count_Master")
            Writer.Statement("self.N_Counts = len(self.Counts)")
            Writer.InitZeros("self.DeltaCounts", "self.N_Counts")
            Writer.Statement("self.MWs = self.MW_Master")
            Writer.BlankLine()

            Writer.Statement("self.Coeff_Complexation = self.Coeff_MolsInCPLXRXN")
            Writer.Statement("self.Coeff_Equilibrium = self.Coeff_MolsInEQMRXN")
            Writer.Statement("self.Coeff_Metabolism = self.Coeff_MolsInMETRXN")
            Writer.BlankLine()

            # Temporary code
            # E coli cell volume: 0.7 um3 (average), which is 7e-16 liters
            Writer.Variable_("self.Vol", 7e-16)  # TO BE REPLACED AND MOVED INTO SIMULATION
            Writer.BlankLine()

            Writer.Statement("self.SwitchToSparseMatrices()")
            Writer.BlankLine()

        with Writer.Statement("def SwitchToSparseMatrices(self):"):
            ListOfMatrixVariablesToSwitchToSparseTensorType = [
                "self.Coeff_Complexation",
                # "self.Coeff_Equilibrium",
                # "self.Coeff_Reaction",
            ]
            for Variable in ListOfMatrixVariablesToSwitchToSparseTensorType:
                Writer.ToSparse_(Variable, Variable)
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def TransposeFreqNCountMatrices(self):"):
            ListOfMatrixVariablesToTranspose = [
                "self.Freq_NTsInChromosomes",
                "self.Freq_NTsInChromosomesReplicating",
                "self.Freq_NTsInChromosomesInGenome",
                "self.Freq_NTsInRNAs",
                "self.Freq_AAsInProteins",
                "self.Count_NTsInChromosomesInGenome",
                "self.Count_NTsInRNAs",
                "self.Count_AAsInProteins",
            ]
            for Variable in ListOfMatrixVariablesToTranspose:
                Writer.Transpose(Variable, Variable)
            Writer.BlankLine()

        with Writer.Statement("def InitializeVariablesForCellProcesses(self):"):

            Writer.Comment__("Popular")
            Writer.BlankLine()
            Writer.Variable_("self.One", 1)
            Writer.Variable_("self.Zero", 0)
            Writer.BlankLine()

            Idx_H2O = Comp.Master.ID2Idx_Master['WATER[c]']
            Idx_Proton = Comp.Master.ID2Idx_Master['PROTON[c]']
            Idx_ATP = Comp.Master.ID2Idx_Master['ATP[c]']
            Idx_ADP = Comp.Master.ID2Idx_Master['ADP[c]']
            Idx_Pi = Comp.Master.ID2Idx_Master['PI[c]']
            Idx_PPi = Comp.Master.ID2Idx_Master['PPI[c]']

            Writer.Variable_("self.Idx_H2O", Idx_H2O)
            Writer.Variable_("self.Idx_Proton", Idx_Proton)
            Writer.Variable_("self.Idx_ATP", Idx_ATP)
            Writer.Variable_("self.Idx_ADP", Idx_ADP)
            Writer.Variable_("self.Idx_Pi", Idx_Pi)
            Writer.Variable_("self.Idx_PPi", Idx_PPi)
            Writer.BlankLine()

            Idx_dNTPs = ProGen.BuildingBlockIdxs('dNTPs')
            Idx_NTPs = ProGen.BuildingBlockIdxs('NTPs')
            Idx_AAs = ProGen.BuildingBlockIdxs('AAs')

            Writer.Variable_("self.Idx_dNTPs", Idx_dNTPs)
            Writer.Variable_("self.Idx_NTPs", Idx_NTPs)
            Writer.Variable_("self.Idx_AAs", Idx_AAs)
            Writer.BlankLine()

            Writer.Comment__("Replication")
            Writer.BlankLine()

            Writer.Variable_("self.Idx_DnaA", 0)
            Writer.Variable_("self.Idx_DnaB", 0)
            Writer.Variable_("self.Idx_DnaC", 0)

            Writer.Variable_("self.Idx_Ch_Original", 0)
            Writer.Variable_("self.Idx_Ch_Replicating", 0)


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
            Writer.BlankLine()
            Writer.Variable_("self.Len_RNAsNascent", 0)
            Writer.BlankLine()

            # Translation
            Writer.Comment__("Translation")
            Writer.BlankLine()
            Writer.Variable_("self.Idx_Ribosome30S", 0)
            Writer.Variable_("self.Idx_Ribosome50S", 0)
            Writer.Variable_("self.Idx_Ribosome70S", 0)
            Writer.Variable_("self.Idx_SelenoCysteineInAAs", 0)
            Writer.Variable_("self.Idx_AAsLocalAssignmentNoSelenoCysteine", 0)
            Writer.BlankLine()

            Writer.Variable_("self.Len_ProteinsNascent", 0)
            Writer.BlankLine()

            # RNA Degradation
            Writer.Comment__("RNA Degradation")
            Writer.BlankLine()
            Writer.Variable_("self.Count_NTsInRNAsCleaved", 0)
            Writer.BlankLine()

            # Protein Degradation
            Writer.Comment__("Protein Degradation")
            Writer.BlankLine()
            Writer.BlankLine()

            # Complexation
            Writer.Comment__("Complexation")
            Writer.BlankLine()
            Writer.Variable_("self.Idx_MolsInCPLXRXN", 0)
            Writer.BlankLine()

            # Equilibrium
            Writer.Comment__("Equilibrium")
            Writer.BlankLine()
            Writer.Variable_("self.Idx_MolsInEQM", 0)
            Writer.BlankLine()

            # Metabolism
            Writer.Comment__("Metabolism")
            Writer.BlankLine()

            Idx_NADH = Comp.Master.ID2Idx_Master['NADH[c]']
            Idx_NADPH = Comp.Master.ID2Idx_Master['NADPH[c]']
            Idx_FADH2 = Comp.Master.ID2Idx_Master['FADH2[p]']

            Writer.Variable_("self.Idx_NADH", Idx_NADH)
            Writer.Variable_("self.Idx_NADPH", Idx_NADPH)
            Writer.Variable_("self.Idx_FADH2", Idx_FADH2)
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
        #     Writer.Concat___("self.MX_Stoichiometries", "self.MX_Stoichiometries", "Stoichiometry")
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def AddToRateMatrix(self, Rate):"):
        #     Writer.Concat___("self.MX_Rates", "self.MX_Rates", "Rate")
        #     Writer.BlankLine()

        with Writer.Statement("def GetCounts(self, MolIdxs):"):
            Writer.Gather___("Counts", "self.Counts", "MolIdxs")
            Writer.Reshape__("Counts", "Counts", -1)
            Writer.ReturnVar("Counts")
            Writer.BlankLine()

        with Writer.Statement("def GetDeltaCounts(self, MolIdxs):"):
            Writer.Gather___("DeltaCounts", "self.DeltaCounts", "MolIdxs")
            Writer.Reshape__("DeltaCounts", "DeltaCounts", -1)
            Writer.ReturnVar("DeltaCounts")
            Writer.BlankLine()

        with Writer.Statement("def AddToDeltaCounts(self, MolIdxs, MolCounts):"):
            Writer.ScatNdAdd("self.DeltaCounts", "self.DeltaCounts", "MolIdxs", "MolCounts")
            Writer.BlankLine()