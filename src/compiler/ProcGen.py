#!/usr/bin/env python3
#
# Copyright 2021,
# Donghoon Lee <dhldjl@gmail.com>,
# Chanhee Park <parkchanhee@gmail.com>, and
# Daehwan Kim <infphilo@gmail.com>,
#
# This file is part of LDP.
#
# LDP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

import numpy as np
import os
import abc
import CodeGen
from lccclass.cellprocess.synthesis import Replication
from lccclass.cellprocess.synthesis import Transcription
from lccclass.cellprocess.synthesis import Translation

class FProcessGenerator():
    def __init__(self):
        self.Dict_CellProcesses = {}
        self.Dict_ProcessReactions = {}

        self.Reactions_Replication = []
        self.Reactions_Transcription = []
        self.Reactions_Translation = []

    def LinkCompilerObj(self, Comp):
        self.Comp = Comp
        self.ImportSaveVariables()

    def ImportSaveVariables(self):
        self.Switch4SaveAllData = self.Comp.Switch4SaveAllData
        self.SavePath = self.Comp.SavePath

    # def PrepareReactionForMatrix(self, Reaction):
    #     # Put all reaction info in the lists
    #     for Type, RXNID, MoleculeIDs, Stoichiometry, Reversibility, Order, Rate, ModulatorIDs, Trigger in Reaction:
    #         MoleculeIDs, Stoichiometry = Utils.RefineBuildingBlocks(MoleculeIDs, Stoichiometry, self.Comp)
    #         return MoleculeIDs, Stoichiometry
    #

    def SetUpProcesses(self):

        self.Dict_CellProcesses['Replication'] = Replication
        self.Reactions_Replication = Replication.SetUpReactions(self)
        self.Dict_ProcessReactions['Replication'] = self.Reactions_Replication

        # self.Reactions_Transcription = Transcription.SetUpReactions(self)
        # self.Dict_Reactions['Transcription'] = self.Reactions_Transcription

        # self.Reactions_Translation = Translation.SetUpReactions(self)
        # self.Dict_Reactions['Translation'] = self.Reactions_Translation

    def SaveProcesses(self):
        # Save all self.Reactions
        if self.Switch4SaveAllData:
            for ProcessName, Reactions in self.Dict_ProcessReactions.items():
                SaveFileName = "%s/%s" % (self.SavePath, ProcessName)
                np.save('%s.npy' % SaveFileName, Reactions)
                print('Process data have been successfully saved. Process: "%s"' % ProcessName)

    def SetUpReaction(self, MolIDs, Coeffs, Rate):
        # Information from L++ code parsed by the compiler:

        MolIDs_Parsed, Coeffs_Parsed = self.RefineBuildingBlocks(MolIDs, Coeffs)
        MolIdxs = self.GetMolIdx(MolIDs_Parsed, self.Comp.Master.ID2Idx_Master)
        Stoichiometry = [MolIDs_Parsed, MolIdxs, Coeffs_Parsed]

        # Rate.append(Rate)

        return [Stoichiometry, Rate]

    # def SetUpReactions(self, Reactions):
    #     for Reaction in Reactions:
    #         Reaction = self.SetUpReaction(Reaction)

    # Reaction Matrix Building functions
    def ParseBuildingBlocks(self, MolGroup, Stoichiometry):
        MolGroupParsed = []
        StoichiometryParsed = []
        if MolGroup == 'dNTP':
            MolGroupParsed = self.Comp.BuildingBlock.Name_dNTPs
            StoichiometryParsed = self.FlatList(
                self.Comp.Chromosome.Freq_NTsInChromosomesInGenome)  # This is a temporary solution for a single chromosome
        elif MolGroup == 'NTP':
            MolGroupParsed = self.Comp.BuildingBlock.Name_NTPs
            StoichiometryParsed = self.Comp.RNA.Freq_NTsInRNAs
        elif MolGroup == 'AA':
            MolGroupParsed = self.Comp.BuildingBlock.Name_AAs
            StoichiometryParsed = self.Comp.Protein.Freq_AAsInProteins
        StoichiometryParsed = [i * Stoichiometry for i in StoichiometryParsed]
        return MolGroupParsed, StoichiometryParsed

    def RefineBuildingBlocks(self, MolIDs, Coeffs, ReactionOrder=None):
        assert len(MolIDs) == len(Coeffs), "#'s of Molecule IDs and Stoichiometries do not match"
        BuildingBlocks = ['dNTP', 'NTP', 'AA']
        # MolTypes = ['Chromosome', 'Gene', 'Promoter', 'RNA', 'Protein', 'Complex', 'Metabolite']
        MolIDs_Refined = []
        Coeffs_Refined = []
        ReactionOrder_Refined = []
        for MolID, Coeff in zip(MolIDs, Coeffs):
            if MolID in self.Comp.Master.ID_Master:
                MolIDs_Refined.append(MolID)
                Coeffs_Refined.append(Coeff)
            elif MolID in BuildingBlocks:
                MolID_BuildingBlocks, Coeff_BuildingBlocks = self.ParseBuildingBlocks(MolID, Coeff)
                for MolID_BuildingBlock in MolID_BuildingBlocks:
                    assert MolID_BuildingBlock in self.Comp.Master.ID_Master, '%s not found in Master ID list' % MolID_BuildingBlock
                MolIDs_Refined.append(MolID_BuildingBlocks)
                Coeffs_Refined.append(Coeff_BuildingBlocks)
            else:
                print('Molecule ID not defined in the organism: %s' % MolID)
        MolIDs_Refined = self.FlatList(MolIDs_Refined)
        Coeffs_Refined = self.FlatList(Coeffs_Refined)
        return MolIDs_Refined, Coeffs_Refined

    def ParseRXNRate(self, Rate):
        # Parse Rxn rate (to be expanded)
        RXNRateParsed = Rate
        return RXNRateParsed

    def FlatList(self, List):
        ListFlattened = []
        for i in List:
            if isinstance(i, str):
                ListFlattened.append(i)
            elif isinstance(i, int):
                ListFlattened.append(i)
            elif isinstance(i, float):
                ListFlattened.append(i)
            else:
                for j in i:
                    ListFlattened.append(j)
        return ListFlattened

    def GetMolIdx(self, Molecules, MolIdxRef):
        MolIdxList = []

        # For ID2Idx cases
        if isinstance(MolIdxRef, dict):
            for Molecule in Molecules:
                MolIdx = MolIdxRef[Molecule]
                MolIdxList.append(MolIdx)
            return MolIdxList
        # For Type2Idx cases
        elif isinstance(MolIdxRef, list):
            for MolIdx, Type in enumerate(MolIdxRef):
                if Type == Molecules:
                    MolIdxList.append(MolIdx)
            return MolIdxList
        else:
            print("Inappropriate reference type used in GetMolIdx function parameter: %s" % MolIdxRef)

    def ReverseNumberSign(self, ListOfNumbers):
        NewListOfNumbers = [-Number for Number in ListOfNumbers]
        return NewListOfNumbers

    def GenerateCellProcessInterface(self, Writer):
        with Writer.Statement("class FCellProcess():"):
            with Writer.Statement("def __init__(self):"):
                Writer.BlankLine()

                Writer.Statement("super().__init__()")
                Writer.BlankLine()

            Writer.AbsMethod()
            with Writer.Statement("def AddToStoichiometryMatrix(self):"):
                Writer.Pass_____()
                Writer.BlankLine()

            Writer.AbsMethod()
            with Writer.Statement("def CalculateRate(self):"):
                Writer.Pass_____()
                Writer.BlankLine()

            Writer.AbsMethod()
            with Writer.Statement("def AddToRateMatrix(self):"):
                Writer.Pass_____()
                Writer.BlankLine()

            Writer.TF_Graph_()
            with Writer.Statement("def ClearRateMatrix(self):"):
                Writer.Statement("self.Cel.ClearRateMatrix()")
                Writer.BlankLine()

            Writer.TF_Graph_()
            with Writer.Statement("def SetUpStoichiometryMatrix(self):"):
                Writer.Statement("self.AddToStoichiometryMatrix()")
                Writer.BlankLine()

            Writer.TF_Graph_()
            with Writer.Statement("def UpdateRates(self):"):
                Writer.Statement("self.ClearRateMatrix()")
                Writer.Statement("self.CalculateRate()")
                Writer.Statement("self.AddToRateMatrix()")
                Writer.BlankLine()

    def GenerateCellProcess(self, Writer, ProcessID, Reactions):
        with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
            with Writer.Statement("def __init__(self):"):
                Writer.BlankLine()

                Writer.Statement("super().__init__()")
                Writer.BlankLine()

            with Writer.Statement("def AddToStoichiometryMatrix(self):"):
                for Reaction in Reactions:
                    Stoichiometry, Rate = Reaction
                    MolIDs, MolIdxs, Coeffs = Stoichiometry

                    # Convert index and stoichiometry to tensor for simulation
                    Writer.Variable_("MolIdxs", MolIdxs)
                    Writer.Variable_("Coeffs", Coeffs)

                    Writer.Comment__("Process Reaction Added:")
                    Writer.Comment__("MolID\t:\tMolIdx\t:\tCoeff)")

                    for MolID, MolIdx, Coeff in zip(MolIDs, MolIdxs, Coeffs):
                        Writer.Comment__("%s\t:\t%s\t:\t%s" % (MolID, MolIdx, Coeff))

                    Writer.BlankLine()

                    # Initialize RXN stoichiometry array for all molecules
                    Writer.InitZeros("StoichiometryArray", self.Comp.Master.NUniq_Master)

                    # Update RXN stoichiometry array with Idx and Stoich participating in the RXN
                    Writer.OperScUpd("StoichiometryArray", "MolIdxs", "Coeffs")
                    Writer.Reshape__("StoichiometryArray", "StoichiometryArray", [1, -1])

                    # Add StoichiometryArray to the Cel.MX_Stoichiometries matrix
                    Writer.Statement("self.Cel.AddToStoichiometryMatrix(StoichiometryArray)")
                    Writer.BlankLine()

            with Writer.Statement("def CalculateRate(self):"):
                Writer.Pass_____()
                Writer.BlankLine()

            with Writer.Statement("def AddToRateMatrix(self):"):
                for Reaction in Reactions:
                    Stoichiometry, Rate = Reaction
                    # Calculate RXN rate
                    Writer.Statement("Rate = self.CalculateRate()")
                    Writer.Statement("Rate = %s" % Rate)  # Temporary code
                    Writer.Reshape__("Rate", "Rate", [1, 1])
                    Writer.BlankLine()

                    # Add RXN to the rate matrix
                    Writer.Statement("self.Cel.AddToRateMatrix(Rate)")
                    Writer.BlankLine()
