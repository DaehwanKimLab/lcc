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
import os, sys
import abc
import CodeGen
from lccclass.cellprocess.signaling import SigmaFactor
from lccclass.cellprocess.synthesis import Replication, Translation, Transcription
from lccclass.cellprocess.degradation import ProteinDegradation, RNADegradation, DNADegradation
from lccclass.cellprocess.conversion import Complexation, Equilibrium
from lccclass.cellprocess.biochemicalreactions import Metabolism
# from lccclass.cellprocess.biochemicalreactions import Transport
from lccclass.cellprocess.biochemicalreactions.pathways import TCACycle
from lccclass.cellprocess.biochemicalreactions.pathways import Magnesium
from lccclass.cellprocess.biochemicalreactions.pathways import CarbonicAnhydrase
from lccclass.cellprocess.division import CellDivision


class FProcessGenerator():
    def __init__(self):
        self.CellProcesses = list()
        self.Dict_CellProcesses = dict()
        # self.Dict_ProcessReactions = dict()
        # self.Dict_ReactionTypeReference = dict()

        # self.Reactions_Replication = list()
        # self.Reactions_Transcription = list()
        # self.Reactions_Translation = list()

        # self.SetUpReactionTypeReference()

    def LinkCompilerObj(self, Comp):
        self.Comp = Comp
        self.ImportSaveVariables()

    def ImportSaveVariables(self):
        self.Switch4SaveAllData = self.Comp.Switch4SaveAllData
        self.SavePath = self.Comp.SavePath

    def SetProcessList(self):
        CellProcessesAvailable = [
            # Signaling
            [SigmaFactor, "SigmaFactor"],

            # Biosynthesis
            [Replication, "Replication"],
            [Transcription, "Transcription"],
            [Translation, "Translation"],

            # Conversion
            [Complexation, "Complexation"],
            # [Equilibrium, "Equilibrium"],

            # Modification

            # Degradation
            [RNADegradation, "RNADegradation"],
            [ProteinDegradation, "ProteinDegradation"],

            # Metabolism
            [Metabolism, "Metabolism"],
            # [TCACycle, "TCACycle"],
            # [Transport, "Transport"],
            # [Magnesium, "Magnesium"],

            # Cell division
            [CellDivision, "CellDivision"],
        ]

        if not self.Comp.UserInput.CellProcesses or "All" in self.Comp.UserInput.CellProcesses:
            self.CellProcesses = list([CellProcess[0] for CellProcess in CellProcessesAvailable])

        else:
            for CellProcess in CellProcessesAvailable:
                if CellProcess[1] in self.Comp.UserInput.CellProcesses:
                    self.CellProcesses.append(CellProcess[0])


    def PrintProcessID(self, ProcessID):
        print("Cellular process written in simulation code: %s" % ProcessID)

    # def PrintCompletion(self, ProcessID):
    #     print("Reactions are successfully set up. Process: %s" % ProcessID)

    # def SetUpReactionTypeReference(self):
    #     self.Dict_ReactionTypeReference = {
    #         'Biochemical Reaction': 'Bch',
    #         'Polymerization': 'Pol'
    #     }

    # def PrepareReactionForMatrix(self, Reaction):
    #     # Put all reaction info in the lists
    #     for Type, RXNID, MoleculeIDs, Stoichiometry, Reversibility, Order, Rate, ModulatorIDs, Trigger in Reaction:
    #         MoleculeIDs, Stoichiometry = Utils.RefineBuildingBlocks(MoleculeIDs, Stoichiometry, self.Comp)
    #         return MoleculeIDs, Stoichiometry
    #

    def SetUpProcesses(self):

        for Process in self.CellProcesses:
            ProcessStr = Process.__name__.split('.')[-1]
            self.Dict_CellProcesses[ProcessStr] = Process
            self.PrintProcessID(ProcessStr)

            # self.Dict_ProcessReactions[ProcessStr] = Process.SetUpReactions(self)
            # self.PrintCompletion(ProcessStr)
        #
        # self.Dict_CellProcesses['Replication'] = Replication
        # self.Reactions_Replication = Replication.SetUpReactions(self)
        # self.Dict_ProcessReactions['Replication'] = self.Reactions_Replication
        #
        # self.Dict_CellProcesses['Transcription'] = Transcription
        # self.Reactions_Transcription = Transcription.SetUpReactions(self)
        # self.Dict_ProcessReactions['Transcription'] = self.Reactions_Transcription

        # self.Reactions_Translation = Translation.SetUpReactions(self)
        # self.Dict_Reactions['Translation'] = self.Reactions_Translation

    def SaveProcesses(self):
        # Save all self.Reactions
        if self.Switch4SaveAllData:
            for ProcessName, Reactions in self.Dict_ProcessReactions.items():
                SaveFileName = "%s/%s" % (self.SavePath, ProcessName)
                np.save('%s.npy' % SaveFileName, Reactions)
                print('Process data have been successfully saved. Process: "%s"' % ProcessName)

    def BuildingBlockIDs(self, BuildingBlocks):
        BuildingBlocksExpanded = None
        if BuildingBlocks == 'dNTPs':
            BuildingBlocksExpanded = self.Comp.BuildingBlock.Name_dNTPs
        elif BuildingBlocks == 'NTPs':
            BuildingBlocksExpanded = self.Comp.BuildingBlock.Name_NTPs
        elif BuildingBlocks == 'AAs':
            BuildingBlocksExpanded = self.Comp.BuildingBlock.Name_AAs
        else:
            BuildingBlocksExpanded = BuildingBlocks
            # print("BuildingBlocks do not exist: %s" % BuildingBlocks)
        return BuildingBlocksExpanded

    def BuildingBlockIdxs(self, BuildingBlocks):
        BuildingBlocksExpanded = self.BuildingBlockIDs(BuildingBlocks)
        BuildingBlocksIdxs = self.GetMolIdx(BuildingBlocksExpanded, self.Comp.Master.ID2Idx_Master)
        return BuildingBlocksIdxs




    def SetUpType(self, Type):
        Type_Parsed = self.Dict_ReactionTypeReference[Type]
        return Type_Parsed

    def SetUpStoichiometry(self, Stoichiometry):
        MolIDs, Coeffs = Stoichiometry
        MolIDs_Parsed, Coeffs_Parsed = self.RefineBuildingBlocks(MolIDs, Coeffs)
        MolIdxs = self.GetMolIdx(MolIDs_Parsed, self.Comp.Master.ID2Idx_Master)
        return [MolIDs_Parsed, MolIdxs, Coeffs_Parsed]

    def SetUpRate(self, Rate):
        Rate_Parsed = Rate
        return Rate_Parsed

    def SetUpTrigger(self, Trigger):
        MolIDs, Thresholds, Conditions = Trigger
        MolIDs_Parsed = MolIDs # Placeholder line
        MolIdxs = self.GetMolIdx(MolIDs, self.Comp.Master.ID2Idx_Master)
        Thresholds_Parsed = Thresholds # Placeholder line
        Conditions_Parsed = Conditions  # Placeholder line
        return [MolIDs_Parsed, MolIdxs, Thresholds_Parsed, Conditions_Parsed]


    def SetUpReaction(self, Reaction):
        '''

        :param Reaction: Type, Rate, Sto, Tri
        :return:
        '''
        Reaction_SetUp = dict()
        Reaction_SetUp['Type'] = self.SetUpType(Reaction['Type'])
        Reaction_SetUp['Stoichiometry'] = self.SetUpStoichiometry(Reaction['Stoichiometry'])
        Reaction_SetUp['Rate'] = self.SetUpRate(Reaction['Rate'])
        # Reaction_SetUp['Trigger'] = self.SetUpTrigger(Reaction['Trigger'])
        return Reaction_SetUp

    # def SetUpReactions(self, Reactions):
    #     for Reaction in Reactions:
    #         Reaction = self.SetUpReaction(Reaction)

    # Reaction Matrix Building functions
    def ParseBuildingBlocks(self, MolGroup, Stoichiometry):
        MolGroupParsed = list()
        StoichiometryParsed = list()
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
        MolIDs_Refined = list()
        Coeffs_Refined = list()
        ReactionOrder_Refined = list()
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
        ListFlattened = list()
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
        MolIdxList = list()

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

    def GetMolIdx_Master(self, Molecules):
        return self.GetMolIdx(Molecules, self.Comp.Master.ID2Idx_Master)

    def ReverseNumberSign(self, ListOfNumbers):
        NewListOfNumbers = [-Number for Number in ListOfNumbers]
        return NewListOfNumbers

    def FindMolID(self, MolName, Type=None, Localization='[c]'):
        # Type can be 'Gene' or 'RNA' or 'Protein' or 'Complex' or 'Metabolite'
        MolID = None
        if Type:
            if MolName in self.Comp.Metabolite.ID_Metabolites:
                MolID = MolName
            elif MolName in self.Comp.BuildingBlock.Name2Key_BuildingBlocks.keys():
                MolID = MolName
            elif Type == 'Gene':
                MolID = self.Comp.Dict4ID_Gene[MolName]
            elif Type == 'RNA':
                MolID = self.Comp.Dict4ID_RNA[MolName]
            elif Type == 'Protein':
                MolID = self.Comp.Dict4ID_Protein[MolName]
            else:
                # print('MolName %s does not match any MolID in Compiler Data' % MolName, sys.exc_info()[0])
                MolID = MolName
        else:
            if MolName in self.Comp.Dict4ID_Gene:
                MolID = self.Comp.Dict4ID_Gene[MolName]
            elif MolName in self.Comp.Dict4ID_RNA:
                MolID = self.Comp.Dict4ID_RNA[MolName]
            elif MolName in self.Comp.Dict4ID_Protein:
                MolID = self.Comp.Dict4ID_Protein[MolName]
            else:
                # print('MolName %s does not match any MolID in Compiler Data' % MolName, sys.exc_info()[0])
                MolID = MolName
        return MolID


    # Generate Cell Process Methods
    def Init_Common(self, Writer):
        with Writer.Statement("def __init__(self, Cel, Cst, Env, Exe):"):
            Writer.Statement("super().__init__(Cel, Cst, Env, Exe)")
            Writer.Statement("self.Init_ProcessSpecificVariables()")
            Writer.BlankLine()



    # def GenerateCellProcessInterface(self, Writer):
    #     with Writer.Statement("class FCellProcess():"):
    #         with Writer.Statement("def __init__(self, Bch, Cel, Cst, Env, Exe, Pol):"):
    #             Writer.LinkClObj('Cel')
    #             Writer.LinkClObj('Cst')
    #             Writer.LinkClObj('Env')
    #
    #             Writer.LinkClObj('Bch')
    #             Writer.LinkClObj('Pol')
    #             Writer.LinkClObj('Exe')
    #             Writer.BlankLine()
    #
    #             Writer.Statement("super().__init__()")
    #             Writer.BlankLine()
    #
    #         Writer.AbsMethod()
    #         with Writer.Statement("def AddToStoichiometryMatrix(self):"):
    #             Writer.Pass_____()
    #             Writer.BlankLine()
    #
    #         Writer.AbsMethod()
    #         with Writer.Statement("def CalculateRate(self):"):
    #             Writer.Pass_____()
    #             Writer.BlankLine()
    #
    #         Writer.AbsMethod()
    #         with Writer.Statement("def AddToRateMatrix(self):"):
    #             Writer.Pass_____()
    #             Writer.BlankLine()
    #
    #         Writer.AbsMethod()
    #         with Writer.Statement("def PrintMolCounts(self):"):
    #             Writer.Pass_____()
    #             Writer.BlankLine()
    #
    #         Writer.TF_Graph_()
    #         with Writer.Statement("def SetUpStoichiometryMatrix(self):"):
    #             Writer.Statement("self.AddToStoichiometryMatrix()")
    #             Writer.BlankLine()
    #
    #         Writer.TF_Graph_()
    #         with Writer.Statement("def UpdateRates(self):"):
    #             Writer.Statement("self.CalculateRate()")
    #             Writer.Statement("self.AddToRateMatrix()")
    #             Writer.BlankLine()

    # def GenerateCellProcess(self, Writer, ProcessID):
    #     with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
    #         with Writer.Statement("def __init__(self, Bch, Cel, Cst, Env, Exe, Pol):"):
    #             Writer.Statement("super().__init__(Bch, Cel, Cst, Env, Exe, Pol)")
    #             Writer.Statement("self.Init_ProcessSpecificVariables()")
    #             Writer.BlankLine()


            # with Writer.Statement("def AddToStoichiometryMatrix(self):"):
            #     for Reaction in Reactions:
            #
            #         if Reaction is None:
            #             Writer.Pass_____()
            #             Writer.BlankLine()
            #             continue
            #
            #         MolIDs, MolIdxs, Coeffs = Reaction['Stoichiometry']
            #
            #         # Convert index and stoichiometry to tensor for simulation
            #         Writer.Comment__("MolIDs = %s" % MolIDs)
            #         Writer.Variable_("MolIdxs", MolIdxs)
            #         Writer.Variable_("Coeffs", Coeffs)
            #         Writer.BlankLine()
            #
            #         # Initialize RXN stoichiometry array for all molecules
            #         Writer.InitZeros("StoichiometryArray", self.Comp.Master.NUniq_Master)
            #
            #         # Update RXN stoichiometry array with Idx and Stoich participating in the RXN
            #         Writer.ScatNdUpd("StoichiometryArray", "MolIdxs", "Coeffs")
            #         Writer.Reshape__("StoichiometryArray", "StoichiometryArray", [1, -1])
            #
            #         # Add StoichiometryArray to the Cel.MX_Stoichiometries matrix
            #         Writer.Statement("self.Cel.AddToStoichiometryMatrix(StoichiometryArray)")
            #         Writer.BlankLine()
            #
            # with Writer.Statement("def CalculateRate(self):"):
            #     Writer.Pass_____()
            #     Writer.BlankLine()
            #
            # with Writer.Statement("def AddToRateMatrix(self):"):
            #     for Reaction in Reactions:
            #
            #         if Reaction is None:
            #             Writer.Pass_____()
            #             Writer.BlankLine()
            #             continue
            #
            #         # Initialize Rate value to 0
            #         Writer.Variable_("Rate", 0)
            #
            #         # Reaction type dependent rate calculation
            #         Type = Reaction['Type']
            #
            #         if Type == 'Pol':
            #             MolIDs, MolIdxs, Thresholds, Conditions = Reaction['Trigger']
            #             for MolID, MolIdx, Threshold, Condition in zip(MolIDs, MolIdxs, Thresholds, Conditions):
            #                 # Generate an evaluation phrase
            #                 Trigger = "self.Cel.Counts[%s] " % MolIdx + Condition + ' ' + Threshold
            #                 with Writer.Statement("if " + Trigger + ":  # Evaluate Molecule ID: \"%s\"" % MolID):
            #
            #                     # Update rate only if this is the first time to determine rate for the current reaction
            #                     with Writer.Statement("if Rate == 0:"):
            #
            #                         Rate_Mean, Rate_SD, Rate_UnitTime = Reaction['Rate']
            #                         Writer.Variable_("Rate_Mean", Rate_Mean)
            #                         Writer.Variable_("Rate_SD", Rate_SD)
            #
            #                         Writer.Statement("self.Pol.LoadRateMean(Rate_Mean)")
            #                         Writer.Statement("self.Pol.LoadRateSD(Rate_SD)")
            #
            #                         Writer.Statement("Rate = self.Pol.DetermineRate()")
            #
            #                 with Writer.Statement("else:"):
            #                     Writer.Variable_("Rate", 0)
            #                 Writer.BlankLine()
            #
            #         elif Type == 'Bch':
            #             # To be implemented
            #             pass
            #
            #         Writer.Cast_____("Rate", "Rate", 'float32')
            #         Writer.Reshape__("Rate", "Rate", [1, 1])
            #         Writer.BlankLine()
            #
            #         # Add RXN to the rate matrix
            #         Writer.Statement("self.Cel.AddToRateMatrix(Rate)")
            #         Writer.BlankLine()
            #
            # with Writer.Statement("def DisplayCounts(self):"):
            #     Dict_MolIdx_ID_Pairs = dict()
            #     for Reaction in Reactions:
            #
            #         if Reaction is None:
            #             Writer.PrintStrg("No reactions found in %s" % ProcessID)
            #             Writer.BlankLine()
            #             continue
            #
            #         MolIDs, MolIdxs, Coeffs = Reaction['Stoichiometry']
            #         for MolID, MolIdx in zip(MolIDs, MolIdxs):
            #             if MolIdx not in Dict_MolIdx_ID_Pairs.keys():
            #                 Dict_MolIdx_ID_Pairs[MolIdx] = MolID
            #
            #     for MolIdx, MolID in sorted(Dict_MolIdx_ID_Pairs.items()):
            #
            #         if Reaction is None:
            #             Writer.Pass_____()
            #             Writer.BlankLine()
            #             continue
            #
            #         Writer.Statement("print('%s Count:\t', self.Cel.Counts[%s])" % (MolID, MolIdx))
            #     Writer.BlankLine()
            #
            #
            # # Class-specific methods
            # if ProcessID == 'Metabolism':
            #     with Writer.Statement("def GetMetaboliteIdxsAndCounts(self):"):
            #         MetaboliteIdxs = list()
            #         MetaboliteCountsInitial = list()
            #         for MolIdx in range(self.Comp.Master.NUniq_Master):
            #             if self.Comp.Master.Type_Master[MolIdx] == 'Metabolite':
            #                 MetaboliteCountsInitial.append(self.Comp.Master.Count_Master[MolIdx])
            #                 MetaboliteIdxs.append(MolIdx)
            #
            #         Writer.Variable_("self.MetaboliteIdxs", MetaboliteIdxs)
            #         Writer.Reshape__("self.MetaboliteIdxs", "self.MetaboliteIdxs", [1, -1])
            #         Writer.Statement("self.Exe.LoadMetaboliteIdxs(self.MetaboliteIdxs)")
            #
            #         Writer.Variable_("self.MetaboliteCountsInitial", MetaboliteCountsInitial)
            #         Writer.Reshape__("self.MetaboliteCountsInitial", "self.MetaboliteCountsInitial", [-1])
            #         Writer.Statement("self.Exe.LoadMetaboliteCountsInitial(self.MetaboliteCountsInitial)")
            #         Writer.BlankLine()
