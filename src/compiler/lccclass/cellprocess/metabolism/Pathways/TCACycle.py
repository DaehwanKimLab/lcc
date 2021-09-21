# TCA Cycle

import numpy as np

'''
'''

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):
    """

    :type Comp: object
    """
    # Parse flat data

    TCARXNs = list()



#
#     MolsInTCARXN =
#
#     self.ID_METRXNs = list()
#     self.ID2Idx_METRXNs = dict()
#     self.Dir_RevMETRXNs = 0  # Rev for Reversibility
#     self.ID_Enzymes4METRXN = list()
#     self.ID_Substrates4METRXN = list()
#     self.ID2Idx_Enzymes4METRXN = dict()
#
#     self.ID_MolsInMETRXN = list()
#     self.ID2Idx_MolsInMETRXN = dict()
#     self.Coeff_MolsInMETRXN = list()
#     self.NUniq_MolsInMETRXN = 0
#
#     self.NUniq_METRXNs = 0
#
#     self.IDModification = dict()
#
#     # Unregistered Metabolites
#     self.ID_UnregisteredFromMETRXN_Metabolites = list()  # Mostly additional metabolites
#     self.MolType_UnregisteredFromMETRXN_Metabolites = 'UnregisteredMetabolitesFromMETRXN'
#     self.Count_UnregisteredFromMETRXN_Metabolites = 0
#     self.MW_UnregisteredFromMETRXN_Metabolites = 0
#     self.NUniq_UnregisteredFromMETRXN_Metabolites = 0
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#     Idx_NADH = Comp.Master.ID2Idx_Master['NADH[c]']
#     Idx_NADPH = Comp.Master.ID2Idx_Master['NADPH[c]']
#     Idx_FADH2 = Comp.Master.ID2Idx_Master['FADH2[p]']
#
#     # For the reactions with kinetic data
#     Idx_SubstratesKINRXN = ProGen.GetMolIdx_Master(Comp.Kinetics.ID_Substrates4KINRXN)
#     Idx_EnzymesKINRXN = ProGen.GetMolIdx_Master(Comp.Kinetics.ID_Enzymes4KINRXN)
#
#     # For the reactions without kinetic data
#     Idx_SubstratesMETRXN = ProGen.GetMolIdx_Master(Comp.Metabolism.ID_Substrates4METRXN)
#     Idx_EnzymesMETRXN = ProGen.GetMolIdx_Master(Comp.Metabolism.ID_Enzymes4METRXN)   # Not for all METRXNs
#
#     # Enzymes available for a subset of metabolism reactions
#     Idx_METRXNs_WithEnzymeWithoutKineticData = list()
#     for EnzymeID in Comp.Metabolism.ID_Enzymes4METRXN:
#         Idx_METRXNs_WithEnzymeWithoutKineticData.append(Comp.Metabolism.ID2Idx_Enzymes4METRXN[EnzymeID])
#
#     Idx_METRXNsWithKineticData = list()
#     for ID_KINRXN in Comp.Kinetics.ID_KINRXNs:
#         Idx_METRXNsWithKineticData.append(Comp.Metabolism.ID2Idx_METRXNs[ID_KINRXN])
#
#     # Kcat and Km default values for reactions without kinetic data have been set in Dataset.py
#     # This serves as a local control in cell.py
#     Kcat_Default = Comp.Kinetics.Const_Kcat_Default
#     Km_Default = Comp.Kinetics.Const_Km_Default
#
#     NUniq_METRXN = Comp.Metabolism.NUniq_METRXNs
#
#     # DEBUG COUNTS:
#     IdxsToDebug = [[14766], [14771], [14785], [14806], [14813], [14814], [14821], [14839], [14842], [14856], [14857], [14859], [14865], [14866], [14869], [14880], [14881], [14887], [14889], [14892], [14923], [14962], [14985], [15012], [15030], [15031], [15051], [15053], [15057], [15060], [15066], [15077], [15081], [15083], [15085], [15105], [15125], [15148], [15162], [15178], [15179], [15183], [15188], [15211], [15246], [15254], [15290], [15360], [15363], [15369], [15384], [15389], [15390], [15399], [15401], [15419], [15426], [15430], [15433], [15435], [15436], [15438], [15439], [15440], [15441], [15445], [15448], [15452], [15457], [15468], [15479], [15481], [15483], [15485], [15503], [15508], [15523], [15524], [15545], [15572], [15574], [15579], [15600], [15606], [15613], [15621], [15623], [15637], [15651], [15674], [15679], [15729], [15745], [15746], [15748], [15755], [15763], [15768], [15820], [15829], [15841], [15844], [15848], [15854], [15878], [15891], [15898], [15903], [15914], [15919], [15922], [15993], [16019], [16069], [16080], [16133], [16303], [16304], [16311]]
#     IDsToDebug = list()
#     for Idx in IdxsToDebug:
#         IDsToDebug.append(Comp.Master.ID_Master[Idx[0]])
#
#     MolsToAnalyze = ['dNTPs', 'NTPs', 'AAs', 'NADH', 'NADPH', 'FADH2']
#
#     # TCA cycle I
#     # ID2ID_BCID2EC_TCACycle =
#

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)
#
#         with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
#             Writer.Variable_("self.Count_MetabolitesInitial", 0)
#
#             Writer.Variable_("self.Idx_Enzymes", 0)
#             Writer.Variable_("self.Idx_Substrates", 0)
#             Writer.Variable_("self.Idx_EnzymesMETRXN", 0)
#             Writer.Variable_("self.Idx_SubstratesMETRXN", 0)
#             Writer.Variable_("self.Idx_METRXNsWithKineticData", 0)
#             Writer.Variable_("self.Idx_METRXNs_WithEnzymeWithoutKineticData", 0)
#             Writer.BlankLine()
#
#             Writer.Variable_("self.Const_Kcat_Default", Kcat_Default)
#             Writer.Variable_("self.Const_Km_Default", Km_Default)
#             Writer.BlankLine()
#
#             Writer.Variable_("self.Coeff_Metabolism_Transposed", 0)
#             Writer.BlankLine()
#
#             Writer.Variable_("self.IdxsToDebug", 0)
#             Writer.BlankLine()
#
#         with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):
#             Writer.Variable_("self.Cel.Idx_NADH", Idx_NADH)
#             Writer.Variable_("self.Cel.Idx_NADPH", Idx_NADPH)
#             Writer.Variable_("self.Cel.Idx_NADPH", Idx_FADH2)
#             Writer.BlankLine()
#
#             Writer.Variable_("self.Idx_EnzymesKINRXN", Idx_EnzymesKINRXN)
#             Writer.Variable_("self.Idx_SubstratesKINRXN", Idx_SubstratesKINRXN)
#             Writer.BlankLine()
#             Writer.Variable_("self.Idx_EnzymesMETRXN", Idx_EnzymesMETRXN)
#             Writer.Variable_("self.Idx_SubstratesMETRXN", Idx_SubstratesMETRXN)
#             Writer.Variable_("self.Idx_METRXNsWithKineticData", Idx_METRXNsWithKineticData)
#             Writer.Variable_("self.Idx_METRXNs_WithEnzymeWithoutKineticData", Idx_METRXNs_WithEnzymeWithoutKineticData)
#             Writer.BlankLine()
#             Writer.Transpose("self.Coeff_Metabolism_Transposed", "self.Cel.Coeff_Metabolism")
#             Writer.Cast_____("self.Coeff_Metabolism_Transposed", "self.Coeff_Metabolism_Transposed", 'float32')
#             Writer.BlankLine()
#
#             Writer.Comment__("For Analysis Only")
#             Writer.Variable_("self.IdxsToDebug", IdxsToDebug)
#             Writer.BlankLine()
#
#             Writer.Comment__("For Replenishing Metabolites")
#             Writer.Statement("self.GetMetaboliteCounts()")
#             Writer.BlankLine()
#
#
#         with Writer.Statement("def GetMetaboliteCounts(self):"):
#             Writer.Gather___("self.Count_MetabolitesInitial", "self.Cel.Count_Master", "self.Cel.Idx_Master_Metabolites")
#             Writer.BlankLine()
#
#         with Writer.Statement("def ReplenishMetabolites(self):"):
#             Writer.ScatNdUpd("self.Cel.Counts", "self.Cel.Counts", "self.Cel.Idx_Master_Metabolites", "self.Count_MetabolitesInitial")
#             Writer.BlankLine()
#
#         with Writer.Statement("def Message_ReplenishMetabolites(self):"):
#             Writer.PrintStrg("\t- Metabolites have been replenished.")
#             Writer.BlankLine()
#
#
#
#         with Writer.Statement("def ExecuteProcess(self):"):
#             if Writer.Switch4Kinetics:
#                 # Calculate reaction rate for the reactions with kinetic data
#                 Writer.Statement("Rate_KINRXN = self.CalculateRateForReactionsWithKineticData()")
#                 Writer.Statement("Rate_METRXN = self.CalculateRateForReactionsWithoutKineticData()")
#                 Writer.BlankLine()
#                 Writer.ScatNdUpd("Rate_METRXN", "Rate_METRXN", "self.Idx_METRXNsWithKineticData", "Rate_KINRXN", Type='float32')
#                 Writer.Reshape__("Rate_METRXN", "Rate_METRXN", [-1, 1])
#                 Writer.MatrixMul("Count_MolsInMETRXN", "self.Coeff_Metabolism_Transposed", "Rate_METRXN")
#                 Writer.RoundInt_("Count_MolsInMETRXN", "Count_MolsInMETRXN")
#                 Writer.BlankLine()
#                 Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_Master_MolsInMETRXN, Count_MolsInMETRXN)")
#                 Writer.BlankLine()
#                 Writer.Comment__("For Analysis Only")
#                 Writer.Statement("self.ForAnalysis(Count_MolsInMETRXN)")
#                 Writer.BlankLine()
#
#             else:
#                 Writer.Pass_____()
#                 Writer.BlankLine()
#
#         with Writer.Statement("def CalculateRateForReactionsWithKineticData(self):"):
#             Writer.Statement("Count_Enzymes = self.GetCounts_Float(self.Idx_EnzymesKINRXN)")
#             Writer.Statement("Count_Substrates = self.GetCounts_Float(self.Idx_SubstratesKINRXN)")
#             Writer.BlankLine()
#             Writer.Statement("Rate_KINRXN = self.CalculateReactionRate(Count_Enzymes, Count_Substrates, self.Cel.Const_Kcat, self.Cel.Const_Km)")
#             Writer.ReturnVar("Rate_KINRXN")
#             Writer.BlankLine()
#
#         # TODO: the reactions without kinetic data need linear algebra
#         with Writer.Statement("def CalculateRateForReactionsWithoutKineticData(self):"):
#             # Currently: calculate reaction rate for the reactions without kinetic data using average values
#             Writer.Statement("Count_Enzymes = self.GetCounts_Float(self.Idx_EnzymesMETRXN)")
#             Writer.ReduceMin("Count_EnzymesAverage", "Count_Enzymes")
#             # Writer.ReduceAvg("Count_EnzymesAverage", "Count_Enzymes")
#
#             Writer.Statement("Count_Substrates = self.GetCounts_Float(self.Idx_SubstratesMETRXN)")
#             Writer.Gather___("Count_Substrates_WithEnzymes", "Count_Substrates", "self.Idx_METRXNs_WithEnzymeWithoutKineticData")
#             Writer.Reshape__("Count_Substrates_WithEnzymes", "Count_Substrates_WithEnzymes", -1)
#             Writer.BlankLine()
#             Writer.Statement("Rate_METRXN_WithoutEnzymes = self.CalculateReactionRate(Count_EnzymesAverage, Count_Substrates, self.Const_Kcat_Default, self.Const_Km_Default)")
#             Writer.Statement("Rate_METRXN_WithEnzymeWithoutKineticData = self.CalculateReactionRate(Count_Enzymes, Count_Substrates_WithEnzymes, self.Const_Kcat_Default, self.Const_Km_Default)")
#
#             Writer.ScatNdUpd("Rate_METRXN", "Rate_METRXN_WithoutEnzymes", "self.Idx_METRXNs_WithEnzymeWithoutKineticData", "Rate_METRXN_WithEnzymeWithoutKineticData", Type='float32')
#             Writer.ReturnVar("Rate_METRXN")
#             Writer.BlankLine()
#
#         with Writer.Statement("def CalculateReactionRate(self, Count_Enzymes, Count_Substrates, Const_Kcat, Const_Km):"):
#             # (kcat * E * S) / (S + kM)
#             Writer.Multiply_("Numerator", "Const_Kcat", "Count_Enzymes")
#             Writer.Multiply_("Numerator", "Numerator", "Count_Substrates")
#             Writer.Add______("Denominator", "Count_Substrates", "Const_Km")
#             Writer.Divide___("Rate_RXN", "Numerator", "Denominator")
#             Writer.ReturnVar("Rate_RXN")
#             Writer.BlankLine()
#
#         with Writer.Statement("def ForAnalysis(self, Count_MolsInMETRXN):"):
#             Writer.ZerosLike("DeltaCountFromMetabolism", "self.Cel.Counts", 'int32')
#             Writer.ScatNdAdd("DeltaCountFromMetabolism", "DeltaCountFromMetabolism", "self.Cel.Idx_Master_MolsInMETRXN", "Count_MolsInMETRXN")
#             for Mol in MolsToAnalyze:
#                 Writer.Gather___("self.DeltaCounts_%s" % Mol, "DeltaCountFromMetabolism", "self.Cel.Idx_%s" % Mol)
#             Writer.Gather___("self.DeltaCounts_IDsToDebug", "DeltaCountFromMetabolism", "self.IdxsToDebug")
#             Writer.BlankLine()
#
#         with Writer.Statement("def ViewProcessSummary(self):"):
#             if Writer.Switch4Kinetics:
#                 Writer.PrintStrg("===== Metabolism ===== ")
#                 for Mol in MolsToAnalyze:
#                     Writer.PrintStVa("deltaCount of %s" % Mol, "self.DeltaCounts_%s" % Mol)
#                 Writer.PrintStVa("deltaCount of IDsToDebug", "self.DeltaCounts_IDsToDebug")
#                 Writer.BlankLine()
#             else:
#                 Writer.Pass_____()
#                 Writer.BlankLine()
#
# # def SetUpReactions(ProGen):
# #     Reactions = list()
# #
# #     Reaction = None
# #
# #     Reaction_SetUp = Reaction
# #     Reactions.append(Reaction_SetUp)
# #
# #     return Reactions
# #
#
