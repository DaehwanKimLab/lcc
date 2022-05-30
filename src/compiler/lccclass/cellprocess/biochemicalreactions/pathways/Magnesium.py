# Magnesium

import numpy as np

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):
    """

    :type Comp: object
    """
    '''
    # Ion regulation
    IonsOfInterest = ['K+[c]', 'CL-[c]', 'MG+2[c]', 'FE+2[c]']  # wcs does not have NA+[c], and current re
    IdxOfTRNRXNs = dict()

    for Ion in IonsOfInterest:
        Idx_Ion = Comp.Transporters.ID2Idx_MolsInTRNRXN[Ion]
        Idx_TRNRXNs = list()
        for i, CoeffArray in enumerate(Comp.Transporters.Coeff_MolsInTRNRXN):
            if CoeffArray[Idx_Ion] != 0:
                Idx_TRNRXNs.append(i)
                # print('%s RXN \t: %s' % (Ion, Comp.Transporters.ID2EQN_TRNRXNs[Comp.Transporters.ID_TRNRXNs[i]]))
        IdxOfTRNRXNs[Ion] = Idx_TRNRXNs
    '''
    # Indices
    Idx_SubstratesInMgRXN = ProGen.GetMolIdx_Master(Comp.Transporters.ID_MolsInMgRXN)
    Idx_MgtA = Comp.Master.ID2Idx_Master[Comp.Transporters.ID_Enzymes4MgRXN[2]]
    Idx_Mg_Cytosol = Comp.Master.ID2Idx_Master['MG+2[c]']
    Idx_Mg_Periplasm = Comp.Master.ID2Idx_Master['MG+2[p]']

    Coeff_MgRXNs = Comp.Transporters.Coeff_MolsInMgRXN
    N_MgRXNs = len(Comp.Transporters.ID_MgRXNs)

    Kcat_Default = 0
    Km_Default = 0

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Function_("Init_ProcessSpecificVariables"):
            # Writer.Variable_("self.Idx_Enzymes", 0)
            # Writer.Variable_("self.Idx_Substrates", 0)
            # Writer.Variable_("self.Idx_EnzymesMgRXN", 0)
            Writer.Variable_("self.Idx_SubstratesMgRXN", 0)
            # Writer.Variable_("self.Idx_MgRXNsWithKineticData", 0)
            Writer.Variable_("self.Idx_MgRXNs_WithEnzymeWithoutKineticData", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Dir_RevMgRXNs", 0)

            Writer.Variable_("self.Const_Kcat_Default", Kcat_Default)
            Writer.Variable_("self.Const_Km_Default", Km_Default)
            Writer.BlankLine()

            Writer.Variable_("self.Coeff_MgRXN_Transposed", 0)
            Writer.Variable_("self.Rate_MgRXN", 0)
            Writer.BlankLine()
            
            Writer.Variable_("self.Bin_ReactantsMgRXN", 0)
            Writer.Variable_("self.Bin_ProductsMgRXN", 0)
            Writer.BlankLine()

        with Writer.Function_("SetUp_ProcessSpecificVariables"):
            Writer.Variable_("self.Idx_SubstratesMgRXN", Idx_SubstratesInMgRXN)
            Writer.BlankLine()

            Writer.Variable_("self.Const_Kcat_Default", Kcat_Default)
            Writer.Variable_("self.Const_Km_Default", Km_Default)
            Writer.BlankLine()

            Writer.Transpose("self.Coeff_MgRXN_Transposed", "self.Cel.Coeff_Magnesium")
            Writer.Cast_____("self.Coeff_MgRXN_Transposed", "self.Coeff_MgRXN_Transposed", 'float32')
            Writer.BlankLine()

            Writer.InitZeros("self.Rate_MgRXN", N_MgRXNs)
            Writer.BlankLine()

            Writer.ConvToBin("self.Bin_ReactantsMgRXN", "self.Bin_ReactantsMgRXN", ">", 0)
            Writer.ConvToBin("self.Bin_ProductsMgRXN", "self.Bin_ProductsMgRXN", ">", 0)
            Writer.BlankLine()

        with Writer.Function_("ExecuteProcess"):
            if Writer.Switch4Kinetics:
                # Calculate reaction rate for the reactions with kinetic data
                Writer.Statement("Rate_KINRXN = self.CalculateRateForReactionsWithKineticData()")
                Writer.Statement("Rate_MgRXN = self.CalculateRateForReactionsWithoutKineticData()")
                Writer.BlankLine()
                Writer.ScatNdUpd("Rate_MgRXN", "Rate_MgRXN", "self.Idx_MgRXNsWithKineticData", "Rate_KINRXN", Type='float32')
                Writer.Reshape__("Rate_MgRXN", "Rate_MgRXN", [-1, 1])
                Writer.MatrixMul("Count_MolsInMgRXN", "self.Coeff_MgRXN_Transposed", "Rate_MgRXN")
                Writer.RoundInt_("Count_MolsInMgRXN", "Count_MolsInMgRXN")
                Writer.BlankLine()
                Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_Master_MolsInMgRXN, Count_MolsInMgRXN)")
                Writer.BlankLine()

            else:
                Writer.Pass_____()
                Writer.BlankLine()

        with Writer.Function_("CalculateRateForReactionsWithKineticData"):
            Writer.Statement("Count_Enzymes = self.GetCounts_Float(self.Idx_EnzymesKINRXN)")
            Writer.Statement("Count_Substrates = self.GetCounts_Float(self.Idx_SubstratesKINRXN)")
            Writer.BlankLine()
            Writer.Statement("Rate_KINRXN = self.CalculateReactionRate(Count_Enzymes, Count_Substrates, self.Cel.Const_Kcat, self.Cel.Const_Km)")
            Writer.ReturnVar("Rate_KINRXN")
            Writer.BlankLine()

        # TODO: the reactions without kinetic data need linear algebra
        with Writer.Function_("CalculateRateForReactionsWithoutKineticData"):
            # Currently: calculate reaction rate for the reactions without kinetic data using average values
            Writer.Statement("Count_Enzymes = self.GetCounts_Float(self.Idx_EnzymesMgRXN)")
            Writer.ReduceMin("Count_EnzymesAverage", "Count_Enzymes")
            # Writer.ReduceAvg("Count_EnzymesAverage", "Count_Enzymes")

            Writer.Statement("Count_Substrates = self.GetCounts_Float(self.Idx_SubstratesMgRXN)")
            Writer.Gather___("Count_Substrates_WithEnzymes", "Count_Substrates", "self.Idx_MgRXNs_WithEnzymeWithoutKineticData")
            Writer.Reshape__("Count_Substrates_WithEnzymes", "Count_Substrates_WithEnzymes", -1)
            Writer.BlankLine()
            Writer.Statement("Rate_MgRXN_WithoutEnzymes = self.CalculateReactionRate(Count_EnzymesAverage, Count_Substrates, self.Const_Kcat_Default, self.Const_Km_Default)")
            Writer.Statement("Rate_MgRXN_WithEnzymeWithoutKineticData = self.CalculateReactionRate(Count_Enzymes, Count_Substrates_WithEnzymes, self.Const_Kcat_Default, self.Const_Km_Default)")

            Writer.ScatNdUpd("Rate_MgRXN", "Rate_MgRXN_WithoutEnzymes", "self.Idx_MgRXNs_WithEnzymeWithoutKineticData", "Rate_MgRXN_WithEnzymeWithoutKineticData", Type='float32')
            Writer.ReturnVar("Rate_MgRXN")
            Writer.BlankLine()

        with Writer.Function_("CalculateReactionRate", "Count_Enzymes", "Count_Substrates", "Const_Kcat", "Const_Km"):
            # (kcat * E * S) / (S + kM)
            Writer.Multiply_("Numerator", "Const_Kcat", "Count_Enzymes")
            Writer.Multiply_("Numerator", "Numerator", "Count_Substrates")
            Writer.Add______("Denominator", "Count_Substrates", "Const_Km")
            Writer.Divide___("Rate_RXN", "Numerator", "Denominator")
            Writer.ReturnVar("Rate_RXN")
            Writer.BlankLine()
        #
        # with Writer.Function_("ForAnalysis", "Count_MolsInMgRXN"):
        #     Writer.ZerosLike("DeltaCountFromMgRXN", "self.Cel.Counts", 'int32')
        #     Writer.ScatNdAdd("DeltaCountFromMgRXN", "DeltaCountFromMgRXN", "self.Cel.Idx_Master_MolsInMgRXN", "Count_MolsInMgRXN")
        #     for Mol in MolsToAnalyze:
        #         Writer.Gather___("self.DeltaCounts_%s" % Mol, "DeltaCountFromMgRXN", "self.Cel.Idx_%s" % Mol)
        #     Writer.Gather___("self.DeltaCounts_IDsToDebug", "DeltaCountFromMgRXN", "self.IdxsToDebug")
        #     Writer.BlankLine()
        #
        # with Writer.Function_("ViewProcessSummary"):
        #     if Writer.Switch4Kinetics:
        #         Writer.PrintStrg("===== MgRXN ===== ")
        #         for Mol in MolsToAnalyze:
        #             Writer.PrintStVa("deltaCount of %s" % Mol, "self.DeltaCounts_%s" % Mol)
        #         Writer.PrintStVa("deltaCount of IDsToDebug", "self.DeltaCounts_IDsToDebug")
        #         Writer.BlankLine()
        #     else:
        #         Writer.Pass_____()
        #         Writer.BlankLine()

# def SetUpReactions(ProGen):
#     Reactions = list()
#
#     Reaction = None
#
#     Reaction_SetUp = Reaction
#     Reactions.append(Reaction_SetUp)
#
#     return Reactions
#

