# Transport

import numpy as np

'''
rate_sugar transporter = ~100/s  (typically kcat 30~300/s)   if protein-coupled
 
'''

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):
    """

    :type Comp: object
    """

    # Indices
    Idx_Substrates4TRNRXN = ProGen.GetMolIdx_Master(Comp.Transporters.ID_MolsInTRNRXN)
    Idx_Enzymes4TRNRXN = ProGen.GetMolIdx_Master(Comp.Transporters.ID_Enzymes4TRNRXN)

    # Coeff_TRNRXNs = Comp.Transporters.Coeff_MolsInTRNRXN
    N_TRNRXNs = len(Comp.Transporters.Coeff_MolsInTRNRXN)

    Rate_Enzymes4TRNRXN = list()

    Type2Rate_Transporter = {
        'Channel'   : 1,
        'Carrier'   : 0,
    }

    for Transporter in Comp.Transporters.ID_Enzymes4TRNRXN:
        Rate_Enzymes4TRNRXN.append(Type2Rate_Transporter[Comp.Transporters.ID2Type_Enzymes4TRNRXN[Transporter]])

    print(Rate_Enzymes4TRNRXN)

    Kcat_Default = 0.1
    Km_Default = 0.0001

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Function_("Init_ProcessSpecificVariables"):
            Writer.Variable_("self.Count_MetabolitesInitial", 0)

            Writer.Variable_("self.Idx_Enzymes", 0)
            Writer.Variable_("self.Idx_Substrates", 0)
            Writer.Variable_("self.Idx_EnzymesMETRXN", 0)
            Writer.Variable_("self.Idx_SubstratesMETRXN", 0)
            Writer.Variable_("self.Idx_METRXNsWithKineticData", 0)
            Writer.Variable_("self.Idx_METRXNs_WithEnzymeWithoutKineticData", 0)
            Writer.BlankLine()

            Writer.Variable_("self.Const_Kcat_Default", Kcat_Default)
            Writer.Variable_("self.Const_Km_Default", Km_Default)
            Writer.BlankLine()

            Writer.Variable_("self.Coeff_Metabolism_Transposed", 0)
            Writer.BlankLine()

            Writer.Variable_("self.IdxsToDebug", 0)
            Writer.BlankLine()

        with Writer.Function_("SetUp_ProcessSpecificVariables"):

            Writer.Variable_("self.Idx_Enzymes4TRNRXN", Idx_Enzymes4TRNRXN)
            Writer.Variable_("self.Idx_Substrates4TRNRXN", Idx_Substrates4TRNRXN)
            Writer.BlankLine()
            # Writer.Variable_("self.Idx_EnzymesMETRXN", Idx_EnzymesMETRXN)
            # Writer.Variable_("self.Idx_SubstratesMETRXN", Idx_SubstratesMETRXN)
            # Writer.Variable_("self.Idx_METRXNsWithKineticData", Idx_METRXNsWithKineticData)
            # Writer.Variable_("self.Idx_METRXNs_WithEnzymeWithoutKineticData", Idx_METRXNs_WithEnzymeWithoutKineticData)
            # Writer.BlankLine()
            Writer.Transpose("self.Coeff_Metabolism_Transposed", "self.Cel.Coeff_Metabolism")
            Writer.Cast_____("self.Coeff_Metabolism_Transposed", "self.Coeff_Metabolism_Transposed", 'float32')
            Writer.BlankLine()

            # Writer.Comment__("For Analysis Only")
            # Writer.Variable_("self.IdxsToDebug", IdxsToDebug)
            # Writer.BlankLine()

            # Writer.Comment__("For Replenishing Metabolites")
            # Writer.Statement("self.GetMetaboliteCounts()")
            # Writer.BlankLine()


        # with Writer.Function_("GetMetaboliteCounts"):
        #     Writer.Gather___("self.Count_MetabolitesInitial", "self.Cel.Count_Master", "self.Cel.Idx_Master_Metabolites")
        #     Writer.BlankLine()
        # 
        # with Writer.Function_("ReplenishMetabolites"):
        #     Writer.ScatNdUpd("self.Cel.Counts", "self.Cel.Counts", "self.Cel.Idx_Master_Metabolites", "self.Count_MetabolitesInitial")
        #     Writer.BlankLine()
        # 
        # with Writer.Function_("Message_ReplenishMetabolites"):
        #     Writer.PrintStrg("\t- Metabolites have been replenished.")
        #     Writer.BlankLine()
        # 

        with Writer.Function_("ExecuteProcess"):
            if Writer.Switch4Kinetics:
                # Calculate reaction rate for the reactions with kinetic data
                Writer.Statement("Rate_TRNRXN = self.CalculateRateForReactionsWithKineticData()")
                Writer.Statement("Rate_TRNRXN = self.CalculateRateForReactionsWithoutKineticData()")
                Writer.BlankLine()
                Writer.ScatNdUpd("Rate_TRNRXN", "Rate_TRNRXN", "self.Idx_TRNRXNsWithKineticData", "Rate_TRNRXN", Type='float32')
                Writer.Reshape__("Rate_TRNRXN", "Rate_TRNRXN", [-1, 1])
                Writer.MatrixMul("Count_MolsInTRNRXN", "self.Coeff_Metabolism_Transposed", "Rate_TRNRXN")
                Writer.RoundInt_("Count_MolsInTRNRXN", "Count_MolsInTRNRXN")
                Writer.BlankLine()
                Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_Master_MolsInTRNRXN, Count_MolsInTRNRXN)")
                Writer.BlankLine()
                Writer.Comment__("For Analysis Only")
                Writer.Statement("self.ForAnalysis(Count_MolsInTRNRXN)")
                Writer.BlankLine()

            else:
                Writer.Pass_____()
                Writer.BlankLine()

        with Writer.Function_("CalculateRateForReactionsWithKineticData"):
            Writer.Statement("Count_Enzymes = self.GetCounts_Float(self.Idx_EnzymesTRNRXN)")
            Writer.Statement("Count_Substrates = self.GetCounts_Float(self.Idx_SubstratesTRNRXN)")
            Writer.BlankLine()
            Writer.Statement("Rate_TRNRXN = self.CalculateReactionRate(Count_Enzymes, Count_Substrates, self.Cel.Const_Kcat, self.Cel.Const_Km)")
            Writer.ReturnVar("Rate_TRNRXN")
            Writer.BlankLine()

        # TODO: the reactions without kinetic data need linear algebra
        with Writer.Function_("CalculateRateForReactionsWithoutKineticData"):
            # Currently: calculate reaction rate for the reactions without kinetic data using average values
            Writer.Statement("Count_Enzymes = self.GetCounts_Float(self.Idx_EnzymesTRNRXN)")
            Writer.ReduceMin("Count_EnzymesAverage", "Count_Enzymes")
            # Writer.ReduceAvg("Count_EnzymesAverage", "Count_Enzymes")

            Writer.Statement("Count_Substrates = self.GetCounts_Float(self.Idx_SubstratesTRNRXN)")
            Writer.Gather___("Count_Substrates_WithEnzymes", "Count_Substrates", "self.Idx_TRNRXNs_WithEnzymeWithoutKineticData")
            Writer.Reshape__("Count_Substrates_WithEnzymes", "Count_Substrates_WithEnzymes", -1)
            Writer.BlankLine()
            Writer.Statement("Rate_TRNRXN_WithoutEnzymes = self.CalculateReactionRate(Count_EnzymesAverage, Count_Substrates, self.Const_Kcat_Default, self.Const_Km_Default)")
            Writer.Statement("Rate_TRNRXN_WithEnzymeWithoutKineticData = self.CalculateReactionRate(Count_Enzymes, Count_Substrates_WithEnzymes, self.Const_Kcat_Default, self.Const_Km_Default)")

            Writer.ScatNdUpd("Rate_TRNRXN", "Rate_TRNRXN_WithoutEnzymes", "self.Idx_TRNRXNs_WithEnzymeWithoutKineticData", "Rate_TRNRXN_WithEnzymeWithoutKineticData", Type='float32')
            Writer.ReturnVar("Rate_TRNRXN")
            Writer.BlankLine()

        with Writer.Function_("CalculateReactionRate", "Count_Enzymes", "Count_Substrates", "Const_Kcat", "Const_Km"):
            # (kcat * E * S) / (S + kM)
            Writer.Multiply_("Numerator", "Const_Kcat", "Count_Enzymes")
            Writer.Multiply_("Numerator", "Numerator", "Count_Substrates")
            Writer.Add______("Denominator", "Count_Substrates", "Const_Km")
            Writer.Divide___("Rate_RXN", "Numerator", "Denominator")
            Writer.ReturnVar("Rate_RXN")
            Writer.BlankLine()

        with Writer.Function_("ForAnalysis", "Count_MolsInTRNRXN"):
            Writer.ZerosLike("DeltaCountFromMetabolism", "self.Cel.Counts", 'int32')
            Writer.ScatNdAdd("DeltaCountFromMetabolism", "DeltaCountFromMetabolism", "self.Cel.Idx_Master_MolsInTRNRXN", "Count_MolsInTRNRXN")
            for Mol in MolsToAnalyze:
                Writer.Gather___("self.DeltaCounts_%s" % Mol, "DeltaCountFromMetabolism", "self.Cel.Idx_%s" % Mol)
            Writer.Gather___("self.DeltaCounts_IDsToDebug", "DeltaCountFromMetabolism", "self.IdxsToDebug")
            Writer.BlankLine()

        with Writer.Function_("ViewProcessSummary"):
            if Writer.Switch4Kinetics:
                Writer.PrintStrg("===== Metabolism ===== ")
                for Mol in MolsToAnalyze:
                    Writer.PrintStVa("deltaCount of %s" % Mol, "self.DeltaCounts_%s" % Mol)
                Writer.PrintStVa("deltaCount of IDsToDebug", "self.DeltaCounts_IDsToDebug")
                Writer.BlankLine()
            else:
                Writer.Pass_____()
                Writer.BlankLine()

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

