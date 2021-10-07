# TCA Cycle

import numpy as np

'''
'''

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):
    """

    :type Comp: object
    """
    # Parse flat data

    Idx_RXNsInTCACycle = list()
    Idx_SubstratesInTCACycle4Count = list()
    Idx_EnzymesInTCACycle4Count = list()

    ID_SubstratesInTCACycle = list()
    ID_EnzymesInTCACycle = list()

    Idx_RXNsInTCACycle_WithEnzymeWithoutKineticData = list()   # Partial list with kinetic data
    Idx_EnzymesInTCACycle4Kinetics = list()

    for ID_RXN in Comp.Metabolism.ID_RXNsInTCACycle:
        Idx = Comp.Metabolism.ID2Idx_METRXNs[ID_RXN]
        Idx_RXNsInTCACycle.append(Idx)

        SubstrateID = Comp.Metabolism.ID_Substrates4METRXN[Idx]
        ID_SubstratesInTCACycle.append(SubstrateID)
        Idx_SubstratesInTCACycle4Count.append(Comp.Master.ID2Idx_Master[SubstrateID])

        EnzymeID = Comp.Metabolism.ID_Enzymes4METRXN[Idx]
        ID_EnzymesInTCACycle.append(EnzymeID)
        Idx_EnzymesInTCACycle4Count.append(Comp.Master.ID2Idx_Master[EnzymeID])

        if ID_RXN in Comp.Kinetics.ID2Idx_RXNID2KINRXN:
            Idx_RXNsInTCACycle_WithEnzymeWithoutKineticData.append(Comp.Metabolism.ID2Idx_METRXNs[ID_RXN])
            Idx_EnzymesInTCACycle4Kinetics.append(Comp.Kinetics.ID2Idx_RXNID2KINRXN[ID_RXN])
        else:
            print('No kinetic data available for TCA Cycle ReactionID: %s' % ID_RXN)

    NUniq_RXNsInTCACycle = Comp.Metabolism.NUniq_RXNsInTCACycle

    ID_TCACycleMetabolitesToReplenish = [
        'ACETYL-COA[c]',
        'NAD[c]',
        'FAD[c]',
        'WATER[c]',
        'GDP[c]',
        'ADP[c]',
    ]

    Idx_TCACycleMetabolitesToReplenish = ProGen.GetMolIdx_Master(ID_TCACycleMetabolitesToReplenish)

    ID_SubstratesInTCACycle4Display = Comp.Metabolism.ID_Substrates4TCARXN
    Idx_SubstratesInTCACycle4Display = ProGen.GetMolIdx_Master(ID_SubstratesInTCACycle4Display)

    ID_EnergyProductsInTCACycle = [
        'CO-A[c]',
        'NADH[c]',
        'ATP[c]',
        'GTP[c]',
    ]

    Idx_EnergyProductsInTCACycle = ProGen.GetMolIdx_Master(ID_EnergyProductsInTCACycle)

    # FROM METABOLISM

    # For the reactions with kinetic data
    Idx_SubstratesKINRXN = ProGen.GetMolIdx_Master(Comp.Kinetics.ID_Substrates4KINRXN)
    Idx_EnzymesKINRXN = ProGen.GetMolIdx_Master(Comp.Kinetics.ID_Enzymes4KINRXN)

    # For the reactions without kinetic data
    Idx_SubstratesMETRXN = ProGen.GetMolIdx_Master(Comp.Metabolism.ID_Substrates4METRXN)
    Idx_EnzymesMETRXN = ProGen.GetMolIdx_Master(Comp.Metabolism.ID_Enzymes4METRXN)   # Not for all METRXNs

    # Enzymes available for a subset of metabolism reactions
    Idx_METRXNs_WithEnzymeWithoutKineticData = list()
    for EnzymeID in Comp.Metabolism.ID_Enzymes4METRXN:
        Idx_METRXNs_WithEnzymeWithoutKineticData.append(Comp.Metabolism.ID2Idx_Enzymes4METRXN[EnzymeID])

    Idx_METRXNsWithKineticData = list()
    for ID_KINRXN in Comp.Kinetics.ID_KINRXNs:
        Idx_METRXNsWithKineticData.append(Comp.Metabolism.ID2Idx_METRXNs[ID_KINRXN])

    # Kcat and Km default values for reactions without kinetic data have been set in Dataset.py
    # This serves as a local control in cell.py
    Kcat_Default = Comp.Kinetics.Const_Kcat_Default
    Km_Default = Comp.Kinetics.Const_Km_Default


    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Function_("Init_ProcessSpecificVariables"):
            Writer.Variable_("self.Idx_TCACycleMetabolitesToReplenish", 0)
            Writer.Variable_("self.Count_TCACycleMetabolitesToReplenish", 0)

            Writer.Variable_("self.Idx_RXNsInTCACycle", 0)
            Writer.Variable_("self.Idx_SubstratesInTCACycle4Count", 0)
            Writer.Variable_("self.Idx_EnzymesInTCACycle4Count", 0)
            Writer.Variable_("self.Idx_RXNsInTCACycle_WithEnzymeWithoutKineticData", 0)
            Writer.Variable_("self.Idx_EnzymesInTCACycle4Kinetics", 0)
            Writer.Variable_("self.Idx_EnergyProductsInTCACycle", 0)
            Writer.BlankLine()

            Writer.Variable_("self.Const_Kcat_Default", Kcat_Default)
            Writer.Variable_("self.Const_Km_Default", Km_Default)
            Writer.BlankLine()

            Writer.Variable_("self.Coeff_TCACycle_Transposed", 0)
            Writer.BlankLine()

            # FROM METABOLISM
            Writer.Variable_("self.Idx_Enzymes", 0)
            Writer.Variable_("self.Idx_Substrates", 0)
            Writer.Variable_("self.Idx_EnzymesMETRXN", 0)
            Writer.Variable_("self.Idx_SubstratesMETRXN", 0)
            Writer.Variable_("self.Idx_METRXNsWithKineticData", 0)
            Writer.Variable_("self.Idx_METRXNs_WithEnzymeWithoutKineticData", 0)
            Writer.BlankLine()

        with Writer.Function_("SetUp_ProcessSpecificVariables"):
            Writer.Variable_("self.Idx_TCACycleMetabolitesToReplenish", Idx_TCACycleMetabolitesToReplenish)

            Writer.Variable_("self.Idx_RXNsInTCACycle", Idx_RXNsInTCACycle)
            Writer.Variable_("self.Idx_SubstratesInTCACycle4Count", Idx_SubstratesInTCACycle4Count)
            Writer.Variable_("self.Idx_EnzymesInTCACycle4Count", Idx_EnzymesInTCACycle4Count)
            Writer.Variable_("self.Idx_RXNsInTCACycle_WithEnzymeWithoutKineticData", Idx_RXNsInTCACycle_WithEnzymeWithoutKineticData)
            Writer.Variable_("self.Idx_EnzymesInTCACycle4Kinetics", Idx_EnzymesInTCACycle4Kinetics)
            Writer.Variable_("self.Idx_EnergyProductsInTCACycle", Idx_EnergyProductsInTCACycle)
            Writer.BlankLine()

            # Transform metabolism coeff matrix to TCA cycle only with the rest zeros
            Writer.Variable_("NUniq_MolsInMETRXN", Comp.Metabolism.NUniq_MolsInMETRXN)
            Writer.Statement("Coeff_TCACycle = tf.reshape(tf.gather(self.Cel.Coeff_Metabolism, self.Idx_RXNsInTCACycle), [-1, NUniq_MolsInMETRXN[0]])")
            Writer.ZerosLike("Zeros_Coeff_Metabolism", "self.Cel.Coeff_Metabolism", 'int32')
            Writer.Statement("Coeff_TCACycleWithZeros = tf.tensor_scatter_nd_update(Zeros_Coeff_Metabolism, tf.reshape(self.Idx_RXNsInTCACycle, [-1, 1]), Coeff_TCACycle)")

            # # Get indices of molecules involved in TCA cycle for count update (may not be necessary)
            # Writer.ConvToBin("Bin_MolsInTCACycle", "Coeff_TCACycle", "!=", 0)
            # Writer.ReduceSum("Bin_MolsInTCACyclePerMol", "Bin_MolsInTCACycle", 1)
            # Writer.Reshape__("Bin_MolsInTCACyclePerMol", "Bin_MolsInTCACyclePerMol", -1)
            # Writer.GenIdx___("Idx_MolsInTCACycle", "Bin_MolsInTCACyclePerMol")

            Writer.Comment__("Replacing the metabolism coefficient matrix with only TCA reactions with values")
            Writer.Transpose("self.Coeff_Metabolism_Transposed", "Coeff_TCACycleWithZeros")
            Writer.Cast_____("self.Coeff_Metabolism_Transposed", "self.Coeff_Metabolism_Transposed", 'float32')
            Writer.BlankLine()

            # FROM METABOLISM
            Writer.Variable_("self.Idx_EnzymesKINRXN", Idx_EnzymesKINRXN)
            Writer.Variable_("self.Idx_SubstratesKINRXN", Idx_SubstratesKINRXN)
            Writer.BlankLine()
            Writer.Variable_("self.Idx_EnzymesMETRXN", Idx_EnzymesMETRXN)
            Writer.Variable_("self.Idx_SubstratesMETRXN", Idx_SubstratesMETRXN)
            Writer.Variable_("self.Idx_METRXNsWithKineticData", Idx_METRXNsWithKineticData)
            Writer.Variable_("self.Idx_METRXNs_WithEnzymeWithoutKineticData", Idx_METRXNs_WithEnzymeWithoutKineticData)

            Writer.BlankLine()
            Writer.Transpose("self.Coeff_Metabolism_Transposed", "self.Cel.Coeff_Metabolism")
            Writer.Cast_____("self.Coeff_Metabolism_Transposed", "self.Coeff_Metabolism_Transposed", 'float32')
            Writer.BlankLine()

            Writer.Comment__("For Replenishing Metabolites")
            Writer.Statement("self.GetMetaboliteCounts()")
            Writer.BlankLine()


        with Writer.Function_("GetMetaboliteCounts"):
            Writer.Gather___("self.Count_TCACycleMetabolitesToReplenish", "self.Cel.Count_Master", "self.Idx_TCACycleMetabolitesToReplenish")
            Writer.BlankLine()

        with Writer.Function_("ReplenishMetabolites"):
            Writer.ScatNdUpd("self.Cel.Counts", "self.Cel.Counts", "self.Idx_TCACycleMetabolitesToReplenish", "self.Count_TCACycleMetabolitesToReplenish")
            Writer.BlankLine()

        with Writer.Function_("Message_ReplenishMetabolites"):
            Writer.PrintStrg("\t- TCA Cycle Input Metabolites have been replenished.")
            Writer.BlankLine()


        with Writer.Function_("ExecuteProcess"):
            if Writer.Switch4Kinetics:
                # Calculate reaction rate for the reactions with kinetic data
                Writer.Statement("Rate_KINRXN = self.CalculateRateForReactionsWithKineticData()")
                Writer.Statement("Rate_METRXN = self.CalculateRateForReactionsWithoutKineticData()")
                Writer.BlankLine()
                Writer.ScatNdUpd("Rate_METRXN", "Rate_METRXN", "self.Idx_METRXNsWithKineticData", "Rate_KINRXN", Type='float32')
                Writer.Reshape__("Rate_METRXN", "Rate_METRXN", [-1, 1])
                Writer.MatrixMul("Count_MolsInMETRXN", "self.Coeff_Metabolism_Transposed", "Rate_METRXN")
                Writer.RoundInt_("Count_MolsInMETRXN", "Count_MolsInMETRXN")
                Writer.BlankLine()
                Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_Master_MolsInMETRXN, Count_MolsInMETRXN)")
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
            Writer.Statement("Count_Enzymes = self.GetCounts_Float(self.Idx_EnzymesMETRXN)")
            Writer.ReduceMin("Count_EnzymesAverage", "Count_Enzymes")
            # Writer.ReduceAvg("Count_EnzymesAverage", "Count_Enzymes")

            Writer.Statement("Count_Substrates = self.GetCounts_Float(self.Idx_SubstratesMETRXN)")
            Writer.Gather___("Count_Substrates_WithEnzymes", "Count_Substrates", "self.Idx_METRXNs_WithEnzymeWithoutKineticData")
            Writer.Reshape__("Count_Substrates_WithEnzymes", "Count_Substrates_WithEnzymes", -1)
            Writer.BlankLine()
            Writer.Statement("Rate_METRXN_WithoutEnzymes = self.CalculateReactionRate(Count_EnzymesAverage, Count_Substrates, self.Const_Kcat_Default, self.Const_Km_Default)")
            Writer.Statement("Rate_METRXN_WithEnzymeWithoutKineticData = self.CalculateReactionRate(Count_Enzymes, Count_Substrates_WithEnzymes, self.Const_Kcat_Default, self.Const_Km_Default)")

            Writer.ScatNdUpd("Rate_METRXN", "Rate_METRXN_WithoutEnzymes", "self.Idx_METRXNs_WithEnzymeWithoutKineticData", "Rate_METRXN_WithEnzymeWithoutKineticData", Type='float32')
            Writer.ReturnVar("Rate_METRXN")
            Writer.BlankLine()

        with Writer.Function_("CalculateReactionRate", "Count_Enzymes", "Count_Substrates", "Const_Kcat", "Const_Km"):
            # (kcat * E * S) / (S + kM)
            Writer.Multiply_("Numerator", "Const_Kcat", "Count_Enzymes")
            Writer.Multiply_("Numerator", "Numerator", "Count_Substrates")
            Writer.Add______("Denominator", "Count_Substrates", "Const_Km")
            Writer.Divide___("Rate_RXN", "Numerator", "Denominator")
            Writer.ReturnVar("Rate_RXN")
            Writer.BlankLine()

        with Writer.Function_("ViewProcessSummary"):
            if Writer.Switch4Kinetics:
                Writer.PrintStrg("===== TCA Cycle ===== ")
                Writer.PrintStrg("Substrates:")
                Writer.PrintStrg(str(ID_SubstratesInTCACycle[:5]))
                Writer.PrintVar_("tf.gather(self.Cel.Counts, self.Idx_SubstratesInTCACycle4Count)[0][:5]")
                Writer.PrintStrg(str(ID_SubstratesInTCACycle[5:]))
                Writer.PrintVar_("tf.gather(self.Cel.Counts, self.Idx_SubstratesInTCACycle4Count)[0][5:]")
                Writer.BlankLine()
                Writer.PrintStrg("Energy Products:")
                Writer.PrintStrg(str(ID_EnergyProductsInTCACycle))
                Writer.PrintVar_("tf.gather(self.Cel.Counts, self.Idx_EnergyProductsInTCACycle)")
                Writer.BlankLine()
            else:
                Writer.Pass_____()
                Writer.BlankLine()
