# DNA polymerase
# RNA polymerase
# Ribosome assembly and disassembly
# check for sigma factor regulation in place?


def Write_CellProcess(Writer, Comp, ProGen, ProcessID):

    # Molecule indices for Molecular IDs
    # Idx_ClpP_SerineProtease = Comp.Master.ID2Idx_Master['EG10158-MONOMER']

    ID_Mols = Comp.Complexation.ID_MolsInComplexations
    Idx_Mols = list()

    for ID_Mol in ID_Mols:
        Idx_Mols.append(Comp.Master.ID2Idx_Master[ID_Mol])

    assert len(ID_Mols) == len(Idx_Mols)
    NUniq_Mols = len(ID_Mols)

    N_CPLXRXNsTotal = len(Comp.Complexation.ID_Complexations)

    # Arbitrarily set number of complexation reactions to run each simulation step
    N_CPLXRXNsToRun = 2

    # Resolution for weight in selecting a complexation reaction
    WeightResolution = 2 ** 30

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.N_CPLXRXNsToRun", 0)
            Writer.Variable_("self.Idx_MolsInComplexation_Local", 0)
            Writer.Variable_("self.Coeff_Complexation_Transposed", 0)
            Writer.BlankLine()

        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):

            # Master indices
            Writer.BlankLine()

            # Local indices
            Writer.Variable_("self.N_CPLXRXNsToRun", N_CPLXRXNsToRun)
            Writer.VarRange_("self.Idx_MolsInComplexation_Local", 0, NUniq_Mols)
            Writer.BlankLine()

            # Boolean and binary matrices for coefficient matrices
            # For all reactants and products
            Writer.NotEqual_("self.Bool_MolsInComplexation", "self.Cel.Coeff_Complexation", 0)

            # For reactants only
            Writer.Less_____("self.Bool_ReactantsInComplexation", "self.Cel.Coeff_Complexation", 0)
            Writer.BoolToBin("self.Bin_ReactantsInComplexation", "self.Bool_ReactantsInComplexation")

            Writer.ToSparse_("self.Bool_ReactantsInComplexation", "self.Bool_ReactantsInComplexation")
            Writer.ToSparse_("self.Bin_ReactantsInComplexation", "self.Bin_ReactantsInComplexation")
            Writer.BlankLine()

            Writer.NonZeros_("self.Count_ReactantsPerReaction", "self.Bool_ReactantsInComplexation", 1)

            Writer.Transpose("self.Coeff_Complexation_Transposed", "self.Cel.Coeff_Complexation")
            Writer.Multiply_("Coeff_Reactants", "self.Bin_ReactantsInComplexation", "self.Cel.Coeff_Complexation")
            Writer.ReduceSum("self.N_CoeffTotalPerReactants", "Coeff_Reactants", 1)

            Writer.BlankLine()

        with Writer.Statement("def ExecuteProcess(self):"):
            Writer.Statement("self.ResetVariables()")
            Writer.Statement("Count_MolsInComplexation = self.GetCounts(self.Cel.Idx_Master_MolsInComplexation)")
            Writer.Statement("Idx_CPLXRXNAvailable, Weight_CPLXRXNAvailable = self.DetermineAvailableCPLXRXN(Count_MolsInComplexation)")
            Writer.Statement("Idx_RndCPLXRXNToRun = self.PickRandomIndexFromPool_Weighted_Global(self.N_CPLXRXNsToRun, Idx_CPLXRXNAvailable, Weight_CPLXRXNAvailable)")
            Writer.Statement("Idx_MolsParticipatingInComplexation, Count_MolsParticipatingInComplexation = self.DetermineAmountOfMoleculesParticipatingInComplexation(Idx_RndCPLXRXNToRun)")
            Writer.Statement("self.AddToDeltaCounts(Idx_MolsParticipatingInComplexation, Count_MolsParticipatingInComplexation)")
            Writer.BlankLine()

        with Writer.Statement("def ResetVariables(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def DetermineAvailableCPLXRXN(self, Count_MolsInComplexation):"):
            Writer.Transpose("Count_MolsInComplexation", "Count_MolsInComplexation")
            Writer.Statement("Bool_ReactionsAvailable = self.DetermineIdxOfAvailableCPLXRXN(Count_MolsInComplexation)")
            Writer.GenIdx___("Idx_CPLXRXNAvailable", "Bool_ReactionsAvailable")
            Writer.Statement("Weight_CPLXRXN = self.DetermineWeightOfAvailableCPLXRXN(Count_MolsInComplexation)")
            Writer.BoolMask_("Weight_CPLXRXNAvailable", "Weight_CPLXRXN", "Bool_ReactionsAvailable")
            Writer.ReturnVar("Idx_CPLXRXNAvailable", "Weight_CPLXRXNAvailable")
            Writer.BlankLine()

        with Writer.Statement("def DetermineIdxOfAvailableCPLXRXN(self, Count_MolsInComplexation):"):
            Writer.Comment__("Determine reactions with sufficient reactant count")
            Writer.Greater__("Bool_CountsSufficient", "Count_MolsInComplexation", "-self.Cel.Coeff_Complexation")
            Writer.LogicAnd_("Bool_ReactantsWithSufficientCount", "self.Bool_ReactantsInComplexation", "Bool_CountsSufficient")
            Writer.NonZeros_("Count_ReactantsSufficientPerReaction", "Bool_ReactantsWithSufficientCount", 1)
            Writer.Equal____("Bool_ReactionsAvailable", "Count_ReactantsSufficientPerReaction", "self.Count_ReactantsPerReaction")
            Writer.ReturnVar("Bool_ReactionsAvailable")
            Writer.BlankLine()

        with Writer.Statement("def DetermineWeightOfAvailableCPLXRXN(self, Count_MolsInComplexation):"):
            Writer.Comment__("Determine counts of reactions with sufficient reactants")
            Writer.Negative_("Coeff_Complexation_ReactantsWithPositiveValues", "self.Cel.Coeff_Complexation")
            Writer.Divide___("Count_ReactionsPossiblePerMolecule", "Count_MolsInComplexation", "Coeff_Complexation_ReactantsWithPositiveValues")
            Writer.RoundInt_("Count_ReactionsPossiblePerMolecule_Rounded", "Count_ReactionsPossiblePerMolecule")
            Writer.Multiply_("Count_ReactionsPossiblePerReactant", "self.Bin_ReactantsInComplexation", "Count_ReactionsPossiblePerMolecule_Rounded")

            Writer.Comment__("Generate relative weight of all complexation reactions in Int32")
            Writer.Statement("Weight_ReactionsPossiblePerReactant_Squeezed = self.SqueezeDistributionRangeZeroAndOne(Count_ReactionsPossiblePerReactant)")
            Writer.Statement("Weight_ReactionsPossiblePerReactant_SqueezedWithZerosToOnes = self.ResetZerosToSpecificValue_Float(Weight_ReactionsPossiblePerReactant_Squeezed, 1)")
            Writer.ReduceMul("Weight_CPLXRXNsPossible", "Weight_ReactionsPossiblePerReactant_SqueezedWithZerosToOnes", 1)
            Writer.Multiply_("Weight_CPLXRXNsPossible_Stretched", "Weight_CPLXRXNsPossible", WeightResolution)
            # Writer.Statement("Weight_CPLXRXNsPossible_Stretched = self.StretchDistributionRangeToNegAndPosValue(Weight_CPLXRXNsPossible, %s)" % WeightResolution)
            Writer.RoundInt_("Weight_CPLXRXNsPossible_StretchedInt", "Weight_CPLXRXNsPossible_Stretched")
            Writer.Comment__("Give a minimum weight to the reactions rounded to 0 by adding 1 to all")
            Writer.Add______("Weight_CPLXRXNsPossible_Adjusted", "Weight_CPLXRXNsPossible_StretchedInt", 1)
            Writer.ReturnVar("Weight_CPLXRXNsPossible_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def DetermineAmountOfMoleculesParticipatingInComplexation(self, Idx_CPLXRXNs):"):
            Writer.Statement("Count_Idx_CPLXRXNs = self.GenerateCountMatrixForSelectedIndices(Idx_CPLXRXNs, %s)" % N_CPLXRXNsTotal)
            Writer.Reshape__("Count_Idx_CPLXRXNs", "Count_Idx_CPLXRXNs", [-1, 1])
            Writer.MatrixMul("Count_MolsParticipatingInComplexation", "self.Coeff_Complexation_Transposed", "Count_Idx_CPLXRXNs")
            Writer.Statement("Count_MolsParticipatingInComplexation_ZerosRemoved = self.RemoveZeroElements(Count_MolsParticipatingInComplexation)")
            Writer.GenIdx___("Idx_MolsParticipatingInComplexation_Local", "Count_MolsParticipatingInComplexation")
            Writer.Statement("Idx_MolsParticipatingInComplexation = self.IdxFromLocalToReference(Idx_MolsParticipatingInComplexation_Local, self.Cel.Idx_Master_MolsInComplexation)")
            Writer.ReturnVar("Idx_MolsParticipatingInComplexation", "Count_MolsParticipatingInComplexation_ZerosRemoved")
            Writer.BlankLine()

        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.PrintStrg("===== Complexation ===== ")
            Writer.PrintStVa("# of Complexation Reactions Occurred",
                             "self.N_CPLXRXNsToRun[0]")
            Writer.BlankLine()
