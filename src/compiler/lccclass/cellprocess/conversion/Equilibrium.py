
def Write_CellProcess(Writer, Comp, ProGen, ProcessID):

    ID_Mols = Comp.Equilibrium.ID_MolsInEQMRXN
    Idx_Mols = list()

    for ID_Mol in ID_Mols:
        Idx_Mols.append(Comp.Master.ID2Idx_Master[ID_Mol])

    assert len(ID_Mols) == len(Idx_Mols)
    NUniq_Mols = len(ID_Mols)

    Rate_EQMRXNAdjustment = 1/10000

    N_EQMRXNsTotal = len(Comp.Equilibrium.ID_EQMRXNs)

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.N_EQMRXNsToRun", 0)
            Writer.Variable_("self.Idx_MolsInEQMRXN_Local", 0)

            Writer.Variable_("self.Bool_AllMolsInEQMRXN", 0)
            Writer.Variable_("self.Bool_ReactantsInEQMRXN", 0)
            Writer.Variable_("self.Bin_ReactantsInEQMRXN", 0)
            Writer.Variable_("self.Bool_ProductsInEQMRXN", 0)
            Writer.Variable_("self.Bin_ProductsInEQMRXN", 0)
            Writer.Variable_("self.Count_ReactantsPerReaction", 0)
            Writer.Variable_("self.Count_ProductsPerReaction", 0)
            Writer.Variable_("self.N_CoeffTotalPerReactant", 0)
            Writer.Variable_("self.N_CoeffTotalPerProduct", 0)

            Writer.Variable_("self.Coeff_EQMRXN_Transposed", 0)
            Writer.BlankLine()
            
        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):

            # Master indices
            Writer.BlankLine()

            # Local indices
            Writer.VarRange_("self.Idx_MolsInEQMRXN_Local", 0, NUniq_Mols)
            Writer.BlankLine()

            # Boolean and binary matrices for coefficient matrices
            # For all reactants and products
            Writer.NotEqual_("self.Bool_AllMolsInEQMRXN", "self.Cel.Coeff_EQMRXN", 0)

            # For reactants only
            Writer.Less_____("self.Bool_ReactantsInEQMRXN", "self.Cel.Coeff_EQMRXN", 0)
            Writer.BoolToBin("self.Bin_ReactantsInEQMRXN", "self.Bool_ReactantsInEQMRXN")

            Writer.ToSparse_("self.Bool_ReactantsInEQMRXN", "self.Bool_ReactantsInEQMRXN")
            Writer.ToSparse_("self.Bin_ReactantsInEQMRXN", "self.Bin_ReactantsInEQMRXN")
            Writer.BlankLine()

            Writer.NonZeros_("self.Count_ReactantsPerReaction", "self.Bool_ReactantsInEQMRXN", 1)

            Writer.Transpose("self.Coeff_EQMRXN_Transposed", "self.Cel.Coeff_EQMRXN")
            Writer.Multiply_("Coeff_Reactants", "self.Bin_ReactantsInEQMRXN", "self.Cel.Coeff_EQMRXN")
            Writer.ReduceSum("self.N_CoeffTotalPerReactants", "Coeff_Reactants", 1)

            Writer.BlankLine()

        with Writer.Statement("def ExecuteProcess(self):"):
            Writer.Statement("self.ResetVariables()")
            Writer.Statement("Count_MolsInEQMRXN = self.GetCounts(self.Cel.Idx_Master_MolsInEQMRXN)")
            Writer.Statement("Rate_EQMRXN = DetermineRateOfEQMRXNs(Count_MolsInEQMRXN)")
            Writer.Statement("Count_MolsParticipatingInEQMRXN = self.DetermineAmountOfMoleculesParticipatingInEQMRXN(Rate_EQMRXN, Count_MolsInEQMRXN)")
            Writer.Statement("self.AddToDeltaCounts(Idx_MolsParticipatingInEQMRXN, Count_MolsParticipatingInEQMRXN)")
            Writer.BlankLine()

        with Writer.Statement("def ResetVariables(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def DetermineRateOfEQMRXNs(self, Count_MolsInEQMRXN):"):
            # Rate is needed only if ODE is used

            # Kc must be known in the first place for the given condition

            # Reactants by boolean then multiply axis
            # Counts ** -Coeff

            # Products by boolean then multiply axis
            # Counts ** Coeff (0 power becomes 1)

            # Qc = Product / Reactants

            # If Qc is > Kc, then the reaction occurs
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def DetermineAmountOfMoleculesParticipatingInEQMRXN(self, Rate_EQMRXN, Count_MolsInEQMRXN):"):
            Writer.MatrixMul("Count_MolsParticipatingInEQMRXN", "self.Coeff_EQMRXN", "Rate_EQMRXN")

            # Get the number of reactions to run


            # Adjust below zeros

            Writer.ReturnVar("Count_MolsParticipatingInEQMRXN_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.PrintStrg("===== Equilibrium Reactions ===== ")
            Writer.PrintStVa("",
                             "")
            Writer.BlankLine()
