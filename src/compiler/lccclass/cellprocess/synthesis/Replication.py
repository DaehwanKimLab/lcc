# dna replication: DNA Polymerase III holozyme,  DNA Polymerase II (back-up)
'''
https://en.wikipedia.org/wiki/DNA_polymerase_III_holoenzyme
The replisome is composed of the following:

2 DNA Pol III enzymes, each comprising α, ε and θ subunits. (It has been proven that there is a third copy of Pol III at the replisome.[1])
the α subunit (encoded by the dnaE gene) has the polymerase activity.
the ε subunit (dnaQ) has 3'→5' exonuclease activity.
the θ subunit (holE) stimulates the ε subunit's proofreading.
2 β units (dnaN) which act as sliding DNA clamps, they keep the polymerase bound to the DNA.
2 τ units (dnaX) which act to dimerize two of the core enzymes (α, ε, and θ subunits).
1 γ unit (also dnaX) which acts as a clamp loader for the lagging strand Okazaki fragments, helping the two β subunits to form a unit and bind to DNA. The γ unit is made up of 5 γ subunits which include 3 γ subunits, 1 δ subunit (holA), and 1 δ' subunit (holB). The δ is involved in copying of the lagging strand.
Χ (holC) and Ψ (holD) which form a 1:1 complex and bind to γ or τ. X can also mediate the switch from RNA primer to DNA.[2]

1000 nucleotides per second
'''

'''
When Chr bp + 1,
+ 2 Pi
- 1 (use ACGT (frequency table for stoichiometry))

rate: + 1000 bp / 1 second
'''

import Utils

# Comp is a short hand for CompilerData
def Write_SynDNA(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FSynDNA(FCellProcess):"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.Idx_RXN", 0)

            Writer.Variable_("self.Idx_Mols", 0)
            Writer.Variable_("self.Stoich", 0)

            Writer.Variable_("self.Rate", 0)

            Writer.BlankLine()
            Writer.Statement("super().__init__()")

            Writer.BlankLine()


        with Writer.Statement("def InitProcess(self):"):
            # Enter elementary reaction equation and rate

            Writer.Statement("self.DefineElementaryRXN()")
            Writer.Statement("self.AddElementaryRXN()")

            Writer.BlankLine()

        Writer.TF_Graph_()
        with Writer.Statement("def LoopProcess(self):"):

            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def DefineElementaryRXN(self):"):
            # Information from L++ code parsed by the compiler:

            RXNEquation = '2 dNTP -> 1 ChromosomeSize + 2 PPi'
            RXNRate = '1000 events per second'

            Writer.Comment__(RXNEquation)
            Writer.Comment__(RXNRate)

            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def AddElementaryRXN(self):"):
            # Use regular expression to parse the reaction stoichiometry info
            # The final product would be the following:

            ID_Consumed = ['dNTP']
            ID_Produced = ['Chromosome1_Rep2', 'PPI[c]']
            Stoich_Consumed = [2]
            Stoich_Produced = [1, 2]
            Rate = 1000.0

            # Parse the molecule inputs for indexing
            ID_Consumed, Stoich_Consumed = Utils.RefineBuildingBlocks(ID_Consumed, Stoich_Consumed, Comp)
            ID_Produced, Stoich_Produced = Utils.RefineBuildingBlocks(ID_Produced, Stoich_Produced, Comp)

            # Correct Stoich_Consumed to negative values
            Stoich_Consumed = Utils.ReverseNumberSign(Stoich_Consumed)

            # Concatenate ID and Stoich of the Consumed and Produced
            ID_All = ID_Consumed + ID_Produced
            Stoich_All = Stoich_Consumed + Stoich_Produced

            # Find Indexes for Molecules
            Idx_All = Utils.GetMolIdx(ID_All, Comp.Master.ID2Idx_Master)

            # Convert index and stoichiometry to tensor for simulation
            Writer.Variable_("Idx", Idx_All)
            Writer.Variable_("Stoich", Stoich_All)

            # Initialize RXN stoichiometry array for all molecules
            Writer.InitZeros("RXNEquationParsed", Comp.Master.NUniq_Master)

            # Update RXN stoichiometry array with Idx and Stoich participating in the RXN
            Writer.OperScUpd("RXNEquationParsed", "Idx", "Stoich")
            Writer.Reshape__("RXNEquationParsed", "RXNEquationParsed", [1, -1])

            # Add RXN Stoichiometry array to the Cel.Stoichs
            Writer.OperCncat("self.Cel.Stoichs", "self.Cel.Stoichs", "RXNEquationParsed", 0)




            # Parse Rxn rate (to be expanded)
            Writer.Variable_("RXNRateParsed", Rate)
            Writer.Reshape__("RXNRateParsed", "RXNRateParsed", [-1, 1])

            # Add RXN Rate to the Cel.Rates
            Writer.OperCncat("self.Cel.Rates", "self.Cel.Rates", "RXNRateParsed", 0)

            # Call GetReactionMolIndex, GetReactionStoich, GetReactionRate methods
            # Call AddToMasterReactionStoichs, AddToMasterReactionRates

            Writer.BlankLine()


        Writer.TF_Graph_()
        with Writer.Statement("def GetReactionRate(self):"):
            # AdjustRate
            Writer.Pass_____()
            Writer.BlankLine()

