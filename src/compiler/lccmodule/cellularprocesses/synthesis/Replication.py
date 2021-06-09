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

# Comp is a short hand for CompilerData
def Write_SynDNA(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FSynDNA(FCellProcess):"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.Idx_Consumed", 0)
            Writer.Variable_("self.Idx_Produced", 0)
            Writer.Variable_("self.Rate", 0)

            Writer.BlankLine()
            Writer.Statement("super().__init__()")

            Writer.BlankLine()

        with Writer.Statement("def InitProcess(self, Cel, Cst, Env):"):
            Writer.Statement("self.ElementaryProcess(Cel, Cst, Env)")

            Writer.BlankLine()

        Writer.TF_Graph_()
        with Writer.Statement("def LoopProcess(self, Cel, Cst, Env, Sim):"):

            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def AddElementaryProcess(self, Cel, Cst, Env):"):


            RXN = '2 dNTP -> 1 ChromosomeSize + 2 PPi'
            Rate = '1000 events per second'
            Type = 'Specific'

            Writer.Comment__(RXN)
            Writer.Comment__(Rate)

            # Use regular expression to parse the reaction stoichiometry info
            # The final product would be the following:

            ID_Consumed = ['dNTP']
            ID_Produced = ['Chromosome1_Rep2', 'PPI[c]']
            Stoich_Consumed = [2]
            Stoich_Produced = [1, 2]
            Rate = 1000

            # Parse the molecule inputs for indexing
            ID_Consumed, Stoich_Consumed = Writer.RefineBBs(ID_Consumed, Stoich_Consumed, Comp)
            ID_Produced, Stoich_Produced = Writer.RefineBBs(ID_Produced, Stoich_Produced, Comp)

            # Find Indexes for Molecules
            Idx_Consumed = Writer.GetMolIdx(ID_Consumed, Comp.Master.ID2Idx_Master)
            Idx_Produced = Writer.GetMolIdx(ID_Produced, Comp.Master.ID2Idx_Master)

            # Convert index and stoichiometry to tensor for simulation
            Writer.Variable_("Idx_Consumed", Idx_Consumed)
            Writer.Variable_("Idx_Produced", Idx_Produced)

            # Initialize RXN stoichiometry array for all molecules
            Writer.InitZeros("Stoich_RXN", Comp.Master.NUniq_Master)

            "Cel.MasterStoich"
            "Cel.MasterRate"

            # Call GetReactionMolIndex, GetReactionStoich, GetReactionRate methods
            # Call AddToMasterReactionStoichs, AddToMasterReactionRates

            Writer.BlankLine()

        with Writer.Statement("def GetReactionMolIndex(self, Cel):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def GetReactionStoich(self, Cel):"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.TF_Graph_()
        with Writer.Statement("def GetReactionRate(self, Cel):"):
            # AdjustRate
            Writer.Pass_____()
            Writer.BlankLine()

