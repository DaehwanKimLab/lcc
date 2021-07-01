# Interface for all lcc modules, a.k.a. cellular processes


# Comp is a short hand for CompilerData
def Write_CellProcess(Writer, ProGen):
    ProGen.GenerateCellProcessInterface(Writer)

        # with Writer.Statement("class FCellProcess():"):
    #     with Writer.Statement("def __init__(self):"):
    #         Writer.BlankLine()
    #
    #         Writer.Statement("super().__init__()")
    #         Writer.BlankLine()
    #
    #     Writer.AbsMethod()
    #     with Writer.Statement("def AddToStoichiometryMatrix(self):"):
    #         Writer.Pass_____()
    #         Writer.BlankLine()
    #
    #     Writer.AbsMethod()
    #     with Writer.Statement("def CalculateRate(self):"):
    #         Writer.Pass_____()
    #         Writer.BlankLine()
    #
    #     Writer.AbsMethod()
    #     with Writer.Statement("def AddToRateMatrix(self):"):
    #         Writer.Pass_____()
    #         Writer.BlankLine()
    #
    #     Writer.TF_Graph_()
    #     with Writer.Statement("def UpdateRates(self):"):
    #         Writer.Statement("self.Cel.ClearRateMatrix()")
    #         Writer.Statement("self.AddToReactionRateMatrix()")
    #         Writer.BlankLine()



"""
Reaction (RXN) Types

== Deterministic (Specific) RXN ==
Deterministic RXNs occur to specific molecules or all members of the molecule group participating in the reaction
Features:
    - RXN Matrix: STATIC
    - Rate Matrix: DYNAMIC
    - Indexing: STATIC
Examples:
    - Metabolic networks
    - Equilibrium reaction

== Stochastic (General) RXN ==
Stochastic RXNs occur to a subset of the molecule group participating in the reaction
RXN matrix must be set up for all possible reactions, where the rate matrix can be controlled
Features:
    - RXN Matrix: STATIC
    - Rate Matrix: DYNAMIC
    - Indexing: DYNAMIC (apply the same rate to all chosen indexes)
Examples:
    - Replication
    - Transcription (NTP Consumption and RNA Production for each RNA)
    - Protein degradation (AA production and Protein_Cleaved count decreasing during Protein degradation)


== Conditional & Deterministic RXN ==
a.k.a. Utility RXN for Deterministic RXN setup by the compiler
Conditional & Deterministic RXNs occur to specific molecules or all members of the molecule group participating in the reaction,
based on a static condition coded by the compiler, lcc.py
Features:
    - Reaction Matrix: STATIC
    - Rate Matrix: DYNAMIC (all or none, based on the condition)
    - Indexing: STATIC
Examples:
    - Gene count during replication based on the DNAP position in the current simulation

== Conditional & Stochastic RXN ==
a.k.a. Utility RXN for Stochastic RXN setup by the compiler
Conditional & Stochastic RXNs occur to a subset of the molecule group participating in the reaction,
based on a dynamic condition coded by the simulator, cell.py
Features:
    - RXN Matrix: STATIC
    - Rate Matrix: DYNAMIC (all or none, based on the condition)
    - Indexing: DYNAMIC to switch on and off Rate matrix
Examples:
    - Transcription (transcription initiation upon sigma factor binding)

"""

# 'Type'
# 'RXNID'
# 'MoleculeIDs'
# 'Stoichiometry'
# 'Reversibility'
# 'Order'
# 'Rate'
# 'ModulatorIDs'
# 'Trigger'

# ReactionType = ['Non-biochemical'] # Consider boolean for biochemical vs. non-biochemical
# Condition = ['DNAPolymerase']  # Enzymes for biochemical reactions
# ReactionTrigger = ['']
# ID_Consumed = ['dNTP']
# ID_Produced = ['Chromosome1_Rep2', 'PPI[c]']
# Stoich_Consumed = [2]
# Stoich_Produced = [1, 2]
# ReactionOrder = [0]
# Rate = [1000.0]