# Interface for all lcc modules, a.k.a. cellular processes

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

# Comp is a short hand for CompilerData
def Write_CellProcess_Init(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FCellProcess():"):
        with Writer.Statement("def __init__(self):"):
            Writer.BlankLine()
            Writer.Statement("super().__init__()")

            Writer.BlankLine()

        # Abstract Methods for CellProcess
        Writer.AbsMethod()
        with Writer.Statement("def InitProcess(self, Cel, Cst, Env):"):
            # Call AddElementaryProcess
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def LoopProcess(self, Cel, Cst, Env, Sim):"):
            # Call UpdateReactionRate
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def AddElementaryProcess(self, Cel, Cst, Env):"):
            # Call GetReactionMolIndex, GetReactionStoich, GetReactionRate methods
            # Call AddToMasterReactionStoichs, AddToMasterReactionRates
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def GetReactionMolIndex(self, Cel):"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def GetReactionStoich(self, Cel):"):
            Writer.Pass_____()
            Writer.BlankLine()

        Writer.AbsMethod()
        with Writer.Statement("def GetReactionRate(self, Cel):"):
            # AdjustRate
            Writer.Pass_____()
            Writer.BlankLine()

        # Methods
        with Writer.Statement("def AdjustReactionRate(self, Cel):"):
            Writer.Pass_____()

            Writer.BlankLine()

        with Writer.Statement("def AddToMasterReactionStoichs(self, Cel):"):
            Writer.Pass_____()

            Writer.BlankLine()

        with Writer.Statement("def AddToMasterReactionRates(self, Cel):"):
            Writer.Pass_____()

            Writer.BlankLine()

        Writer.TF_Graph_()
        with Writer.Statement("def UpdateMasterReactionRates(self, Cel):"):
            Writer.Pass_____()

            Writer.BlankLine()

