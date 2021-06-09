# dna damage: Use DNA polymerase V (and DNA polymerase I) counts for repair
# restriction enzymes for surveillance

'''
https://en.wikipedia.org/wiki/DNA_polymerase_V#Function
SOS Response
SOS response in E. coli attempts to alleviate the effect of a damaging stress in the cell. The role of Pol V in SOS response triggered by UV-radiation is described as follows:

Pol III stalls at lesion site.
DNA replication helicase DnaB continues to expand the replication fork creating single stranded DNA (ssDNA) segments ahead of from the lesion.
ssDNA binding proteins (SSBs) stabilize ssDNA.
RecA recruited and loaded onto ssDNA by RecFOR replacing SSBs. Formation of RecA nucleoprotein filament (RecA*).
RecA functions through mediator proteins to activate Pol V (see Regulation).
Pol V accesses 3'-OH of nascent DNA strand and extends strand past the lesion site.
Pol III resumes elongation.[8]
'''

def Write_DegDNA(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FDegDNA():"):
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