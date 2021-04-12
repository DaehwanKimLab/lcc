import numpy as np
# import tensorflow as tf

# Cell Matrix class

def Write_CellMX_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FCellMX():"):
        with Writer.Statement("def __init__(self):"):

            # Single number variables (subject to change)
            Writer.Statement("# Single number variables (subject to change)")
            Writer.Variable_("self.NumberOfUniqueRNA", 0)
            Writer.Variable_("self.NumberOfUniquemRNA", 0)
            Writer.Variable_("self.NumberOfUniquetRNA", 0)
            Writer.Variable_("self.NumberOfUniquerRNA", 0)
            Writer.Variable_("self.NumberOfUniquemiscRNA", 0)
            Writer.Variable_("self.ActiveRNAPCount", 0)
            Writer.Variable_("self.ActiveRNAPAvailCount", 0)
            Writer.Variable_("self.ElongationRate", 0)
            Writer.Variable_("self.TCSModel", 0)
            Writer.Variable_("self.ActiveEndoRNaseCount", 0)
            Writer.Variable_("self.ActiveEndoRNase4mRNACount", 0)
            Writer.Variable_("self.ActiveEndoRNase4tRNACount", 0)
            Writer.Variable_("self.ActiveEndoRNase4rRNACount", 0)
            Writer.Variable_("self.ActiveEndoRNase4miscRNACount", 0)
            Writer.Variable_("self.RNADegRate", 0)
            Writer.BlankLine()

            # Matrices
            Writer.Statement("# Matrices")
            Writer.BlankLine()
            # Index matrices
            Writer.Statement("# Index matrices")
            Writer.Variable_("self.NTConcsIndexTF", 0)
            Writer.Variable_("self.TCSMolIndexTF", 0)
            Writer.Variable_("self.RNAIndex4AllRNATF", 0)
            Writer.Variable_("self.RNAIndex4mRNATF", 0)
            Writer.Variable_("self.RNAIndex4tRNATF", 0)
            Writer.Variable_("self.RNAIndex4rRNATF", 0)
            Writer.Variable_("self.RNAIndex4miscRNATF", 0)
            Writer.Variable_("self.TranscriptIndexTF", 0)

            Writer.BlankLine()

            # Data matrices
            Writer.Statement("# Data matrices")

            Writer.Variable_("self.MetaboliteConcsTF", 0)
            Writer.Variable_("self.MetaboliteConcs", 0)
            Writer.Variable_("self.MolCountsTF", 0)
            Writer.Variable_("self.TCSODETimeStepTF", 0)
            Writer.Variable_("self.RNACountsTF", 0)
            Writer.Variable_("self.RNAPDurationsTF", 0)
            Writer.Variable_("self.ElongCompletionDurationTF", 0)
            Writer.Variable_("self.RNAPPerTranscriptTF", 0)
            Writer.Variable_("self.RNADeg_Transcript", [])
            Writer.Variable_("self.RNALengthsTF", 0)
            Writer.Variable_("self.RNANTFreqsTF", 0)
            Writer.Variable_("self.TranscriptCountsTF", 0)
            Writer.BlankLine()

            # Temporary - these empty lists do not seem to show up? currently moved to lcc.py
            Writer.Statement("# Temporary variables")
            Writer.Variable_("self.SimStep", [])
            Writer.Variable_("self.TE_ACGU", [])
            Writer.Variable_("self.RNADeg_Transcript", [])
            Writer.BlankLine()





