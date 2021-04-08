import numpy as np
# import tensorflow as tf

# Cell Matrix class

def Write_CellMX_Init(writer):
    writer.BlankLine()
    with writer.Statement("class FCellMX():"):
        with writer.Statement("def __init__(self):"):

            # Single number variables (subject to change)
            writer.Statement("# Single number variables (subject to change)")
            writer.Variable_("self.NumberOfUniqueTranscripts", 0)
            writer.Variable_("self.ActiveRNAPCount", 0)
            writer.Variable_("self.ElongationRate", 0)
            writer.Variable_("self.TCSModel", 0)
            writer.Variable_("self.ActiveEndonucleaseCount", 0)
            writer.Variable_("self.RNADegRate", 0)
            writer.BlankLine()

            # Matrices
            writer.Statement("# Matrices")
            writer.Variable_("self.TranscriptNTFreqsTF", 0)
            writer.Variable_("self.MetaboliteConcsTF", 0)
            writer.Variable_("self.MetaboliteConcs", 0)
            writer.Variable_("self.NTConcsIndexTF", 0)
            writer.Variable_("self.MolCountsTF", 0)
            writer.Variable_("self.TCSMolIndexTF", 0)
            writer.Variable_("self.TCSODETimeStepTF", 0)
            writer.Variable_("self.RNACountsTF", 0)
            writer.Variable_("self.TranscriptCountsTF", 0)
            writer.Variable_("self.TranscriptIndexTF", 0)
            writer.BlankLine()

            # Temporary - these empty lists do not seem to show up? currently moved to lcc.py
            writer.Statement("# Temporary variables")
            writer.Variable_("self.SimStep", [])
            writer.Variable_("self.TE_ACGU", [])
            writer.Variable_("self.RNADeg_Transcript", [])
            writer.BlankLine()





