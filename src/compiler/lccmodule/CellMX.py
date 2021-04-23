import numpy as np
# import tensorflow as tf

# Cell Matrix class

def Write_CellMX_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FCellMX():"):
        with Writer.Statement("def __init__(self):"):

            # Single number variables (subject to change)
            Writer.Statement("# Single number variables (subject to change)")
            Writer.Variable_("self.CellVol", 0)
            Writer.Variable_("self.NumberOfUniqueRNAs", 0)
            Writer.Variable_("self.NumberOfUniquemRNAs", 0)
            Writer.Variable_("self.NumberOfUniquetRNAs", 0)
            Writer.Variable_("self.NumberOfUniquerRNAs", 0)
            Writer.Variable_("self.NumberOfUniquemiscRNAs", 0)
            Writer.Variable_("self.ActiveRNAPCount", 0)
            Writer.Variable_("self.ActiveRNAPAvailCount", 0)
            Writer.Variable_("self.RNAPElongationRate", 0)
            Writer.Variable_("self.TCSModel", 0)
            Writer.Variable_("self.ActiveEndoRNaseCount", 0)
            Writer.Variable_("self.ActiveEndoRNase4mRNACount", 0)
            Writer.Variable_("self.ActiveEndoRNase4tRNACount", 0)
            Writer.Variable_("self.ActiveEndoRNase4rRNACount", 0)
            Writer.Variable_("self.ActiveEndoRNase4miscRNACount", 0)
            Writer.Variable_("self.RNADegRate", 0)

            Writer.Variable_("self.NumberOfUniqueProts", 0)
            Writer.Variable_("self.ActiveRibosomeCount", 0)
            Writer.Variable_("self.RibosomeElongationRate", 0)
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
            Writer.Variable_("self.TECompletionIndexTF", 0)
            Writer.Variable_("self.ExoRNaseTargetIndexTF", 0)
            Writer.Variable_("self.ProtIndexTF", 0)
            Writer.Variable_("self.PECompletionIndexTF", 0)

            Writer.BlankLine()

            # Data matrices
            Writer.Statement("# Data matrices")
            Writer.Variable_("self.MetaboliteConcsTF", 0)
            Writer.Variable_("self.MetaboliteConcs", 0)
            Writer.Variable_("self.MetaboliteConcsResetTF", 0)
            Writer.Variable_("self.MetaboliteMWsTF", 0)
            Writer.Variable_("self.MetaboliteCountsTF", 0)
            Writer.Variable_("self.MolCountsTF", 0)
            Writer.Variable_("self.TCSODETimeStepTF", 0)
            Writer.Variable_("self.RNACountsTF", 0)
            Writer.Variable_("self.RNAPDurationsTF", 0)
            Writer.Variable_("self.TECompletionDurationTF", 0)
            Writer.Variable_("self.RNAPPerTranscriptTF", 0)
            Writer.Variable_("self.RNADeg_Transcript", [])
            Writer.Variable_("self.RNALengthsTF", 0)
            Writer.Variable_("self.RNANTFreqsTF", 0)
            Writer.Variable_("self.TranscriptCountsTF", 0)
            Writer.Variable_("self.RNACleavedCountsTF", 0)
            Writer.Variable_("self.ProtAAFreqsTF", 0)
            Writer.Variable_("self.AAConcsTF", 0)
            Writer.Variable_("self.ProtCountsTF", 0)
            Writer.Variable_("self.RibosomeDurationsTF", 0)
            Writer.Variable_("self.PECompletionDurationTF", 0)
            Writer.Variable_("self.RibosomePerPolypeptide", 0)
            Writer.Variable_("self.ProtDeg_Transcript", [])
            Writer.Variable_("self.ProtLengthsTF", 0)
            Writer.Variable_("self.ProtAAFreqsTF", 0)
            Writer.Variable_("self.ProtCleavedCountsTF", 0)
            Writer.BlankLine()

            # Temporary
            Writer.Statement("# Temporary variables for misc processes, i.e 2D visualization")
            Writer.Variable_("self.SimStep", [])
            Writer.Variable_("self.TE_ACGU", [])
            Writer.Variable_("self.RNADeg_Transcript", [])
            Writer.Variable_("self.PE_AAs", [])
            Writer.BlankLine()





