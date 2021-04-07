import numpy as np


# Transcript Elongation.

def Write_TE_Init(writer, CompilerData):
    writer.BlankLine()
    with writer.Statement("def TE_Init():"):
        # Matrices for Transcript Elongation
        np.save("TranscriptNTFreqs", CompilerData.TranscriptNTFreqs)
        writer.Statement("# Matrices for Transcript Elongation")

        # NT frequency table for transcripts
        writer.Statement("# Fetch NT frequency table for transcripts")
        writer.Statement("TranscriptNTFreqs = np.load(\"TranscriptNTFreqs.npy\").astype('float32')")
        writer.Statement("CellMX.TranscriptNTFreqsTF = tf.convert_to_tensor(TranscriptNTFreqs)")
        writer.Statement("CellMX.TranscriptNTFreqsTF = tf.transpose(TranscriptNTFreqs)")
        writer.DebugPVar("TranscriptNTFreqs")
        writer.DebugPVar("CellMX.TranscriptNTFreqsTF")
        writer.Statement("CellMX.NumberOfUniqueTranscripts = len(TranscriptNTFreqs)")
        writer.BlankLine()

        # NT counts per transcript
        NTIndexList = list()
        writer.Statement("# Fetch NT concentration")
        ACGU = ["ATP", "CTP", "GTP", "UTP"]
        writer.Statement("NTConcs = np.zeros(" + str(len(ACGU)) + ").astype('float32')")
        for i, NTName in enumerate(ACGU):
            NTIndex = CompilerData.MetaboliteName2Index[NTName]
            NTIndexList.append(int(NTIndex))
            writer.Statement("NTConcs[%d] = CellMX.MetaboliteConcs[%d] # %s" % (i, NTIndex, NTName))
        writer.Statement("NTConcsTF = tf.convert_to_tensor(NTConcs)")
        writer.Statement("CellMX.NTConcsIndexTF = tf.reshape(tf.constant(" + str(NTIndexList) + "), [4, -1])")
        writer.DebugPVar("NTConcs")
        writer.DebugPVar("NTConcsTF")
        writer.DebugPVar("CellMX.NTConcsIndexTF")
        writer.BlankLine()

        writer.Statement("# Determine elongation rate")
        writer.Variable_("CellMX.ElongationRate", 10) # TO BE REPLACED AND MOVED INTO SIMULATION
        writer.BlankLine()

        writer.Statement("# Determine active RNAP count")
        writer.Variable_("CellMX.ActiveRNAPCount", 829) # TO BE REPLACED AND MOVED INTO SIMULATION
        writer.BlankLine()
    return

def Write_TE_Loop(writer):
    writer.BlankLine()
    with writer.Statement("def TE_Loop():"):
        writer.Statement("# Transcript Elongation (TE)")
        writer.BlankLine()

        # TE - Allocate RNAP to transcript
        writer.Statement("# TE - Allocate RNAP to transcript")
        writer.Statement("RNAPPerTranscriptTF = tf.zeros(CellMX.NumberOfUniqueTranscripts)")
        with writer.Statement("for RNAPPosition in range(CellMX.ActiveRNAPCount):"):
            writer.Statement(
                "RNAPPosition = tf.random.uniform(shape=[1,1], minval=1, maxval=CellMX.NumberOfUniqueTranscripts, dtype='int32')")
            writer.Statement("RNAPPerTranscriptTF = tf.tensor_scatter_nd_add(RNAPPerTranscriptTF, RNAPPosition, OneTF)")
        writer.DebugAsrt("tf.math.reduce_sum(RNAPPerTranscriptTF) == CellMX.ActiveRNAPCount",
                         'Active RNAP is not properly allocated')
        writer.Statement("RNAPPerTranscriptTF = tf.reshape(RNAPPerTranscriptTF, [-1, 1])")
        writer.DebugPVar("RNAPPerTranscriptTF")
        writer.BlankLine()

        # TE - Determine NT consumption
        writer.Statement("# TE - Determine NT consumption")
        writer.Statement(
            "DeltaNTCountsTF = tf.linalg.matmul(CellMX.TranscriptNTFreqsTF, RNAPPerTranscriptTF) * CellMX.ElongationRate")
        writer.Statement("DeltaNTCountsTF = tf.reshape(DeltaNTCountsTF, -1)")
        writer.Statement("DeltaNTConcsTF = DeltaNTCountsTF / (CellVol * AvogadroNum)")  # final unit: mol/L
        writer.DebugVari("NTConcsAvailTF", "tf.gather(CellMX.MetaboliteConcsTF, CellMX.NTConcsIndexTF)")
        writer.DebugPVar("DeltaNTCountsTF")
        writer.DebugPVar("DeltaNTConcsTF")
        writer.DebugPVar("NTConcsAvailTF")
        writer.DebugSTMT(
            "tf.debugging.assert_positive(NTConcsAvailTF - DeltaNTConcsTF), 'The cell is running out of NTs'")
        writer.BlankLine()

        # TE - Update NT concs
        writer.Statement("# TE - Update NT counts")
        writer.Statement(
            "CellMX.MetaboliteConcsTF = tf.tensor_scatter_nd_sub(CellMX.MetaboliteConcsTF, CellMX.NTConcsIndexTF, DeltaNTConcsTF)")
        writer.DebugVari("NTConcsNewTF", "tf.gather(CellMX.MetaboliteConcsTF, CellMX.NTConcsIndexTF)")
        writer.DebugPVar("NTConcsNewTF")
        writer.DebugSTMT(
            "tf.debugging.assert_none_equal(NTConcsNewTF, NTConcsAvailTF), 'NT consumption is not properly applied'")
        writer.PrintVari("NTConcsNewTF")
        # writer.Statement("TE_Y.append(tf.reshape(NTConcsNewTF, -1).numpy())")
        writer.BlankLine()

        # TE - Update Transcript counts - TO BE IMPLEMENTED
    return