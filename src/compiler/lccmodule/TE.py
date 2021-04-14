import numpy as np


# Transcript Elongation.

def Write_TE_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def TE_Init():"):
        # Matrices for Transcript Elongation
        np.save("RNANTFreqs", CompilerData.RNANTFreqs)
        np.save("RNALengths", CompilerData.RNALengths)
        
        Writer.Statement("# Matrices for Transcript Elongation")
        # NT frequency table for RNA
        Writer.Statement("# Fetch NT frequency table for transcripts")
        Writer.Statement("RNANTFreqs = np.load(\"RNANTFreqs.npy\").astype('float32')")
        Writer.Statement("CellMX.RNANTFreqsTF = tf.convert_to_tensor(RNANTFreqs)")
        Writer.Statement("CellMX.RNANTFreqsTF = tf.transpose(RNANTFreqs)")
        Writer.DebugPVar("RNANTFreqs")
        Writer.DebugPVar("CellMX.RNANTFreqsTF")
        Writer.Statement("CellMX.NumberOfUniqueRNA = len(RNANTFreqs)")
        Writer.BlankLine()
        # RNA lengths
        Writer.Statement("# Fetch RNA lengths")
        Writer.Statement("RNALengths = np.load(\"RNALengths.npy\").astype('float32')")
        Writer.Statement("CellMX.RNALengthsTF = tf.convert_to_tensor(RNALengths)")
        Writer.BlankLine()

        # NT counts per transcript
        NTIndexList = list()
        Writer.Statement("# Fetch NT concentration")
        ACGU = ["ATP", "CTP", "GTP", "UTP"]
        Writer.Statement("NTConcs = np.zeros(" + str(len(ACGU)) + ").astype('float32')")
        for i, NTName in enumerate(ACGU):
            NTIndex = CompilerData.MetaboliteName2Index[NTName]
            NTIndexList.append(int(NTIndex))
            Writer.Statement("NTConcs[%d] = CellMX.MetaboliteConcs[%d] # %s" % (i, NTIndex, NTName))
        Writer.Statement("NTConcsTF = tf.convert_to_tensor(NTConcs)")
        Writer.Statement("CellMX.NTConcsIndexTF = tf.reshape(tf.constant(" + str(NTIndexList) + "), [4, -1])")
        Writer.DebugPVar("CellMX.NTConcsIndexTF")
        Writer.BlankLine()

        # Determine elongation rate
        Writer.Statement("# Determine elongation rate")
        Writer.Variable_("CellMX.ElongationRate", 10) # TO BE REPLACED AND MOVED INTO SIMULATION
        Writer.BlankLine()

        # Determine active RNAP count
        Writer.Statement("# Determine active RNAP count")
        Writer.Variable_("CellMX.ActiveRNAPCount", 829) # TO BE REPLACED AND MOVED INTO SIMULATION
        Writer.Statement("CellMX.ActiveRNAPAvailCount = CellMX.ActiveRNAPCount") # Initialize
        Writer.Statement("CellMX.RNAPPerTranscriptTF = tf.zeros(CellMX.NumberOfUniqueRNA, dtype='int32')")
        Writer.Statement("CellMX.RNAPDurationsTF = tf.zeros(CellMX.NumberOfUniqueRNA, dtype='int32')")
        Writer.BlankLine()

        # Index for all unique RNA
        Writer.Statement("# Index for all unique RNA")
        Writer.Statement(
            "CellMX.RNAIndex4AllRNATF = tf.convert_to_tensor(list(range(CellMX.NumberOfUniqueRNA)))")
        Writer.Statement(
            "CellMX.RNAIndex4AllRNATF = tf.reshape(CellMX.RNAIndex4AllRNATF, [-1, 1])")
        Writer.BlankLine()

        # Elongation completion duration for each RNA
        Writer.Statement("# Elongation completion duration for each RNA")
        Writer.Statement("ElongCompletionDuration = RNALengths // CellMX.ElongationRate")
        Writer.Statement("CellMX.ElongCompletionDurationTF = tf.convert_to_tensor(ElongCompletionDuration, dtype='int32')")
        Writer.Statement("CellMX.ElongCompletionDurationTF = tf.reshape(CellMX.ElongCompletionDurationTF, -1)")
        Writer.BlankLine()

    return

def Write_TE_Loop(Writer):
    Writer.BlankLine()
    with Writer.Statement("def TE_Loop():"):
        Writer.Statement("# Transcript Elongation (TE)")
        Writer.PrintStrg("# Transcript Elongation Loop outputs:")
        Writer.BlankLine()

        # TE - Allocate RNAP to transcript
        Writer.Statement("# TE - Allocate RNAP to transcript")
        with Writer.Statement("if CellMX.ActiveRNAPAvailCount > 0:"):
            Writer.RndIncrmt("CellMX.RNAPPerTranscriptTF", "CellMX.ActiveRNAPAvailCount", "CellMX.RNAIndex4AllRNATF", "1")
            # with Writer.Statement("for RNAPPosition in range(CellMX.ActiveRNAPAvailCount):"):
            #     Writer.Statement(
            #         "RNAPPosition = tf.random.uniform(shape=[1,1], minval=0, maxval=CellMX.NumberOfUniqueRNA, dtype='int32')")
            #     Writer.Statement("CellMX.RNAPPerTranscriptTF = tf.tensor_scatter_nd_add(CellMX.RNAPPerTranscriptTF, RNAPPosition, OneTF)")
            Writer.Statement("CellMX.ActiveRNAPAvailCount = 0")
            Writer.DebugAsrt("tf.math.reduce_sum(CellMX.RNAPPerTranscriptTF) == CellMX.ActiveRNAPCount",
                             'Active RNAPs are not properly allocated')
            Writer.DebugPVar("CellMX.RNAPPerTranscriptTF")
            Writer.BlankLine()

        # TE - Determine NT consumption
        Writer.Statement("# TE - Determine NT consumption")
        Writer.Statement("RNAPPerTranscriptTF_Float = tf.cast(tf.reshape(CellMX.RNAPPerTranscriptTF, [-1, 1]), dtype='float32')")
        Writer.Statement(
            "DeltaNTCountsTF = tf.linalg.matmul(CellMX.RNANTFreqsTF, RNAPPerTranscriptTF_Float) * CellMX.ElongationRate")
        Writer.Statement("DeltaNTCountsTF = tf.reshape(DeltaNTCountsTF, -1)")
        Writer.Statement("DeltaNTConcsTF = DeltaNTCountsTF / (CellVol * AvogadroNum)")  # final unit: mol/L
        Writer.DebugVari("NTConcsAvailTF", "tf.gather(CellMX.MetaboliteConcsTF, CellMX.NTConcsIndexTF)")
        Writer.DebugPVar("DeltaNTCountsTF")
        Writer.DebugPVar("DeltaNTConcsTF")
        Writer.DebugPVar("NTConcsAvailTF")
        Writer.DebugSTMT(
            "tf.debugging.assert_positive(NTConcsAvailTF - DeltaNTConcsTF), 'The cell is running out of NTs'")
        Writer.BlankLine()

        # TE - Determine RNA production and RNAP release
        Writer.Statement("# TE - Determine RNA production and RNAP release")
        Writer.Statement(
            "CellMX.RNAPDurationsTF = tf.tensor_scatter_nd_add(CellMX.RNAPDurationsTF, CellMX.RNAIndex4AllRNATF, tf.reshape(CellMX.RNAPPerTranscriptTF, -1))")
        Writer.Statement("CellMX.RNAPDurationsTF = tf.reshape(CellMX.RNAPDurationsTF, -1)")
        Writer.BlankLine()

        Writer.Statement("# Determine Index for Elongation Completion")
        Writer.Statement("ElongCompletionIndexTF = tf.where(tf.math.greater_equal(CellMX.RNAPDurationsTF, CellMX.ElongCompletionDurationTF))")
        Writer.Statement("ElongCompletionIndexTF = tf.reshape(ElongCompletionIndexTF, [-1, 1])")
        Writer.BlankLine()

        # TE - Update NT concs
        Writer.Statement("# TE - Update NT counts")
        Writer.Statement(
            "CellMX.MetaboliteConcsTF = tf.tensor_scatter_nd_sub(CellMX.MetaboliteConcsTF, CellMX.NTConcsIndexTF, DeltaNTConcsTF)")
        Writer.DebugVari("NTConcsNewTF", "tf.reshape(tf.gather(CellMX.MetaboliteConcsTF, CellMX.NTConcsIndexTF), [-1])")
        Writer.DebugSTMT(
            "tf.debugging.assert_none_equal(NTConcsNewTF, NTConcsAvailTF), 'NT consumption is not properly applied'")
        Writer.PrintVari("NTConcsNewTF")
        Writer.BlankLine()

        # If elongation was complete, update RNAP duration on RNA, RNA Counts, endoRNase counts
        Writer.Statement("# Updates upon elongation completion")
        # Writer.PrintVari("ElongCompletionIndexTF")
        with Writer.Statement("if tf.math.count_nonzero(ElongCompletionIndexTF):"):
            # Update RNAP duration on RNA
            Writer.Statement("# Update RNAP duration on RNA")
            Writer.Statement("DeltaDurationsTF = tf.gather(CellMX.ElongCompletionDurationTF, ElongCompletionIndexTF)")
            Writer.Statement("DeltaDurationsTF = tf.reshape(DeltaDurationsTF, -1)")
            # Writer.PrintStrg("RNAPDuration with Elongation completion index before update:")
            # Writer.PrintVari("tf.reshape(tf.gather(CellMX.RNAPDurationsTF, ElongCompletionIndexTF), -1)")
            Writer.Statement("CellMX.RNAPDurationsTF = tf.tensor_scatter_nd_sub(CellMX.RNAPDurationsTF, ElongCompletionIndexTF, DeltaDurationsTF)")
            # Writer.PrintStrg("RNAPDuration with Elongation completion index after update:")
            # Writer.PrintVari("tf.reshape(tf.gather(CellMX.RNAPDurationsTF, ElongCompletionIndexTF), -1)")
            Writer.BlankLine()

            # Update RNAPPerTranscript
            Writer.Statement("# Update RNAPPerTranscript")
            # Writer.PrintStrg("RNAPPerTranscript with Elongation completion index before update:")
            # Writer.PrintVari("tf.reshape(tf.gather(CellMX.RNAPPerTranscriptTF, ElongCompletionIndexTF), -1)")
            Writer.InitArrayWithOne("ElongCompletionOnesTF", "tf.size(ElongCompletionIndexTF)", 'int32')
            Writer.Statement("CellMX.RNAPPerTranscriptTF = tf.tensor_scatter_nd_sub(CellMX.RNAPPerTranscriptTF, ElongCompletionIndexTF, ElongCompletionOnesTF)")
            # Writer.PrintStrg("RNAPPerTranscript with Elongation completion index after update:")
            # Writer.PrintVari("tf.reshape(tf.gather(CellMX.RNAPPerTranscriptTF, ElongCompletionIndexTF), -1)")
            Writer.BlankLine()

            # TE - Update RNA counts
            Writer.Statement("# TE - Update RNA counts")
            Writer.Statement(
                "CellMX.RNACountsTF = tf.tensor_scatter_nd_add(CellMX.RNACountsTF, ElongCompletionIndexTF, ElongCompletionOnesTF)")

            Writer.BlankLine()

            # TE - Update free active RNAP counts
            Writer.Statement("# TE - Update free active RNAP counts")
            Writer.Statement("CellMX.ActiveRNAPAvailCount = tf.size(ElongCompletionIndexTF)")
            Writer.BlankLine()

        Writer.PrintVari("CellMX.ElongCompletionDurationTF[:20]")
        Writer.PrintVari("CellMX.RNAPDurationsTF[:20]")
        Writer.PrintVari("CellMX.RNACountsTF[:20]")
        Writer.PrintVari("CellMX.ActiveRNAPAvailCount")

        # temporary visualization code
        Writer.Statement("CellMX.TE_ACGU.append(tf.reshape(NTConcsNewTF, -1).numpy())")


        Writer.BlankLine()

        # TE - Update Transcript counts - TO BE IMPLEMENTED
    return