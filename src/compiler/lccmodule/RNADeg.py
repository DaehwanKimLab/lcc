import numpy as np


# RNA Degradation

def Write_RNADeg_Init(writer, CompilerData):
    writer.BlankLine()
    with writer.Statement("def RNADeg_Init():"):
        # Matrices for RNA Degradation

        # Set up matrices for RNA Degradation
        # Currently RNA is general, Transcript is specific to the process
        TranscriptIndexList = []
        TranscriptIDs = np.load('RNAIDs.npy') # can be a specific set to be loaded, i.e. types of RNA?
        writer.Statement("TranscriptCounts = np.zeros(" + str(len(TranscriptIDs)) + ").astype('int32')")

        # for i, TranscriptID in enumerate(TranscriptIDs):
        #     TranscriptIndex = CompilerData.RNAID2Index[TranscriptID]
        #     TranscriptIndexList.append(int(TranscriptIndex))
        #     writer.Statement("TranscriptCounts[%d] = RNACounts[%d] # %s" % (i, TranscriptIndex, TranscriptID))

        # Temporary replacement for loop in life code
        with writer.Statement("for i in range(len(TranscriptCounts)):"):
            writer.Statement("TranscriptCounts[i] = RNACounts[i]")
        for i in range(len(TranscriptIDs)):
            TranscriptIndexList.append(i)

        writer.Statement("TranscriptCountsTF = tf.convert_to_tensor(TranscriptCounts)")
        writer.Statement("CellMX.TranscriptIndexTF = tf.reshape(tf.constant(" + str(TranscriptIndexList) + ", dtype='int32'), (-1, 1))")
        writer.DebugPVar("TranscriptCountsTF")
        writer.DebugPVar("CellMX.TranscriptIndexTF")
        writer.BlankLine()

        # writer.Statement("# RNADeg - Determine RNA Deg rate")
        # writer.Variable_("CellMX.RNADegRate", 10) # TO BE REPLACED AND MOVED INTO SIMULATION
        # writer.BlankLine()

        writer.Statement("# RNADeg - Determine active endonuclease count")
        writer.Variable_("CellMX.ActiveEndonucleaseCount", 500) # TO BE REPLACED AND MOVED INTO SIMULATION
        writer.BlankLine()


    return

def Write_RNADeg_Loop(writer):
    writer.BlankLine()
    with writer.Statement("def RNADeg_Loop():"):
        writer.Statement("# RNA Degradation (RNADeg)")
        writer.BlankLine()

        # RNA Deg - Allocate endonuclease to transcript
        writer.Statement("# RNADeg - Allocate endonuclease to transcript")
        writer.Statement("EndonucleasePerTranscriptTF = tf.zeros(CellMX.NumberOfUniqueTranscripts, dtype='int32')")
        with writer.Statement("for EndonucleasePosition in range(CellMX.ActiveEndonucleaseCount):"):
            writer.Statement(
                "EndonucleasePosition = tf.random.uniform(shape=[1,1], minval=0, maxval=CellMX.NumberOfUniqueTranscripts, dtype='int32')")
            writer.Statement("EndonucleasePerTranscriptTF = tf.tensor_scatter_nd_add(EndonucleasePerTranscriptTF, EndonucleasePosition, OneTF)")
        writer.DebugAsrt("tf.math.reduce_sum(EndonucleasePerTranscriptTF) == CellMX.ActiveEndonucleaseCount",
                         'Active Endonucleases are not properly allocated')
        writer.Statement("EndonucleasePerTranscriptTF = tf.reshape(EndonucleasePerTranscriptTF, [-1])")
        writer.BlankLine()

        # RNADeg - Determine transcript degradation count
        writer.Statement("# Update RNA Deg Transcript counts")
        writer.Statement("CellMX.RNACountsTF = tf.tensor_scatter_nd_sub(CellMX.RNACountsTF, CellMX.TranscriptIndexTF, EndonucleasePerTranscriptTF)")
        # writer.DebugSTMT("RNACountsUpdated = tf.gather(CellMX.RNACountsTF, CellMX.TranscriptIndexTF)")
        # writer.DebugSTMT("tf.assert_equal(TCSMolCountsUpdated, TCSMolCountsNewTF, 'TCSMolCounts is not updated')")
        writer.PrintVari("CellMX.RNACountsTF[:10]")
        writer.BlankLine()

        # temporary visualization code
        writer.Statement("CellMX.RNADeg_Transcript.append(tf.reshape(CellMX.RNACountsTF, -1).numpy())")

        writer.BlankLine()

    return