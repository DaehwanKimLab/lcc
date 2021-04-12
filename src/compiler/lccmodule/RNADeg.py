import numpy as np


# RNA Degradation

def Write_RNADeg_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def RNADeg_Init():"):
        Writer.Statement("# RNA Degradation (RNADeg)")
        # Matrices for RNA Degradation

        # Set up matrices for RNA Degradation
        # Currently RNA is general, Transcript is specific to the process
        TranscriptIndexList = []
        TranscriptIDs = np.load('RNAIDs.npy') # can be a specific set to be loaded, i.e. types of RNA?
        Writer.Statement("TranscriptCounts = np.zeros(" + str(len(TranscriptIDs)) + ").astype('int32')")

        # for i, TranscriptID in enumerate(TranscriptIDs):
        #     TranscriptIndex = CompilerData.RNAID2Index[TranscriptID]
        #     TranscriptIndexList.append(int(TranscriptIndex))
        #     Writer.Statement("TranscriptCounts[%d] = RNACounts[%d] # %s" % (i, TranscriptIndex, TranscriptID))

        # Temporary replacement for loop in life code
        with Writer.Statement("for i in range(len(TranscriptCounts)):"):
            Writer.Statement("TranscriptCounts[i] = RNACounts[i]")
        for i in range(len(TranscriptIDs)):
            TranscriptIndexList.append(i)

        Writer.Statement("TranscriptCountsTF = tf.convert_to_tensor(TranscriptCounts)")
        Writer.Statement("CellMX.TranscriptIndexTF = tf.reshape(tf.constant(" + str(TranscriptIndexList) + ", dtype='int32'), [-1, 1])")
        Writer.DebugPVar("TranscriptCountsTF")
        Writer.DebugPVar("CellMX.TranscriptIndexTF")
        Writer.BlankLine()

        # Writer.Statement("# RNADeg - Determine RNA Deg rate")
        # Writer.Variable_("CellMX.RNADegRate", 10) # TO BE REPLACED AND MOVED INTO SIMULATION
        # Writer.BlankLine()

        # RNADeg - Determine active endoRNase count - to be revised with EndoRNase with specific RNA type targets
        Writer.Statement("# RNADeg - Determine active endoRNase count")

        # EndoRNaseIndexList = []
        # EndoRNaseIDs = np.load('EndoRNaseIDs.npy') # can be a specific set to be loaded, i.e. mRNA or tRNA or rRNA
        # Writer.Statement("TranscriptCounts = np.zeros(" + str(len(TranscriptIDs)) + ").astype('int32')")
        # # for i, EndoRNaseID in enumerate(EndoRNaseIDs):
        #     EndoRNaseIndex = CompilerData.EndoRNaseID2Index[EndoRNaseID]
        #     TranscriptIndexList.append(int(EndoRNaseIndex)
        #     Writer.Statement("EndoRNaseCounts[%d] = ProteinCounts[%d] # %s" % (i, EndoRNaseIndex, EndoRNaseID))

        # EndoRNase4mRNAIndexList = []
        # EndoRNase4tRNAIndexList = []
        # EndoRNase4rRNAIndexList = []
        # EndoRNase4miscRNAIndexList = []

        EndoRNase4mRNACount = 1000
        EndoRNase4tRNACount = 100
        EndoRNase4rRNACount = 200
        EndoRNase4miscRNACount = 50

        EndoRNase4mRNAAvailCount = 1000
        EndoRNase4tRNAAvailCount = 100
        EndoRNase4rRNAAvailCount = 200
        EndoRNase4miscRNAAvailCount = 50

        Writer.Variable_("CellMX.ActiveEndoRNaseCount", 100) # TO BE REPLACED AND MOVED INTO SIMULATION
        Writer.Statement("CellMX.ActiveEndoRNaseAvailCount = CellMX.ActiveEndoRNaseCount")

        Writer.BlankLine()


    return

def Write_RNADeg_Loop(Writer):
    Writer.BlankLine()
    with Writer.Statement("def RNADeg_Loop():"):
        Writer.Statement("# RNA Degradation (RNADeg)")
        Writer.PrintStrg("# RNA Degradation Loop outputs:")
        Writer.BlankLine()

        # RNA Deg - Allocate endoRNase to transcript
        Writer.Statement("# RNADeg - Allocate endoRNase to transcript")
        Writer.Statement("EndoRNasePerTranscriptTF = tf.zeros(CellMX.NumberOfUniqueRNA, dtype='int32')") # Initialize new Endonuclease targets
        with Writer.Statement("for EndoRNasePosition in range(CellMX.ActiveEndoRNaseAvailCount):"):
            Writer.Statement(
                "EndoRNasePosition = tf.random.uniform(shape=[1,1], minval=0, maxval=CellMX.NumberOfUniqueRNA, dtype='int32')")
            Writer.Statement("EndoRNasePerTranscriptTF = tf.tensor_scatter_nd_add(EndoRNasePerTranscriptTF, EndoRNasePosition, OneTF)")
        Writer.DebugAsrt("tf.math.reduce_sum(EndoRNasePerTranscriptTF) == CellMX.ActiveEndoRNaseAvailCount",
                         'Active EndoRNases are not properly allocated')
        Writer.Statement("EndoRNasePerTranscriptTF = tf.reshape(EndoRNasePerTranscriptTF, -1)")
        Writer.BlankLine()

        # RNADeg - Determine transcript degradation count
        Writer.Statement("# Update RNA Deg Transcript counts")
        Writer.Statement("CellMX.RNACountsTF = tf.tensor_scatter_nd_sub(CellMX.RNACountsTF, CellMX.TranscriptIndexTF, EndoRNasePerTranscriptTF)")
        # Writer.DebugSTMT("RNACountsUpdated = tf.gather(CellMX.RNACountsTF, CellMX.TranscriptIndexTF)")
        # Writer.DebugSTMT("tf.assert_equal(TCSMolCountsUpdated, TCSMolCountsNewTF, 'TCSMolCounts is not updated')")
        Writer.PrintVari("CellMX.RNACountsTF[:20]")
        Writer.BlankLine()

        # RNADeg - Release NTPs from Exonuclease activity


        # temporary visualization code
        Writer.Statement("CellMX.RNADeg_Transcript.append(tf.reshape(CellMX.RNACountsTF, -1).numpy())")

        Writer.BlankLine()

    return