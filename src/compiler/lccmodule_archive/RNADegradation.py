import numpy as np


# RNA Degradation

def Write_RNADeg_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def RNADeg_Init():"):
        Writer.Statement("# RNA Degradation (RNADeg)")
        # Matrices for RNA Degradation

        # Set up matrices for RNA Degradation
        # # Currently RNA is general, Transcript is specific to the process
        # TranscriptIndexList = list()
        # TranscriptIDs = np.load('RNAIDs.npy') # can be a specific set to be loaded, i.e. types of RNA?
        # Writer.Statement("TranscriptCounts = np.zeros(" + str(len(TranscriptIDs)) + ").astype('int32')")
        #
        # # for i, TranscriptID in enumerate(TranscriptIDs):
        # #     TranscriptIndex = CompilerData.RNAID2Index[TranscriptID]
        # #     TranscriptIndexList.append(int(TranscriptIndex))
        # #     Writer.Statement("TranscriptCounts[%d] = RNACounts[%d] # %s" % (i, TranscriptIndex, TranscriptID))
        #
        # # Temporary replacement for loop in life code
        # with Writer.Statement("for i in range(len(TranscriptCounts)):"):
        #     Writer.Statement("TranscriptCounts[i] = RNACounts[i]")
        # for i in range(len(TranscriptIDs)):
        #     TranscriptIndexList.append(i)
        #
        # Writer.Statement("TranscriptCountsTF = tf.convert_to_tensor(TranscriptCounts)")
        # Writer.Statement("CellMX.TranscriptIndexTF = tf.reshape(tf.constant(" + str(TranscriptIndexList) + ", dtype='int32'), [-1, 1])")
        # Writer.DebugPVar("TranscriptCountsTF")
        # Writer.DebugPVar("CellMX.TranscriptIndexTF")
        # Writer.BlankLine()

        # RNADeg - Determine active endoRNase count - to be revised with EndoRNase with specific RNA type targets
        Writer.Statement("# RNADeg - Determine active endoRNase count")

        # EndoRNaseIndexList = list()
        # EndoRNaseIDs = np.load('EndoRNaseIDs.npy') # can be a specific set to be loaded, i.e. mRNA or tRNA or rRNA
        # Writer.Statement("TranscriptCounts = np.zeros(" + str(len(TranscriptIDs)) + ").astype('int32')")
        # # for i, EndoRNaseID in enumerate(EndoRNaseIDs):
        #     EndoRNaseIndex = CompilerData.EndoRNaseID2Index[EndoRNaseID]
        #     TranscriptIndexList.append(int(EndoRNaseIndex)
        #     Writer.Statement("EndoRNaseCounts[%d] = ProteinCounts[%d] # %s" % (i, EndoRNaseIndex, EndoRNaseID))

        # EndoRNase4mRNAIndexList = list()
        # EndoRNase4tRNAIndexList = list()
        # EndoRNase4rRNAIndexList = list()
        # EndoRNase4miscRNAIndexList = list()

        # EndoRNase4mRNACount = 1000
        # EndoRNase4tRNACount = 100
        # EndoRNase4rRNACount = 200
        # EndoRNase4miscRNACount = 50

        Writer.Variable_("CellMX.ActiveEndoRNase4mRNAAvailCount", 500) # TO BE REPLACED
        Writer.Variable_("CellMX.ActiveEndoRNase4tRNAAvailCount", 0) # TO BE REPLACED
        Writer.Variable_("CellMX.ActiveEndoRNase4rRNAAvailCount", 0) # TO BE REPLACED
        Writer.Variable_("CellMX.ActiveEndoRNase4miscRNAAvailCount", 0) # TO BE REPLACED

        Writer.Statement("CellMX.ActiveEndoRNaseCount += CellMX.ActiveEndoRNase4mRNAAvailCount")
        Writer.Statement("CellMX.ActiveEndoRNaseCount += CellMX.ActiveEndoRNase4tRNAAvailCount")
        Writer.Statement("CellMX.ActiveEndoRNaseCount += CellMX.ActiveEndoRNase4rRNAAvailCount")
        Writer.Statement("CellMX.ActiveEndoRNaseCount += CellMX.ActiveEndoRNase4miscRNAAvailCount")

        Writer.Variable_("CellMX.ActiveExoRNaseCount", 10)

        Writer.InitArrayWithZero("CellMX.ExoRNaseTargetIndexTF", 0, 'int32')

        Writer.BlankLine()


    return

# This code needs to be improved to proportionally predict the RNA Degradation target for its abundance
def Write_RNADeg_Loop(Writer):
    Writer.BlankLine()
    with Writer.Statement("def RNADeg_Loop():"):
        Writer.Statement("# RNA Degradation (RNADeg)")
        Writer.PrintStrg("# RNA Degradation Loop outputs:")
        Writer.BlankLine()

        # RNA Deg - Allocate endoRNase to specific substrates
        Writer.Statement("EndoRNasePerTranscriptTF = tf.zeros(CellMX.NumberOfUniqueRNAs, dtype='int32')")
        FourRNATypes = ['mRNA', 'tRNA', 'rRNA', 'miscRNA']
        for RNAType in FourRNATypes:
            Writer.RndIncrmt("EndoRNasePerTranscriptTF", "CellMX.ActiveEndoRNase4" + RNAType + "AvailCount", "CellMX.RNAIndex4" + RNAType + "TF", "1")
            Writer.BlankLine()
        Writer.DebugAsrt("tf.math.reduce_sum(EndoRNasePerTranscriptTF) == CellMX.ActiveEndoRNaseCount",
                         'Active EndoRNases are not properly allocated')
        Writer.Statement("EndoRNasePerTranscriptTF = tf.reshape(EndoRNasePerTranscriptTF, -1)")
        Writer.BlankLine()

        # RNADeg - Update RNA Deg Transcript counts and RNA cleaved counts
        Writer.Statement("# Update RNA Deg Transcript counts")
        Writer.Statement("CellMX.RNACountsTF = tf.tensor_scatter_nd_sub(CellMX.RNACountsTF, CellMX.RNAIndex4AllRNATF, EndoRNasePerTranscriptTF)")
        # Writer.DebugSTMT("RNACountsUpdated = tf.gather(CellMX.RNACountsTF, CellMX.TranscriptIndexTF)")
        # Writer.DebugSTMT("tf.assert_equal(TCSMolCountsUpdated, TCSMolCountsNewTF, 'TCSMolCounts is not updated')")
        Writer.PrintVari("CellMX.RNACountsTF[:20]")
        Writer.BlankLine()

        # RNADeg - Add RNA cleaved counts
        Writer.Statement("# RNADeg - Add RNA cleaved counts")
        Writer.Statement("CellMX.RNACleavedCountsTF += EndoRNasePerTranscriptTF")
        Writer.PrintVari("CellMX.RNACleavedCountsTF[:20]")
        Writer.BlankLine()

        # RNADeg - Find index for all RNA cleaved thus far
        Writer.Statement("# RNADeg - Find index for all RNA cleaved thus far")
        Writer.Statement("CellMX.ExoRNaseTargetIndexTF = tf.where(tf.math.greater(CellMX.RNACleavedCountsTF, 0))")
        Writer.Statement("CellMX.ExoRNaseTargetIndexTF = tf.cast(CellMX.ExoRNaseTargetIndexTF, dtype='int32')")
        Writer.Statement("CellMX.ExoRNaseTargetIndexTF = tf.reshape(CellMX.ExoRNaseTargetIndexTF, -1)")
        Writer.PrintVari("CellMX.ExoRNaseTargetIndexTF")
        Writer.BlankLine()

        # RNADeg - Allocate ExoRNases to the EndoRNase-processed RNAs
        Writer.RndIncrmt("CellMX.RNACleavedCountsTF", "CellMX.ActiveExoRNaseCount", "CellMX.ExoRNaseTargetIndexTF", "-1")
        # Writer.WgtIncrmt("CellMX.RNACleavedCountsTF", "CellMX.ActiveExoRNaseCount", "CellMX.ExoRNaseTargetIndexTF", "-1")
        Writer.DebugSTMT("tf.debugging.assert_non_negative(CellMX.RNACleavedCountsTF, message='RNA Cleaved Counts became negative')")
        Writer.PrintVari("CellMX.RNACleavedCountsTF[:20]")
        Writer.BlankLine()


        # temporary visualization code
        # Writer.Statement("CellMX.RNADeg_Transcript.append(tf.reshape(CellMX.RNACountsTF, -1).numpy())")

        Writer.BlankLine()

    return