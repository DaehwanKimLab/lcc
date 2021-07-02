import numpy as np


# Two Component Systems

def Write_TCS_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def TCS_Init():"):
        Writer.Statement("# Two Component Systems (TCS)")
        Writer.BlankLine()

        # Set up matrices for Two Component Systems
        Writer.Statement("# Two Component Systems model")
        Writer.Statement("CellMX.TCSModel = tf.keras.models.load_model(LCCDataPath + '/two_component.h5')")
        Writer.Variable_("TCSODETimeStep", 1)
        Writer.Statement("CellMX.TCSODETimeStepTF = tf.reshape(tf.constant([TCSODETimeStep], dtype='float32'), [1, 1])")
        Writer.DebugPVar("TCSODDETimeStepTF")
        Writer.BlankLine()

        # Set up matrices for Two Component Systems
        TCSMolIndexList = list()
        TCSMolNames = np.load('TCSMolNames.npy')
        Writer.Statement("TCSMolCounts = np.zeros(" + str(len(TCSMolNames)) + ").astype('float32')")

        # for i, TCSMolName in enumerate(TCSMolNames):
        #     TCSMolIndex = CompilerData.MolName2Index[TCSMolName]
        #     TCSMolIndexList.append(int(TCSMolIndex))
        #     Writer.Statement("TCSMolCounts[%d] = MolCounts[%d] # %s" % (i, TCSMolIndex, TCSMolName))

        # Temporary replacement for loop in life code
        with Writer.Statement("for i in range(len(TCSMolCounts)):"):
            Writer.Statement("TCSMolCounts[i] = MolCounts[i]")
        for i in range(len(TCSMolNames)):
            TCSMolIndexList.append(i)

        Writer.Statement("TCSMolCountsTF = tf.convert_to_tensor(TCSMolCounts)")
        Writer.Statement("CellMX.TCSMolIndexTF = tf.reshape(tf.constant(" + str(TCSMolIndexList) + "), [-1, 1])")
        Writer.DebugPVar("TCSMolCountsTF")
        Writer.DebugPVar("CellMX.TCSMolIndexTF")
        Writer.BlankLine()
    return

def Write_TCS_Loop(Writer):
    Writer.BlankLine()
    with Writer.Statement("def TCS_Loop():"):
        # Two Component Systems code (TCS)
        Writer.Statement("# Two Component Systems code (TCS)")
        Writer.PrintStrg("# Two Component Systems Loop outputs:")

        # TCS - Run machine learned model
        Writer.Statement("# TCS - Run machine learned model")
        Writer.OperGathr("TCSMolCountsTF", "CellMX.MolCountsTF", "CellMX.TCSMolIndexTF")
        Writer.Statement("TCSMolConcsTF = TCSMolCountsTF / (CellMX.CellVol * AvogadroNum)")  # TCSMolConcsTF == y_init
        Writer.Statement("TCSModelInput = tf.concat([TCSMolConcsTF, CellMX.TCSODETimeStepTF], axis=0)")
        Writer.Statement("TCSModelInput = tf.reshape(TCSModelInput, [1, -1])")
        Writer.Statement("TCSMolConcsNewTF = CellMX.TCSModel(TCSModelInput)[0, :]")
        Writer.BlankLine()

        # TCS - Replace values < 0 to 0
        Writer.Statement("# TCS - Replace values < 0 to 0")
        Writer.Statement("TCSMolConcsNewZeroIndexTF = tf.where(tf.less(TCSMolConcsNewTF, 0))")
        Writer.Statement("TCSMolConcsNewReplaceTF = tf.zeros(tf.shape(TCSMolConcsNewZeroIndexTF)[0])")
        Writer.Statement(
            "TCSMolConcsNewTF = tf.tensor_scatter_nd_update(TCSMolConcsNewTF, TCSMolConcsNewZeroIndexTF, TCSMolConcsNewReplaceTF)")
        Writer.DebugPVar("TCSMolConcsNewZeroIndexTF")
        Writer.DebugPVar("TCSMolConcsNewReplaceTF")
        Writer.DebugPVar("TCSMolConcsNewTF")
        Writer.DebugSTMT(
            "tf.debugging.assert_non_negative(TCSMolConcsNewTF, 'TCSMolConcsNewTF contains a negative value(s)')")

        Writer.DebugPVar("TCSMolCountsTF")
        Writer.DebugPVar("TCSMolConcsTF")
        Writer.DebugPVar("TCSMolConcsNewTF")
        Writer.DebugPVar("TCSMolCountsNewTF")
        Writer.BlankLine()

        # TCS - Update TCS molecule counts
        Writer.Statement("# Update two component systems molecule counts")
        Writer.Statement("TCSMolCountsNewTF = TCSMolConcsNewTF * (CellMX.CellVol * AvogadroNum)")
        Writer.Statement("CellMX.MolCountsTF = tf.tensor_scatter_nd_update(CellMX.MolCountsTF, CellMX.TCSMolIndexTF, TCSMolCountsNewTF)")
        # Writer.DebugSTMT("TCSMolCountsUpdated = tf.gather(CellMX.MolCountsTF, CellMX.TCSMolIndexTF)")
        # Writer.DebugSTMT("tf.assert_equal(TCSMolCountsUpdated, TCSMolCountsNewTF, 'TCSMolCounts is not updated')")
        Writer.PrintVari("TCSMolCountsNewTF[:10]")
        Writer.BlankLine()
    return

def Write_TCS_Body(Writer, CompilerData):
    # Write init
    # loop






    return