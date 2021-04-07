import numpy as np


# Two Component Systems

def Write_TCS_Init(writer, CompilerData):
    writer.BlankLine()
    with writer.Statement("def TCS_Init():"):

        # Set up matrices for Two Component Systems
        writer.Statement("# Two Component Systems model")
        writer.Statement("CellMX.TCSModel = tf.keras.models.load_model(LCCDataPath + '/two_component.h5')")
        writer.Variable_("TCSODETimeStep", 1)
        writer.Statement("CellMX.TCSODETimeStepTF = tf.reshape(tf.constant([TCSODETimeStep], dtype='float32'), (1, 1))")
        writer.DebugPVar("TCSODDETimeStepTF")
        writer.BlankLine()

        # Set up matrices for Two Component Systems
        TCSMolIndexList = []
        TCSMolNames = np.load('TCSMolNames.npy')
        writer.Statement("TCSMolCounts = np.zeros(" + str(len(TCSMolNames)) + ").astype('float32')")
        for i, TCSMolName in enumerate(TCSMolNames):
            TCSMolIndex = CompilerData.MolName2Index[TCSMolName]
            TCSMolIndexList.append(int(TCSMolIndex))
            writer.Statement("TCSMolCounts[%d] = MolCounts[%d] # %s" % (i, TCSMolIndex, TCSMolName))
        writer.Statement("TCSMolCountsTF = tf.convert_to_tensor(TCSMolCounts)")
        writer.Statement("CellMX.TCSMolIndexTF = tf.reshape(tf.constant(" + str(TCSMolIndexList) + "), (-1, 1))")
        writer.DebugPVar("TCSMolCountsTF")
        writer.DebugPVar("CellMX.TCSMolIndexTF")
        writer.BlankLine()
    return

def Write_TCS_Loop(writer):
    writer.BlankLine()
    with writer.Statement("def TCS_Loop():"):
        # Two Component Systems code (TCS)
        writer.Statement("# Two Component Systems code (TCS)")

        # TCS - Run machine learned model
        writer.Statement("# TCS - Run machine learned model")
        writer.Statement("TCSMolCountsTF = tf.gather(CellMX.MolCountsTF, CellMX.TCSMolIndexTF)")
        writer.Statement("TCSMolConcsTF = TCSMolCountsTF / (CellVol * AvogadroNum)")  # TCSMolConcsTF == y_init
        writer.Statement("TCSModelInput = tf.concat([TCSMolConcsTF, CellMX.TCSODETimeStepTF], axis=0)")
        writer.Statement("TCSMolConcsNewTF = CellMX.TCSModel.predict(tf.reshape(TCSModelInput, (1, -1)))[0, :]")
        writer.BlankLine()

        # TCS - Replace values < 0 to 0
        writer.Statement("# TCS - Replace values < 0 to 0")
        writer.Statement("TCSMolConcsNewZeroIndexTF = tf.where(tf.less(TCSMolConcsNewTF, 0))")
        writer.Statement("TCSMolConcsNewReplaceTF = tf.zeros(TCSMolConcsNewZeroIndexTF.shape[0])")
        writer.Statement(
            "TCSMolConcsNewTF = tf.tensor_scatter_nd_update(TCSMolConcsNewTF, TCSMolConcsNewZeroIndexTF, TCSMolConcsNewReplaceTF)")
        writer.DebugPVar("TCSMolConcsNewZeroIndexTF")
        writer.DebugPVar("TCSMolConcsNewReplaceTF")
        writer.DebugPVar("TCSMolConcsNewTF")
        writer.DebugSTMT(
            "tf.debugging.assert_non_negative(TCSMolConcsNewTF, 'TCSMolConcsNewTF contains a negative value(s)')")

        writer.DebugPVar("TCSMolCountsTF")
        writer.DebugPVar("TCSMolConcsTF")
        writer.DebugPVar("TCSMolConcsNewTF")
        writer.DebugPVar("TCSMolCountsNewTF")
        writer.BlankLine()

        # TCS - Update TCS molecule counts
        writer.Statement("# Update two component systems molecule counts")
        writer.Statement("TCSMolCountsNewTF = TCSMolConcsNewTF * (CellVol * AvogadroNum)")
        writer.Statement("CellMX.MolCountsTF = tf.tensor_scatter_nd_update(CellMX.MolCountsTF, CellMX.TCSMolIndexTF, TCSMolCountsNewTF)")
        writer.DebugSTMT("TCSMolCountsUpdated = tf.gather(CellMX.MolCountsTF, CellMX.TCSMolIndexTF)")
        # writer.DebugSTMT("tf.assert_equal(TCSMolCountsUpdated, TCSMolCountsNewTF, 'TCSMolCounts is not updated')")
        writer.PrintVari("TCSMolCountsNewTF")
        writer.BlankLine()
    return

def Write_TCS_Body(writer, CompilerData):
    # Write init
    # loop






    return