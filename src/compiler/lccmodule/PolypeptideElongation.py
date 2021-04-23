import numpy as np


# Polypeptide Elongation.

def Write_PE_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def PE_Init():"):
        # Matrices for Polypeptide Elongation

        Writer.Statement("# Matrices for Polypeptide Elongation")
        # AA frequency table for proteins
        Writer.Statement("# Fetch AA frequency table for Polypeptides")
        Writer.Statement("ProtAAFreqs = np.load(\"ProtAAFreqs.npy\").astype('float32')")
        Writer.Statement("CellMX.ProtAAFreqsTF = tf.convert_to_tensor(ProtAAFreqs)")
        Writer.Statement("CellMX.ProtAAFreqsTF = tf.transpose(ProtAAFreqs)")
        Writer.DebugPVar("ProtAAFreqs")
        Writer.DebugPVar("CellMX.ProtAAFreqsTF")
        Writer.Statement("CellMX.NumberOfUniqueProt = len(ProtAAFreqs)")
        Writer.BlankLine()

        # AA counts per polypeptide
        AAIndexList = list()
        Writer.Statement("# Fetch AA concentration")
        Writer.Statement("AACounts = np.zeros(" + str(len(CompilerData.AAsWithLoc)) + ").astype('float32')")
        for i, AAName in enumerate(CompilerData.AAsWithLoc):
            assert AAName in (CompilerData.MetaboliteName2CountIndex or CompilerData.MetaboliteName2ConcIndex), "" + str(AAName) + " not found in CellMX.MetaboliteCounts"
            # Gly[c] has been arbitrarily set to 100000 in DL_0_Metabolism_BulkMolecules.tsv
            # L-SELENOCYSTEINE[c] has been arbitrarily set to 0 in DL_0_Metabolism_BulkMolecules.tsv
            if AAName in CompilerData.MetaboliteName2CountIndex:
                AAIndex = CompilerData.MetaboliteName2CountIndex[AAName]
                AAIndexList.append(int(AAIndex))
                Writer.Statement("AACounts[%d] = CellMX.MetaboliteCounts[%d] # %s" % (i, AAIndex, AAName))
            else:
                AAName = CompilerData.AAs[i]
                AAIndex = CompilerData.MetaboliteName2ConcIndex[AAName]
                Writer.Statement("AACounts[%d] = CellMX.MetaboliteCounts[%d] * CellMX.CellVol * AvogadroNum # %s" % (i, AAIndex, AAName))
        Writer.Statement("AAConcs = AACounts / CellMX.CellVol / AvogadroNum")
        Writer.Statement("AAConcsTF = tf.convert_to_tensor(AAConcs)")
        Writer.Statement("CellMX.AAConcsIndexTF = tf.constant(" + str(AAIndexList) + ")")
        Writer.DebugPVar("CellMX.AAConcsIndexTF")
        Writer.BlankLine()

        # Determine elongation rate
        Writer.Statement("# Determine elongation rate")
        Writer.Variable_("CellMX.RibosomeElongationRate", 5) # TO BE REPLACED AND MOVED INTO SIMULATION
        Writer.BlankLine()

        # Determine active Ribosome count and Ribosome-related matrices
        Writer.Statement("# Determine active Ribosome count")
        Writer.Variable_("CellMX.ActiveRibosomeCount", 13408) # TO BE REPLACED AND MOVED INTO SIMULATION
        Writer.Statement("CellMX.ActiveRibosomeAvailCount = CellMX.ActiveRibosomeCount") # Initialize
        Writer.Statement("CellMX.RibosomePerPolypeptideTF = tf.zeros(CellMX.NumberOfUniqueProts, dtype='int32')")
        Writer.Statement("CellMX.RibosomeDurationsTF = tf.zeros(CellMX.NumberOfUniqueProts, dtype='int32')")
        Writer.BlankLine()

        # Elongation completion duration for each polypeptide
        Writer.Statement("# Elongation completion duration for each protein")
        Writer.Statement("ElongCompletionDuration = CellMX.ProtLengthsTF // CellMX.RibosomeElongationRate")
        Writer.Statement("CellMX.PECompletionDurationTF = tf.convert_to_tensor(ElongCompletionDuration, dtype='int32')")
        Writer.BlankLine()

    return

def Write_PE_Loop(Writer):
    Writer.BlankLine()
    with Writer.Statement("def PE_Loop():"):
        Writer.Statement("# Polypeptide Elongation (PE)")
        Writer.PrintStrg("# Polypeptide Elongation Loop outputs:")
        Writer.BlankLine()

        # PE - Allocate Ribosomes to polypeptide
        Writer.Statement("# PE - Allocate Ribosomes to Polypeptide")
        with Writer.Statement("if CellMX.ActiveRibosomeAvailCount > 0:"):
            Writer.RndIncrmt("CellMX.RibosomePerPolypeptideTF", "CellMX.ActiveRibosomeAvailCount", "CellMX.ProtIndexTF", "1")
            Writer.Statement("CellMX.ActiveRibosomeAvailCount = 0")
            Writer.DebugAsrt("tf.math.reduce_sum(CellMX.RibosomePerPolypeptideTF) == CellMX.ActiveRibosomeCount",
                             'Active Ribosomes are not properly allocated')
            Writer.DebugPVar("CellMX.RibosomePerPolypeptideTF")
            Writer.BlankLine()

        # PE - Determine AA consumption
        Writer.Statement("# PE - Determine AA consumption")
        Writer.Statement("RibosomePerPolypeptideTF_Float = tf.cast(tf.reshape(CellMX.RibosomePerPolypeptideTF, [-1, 1]), dtype='float32')")
        Writer.OperMXMul("DeltaAACountsTF", "CellMX.ProtAAFreqsTF", "RibosomePerPolypeptideTF_Float")
        Writer.Statement("DeltaAACountsTF *= CellMX.RibosomeElongationRate")
        Writer.Statement("DeltaAACountsTF = tf.reshape(DeltaAACountsTF, -1)")
        Writer.Statement("DeltaAAConcsTF = DeltaAACountsTF / (CellMX.CellVol * AvogadroNum)")  # final unit: mol/L

        Writer.DebugVari("AAConcsAvailTF", "tf.gather(CellMX.MetaboliteConcsTF, CellMX.AAConcsIndexTF)")
        Writer.DebugPVar("DeltaAACountsTF")
        Writer.DebugPVar("DeltaAAConcsTF")
        Writer.DebugPVar("AAConcsAvailTF")
        Writer.DebugSTMT(
            "tf.debugging.assert_positive(AAConcsAvailTF - DeltaAAConcsTF), 'The cell is running out of AAs'")
        Writer.BlankLine()

        # PE - Determine Prot production and Ribosome release
        Writer.Statement("# PE - Determine Prot production and Ribosome release")
        Writer.OperScAdd("CellMX.RibosomeDurationsTF", "CellMX.ProtIndexTF", "CellMX.RibosomePerPolypeptideTF")
        Writer.Statement("CellMX.RibosomeDurationsTF = tf.reshape(CellMX.RibosomeDurationsTF, -1)")
        Writer.BlankLine()

        Writer.Statement("# Determine Index for Elongation Completion")
        Writer.Statement("CellMX.PECompletionIndexTF = tf.where(tf.math.greater_equal(CellMX.RibosomeDurationsTF, CellMX.PECompletionDurationTF))")
        Writer.Statement("CellMX.PECompletionIndexTF = tf.reshape(CellMX.PECompletionIndexTF, [-1, 1])")
        Writer.BlankLine()

        # PE - Update AA concs
        Writer.Statement("# PE - Update AA counts")
        Writer.OperScSub("CellMX.MetaboliteConcsTF", "CellMX.AAConcsIndexTF", "DeltaAAConcsTF")
        Writer.DebugVari("AAConcsNewTF", "tf.reshape(tf.gather(CellMX.MetaboliteConcsTF, CellMX.AAConcsIndexTF), [-1])")
        Writer.DebugSTMT(
            "tf.debugging.assert_none_equal(AAConcsNewTF, AAConcsAvailTF), 'AA consumption is not properly applied'")
        Writer.PrintVari("AAConcsNewTF")
        Writer.BlankLine()

        # If elongation was complete, update Ribosome duration on Polypeptides, Prot Counts, free Ribosome counts
        Writer.Statement("# Updates upon elongation completion")
        # Writer.PrintVari("CellMX.PECompletionIndexTF")
        with Writer.Statement("if tf.math.count_nonzero(CellMX.PECompletionIndexTF):"):
            # Update Ribosome duration on Polypeptide
            Writer.Statement("# Update Ribosome duration on Polypeptide")
            Writer.OperGathr("DeltaDurationsTF", "CellMX.PECompletionDurationTF", "CellMX.PECompletionIndexTF")
            Writer.Statement("DeltaDurationsTF = tf.reshape(DeltaDurationsTF, -1)")
            Writer.PrintStrg("RibosomeDuration with Elongation completion index before update:")
            Writer.PrintVari("tf.reshape(tf.gather(CellMX.RibosomeDurationsTF, CellMX.PECompletionIndexTF), -1)")
            Writer.OperScSub("CellMX.RibosomeDurationsTF", "CellMX.PECompletionIndexTF", "DeltaDurationsTF")
            Writer.PrintStrg("RibosomeDuration with Elongation completion index after update:")
            Writer.PrintVari("tf.reshape(tf.gather(CellMX.RibosomeDurationsTF, CellMX.PECompletionIndexTF), -1)")
            Writer.BlankLine()

            # Update RibosomePerPolypeptide
            Writer.Statement("# Update RibosomePerPolypeptide")
            Writer.PrintStrg("RibosomePerPolypeptide with Elongation completion index before update:")
            Writer.PrintVari("tf.reshape(tf.gather(CellMX.RibosomePerPolypeptideTF, CellMX.PECompletionIndexTF), -1)")
            Writer.InitArrayWithOne("ElongCompletionOnesTF", "tf.size(CellMX.PECompletionIndexTF)", 'int32')
            Writer.OperScSub("CellMX.RibosomePerPolypeptideTF", "CellMX.PECompletionIndexTF", "ElongCompletionOnesTF")
            Writer.PrintStrg("RibosomePerPolypeptide with Elongation completion index after update:")
            Writer.PrintVari("tf.reshape(tf.gather(CellMX.RibosomePerPolypeptideTF, CellMX.PECompletionIndexTF), -1)")
            Writer.BlankLine()

            # PE - Update Protein counts
            Writer.Statement("# PE - Update Protein counts")
            Writer.OperScAdd("CellMX.ProtCountsTF", "CellMX.PECompletionIndexTF", "ElongCompletionOnesTF")
            Writer.BlankLine()

            # PE - Update free active Ribosome counts
            Writer.Statement("# PE - Update free active Ribosome counts")
            Writer.Statement("CellMX.ActiveRibosomeAvailCount = tf.size(CellMX.PECompletionIndexTF)")
            Writer.BlankLine()

        Writer.PrintVari("CellMX.PECompletionDurationTF[:20]")
        Writer.PrintVari("CellMX.RibosomeDurationsTF[:20]")
        Writer.PrintVari("CellMX.ProtCountsTF[:20]")
        Writer.PrintVari("CellMX.ActiveRibosomeAvailCount")

        # temporary visualization code
        Writer.Statement("CellMX.PE_AAs.append(tf.reshape(AAConcsNewTF, -1).numpy())")

        Writer.BlankLine()

        # PE - Update Polypeptide counts - TO BE IMPLEMENTED
    return