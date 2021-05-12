import numpy as np


# Metabolism

def Write_Metab_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def Metab_Init():"):
        Writer.Statement("# Metabolism (Metab)")
        # Matrices for Metabolism

        # Initialize Metabolite Concentrations
        Writer.Statement("# Metab - Initialize metabolite concentrations")
        Writer.Statement("CellMX.MetaboliteConcsResetTF = CellMX.MetaboliteConcsTF")

        # Set up Molecular Weight Matrix
        # Writer.Statement("# Metab - Set up Molecular Weight Matrix for all metabolites for Metab function")
        # MetaboliteMWs4All = np.load('MetaboliteMWs.npy')
        # MetaboliteIndexList4Metab = []
        # Writer.Statement("MetaboliteMWs = np.zeros(" + str(len(CompilerData.MetaboliteNames4Conc)) + ").astype('float32')")
        # for i, Name in enumerate(CompilerData.MetaboliteNames4Conc):
        #     MetaboliteIndex = CompilerData.MetaboliteName2MWIndex[Name]
        #     MetaboliteIndexList4Metab.append(int(MetaboliteIndex))
        #     Writer.Statement("MetaboliteMWs[%d] = %s # %s" % (i, MetaboliteMWs4All[MetaboliteIndex], Name))
        # Writer.Statement("CellMX.MetaboliteMWsTF = tf.convert_to_tensor(MetaboliteMWs)")
        # Writer.DebugPVar("CellMX.MetaboliteMWsTF")
        # Writer.BlankLine()

        # Convert metabolite concentrations into metabolite counts
        # Writer.Statement("Convert metabolite concentrations into metabolite counts")
        # Writer.Statement("CellMX.MetaboliteMWsTF = tf.transpose(CellMX.MetaboliteMWsTF")
        # Writer.Statement("CellMX.MetaboliteCountsTF = tf.matmul(CellMX.MetaboliteConcsTF, CellMX.MetaboliteMWsTF)")
        # Writer.Statement("CellMX.MetaboliteCountsTF = tf.reshape(CellMX.MetaboliteCountsTF, -1)")
        # Writer.PrintVari("CellMX.MetaboliteCountsTF")

        Writer.BlankLine()
    return

# This code needs to be improved to proportionally predict the Metabolism target for its abundance

def Write_Metab_Loop(Writer):
    Writer.BlankLine()
    with Writer.Statement("def Metab_Loop():"):
        Writer.Statement("# Metabolism (Metab)")
        Writer.PrintStrg("# Metabolism Loop outputs:")
        Writer.BlankLine()

        # Metab - Replenish all metabolites (temporary)
        Writer.Statement("# Metab - Replenish all metabolites (temporary)")
        Writer.Statement("CellMX.MetaboliteConcsTF = CellMX.MetaboliteConcsResetTF")
        Writer.PrintVari("CellMX.MetaboliteConcsTF")


        Writer.BlankLine()

    return

