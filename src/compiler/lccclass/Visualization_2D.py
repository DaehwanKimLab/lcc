import os
import datetime

# Comp is a short hand for CompilerData
def Write_DataExport(Writer, Comp):
    Writer.BlankLine()

    SavePath = os.path.realpath(Comp.SavePath)
    DateTime = datetime.datetime.now().strftime('%Y%m%d-%H%M')
    SaveFileName = SavePath + '/DL_EcoliSimulation_%s_Transcription.csv' % (DateTime)

    # Temporary export

    GeneNames4mRNAsToTrack_Random10 = [
        'alaS',  # Amino acyl tRNA synthesis
        # 'rplC',  # Ribosome
        'def',   # Translation
        'groS',  # Protein folding
        'dnaA',  # Replication
        'ftsZ',  # Cell division
        'nusB',  # Transcription factor
        'pyrH',  # NTP biosynthesis
        'nadE',  # NAD biosynthesis
        'murC'   # Peptidoglycan biosynthesis
    ]
    Idx_Random10mRNAs = list()
    for GeneName in GeneNames4mRNAsToTrack_Random10:
        Idx_Random10mRNAs.append(Comp.Master.ID2Idx_Master[Comp.Master.ID2ID_Gene2RNA_Master[Comp.Gene.Sym2ID_Genes[GeneName]]])

    with Writer.Statement("class FSimDataExport():"):
        with Writer.Function_("__init__"):
            Writer.Statement("# Define temporary variables for visualization purposes")
            Writer.BlankLine()

        with Writer.Function_("ExportData", "X", "Y", "XLabel", "YLabel", "Legends", "Title"):
            # Temporary export code
            with Writer.Statement("with open('%s', 'w', newline='') as SaveFile:" % SaveFileName):
                Header = list(["# of active RNAPs", "# of RNAPs Newly Bound To Gene", "# of Nascent RNAs Elongated", "# of New RNAs Generated", "# of ATP consumption", "# of CTP consumption", "# of GTP consumption", "# of UTP consumption"])
                Header += GeneNames4mRNAsToTrack_Random10
                Writer.Statement("Header = %s" % str(Header))
                Writer.Statement("Writer_Save = csv.writer(SaveFile)")
                Writer.Statement("Writer_Save.writerow(list(np.array(Header)))")
            Writer.BlankLine()

            # Temporary export
            Writer.Statement("Array = [self.Count_ElongatedTotal, tf.shape(self.Idx_RndRNAsNascent)[0], self.Count_ElongatedTotal, self.Count_ElongationCompletionTotal, self.Count_NTPConsumption[0], self.Count_NTPConsumption[1], self.Count_NTPConsumption[2], self.Count_NTPConsumption[3]]")
            for Idx_RNA in Idx_Random10mRNAs:
                Writer.Statement("Array.append(self.Cel.Counts[%s])" % Idx_RNA)
            with Writer.Statement("with open('%s', 'a+', newline='') as SaveFile:" % SaveFileName):
                Writer.Statement("Writer_Save = csv.writer(SaveFile)")
                Writer.Statement("Writer_Save.writerow(np.array(Array))")
            Writer.BlankLine()
