# Includes all scalar data

def Write_Class_CellState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FCellState():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Variable_("self.Species", None) # Exact or closest species
            Writer.Variable_("self.CellID", 0)
            Writer.Variable_("self.CellVol", 0)
            Writer.BlankLine()

            Writer.Statement("# Number of unique molecule species")
            Writer.Variable_("self.NUniq_Genome", 0)
            Writer.Variable_("self.NUniq_Genes", 0)
            Writer.Variable_("self.NUniq_Promoters", 0)
            Writer.Variable_("self.NUniq_RNAs", 0)
            Writer.Variable_("self.NUniq_mRNAs", 0)
            Writer.Variable_("self.NUniq_tRNAs", 0)
            Writer.Variable_("self.NUniq_rRNAs", 0)
            Writer.Variable_("self.NUniq_miscRNAs", 0)
            Writer.Variable_("self.NUniq_Proteins", 0)
            Writer.Variable_("self.NUniq_CPLXs", 0)
            Writer.Variable_("self.NUniq_Lipids", 0)
            Writer.BlankLine()

            # May be replaced with transcript-specific rates (arrays) if known
            Writer.Statement("# Rates of processes")
            Writer.Variable_("self.Rate_DNAReplication", 0)
            Writer.Variable_("self.Rate_RNASynthesis", 0)
            Writer.Variable_("self.Rate_RNADegradation", 0)
            Writer.Variable_("self.Rate_ProteinSynthesis", 0)
            Writer.Variable_("self.Rate_ProteinDegradation", 0)
            Writer.Variable_("self.Rate_RNAHalfLives", 0)
            Writer.BlankLine()

            Writer.Statement("# Number of all molecules present")
            Writer.Variable_("self.N_DNAPActive", 0)
            Writer.Variable_("self.N_DNAPActiveAvail", 0)
            Writer.Variable_("self.N_RNAPActive", 0)
            Writer.Variable_("self.N_RNAPActiveAvail", 0)
            Writer.Variable_("self.N_EndoRNaseActive", 0)
            Writer.Variable_("self.N_EndoRNaseActive4mRNA", 0)
            Writer.Variable_("self.N_EndoRNaseActive4tRNA", 0)
            Writer.Variable_("self.N_EndoRNaseActive4rRNA", 0)
            Writer.Variable_("self.N_EndoRNaseActive4miscRNA", 0)
            Writer.Variable_("self.N_EndoRNaseActiveAvail4RNA", 0)
            Writer.Variable_("self.N_EndoRNaseActiveAvail4mRNA", 0)
            Writer.Variable_("self.N_EndoRNaseActiveAvail4tRNA", 0)
            Writer.Variable_("self.N_EndoRNaseActiveAvail4rRNA", 0)
            Writer.Variable_("self.N_EndoRNaseActiveAvail4miscRNA", 0)
            Writer.Variable_("self.N_RibosomeActive", 0)
            Writer.Variable_("self.N_RibosomeActiveAvail", 0)
            Writer.Variable_("self.N_ProteaseActive", 0)
            Writer.Variable_("self.N_ProteaseActiveAvail", 0)
            Writer.BlankLine()

            Writer.Statement("# Models")
            Writer.Variable_("self.TCSModel", 0)
            Writer.Variable_("self.TCSODETimeStep", 0)
            Writer.BlankLine()

            # Arrays
            Writer.Statement("# Total lengths")
            Writer.Variable_("self.Len_Genome", 0)
            Writer.Variable_("self.Len_Genes", 0)
            Writer.Variable_("self.Len_RNAs", 0)
            Writer.Variable_("self.Len_Proteins", 0)
            Writer.Variable_("self.Len_Metabolites", 0)
            Writer.BlankLine()

            Writer.Statement("# Molecular Weights")
            Writer.Variable_("self.MW_DNAs", 0)
            Writer.Variable_("self.MW_RNAs", 0)
            Writer.Variable_("self.MW_Proteins", 0)
            Writer.Variable_("self.MW_Metabolites", 0)
            Writer.BlankLine()

            Writer.Statement("# Indexes")
            Writer.Variable_("self.Index_Genes", 0)
            Writer.Variable_("self.Index_Promoters", 0)
            Writer.Variable_("self.Index_RNAs", 0)
            Writer.Variable_("self.Index_mRNAs", 0)
            Writer.Variable_("self.Index_tRNAs", 0)
            Writer.Variable_("self.Index_rRNAs", 0)
            Writer.Variable_("self.Index_miscRNAs", 0)
            Writer.Variable_("self.Index_RNAsInitiated", 0)
            Writer.Variable_("self.Index_RNAsCompleted", 0)
            Writer.Variable_("self.Index_RNAsCleaved", 0)
            Writer.Variable_("self.Index_ExoRNaseTargets", 0)
            Writer.Variable_("self.Index_Proteins", 0)
            Writer.Variable_("self.Index_ProteinsInitiated", 0)
            Writer.Variable_("self.Index_ProteinsCompleted", 0)
            Writer.Variable_("self.Index_ProteinsCleaved", 0)
            Writer.Variable_("self.Index_Metabolites", 0)

            Writer.Variable_("self.Index_NTConcs", 0)
            Writer.Variable_("self.Index_TCSMol", 0)
            Writer.BlankLine()

            Writer.Statement("# Frequency tables")
            Writer.Variable_("self.Freq_dNTPsInGenome", 0)
            Writer.Variable_("self.Freq_NTPsInRNAs", 0)
            Writer.Variable_("self.Freq_AAinProteins", 0)
            Writer.BlankLine()

            Writer.Statement("# Durations")
            Writer.Variable_("self.Dur_RNAPsOngoing", 0)
            Writer.Variable_("self.Dur_RNAPsComplete", 0)
            Writer.Variable_("self.Dur_RibosomesOngoing", 0)
            Writer.Variable_("self.Dur_RibosomesComplete", 0)
            Writer.BlankLine()

            Writer.Statement("# Distribution per target")
            Writer.Variable_("self.Dist_DNAPsOnGenome", 0)
            Writer.Variable_("self.Dist_RNAPsOnGene", 0)
            Writer.Variable_("self.Dist_ExoRNasesOnDNA", 0)
            Writer.Variable_("self.Dist_RibosomesOnRNA", 0)
            Writer.Variable_("self.Dist_ProteasesOnProtein", 0)
            Writer.BlankLine()

            Writer.Statement("# Counts of molecules")
            Writer.Variable_("self.Count_Metabolites", 0)
            Writer.Variable_("self.Count_RNAs", 0)
            Writer.Variable_("self.Count_RNAsCompleted", 0)
            Writer.Variable_("self.Count_RNAsCleaved", 0)
            Writer.Variable_("self.Count_Proteins", 0)
            Writer.Variable_("self.Count_ProteinsCompleted", 0)
            Writer.Variable_("self.Count_ProteinsCleaved", 0)
            Writer.Variable_("self.Count_TCSMols", 0)
            Writer.BlankLine()

            Writer.Statement("# LocalCounts for applicable molecular species")
            Writer.Statement("# LocalCount_Proteins")
            Writer.Statement("# LocalCount_Lipids")

            Writer.Statement("# LocalRatio for applicable molecular species")
            Writer.Statement("# LocalRatio_Proteins")
            Writer.Statement("# LocalRatio_Lipids")


            Writer.Statement("# Concentrations of molecules")
            Writer.Variable_("self.Concs_Metabolites", 0)
            Writer.Variable_("self.Concs_MetabolitesReset", 0)

            Writer.BlankLine()

        # Load from CompilerData
        with Writer.Statement("def Initialize_Class_CellState(self):"):

            Writer.Statement("# Load compiler data for metabolites")
            Writer.LoadNP2TF("self.MW_Metabolites", "MetaboliteMWs.npy", 'float32')
            Writer.LoadNP2TF("self.Counts_Metabolites", "MetaboliteCounts.npy", 'int32')

            Writer.Statement("# Load compiler data for RNAs")

            Writer.Variable_("self.NUniq_RNAs", 0)
            Writer.Variable_("self.NUniq_mRNAs", 0)
            Writer.Variable_("self.NUniq_tRNAs", 0)
            Writer.Variable_("self.NUniq_rRNAs", 0)
            Writer.Variable_("self.NUniq_miscRNAs", 0)

            Writer.LoadNP2TF("self.MW_RNAs", "RNAMWs.npy", 'float32')
            Writer.LoadNP2TF("self.Counts_RNAs", "RNACounts.npy", 'int32')
            Writer.LoadNP2TF("self.Rate_RNAHalfLives", "RNAHalfLives.npy", 'float32')
            Writer.LoadNP2TF("self.Len_RNAs", "RNALengths.npy", 'int32')
            Writer.LoadNP2TF("self.Count_RNAs", "RNANTCounts.npy", 'int32')
            Writer.LoadNP2TF("self.Freq_NTPsInRNAs", "RNANTFreqs.npy", 'int32')
            Writer.LoadNP2TF("self.Index_RNAs", "RNATypeIndex4AllRNA.npy", 'int32')
            Writer.LoadNP2TF("self.Index_mRNAs", "RNATypeIndex4mRNA.npy", 'int32')
            Writer.LoadNP2TF("self.Index_tRNAs", "RNATypeIndex4tRNA.npy", 'int32')
            Writer.LoadNP2TF("self.Index_rRNAs", "RNATypeIndex4rRNA.npy", 'int32')
            Writer.LoadNP2TF("self.Index_miscRNAs", "RNATypeIndex4miscRNA.npy", 'int32')



'''
            # Index for all unique RNA
            Writer.Statement("# Index for all unique RNA")
            Writer.Statement(
                "CellMX.RNAIndex4AllRNATF = tf.convert_to_tensor(list(range(CellMX.NumberOfUniqueRNAs)))")
            Writer.BlankLine()

            # Load all RNA lengths
            Writer.Statement("RNALengths = np.load(\"RNALengths.npy\").astype('int32')")
            Writer.Statement("CellMX.RNALengthsTF = tf.convert_to_tensor(RNALengths)")
            Writer.BlankLine()

            # Load all Protein counts (placeholder)
            Writer.Statement("# Protein Counts")
            Writer.Variable_("DefaultCount", 100)
            Writer.Statement("ProtCounts = np.ones(" + str(len(CompilerData.ProtIDs)) + ").astype('int32') * DefaultCount")
            Writer.Statement("CellMX.ProtCountsTF = tf.convert_to_tensor(ProtCounts)")
            Writer.DebugPVar("CellMX.ProtCountsTF")
            Writer.BlankLine()

            # Load the number of unique Proteins
            Writer.Statement("CellMX.NumberOfUniqueProts = len(ProtCounts)")
            Writer.BlankLine()

            # Index for all unique Protein
            Writer.Statement("# Index for all unique Prot")
            Writer.Statement(
                "CellMX.ProtIndexTF = tf.convert_to_tensor(list(range(CellMX.NumberOfUniqueProts)))")
            Writer.BlankLine()

            # Load all Protein lengths
            Writer.Statement("ProtLengths = np.load(\"ProtLengths.npy\").astype('int32')")
            Writer.Statement("CellMX.ProtLengthsTF = tf.convert_to_tensor(ProtLengths)")
            Writer.BlankLine()
'''
