import os, sys
import csv
import numpy as np
import ast
import re
from abc import abstractmethod

class FDataset():
    def __init__(self):
        self.Switch4DebugDataset = False


    @abstractmethod
    def SetUpData(self, Dataset, MasterDataset = None):
        pass

    @abstractmethod
    def SaveData(self, SavePath):
        pass

    def AddLocalizationToMolName(self, VariableName, Localization):
        assert VariableName[-3] != '[' and VariableName[-1] != ']', 'Variable already has a localization info to be added'
        assert Localization[0] == '[' and Localization[2] == ']', 'Localization String must be in the format of "[x]"'
        return VariableName + Localization

    def RemoveLocalizationFromMolName(self, VariableName):
        assert VariableName[-3] == '[' and VariableName[-1] == ']', 'Variable has no localization info to remove'
        return VariableName[:-3]



    def AddToMaster(self, MasterDataset, IDs, MolType, Counts, MWs, Localizations='None'):
        # self.AddToMaster(self.ID_Metabolites, self.Count_Metabolites, self.MW_Metabolites)

        for ID in IDs:
            MasterDataset.ID2Idx_Master[ID] = len(MasterDataset.ID_Master)
            MasterDataset.ID_Master.append(ID)
            MasterDataset.ID2Type_Master[ID] = MolType
            MasterDataset.Type_Master.append(MolType)


        MasterDataset.Count_Master = np.append(MasterDataset.Count_Master, Counts)
        MasterDataset.MW_Master = np.append(MasterDataset.MW_Master, MWs)
        # if Localizations:
        #     self.MasterLocalizations = np.append(self.MasterLocalizations, Localizations)

        if MasterDataset.Switch4DebugMasterDataset:
            for Key, Value in self.__dict__.items():
                if (Value == IDs) and (Key[:6] != 'Master'):
                    print('MasterDataset addition completed: Class "%s"' % Key[3:])

                if Key[:6] == 'Master':
                    LengthOfValue = len(Value)
                    print('The length of "%s": %d' % (Key, LengthOfValue))

        MasterDataset.NUniq_Master = len(MasterDataset.ID_Master)

        assert len(MasterDataset.ID_Master) == len(MasterDataset.ID2Idx_Master) and len(MasterDataset.ID_Master) == len(MasterDataset.Count_Master) and len(MasterDataset.ID_Master) == len(MasterDataset.MW_Master)

    def AppendStr(self, List, Str):
        # Appends Str to the List
        ListStr = []
        for i in List:
            i += Str
            ListStr.append(i)
        return ListStr

class FMetabolite(FDataset):
    def __init__(self):
        self.MolType_Metabolites = 'Metabolite'
        self.ID2MW_Metabolites_Temp = {}
        self.ID_Metabolites = []
        self.ID2Idx_Metabolites = {}
        self.Count_Metabolites = 0
        self.MW_Metabolites = 0
        self.NUniq_Metabolites = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        MW_Metabolites = Dataset['metabolites.tsv']
        MW_Metabolites = Dataset['metabolites.tsv']
        WaterMW = Dataset['water.tsv'][0]
        MW_Metabolites.append(WaterMW)

        for Value in MW_Metabolites:
            Name, MW, Localization = Value
            assert Name not in self.ID2MW_Metabolites_Temp
            self.ID2MW_Metabolites_Temp[Name] = MW

        Count_Metabolites = Dataset['Counts4AllMetabolites_DL_2_totalIs_0_0_Metabolism_BulkMolecules.tsv']
        self.NUniq_Metabolites = len(Count_Metabolites)
        self.Count_Metabolites = np.zeros(self.NUniq_Metabolites)
        self.MW_Metabolites = np.zeros(self.NUniq_Metabolites)

        for i, Value in enumerate(Count_Metabolites):
            Name, Count, Request = Value
            NameLocRemoved = self.RemoveLocalizationFromMolName(Name)
            assert Name not in self.ID2Idx_Metabolites
            self.ID2Idx_Metabolites[Name] = len(self.ID_Metabolites)  # = i
            self.ID_Metabolites.append(Name)
            self.Count_Metabolites[i] = Count
            self.MW_Metabolites[i] = self.ID2MW_Metabolites_Temp[NameLocRemoved]
            assert self.MW_Metabolites[i] != 0
        
        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.ID_Metabolites, self.MolType_Metabolites, self.Count_Metabolites, self.MW_Metabolites)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FChromosome(FDataset):
    def __init__(self):
        self.MolType_Chromosomes = 'Chromosome'
        self.Len_ChromosomesInGenome = 0
        self.Count_NTsInChromosomesInGenome = []
        self.Freq_NTsInChromosomesInGenome = []
        self.Genome = None

        self.NMax_Chromosomes = 0
        self.ID_Chromosomes = []
        self.Count_BasePairsInChromosomes = 0
        self.MW_DNABasePairs = 0
        self.Freq_NTsInChromosomes = []
        self.NUniq_ChromosomesInGenome = 0

        super().__init__()

    def SetUpData(self, Dataset, MasterDataset = None):
        Genome = Dataset['EscherichiaColi.fasta']
        # Genome is a dictionary where its key is 'Ch#' and its value is its corresponding DNA seq.

        self.Genome = Genome
        self.NUniq_ChromosomesInGenome = len(self.Genome)
        self.Len_ChromosomesInGenome = np.zeros(self.NUniq_ChromosomesInGenome)

        for i, (ChromosomeName, ChromosomeSeq) in enumerate(self.Genome.items()):
            self.Len_ChromosomesInGenome[i] = len(ChromosomeSeq)

            # ChromosomeNTFreq
            NTCounts = np.array([ChromosomeSeq.count("A"), ChromosomeSeq.count("C"), ChromosomeSeq.count("G"), ChromosomeSeq.count("T")])
            if np.sum(NTCounts) != self.Len_ChromosomesInGenome[i]:
                print("WARNING: RNA seq may contain non-ACGT text.", file=sys.stderr)
            NTFreq = NTCounts / self.Len_ChromosomesInGenome[i]
            self.Count_NTsInChromosomesInGenome.append(NTCounts)
            self.Freq_NTsInChromosomesInGenome.append(NTFreq)

        # Set up an arbitrary max number to keep NT Count arrays for all chromosomes
        MaxConcurrentReplications = 1
        MaxCopyNumberOfEachChromosome = 2 ** MaxConcurrentReplications
        self.NMax_Chromosomes = self.NUniq_ChromosomesInGenome * MaxCopyNumberOfEachChromosome
        self.Count_BasePairsInChromosomes = np.zeros(self.NMax_Chromosomes)
        self.MW_DNABasePairs = np.zeros(self.NMax_Chromosomes)
        DNABasePairMW = 650 # average molecular weight of one base pair: 650 Daltons, which is 650 g/mol

        for i in range(self.NMax_Chromosomes): # i = all chromosomes (original and replicates)
            for j in range(self.NUniq_ChromosomesInGenome): # j is the reference to the original set of chromosomes
                ChromosomeNumber = j + 1
                ChromosomeReplication = i + 1
                ChromosomeID = 'Chromosome%d_Rep%d' % (ChromosomeNumber, ChromosomeReplication)
                self.ID_Chromosomes.append(ChromosomeID)
                if i == 0:
                    self.Count_BasePairsInChromosomes[j] = self.Len_ChromosomesInGenome[j]
                self.MW_DNABasePairs[i] = DNABasePairMW
                self.Freq_NTsInChromosomes.append(self.Freq_NTsInChromosomesInGenome[j])

        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.ID_Chromosomes, self.MolType_Chromosomes, self.Count_BasePairsInChromosomes, self.MW_DNABasePairs)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FGene(FDataset):
    def __init__(self):
        self.MolType_Genes = 'Gene'

        self.Len_Genes = 0
        self.Coord_Genes = 0 # Coord for Coordinate
        self.MW_Genes = 0
        self.Dir_Genes = 0 # Dir for Direction
        self.Name_Genes = []
        self.Sym_Genes = [] # Sym for Symbols
        self.ID_Genes = []
        self.Seq_Genes = []
        self.Type_Genes = []
        self.ID2Idx_Genes = {}
        self.Sym2Idx_Genes = {}
        self.Name2ID_Genes = {}
        self.Count_Genes = 0
        self.NUniq_Genes = 0

        super().__init__()

    def SetUpData(self, Dataset, MasterDataset = None):
        Genes = Dataset['genes.tsv']
        self.NUniq_Genes = len(Genes)
        self.Len_Genes = np.zeros(self.NUniq_Genes)
        self.Coord_Genes = np.zeros(self.NUniq_Genes)
        self.MW_Genes = np.zeros(self.NUniq_Genes) # No molecular weights to be calculated
        self.Dir_Genes = np.zeros(self.NUniq_Genes)
        self.Count_Genes = np.ones(self.NUniq_Genes) # The count may be adjusted according to the starting replication state of the chromosomes.
        self.Name_Genes = []
        self.Sym_Genes = []
        self.ID_Genes = []
        self.Type_Genes = []
        self.Name2ID_Genes = {}
        self.ID2Idx_Genes = {}
        self.Sym2Idx_Genes = {}

        DirectionBinaryDict = {}
        DirectionBinaryDict['+'] = 1
        DirectionBinaryDict['-'] = 0

        for i, Value in enumerate(Genes):
            Length, Name, Seq, RNAID, Coordinate, Direction, Symbol, Type, GeneID, MonomerID = Value
            # assert Name not in self.Name2ID_Genes, '"%s" has a redundant name in GeneName2ID dictionary' % Name
            # self.Name_Genes.append(Name)
            # self.Name2ID_Genes[Name] = GeneID
            self.Sym_Genes.append(Symbol)
            self.Sym2Idx_Genes[Symbol] = len(self.ID_Genes)
            assert GeneID not in self.ID2Idx_Genes
            self.ID_Genes.append(GeneID)
            self.ID2Idx_Genes[GeneID] = len(self.Seq_Genes)
            self.Seq_Genes.append(Seq)
            self.Len_Genes[i] = (len(Seq))
            self.Coord_Genes[i] = Coordinate
            DirectionBinary = DirectionBinaryDict[Direction]
            self.Dir_Genes[i] = DirectionBinary
            self.Type_Genes.append(Type)

        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.ID_Genes, self.MolType_Genes, self.Count_Genes, self.MW_Genes)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FPromoter(FDataset):
    def __init__(self):
        self.MolType_Promoters = 'Promoter'
        self.Coord_Promoters = 0
        self.Dir_Promoters = 0
        self.ID_Promoters = []
        self.ID2Idx_Promoters = {}
        self.Sym_PromoterTargetGenes = []
        self.NUniq_Promoters = 0

        super().__init__()

    def SetUpData(self, Dataset, MasterDataset = None):
        Promoters = Dataset['promoters.tsv']

        self.NUniq_Promoters = len(Promoters)
        self.Coord_Promoters = np.zeros(self.NUniq_Promoters)
        self.Dir_Promoters = np.zeros(self.NUniq_Promoters)
        IrregularNames = 0

        DirectionBinaryDict = {}
        DirectionBinaryDict['+'] = 1
        DirectionBinaryDict['-'] = 0

        for i, Value in enumerate(Promoters):
            Position, Direction, PromoterID, Name = Value
            self.ID_Promoters.append(PromoterID)
            self.ID2Idx_Promoters = len(self.Dir_Promoters)  # i
            DirectionBinary = DirectionBinaryDict[Direction]
            self.Dir_Promoters[i] = DirectionBinary
            self.Coord_Promoters[i] = Position
            # Parse Name to identify promoter's target gene. Cases to cover:
            # "grxD"
            # "cadBp"
            # "exutP1"
            # "glgCp3"
            # "insA-5p"
            # "p_WC_ypeC"
            # "p_WC_gfcE-etp-etk"
            # "p_WC_insC-4D-4-ygeONM-insCD-4""

            if self.Switch4DebugDataset:
                if 'ins' in Name:
                    print(Name)
                if len(Name.split('-')) > 1 and not Name[:5] == 'p_WC_':
                    print("Promoter '%s' contains '-' in XXXp promoter name format: %s" % (PromoterID, Name))
                if len(Name.split('-')) > 2:
                    print("Promoter '%s' contains more than two '-' in the promoter name: %s" % (PromoterID, Name))

            TargetGeneSymbol = None
            TargetGeneSymbolsOnly = []
            if Name[:5] == 'p_WC_':
                Name = Name[5:]
                N_SplitName = len(Name.split('-'))
                if N_SplitName == 1:
                    TargetGeneSymbol = Name

                # DL: INCOMPLETE Exception handling:
                    # p_WC_insA-1B-1-insAB-1
                    # p_WC_insA-2B-2-afuBC-insAB-2
                    # p_WC_insA-3B-3-insAB-3
                    # p_WC_insA-4B-4-insAB-4
                    # p_WC_insC-1D-1-insCD-1
                    # p_WC_insC-2D-2-insCD-2
                    # p_WC_insC-3D-3-insCD-3
                    # p_WC_insC-4D-4-ygeONM-insCD-4
                    # p_WC_insC-6D-6-insCD-6
                    # p_WC_insE-2F-2-insEF-2
                    # p_WC_insE-3F-3-insEF-3
                    # p_WC_insE-4F-4-insEF-4
                    # p_WC_insE-5F-5-insEF-5
                    # p_WC_insO-2-yjhWV
                    # p_WC_eaeH-insE-1F-1-insEF-1
                # elif N_SplitName == 2:
                #
                #     TargetGeneSymbol =
                #
                #     TargetGeneSymbol = []
                else:
                    for Item in Name.split('-'):
                        TargetGeneSymbolsOnly.append(Item)
                    TargetGeneSymbol = [TargetGeneSymbolsOnly]
            else:
                for i in (1, 2, 3):
                    if Name[-i] == ('p' or 'P'):
                        TargetGeneSymbol = Name[:-i]
                        continue

            if not TargetGeneSymbol:
                if self.Switch4DebugDataset:
                    print("Promoter %s contains an irregular name format to parse: %s" % (PromoterID, Name))
                IrregularNames += 1
                TargetGeneSymbol = Name
            # Check if the gene is expressed by multiple promoters
            if self.Switch4DebugDataset:
                if TargetGeneSymbol in self.Sym_PromoterTargetGenes:
                    print('The gene symbol %s is expressed by multiple promoters' % TargetGeneSymbol)
            self.Sym_PromoterTargetGenes.append(TargetGeneSymbol)
        if self.Switch4DebugDataset:
            print("# of Irregular Names: ", IrregularNames)
            # print("Names containing '-#' pattern:", NamesWithDashDigit)
        # elif Name[-1] == ('p' or 'P'):
        #         assert len(Name) <= 7, Name
        #         TargetGeneSymbol = Name[:-1]
        #     elif Name[-2] == ('p' or 'P'):
        #         assert len(Name) <= 8, Name
        #         if Name[-1].isdigit():
        #             TargetGeneSymbol = Name[:-2]
        #             print("Promoter '%s' contains an irregular naming style: %s" % (PromoterID, Name))

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FRNA(FDataset):
    def __init__(self):
        self.MolType_RNAs = 'RNA'
        self.MolType_RNAsNascent = 'RNA_Nascent'
        self.MolType_RNAsCleaved = 'RNA_Cleaved'

        self.ID2Count_RNAs_Temp = {}
        self.Name_RNAs = []
        self.Name2ID_RNAs = {}
        self.ID2Idx_RNAs = {}
        self.ID_RNAs = []
        self.HalfLife_RNAs = 0
        self.Seq_RNAs = []
        self.Len_RNAs = 0
        self.Type_RNAs = []
        self.MW_RNAs = 0
        self.Count_NTsInRNAs = []
        self.Freq_NTsInRNAs = []
        self.Count_RNAs = 0

        self.NUniq_RNAs = 0
        self.NUniq_mRNAs = 0
        self.NUniq_tRNAs = 0
        self.NUniq_rRNAs = 0
        self.NUniq_miscRNAs = 0

        self.ID2Type_RNA = {}

        self.Idx_AllRNAs = []
        self.Idx_mRNA = []
        self.Idx_tRNA = []
        self.Idx_rRNA = []
        self.Idx_miscRNA = []

        self.ID_RNAsNascent = []
        self.Count_RNAsNascent = 0
        self.MW_RNAsNascent = 0

        self.ID_RNAsCleaved = []
        self.Count_RNAsCleaved = 0
        self.MW_RNAsCleaved = 0

        super().__init__()

    def SetUpData(self, Dataset, MasterDataset = None):
        RNACounts = Dataset['Counts4AllRNAs_DL_19_totalIs_17_3_TranscriptElongation_BulkMolecules.tsv']
        for Value in RNACounts:
            ID, Count, Request = Value
            assert ID not in self.ID2Count_RNAs_Temp
            self.RemoveLocalizationFromMolName(ID)
            ID = self.RemoveLocalizationFromMolName(ID)
            self.ID2Count_RNAs_Temp[ID] = int(Count)

        RNAs = Dataset['rnas.tsv']
        self.NUniq_RNAs = len(RNAs)
        self.HalfLife_RNAs = np.zeros(self.NUniq_RNAs)
        self.Len_RNAs = np.zeros(self.NUniq_RNAs)
        self.MW_RNAs = np.zeros(self.NUniq_RNAs)
        self.Count_RNAs = np.zeros(self.NUniq_RNAs)

        for i, Value in enumerate(RNAs):
            HalfLife, Name, Seq, Type, ModifiedForms, MonomerID, Comments, MW, Location, NTCount, RNAID, GeneID, MicArrExp = Value
            assert Name not in self.ID2Idx_RNAs
            self.Name_RNAs.append(Name)
            self.Name2ID_RNAs[Name] = RNAID
            self.ID2Idx_RNAs[RNAID] = len(self.ID_RNAs)
            self.ID_RNAs.append(RNAID)
            self.HalfLife_RNAs[i] = HalfLife
            self.Seq_RNAs.append(Seq)
            self.Len_RNAs[i] = (len(Seq))
            self.Type_RNAs.append(Type)
            self.MW_RNAs[i] = max(np.fromstring(MW[1:-1], dtype='float32', sep=','))
            assert self.MW_RNAs[i] != 0
            self.Count_NTsInRNAs.append(NTCount)
            self.Count_RNAs[i] = self.ID2Count_RNAs_Temp[self.Name2ID_RNAs[Name]]
            self.Idx_AllRNAs.append(i)

            if Type == 'mRNA':
                self.NUniq_mRNAs += 1
                self.Idx_mRNA.append(i)
            elif Type == 'tRNA':
                self.NUniq_tRNAs += 1
                self.Idx_tRNA.append(i)
            elif Type == 'rRNA':
                self.NUniq_rRNAs += 1
                self.Idx_rRNA.append(i)
            elif Type == 'miscRNA':
                self.NUniq_miscRNAs += 1
                self.Idx_miscRNA.append(i)
            else:
                print('Warning: Unaccounted RNA type detected: ', Type)
            self.ID2Type_RNA[RNAID] = Type

            # RNANTFreq
            NTTotalCount = len(Seq)
            NTCounts = np.array([Seq.count("A"), Seq.count("C"), Seq.count("G"), Seq.count("U")])
            if np.sum(NTCounts) != NTTotalCount:
                print("WARNING: RNA seq may contain non-ACGU text.", file=sys.stderr)
            NTFreq = NTCounts / NTTotalCount
            self.Count_NTsInRNAs[i] = NTCounts
            self.Freq_NTsInRNAs = NTFreq

        self.ID_RNAsNascent = self.AppendStr(self.ID_RNAs, '_Nas')
        self.Count_RNAsNascent = np.zeros(self.NUniq_RNAs)
        self.MW_RNAsNascent = self.MW_RNAs / 2 # On average assuming transcription runs completely

        self.ID_RNAsCleaved = self.AppendStr(self.ID_RNAs, '_Cleaved')
        self.Count_RNAsCleaved = np.zeros(self.NUniq_RNAs)
        self.MW_RNAsCleaved = self.MW_RNAs

        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.ID_RNAs, self.MolType_RNAs, self.Count_RNAs, self.MW_RNAs)
        self.AddToMaster(MasterDataset, self.ID_RNAsNascent, self.MolType_RNAsNascent, self.Count_RNAsNascent, self.MW_RNAsNascent)
        self.AddToMaster(MasterDataset, self.ID_RNAsCleaved, self.MolType_RNAsCleaved, self.Count_RNAsCleaved, self.MW_RNAsCleaved)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FProtein(FDataset):
    def __init__(self):
        self.MolType_Proteins = 'Protein'
        self.MolType_Proteins_Nascent = 'Protein_Nascent'
        self.MolType_Proteins_Cleaved = 'Protein_Cleaved'

        self.ID2Count_Proteins_Temp = {}
        self.Name_Proteins = []
        self.Name2ID_Proteins = {}
        self.ID2Idx_Proteins = {}
        self.ID_Proteins = []
        self.Seq_Proteins = []
        self.Len_Proteins = 0
        self.MW_Proteins = 0
        self.Loc_Proteins = []
        self.Count_AAsInProteins = []
        self.Count_Proteins = 0
        self.Freq_AAsInProteins = []
        self.ID2ID_Protein2Gene = {}
        self._GeneIDs = []
        self.ID2ID_Protein2RNA = {}
        self._RNAIDs = []

        self.NUniq_Proteins = 0
        
        self.ID_ProteinsNascent = []
        self.Count_ProteinsNascent = 0
        self.MW_ProteinsNascent = 0

        self.ID_ProteinsCleaved = []
        self.Count_ProteinsCleaved = 0
        self.MW_ProteinsCleaved = 0

        super().__init__()

        return

    # Protein Monomers
    def SetUpData(self, Dataset, MasterDataset = None):
        ProteinCounts = Dataset['Counts4AllProteinMonomers_DL_28_totalIs_26_5_PolypeptideElongation_BulkMolecules.tsv']
        for Value in ProteinCounts:
            ID, Count, Request = Value
            assert ID not in self.ID2Count_Proteins_Temp
            # self.RemoveLocalizationFromMolID(ID)
            ID = self.RemoveLocalizationFromMolName(ID)
            self.ID2Count_Proteins_Temp[ID] = int(Count)

        Proteins = Dataset['proteins.tsv']
        self.NUniq_Proteins = len(Proteins)
        self.Len_Proteins = np.zeros(self.NUniq_Proteins)
        self.MW_Proteins = np.zeros(self.NUniq_Proteins)
        self.Count_Proteins = np.zeros(self.NUniq_Proteins)

        for i, Value in enumerate(Proteins):
            AACount, Name, Seq, Comments, CodingRNASeq, MW, Location, RNAID, ProtMonomerID, GeneID = Value
            self.Name_Proteins.append(Name)
            self.Name2ID_Proteins[Name] = ProtMonomerID
            self.ID2Idx_Proteins[ProtMonomerID] = len(self.ID_Proteins)
            self.ID_Proteins.append(ProtMonomerID)
            self.Seq_Proteins.append(Seq)
            self.Len_Proteins[i] = len(Seq)
            self.MW_Proteins[i] = max(np.fromstring(MW[1:-1], dtype='float32', sep=','))
            assert self.MW_Proteins[i] != 0
            self.Loc_Proteins.append(Location)
            AACount = ast.literal_eval(AACount)
            self.Count_AAsInProteins.append(AACount)
            self.Count_Proteins[i] = self.ID2Count_Proteins_Temp[ProtMonomerID]
            self.ID2ID_Protein2RNA[ProtMonomerID] = RNAID
            # ProteinAAFreq
            AATotalCount = len(Seq)
            if AATotalCount != self.Len_Proteins[i]:
                print("WARNING: Protein length data may contain error.", file=sys.stderr)
            AACounts = np.array(AACount)
            AAFreq = AACounts / AATotalCount
            self.Freq_AAsInProteins.append(AAFreq)

        self.ID_ProteinsNascent = self.AppendStr(self.ID_Proteins, '_Nas')
        self.Count_ProteinsNascent = np.zeros(self.NUniq_Proteins)
        self.MW_ProteinsNascent = self.MW_Proteins / 2 # On average assuming transcription runs completely

        self.ID_ProteinsCleaved = self.AppendStr(self.ID_Proteins, '_Cleaved')
        self.Count_ProteinsCleaved = np.zeros(self.NUniq_Proteins)
        self.MW_ProteinsCleaved = self.MW_Proteins

        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.ID_Proteins, self.MolType_Proteins, self.Count_Proteins, self.MW_Proteins)
        self.AddToMaster(MasterDataset, self.ID_ProteinsNascent, self.MolType_Proteins_Nascent, self.Count_ProteinsNascent, self.MW_ProteinsNascent)
        self.AddToMaster(MasterDataset, self.ID_ProteinsCleaved, self.MolType_Proteins_Cleaved, self.Count_ProteinsCleaved, self.MW_ProteinsCleaved)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FComplex(FDataset):
    def __init__(self):
        self.MolType_Complexes = 'Complex'

        self.ID_Proteins = []
        self.ID2Count_Proteins_Temp = {}

        self.ID2Count_Complexes_Temp = {}
        self.IDWithoutLoc2Count_Complexes_Temp = {}
        self.Name_Complexes = []
        self.Name2ID_Complexes = {}
        self.ID_Complexes = []
        self.ID2Idx_Complexes = {}
        self.Count_Complexes = []
        self.MW_Complexes = []
        self.Loc_Complexes = []

        self.NUniq_Complexes = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        Proteins = Dataset['proteins.tsv']
        for Value in Proteins:
            AACount, Name, Seq, Comments, CodingRNASeq, MW, Location, RNAID, ProtMonomerID, GeneID = Value
            self.ID_Proteins.append(ProtMonomerID)

        ProteinCounts = Dataset['Counts4AllProteinMonomers_DL_28_totalIs_26_5_PolypeptideElongation_BulkMolecules.tsv']
        for Value in ProteinCounts:
            ID, Count, Request = Value
            assert ID not in self.ID2Count_Proteins_Temp
            ID = self.RemoveLocalizationFromMolName(ID)
            self.ID2Count_Proteins_Temp[ID] = int(Count)

        ComplexCounts = Dataset['Counts4AllProteinMonomersAndCPLXsForComplexation_DL_44_totalIs_42_8_Complexation_BulkMolecules.tsv']
        for Value in ComplexCounts:
            ID, Count, Request = Value
            IDWithoutLocalization = self.RemoveLocalizationFromMolName(ID)
            assert ID not in self.ID2Count_Complexes_Temp
            assert IDWithoutLocalization not in self.IDWithoutLoc2Count_Complexes_Temp
            if (ID or IDWithoutLocalization) in self.ID_Proteins:
                continue
            # self.RemoveLocalizationFromMolName(ID)
            self.ID2Count_Complexes_Temp[ID] = int(Count)
            self.IDWithoutLoc2Count_Complexes_Temp[IDWithoutLocalization] = int(Count)

        Complexes = Dataset['proteinComplexes_large.tsv'] + Dataset['proteinComplexes.tsv']

        NoCountInfo = 0
        # Determine the number of unique complex entries from two combined datasets
        for Value in Complexes:
            Name, Comments, MW, Location, RxnID, ComplexID = Value
            if ComplexID in self.ID_Complexes:
                continue
            self.Name_Complexes.append(Name)
            self.Name2ID_Complexes[Name] = ComplexID
            self.ID2Idx_Complexes[ComplexID] = len(self.ID_Complexes)
            self.ID_Complexes.append(ComplexID)
            self.MW_Complexes.append(max(np.fromstring(MW[1:-1], dtype='float32', sep=',')))
            assert self.MW_Complexes[-1] != 0
            self.Loc_Complexes.append(Location)
            if ComplexID in self.IDWithoutLoc2Count_Complexes_Temp:
                Count = self.IDWithoutLoc2Count_Complexes_Temp[ComplexID]
                self.Count_Complexes.append(Count)
            else:
                if ComplexID not in self.ID2Count_Proteins_Temp:
                    NoCountInfo += 1
                    # print('%s is not found in either Protein or Complex IDs' % ComplexID)
                    self.Count_Complexes.append(0)

        print('# of Complexes missing Count info from wcm:', NoCountInfo)
        self.NUniq_Complexes = len(self.ID_Complexes)

        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.ID_Complexes, self.MolType_Complexes, self.Count_Complexes, self.MW_Complexes)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FReaction(FDataset):
    def __init__(self):
        self.MolType_RXNs = 'RXN'

        self.ID_RXNs = []
        self.ID2Idx_RXNs = {}
        self.Stoich_MolStoichDictInRXN = []
        self.ID_MolsInRXNs = []
        self.Stoich_MolsInRXNs = [] # Stoich for Stoichiometry
        self.Dir_RevRXNs = 0 # Rev for Reversibility
        self.ID_Enzs4RXNs = []
        # self.ID_SubstratesInRXNs = []
        # self.Stoich_SubstratesInRXNs = []
        # self.ID_ProductsInRXNs = []
        # self.Stoich_ProductsInRXNs = []

        self.NUniq_RXNs = 0

        super().__init__()

    # Reactions and kinetic parameters (K values)
    def SetUpData(self, Dataset, MasterDataset = None):
        RXNs = Dataset['reactions.tsv']

        self.NUniq_RXNs = len(RXNs)
        self.Dir_RevRXNs = np.zeros(self.NUniq_RXNs)

        ReversibilityBinaryDict = {}
        ReversibilityBinaryDict['true'] = 1
        ReversibilityBinaryDict['false'] = 0

        for i, Value in enumerate(RXNs):
            RXNID, RXNStoichiometry, RXNReversibility, RXNEnzymeID = Value
            assert RXNID not in self.ID2Idx_RXNs
            self.ID2Idx_RXNs[RXNID] = len(self.ID_RXNs) # = i
            self.ID_RXNs.append(RXNID)
            RXNReversibilityBinary = ReversibilityBinaryDict[RXNReversibility]
            self.Dir_RevRXNs[i] = RXNReversibilityBinary
            if RXNEnzymeID.strip('[').strip(']'):
                RXNEnzymeID = ast.literal_eval(RXNEnzymeID.strip('[').strip(']'))
                if type(RXNEnzymeID) == tuple:
                    for ID in RXNEnzymeID:
                        self.ID_Enzs4RXNs.append(ID)
                else:
                    self.ID_Enzs4RXNs.append(RXNEnzymeID)
            MolsInRXN = []
            StoichInRXN = []
            # Substrates = []
            # SubstrateStoichs = []
            # Products = []
            # ProductStoichs = []
            RXNStoichiometry = ast.literal_eval(RXNStoichiometry)
            for Mol, Stoich in RXNStoichiometry.items(): # RXNStoichiometry is a string dictionary
                MolsInRXN.append(Mol)
                StoichInRXN.append(Stoich)
            #     Mol = self.RemoveLocalizationFromMolName(Mol)
            #     if Stoich < 0:
            #         Substrates.append(Mol)
            #         SubstrateStoichs.append(-Stoich)
            #     if Stoich > 0:
            #         Products.append(Mol)
            #         ProductStoichs.append(Stoich)
            # assert SubstrateStoichs or ProductStoichs > 0, '%s: The exported Stoich value has to be greater than 0' % RXNID
            self.Stoich_MolStoichDictInRXN.append(RXNStoichiometry)
            self.ID_MolsInRXNs.append(MolsInRXN)
            self.Stoich_MolsInRXNs.append(StoichInRXN)
            # self.ID_SubstratesInRXNs.append([Substrates])
            # self.Stoich_SubstratesInRXNs.append([SubstrateStoichs])
            # self.ID_ProductsInRXNs.append([Products])
            # self.Stoich_ProductsInRXNs.append([ProductStoichs])

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FEnzyme(FDataset):
    def __init__(self): # MasterLocalizations
        self.MolType_Kinetics = 'Kinetics'

        self.KinRXNID2KineticData = {}
        self.KinRXNClassIDs = []
        self.KinEnzymeIDs = []
        self.KinSubstrateIDs = []
        self.Temp_EnzKinetics = 0
        self.Kcat_EnzKinetics = 0
        self.Km_EnzKinetics = []
        self.Ki_EnzKinetics = 0

        self.NUniq_EnzKinetics = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        # This kinetic table has no info on pH condition
        EnzymeKinetics = Dataset['enzymeKinetics.tsv']
        self.NUniq_EnzKinetics = len(EnzymeKinetics)
        self.Temp_EnzKinetics = np.zeros(self.NUniq_EnzKinetics)
        self.Kcat_EnzKinetics = np.zeros(self.NUniq_EnzKinetics)
        self.Ki_EnzKinetics = np.zeros(self.NUniq_EnzKinetics)
        N_TempAbscent = 0
        N_KcatValueAbscent = 0
        N_KmValueAbscent = 0
        N_KiValueAbscent = 0
        N_KmValueUnmatched = 0

        for i, Value in enumerate(EnzymeKinetics):
            PubmedID, Temp, Checked, Notes, Exclude, ConcentrationSubstrates, ReactionClassID, ReactionID, EnzymeIDs,\
            SubstrateIDs, RateEquationType, Direction, kcat, kM, kI, CustomRateEquation, CustomParameters,\
            CustomParameterConstants, CustomParameterConstantValues, CustomParameterVariables, CheckedBy = Value
            self.KinRXNID2KineticData[ReactionID] = len(self.KinEnzymeIDs)
            self.KinRXNClassIDs.append(ReactionClassID)
            self.KinEnzymeIDs.append(EnzymeIDs)
            # Exception handling: Substrates may be more than one. Same applies to Km values
            SubstrateIDs_Rebuilt = []
            for Substrate in SubstrateIDs.split(','):
                SubstrateIDs_Rebuilt.append(Substrate)
            self.KinSubstrateIDs.append(SubstrateIDs_Rebuilt)

            # Exception handling: Temp and K values are missing
            Temp = Temp.replace('[', '').replace(']', '')
            if Temp:
                self.Temp_EnzKinetics[i] = float(Temp)
            else:
                N_TempAbscent += 1
                if self.Switch4DebugDataset:
                    print('TempAbscent:', ReactionID)

            kcat = kcat.replace(',', '').replace('[', '').replace(']', '')
            if kcat:
                self.Kcat_EnzKinetics[i] = float(kcat)
            else:
                N_KcatValueAbscent += 1
                if self.Switch4DebugDataset:
                    print('KcatValueAbscent:', ReactionID)

            # In case there were more than one substrate
            kM = kM.replace('[', '').replace(']', '')
            N_SubstrateIDs = len(SubstrateIDs_Rebuilt)
            if N_SubstrateIDs > 1:
                if kM.count(',') >= 1:
                    kM_Rebuilt = np.zeros(N_SubstrateIDs)
                    for j, kM_value in enumerate(kM.split(',')):
                        kM_Rebuilt[j] = float(kM_value)
                    self.Km_EnzKinetics.append(kM_Rebuilt)
                else:
                    N_KmValueUnmatched += 1
                    if self.Switch4DebugDataset:
                        print(ReactionID, ': The number of kM values do not match the number of substrates')
                    self.Km_EnzKinetics.append(0.0)
            else:
                if kM:
                    self.Km_EnzKinetics.append(float(kM))
                else:
                    self.Km_EnzKinetics.append(0.0)
                    N_KmValueAbscent += 1
                    if self.Switch4DebugDataset:
                        print('KmValueAbscent:', ReactionID)

            kI = kI.replace(',', '').replace('[', '').replace(']', '')
            if kI:
                self.Ki_EnzKinetics[i] = float(kI)
            else:
                N_KiValueAbscent += 1
                if self.Switch4DebugDataset:
                    print('KiValueAbscent:', ReactionID)
        Km_ExportShort = self.NUniq_EnzKinetics - len(self.Km_EnzKinetics)
        assert len(self.Km_EnzKinetics) == self.NUniq_EnzKinetics, '%s Km exporting fell short' % Km_ExportShort
        if self.Switch4DebugDataset:
            print('N_Total:', self.NUniq_EnzKinetics)
            print('N_TempAbscent:', N_TempAbscent)
            print('N_KcatValueAbscent:', N_KcatValueAbscent)
            print('N_KmValueAbscent:', N_KmValueAbscent)
            print('N_KiValueAbscent:', N_KiValueAbscent)
            print('N_KmValueUnmatched:', N_KmValueUnmatched)

      # # Adjust kcat for temperature
        # temperature = constraint["Temp"]
        #
        # # If temperature not reported, assume 25 C
        # if type(temperature) == str:
        #     temperature = 25
        # constraint["kcatAdjusted"] = 2 ** ((37. - temperature) / 10.) * constraint["kcat"]

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FCompartment(FDataset):
    def __init__(self): # MasterLocalizations

        self.ID_Compartments = []
        self.ID2Key_Compartments = {}
        self.Key_Compartments = []
        self.Key2Idx_Compartments = {}
        self.NUniq_Compartments = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        Compartments = Dataset['compartments.tsv']
        self.NUniq_Compartments = len(Compartments)
        self.Idx_Compartments = np.zeros(self.NUniq_Compartments)

        for i, Value in enumerate(Compartments):
            Abbrev, ID = Value
            self.Key2Idx_Compartments[Abbrev] = len(self.Key_Compartments)
            self.Idx_Compartments[i] = len(self.Key_Compartments)
            self.Key_Compartments.append(Abbrev)
            self.ID_Compartments.append(ID)
            self.ID2Key_Compartments[ID] = Abbrev

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FBuildingBlock(FDataset):
    def __init__(self):
        # Lists of molecular specie
        self.Name_dNTPs = []
        self.Name_NTPs = []
        self.Name_AAs = []
        self.Key_dNTPs = []
        self.Key_NTPs = []
        self.Key_AAs = []
        self.NUniq_dNTPs = 0
        self.NUniq_NTPs = 0
        self.NUniq_AAs = 0

        super().__init__()

    def SetUpData(self, Dataset, MasterDataset = None):
        dNTPs = Dataset['dntps.txt']
        for BuildingBlock in dNTPs:
            BuildingBlock = self.AddLocalizationToMolName(BuildingBlock, '[c]')
            self.Name_dNTPs.append(BuildingBlock)
        self.Key_dNTPs = ['A', 'C', 'G', 'T']
        self.NUniq_dNTPs = len(self.Name_dNTPs)

        NTPs = Dataset['ntps.txt']
        for BuildingBlock in NTPs:
            BuildingBlock = self.AddLocalizationToMolName(BuildingBlock, '[c]')
            self.Name_NTPs.append(BuildingBlock)
        self.Key_NTPs = ['A', 'C', 'G', 'U']
        self.NUniq_NTPs = len(self.Name_NTPs)

        AAs = Dataset['amino_acids.txt']
        for BuildingBlock in AAs:
            BuildingBlock = self.AddLocalizationToMolName(BuildingBlock, '[c]')
            self.Name_AAs.append(BuildingBlock)
        self.NUniq_AAs = len(self.Name_AAs)

        AAKeys = Dataset['amino_acid_keys.txt']
        for BuildingBlock in AAKeys:
            self.Key_AAs.append(BuildingBlock)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)

# Master Dataset is an exception to the SetUpData method
class FMaster():
    def __init__(self):
        self.Switch4DebugMasterDataset = True

        self.ID2Idx_Master = {}
        self.Type_Master = []
        self.ID2Type_Master = {}
        self.ID_Master = []
        self.Count_Master = []
        self.MW_Master = []
        self.NUniq_Master = 0

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)
