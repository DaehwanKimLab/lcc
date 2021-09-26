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
    def SetUpData(self, Dataset, Comp):
        pass

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)

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
            # assert ID not in MasterDataset.ID_Master, 'Redundant ID in MasterDataset.ID_Master: %s' % ID
            # assert ID not in MasterDataset.ID2Idx_Master, 'Redundant ID in MasterDataset.ID2Idx_Master: %s' % ID
            # assert ID not in MasterDataset.ID2Type_Master, 'Redundant ID in MasterDataset.ID2Type_Master: %s' % ID
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
                    print('Master Dataset | ', 'addition completed: Class "%s"' % Key[3:])

                if Key[:6] == 'Master':
                    LengthOfValue = len(Value)
                    print('Master Dataset | ', 'The length of "%s": %d' % (Key, LengthOfValue))

        MasterDataset.NUniq_Master = len(MasterDataset.ID_Master)

        assert len(MasterDataset.ID_Master) == len(MasterDataset.ID2Idx_Master) and len(MasterDataset.ID_Master) == len(MasterDataset.Count_Master) and len(MasterDataset.ID_Master) == len(MasterDataset.MW_Master)

    def AppendStr(self, List, Str):
        # Appends Str to the List
        ListStr = list()
        for i in List:
            i += Str
            ListStr.append(i)
        return ListStr

    def SetUpMetaboliteIDRef(self, Dataset):
        # To handle unregistered metabolites
        self.ID_MetabolitesRef = list()

        Metabolites = Dataset['Counts4AllMetabolites_DL_2_totalIs_0_0_Metabolism_BulkMolecules.tsv']
        for Metabolite in Metabolites:
            self.ID_MetabolitesRef.append(Metabolite[0][:-3])

    def CheckIDForModification(self, ID, IDModificationRef):
        if ID in IDModificationRef:
            return IDModificationRef[ID]
        else:
            # if ID.find('+2') >= 0:
            #     ID += '[c]'
            if ID in self.ID_MetabolitesRef:
                ID += '[c]'
            return ID

class FMetabolite(FDataset):
    def __init__(self):
        self.MolType_Metabolites = 'Metabolite'
        self.ID2MW_Metabolites_Temp = dict()
        self.ID_Metabolites = list()
        self.ID2Idx_Metabolites = dict()
        self.Count_Metabolites = 0
        self.MW_Metabolites = 0
        self.NUniq_Metabolites = 0
        self.SpecieIdx_Metabolites = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, Comp):
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
            self.MW_Metabolites[i] = self.ID2MW_Metabolites_Temp[NameLocRemoved]
            self.Count_Metabolites[i] = Count
            assert self.MW_Metabolites[i] != 0

            # Temporary correction for water count(out of int32 range)
            if (int(Count) < 0) or (int(Count) > 2147483647):
                Count_Temporary_Int32 = int(1e7)
                self.Count_Metabolites[i] = Count_Temporary_Int32
                if self.Switch4DebugDataset:
                    print("Metabolites | The Count of %s is %s, which is out of int32 range (> 2147483647). It has been temporarily replaced with the value of %s" % (Name, Count, Count_Temporary_Int32))

        # Add to the Master Dataset
        self.AddToMaster(Comp.Master, self.ID_Metabolites, self.MolType_Metabolites, self.Count_Metabolites, self.MW_Metabolites)


class FChromosome(FDataset):
    def __init__(self):
        self.MolType_Chromosomes = 'Chromosome'
        self.Len_ChromosomesInGenome = 0
        self.Count_NTsInChromosomesInGenome = list()
        self.Freq_NTsInChromosomesInGenome = list()
        self.Genome = None

        # self.N_InitialCopyPerChromosome = 0
        # self.N_FullChromosomesPerReplication = 0
        # self.N_PartialChromosomesPerReplication = 0
        self.N_PossibleRoundsOfConcurrentReplications = 0
        # self.ReplicatingChromosomeAttributes = list()

        # self.N_FullChromosomes = 0
        # self.N_PartialChromosomes = 0
        # self.N_AllChromosomes = 0

        self.ID_Chromosomes = list()
        self.Count_Chromosomes = 0
        self.MW_DNABasePairs = 0
        self.Freq_NTsInChromosomes = list()
        self.Freq_NTsInChromosomesReplicating = list()
        self.NUniq_ChromosomesInGenome = 0
        self.NMax_Chromosomes = 0

        super().__init__()

    def SetUpData(self, Dataset, Comp):
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
                print("Genome Dataset | ', 'WARNING: RNA seq may contain non-ACGT text.", file=sys.stderr)
            NTFreq = NTCounts / self.Len_ChromosomesInGenome[i]
            self.Count_NTsInChromosomesInGenome.append(NTCounts)
            self.Freq_NTsInChromosomesInGenome.append(NTFreq)

        # Set up an arbitrary max number to keep NT Count arrays for all chromosomes.
        self.N_PossibleRoundsOfConcurrentReplications = 2  # n
        RepFactor = 2 ** self.N_PossibleRoundsOfConcurrentReplications
        self.NMax_Chromosomes = self.NUniq_ChromosomesInGenome * RepFactor

        self.Count_Chromosomes = np.zeros(self.NMax_Chromosomes)
        self.MW_DNABasePairs = np.zeros(self.NMax_Chromosomes)
        DNABasePairMW = 650 # average molecular weight of one base pair: 650 Daltons, which is 650 g/mol

        ChromosomeCount = 0
        for ReplicationRoundID in range(RepFactor):

            for ChromosomeNumber in range(self.NUniq_ChromosomesInGenome):
                if ReplicationRoundID == 0:
                    ChromosomeID = 'Ch%d_Original' % (ChromosomeNumber)
                    self.ID_Chromosomes.append(ChromosomeID)
                    self.Count_Chromosomes[ChromosomeCount] = 1
                    self.MW_DNABasePairs[ChromosomeCount] = self.Len_ChromosomesInGenome[
                                                                ChromosomeNumber] * DNABasePairMW
                    self.Freq_NTsInChromosomes.append(self.Freq_NTsInChromosomesInGenome[ChromosomeNumber])
                    ChromosomeCount += 1
                    continue

                else:
                    ChromosomeID = 'Ch%d_Replicating_Round%d' % (ChromosomeNumber, ReplicationRoundID)
                    self.ID_Chromosomes.append(ChromosomeID)
                    self.MW_DNABasePairs[ChromosomeCount] = self.Len_ChromosomesInGenome[ChromosomeNumber] * DNABasePairMW / 2  # average length
                    self.Freq_NTsInChromosomes.append(self.Freq_NTsInChromosomesInGenome[ChromosomeNumber])
                    self.Freq_NTsInChromosomesReplicating.append(self.Freq_NTsInChromosomesInGenome[ChromosomeNumber])

                    ChromosomeCount += 1

                # ChromosomeState = 'Replicating'
                # ChromosomeID = 'Ch%d_' % (ChromosomeNumber, ChromosomeRepNumber, ChromosomeState)
                # self.ID_Chromosomes.append(ChromosomeID)
                # self.Freq_NTsInChromosomes.append(self.Freq_NTsInChromosomesInGenome[i])  # Just a filler for the list
                # ChromosomeCount += 1
                #
                # for Attribute in self.ReplicatingChromosomeAttributes:
                #     ChromosomeState = 'Partial_' + Attribute
                #     ChromosomeID = 'Ch%d_Rep%d_%s' % (ChromosomeNumber, ChromosomeRepNumber, ChromosomeState)
                #     self.ID_Chromosomes.append(ChromosomeID)
                #     self.MW_DNABasePairs[ChromosomeCount] = DNABasePairMW
                #     self.Freq_NTsInChromosomes.append(self.Freq_NTsInChromosomesInGenome[i])
                #     ChromosomeCount += 1

        # Add to the Master Dataset
        self.AddToMaster(Comp.Master, self.ID_Chromosomes, self.MolType_Chromosomes, self.Count_Chromosomes, self.MW_DNABasePairs)


class FGene(FDataset):
    def __init__(self):
        self.MolType_Genes = 'Gene'

        self.Len_Genes = 0
        self.Coord_Genes = 0 # Coord for Coordinate
        self.MW_Genes = 0
        self.Dir_Genes = 0 # Dir for Direction
        self.Name_Genes = list()
        self.Sym_Genes = list() # Sym for Symbols
        self.ID_Genes = list()
        self.Seq_Genes = list()
        self.Type_Genes = list()
        self.ID2Idx_Genes = dict()
        self.Sym2Idx_Genes = dict()
        self.Sym2ID_Genes = dict()
        self.Name2ID_Genes = dict()
        self.Count_Genes = 0
        self.NUniq_Genes = 0

        self.Coord_Genes_Reindexed = 0
        self.Coord_Genes_Reindexed_Leftward = 0
        self.Coord_Genes_Reindexed_Rightward = 0


        super().__init__()

    def SetUpData(self, Dataset, Comp):
        Genes = Dataset['genes.tsv']
        self.NUniq_Genes = len(Genes)
        self.Len_Genes = np.zeros(self.NUniq_Genes)
        self.Coord_Genes = np.zeros(self.NUniq_Genes)
        self.MW_Genes = np.zeros(self.NUniq_Genes) # No molecular weights to be calculated
        self.Dir_Genes = np.zeros(self.NUniq_Genes)
        self.Count_Genes = np.ones(self.NUniq_Genes) # The count may be adjusted according to the starting replication state of the chromosomes.
        self.Name_Genes = list()
        self.Sym_Genes = list()
        self.ID_Genes = list()
        self.Type_Genes = list()
        self.Name2ID_Genes = dict()
        self.ID2Idx_Genes = dict()
        self.Sym2Idx_Genes = dict()

        DirectionBinaryDict = dict()
        DirectionBinaryDict['+'] = 1
        DirectionBinaryDict['-'] = 0

        for i, Value in enumerate(Genes):
            Length, Name, Seq, RNAID, Coordinate, Direction, Symbol, Type, GeneID, MonomerID = Value
            # assert Name not in self.Name2ID_Genes, '"%s" has a redundant name in GeneName2ID dictionary' % Name
            # self.Name_Genes.append(Name)
            # self.Name2ID_Genes[Name] = GeneID
            self.Sym_Genes.append(Symbol)
            self.Sym2ID_Genes[Symbol] = GeneID
            self.Sym2Idx_Genes[Symbol] = len(self.ID_Genes)
            assert GeneID not in self.ID2Idx_Genes
            self.ID_Genes.append(GeneID)
            self.ID2Idx_Genes[GeneID] = len(self.Seq_Genes)
            self.Seq_Genes.append(Seq)
            self.Len_Genes[i] = (len(Seq))
            self.Coord_Genes[i] = int(Coordinate)
            DirectionBinary = DirectionBinaryDict[Direction]
            self.Dir_Genes[i] = DirectionBinary
            self.Type_Genes.append(Type)

        # Add to the Master Dataset
        self.AddToMaster(Comp.Master, self.ID_Genes, self.MolType_Genes, self.Count_Genes, self.MW_Genes)

        # Reindexing coordinates according to the center of OriC as 0
        Genome = Dataset['EscherichiaColi.fasta']
        Len_Genome = len(Genome['Ch1'])

        Coord_Ori_Start = 3925744
        Len_Ori = 232
        Coord_Ori = Coord_Ori_Start + int(Len_Ori / 2)

        self.Coord_Genes_Reindexed = np.zeros(self.NUniq_Genes)
        for i, Coord in enumerate(self.Coord_Genes):
            if (Coord >= Coord_Ori) and (Coord <= Len_Genome):
                self.Coord_Genes_Reindexed[i] = Coord - Coord_Ori
            elif (Coord >= 0) and (Coord <= Coord_Ori):
                self.Coord_Genes_Reindexed[i] = Coord + (Len_Genome - Coord_Ori)
            else:
                print("Genome Dataset | ', 'WARNING: Unidentifiable coordinate range", file=sys.stderr)
        assert max(self.Coord_Genes_Reindexed) < Len_Genome, 'Reindexing failed'

        # Buffering and organizing coordinates for bidirectional scanning during replication

        Len_Buffer4Promoter = 50

        self.Coord_Genes_Reindexed_Leftward = np.zeros(self.NUniq_Genes)
        self.Coord_Genes_Reindexed_Rightward = np.zeros(self.NUniq_Genes)
        for i, Coord in enumerate(self.Coord_Genes_Reindexed):
            if self.Dir_Genes[i]:
                self.Coord_Genes_Reindexed_Leftward[i] = (Len_Genome - Coord) + Len_Buffer4Promoter
                self.Coord_Genes_Reindexed_Rightward[i] = Coord
            else:
                self.Coord_Genes_Reindexed_Leftward[i] = (Len_Genome - Coord)
                self.Coord_Genes_Reindexed_Rightward[i] = Coord + Len_Buffer4Promoter

        assert max(self.Coord_Genes_Reindexed_Leftward) < Len_Genome, 'Reindexing failed'
        assert max(self.Coord_Genes_Reindexed_Rightward) < Len_Genome, 'Reindexing failed'


class FPromoter(FDataset):
    def __init__(self):
        self.MolType_Promoters = 'Promoter'
        self.Coord_Promoters = 0
        self.Dir_Promoters = 0
        self.ID_Promoters = list()
        self.ID2Idx_Promoters = dict()
        self.Sym_PromoterTargetGenes = list()
        self.NUniq_Promoters = 0

        super().__init__()

    def SetUpData(self, Dataset, Comp):
        Promoters = Dataset['promoters.tsv']

        self.NUniq_Promoters = len(Promoters)
        self.Coord_Promoters = np.zeros(self.NUniq_Promoters)
        self.Dir_Promoters = np.zeros(self.NUniq_Promoters)
        IrregularNames = 0

        DirectionBinaryDict = dict()
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
                    print("Promoter Dataset | ', 'Promoter '%s' contains '-' in XXXp promoter name format: %s" % (PromoterID, Name))
                if len(Name.split('-')) > 2:
                    print("Promoter Dataset | ', 'Promoter '%s' contains more than two '-' in the promoter name: %s" % (PromoterID, Name))

            TargetGeneSymbol = None
            TargetGeneSymbolsOnly = list()
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
                #     TargetGeneSymbol = list()
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
                    print("Promoter Dataset | ', 'Promoter %s contains an irregular name format to parse: %s" % (PromoterID, Name))
                IrregularNames += 1
                TargetGeneSymbol = Name
            # Check if the gene is expressed by multiple promoters
            if self.Switch4DebugDataset:
                if TargetGeneSymbol in self.Sym_PromoterTargetGenes:
                    print('Promoter Dataset | ', 'The gene symbol %s is expressed by multiple promoters' % TargetGeneSymbol)
            self.Sym_PromoterTargetGenes.append(TargetGeneSymbol)
        if self.Switch4DebugDataset:
            print("Promoter Dataset | ', '# of Irregular Names: ", IrregularNames)
            # print("Names containing '-#' pattern:", NamesWithDashDigit)
        # elif Name[-1] == ('p' or 'P'):
        #         assert len(Name) <= 7, Name
        #         TargetGeneSymbol = Name[:-1]
        #     elif Name[-2] == ('p' or 'P'):
        #         assert len(Name) <= 8, Name
        #         if Name[-1].isdigit():
        #             TargetGeneSymbol = Name[:-2]
        #             print("Promoter '%s' contains an irregular naming style: %s" % (PromoterID, Name))


class FRNA(FDataset):
    def __init__(self):
        self.MolType_RNAs = 'RNA'

        self.ID2Count_RNAs_Temp = dict()
        self.Name_RNAs = list()
        self.Name2ID_RNAs = dict()
        self.ID2Idx_RNAs = dict()
        self.ID_RNAs = list()
        self.HalfLife_RNAs = 0
        self.Seq_RNAs = list()
        self.Len_RNAs = 0
        self.Type_RNAs = list()
        self.MW_RNAs = 0
        self.Count_NTsInRNAs = list()
        self.Freq_NTsInRNAs = list()
        self.Count_RNAs = 0

        self.NUniq_RNAs = 0
        self.NUniq_mRNAs = 0
        self.NUniq_tRNAs = 0
        self.NUniq_rRNAs = 0
        self.NUniq_miscRNAs = 0

        self.ID2Type_RNA = dict()

        self.Idx_AllRNAs = list()
        self.Idx_mRNA = list()
        self.Idx_tRNA = list()
        self.Idx_rRNA = list()
        self.Idx_miscRNA = list()

        self.ID2ID_Gene2RNA = dict()

        super().__init__()

    def SetUpData(self, Dataset, Comp):
        RNACounts = Dataset['Counts4AllRNAs_DL_19_totalIs_17_3_TranscriptElongation_BulkMolecules.tsv']
        for Value in RNACounts:
            ID, Count, Request = Value
            assert ID not in self.ID2Count_RNAs_Temp
            self.RemoveLocalizationFromMolName(ID)
            ID = self.RemoveLocalizationFromMolName(ID)
            self.ID2Count_RNAs_Temp[ID] = int(Count)

        RNAs = Dataset['rnas_Sorted_DL.tsv']  # sorted by GeneID
        # RNAs = Dataset['rnas.tsv']
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
            self.ID2ID_Gene2RNA[GeneID] = RNAID

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
                print('RNA Dataset | ', 'Warning: Unaccounted RNA type detected: ', Type)
            self.ID2Type_RNA[RNAID] = Type

            # RNANTFreq
            NTTotalCount = len(Seq)
            NTCounts = np.array([Seq.count("A"), Seq.count("C"), Seq.count("G"), Seq.count("U")])
            if np.sum(NTCounts) != NTTotalCount:
                print("RNA Dataset | ', 'WARNING: RNA seq may contain non-ACGU text.", file=sys.stderr)
            NTFreq = NTCounts / NTTotalCount
            self.Count_NTsInRNAs[i] = NTCounts
            self.Freq_NTsInRNAs.append(NTFreq)

        # Add to the Master Dataset
        self.AddToMaster(Comp.Master, self.ID_RNAs, self.MolType_RNAs, self.Count_RNAs, self.MW_RNAs)

        # Add a single variable to the Master Dataset
        Comp.Master.ID2ID_Gene2RNA_Master = self.ID2ID_Gene2RNA


class FProtein(FDataset):
    def __init__(self):
        self.MolType_Proteins = 'Protein'

        self.ID2Count_Proteins_Temp = dict()
        self.Name_Proteins = list()
        self.Name2ID_Proteins = dict()
        self.ID2Idx_Proteins = dict()
        self.ID_Proteins = list()
        self.Seq_Proteins = list()
        self.Len_Proteins = 0
        self.MW_Proteins = 0
        self.Loc_Proteins = list()
        self.Count_AAsInProteins = list()
        self.Count_Proteins = 0
        self.Freq_AAsInProteins = list()

        self.ID2ID_Protein2Gene = dict()
        self.ID2ID_Gene2Protein = dict()
        self.ID2ID_Protein2RNA = dict()
        self.ID2ID_RNA2Protein = dict()

        self.NUniq_Proteins = 0

        self.ID2ID_mRNA2Protein = dict()

        super().__init__()

        return

    # Protein Monomers
    def SetUpData(self, Dataset, Comp):
        ProteinCounts = Dataset['Counts4AllProteinMonomers_DL_28_totalIs_26_5_PolypeptideElongation_BulkMolecules.tsv']
        for Value in ProteinCounts:
            ID, Count, Request = Value
            assert ID not in self.ID2Count_Proteins_Temp
            # self.RemoveLocalizationFromMolID(ID)
            ID = self.RemoveLocalizationFromMolName(ID)
            self.ID2Count_Proteins_Temp[ID] = int(Count)

        Proteins = Dataset['proteins_Sorted_DL.tsv']  # Sorted by GeneID
        # Proteins = Dataset['proteins.tsv']
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
            self.ID2ID_Protein2Gene[ProtMonomerID] = GeneID
            self.ID2ID_Gene2Protein[GeneID] = ProtMonomerID
            self.ID2ID_Protein2RNA[ProtMonomerID] = RNAID
            self.ID2ID_RNA2Protein[RNAID] = ProtMonomerID
            self.ID2ID_mRNA2Protein[RNAID] = ProtMonomerID
            # ProteinAAFreq
            AATotalCount = len(Seq)
            if AATotalCount != self.Len_Proteins[i]:
                print("Protein Dataset | ', 'WARNING: Protein length data may contain error.", file=sys.stderr)
            AACounts = np.array(AACount)
            AAFreq = AACounts / AATotalCount
            self.Freq_AAsInProteins.append(AAFreq)

        # Add to the Master Dataset
        self.AddToMaster(Comp.Master, self.ID_Proteins, self.MolType_Proteins, self.Count_Proteins, self.MW_Proteins)

        # Add a single variable to the Master Dataset
        Comp.Master.ID2ID_mRNA2Protein_Master = self.ID2ID_mRNA2Protein

        # # Search for protein names containing 'protease' or 'peptidase'
        # for i, Value in enumerate(Proteins):
        #     AACount, Name, Seq, Comments, CodingRNASeq, MW, Location, RNAID, ProtMonomerID, GeneID = Value
        #     if 'endoribonuclease' in Name:
        #         print('GeneID: ',GeneID,'\tProtID: ',ProtMonomerID, '\tName: ', Name)
        #     if 'RNA helicase' in Name:
        #         print('GeneID: ',GeneID,'\tProtID: ',ProtMonomerID, '\tName: ', Name)
        #     if 'EG10844-MONOMER' in ProtMonomerID:
        #         print('GeneID: ',GeneID,'\tProtID: ',ProtMonomerID, '\tName: ', Name)


class FComplex(FDataset):
    def __init__(self):
        self.MolType_Complexes = 'Complex'

        self.ID_Proteins = list()
        self.ID2Count_Proteins_Temp = dict()

        self.ID2Count_Complexes_Temp = dict()
        self.IDWithoutLoc2Count_Complexes_Temp = dict()
        self.Name_Complexes = list()
        self.Name2ID_Complexes = dict()
        self.ID_Complexes = list()
        self.ID2Idx_Complexes = dict()
        self.Count_Complexes = list()
        self.MW_Complexes = list()
        self.Loc_Complexes = list()

        self.NUniq_Complexes = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, Comp):
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

        if self.Switch4DebugDataset:
            print('Complexes Dataset | ', '# of Complexes missing Count info from wcm:', NoCountInfo)
        self.NUniq_Complexes = len(self.ID_Complexes)

        # Add to the Master Dataset
        self.AddToMaster(Comp.Master, self.ID_Complexes, self.MolType_Complexes, self.Count_Complexes, self.MW_Complexes)


class FComplexation(FDataset):
    def __init__(self):
        self.MolType_Complexations = 'Complexation'

        self.ID_Complexations = list()
        self.ID2Idx_Complexations = dict()
        self.Rate_Complexations = list()
        self.Dir_RevComplexations = 0  # Rev for Reversibility

        self.ID_MolsInCPLXRXN = list()
        self.ID2Idx_MolsInCPLXRXN = dict()
        self.Coeff_MolsInCPLXRXN = list()
        self.NUniq_MolsInCPLXRXN = 0

        self.NUniq_Complexations = 0

        self.IDModification = dict()

        self.ID_UnregisteredFromCPLXRXN = list()   # Mostly additional metabolites
        self.MolType_UnregisteredFromCPLXRXN = 'MetabolitesUnregisteredFromCPLXRXN'
        self.Count_UnregisteredFromCPLXRXN = 0
        self.MW_UnregisteredFromCPLXRXN = 0
        self.NUniq_UnregisteredFromCPLXRXN = 0

        super().__init__()

    # Reactions and kinetic parameters (K values)
    def setupdata(self, Dataset, Comp):
        self.SetUpMetaboliteIDRef(Dataset)

        Complexations = Dataset['complexationReactions_large.tsv']

        # Remove molecules not listed in other flat data (mostly copied and pasted from wcs)

        NONSPECIFIC_METABOLITES = [
            "ALLOSE",
            "N-ACETYL-D-GLUCOSAMINE-6-P",
            'RIBOSE',
            'ACYL-COA',
            'N-ACETYLNEURAMINATE',
            'FRUCTURONATE',
            'ALLANTOIN',  # CPLX0-8071 (also glyxoylate)
        ]

        UNIDENTIFIED_GENE = [
            "TRANSENOYLCOARED-MONOMER",

            "GLUTAMINA-MONOMER",
            'GLUTAMINA-CPLX',

            "GLUTAMINB-MONOMER",
            'GLUTAMINB-CPLX',

            'TRANSENOYLCOARED-MONOMER',
            'TRANSENOYLCOARED-CPLX',

            'NQOR-MONOMER',
            'NQOR-CPLX',
        ]

        MISC_OR_UNEXPLAINED = [
            "ACETYL-COA-CARBOXYLMULTI-CPLX",  # biotinylated, many subunits
            'GCVMULTI-CPLX',  # many subunits
            'CPLX0-7754',  # modified form for one subunit (phos'd)
            'CPLX0-7849',  # proenzymes and modified forms, many forms
            'EG10245-MONOMER',  # DNA poly III
            'modified-charged-selC-tRNA',  #

            # biotin carrier protein
            'BCCP-BIOTIN',
            'BCCP-CPLX'

            # FadR plus CoA, met id not recognized
            "MONOMER-51-CPD-18346",
            'MONOMER-51-CPD-10269',

            # enterobactin multicomplex (also see modified forms)
            'ENTB-CPLX',
            'HOL-ENTB',
            'ENTMULTI-CPLX',

            # no idea
            'CPD-10269',
            'CPD-18346',

            # RcsB phosphorylated transcription factor, many forms
            'PHOSPHO-RCSB',
            'CPLX0-7884',
            'CPLX0-7978',

            # phos'd TF
            'CPLX0-7795',
            'PROTEIN-NRIP',
            'PHOSPHO-UHPA',
            'PHOSPHO-KDPE',
            'CPLX0-7721',
            'MONOMER0-4198',
            'CPLX0-7748',
            'PHOSPHO-CPXR',

            # weird orphan metabolite, might actually be a protein
            'CPD0-2342',

            # modified form?
            'LIPOYL-GCVH',
        ]

        DISABLED = [  # complexes we don't form for modeling reasons
            "CPLX0-3964",  # full ribosome

            # 50S subcomplex (can't procede to full complex)
            "CPLX0-3956",
            "CPLX0-3955",

            # non-apo RNA polymerase complexes
            "CPLX0-221",
            "RNAP54-CPLX",
            "RNAPS-CPLX",
            "RNAP32-CPLX",
            "RNAPE-CPLX",
            "CPLX0-222",
            "RNAP70-CPLX",
        ]

        SUBUNIT_PROENZYME = [
            'CPLX0-7885',
            # 'MONOMER0-4195': 'EG10374-MONOMER',
            # 'MONOMER0-4196': 'EG10374-MONOMER',
            #
            # 'SAMDC-ALPHA-MONOMER': 'SPED-MONOMER',
            # 'SAMDC-BETA-MONOMER': 'SPED-MONOMER',

            'PHOSPHASERDECARB-DIMER',
            'PHOSPHASERDECARB-CPLX',
            # 'PHOSPHASERDECARB-ALPHA-MONOMER': 'PSD-MONOMER',
            # 'PHOSPHASERDECARB-BETA-MONOMER': 'PSD-MONOMER',

            'CPLX0-263',
            # 'MONOMER0-2': 'EG12407-MONOMER',
            # 'MONOMER0-3': 'EG12407-MONOMER',

            'CPLX0-2901',
            # 'MONOMER0-1842': 'ASPDECARBOX-MONOMER',
            # 'MONOMER0-1843': 'ASPDECARBOX-MONOMER',
        ]

        DL_Addition = [
            'PHOSPHO-OMPR-MONOMER',  # Unable to map
            'SAMDECARB-CPLX',

            # # Unimplemented metabolites in the system
            # 'CU+',
            # 'CPD-207',
            # 'CPD0-881',
            # 'PUTRESCINE',
            # 'CPD-15818',
            # 'CPD0-1110',
            # 'CPD0-1108',
            # 'CPD-10330',
            # 'CPD-12537',
            # 'HYPOXANTHINE',
            # 'FUCULOSE-1P',
            # 'ANTIMONITE',
            # 'GLYCEROL',   # strangely not implemented
            # '2-KETO-3-DEOXY-6-P-GLUCONATE',
            # 'FRU1P',
            # 'BIO-5-AMP',
            # 'ALLOLACTOSE',
            # 'CPD-622',
            # 'CAMP',
            # 'XANTHOSINE',
            # 'NITRIC-OXIDE',   # strangely not implemented
            # 'GALACTOSE',   # strangely not implemented
            # 'NA+',   # strangely not implemented
            # 'CPD-3',

        ]

        BlackList = NONSPECIFIC_METABOLITES + UNIDENTIFIED_GENE + MISC_OR_UNEXPLAINED + DISABLED + SUBUNIT_PROENZYME + DL_Addition

        self.IDModification = {
            'ENTF-PANT': 'ENTF-MONOMER',
            'PHOSPHO-OMPR': 'OMPR-MONOMER',
            'G7678-MONOMER': 'MONOMER0-2341',
            'EG10245-MONOMER': 'MONOMER0-2383',   # DNA polymerase III, &gamma; subunit
            # 'MONOMER0-4195': 'EG10374-MONOMER',   # CPLX0-7885
            # 'MONOMER0-4196': 'EG10374-MONOMER',   # CPLX0-7885
            # 'SAMDC-ALPHA-MONOMER': 'SPED-MONOMER',
            # 'SAMDC-BETA-MONOMER': 'SPED-MONOMER',
            # 'PHOSPHASERDECARB-ALPHA-MONOMER': 'PSD-MONOMER',   # PHOSPHASERDECARB-DIMER, PHOSPHASERDECARB-CPLX
            # 'PHOSPHASERDECARB-BETA-MONOMER': 'PSD-MONOMER',   # PHOSPHASERDECARB-DIMER, PHOSPHASERDECARB-CPLX
            # 'MONOMER0-2': 'EG12407-MONOMER',   # CPLX0-263
            # 'MONOMER0-3': 'EG12407-MONOMER',   # CPLX0-263
            # 'MONOMER0-1842': 'ASPDECARBOX-MONOMER',   # CPLX0-2901
            # 'MONOMER0-1843': 'ASPDECARBOX-MONOMER',   # CPLX0-2901
        }

        Complexations_Filtered = list()
        Complexations_Removed = list()
        N_Filtered = 0
        for Value in Complexations:
            Filter = True
            for Mol in BlackList:
                if Value[1].find(Mol) >= 0:
                    Filter = False
                    N_Filtered += 1
                    Complexations_Removed.append(Value)
                    break
            if Filter:
                Complexations_Filtered.append(Value)
        assert len(Complexations) == (len(Complexations_Filtered) + len(Complexations_Removed))

        Complexations = Complexations_Filtered

        self.NUniq_Complexations = len(Complexations)
        self.Dir_RevComplexations = np.zeros(self.NUniq_Complexations)

        # Build indices for stoichiometry matrix
        N_Mol_Total = 0
        N_Mol_Repeated = 0
        Data_NotInMasterDataset_NotAddedToUnregistered_Temp = list()
        for Value in Complexations:
            ComplexationProcess, ComplexationStoichiometry, ComplexationID, ComplexationDirection = Value
            for MoleculeInfo in ComplexationStoichiometry.strip().strip('"').strip('[').strip(']').strip('{').split('},'):
                for LabelDataPair in MoleculeInfo.split(','):
                    Label, Data = LabelDataPair.strip().split(': ')
                    if Label == '"molecule"':
                        N_Mol_Total += 1
                        Data = self.CheckIDForModification(Data.strip('"'), self.IDModification)

                        # Check for unregistered molecules
                        if Data not in Comp.Master.ID_Master:
                            Data = self.AddLocalizationToMolName(Data, '[c]')
                            if Data not in Comp.Master.ID_Master:
                                if Data not in self.ID_UnregisteredFromCPLXRXN:
                                    self.ID_UnregisteredFromCPLXRXN.append(Data)
                                else:
                                    # Accumulate the original data without localization tag
                                    Data_NotInMasterDataset_NotAddedToUnregistered_Temp.append(Data)

                        # Data at this point has localization tag
                        if Data not in self.ID_MolsInCPLXRXN:
                            self.ID2Idx_MolsInCPLXRXN[Data] = len(self.ID_MolsInCPLXRXN)
                            self.ID_MolsInCPLXRXN.append(Data)
                        else:
                            N_Mol_Repeated += 1

        self.NUniq_UnregisteredFromCPLXRXN = len(self.ID_UnregisteredFromCPLXRXN)
        self.Count_UnregisteredFromCPLXRXN = np.zeros(self.NUniq_UnregisteredFromCPLXRXN)
        self.MW_UnregisteredFromCPLXRXN = np.zeros(self.NUniq_UnregisteredFromCPLXRXN)   # TODO: get values from database

        self.AddToMaster(Comp.Master, self.ID_UnregisteredFromCPLXRXN, self.MolType_UnregisteredFromCPLXRXN, self.Count_UnregisteredFromCPLXRXN, self.MW_UnregisteredFromCPLXRXN)

        self.NMax_ComplexationMolParts = len(self.ID_MolsInCPLXRXN)

        for i, Value in enumerate(Complexations):
            ComplexationProcess, ComplexationStoichiometry, ComplexationID, ComplexationDirection = Value
            assert ComplexationID not in self.ID_Complexations, '%s' % ComplexationID
            self.ID_Complexations.append(ComplexationID)
            self.ID2Idx_Complexations[ComplexationID] = len(self.ID_Complexations) # = i
            self.Dir_RevComplexations[i] = ComplexationDirection

            CoeffArray = np.zeros(self.NMax_ComplexationMolParts)
            for MoleculeInfo in ComplexationStoichiometry.strip().strip('"').strip('[').strip(']').strip('{').split('},'):
                Coeff = None
                Idx = None
                for LabelDataPair in MoleculeInfo.split(','):
                    Label, Data = LabelDataPair.strip().split(': ')
                    if Label == '"coeff"':
                        Coeff = Data
                    elif Label == '"molecule"':
                        MolID = self.CheckIDForModification(Data.strip('"'), self.IDModification)
                        if MolID in self.ID_MetabolitesRef:
                            MolID = self.AddLocalizationToMolName(MolID, '[c]')
                        if MolID not in self.ID2Idx_MolsInCPLXRXN:
                            MolID = self.AddLocalizationToMolName(MolID, '[c]')
                        Idx = self.ID2Idx_MolsInCPLXRXN[MolID]
                CoeffArray[Idx] = Coeff

            self.Coeff_MolsInCPLXRXN.append(CoeffArray)
        self.NUniq_MolsInCPLXRXN = len(self.ID_MolsInCPLXRXN)


class FEquilibrium(FDataset):
    def __init__(self):
        self.MolType_EQMRXNs = 'EQMRXN'

        self.ID_EQMRXNs = list()
        self.ID2Idx_EQMRXNs = dict()
        self.Rate_EQMRXNForward = 0
        self.Rate_EQMRXNReverse = 0
        self.Dir_RevEQMRXNs = 0  # Rev for Reversibility


        self.ID_MolsInEQMRXN = list()
        self.ID2Idx_MolsInEQMRXN = dict()
        self.Coeff_MolsInEQMRXN = list()
        self.NUniq_MolsInEQMRXN = 0

        self.NUniq_EQMRXNs = 0

        self.IDModification = dict()

        # To handle unregistered metabolites
        self.ID_Metabolites = list()

        self.ID_UnregisteredFromEQMRXN = list()   # Mostly additional metabolites
        self.MolType_UnregisteredFromEQMRXN = 'MetabolitesUnregisteredFromEQMRXN'
        self.Count_UnregisteredFromEQMRXN = 0
        self.MW_UnregisteredFromEQMRXN = 0
        self.NUniq_UnregisteredFromEQMRXN = 0

        super().__init__()

    # Reactions and kinetic parameters (K values)
    def setupdata(self, Dataset, Comp):
        self.SetUpMetaboliteIDRef(Dataset)
        # Metabolites = Dataset['Counts4AllMetabolites_DL_2_totalIs_0_0_Metabolism_BulkMolecules.tsv']
        # for Metabolite in Metabolites:
        #     self.ID_Metabolites.append(Metabolite[0][:-3])
        #
        # # This requires complexation reaction flatdata to be added to the master dataset before equilibrium flatdata
        # self.ID_Metabolites.append(MasterDataset.ID_UnregisteredFromCPLXRXN)
        #
        # EQMRXNCounts = Dataset['Counts4AllMols4Equilibrium_DL_47_totalIs_45_10_Equilibrium_BulkMolecules.tsv']
        # for EQMRXNCount in EQMRXNCounts:
        #     MolName, MolTotalCounts, MolRequestedCounts = EQMRXNCount
        #     self.ID_
        #
        #     Process, Stoichiometry, ID, Direction, ForwardRate, ReverseRate, OriginalReverseRate = EQMRXNs


        EQMRXNs = Dataset['equilibriumReactions.tsv']

        BlackListed_EQMRXNID = [
            "CPLX0-3964",  # Full ribosome--this gets formed in the simulation
            "RNAP32-CPLX",  # Full RNA Polymerase--simulation only uses apo form
            "RNAP54-CPLX",  # Full RNA Polymerase--simulation only uses apo form
            "RNAP70-CPLX",  # Full RNA Polymerase--simulation only uses apo form
            "RNAPE-CPLX",  # Full RNA Polymerase--simulation only uses apo form
            "RNAPS-CPLX",  # Full RNA Polymerase--simulation only uses apo form
            "CPLX0-221",  # Full RNA Polymerase--simulation only uses apo form
            "CPLX0-222",  # Full RNA Polymerase--simulation only uses apo form
        ]

        EQMRXNs_Filtered = list()
        EQMRXNs_Removed = list()
        N_Filtered = 0
        for Value in EQMRXNs:
            Filter = True
            for EQMRXNID in BlackListed_EQMRXNID:
                if Value[2].find(EQMRXNID) >= 0:
                    Filter = False
                    N_Filtered += 1
                    EQMRXNs_Removed.append(Value)
                    break
            if Filter:
                EQMRXNs_Filtered.append(Value)
        assert len(EQMRXNs) == (len(EQMRXNs_Filtered) + len(EQMRXNs_Removed))

        EQMRXNs = EQMRXNs_Filtered

        self.NUniq_EQMRXNs = len(EQMRXNs)
        self.Dir_RevEQMRXNs = np.zeros(self.NUniq_EQMRXNs)
        self.Rate_EQMRXNForward = np.zeros(self.NUniq_EQMRXNs)
        self.Rate_EQMRXNReverse = np.zeros(self.NUniq_EQMRXNs)

        # Build indices for stoichiometry matrix
        N_Mol_Total = 0
        N_Mol_Repeated = 0
        Data_NotInMasterDataset_NotAddedToUnregistered_Temp = list()
        for Value in EQMRXNs:
            Process, Stoichiometry, EQMRXNID, EQMRXNDirection, ForwardRate, ReverseRate, OriginalReverseRate = Value
            for MoleculeInfo in Stoichiometry.strip().strip('"').strip('[').strip(']').strip('{').split('},'):
                for LabelDataPair in MoleculeInfo.split(','):
                    Label, Data = LabelDataPair.strip().split(': ')
                    if Label == '"molecule"':
                        N_Mol_Total += 1
                        Data = self.CheckIDForModification(Data.strip('"'), self.IDModification)

                        # Check for unregistered molecules
                        if Data not in Comp.Master.ID_Master:
                            Data = self.AddLocalizationToMolName(Data, '[c]')
                            if Data not in Comp.Master.ID_Master:
                                if Data not in self.ID_UnregisteredFromEQMXRXN:
                                    self.ID_UnregisteredFromEQMRXN.append(Data)
                                else:
                                    # Accumulate the original data without localization tag
                                    Data_NotInMasterDataset_NotAddedToUnregistered_Temp.append(Data)

                        # Data at this point has localization tag
                        if Data not in self.ID_MolsInEQMRXN:
                            self.ID2Idx_MolsInEQMRXN[Data] = len(self.ID_MolsInEQMRXN)
                            self.ID_MolsInEQMRXN.append(Data)
                        else:
                            N_Mol_Repeated += 1


        self.NUniq_UnregisteredFromEQMRXN = len(self.ID_UnregisteredFromEQMRXN)
        self.Count_UnregisteredFromEQMRXN = np.zeros(self.NUniq_UnregisteredFromEQMRXN)
        self.MW_UnregisteredFromEQMRXN = np.zeros(self.NUniq_UnregisteredFromEQMRXN)   # TODO: get values from database

        self.AddToMaster(Comp.Master, self.ID_UnregisteredFromEQMRXN, self.MolType_UnregisteredFromEQMRXN, self.Count_UnregisteredFromEQMRXN, self.MW_UnregisteredFromEQMRXN)

        self.NMax_EQMRXNMolParts = len(self.ID_MolsInEQMRXN)

        for i, Value in enumerate(EQMRXNs):
            Process, Stoichiometry, EQMRXNID, EQMRXNDirection, ForwardRate, ReverseRate, OriginalReverseRate = Value
            assert EQMRXNID not in self.ID_EQMRXNs, '%s' % EQMRXNID
            self.ID_EQMRXNs.append(EQMRXNID)
            self.ID2Idx_EQMRXNs[EQMRXNID] = len(self.ID_EQMRXNs) # = i
            self.Rate_EQMRXNForward[i] = ForwardRate
            self.Rate_EQMRXNReverse[i] = ReverseRate
            self.Dir_RevEQMRXNs[i] = EQMRXNDirection

            CoeffArray = np.zeros(self.NMax_EQMRXNMolParts)
            for MoleculeInfo in Stoichiometry.strip().strip('"').strip('[').strip(']').strip('{').split('},'):
                Coeff = None
                Idx = None
                for LabelDataPair in MoleculeInfo.split(','):
                    Label, Data = LabelDataPair.strip().split(': ')
                    if Label == '"coeff"':
                        Coeff = Data
                    elif Label == '"molecule"':
                        MolID = self.CheckIDForModification(Data.strip('"'), self.IDModification)
                        Idx = self.ID2Idx_MolsInEQMRXN[MolID]

                CoeffArray[Idx] = Coeff

            self.Coeff_MolsInEQMRXN.append(CoeffArray)
        self.NUniq_MolsInEQMRXN = len(self.ID_MolsInEQMRXN)


class FMetabolism(FDataset):
    def __init__(self):
        self.MolType_METRXNs = 'METRXN'

        self.ID_METRXNs = list()
        self.ID2Idx_METRXNs = dict()
        self.Dir_RevMETRXNs = list()   # Rev for Reversibility
        self.ID_Enzymes4METRXN = list()
        self.ID_Substrates4METRXN = list()
        self.ID2Idx_Enzymes4METRXN = dict()

        self.ID_MolsInMETRXN = list()
        self.ID2Idx_MolsInMETRXN = dict()
        self.Coeff_MolsInMETRXN = list()
        self.NUniq_MolsInMETRXN = 0

        self.NUniq_METRXNs = 0

        self.IDModification = dict()

        # Unregistered Metabolites
        self.ID_UnregisteredFromMETRXN_Metabolites = list()   # Mostly additional metabolites
        self.MolType_UnregisteredFromMETRXN_Metabolites = 'UnregisteredMetabolitesFromMETRXN'
        self.Count_UnregisteredFromMETRXN_Metabolites = 0
        self.MW_UnregisteredFromMETRXN_Metabolites = 0
        self.NUniq_UnregisteredFromMETRXN_Metabolites = 0

        # Unregistered Enzymes
        self.ID_UnregisteredFromMETRXN_Enzymes = list()   # Mostly additional Enzymes
        self.MolType_UnregisteredFromMETRXN_Enzymes = 'UnregisteredEnzymesFromMETRXN'
        self.Count_UnregisteredFromMETRXN_Enzymes = 0
        self.MW_UnregisteredFromMETRXN_Enzymes = 0
        self.NUniq_UnregisteredFromMETRXN_Enzymes = 0

        # TCA Cycle
        self.ID_RXNsInTCACycle = list()
        self.ID_Substrates4TCACycle = list()
        self.ID_Enzymes4TCACycle = list()
        self.NUniq_RXNsInTCACycle = 0
        self.Dir_RevMETRXNs = list()  # Rev for Reversibility

        super().__init__()

    # Reactions and kinetic parameters (K values)
    def SetUpData(self, Dataset, Comp):
        self.SetUpMetaboliteIDRef(Dataset)

        METRXNs = Dataset['reactions.tsv']

        Enzymes_Blacklisted = [
            'METHYLGLYREDUCT-MONOMER',
            'ACETYL-COA-CARBOXYLMULTI-CPLX',
            'CPLX0-2901',
            'DIHYDRONEOPTERIN-MONO-P-DEPHOS-MONOMER',
            'ENTF-PANT',
            'ENTMULTI-CPLX',
            'CPLX0-8212',
            'GCVMULTI-CPLX',
            'CPLX0-8205',
            'GUANYLCYC-MONOMER',
            'HYDGLUTSYN-MONOMER',
            'ENTB-CPLX',
            'MANNKIN-MONOMER',
            'PHOSPHO-CHEB',
            'METHYLMALONYL-COA-EPIM-MONOMER',
            'NQOR-CPLX',
            'NMNNUCLEOSID-MONOMER',
            'PHOSPHASERDECARB-CPLX',
            'PROPIONYL-COA-CARBOXY-MONOMER',
            'PYROXALTRANSAM-MONOMER',
            'PYRUFLAVREDUCT-MONOMER',
            'CPLX0-7885',
            'NADNUCLEOSID-MONOMER',
            'ENTF-PANT',
            'NADPPHOSPHAT-MONOMER',
            'GCVMULTI-CPLX',
            'MONOMER0-702',
            'MONOMER0-2838',
            'PYRDAMPTRANS-MONOMER',
            'SAMDECARB-CPLX,',
            'TRANSENOYLCOARED-CPLX',
        ]

        self.NUniq_METRXNs = len(METRXNs)

        ReversibilityBinaryDict = dict()
        ReversibilityBinaryDict['true'] = 1
        ReversibilityBinaryDict['false'] = 0

        # Build indices for stoichiometry matrix
        N_Mol_Total = 0
        N_Mol_Repeated = 0
        for Value in METRXNs:
            RXNID, Stoichiometry, RXNReversibility, RXNEnzymeID = Value
            for MoleculeInfo in Stoichiometry.strip().strip('"').strip('[').strip(']').strip('{').split('},'):
                for MolIDCoeffPair in MoleculeInfo.strip().split(','):
                    MolID, Coeff = MolIDCoeffPair.strip().split(': ')
                    MolID = MolID[1:-1]
                    N_Mol_Total += 1
                    if MolID not in self.ID_MolsInMETRXN:
                        self.ID2Idx_MolsInMETRXN[MolID] = len(self.ID_MolsInMETRXN)
                        self.ID_MolsInMETRXN.append(MolID)
                    else:
                        N_Mol_Repeated += 1

                    if MolID not in Comp.Master.ID_Master:
                        if MolID not in self.ID_UnregisteredFromMETRXN_Metabolites:
                            self.ID_UnregisteredFromMETRXN_Metabolites.append(MolID)

            RXNEnzymeIDs = ast.literal_eval(RXNEnzymeID)
            for EnzymeID in RXNEnzymeIDs:
                if EnzymeID not in Comp.Master.ID_Master:
                    if EnzymeID not in self.ID_UnregisteredFromMETRXN_Enzymes:
                        self.ID_UnregisteredFromMETRXN_Enzymes.append(EnzymeID)

        self.NUniq_UnregisteredFromMETRXN_Metabolites = len(self.ID_UnregisteredFromMETRXN_Metabolites)
        self.Count_UnregisteredFromMETRXN_Metabolites = np.zeros(self.NUniq_UnregisteredFromMETRXN_Metabolites)
        self.MW_UnregisteredFromMETRXN_Metabolites = np.zeros(self.NUniq_UnregisteredFromMETRXN_Metabolites)   # TODO: get values from database
        
        self.AddToMaster(Comp.Master, self.ID_UnregisteredFromMETRXN_Metabolites, self.MolType_UnregisteredFromMETRXN_Metabolites, self.Count_UnregisteredFromMETRXN_Metabolites, self.MW_UnregisteredFromMETRXN_Metabolites)

        self.NUniq_UnregisteredFromMETRXN_Enzymes = len(self.ID_UnregisteredFromMETRXN_Enzymes)
        self.Count_UnregisteredFromMETRXN_Enzymes = np.zeros(self.NUniq_UnregisteredFromMETRXN_Enzymes)
        self.MW_UnregisteredFromMETRXN_Enzymes = np.zeros(self.NUniq_UnregisteredFromMETRXN_Enzymes)   # TODO: get values from database
        
        self.AddToMaster(Comp.Master, self.ID_UnregisteredFromMETRXN_Enzymes, self.MolType_UnregisteredFromMETRXN_Enzymes, self.Count_UnregisteredFromMETRXN_Enzymes, self.MW_UnregisteredFromMETRXN_Enzymes)

        self.NMax_METRXNMolParts = len(self.ID_MolsInMETRXN)

        SubstratesToSkip = [
            "WATER",
            "CARBON-DIOXIDE",
            "PROTON",
            "FAD",
            "NAD",
            "NADP",
            "NADH",
            "NADPH",
            "FADH2",
            "CO-A",
            "ATP",
        ]

        LocalizationTags = ['[c]', '[p]']

        SubstratesToSkipWithTags = list()
        for Substrate in SubstratesToSkip:
            for Tag in LocalizationTags:
                SubstrateWithTag = Substrate + Tag
                SubstratesToSkipWithTags.append(SubstrateWithTag)

        for i, Value in enumerate(METRXNs):
            RXNID, Stoichiometry, RXNReversibility, RXNEnzymeID = Value
            assert RXNID not in self.ID2Idx_METRXNs
            self.ID2Idx_METRXNs[RXNID] = len(self.ID_METRXNs) # = i
            self.ID_METRXNs.append(RXNID)
            RXNReversibilityBinary = ReversibilityBinaryDict[RXNReversibility]
            self.Dir_RevMETRXNs.append(RXNReversibilityBinary)

            if RXNEnzymeID.strip('[').strip(']'):
                RXNEnzymeID = ast.literal_eval(RXNEnzymeID.strip('[').strip(']'))
                # TODO: Only one enzyme is taken for now
                EnzymeID = None
                if isinstance(RXNEnzymeID, tuple):
                    EnzymeID = RXNEnzymeID[0]
                    # for ID in RXNEnzymeID:
                    #     self.ID_Enzymes4METRXN.append(ID)
                else:
                    EnzymeID = RXNEnzymeID
                if EnzymeID not in Comp.Master.ID_Master:
                    print('Metabolism Dataset | ' + 'EnzymeID not found in ID_Master: %s' % EnzymeID)
                assert EnzymeID in Comp.Master.ID_Master, 'Metabolism Dataset | ' + 'EnzymeID not found in ID_Master: %s' % EnzymeID
                self.ID2Idx_Enzymes4METRXN[EnzymeID] = len(self.ID_METRXNs)
                self.ID_Enzymes4METRXN.append(EnzymeID)

            CoeffArray = np.zeros(self.NMax_METRXNMolParts)
            Substrate = None
            for MoleculeInfo in Stoichiometry.strip().strip('"').strip('[').strip(']').strip('{').split('},'):
                for MolIDCoeffPair in MoleculeInfo.split(','):
                    MolID, Coeff = MolIDCoeffPair.strip().strip('}').split(': ')
                    MolID = MolID[1:-1]
                    N_Mol_Total += 1
                    Idx = self.ID2Idx_MolsInMETRXN[MolID]
                    CoeffArray[Idx] = Coeff

                    # Substrate
                    if int(Coeff) < 0:
                        if Substrate == None:
                            Substrate = MolID
                        elif Substrate in SubstratesToSkipWithTags:
                            Substrate = MolID
                        elif Substrate.__contains__('+'):
                            Substrate = MolID

            self.ID_Substrates4METRXN.append(Substrate)   # Only one substrate taken from one reaction
            self.Coeff_MolsInMETRXN.append(CoeffArray)
        self.NUniq_MolsInMETRXN = len(self.ID_MolsInMETRXN)


        TCACycle = Dataset['DL_TCAcycleI_prokaryotic.tsv']

        for i, RXN in enumerate(TCACycle):
            ID_RXN, ID_Substrates, ID_Enzymes = RXN
            assert ID_RXN in self.ID_METRXNs, 'Metabolism Dataset | ' + 'Reaction ID not found in ID_METRXNs: %s' % ID_RXN
            self.ID_RXNsInTCACycle.append(ID_RXN)
            for ID_Substrate in ID_Substrates.split('//'):
                assert ID_Substrate in Comp.Master.ID_Master, 'Metabolism Dataset | ' + 'Substrate ID not found in ID_Master: %s' % ID_Substrate
            for ID_Enzyme in ID_Enzymes.split('//'):
                assert ID_Enzyme in Comp.Master.ID_Master, 'Metabolism Dataset | ' + 'Enzyme ID not found in ID_Master: %s' % ID_Enzyme

        self.NUniq_RXNsInTCACycle = len(self.ID_RXNsInTCACycle)


class FKinetics(FDataset):
    def __init__(self): # MasterLocalizations
        self.MolType_Kinetics = 'Kinetics'

        self.ID_KINRXNs = list()
        self.ID2Idx_RXNID2KIN = dict()
        self.ID_Enzymes4KINRXN = list()
        self.ID_Substrates4KINRXN = list()
        self.Const_Temp = 0
        self.Const_Kcat = 0
        self.Const_Km = 0
        self.Const_Ki = 0
        self.Const_Kcat_Default = 0
        self.Const_Km_Default = 0

        self.ID_UnregisteredFromKINRXN = list()   # Mostly additional metabolites
        self.MolType_UnregisteredFromKINRXN = 'MetabolitesUnregisteredFromKINRXN'
        self.Count_UnregisteredFromKINRXN = 0
        self.MW_UnregisteredFromKINRXN = 0
        self.NUniq_UnregisteredFromKINRXN = 0

        self.NUniq_KINRXNs = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, Comp):
        # This kinetic table has no info on pH condition
        Kinetics = Dataset['enzymeKinetics_Sorted_DL.tsv']

        METRXNs = Dataset['reactions.tsv']
        ID_METRXNs = list()
        for METRXN in METRXNs:
            ID_METRXNs.append(METRXN[0])   # METRXNID

        Kinetics_Filtered = list()

        N_RXNsNotInMETRXNID = 0
        N_RXNsWithMoreThanOneSubstrate = 0
        N_SubstratesNotInMasterID = 0
        N_RXNsRedundantKINRXN = 0

        # Flatdata quality check for unregistered data 
        for Value in Kinetics:
            PubmedID, Temp, Checked, Notes, Exclude, ConcentrationSubstrates, ReactionClassID, ReactionID, EnzymeIDs,\
            SubstrateIDs, RateEquationType, Direction, kcat, kM, kI, CustomRateEquation, CustomParameters,\
            CustomParameterConstants, CustomParameterConstantValues, CustomParameterVariables, CheckedBy = Value

            SubstrateIDs = ast.literal_eval(SubstrateIDs)

            # Check unregistered RXN ID
            if ReactionID not in ID_METRXNs:
                N_RXNsNotInMETRXNID += 1
                print('Kinetics Dataset | ', 'ReactionID not in ID_METRXNs: ', ReactionID)

            # Check and add unregistered molecule IDs
            if EnzymeIDs not in Comp.Master.ID_Master:   # All EnzymeIDs contain only one EnzymeID
                self.ID_UnregisteredFromKINRXN.append(EnzymeIDs)
                print('Kinetics Dataset | ', 'EnzymeID not in ID_Master: ', EnzymeIDs)

            for i, SubstrateID in enumerate(SubstrateIDs):
                SubstrateID += '[c]'
                if SubstrateID not in Comp.Master.ID_Master:
                    self.ID_UnregisteredFromKINRXN.append(SubstrateID)
                    N_SubstratesNotInMasterID += 1
                    if self.Switch4DebugDataset:
                        print('Kinetics Dataset | ', 'SubstrateID not in ID_Master: ', SubstrateID)
                if i == 1:
                    N_RXNsWithMoreThanOneSubstrate += 1

            # Load data with the filters below
            # Exclude for repeated Reaction ID
            if ReactionID in self.ID_KINRXNs:
                N_RXNsRedundantKINRXN += 1
                continue

            # Exclude the data requiring custom rate equations
            elif RateEquationType != 'standard':
                continue

            else:
                self.ID_KINRXNs.append(ReactionID)
                Kinetics_Filtered.append(Value)

        if self.Switch4DebugDataset:
            print('Kinetics Dataset | ', '# of RXNs not found in ID_METRXNIDs', N_RXNsNotInMETRXNID)
            print('Kinetics Dataset | ', '# of RXNs with multiple substrate IDs: ', N_RXNsWithMoreThanOneSubstrate)
            print('Kinetics Dataset | ', '# of Substrate IDs not found in ID_Master', N_SubstratesNotInMasterID)
            print('Kinetics Dataset | ', '# of RXNs repeated in ID_KMaster', N_SubstratesNotInMasterID)

        self.NUniq_UnregisteredFromKINRXN = len(self.ID_UnregisteredFromKINRXN)
        self.Count_UnregisteredFromKINRXN = np.zeros(self.NUniq_UnregisteredFromKINRXN)
        self.MW_UnregisteredFromKINRXN = np.zeros(self.NUniq_UnregisteredFromKINRXN)  # TODO: get values from database

        self.AddToMaster(Comp.Master, self.ID_UnregisteredFromKINRXN, self.MolType_UnregisteredFromKINRXN,
                         self.Count_UnregisteredFromKINRXN, self.MW_UnregisteredFromKINRXN)

        # Data extracting
        Kinetics = Kinetics_Filtered

        self.NUniq_KINRXNs = len(Kinetics)
        self.Const_Temp = np.zeros(self.NUniq_KINRXNs)
        self.Const_Kcat = np.zeros(self.NUniq_KINRXNs)
        self.Const_Km = np.zeros(self.NUniq_KINRXNs)
        self.Const_Ki = np.zeros(self.NUniq_KINRXNs)

        # TODO: Reactions with multiple entries and multiple substrate IDs need to be fully compiled. Currently only one adjusted values have been filtered and loaded.
        for i, Value in enumerate(Kinetics):
            PubmedID, Temp, Checked, Notes, Exclude, ConcentrationSubstrates, ReactionClassID, ReactionID, EnzymeIDs, \
            SubstrateIDs, RateEquationType, Direction, kcat, kM, kI, CustomRateEquation, CustomParameters, \
            CustomParameterConstants, CustomParameterConstantValues, CustomParameterVariables, CheckedBy = Value

            # TODO: Currently only accepting a single substrate from the flatdata
            SubstrateIDs = ast.literal_eval(SubstrateIDs)
            SubstrateID_Single = str(SubstrateIDs[0]) + '[c]'

            if Temp == '':
                Temp = 25
            else:
                Temp = int(Temp)

            # Undetermined kM values are arbitrarily set for now
            self.Const_Kcat_Default = 10 ** -4
            self.Const_Km_Default = 10 ** -1

            # TODO: Currently only accepting a single kM from the flatdata,
            kcat_Adjusted = self.ReduceListWithMoreThanTwoElementsToOne(self.ReplaceEmptyListWithValue(ast.literal_eval(kcat), 0))
            kM_Adjusted = self.ReduceListWithMoreThanTwoElementsToOne(self.ReplaceEmptyListWithValue(ast.literal_eval(kM), self.Const_Km_Default))
            kI_Adjusted = self.ReduceListWithMoreThanTwoElementsToOne(self.ReplaceEmptyListWithValue(ast.literal_eval(kI), 0))

            # Adjust kcat based on temperature
            if Temp != 37:
                kcat_Adjusted = (2 ** ((37 - Temp) / 10)) * kcat_Adjusted

            if self.Switch4DebugDataset:
                print('Kinetics Dataset | ', '', 'Temp:', Temp)
                print('Kinetics Dataset | ', '', 'kcat:', kcat, '|', 'kcat_Adjusted', kcat_Adjusted)
                print('Kinetics Dataset | ', '', 'kM:', kM, '|', 'kM_Adjusted', kM_Adjusted)
                print('Kinetics Dataset | ', '', 'kI:', kI, '|', 'kI_Adjusted', kI_Adjusted)

            assert ReactionID in self.ID_KINRXNs
            self.ID2Idx_RXNID2KIN[ReactionID] = len(self.ID2Idx_RXNID2KIN)
            self.ID_Enzymes4KINRXN.append(EnzymeIDs)
            self.ID_Substrates4KINRXN.append(SubstrateID_Single)
            self.Const_Temp[i] = 37
            self.Const_Kcat[i] = kcat_Adjusted
            self.Const_Km[i] = kM_Adjusted
            self.Const_Ki[i] = kI_Adjusted


    def ReplaceEmptyListWithValue(self, Variable, Value):
        if Variable == list():
            Variable = Value
        elif Variable == str():
            Variable = Value
        return Variable

    def ReduceListWithMoreThanTwoElementsToOne(self, Variable):
        if isinstance(Variable, list):
            Variable = Variable[0]
        return Variable

class FCompartment(FDataset):
    def __init__(self): # MasterLocalizations

        self.ID_Compartments = list()
        self.ID2Key_Compartments = dict()
        self.Key_Compartments = list()
        self.Key2Idx_Compartments = dict()
        self.NUniq_Compartments = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, Comp):
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


class FBuildingBlock(FDataset):
    def __init__(self):
        # Lists of molecular specie
        self.Name_dNTPs = list()
        self.Name_NTPs = list()
        self.Name_AAs = list()
        self.Key_dNTPs = list()
        self.Key_NTPs = list()
        self.Key_AAs = list()
        self.NUniq_dNTPs = 0
        self.NUniq_NTPs = 0
        self.NUniq_AAs = 0

        self.Name2Key_BuildingBlocks = dict()

        super().__init__()

    def SetUpData(self, Dataset, Comp):
        dNTPs = Dataset['dntps.txt']
        for BuildingBlock in dNTPs:
            BuildingBlock = self.AddLocalizationToMolName(BuildingBlock, '[c]')
            self.Name_dNTPs.append(BuildingBlock)
        self.Key_dNTPs = ['A', 'C', 'G', 'T']
        self.NUniq_dNTPs = len(self.Name_dNTPs)
        self.Name2Key_BuildingBlocks['dNTPs'] = self.Key_dNTPs

        NTPs = Dataset['ntps.txt']
        for BuildingBlock in NTPs:
            BuildingBlock = self.AddLocalizationToMolName(BuildingBlock, '[c]')
            self.Name_NTPs.append(BuildingBlock)
        self.Key_NTPs = ['A', 'C', 'G', 'U']
        self.NUniq_NTPs = len(self.Name_NTPs)
        self.Name2Key_BuildingBlocks['NTPs'] = self.Key_NTPs

        AAs = Dataset['amino_acids.txt']
        for BuildingBlock in AAs:
            BuildingBlock = self.AddLocalizationToMolName(BuildingBlock, '[c]')
            self.Name_AAs.append(BuildingBlock)
        self.NUniq_AAs = len(self.Name_AAs)

        AAKeys = Dataset['amino_acid_keys.txt']
        for BuildingBlock in AAKeys:
            self.Key_AAs.append(BuildingBlock)
        self.Name2Key_BuildingBlocks['AAs'] = self.Key_AAs


class FUserInput(FDataset):
    def __init__(self):
        self.CellProcesses = list()

        super().__init__()  # MasterLocalizations

    def SetUpData(self, Dataset, Comp):
        CellProcesses = Dataset['process_module.tsv']
        InputReactions = Dataset['reaction_data.tsv']
        # MoleculeData = Dataset['molecule_data.tsv']

        for CellProcess in CellProcesses:
            self.CellProcesses.append(CellProcess[0])

        ReversibilityBinaryDict = dict()
        ReversibilityBinaryDict['true'] = 1
        ReversibilityBinaryDict['false'] = 0

        SubstratesToSkip = [
            "WATER",
            "CARBON-DIOXIDE",
            "PROTON",
            "FAD",
            "NAD",
            "NADP",
            "NADH",
            "NADPH",
            "FADH2",
            "CO-A",
            "ATP",
        ]

        LocalizationTags = ['[c]', '[p]']

        SubstratesToSkipWithTags = list()
        for Substrate in SubstratesToSkip:
            for Tag in LocalizationTags:
                SubstrateWithTag = Substrate + Tag
                SubstratesToSkipWithTags.append(SubstrateWithTag)

        # TODO: Allow unregistered molecules in the future
        # Overwrite reaction information
        if InputReactions:
            for InputReaction in InputReactions:
                RXNID, Stoichiometry, RXNReversibility, RXNEnzymeID = InputReaction

                EnzymeID = None

                if RXNEnzymeID.strip('[').strip(']'):
                    RXNEnzymeID = ast.literal_eval(RXNEnzymeID.strip('[').strip(']'))
                    # TODO: Only one enzyme is taken for now
                    if isinstance(RXNEnzymeID, tuple):
                        EnzymeID = RXNEnzymeID[0]
                        # for ID in RXNEnzymeID:
                        #     self.ID_Enzymes4METRXN.append(ID)
                    else:
                        EnzymeID = RXNEnzymeID
                    if EnzymeID not in Comp.Master.ID_Master:
                        print('Metabolism Dataset | ' + 'EnzymeID not found in ID_Master: %s' % EnzymeID)
                    assert EnzymeID in Comp.Master.ID_Master, 'Metabolism Dataset | ' + 'EnzymeID not found in ID_Master: %s' % EnzymeID

                CoeffArray = np.zeros(Comp.Metabolism.NMax_METRXNMolParts)
                Substrate = None

                for MoleculeInfo in Stoichiometry.strip().strip('"').strip('[').strip(']').strip('{').split('},'):
                    for MolIDCoeffPair in MoleculeInfo.split(','):
                        MolID, Coeff = MolIDCoeffPair.strip().strip('}').split(': ')
                        MolID = MolID[1:-1]
                        Idx = Comp.Metabolism.ID2Idx_MolsInMETRXN[MolID]
                        CoeffArray[Idx] = Coeff

                        # Substrate
                        if int(Coeff) < 0:
                            if Substrate == None:
                                Substrate = MolID
                            elif Substrate in SubstratesToSkipWithTags:
                                Substrate = MolID
                            elif Substrate.__contains__('+'):
                                Substrate = MolID

                if RXNID in Comp.Metabolism.ID_METRXNs:
                    Idx_RXN = Comp.Metabolism.ID2Idx_METRXNs[RXNID]
                    Comp.Metabolism.Dir_RevMETRXNs[Idx_RXN] = ReversibilityBinaryDict[RXNReversibility]

                    Comp.Metabolism.ID_Enzymes4METRXN[Idx_RXN] = EnzymeID
                    Comp.Metabolism.ID_Substrates4METRXN[Idx_RXN] = Substrate
                    Comp.Metabolism.Coeff_MolsInMETRXN[Idx_RXN] = CoeffArray

                else:
                    Comp.Metabolism.ID2Idx_METRXNs[RXNID] = len(Comp.Metabolism.ID_METRXNs)
                    Comp.Metabolism.ID_METRXNs.append(RXNID)
                    Comp.Metabolism.Dir_RevMETRXNs.append(ReversibilityBinaryDict[RXNReversibility])

                    Comp.Metabolism.ID_Enzymes4METRXN.append(EnzymeID)
                    Comp.Metabolism.ID_Substrates4METRXN.append(Substrate)
                    Comp.Metabolism.Coeff_MolsInMETRXN.append(CoeffArray)


# Master Dataset is an exception to the SetUpData method
class FMaster():
    def __init__(self):
        self.Switch4DebugMasterDataset = False

        self.ID2Idx_Master = dict()
        self.Type_Master = list()
        self.ID2Type_Master = dict()
        self.ID_Master = list()
        self.Count_Master = list()
        self.MW_Master = list()
        self.NUniq_Master = 0
        self.ID2ID_Gene2RNA_Master = dict()
        self.ID2ID_mRNA2Protein_Master = dict()

        self.Idx2Idx_Gene2RNA_Master = dict()
        self.Idx2Idx_mRNA2ProteinMonomer_Master = dict()
        
        # MasterIndices
        self.Idx_Master_Chromosomes = list()
        self.Idx_Master_Genes = list()
        self.Idx_Master_Promoters = list()
        self.Idx_Master_RNAs = list()
        self.Idx_Master_mRNAs = list()
        self.Idx_Master_rRNAs = list()
        self.Idx_Master_tRNAs = list()
        self.Idx_Master_miscRNAs = list()
        self.Idx_Master_Proteins = list()
        self.Idx_Master_Complexes = list()
        self.Idx_Master_Metabolites = list()
        self.Idx_Master_UnregisteredFromCPLXRXN = list()
        self.Idx_Master_UnregisteredFromEQMRXN = list()
        self.Idx_Master_UnregisteredFromMETRXN_Metabolites = list()
        self.Idx_Master_UnregisteredFromMETRXN_Enzymes = list()

        self.Idx_Master_MolsInCPLXRXN = list()
        self.Idx_Master_MolsInEQMRXN = list()
        self.Idx_Master_MolsInMETRXN = list()

    def SetUpData(self, Comp):
        self.SetUpMasterIndices(Comp)
        self.SetUpIdx2IdxMappingDict()

    def SetUpMasterIndices(self, Comp):
        self.Idx_Master_Chromosomes = self.GetMolIdx(Comp.Chromosome.ID_Chromosomes, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_Chromosomes) == Comp.Chromosome.NMax_Chromosomes

        self.Idx_Master_Genes = self.GetMolIdx(Comp.Gene.ID_Genes, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_Genes) == Comp.Gene.NUniq_Genes

        # self.Idx_Master_Promoters = self.GetMolIdx(Comp.Promoter.ID_Promoters, Comp.Master.ID2Idx_Master)

        self.Idx_Master_RNAs = self.GetMolIdx(Comp.RNA.ID_RNAs, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_RNAs) == Comp.RNA.NUniq_RNAs

        self.Idx_Master_Proteins = self.GetMolIdx(Comp.Protein.ID_Proteins, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_Proteins) == Comp.Protein.NUniq_Proteins

        self.Idx_Master_Complexes = self.GetMolIdx(Comp.Complex.ID_Complexes, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_Complexes) == Comp.Complex.NUniq_Complexes

        self.Idx_Master_Metabolites = self.GetMolIdx(Comp.Metabolite.ID_Metabolites, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_Metabolites) == Comp.Metabolite.NUniq_Metabolites

        self.Idx_Master_UnregisteredFromCPLXRXN = self.GetMolIdx(Comp.Complexation.ID_UnregisteredFromCPLXRXN, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_UnregisteredFromCPLXRXN) == Comp.Complexation.NUniq_UnregisteredFromCPLXRXN

        self.Idx_Master_UnregisteredFromEQMRXN = self.GetMolIdx(Comp.Equilibrium.ID_UnregisteredFromEQMRXN, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_UnregisteredFromEQMRXN) == Comp.Equilibrium.NUniq_UnregisteredFromEQMRXN

        self.Idx_Master_UnregisteredFromMETRXN_Metabolites = self.GetMolIdx(Comp.Metabolism.ID_UnregisteredFromMETRXN_Metabolites, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_UnregisteredFromMETRXN_Metabolites) == Comp.Metabolism.NUniq_UnregisteredFromMETRXN_Metabolites

        self.Idx_Master_UnregisteredFromMETRXN_Enzymes = self.GetMolIdx(Comp.Metabolism.ID_UnregisteredFromMETRXN_Enzymes, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_UnregisteredFromMETRXN_Enzymes) == Comp.Metabolism.NUniq_UnregisteredFromMETRXN_Enzymes

        self.Idx_Master_MolsInCPLXRXN = self.GetMolIdx(Comp.Complexation.ID_MolsInCPLXRXN, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_MolsInCPLXRXN) == Comp.Complexation.NUniq_MolsInCPLXRXN

        self.Idx_Master_MolsInEQMRXN = self.GetMolIdx(Comp.Equilibrium.ID_MolsInEQMRXN, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_MolsInEQMRXN) == Comp.Equilibrium.NUniq_MolsInEQMRXN

        self.Idx_Master_MolsInMETRXN = self.GetMolIdx(Comp.Metabolism.ID_MolsInMETRXN, Comp.Master.ID2Idx_Master)
        assert len(self.Idx_Master_MolsInMETRXN) == Comp.Metabolism.NUniq_MolsInMETRXN

        for i, Idx_Master in enumerate(self.Idx_Master_RNAs):
            if i in Comp.RNA.Idx_mRNA:
                self.Idx_Master_mRNAs.append(Idx_Master)
            elif i in Comp.RNA.Idx_rRNA:
                self.Idx_Master_rRNAs.append(Idx_Master)
            elif i in Comp.RNA.Idx_tRNA:
                self.Idx_Master_tRNAs.append(Idx_Master)
            elif i in Comp.RNA.Idx_miscRNA:
                self.Idx_Master_miscRNAs.append(Idx_Master)
            else:
                print("Master Dataset | ', 'Unmapped RNA type found in 'self.Idx_Master_RNA'")

        assert len(self.Idx_Master_mRNAs) == Comp.RNA.NUniq_mRNAs
        assert len(self.Idx_Master_rRNAs) == Comp.RNA.NUniq_rRNAs
        assert len(self.Idx_Master_tRNAs) == Comp.RNA.NUniq_tRNAs
        assert len(self.Idx_Master_miscRNAs) == Comp.RNA.NUniq_miscRNAs

        assert len(self.Idx_Master_Genes) == len(self.Idx_Master_RNAs)
        assert len(self.Idx_Master_mRNAs) == len(self.Idx_Master_Proteins)

    def SetUpIdx2IdxMappingDict(self):
        # Index to Index mapping dictionaries
        for i, ID in enumerate(self.ID_Master):
            if self.ID2Type_Master[ID] == 'Gene':
                GeneID = ID
                RNAID = self.ID2ID_Gene2RNA_Master[GeneID]
                GeneIdx = self.ID2Idx_Master[GeneID]
                RNAIdx = self.ID2Idx_Master[RNAID]
                self.Idx2Idx_Gene2RNA_Master['%s' % str(GeneIdx)] = RNAIdx
            elif self.ID2Type_Master[ID] == 'RNA':
                RNAID = ID
                if RNAID not in self.ID2ID_mRNA2Protein_Master.keys():
                    continue
                else:
                    ProteinID = self.ID2ID_mRNA2Protein_Master[RNAID]
                    RNAIdx = self.ID2Idx_Master[RNAID]
                    ProteinIdx = self.ID2Idx_Master[ProteinID]
                    self.Idx2Idx_mRNA2ProteinMonomer_Master['%s' % str(RNAIdx)] = ProteinIdx

    def GetMolIdx(self, Molecules, MolIdxRef):
        MolIdxList = list()
        for Molecule in Molecules:
            MolIdx = MolIdxRef[Molecule]
            MolIdxList.append(MolIdx)
        return MolIdxList

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)