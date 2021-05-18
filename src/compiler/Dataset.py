import os, sys
import csv
import numpy as np
import ast
import re
from abc import abstractmethod

class FDataset():
    def __init__(self):
        self.Switch4DebugDataset = False
        pass

    @abstractmethod
    def SetUpData(self, Dataset, MasterDataset = None):
        pass

    @abstractmethod
    def SaveData(self, SavePath):
        pass

    def RemoveLocalizationFromMolName(self, VariableName):
        return VariableName[:-3]

    def AddToMaster(self, MasterDataset, IDs, Counts, MWs, Localizations='None'):
        # self.AddToMaster(self.MetaboliteIDs, self.MetaboliteCounts, self.MetaboliteMWs)

        for ID in IDs:
            MasterDataset.MasterID2Index[ID] = len(IDs)
            MasterDataset.MasterIDs.append(ID)

        MasterDataset.MasterCounts = np.append(MasterDataset.MasterCounts, Counts)
        MasterDataset.MasterMWs = np.append(MasterDataset.MasterMWs, MWs)
        # if Localizations:
        #     self.MasterLocalizations = np.append(self.MasterLocalizations, Localizations)

        if MasterDataset.Switch4DebugMasterDataset:
            for Key, Value in self.__dict__.items():
                if (Value == IDs) and (Key[:6] != 'Master'):
                    print('MasterDataset addition completed: Class "%s"' % Key[:-3])

                if Key[:6] == 'Master':
                    LengthOfValue = len(Value)
                    print('The length of "%s": %d' % (Key, LengthOfValue))

        assert len(MasterDataset.MasterIDs) == len(MasterDataset.MasterID2Index) and len(MasterDataset.MasterIDs) == len(MasterDataset.MasterCounts) and len(MasterDataset.MasterIDs) == len(MasterDataset.MasterMWs)


class FMetabolite(FDataset):
    def __init__(self):
        self.MetaboliteID2MWTemp = {}
        self.MetaboliteIDs = []
        self.MetaboliteID2Index = {}
        self.MetaboliteCounts = 0
        self.MetaboliteMWs = 0
        self.N_Metabolites = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        MetaboliteMWs = Dataset['metabolites.tsv']
        WaterMW = Dataset['water.tsv'][0]
        MetaboliteMWs.append(WaterMW)

        for Value in MetaboliteMWs:
            Name, MW, Localization = Value
            assert Name not in self.MetaboliteID2MWTemp
            self.MetaboliteID2MWTemp[Name] = MW

        MetaboliteCounts = Dataset['Counts4AllMetabolites_DL_2_totalIs_0_0_Metabolism_BulkMolecules.tsv']
        self.N_Metabolites = len(MetaboliteCounts)
        self.MetaboliteCounts = np.zeros(self.N_Metabolites)
        self.MetaboliteMWs = np.zeros(self.N_Metabolites)

        for i, Value in enumerate(MetaboliteCounts):
            Name, Count, Request = Value
            NameLocRemoved = self.RemoveLocalizationFromMolName(Name)
            assert Name not in self.MetaboliteID2Index
            self.MetaboliteID2Index[Name] = len(self.MetaboliteIDs)  # = i
            self.MetaboliteIDs.append(Name)
            self.MetaboliteCounts[i] = Count
            self.MetaboliteMWs[i] = self.MetaboliteID2MWTemp[NameLocRemoved]
            assert self.MetaboliteMWs[i] != 0
        
        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.MetaboliteIDs, self.MetaboliteCounts, self.MetaboliteMWs)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FChromosome(FDataset):
    def __init__(self): # MasterLocalizations
        self.ChromosomeLengthsInGenome = 0
        self.Genome = None

        self.NMax_Chromosomes = 0
        self.ChromosomeIDs = []
        self.ChromosomeBPCounts = 0
        self.DNABasePairMWs = 0
        self.N_ChromosomesInGenome = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        Genome = Dataset['EscherichiaColi.fasta']
        # Genome is a dictionary where key is 'Ch#' and value is its seq.

        self.Genome = Genome
        self.N_ChromosomesInGenome = len(self.Genome)

        self.ChromosomeLengthsInGenome = np.zeros(self.N_ChromosomesInGenome)
        for i, (ChromosomeName, ChromosomeSeq) in enumerate(self.Genome.items()):
            self.ChromosomeLengthsInGenome[i] = len(ChromosomeSeq)

        # Set up an arbitrary max number to keep NT Count arrays for all chromosomes
        MaxCopyNumberOfEachChromosome = 5
        self.NMax_Chromosomes = self.N_ChromosomesInGenome * MaxCopyNumberOfEachChromosome
        self.ChromosomeLengths = np.zeros(self.NMax_Chromosomes)
        self.ChromosomeBPCounts = np.zeros(self.NMax_Chromosomes)
        self.DNABasePairMWs = np.zeros(self.NMax_Chromosomes)
        DNABasePairMW = 650

        for i in range(self.NMax_Chromosomes):
            for j in range(self.N_ChromosomesInGenome):
                ChromosomeNumber = j + 1
                ChromosomeReplication = i + 1
                ChromosomeID = 'Ch%d_%d' % (ChromosomeNumber, ChromosomeReplication)
                self.ChromosomeIDs.append(ChromosomeID)
                if i == 0:
                    self.ChromosomeBPCounts[j] = self.ChromosomeLengthsInGenome
                self.DNABasePairMWs[i] = DNABasePairMW
        # average molecular weight of one base pair: 650 Daltons, which is 650 g/mol

        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.ChromosomeIDs, self.ChromosomeBPCounts, self.DNABasePairMWs)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FGene(FDataset):
    def __init__(self): # MasterLocalizations
        self.GeneLengths = 0
        self.GeneCoordinates = 0
        self.GeneDirections = []
        self.GeneNames = []
        self.GeneSymbols = []
        self.GeneIDs = []
        self.GeneSeqs = []
        self.GeneTypes = []
        self.GeneID2Index = {}
        self.GeneSymbol2Index = {}
        self.GeneName2ID = {}
        self.GeneCounts = 0
        self.N_Genes = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        Genes = Dataset['genes.tsv']
        self.N_Genes = len(Genes)
        self.GeneLengths = np.zeros(self.N_Genes)
        self.GeneCoordinates = np.zeros(self.N_Genes)
        self.GeneNames = []
        self.GeneSymbols = []
        self.GeneIDs = []
        self.GeneTypes = []
        self.GeneName2ID = {}
        self.GeneID2Index = {}
        self.GeneSymbol2Index = {}
        
        for i, Value in enumerate(Genes):
            Length, Name, Seq, RNAID, Coordinate, Direction, Symbol, Type, GeneID, MonomerID = Value
            # assert Name not in self.GeneName2ID, '"%s" has a redundant name in GeneName2ID dictionary' % Name
            # self.GeneNames.append(Name)
            # self.GeneName2ID[Name] = GeneID
            self.GeneSymbols.append(Symbol)
            self.GeneSymbol2Index[Symbol] = len(self.GeneIDs)
            assert GeneID not in self.GeneID2Index
            self.GeneIDs.append(GeneID)
            self.GeneID2Index[GeneID] = len(self.GeneSeqs)
            self.GeneSeqs.append(Seq)
            self.GeneLengths[i] = (len(Seq))
            self.GeneCoordinates[i] = Coordinate
            self.GeneDirections.append(Direction)
            self.GeneTypes.append(Type)

        self.GeneCounts = np.ones(self.N_Genes)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FPromoter(FDataset):
    def __init__(self): # MasterLocalizations
        self.PromoterCoordinates = 0
        self.PromoterDirections = []
        self.PromoterIDs = []
        self.PromoterID2Index = {}
        self.PromoterTargetGenesSymbols = []
        self.N_Promoters = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        Promoters = Dataset['promoters.tsv']

        self.N_Promoters = len(Promoters)
        self.PromoterCoordinates = np.zeros(self.N_Promoters)
        IrregularNames = 0

        for i, Value in enumerate(Promoters):
            Position, Direction, PromoterID, Name = Value
            self.PromoterIDs.append(PromoterID)
            self.PromoterID2Index = len(self.PromoterDirections)  # i
            self.PromoterDirections.append(Direction)
            self.PromoterCoordinates[i] = Position
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
                if TargetGeneSymbol in self.PromoterTargetGenesSymbols:
                    print('The gene symbol %s is expressed by multiple promoters' % TargetGeneSymbol)
            self.PromoterTargetGenesSymbols.append(TargetGeneSymbol)
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
    def __init__(self): # MasterLocalizations
        self.RNAName2CountTemp = {}
        self.RNANames = []
        self.RNAName2ID = {}
        self.RNAID2Index = {}
        self.RNAIDs = []
        self.RNAHalfLives = 0
        self.RNASeqs = []
        self.RNALengths = 0
        self.RNATypes = []
        self.RNAMWs = 0
        self.RNANTCounts = []
        self.RNANTFreqs = []

        self.N_RNAs = 0
        self.N_mRNAs = 0
        self.N_tRNAs = 0
        self.N_rRNAs = 0
        self.N_miscRNAs = 0

        self.RNATypeIndex4AllRNA = []
        self.RNATypeIndex4mRNA = []
        self.RNATypeIndex4tRNA = []
        self.RNATypeIndex4rRNA = []
        self.RNATypeIndex4miscRNA = []

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        RNACounts = Dataset['Counts4AllRNAs_DL_19_totalIs_17_3_TranscriptElongation_BulkMolecules.tsv']
        for Value in RNACounts:
            Name, Count, Request = Value
            assert Name not in self.RNAName2CountTemp
            self.RemoveLocalizationFromMolName(Name)
            Name = self.RemoveLocalizationFromMolName(Name)
            self.RNAName2CountTemp[Name] = int(Count)

        RNAs = Dataset['rnas.tsv']
        self.N_RNAs = len(RNAs)
        self.RNAHalfLives = np.zeros(self.N_RNAs)
        self.RNALengths = np.zeros(self.N_RNAs)
        self.RNAMWs = np.zeros(self.N_RNAs)
        self.RNACounts = np.zeros(self.N_RNAs)

        for i, Value in enumerate(RNAs):
            HalfLife, Name, Seq, Type, ModifiedForms, MonomerID, Comments, MW, Location, NTCount, RNAID, GeneID, MicArrExp = Value
            assert Name not in self.RNAID2Index
            self.RNANames.append(Name)
            self.RNAName2ID[Name] = RNAID
            self.RNAID2Index[RNAID] = len(self.RNAIDs)
            self.RNAIDs.append(RNAID)
            self.RNAHalfLives[i] = HalfLife
            self.RNASeqs.append(Seq)
            self.RNALengths[i] = (len(Seq))
            self.RNATypes.append(Type)
            self.RNAMWs[i] = max(np.fromstring(MW[1:-1], dtype='float32', sep=','))
            assert self.RNAMWs[i] != 0
            self.RNANTCounts.append(NTCount)
            self.RNACounts[i] = self.RNAName2CountTemp[self.RNAName2ID[Name]]
            self.RNATypeIndex4AllRNA.append(i)
            if Type == 'mRNA':
                self.RNATypeIndex4mRNA.append(i)
            elif Type == 'tRNA':
                self.RNATypeIndex4tRNA.append(i)
            elif Type == 'rRNA':
                self.RNATypeIndex4rRNA.append(i)
            elif Type == 'miscRNA':
                self.RNATypeIndex4miscRNA.append(i)
            else:
                print('Warning: Unaccounted RNA type detected: ', Type)
            # RNANTFreq
            NTTotalCount = len(Seq)
            NTCounts = np.array([Seq.count("A"), Seq.count("C"), Seq.count("G"), Seq.count("U")])
            if np.sum(NTCounts) != NTTotalCount:
                print("WARNING: RNA seq may contain non-ACGU text.", file=sys.stderr)
            NTFreq = NTCounts / NTTotalCount
            self.RNANTCounts[i] = NTCounts
            self.RNANTFreqs = NTFreq

        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.RNAIDs, self.RNACounts, self.RNAMWs)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FProtein(FDataset):
    def __init__(self): # MasterLocalizations
        self.ProteinName2CountTemp = {}
        self.ProteinNames = []
        self.ProteinName2ID = {}
        self.ProteinID2Index = {}
        self.ProteinIDs = []
        self.ProteinSeqs = []
        self.ProteinLengths = 0
        self.ProteinMWs = 0
        self.ProteinLocations = []
        self.ProteinAACounts = []
        self.ProteinCounts = 0
        self.ProteinAAFreqs = []
        self.ProteinID2GeneID = {}
        self._GeneIDs = []
        self.ProteinID2RNAID = {}
        self._RNAIDs = []

        self.N_Proteins = 0

        super().__init__() # MasterLocalizations

        return

    # Protein Monomers
    def SetUpData(self, Dataset, MasterDataset = None):
        ProteinCounts = Dataset['Counts4AllProteinMonomers_DL_28_totalIs_26_5_PolypeptideElongation_BulkMolecules.tsv']
        for Value in ProteinCounts:
            Name, Count, Request = Value
            assert Name not in self.ProteinName2CountTemp
            # self.RemoveLocalizationFromMolName(Name)
            Name = self.RemoveLocalizationFromMolName(Name)
            self.ProteinName2CountTemp[Name] = int(Count)

        Proteins = Dataset['proteins.tsv']
        self.N_Proteins = len(Proteins)
        self.ProteinLengths = np.zeros(self.N_Proteins)
        self.ProteinMWs = np.zeros(self.N_Proteins)
        self.ProteinCounts = np.zeros(self.N_Proteins)

        for i, Value in enumerate(Proteins):
            AACount, Name, Seq, Comments, CodingRNASeq, MW, Location, RNAID, ProtMonomerID, GeneID = Value
            self.ProteinNames.append(Name)
            self.ProteinName2ID[Name] = ProtMonomerID
            self.ProteinID2Index[ProtMonomerID] = len(self.ProteinIDs)
            self.ProteinIDs.append(ProtMonomerID)
            self.ProteinSeqs.append(Seq)
            self.ProteinLengths[i] = len(Seq)
            self.ProteinMWs[i] = max(np.fromstring(MW[1:-1], dtype='float32', sep=','))
            assert self.ProteinMWs[i] != 0
            self.ProteinLocations.append(Location)
            self.ProteinAACounts.append(AACount)
            self.ProteinCounts[i] = self.ProteinName2CountTemp[self.ProteinName2ID[Name]]
            self.ProteinID2RNAID[ProtMonomerID] = RNAID
            # ProteinAAFreq
            AATotalCount = len(Seq)
            if AATotalCount != self.ProteinLengths[i]:
                print("WARNING: Protein length data may contain error.", file=sys.stderr)
            AACounts = np.array(list(map(int, AACount[1:-1].split(','))))
            AAFreq = AACounts / AATotalCount
            self.ProteinAAFreqs.append(AAFreq)
        
        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.ProteinIDs, self.ProteinCounts, self.ProteinMWs)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FComplex(FDataset):
    def __init__(self): # MasterLocalizations
        self.ProteinIDs = []

        self.ComplexName2CountTemp = {}
        self.ComplexNameWithoutLocalization2CountTemp = {}
        self.ComplexNames = []
        self.ComplexName2ID = {}
        self.ComplexIDs = []
        self.ComplexID2Index = {}
        self.ComplexCounts = 0
        self.ComplexMWs = 0
        self.ComplexLocations = []

        self.N_Complexes = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        Proteins = Dataset['proteins.tsv']
        for Value in Proteins:
            AACount, Name, Seq, Comments, CodingRNASeq, MW, Location, RNAID, ProtMonomerID, GeneID = Value
            self.ProteinIDs.append(ProtMonomerID)

        ComplexCounts = Dataset['Counts4AllProteinMonomersAndCPLXsForComplexation_DL_44_totalIs_42_8_Complexation_BulkMolecules.tsv']
        for Value in ComplexCounts:
            Name, Count, Request = Value
            NameWithoutLocalization = self.RemoveLocalizationFromMolName(Name)
            assert Name not in self.ComplexName2CountTemp
            assert NameWithoutLocalization not in self.ComplexName2CountTemp
            if (Name or NameWithoutLocalization) in self.ProteinIDs:
                continue
            # self.RemoveLocalizationFromMolName(Name)
            self.ComplexName2CountTemp[Name] = int(Count)
            self.ComplexNameWithoutLocalization2CountTemp[NameWithoutLocalization] = int(Count)

        Complexes = Dataset['proteinComplexes_large.tsv'] + Dataset['proteinComplexes.tsv']
        for Value in Complexes:
            Name, Comments, MW, Location, RxnID, ComplexID = Value
            if ComplexID in self.ComplexIDs:
                continue
            self.ComplexNames.append(Name)
            self.ComplexName2ID[Name] = ComplexID
            self.ComplexID2Index[ComplexID] = len(self.ComplexIDs)
            self.ComplexIDs.append(ComplexID)
        self.N_Complexes = len(self.ComplexIDs)
        self.ComplexMWs = np.zeros(self.N_Complexes)
        self.ComplexCounts = np.zeros(self.N_Complexes)

        for i, Value in enumerate(Complexes):
            Name, Comments, MW, Location, RxnID, ComplexID = Value
            if ComplexID in self.ComplexIDs:
                continue
            self.ComplexMWs[i] = max(np.fromstring(MW[1:-1], dtype='float32', sep=','))
            assert self.ComplexMWs[i] != 0
            self.ComplexLocations.append(Location)
            if ComplexID in self.ComplexNameWithoutLocalization2CountTemp:
                Count = self.ComplexNameWithoutLocalization2CountTemp[ComplexID]
                self.ComplexCounts[i] = Count
            else:
                if ComplexID not in self.ProteinIDs:
                    print('%s is not found in either Protein or Complex IDs' % ComplexID)
                # assert ComplexID in self.ProteinIDs, '%s is not found in either Protein or Complex IDs' % ComplexID

        # Add to the Master Dataset
        self.AddToMaster(MasterDataset, self.ComplexIDs, self.ComplexCounts, self.ComplexMWs)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FReaction(FDataset):
    def __init__(self): # MasterLocalizations

        self.RXNIDs = []
        self.RXNID2Index = {}
        self.RXNStoichiometries = []
        self.RXNReversibilities = []
        self.RXNEnzymeIDs = []
        self.RXNSubstrates = []
        self.RXNSubstrateStoichs = []
        self.RXNProducts = []
        self.RXNProductStoichs = []

        self.N_RXNs = 0

        super().__init__() # MasterLocalizations

    # Reactions and kinetic parameters (K values)
    def SetUpData(self, Dataset, MasterDataset = None):
        RXNs = Dataset['reactions.tsv']
        for i, Value in enumerate(RXNs):
            RXNID, RXNStoichiometry, RXNReversibility, RXNEnzymeID = Value
            assert RXNID not in self.RXNID2Index
            self.RXNID2Index[RXNID] = len(self.RXNIDs) # = i
            self.RXNIDs.append(RXNID)
            self.RXNStoichiometries.append(RXNStoichiometry)
            self.RXNReversibilities.append(RXNReversibility)
            self.RXNEnzymeIDs.append(RXNEnzymeID)
            Substrates = []
            SubstrateStoichs = []
            Products = []
            ProductStoichs = []
            RXNStoichiometry = ast.literal_eval(RXNStoichiometry)
            for Mol, Stoich in RXNStoichiometry.items(): # RXNStoichiometry is a string dictionary
                self.RemoveLocalizationFromMolName(Mol)
                if Stoich < 0:
                    Substrates.append(Mol)
                    SubstrateStoichs.append(-Stoich)
                if Stoich > 0:
                    Products.append(Mol)
                    ProductStoichs.append(Stoich)
            assert SubstrateStoichs or ProductStoichs > 0, '%s: The exported Stoich value has to be greater than 0' % RXNID
            self.RXNSubstrates.append([Substrates])
            self.RXNSubstrateStoichs.append([SubstrateStoichs])
            self.RXNProducts.append([Products])
            self.RXNProductStoichs.append([ProductStoichs])

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FEnzyme(FDataset):
    def __init__(self): # MasterLocalizations

        self.KinRXNID2KineticData = {}
        self.KinRXNClassIDs = []
        self.KinEnzymeIDs = []
        self.KinSubstrateIDs = []
        self.KinTemperature = 0
        self.KinKcats = 0
        self.KinKms = []
        self.KinKis = 0

        self.N_KinRXNs = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        # This kinetic table has no info on pH condition
        Kinetics = Dataset['enzymeKinetics.tsv']
        self.N_KinRXNs = len(Kinetics)
        self.KinTemperature = np.zeros(self.N_KinRXNs)
        self.KinKcats = np.zeros(self.N_KinRXNs)
        self.KinKis = np.zeros(self.N_KinRXNs)
        N_TempAbscent = 0
        N_KcatValueAbscent = 0
        N_KmValueAbscent = 0
        N_KiValueAbscent = 0
        N_KmValueUnmatched = 0

        for i, Value in enumerate(Kinetics):
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
                self.KinTemperature[i] = float(Temp)
            else:
                N_TempAbscent += 1
                if self.Switch4Debug:
                    print('TempAbscent:', ReactionID)

            kcat = kcat.replace(',', '').replace('[', '').replace(']', '')
            if kcat:
                self.KinKcats[i] = float(kcat)
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
                    self.KinKms.append(kM_Rebuilt)
                else:
                    N_KmValueUnmatched += 1
                    if self.Switch4DebugDataset:
                        print(ReactionID, ': The number of kM values do not match the number of substrates')
                    self.KinKms.append(0.0)
            else:
                if kM:
                    self.KinKms.append(float(kM))
                else:
                    self.KinKms.append(0.0)
                    N_KmValueAbscent += 1
                    if self.Switch4DebugDataset:
                        print('KmValueAbscent:', ReactionID)

            kI = kI.replace(',', '').replace('[', '').replace(']', '')
            if kI:
                self.KinKis[i] = float(kI)
            else:
                N_KiValueAbscent += 1
                if self.Switch4DebugDataset:
                    print('KiValueAbscent:', ReactionID)
        Km_ExportShort = self.N_KinRXNs - len(self.KinKms)
        assert len(self.KinKms) == self.N_KinRXNs, '%s Km exporting fell short' % Km_ExportShort
        if self.Switch4DebugDataset:
            print('N_Total:', self.N_KinRXNs)
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

        self.CompartmentIDs = []
        self.CompartmentKeys = []
        self.CompartmentKey2Index = {}
        self.N_Compartments = 0

        super().__init__() # MasterLocalizations

    def SetUpData(self, Dataset, MasterDataset = None):
        Compartments = Dataset['compartments.tsv']
        self.N_Compartments = len(Compartments)
        self.CompartmentIndexes = np.zeros(self.N_Compartments)

        for i, Value in enumerate(Compartments):
            Abbrev, ID = Value
            self.CompartmentKey2Index[Abbrev] = len(self.CompartmentKeys)
            self.CompartmentKeys.append(Abbrev)
            self.CompartmentIDs.append(ID)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FBuildingBlock(FDataset):
    def __init__(self):
        # Lists of molecular specie
        self.dNTPs = []
        self.NTPs = []
        self.AAs = []
        self.dNTPKeys = []
        self.NTPKeys = []
        self.AAKeys = []
        self.N_dNTPs = 0
        self.N_NTPs = 0
        self.N_AAs = 0

        super().__init__()

    def SetUpData(self, Dataset, MasterDataset = None):
        dNTPs = Dataset['dntps.txt']
        for BuildingBlock in dNTPs:
            self.dNTPs.append(BuildingBlock)
        self.dNTPKeys = ['A', 'C', 'G', 'T']
        self.N_dNTPs = len(self.dNTPs)

        NTPs = Dataset['ntps.txt']
        for BuildingBlock in NTPs:
            self.NTPs.append(BuildingBlock)
        self.NTPKeys = ['A', 'C', 'G', 'U']
        self.N_NTPs = len(self.NTPs)

        AAs = Dataset['amino_acids.txt']
        for BuildingBlock in AAs:
            self.AAs.append(BuildingBlock)
        self.N_AAs = len(self.AAs)

        AAKeys = Dataset['amino_acid_keys.txt']
        for BuildingBlock in AAKeys:
            self.AAKeys.append(BuildingBlock)

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)


class FMaster():
    def __init__(self):
        self.Switch4DebugMasterDataset = True

        self.MasterID2Index = {}
        self.MasterIDs = []
        self.MasterCounts = []
        self.MasterMWs = []

    def SaveData(self, SavePath):
        for Key, Value in self.__dict__.items():
            SaveFileName = "%s/%s" % (SavePath, Key)
            np.save('%s.npy' % SaveFileName, Value)
