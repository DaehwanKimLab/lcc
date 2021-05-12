import os, sys
import csv
import numpy as np
import ast
import re

class FCompilerData:
    def __init__(self):
        self.Switch4DebugCompilerData = True
        self.Switch4SaveAllData = True

        # Genome sequence
        self.InputGenomeSeq = []
        self.GenomeLengths = 0
        self.DNABasePairMW = 0

        # Lists of molecular species
        self.dNTPs = []
        self.NTPs = []
        self.AAs = []
        self.dNTPKeys = []
        self.NTPKeys = []
        self.AAKeys = []
        self.Compartments = []

        # Numbers of molecule species
        self.N_dNTPs = 0
        self.N_NTPs = 0
        self.N_AAs = 0
        self.N_Compartments = 0
        self.N_Chromosomes = 0
        self.N_Promoters = 0
        self.N_Genes = 0
        self.N_Metabolites = 0
        self.N_RNAs = 0
        self.N_mRNAs = 0
        self.N_tRNAs = 0
        self.N_rRNAs = 0
        self.N_miscRNAs = 0
        self.N_Proteins = 0
        self.N_Complexes = 0
        self.N_RXNs = 0
        self.N_KinRXNs = 0
        self.N_TCSMolecules = 0

        # Metabolites
        self.MetaboliteName2MWTemp = {}
        self.MetaboliteNames = []
        self.MetaboliteName2Index = {}
        self.MetaboliteCounts = 0
        self.MetaboliteMWs = 0

        # Chromosomes
        self.ChromosomeLengths = 0
        self.NMax_Chromosomes = 0
        self.ChromosomeIDs = []
        self.ChromosomeBPCounts = 0
        self.DNABasePairMWs = 0

        # Genes
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

        # Promoters
        self.PromoterCoordinates = 0
        self.PromoterDirections = []
        self.PromoterIDs = []
        self.PromoterID2Index = {}
        self.PromoterTargetGenesSymbols = []

        # RNAs
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

        self.RNATypeIndex4AllRNA = []
        self.RNATypeIndex4mRNA = []
        self.RNATypeIndex4tRNA = []
        self.RNATypeIndex4rRNA = []
        self.RNATypeIndex4miscRNA = []

        # Proteins
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

        # Complexes
        self.ComplexName2CountTemp = {}
        self.ComplexNameWithoutLocalization2CountTemp = {}
        self.ComplexNames = []
        self.ComplexName2ID = {}
        self.ComplexIDs = []
        self.ComplexID2Index = {}
        self.ComplexCounts = 0
        self.ComplexLengths = 0
        self.ComplexMWs = 0
        self.ComplexLocations = []

        # Reactions
        self.RXNIDs = []
        self.RXNID2Index = {}
        self.RXNStoichiometries = []
        self.RXNReversibilities = []
        self.RXNEnzymeIDs = []
        self.RXNSubstrates = []
        self.RXNSubstrateStoichs = []
        self.RXNProducts = []
        self.RXNProductStoichs = []

        self.KinRXNID2KineticData = {}
        self.KinRXNClassIDs = []
        self.KinEnzymeIDs = []
        self.KinSubstrateIDs = []
        self.KinTemperature = 0
        self.KinKcats = 0
        self.KinKms = []
        self.KinKis = 0

        # self.TCSRXNs = []
        self.TCSMolNames = []
        self.TCSMolIndexes = 0

        # Compartments
        self.CompartmentKeys = []
        self.CompartmentIDs = []
        self.CompartmentKey2Index = {}
        self.CompartmentIndexes = 0

        # Master Look Up Table
        self.MasterID2Index = {}
        self.MasterIDs = []
        self.MasterCounts = []
        self.MasterMWs = []
        # self.MasterLocalizations = 0

        # Data Path
        self.DataPath = None
        self.SavePath = None

    def SetDataPath(self, InDataPath):
        self.DataPath = InDataPath

    def GetDataPath(self):
        return self.DataPath

    def SetSavePath(self, InDataPath):
        self.SavePath = InDataPath

    def GetSavePath(self):
        return self.SavePath

    def LoadData(self, data_dir):
        dataset = dict()
        def parse_tsv(fpath, fname):
            fullpath = fpath + '/' + fname
            #print(fname)

            with open(fullpath) as fp:
                csv_reader = csv.reader(fp, delimiter = '\t')
                list_of_rows = list(csv_reader)

                dataset[fname] = list_of_rows[1:]

        def dump_dataset():
            for key, value in dataset.items():
                print(key, len(value))

        for fname in os.listdir(data_dir):
            if fname.endswith('.tsv'):
                parse_tsv(data_dir, fname)

        for fname in os.listdir(data_dir + '/wcs_simdata'):
            if fname.endswith('.tsv'):
                parse_tsv(data_dir + '/wcs_simdata', fname)

        if self.Switch4DebugCompilerData:
            dump_dataset()
        return dataset

    def RemoveLocalizationFromMolName(self, VariableName):
        return VariableName[:-3]

    def GetVarName(self, VariableName):
        for Key, Value in self.__dict__.items():
            if id(VariableName) == id(Value):
                return Key

    def SaveToNpy(self, VariableName):
        if self.Switch4SaveAllData:
            NameStr = self.GetVarName(VariableName)
            np.save("%s.npy" % NameStr, VariableName)
        return

    def SetUpData(self, Dataset):
        self.SetUpData_Metabolites(Dataset)
        self.SetUpData_Genomes()
        self.SetUpData_Chromosomes(Dataset)
        self.SetUpData_Genes(Dataset)
        self.SetUpData_Promoters(Dataset)
        self.SetUpData_RNAs(Dataset)
        self.SetUpData_Proteins(Dataset)
        self.SetUpData_Complexes(Dataset)
        self.SetUpData_RXNs(Dataset)
        self.SetUpData_TwoComponentSystems(Dataset)
        self.SetUpData_Compartments(Dataset)
        self.SetUpData_BuildingBlocks()
        self.SetUpData_MasterDataset()
        return

    def SetUpData_Metabolites(self, Dataset):
        MetaboliteMWs = Dataset['metabolites.tsv']
        WaterMW = Dataset['water.tsv'][0]
        MetaboliteMWs.append(WaterMW)
        for Value in MetaboliteMWs:
            Name, MW, Localization = Value
            assert Name not in self.MetaboliteName2MWTemp
            self.MetaboliteName2MWTemp[Name] = MW

        MetaboliteCounts = Dataset['Counts4AllMetabolites_DL_2_totalIs_0_0_Metabolism_BulkMolecules.tsv']
        self.N_Metabolites = len(MetaboliteCounts)
        self.MetaboliteCounts = np.zeros(self.N_Metabolites)
        self.MetaboliteMWs = np.zeros(self.N_Metabolites)
        for i, Value in enumerate(MetaboliteCounts):
            Name, Count, Request = Value
            NameLocRemoved = self.RemoveLocalizationFromMolName(Name)
            assert Name not in self.MetaboliteName2Index
            self.MetaboliteName2Index[Name] = len(self.MetaboliteNames)  # = i
            self.MetaboliteNames.append(Name)
            self.MetaboliteCounts[i] = Count
            self.MetaboliteMWs[i] = self.MetaboliteName2MWTemp[NameLocRemoved]
            assert self.MetaboliteMWs[i] != 0

        # self.SaveToNpy(self.MetaboliteNames)
        self.SaveToNpy(self.MetaboliteCounts)
        self.SaveToNpy(self.MetaboliteMWs)

        return

    # Genomes
    def SetUpData_Genomes(self):
        self.GenomeLengths = len(self.InputGenomeSeq)
        # self.N_AsInGenome = self.InputGenomeSeq.count('A')
        # self.N_CsInGenome = self.InputGenomeSeq.count('C')
        # self.N_GsInGenome = self.InputGenomeSeq.count('G')
        # self.N_TsInGenome = self.InputGenomeSeq.count('T')
        return

    def SetUpData_Chromosomes(self, Dataset):
        # Number of Chromosomes by Input Genome parsing (default = 1). Read from FASTA
        self.N_Chromosomes = 1
        # Set up an arbitrary max number to keep NT Count arrays for all chromosomes
        MaxCopyNumberOfEachChromosome = 3
        self.NMax_Chromosomes = self.N_Chromosomes * MaxCopyNumberOfEachChromosome
        self.ChromosomeLengths = np.zeros(self.NMax_Chromosomes)
        self.ChromosomeBPCounts = np.zeros(self.NMax_Chromosomes)
        self.DNABasePairMWs = np.zeros(self.NMax_Chromosomes)
        DNABasePairMW = 650
        for i in range(self.NMax_Chromosomes):
            for j in range(self.N_Chromosomes):
                ChromosomeNumber = j + 1
                ChromosomeReplication = i + 1
                ChromosomeID = 'Ch%d_%d' % (ChromosomeNumber, ChromosomeReplication)
                self.ChromosomeIDs.append(ChromosomeID)
                if i == 0:
                    self.ChromosomeBPCounts[j] = self.GenomeLengths
                self.DNABasePairMWs[i] = DNABasePairMW
        # average molecular weight of one base pair: 650 Daltons, which is 650 g/mol
        return

    def SetUpData_Genes(self, Dataset):
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
        self.GeneCounts = np.ones(self.N_Chromosomes)

    # Promoters
    def SetUpData_Promoters(self, Dataset):
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

            if self.Switch4DebugCompilerData:
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
                    print("Promoter %s contains an irregular name format to parse: %s" % (PromoterID, Name))
                    IrregularNames += 1
                    TargetGeneSymbol = Name
            # Check if the gene is expressed by multiple promoters
            if self.Switch4DebugCompilerData:
                if TargetGeneSymbol in self.PromoterTargetGenesSymbols:
                    print('The gene symbol %s is expressed by multiple promoters' % TargetGeneSymbol)
            self.PromoterTargetGenesSymbols.append(TargetGeneSymbol)
        if self.Switch4DebugCompilerData:
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
        return

    # RNAs
    def SetUpData_RNAs(self, Dataset):
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

        return

    # Protein Monomers
    def SetUpData_Proteins(self, Dataset):
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
        return


    # Complexes
    def SetUpData_Complexes(self, Dataset):
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

        Complexes = Dataset['proteinComplexes_large.tsv']
        self.N_Complexes = len(Complexes)
        self.ComplexLengths = np.zeros(self.N_Complexes)
        self.ComplexMWs = np.zeros(self.N_Complexes)
        self.ComplexCounts = np.zeros(self.N_Complexes)
        for i, Value in enumerate(Complexes):
            Name, Comments, MW, Location, RxnID, ComplexID = Value
            self.ComplexNames.append(Name)
            self.ComplexName2ID[Name] = ComplexID
            self.ComplexID2Index[ComplexID] = len(self.ComplexIDs)
            self.ComplexIDs.append(ComplexID)
            self.ComplexMWs[i] = max(np.fromstring(MW[1:-1], dtype='float32', sep=','))
            assert self.ComplexMWs[i] != 0
            self.ComplexLocations.append(Location)
            if ComplexID in self.ComplexNameWithoutLocalization2CountTemp:
                Count = self.ComplexNameWithoutLocalization2CountTemp[ComplexID]
                self.ComplexCounts[i] = Count

    # Reactions and kinetic parameters (K values)
    def SetUpData_RXNs(self, Dataset):
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
            self.KinRXNID2KineticData[ReactionID] = len(self.RXNEnzymeIDs)
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
                if self.Switch4DebugCompilerData:
                    print('TempAbscent:', ReactionID)

            kcat = kcat.replace(',', '').replace('[', '').replace(']', '')
            if kcat:
                self.KinKcats[i] = float(kcat)
            else:
                N_KcatValueAbscent += 1
                if self.Switch4DebugCompilerData:
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
                    if self.Switch4DebugCompilerData:
                        print(ReactionID, ': The number of kM values do not match the number of substrates')
                    self.KinKms.append(0.0)
            else:
                if kM:
                    self.KinKms.append(float(kM))
                else:
                    self.KinKms.append(0.0)
                    N_KmValueAbscent += 1
                    if self.Switch4DebugCompilerData:
                        print('KmValueAbscent:', ReactionID)

            kI = kI.replace(',', '').replace('[', '').replace(']', '')
            if kI:
                self.KinKis[i] = float(kI)
            else:
                N_KiValueAbscent += 1
                if self.Switch4DebugCompilerData:
                    print('KiValueAbscent:', ReactionID)
        Km_ExportShort = self.N_KinRXNs - len(self.KinKms)
        assert len(self.KinKms) == self.N_KinRXNs, '%s Km exporting fell short' % Km_ExportShort
        if self.Switch4DebugCompilerData:
            print('N_Total:', self.N_KinRXNs)
            print('N_TempAbscent:', N_TempAbscent)
            print('N_KcatValueAbscent:', N_KcatValueAbscent)
            print('N_KmValueAbscent:', N_KmValueAbscent)
            print('N_KiValueAbscent:', N_KiValueAbscent)
            print('N_KmValueUnmatched:', N_KmValueUnmatched)

        self.SaveToNpy(self.KinRXNClassIDs)
        self.SaveToNpy(self.KinEnzymeIDs)
        self.SaveToNpy(self.KinSubstrateIDs)
        self.SaveToNpy(self.KinTemperature)
        self.SaveToNpy(self.KinKcats)
        self.SaveToNpy(self.KinKms)
        self.SaveToNpy(self.KinKis)
        self.SaveToNpy(self.N_KinRXNs)

        # # Adjust kcat for temperature
        # temperature = constraint["Temp"]
        #
        # # If temperature not reported, assume 25 C
        # if type(temperature) == str:
        #     temperature = 25
        # constraint["kcatAdjusted"] = 2 ** ((37. - temperature) / 10.) * constraint["kcat"]

    def SetUpData_TwoComponentSystems(self, Dataset): # TCS for short
        # TwoComponentSystems = Dataset['twoComponentSystems.tsv']
        TwoComponentSystems = Dataset['TwoComponentSystemsTemporary_DL.tsv'] # temporary data table
        # Get the ID and Index list of molecule involved in TCS
        self.N_TCSMolecules = len(TwoComponentSystems)
        self.TCSMolIndexes = np.zeros(self.N_TCSMolecules)
        for i, Value in enumerate(TwoComponentSystems):
            Name, Count = Value
            self.TCSMolNames.append(Name)
            # self.TCSMolIndexes[i] = self.ProteinID2Index[Name]
        # Get the index list of TCS molecules
        return

    # CompartmentList
    def SetUpData_Compartments(self, Dataset):
        Compartments = Dataset['compartments.tsv']
        self.N_Compartments = len(Compartments)
        self.CompartmentIndexes = np.zeros(self.N_Compartments)
        for i, Value in enumerate(Compartments):
            Abbrev, ID = Value
            self.CompartmentKey2Index[Abbrev] = len(self.CompartmentKeys)
            self.CompartmentKeys.append(Abbrev)
            self.CompartmentIDs.append(ID)
        return

    # BuildingBlockLists
    def SetUpData_BuildingBlocks(self):
        TxtFilePath = self.DataPath + '/'
        with open(TxtFilePath + "dntps.txt", 'r') as OpenFile:
            for line in OpenFile:
                self.dNTPs.append(line[:-1])
        assert '\n' not in self.dNTPs
        self.dNTPKeys = ['A', 'C', 'G', 'T']
        self.N_dNTPs = len(self.dNTPs)

        with open(TxtFilePath + "ntps.txt", 'r') as OpenFile:
            for line in OpenFile:
                self.NTPs.append(line[:-1])
        assert '\n' not in self.NTPs
        self.NTPKeys = ['A', 'C', 'G', 'U']
        self.N_NTPs = len(self.NTPs)

        with open(TxtFilePath + "amino_acids.txt", 'r') as OpenFile:
            for line in OpenFile:
                self.AAs.append(line[:-1])
        assert '\n' not in self.AAs
        self.N_AAs = len(self.AAs)

        with open(TxtFilePath + "amino_acid_keys.txt", 'r') as OpenFile:
            for line in OpenFile:
                self.AAKeys.append(line[:-1])
        assert '\n' not in self.AAKeys
        return

    def SaveAllCompilerData(self):
        if self.Switch4SaveAllData:
            SavePath = 'lccsave'
            self.SetSavePath(SavePath)
            for Key, Value in self.__dict__.items():
                SaveFileName = "%s/%s" % (self.SavePath, Key)
                np.save('%s.npy' % SaveFileName, Value)
        return

    def AddToMaster(self, IDs, Counts, MWs, Localizations='None'):
        for ID in IDs:
            self.MasterID2Index[ID] = len(IDs)
            self.MasterIDs.append(ID)
        self.MasterCounts = np.append(self.MasterCounts, Counts)
        self.MasterMWs = np.append(self.MasterMWs, MWs)
        # if Localizations:
        #     self.MasterLocalizations = np.append(self.MasterLocalizations, Localizations)
        if self.Switch4DebugCompilerData:
            for Key, Value in self.__dict__.items():
                if Key[:6] == 'Master':
                    LengthOfValue = len(Value)
                    print('The length of "%s": %d' % (Key, LengthOfValue))
        assert len(self.MasterIDs) == len(self.MasterID2Index) and len(self.MasterIDs) == len(self.MasterCounts) and len(self.MasterIDs) == len(self.MasterMWs)

    def SetUpData_MasterDataset(self):
        # Concatenate IDs, Counts, MWs
        self.AddToMaster(self.MetaboliteNames, self.MetaboliteCounts, self.MetaboliteMWs)
        self.AddToMaster(self.ChromosomeIDs, self.ChromosomeBPCounts, self.DNABasePairMWs)
        self.AddToMaster(self.RNAIDs, self.RNACounts, self.RNAMWs)
        self.AddToMaster(self.ProteinIDs, self.ProteinCounts, self.ProteinMWs)
        self.AddToMaster(self.ComplexIDs, self.ComplexCounts, self.ComplexMWs)

        return