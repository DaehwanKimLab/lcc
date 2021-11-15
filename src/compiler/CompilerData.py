import os, sys
import csv
import numpy as np
import ast
from Dataset import *


class FCompilerData:
    def __init__(self):
        self.InputGenomeSeq = list()

        self.Switch4DebugCompilerData = False
        self.Switch4SaveAllData = True
        self.Switch4SaveProcessedData = False

        self.Dict_DataClass = dict()

        self.Dict4ID_Gene = dict()
        self.Dict4ID_RNA = dict()
        self.Dict4ID_Protein = dict()

        # Data classes
        self.Genome = None
        self.Metabolite = None
        self.Chromosome = None
        self.Gene = None
        self.Promoter = None
        self.RNA = None
        self.Protein = None
        self.Complex = None
        self.RXN = None
        self.Compartment = None
        self.BuildingBlock = None
        self.Master = None

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

    def LoadRawData(self, data_dir):
        dataset = dict()
        def parse_tsv(fpath, fname):
            fullpath = fpath + '/' + fname
            #print(fname)

            with open(fullpath) as fp:
                csv_reader = csv.reader(fp, delimiter = '\t')
                list_of_rows = list(csv_reader)

                dataset[fname] = list_of_rows[1:]

        def parse_txt(fpath, fname):
            fullpath = fpath + '/' + fname
            #print(fname)

            with open(fullpath, 'r') as fp:
                list_of_rows = fp.read().splitlines()

                dataset[fname] = list_of_rows[:]

        def parse_fasta(fpath, fname):
            fullpath = fpath + '/' + fname
            fasta = dict()
            #print(fname)

            with open(fullpath, 'r') as fp:
                sequences = fp.read().split('>')[1:]
                for sequence in sequences:
                    chr = 'Ch%s' % (len(fasta) + 1)
                    list_of_rows = sequence.split('\n')
                    fasta[chr] = ''.join(list_of_rows[1:])
                dataset[fname] = fasta

        def dump_dataset():
            for key, value in dataset.items():
                print(key, len(value))

        for fname in os.listdir(data_dir):
            if fname.endswith('.tsv'):
                parse_tsv(data_dir, fname)

            if fname.endswith('.txt'):
                parse_txt(data_dir, fname)

            if fname.endswith('.fasta'):
                parse_fasta(data_dir, fname)

        for fname in os.listdir(data_dir + '/wcs_simdata'):
            if fname.endswith('.tsv'):
                parse_tsv(data_dir + '/wcs_simdata', fname)

        for fname in os.listdir(data_dir + '/intermediate'):
            if fname.endswith('.tsv'):
                parse_tsv(data_dir + '/intermediate', fname)

        if self.Switch4DebugCompilerData:
            dump_dataset()
        return dataset

    def InitializeCompilerData(self, Dataset):
        self.BuildingBlock = FBuildingBlock()
        self.Dict_DataClass['BuildingBlock'] = self.BuildingBlock

        self.Compartment = FCompartment()
        self.Dict_DataClass['Compartment'] = self.Compartment

        self.Chromosome = FChromosome()
        self.Dict_DataClass['Chromosome'] = self.Chromosome

        self.Gene = FGene()
        self.Dict_DataClass['Gene'] = self.Gene

        self.Promoter = FPromoter()
        self.Dict_DataClass['Promoter'] = self.Promoter

        self.RNA = FRNA()
        self.Dict_DataClass['RNA'] = self.RNA

        self.Protein = FProtein()
        self.Dict_DataClass['Protein'] = self.Protein

        self.Complex = FComplex()
        self.Dict_DataClass['Complex'] = self.Complex

        self.Metabolite = FMetabolite()
        self.Dict_DataClass['Metabolite'] = self.Metabolite

        self.Complexation = FComplexation()
        self.Dict_DataClass['Complexation'] = self.Complexation

        self.Equilibrium = FEquilibrium()
        self.Dict_DataClass['Equilibrium'] = self.Equilibrium

        self.Metabolism = FMetabolism()
        self.Dict_DataClass['Metabolism'] = self.Metabolism

        self.Kinetics = FKinetics()
        self.Dict_DataClass['Kinetics'] = self.Kinetics

        self.Transporters = FTransporter()
        self.Dict_DataClass['Transporters'] = self.Transporters

        self.UserInput = FUserInput()
        self.Dict_DataClass['UserInput'] = self.UserInput

        self.Master = FMaster()
        self.Dict_DataClass['Master'] = self.Master

        print('lcc compiler data have been initialized. [# of data classes: %s]' % len(self.Dict_DataClass))


    def SetUpCompilerData(self, Dataset):
        for DataClassName, DataClassDataset in self.Dict_DataClass.items():
            DataClassDataset.SetUpData(Dataset, self)
            print('lcc compiler data have been set up. Class: "%s"' % DataClassName)

        self.SetUpDict4FindMolID()


    def SaveCompilerData(self):
        # Save all compiler data
        if self.Switch4SaveAllData:
            for DataClassName, DataClassDataset in self.Dict_DataClass.items():
                DataClassDataset.SaveData(self.SavePath)
                print('lcc compiler data have been saved. Class: "%s"' % DataClassName)
        # Save processed compiler data only used in simulation (not implemented yet)
        # elif self.Switch4SaveProcessedData:
        #     for DataClassName, DataClassDataset in self.Dict_DataClass.items():
        #         DataClassDataset.SaveData(self.SavePath)
        #         print('Processed compiler data have been saved. Class: "%s"' % DataClassName)
        else:
            print('lcc compiler data save option is off')

    def SetUpDict4FindMolID(self):

        ListOfDictsToMerge = [
            self.Dict4ID_Gene,
            self.Dict4ID_RNA,
            self.Dict4ID_Protein,
        ]

        ListOfDicts4Gene = [
            self.Gene.Name2ID_Genes,
            self.Gene.Sym2ID_Genes,
            self.Protein.ID2ID_Protein2Gene,
            self.RNA.ID2ID_RNA2Gene,
        ]

        ListOfDicts4RNA = [
            self.RNA.Name2ID_RNAs,
            self.RNA.ID2ID_Gene2RNA,
            self.Protein.ID2ID_Protein2RNA,
            self.Master.Name2ID_Gene2RNA_Master,
            self.Master.Sym2ID_Gene2RNA_Master,
        ]

        ListOfDicts4Protein = [
            self.Protein.Name2ID_Proteins,
            self.Protein.ID2ID_Gene2Protein,
            self.Protein.ID2ID_RNA2Protein,
            self.Master.Name2ID_Gene2Protein_Master,
            self.Master.Sym2ID_Gene2Protein_Master,
        ]

        ListOfDicts4Type = [
            ListOfDicts4Gene,
            ListOfDicts4RNA,
            ListOfDicts4Protein,
        ]

        for DictToMerge, Dicts4Type in zip(ListOfDictsToMerge, ListOfDicts4Type):
            for Dict4Type in Dicts4Type:
                DictToMerge.update(Dict4Type)
