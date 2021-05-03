#!/usr/bin/env python3
#
# Copyright 2021,
# Donghoon Lee <dhldjl@gmail.com>,
# Chanhee Park <parkchanhee@gmail.com>, and
# Daehwan Kim <infphilo@gmail.com>,
#
# This file is part of LDP.
#
# LDP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

import os, sys
import numpy as np
import csv
# import tensorflow as tf
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import CodeGen
from CodeGen import Target
import inspect
from lccmodule import CellMX
from lccmodule import TwoComponentSystems
from lccmodule import TranscriptElongation
from lccmodule import RNADegradation
from lccmodule import Metabolism
from lccmodule import PolypeptideElongation

LCC_VERSION = "0.1"

def lcc_dummy():
    pass
LCC_PATH = os.path.dirname(os.path.realpath(inspect.getsourcefile(lcc_dummy)))


# Move out to a file
class FCompilerData:
    def __init__(self):
        self.dNTPs = []
        self.NTPs = []
        self.AAs = []
        self.dNTPsWithLoc = []
        self.NTPsWithLoc = []
        self.AAsWithLoc = []
        self.dNTPKeys = []
        self.NTPKeys = []
        self.AAKeys = []


        self.MetaboliteNames4Conc = []
        self.MetaboliteName2ConcIndex = {}
        self.MetaboliteConcs = None
        self.MetaboliteNames4MW = []
        self.MetaboliteName2MWIndex = {}
        self.MetaboliteMWs = None
        self.MetaboliteNames4Counts = []
        self.MetaboliteName2CountIndex = {}
        self.MetaboliteCounts = None


        self.RXNIDs = []
        self.RXNID2Index = {}
        self.RXNStoichiometries = []
        self.RXNReversibilities = []
        self.RXNEnzymes = []

        # self.TCSRXNs = []
        self.MolNames = []
        self.MolName2Index = {}
        self.MolCounts = None
        self.MolConcs = None

        self.RNAName2Index = {}
        self.RNANames = []
        self.RNAMonomerID2Index = {}
        self.RNAMonomerIDs = []
        self.RNAID2Index = {}
        self.RNAIDs = []
        self.RNAHalfLives = []
        self.RNASeqs = []
        self.RNALengths = []
        self.RNATypes = []
        self.RNAMolWeights = []
        self.RNANTCounts = []
        self.RNANTFreqs = []

        self.RNATypeIndex4AllRNA = []
        self.RNATypeIndex4mRNA = []
        self.RNATypeIndex4tRNA = []
        self.RNATypeIndex4rRNA = []
        self.RNATypeIndex4miscRNA = []

        self.ProtName2Index = {}
        self.ProtNames = []
        self.ProtID2Index = {}
        self.ProtIDs = []
        self.ProtSeqs = []
        self.ProtLengths = []
        self.ProtMolWeights = []
        self.ProtLocations = []
        self.ProtAACounts = []
        self.ProtAAFreqs = []
        self.ProtID2GeneID = {}
        self._GeneIDs = []
        self.ProtID2RNAID = {}
        self._RNAIDs = []


        # Data Path
        self.DataPath = LCC_PATH

    def SetDataPath(self, InDataPath):
        self.DataPath = InDataPath

    def GetDataPath(self):
        return self.DataPath

def WriteLicense(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)
    for line in open(LCC_PATH + "/LICENSE.input"):
        line = line.strip()
        Writer.Statement(line)
    Writer.SetIndentLevel(tmpLevel)


def WriteImport(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)

    Writer.Statement("import os, sys")
    Writer.Statement("import numpy as np")
    Writer.Statement("import tensorflow as tf")
    Writer.Statement("import matplotlib.pyplot as plt")
    Writer.Statement("from datetime import datetime")

    Writer.Statement("from argparse import ArgumentParser, FileType")

    Writer.SetIndentLevel(tmpLevel)


def LoadData(data_dir):
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

    for fname in os.listdir(data_dir + '/MolCount_Sim1'):
        if fname.endswith('.tsv'):
            parse_tsv(data_dir + '/MolCount_Sim1', fname)


    dump_dataset()
    return dataset


def SetUpCompilerData(Dataset, CompilerData):

    def SetUpDataIndexes():
        # Saving npy files may be done in a separate temporary compiler data output folder

        # This table is to be discarded and completely replaced by MetaboliteCounts Dataset
        # MetaboliteMWs = Dataset['metabolites.tsv']
        # CompilerData.MetaboliteMWs = np.zeros(len(MetaboliteMWs))
        # for i, Value in enumerate(MetaboliteMWs):
        #     Name, MW, Localization = Value
        #     assert Name not in CompilerData.MetaboliteName2MWIndex
        #     CompilerData.MetaboliteName2MWIndex[Name] = len(CompilerData.MetaboliteNames4MW) # = i
        #     CompilerData.MetaboliteNames4MW.append(Name)
        #     CompilerData.MetaboliteMWs[i] = MW
        # np.save("MetaboliteMWs", CompilerData.MetaboliteMWs)

        MetaboliteConcs = Dataset['metaboliteConcentrations.tsv']
        CompilerData.MetaboliteConcs = np.zeros(len(MetaboliteConcs))
        for i, Value in enumerate(MetaboliteConcs):
            Name, Conc = Value
            assert Name not in CompilerData.MetaboliteName2ConcIndex
            CompilerData.MetaboliteName2ConcIndex[Name] = len(CompilerData.MetaboliteNames4Conc) # = i
            CompilerData.MetaboliteNames4Conc.append(Name)
            CompilerData.MetaboliteConcs[i] = Conc
        np.save("MetaboliteConcs", CompilerData.MetaboliteConcs)

        MetaboliteCounts = Dataset['DL_0_Metabolism_BulkMolecules.tsv']
        CompilerData.MetaboliteCounts = np.zeros(len(MetaboliteCounts))
        for i, Value in enumerate(MetaboliteCounts):
            Name, Count = Value
            assert Name not in CompilerData.MetaboliteName2ConcIndex
            CompilerData.MetaboliteName2CountIndex[Name] = len(CompilerData.MetaboliteNames4Counts)  # = i
            CompilerData.MetaboliteNames4Counts.append(Name)
            CompilerData.MetaboliteCounts[i] = Count
        np.save("MetaboliteCounts", CompilerData.MetaboliteCounts)

        RXNs = Dataset['reactions.tsv']
        for i, Value in enumerate(RXNs):
            RXNID, RXNStoichiometry, RXNReversibility, RXNEnzyme = Value
            assert RXNID not in CompilerData.RXNID2Index
            CompilerData.RXNID2Index[RXNID] = len(CompilerData.RXNIDs) # = i
            CompilerData.RXNIDs.append(RXNID)
            CompilerData.RXNStoichiometries.append(RXNStoichiometry)
            CompilerData.RXNReversibilities.append(RXNReversibility)
            CompilerData.RXNEnzymes.append(RXNEnzyme)

        # TwoComponentSystems = Dataset['twoComponentSystems.tsv']
        MolCounts = Dataset['TwoComponentSystemsTemporary_DL.tsv'] # temporary data table
        CompilerData.MolCounts = np.zeros(len(MolCounts))
        for i, Value in enumerate(MolCounts):
            MolName, MolCount = Value
            CompilerData.MolName2Index[MolName] = len(CompilerData.MolNames) # = i
            CompilerData.MolNames.append(MolName)
            CompilerData.MolCounts[i] = MolCount
        np.save("MolCounts.npy", CompilerData.MolCounts)

        # Temporary command to save a list of participating molecules in TCS
        np.save('TCSMolNames.npy', CompilerData.MolNames)

        RNAs = Dataset['rnas.tsv']
        for i, Value in enumerate(RNAs):
            HalfLife, Name, Seq, Type, ModifiedForms, MonomerID, Comments, MolWeight, Location, NTCount, ID, GeneID, MicArrExp = Value
            assert Name not in CompilerData.RNAID2Index
            CompilerData.RNAName2Index[Name] = len(CompilerData.RNANames) # = i
            CompilerData.RNANames.append(Name)
            CompilerData.RNAMonomerID2Index[MonomerID] = len(CompilerData.RNAMonomerIDs)
            CompilerData.RNAMonomerIDs.append(MonomerID)
            CompilerData.RNAID2Index[ID] = len(CompilerData.RNAIDs)
            CompilerData.RNAIDs.append(ID)
            CompilerData.RNAHalfLives.append(HalfLife)
            CompilerData.RNASeqs.append(Seq)
            CompilerData.RNALengths.append(len(Seq))
            CompilerData.RNATypes.append(Type)
            CompilerData.RNAMolWeights.append(MolWeight[1:-1].split(',')[5])
            CompilerData.RNANTCounts.append(NTCount)
            CompilerData.RNATypeIndex4AllRNA.append(i)
            if Type == 'mRNA':
                CompilerData.RNATypeIndex4mRNA.append(i)
            elif Type == 'tRNA':
                CompilerData.RNATypeIndex4tRNA.append(i)
            elif Type == 'rRNA':
                CompilerData.RNATypeIndex4rRNA.append(i)
            elif Type == 'miscRNA':
                CompilerData.RNATypeIndex4miscRNA.append(i)
            else:
                print('Warning: Unaccounted RNA type detected: ', Type)
        np.save('RNAIDs.npy', CompilerData.RNAIDs)
        np.save('RNATypeIndex4AllRNA.npy', CompilerData.RNATypeIndex4AllRNA)
        np.save('RNATypeIndex4mRNA.npy', CompilerData.RNATypeIndex4mRNA)
        np.save('RNATypeIndex4tRNA.npy', CompilerData.RNATypeIndex4tRNA)
        np.save('RNATypeIndex4rRNA.npy', CompilerData.RNATypeIndex4rRNA)
        np.save('RNATypeIndex4miscRNA.npy', CompilerData.RNATypeIndex4miscRNA)

        Proteins = Dataset['proteins.tsv']
        SameProtName = 0
        SameProtNameSeq = 0
        SameProtNameSeqLocation = 0
        SameProtNameLocation = 0
        SameProtID = 0
        SameGeneID = 0
        SameRNAID = 0
        for i, Value in enumerate(Proteins):
            AACount, Name, Seq, Comments, CodingRNASeq, MolWeight, Location, RNAID, ProtID, GeneID = Value
            if Name in CompilerData.ProtName2Index:
                SameProtName += 1
                # print(Name)
                # print(Value)
                if Seq in CompilerData.ProtSeqs:
                    SameProtNameSeq += 1
                    # print(Name)
                    if Location == CompilerData.ProtLocations[CompilerData.ProtName2Index[Name]]:
                        SameProtNameSeqLocation += 1
                        # print("Same Protein Name, Seq, Location", Name)
                if Location == CompilerData.ProtLocations[CompilerData.ProtName2Index[Name]]:
                    SameProtNameLocation += 1
            # assert Name not in CompilerData.ProtName2Index
            CompilerData.ProtName2Index[Name] = len(CompilerData.ProtNames) # = i
            CompilerData.ProtNames.append(Name)
            if ProtID in CompilerData.ProtIDs:
                SameProtID += 1
                # print("Same ProtID:", ProtID)
            CompilerData.ProtID2Index[ProtID] = len(CompilerData.ProtIDs)
            CompilerData.ProtIDs.append(ProtID)
            CompilerData.ProtSeqs.append(Seq)
            CompilerData.ProtLengths.append(len(Seq))
            CompilerData.ProtMolWeights.append(MolWeight[1:-1].split(',')[6])
            CompilerData.ProtLocations.append(Location)
            CompilerData.ProtAACounts.append(AACount)
            if GeneID in CompilerData._GeneIDs:
                SameGeneID += 1
            CompilerData.ProtID2GeneID[ProtID] = GeneID
            if RNAID in CompilerData._RNAIDs:
                SameRNAID += 1
            CompilerData.ProtID2RNAID[ProtID] = RNAID
        np.save('ProtIDs.npy', CompilerData.ProtIDs)
        np.save('ProtMolWeights.npy', CompilerData.ProtMolWeights)
        np.save('ProtLocations.npy', CompilerData.ProtLocations)
        np.save('PRotAACounts.npy', CompilerData.ProtAACounts)

        # print('Same Name Count: ', SameProtName)
        # print('Same Name and Seq Count: ', SameProtNameSeq)
        # print('Same Name and Location Count: ', SameProtNameLocation)
        # print('Same Name, Seq and Location Count: ', SameProtNameSeqLocation)
        # print('Same Protein ID Count: ', SameProtID)
        # print('Same RNA ID Count: ', SameRNAID)
        # print('Same Gene ID Count: ', SameGeneID)
        return

    # Add [c] for each for reading in for count calling
    def GetBuildingBlockLists():
        # BuildingBlocks4TxtFiles = ['dntps', 'ntps', 'amino_acids']
        # BuildingBlocks4Variables = ['dNTPs', 'NTPs', 'AAs']
        TxtFilePath = CompilerData.GetDataPath() + '/'
        with open(TxtFilePath + "dntps.txt", 'r') as OpenFile:
            for line in OpenFile:
                CompilerData.dNTPs.append(line[:-1])
                CompilerData.dNTPsWithLoc.append(line[:-1] + '[c]')
        assert '\n' not in CompilerData.dNTPs
        CompilerData.dNTPKeys = ['A', 'C', 'G', 'T']

        with open(TxtFilePath + "ntps.txt", 'r') as OpenFile:
            for line in OpenFile:
                CompilerData.NTPs.append(line[:-1])
                CompilerData.NTPsWithLoc.append(line[:-1] + '[c]')
        assert '\n' not in CompilerData.NTPs
        CompilerData.NTPKeys = ['A', 'C', 'G', 'U']

        with open(TxtFilePath + "amino_acids.txt", 'r') as OpenFile:
            for line in OpenFile:
                CompilerData.AAs.append(line[:-1])
                CompilerData.AAsWithLoc.append(line[:-1] + '[c]')
        assert '\n' not in CompilerData.AAs

        with open(TxtFilePath + "amino_acid_keys.txt", 'r') as OpenFile:
            for line in OpenFile:
                CompilerData.AAKeys.append(line[:-1])
        assert '\n' not in CompilerData.AAKeys



    def GetMXRNANTFreq():
        for RNASeq in CompilerData.RNASeqs:
            NTTotalCount = len(RNASeq)
            NTCounts = np.array([RNASeq.count("A"), RNASeq.count("C"), RNASeq.count("G"), RNASeq.count("U")])
            if np.sum(NTCounts) != NTTotalCount:
                print("WARNING: RNA seq may contain non-ACGU text.", file=sys.stderr)
            NTFreq = NTCounts / NTTotalCount
            CompilerData.RNANTCounts.append(NTCounts)
            CompilerData.RNANTFreqs.append(NTFreq)
        np.save("RNANTFreqs.npy", CompilerData.RNANTFreqs)
        np.save("RNALengths.npy", CompilerData.RNALengths)

    def GetMXProtAAFreq():
        for ProtSeq, LineForAACounts in zip(CompilerData.ProtSeqs, CompilerData.ProtAACounts):
            AATotalCount = len(ProtSeq)
            if AATotalCount == CompilerData.ProtLengths:
                print("WARNING: Protein length data may contain error.", file=sys.stderr)
            AACounts = np.array(list(map(int, LineForAACounts[1:-1].split(','))))
            AAFreq = AACounts / AATotalCount
            CompilerData.ProtAAFreqs.append(AAFreq)
        np.save("ProtAAFreqs.npy", CompilerData.ProtAAFreqs)
        np.save("ProtLengths.npy", CompilerData.ProtLengths)

    SetUpDataIndexes()
    GetBuildingBlockLists()
    GetMXRNANTFreq()
    GetMXProtAAFreq()
    # CompilerData.TCSMolNames = TwoComponentSystems()


def WriteBody(Writer, CompilerData):

    Writer.DebugPrintSwitch = False # ON/OFF switch to print data being processed
    Writer.DebugAssertSwitch = True


    Writer.Variable_('LCCDataPath', "\"" + CompilerData.GetDataPath() + "\"")
    Lines = [
        "def main(GenomeFileName, verbose):",
    ]
    for Line in Lines:
        Writer.Statement(Line)

    with Writer:
        # Define simulation parameters
        Writer.Statement("# Define simulation parameters")
        Writer.Variable_("CellCycles", 1)
        Writer.Variable_("SimulationSteps", 100)

        Writer.BlankLine()

        Writer.Statement("# This is a numpy code", TargetCode=Target.Numpy)
        Writer.Statement("# This is a tensorflow code", TargetCode=Target.TensorFlow)

        # Define key constants
        Writer.Statement("# Define key constants")
        Writer.Variable_("AvogadroNum", 6.022141527E23)

        # Define accessory variables - TF version only
        Writer.InitArrayWithOne('OneTF', 1, 'int32')
        Writer.InitArrayWithZero('ZeroTF', 1)
        Writer.BlankLine()

        # GLOBAL INIT?

        # Load CellMX
        CellMX.Write_CellMX_Init(Writer)
        Writer.Statement("CellMX = FCellMX()")
        Writer.BlankLine()

        # Temporary E coli cell volume: 0.7 um3 (average), which is 7e-16 liters
        Writer.Variable_("CellMX.CellVol", 7e-16)  # TO BE REPLACED AND MOVED INTO SIMULATION

        # Set up commonly used Cell State and index matrices - Make a separate py file
        Writer.Statement("# Load metabolite concentration table")
        Writer.Statement("CellMX.MetaboliteConcs = np.load(\"MetaboliteConcs.npy\").astype('float32')")
        Writer.Statement("CellMX.MetaboliteConcsTF = tf.convert_to_tensor(CellMX.MetaboliteConcs)")
        Writer.DebugPVar("CellMX.MetaboliteConcs")
        Writer.DebugPVar("CellMX.MetaboliteConcsTF")
        Writer.BlankLine()

        Writer.Statement("# Load metabolite counts table")
        Writer.Statement("CellMX.MetaboliteCounts = np.load(\"MetaboliteCounts.npy\").astype('float32')")
        Writer.Statement("CellMX.MetaboliteCountsTF = tf.convert_to_tensor(CellMX.MetaboliteCounts)")
        Writer.DebugPVar("CellMX.MetaboliteCounts")
        Writer.DebugPVar("CellMX.MetaboliteCountsTF")

        Writer.Statement("# Load molecule count table")
        Writer.Statement("MolCounts = np.load(\"MolCounts.npy\").astype('float32')")
        Writer.Statement("CellMX.MolCountsTF = tf.convert_to_tensor(MolCounts)")
        Writer.DebugPVar("MolCounts")
        Writer.DebugPVar("CellMX.MolCountsTF")
        Writer.BlankLine()

        # Load all RNA counts (placeholder)
        Writer.Statement("# RNA Counts")
        Writer.Variable_("DefaultCount", 100)
        Writer.Statement("RNACounts = np.ones(" + str(len(CompilerData.RNAIDs)) + ").astype('int32') * DefaultCount")
        Writer.Statement("CellMX.RNACountsTF = tf.convert_to_tensor(RNACounts)")
        Writer.DebugPVar("CellMX.RNACountsTF")
        Writer.BlankLine()

        # Load all RNA Types indexes
        Writer.Statement("# Indices for RNA Types")
        Writer.Statement("CellMX.RNAIndex4mRNATF = np.load(\"RNATypeIndex4mRNA.npy\").astype('int32')")
        Writer.Statement("CellMX.RNAIndex4tRNATF = np.load(\"RNATypeIndex4tRNA.npy\").astype('int32')")
        Writer.Statement("CellMX.RNAIndex4rRNATF = np.load(\"RNATypeIndex4rRNA.npy\").astype('int32')")
        Writer.Statement("CellMX.RNAIndex4miscRNATF = np.load(\"RNATypeIndex4miscRNA.npy\").astype('int32')")
        Writer.BlankLine()

        # Load the number of unique RNAs
        Writer.Statement("CellMX.NumberOfUniqueRNAs = len(RNACounts)")
        Writer.Statement("CellMX.NumberOfUniquemRNAs = len(CellMX.RNAIndex4mRNATF)")
        Writer.Statement("CellMX.NumberOfUniquetRNAs = len(CellMX.RNAIndex4tRNATF)")
        Writer.Statement("CellMX.NumberOfUniquerRNAs = len(CellMX.RNAIndex4rRNATF)")
        Writer.Statement("CellMX.NumberOfUniquemiscRNAs = len(CellMX.RNAIndex4miscRNATF)")
        Writer.BlankLine()

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


        # Load initialization functions for each process
        TranscriptElongation.Write_TE_Init(Writer, CompilerData)
        TwoComponentSystems.Write_TCS_Init(Writer, CompilerData)
        RNADegradation.Write_RNADeg_Init(Writer, CompilerData)
        PolypeptideElongation.Write_PE_Init(Writer, CompilerData)
        Metabolism.Write_Metab_Init(Writer, CompilerData)


        # Load loop functions for each process
        TranscriptElongation.Write_TE_Loop(Writer)
        TwoComponentSystems.Write_TCS_Loop(Writer)
        RNADegradation.Write_RNADeg_Loop(Writer)
        PolypeptideElongation.Write_PE_Loop(Writer)
        Metabolism.Write_Metab_Loop(Writer)

        # Run initialization functions for each process
        Writer.Statement("TE_Init()")
        Writer.Statement("TCS_Init()")
        Writer.Statement("RNADeg_Init()")
        Writer.Statement("PE_Init()")
        Writer.Statement("Metab_Init()")


        Writer.BlankLine()

        # Run simulation
        Writer.Statement("# Run simulation")
        Writer.PrintStrg("Simulation begins...")
        with Writer.Statement("for SimulationStep in range(SimulationSteps):"):
            Writer.PrintStrg('=============================================')
            Writer.PrintStVa('SimulationStep: ', "SimulationStep + 1")
            Writer.Statement("CellMX.SimStep.append(SimulationStep + 1)")
            Writer.BlankLine()

            # Run loop functions for each process
            Writer.Statement("TE_Loop()")
            Writer.Statement("TCS_Loop()")
            Writer.Statement("RNADeg_Loop()")
            Writer.Statement("PE_Loop()")
            Writer.Statement("Metab_Loop()")
            Writer.BlankLine()

        # Temporary TE visualization code
        Writer.Statement("# Temporary TE visualization code")
        Writer.Statement("fig, ax = plt.subplots()")
        # Writer.Statement("ax.plot(CellMX.SimStep, CellMX.TE_ACGU)")
        Writer.Statement("lines = ax.plot(CellMX.SimStep, CellMX.TE_ACGU)")
        Writer.Statement("labels = ['ATP', 'CTP', 'GTP', 'UTP']")
        Writer.Statement("ax.legend(lines, labels)")
        Writer.Statement("ax.set(xlabel='SimStep', ylabel='Concentration (M)', title='NTP level')")
        Writer.Statement("ax.grid()")
        Writer.Statement("plt.show()")

        Writer.BlankLine()
        # End of simulation

        # Print input genome
        with Writer.Statement("if GenomeFileName != \"\":"):
            Writer.Statement("GenomeFile = open(GenomeFileName, 'w')")
            Writer.Statement("InputGenomeFile = open('cell.fa')")
            with Writer.Statement("for Line in InputGenomeFile:"):
                Writer.Statement("Line = Line.strip()")
                Writer.Statement("print(Line, file=GenomeFile)")
            Writer.Statement("GenomeFile.close()")

def WriteMain(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)

    Writer.BlankLine()
    with Writer.Statement("if __name__ == '__main__':"):
        Writer.Statement("parser = ArgumentParser(description='')")
        Writer.Statement("parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='also print some statistics to stderr')")
        Writer.Statement("parser.add_argument('-g', '--genome', dest='GenomeFileName', type=str, default='', help='')")
        Writer.Statement("parser.add_argument('-e', '--eager', dest='TFEagerMode', action='store_true', default=False, help='Running in the eager mode')")
        Writer.Statement("parser.add_argument('--trace-graph', dest='TFTraceGraph', action='store_true', default=False, help='Tracing the TensorFlow Graph')")
        Writer.BlankLine()
        Writer.Statement("args = parser.parse_args()")
        Writer.BlankLine()
        Writer.Statement("print('TF Eager Mode:', args.TFEagerMode)", TargetCode=Target.TensorFlow)
        Writer.Statement("tf.config.run_functions_eagerly(args.TFEagerMode)", TargetCode=Target.TensorFlow)
        Writer.BlankLine()
        with Writer.Statement("if args.TFTraceGraph:"):
            Writer.Statement("tf.summary.trace_on(graph=True)  # , profiler=True)")
        Writer.BlankLine()
        Writer.Statement("main(args.GenomeFileName, args.verbose)")
        Writer.BlankLine()

        with Writer.Statement("if args.TFTraceGraph:"):
            Writer.Statement("timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')")
            Writer.Statement("tf_logdir = 'logs/graph/%s' % timestamp")
            Writer.Statement("summary_writer = tf.summary.create_file_writer(tf_logdir)")
            with Writer.Statement("with summary_writer.as_default():"):
                Writer.Statement("tf.summary.trace_export(name='lcc_trace', step=0, profiler_outdir=tf_logdir)")
        Writer.BlankLine()

    Writer.SetIndentLevel(tmpLevel)


def NewLine(code_file):
    print("\t", file=code_file)


def Parse(CodeFileNames):
    Result = {}

    for CodeFileName in CodeFileNames:
        CodeFile = open(CodeFileName)
        for line in CodeFile:
            line = line.strip()
            if line.startswith('#'):
                continue

            if line.startswith('TemplateOrganism'):
                Organism = line.split(':')[1].strip()
                Result['Organism'] = Organism
            elif line.startswith('RNATranscription'):
                if 'Process' not in Result:
                    Result['Process'] = [line]
                else:
                    Result['Process'].append(line)

    return Result

def CompileToGenome(GenomeFileName,
                    DataDir):
    GeneFileName = DataDir + '/' + 'genes.tsv'
    assert os.path.exists(GeneFileName)

    GenomeFile = open(GenomeFileName, 'w')
    print(">E. coli", file=GenomeFile)

    Seq = ""
    fp = open(GeneFileName)

    # skip first line
    fp.readline()
    for Line in fp:
        Line = Line.strip()
        Fields = Line.split('\t')
        if len(Fields) == 3:
            continue

        TmpSeq = Fields[2][1:-1]
        Seq += TmpSeq

    fp.close()
    print(Seq, file=GenomeFile)
    GenomeFile.close()


def Compile(CodeFileNames,
            PrefixName,
            DataDir,
            TargetCodeModel,
            Verbose):

    OutputCodeName = PrefixName + ".py"

    CompilerData = FCompilerData()
    CompilerData.SetDataPath(os.path.realpath(DataDir))
    CodeInfo = Parse(CodeFileNames)

    GenomeFileName = PrefixName + ".ecoli.fa"
    if 'Organism' in CodeInfo:
        OrganismName = CodeInfo['Organism'].replace(' ', '_')
        GenomeFileName = PrefixName + "." + OrganismName + ".fa"

    CompileToGenome(GenomeFileName, DataDir)

    Dataset = LoadData(DataDir)
    SetUpCompilerData(Dataset, CompilerData)
    OutputFile = open(OutputCodeName, 'w')

    # TODO: use factory
    if TargetCodeModel == Target.TensorFlow:
        Writer = CodeGen.TFCodeWriter(OutputFile, 0)
    else:
        Writer = CodeGen.NumpyCodeWriter(OutputFile, 0)

    WriteLicense(Writer)
    WriteImport(Writer); Writer.BlankLine()
    WriteBody(Writer, CompilerData); Writer.BlankLine()
    WriteMain(Writer)

    OutputFile.close()


"""
"""
# PyCharm: set parameters configuration to "-L ../../data ecoli.lpp"
if __name__ == '__main__':

    version_str = 'lcc version {version}'.format(version=LCC_VERSION)

    parser = ArgumentParser(
        description='')
    parser.add_argument('-L',
                        dest='data_dir',
                        type=str,
                        help='Library/Data directory')
    parser.add_argument('-o', '--out-file',
                        dest='OutputFileName',
                        type=str,
                        default='cell',
                        help='Output code file')
    parser.add_argument('-march',
                        dest='arch',
                        type=str,
                        default='tensorflow',
                        choices=['tensorflow', 'tf', 'numpy', 'np'],
                        help='Specify output code type')
    parser.add_argument('-V', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument('infiles', 
                        type=str, nargs='+',
                        help='Code Files')
    parser.add_argument('-v', '--version',
                        action='version',
                        version=version_str)
    args = parser.parse_args()
    
    if not args.data_dir:
        parser.print_help()
        exit(1)

    for codefile in args.infiles:
        if not os.path.exists(codefile):
            print("Error: %s doesn't exist" % codefile)
            exit(1)

    print(args.arch)
    try:
        TargetCodeModel = Target.from_str(args.arch)
    except NotImplemented:
        parser.print_help()

    Compile(args.infiles,
            args.OutputFileName,
            args.data_dir,
            TargetCodeModel,
            args.verbose)
