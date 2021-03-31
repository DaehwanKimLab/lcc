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
# import matplotlib.pyplot as plt
import datetime
from argparse import ArgumentParser, FileType
import codegen


class FCompilerData:
    def __init__(self):
        self.MetaboliteNames = []
        self.MetaboliteName2Index = {}
        self.MetaboliteConcs = None

        self.TranscriptNTFreqs = None

        self.RXNIDs = []
        self.RXNID2Index = {}
        self.RXNStoichiometry = []
        self.RXNReversibility = []
        self.RXNEnzyme = []

        # self.TCSRXNs = []
        self.TCSMolNames = []
        self.TCSMolName2Index = {}
        self.TCSMolConcs = None


def WriteLicense(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)
    for line in open("LICENSE.input"):
        line = line.strip()
        Writer.WriteStatement(line)
    Writer.SetIndentLevel(tmpLevel)


def WriteImport(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)

    Writer.WriteStatement("import os, sys")
    Writer.WriteStatement("import numpy as np")
    Writer.WriteStatement("import tensorflow as tf")
    Writer.WriteStatement("from argparse import ArgumentParser, FileType")

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

    dump_dataset()
    return dataset


def SetUpMatrix(Dataset, CompilerData):

    def SetUpDataIndexes():
        Metabolites = Dataset['metaboliteConcentrations.tsv']
        CompilerData.MetaboliteConcs = np.zeros(len(Metabolites))
        for i, Value in enumerate(Metabolites):
            Name, Conc = Value
            assert Name not in CompilerData.MetaboliteName2Index
            CompilerData.MetaboliteName2Index[Name] = len(CompilerData.MetaboliteNames) # = i
            CompilerData.MetaboliteNames.append(Name)
            CompilerData.MetaboliteConcs[i] = Conc

        RXNs = Dataset['reactions.tsv']
        for i, Value in enumerate(RXNs):
            RXNID, RXNStoichiometry, RXNReversibility, RXNEnzyme = Value
            assert RXNID not in CompilerData.RXNID2Index
            CompilerData.RXNID2Index[RXNID] = len(CompilerData.RXNIDs) # = i
            CompilerData.RXNIDs.append(RXNID)
            CompilerData.RXNStoichiometry.append(RXNStoichiometry)
            CompilerData.RXNReversibility.append(RXNReversibility)
            CompilerData.RXNEnzyme.append(RXNEnzyme)

        # TwoComponentSystems = Dataset['twoComponentSystems.tsv']
        TwoComponentSystemsTEMP = Dataset['TwoComponentSystemsTemporary_DL.tsv'] # temporary datatable
        CompilerData.TCSMolCounts = np.zeros(len(TwoComponentSystemsTEMP))
        for i, Value in enumerate(TwoComponentSystemsTEMP):
            TCSMolName, TCSMolCount = Value
            CompilerData.TCSMolName2Index[TCSMolName] = len(CompilerData.TCSMolNames) # = i
            CompilerData.TCSMolNames.append(TCSMolName)
            CompilerData.TCSMolCounts[i] = TCSMolCount


    def TranscriptElongation():

        def GetMatrixRNANTFreq():
            NTFreqTable = list()
            # num_rnas = len(dataset['rnas.tsv'])
            RNAs = Dataset['rnas.tsv']
            for RNA in RNAs:
                Seq = RNA[2]
                NTTotalCount = len(Seq)
                NTCounts = np.array([Seq.count("A"), Seq.count("C"), Seq.count("G"), Seq.count("U")])
                if np.sum(NTCounts) != NTTotalCount:
                    print("WARNING: RNA seq may contain non-ACGT text.", file=sys.stderr)
                NTFreq = NTCounts / NTTotalCount
                NTFreqTable.append(NTFreq)
            return NTFreqTable




        # def GetMatrixNtFlux():
        #     metabolites = dataset("metabolites.tsv")
        #     RearrangeLstToDict(metabolites, 1)
        #
        #     return acgu_flux
        #

        TranscriptNTFreqs = GetMatrixRNANTFreq()
        return TranscriptNTFreqs

    # def TwoComponentSystems():
    #
    #
    #
    #     # assert TCSMolNames in
    #     return TCSMolNames
        # mtrx_active_RNAP = GetMatrixActiveRNAP()
        #
        # mtrx_Nt_flux = GetMatrixNTFlux()

    SetUpDataIndexes()
    CompilerData.TranscriptNTFreqs = TranscriptElongation()
    # CompilerData.TCSMolNames = TwoComponentSystems()


def WriteBody(writer, CompilerData):

    # ON/OFF switch to print data being processed
    ShowData = False

    Lines = [
        "def main(GenomeFileName, verbose):",
    ]
    for Line in Lines:
        writer.WriteStatement(Line)

    with writer:
        # Define simulation parameters and global constants
        writer.WriteStatement("# Define simulation parameters")
        writer.WriteVariable("CellCycles", 1)
        writer.WriteVariable("SimulationSteps", 100)
        writer.WriteBlankLine()

        # Define key constants
        writer.WriteStatement("# Define key constants")
        writer.WriteVariable("AvogadroNum", 6.022141527E23)
        # Average E coli cell volume: 0.7 um3, which is 7e-16 liters
        writer.WriteVariable("CellVol", 7e-16)  # TO BE REPLACED AND MOVED INTO SIMULATION

        # Define accessory variables - TF version only
        writer.InitArrayWithOne('OneTF', [1])
        writer.WriteBlankLine()

        # Set up static tables for transcript elongation
        np.save("MetaboliteConcs", CompilerData.MetaboliteConcs)
        writer.WriteStatement("# Load metabolite concentration table")
        writer.WriteStatement("MetaboliteConcs = np.load(\"MetaboliteConcs.npy\").astype('float32')")
        writer.WriteStatement("MetaboliteConcsTF = tf.convert_to_tensor(MetaboliteConcs)")
        if ShowData:
            writer.WriteStatement("print(MetaboliteConcs)")
            writer.WriteStatement("print(MetaboliteConcsTF)")
        writer.WriteBlankLine()

        np.save("TranscriptNTFreqs", CompilerData.TranscriptNTFreqs)
        writer.WriteStatement("# Load NT frequency table for transcripts")
        writer.WriteStatement("TranscriptNTFreqs = np.load(\"TranscriptNTFreqs.npy\").astype('float32')")
        writer.WriteStatement("TranscriptNTFreqsTF = tf.convert_to_tensor(TranscriptNTFreqs)")
        writer.WriteStatement("TranscriptNTFreqsTF = tf.transpose(TranscriptNTFreqs)")
        if ShowData:
            writer.WriteStatement("print(TranscriptNTFreqs)")
            writer.WriteStatement("print(TranscriptNTFreqsTF)")
        writer.WriteStatement("NumberOfUniqueTranscripts = len(TranscriptNTFreqs)")
        writer.WriteBlankLine()

        NTIndexList = list()
        writer.WriteStatement("# Fetch NT counts")
        # writer.InitArrayWithZero('NTCounts', [4], 'float32')
        writer.WriteStatement("NTCounts = np.zeros(4).astype('float32')")
        for i, NTName in enumerate(["ATP", "CTP", "GTP", "UTP"]):
            NTIndex = CompilerData.MetaboliteName2Index[NTName]
            NTIndexList.append(int(NTIndex))
            writer.WriteStatement("NTCounts[%d] = MetaboliteConcs[%d] # %s" % (i, NTIndex, NTName))
        writer.WriteStatement("NTCounts = tf.convert_to_tensor(NTCounts)")
        writer.WriteStatement("NTConcsIndexTF = tf.reshape(tf.constant(" + str(NTIndexList) + "), [4, -1])")
        if ShowData:
            writer.WriteStatement("print('NTConcsIndexTF = ', NTConcsIndexTF)")
            writer.WriteStatement("print(\"NTCounts =\", NTCounts)")
        writer.WriteBlankLine()

        writer.WriteStatement("# Determine elongation rate")
        writer.WriteVariable("ElongationRate", 10) # TO BE REPLACED AND MOVED INTO SIMULATION
        writer.WriteBlankLine()

        writer.WriteStatement("# Determine active RNAP count")
        writer.WriteVariable("ActiveRNAPCount", 829) # TO BE REPLACED AND MOVED INTO SIMULATION
        writer.WriteBlankLine()
        writer.WriteBlankLine()


        # Run simulation
        writer.WriteStatement("# Run simulation")
        with writer.WriteStatement("for SimulationStep in range(SimulationSteps):"):
            writer.WriteStatement("print('SimulationStep: ', SimulationStep + 1)")

            # Transcript elongation code
            writer.WriteStatement("# Transcript elongation code")
            # writer.WriteStatement("RNAPPerTranscriptTF = tf.zeros(NumberOfUniqueTranscripts)")
            writer.InitArrayWithZero("RNAPPerTranscriptTF", "[NumberOfUniqueTranscripts]")

            writer.WriteBlankLine()

            # Allocate RNAP to transcript
            writer.WriteStatement("# Allocate RNAP to transcript")
            with writer.WriteStatement("for position in range(ActiveRNAPCount):"):
                writer.WriteStatement("position = tf.random.uniform(shape=[1,1], minval=1, maxval=NumberOfUniqueTranscripts, dtype='int32')")
                writer.WriteStatement("RNAPPerTranscriptTF = tf.tensor_scatter_nd_add(RNAPPerTranscriptTF, position, OneTF)")
                writer.WriteBlankLine()
            with writer.WriteStatement("if tf.math.count_nonzero(RNAPPerTranscriptTF) == 0:"):
                writer.WriteStatement("print('WARNING: There is no RNAP on RNA.', file=sys.stderr)")
            # writer.WriteStatement("RNAPPerTranscriptTF = tf.reshape(RNAPPerTranscriptTF, [-1, 1])")
            writer.Reshape("RNAPPerTranscriptTF", "RNAPPerTranscriptTF", [-1, 1])
            if ShowData:
                writer.WriteStatement("print(RNAPPerTranscriptTF)")
            writer.WriteBlankLine()

            # Determine NT consumption
            writer.WriteStatement("# Determine NT consumption")
            writer.WriteStatement("DeltaNTCounts = tf.linalg.matmul(TranscriptNTFreqsTF, RNAPPerTranscriptTF) * ElongationRate")
            # writer.WriteStatement("DeltaNTCounts = tf.reshape(DeltaNTCounts, -1)")
            writer.Reshape("DeltaNTCounts", "DeltaNTCounts", -1)
            if ShowData:
                writer.WriteStatement("print(\"DeltaNTCounts:\", DeltaNTCounts)")
            writer.WriteStatement("DeltaNTCounts /= AvogadroNum")
            writer.WriteStatement("DeltaNTCounts /= CellVol")
            writer.WriteStatement("print('Available ACGU (mol)', tf.gather(MetaboliteConcsTF, NTConcsIndexTF))")
            writer.WriteStatement("print(\"DeltaNTCounts (mol):\", DeltaNTCounts)")
            writer.WriteBlankLine()

            # Update NT counts
            writer.WriteStatement("# Update NT counts")
            writer.WriteStatement("MetaboliteConcsTF = tf.tensor_scatter_nd_sub(MetaboliteConcsTF, NTConcsIndexTF, DeltaNTCounts)")
            if ShowData:
                writer.WriteStatement("print('After ACGU', tf.gather(MetaboliteConcsTF, NTConcsIndexTF))")
            writer.WriteStatement("NTCounts -= DeltaNTCounts")
            writer.WriteStatement("print(\"After one simulation unit,\")")
            writer.WriteStatement("print(\"\tNTCounts =\", NTCounts)")
            writer.WriteBlankLine()

            # Update Transcript counts - TO BE IMPLEMENTED

        writer.WriteBlankLine()
        # End of simulation

        # Print input genome
        with writer.WriteStatement("if GenomeFileName != \"\":"):
            writer.WriteStatement("GenomeFile = open(GenomeFileName, 'w')")
            writer.WriteStatement("InputGenomeFile = open('cell.fa')")
            with writer.WriteStatement("for Line in InputGenomeFile:"):
                writer.WriteStatement("Line = Line.strip()")
                writer.WriteStatement("print(Line, file=GenomeFile)")
            writer.WriteStatement("GenomeFile.close()")


def WriteMain(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)

    Writer.WriteBlankLine()
    with Writer.WriteStatement("if __name__ == '__main__':"):
        Writer.WriteStatement("parser = ArgumentParser(description='')")
        Writer.WriteStatement("parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='also print some statistics to stderr')")
        Writer.WriteStatement("parser.add_argument('-g', '--genome', dest='GenomeFileName', type=str, default='', help='')")
        Writer.WriteBlankLine()
        Writer.WriteStatement("args = parser.parse_args()")
        Writer.WriteStatement("main(args.GenomeFileName, args.verbose)")

    Writer.SetIndentLevel(tmpLevel)


def NewLine(code_file):
    print("\t", file=code_file)


def Parse(CodeFileName):
    Result = {}

    CodeFile = open(CodeFileName)
    for line in CodeFile:
        line = line.strip()
        if line.startswith('#'):
            continue

        if line.startswith('TemplateOrganism'):
            Organism = line.split(':')[1]
            Result['Organism'] = Organism
        elif line.startswith('RNATranscription'):
            if 'Process' not in Result:
                Result['Process'] = [line]
            else:
                Result['Process'].append(line)


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


def Compile(CodeFileName,
            OutputFileName,
            DataDir,
            Verbose):

    PrefixName = 'cell'

    CompilerData = FCompilerData()
    CodeInfo = Parse(CodeFileName)

    CompileToGenome(PrefixName + '.fa', DataDir)

    Dataset = LoadData(DataDir)
    SetUpMatrix(Dataset, CompilerData)
    OutputFile = open(OutputFileName, 'w')

    Writer = codegen.TFCodeWriter(OutputFile, 0)
    # Writer = codegen.NumpyCodeWriter(OutputFile, 0)

    WriteLicense(Writer)
    WriteImport(Writer); Writer.WriteBlankLine()
    WriteBody(Writer, CompilerData); Writer.WriteBlankLine()
    WriteMain(Writer)

    OutputFile.close()

    
"""
"""
# PyCharm: set parameters configuration to "-d ../../data"
if __name__ == '__main__':
    parser = ArgumentParser(
        description='')
    parser.add_argument('-d', '--data',
                        dest='data_dir',
                        type=str,
                        help='Data directory')
    parser.add_argument('-c', '--code',
                        dest='CodeFileName',
                        type=str,
                        help='life source code filename')
    parser.add_argument('-o', '--out-file',
                        dest='OutputFileName',
                        type=str,
                        default='cell.py',
                        help='Output code file')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    if not args.data_dir:
        parser.print_help()
        exit(1)

    if not os.path.exists(args.CodeFileName):
        print("Error: %s doesn't exist" % args.CodeFileName)
        exit(1)

    Compile(args.CodeFileName,
            args.OutputFileName,
            args.data_dir,
            args.verbose)
