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
from codegen import Target
import inspect

LCC_VERSION = "0.1"

def lcc_dummy():
    pass
LCC_PATH = os.path.dirname(os.path.realpath(inspect.getsourcefile(lcc_dummy)))

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
        self.MolNames = []
        self.MolName2Index = {}
        self.MolCounts = None
        self.MolConcs = None

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
        MolCounts = Dataset['TwoComponentSystemsTemporary_DL.tsv'] # temporary data table
        CompilerData.MolCounts = np.zeros(len(MolCounts))
        for i, Value in enumerate(MolCounts):
            MolName, MolCount = Value
            CompilerData.MolName2Index[MolName] = len(CompilerData.MolNames) # = i
            CompilerData.MolNames.append(MolName)
            CompilerData.MolCounts[i] = MolCount

        # Temporary command to save a list of participating molecules in TCS
        np.save('TCSMolNames.npy', CompilerData.MolNames)

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

    writer.DebugPrintSwitch = False # ON/OFF switch to print data being processed
    writer.DebugAssertSwitch = True


    writer.WriteVariable('LCCDataPath', "\"" + CompilerData.GetDataPath() + "\"")
    Lines = [
        "def main(GenomeFileName, verbose):",
    ]
    for Line in Lines:
        writer.WriteStatement(Line)

    with writer:
        # Define simulation parameters
        writer.WriteStatement("# Define simulation parameters")
        writer.WriteVariable("CellCycles", 1)
        writer.WriteVariable("SimulationSteps", 100)
        writer.WriteVariable("TCSODETimeStep", 1)
        writer.WriteBlankLine()

        writer.WriteStatement("# This is a numpy code", TargetCode=Target.Numpy)
        writer.WriteStatement("# This is a tensorflow code", TargetCode=Target.TensorFlow)

        # Define key constants
        writer.WriteStatement("# Define key constants")
        writer.WriteVariable("AvogadroNum", 6.022141527E23)
        # Average E coli cell volume: 0.7 um3, which is 7e-16 liters
        writer.WriteVariable("CellVol", 7e-16)  # TO BE REPLACED AND MOVED INTO SIMULATION

        # Define accessory variables - TF version only
        writer.InitArrayWithOne('OneTF', [1])
        writer.InitArrayWithZero('ZeroTF', [1])
        writer.WriteBlankLine()

        # Save and load all Cell State matrices
        np.save("MetaboliteConcs", CompilerData.MetaboliteConcs)
        writer.WriteStatement("# Load metabolite concentration table")
        writer.WriteStatement("MetaboliteConcs = np.load(\"MetaboliteConcs.npy\").astype('float32')")
        writer.WriteStatement("MetaboliteConcsTF = tf.convert_to_tensor(MetaboliteConcs)")
        writer.WriteDebugPrintVar("MetaboliteConcs")
        writer.WriteDebugPrintVar("MetaboliteConcsTF")
        writer.WriteBlankLine()

        # Matrices for Transcript Elongation
        np.save("TranscriptNTFreqs", CompilerData.TranscriptNTFreqs)
        writer.WriteStatement("# Matrices for Transcript Elongation")

        # NT frequency table for transcripts
        writer.WriteStatement("# Fetch NT frequency table for transcripts")
        writer.WriteStatement("TranscriptNTFreqs = np.load(\"TranscriptNTFreqs.npy\").astype('float32')")
        writer.WriteStatement("TranscriptNTFreqsTF = tf.convert_to_tensor(TranscriptNTFreqs)")
        writer.WriteStatement("TranscriptNTFreqsTF = tf.transpose(TranscriptNTFreqs)")
        writer.WriteDebugPrintVar("TranscriptNTFreqs")
        writer.WriteDebugPrintVar("TranscriptNTFreqsTF")
        writer.WriteStatement("NumberOfUniqueTranscripts = len(TranscriptNTFreqs)")
        writer.WriteBlankLine()

        # NT counts per transcript
        NTIndexList = list()
        writer.WriteStatement("# Fetch NT concentration")
        ACGU = ["ATP", "CTP", "GTP", "UTP"]
        writer.WriteStatement("NTConcs = np.zeros(" + str(len(ACGU)) + ").astype('float32')")
        for i, NTName in enumerate(ACGU):
            NTIndex = CompilerData.MetaboliteName2Index[NTName]
            NTIndexList.append(int(NTIndex))
            writer.WriteStatement("NTConcs[%d] = MetaboliteConcs[%d] # %s" % (i, NTIndex, NTName))
        writer.WriteStatement("NTConcsTF = tf.convert_to_tensor(NTConcs)")
        writer.WriteStatement("NTConcsIndexTF = tf.reshape(tf.constant(" + str(NTIndexList) + "), [4, -1])")
        writer.WriteDebugPrintVar("NTConcs")
        writer.WriteDebugPrintVar("NTConcsTF")
        writer.WriteDebugPrintVar("NTConcsIndexTF")
        writer.WriteBlankLine()

        writer.WriteStatement("# Determine elongation rate")
        writer.WriteVariable("ElongationRate", 10) # TO BE REPLACED AND MOVED INTO SIMULATION
        writer.WriteBlankLine()

        writer.WriteStatement("# Determine active RNAP count")
        writer.WriteVariable("ActiveRNAPCount", 829) # TO BE REPLACED AND MOVED INTO SIMULATION
        writer.WriteBlankLine()

        # Set up matrices for Two Component Systems
        writer.WriteStatement("# Two Component Systems model")
        writer.WriteStatement("TCSModel = tf.keras.models.load_model(LCCDataPath + '/two_component.h5')")
        writer.WriteStatement("TCSODETimeStepTF = tf.reshape(tf.constant([TCSODETimeStep], dtype='float32'), (1, 1))")
        writer.WriteDebugPrintVar("TCSODDETimeStepTF")
        writer.WriteBlankLine()

        # Save and load all molecule count
        np.save("MolCounts.npy", CompilerData.MolCounts) # Maybe provided in another folder later
        writer.WriteStatement("# Load molecule count table")
        writer.WriteStatement("MolCounts = np.load(\"MolCounts.npy\").astype('float32')")
        writer.WriteStatement("MolCountsTF = tf.convert_to_tensor(MolCounts)")
        writer.WriteDebugPrintVar("MolCounts")
        writer.WriteDebugPrintVar("MolCountsTF")
        writer.WriteBlankLine()

        # Matrices for Two Component Systems
        TCSMolIndexList = []
        TCSMolNames = np.load('TCSMolNames.npy')
        writer.WriteStatement("TCSMolCounts = np.zeros(" + str(len(TCSMolNames)) + ").astype('float32')")
        for i, TCSMolName in enumerate(TCSMolNames):
            TCSMolIndex = CompilerData.MolName2Index[TCSMolName]
            TCSMolIndexList.append(int(TCSMolIndex))
            writer.WriteStatement("TCSMolCounts[%d] = MolCounts[%d] # %s" % (i, TCSMolIndex, TCSMolName))
        writer.WriteStatement("TCSMolCountsTF = tf.convert_to_tensor(TCSMolCounts)")
        writer.WriteStatement("TCSMolIndexTF = tf.reshape(tf.constant(" + str(TCSMolIndexList) + "), (-1, 1))")
        writer.WriteDebugPrintVar("TCSMolCountsTF")
        writer.WriteDebugPrintVar("TCSMolIndexTF")
        writer.WriteBlankLine()

        # Run simulation
        writer.WriteStatement("# Run simulation")
        writer.WritePrintStr("Simulation begins...")
        with writer.WriteStatement("for SimulationStep in range(SimulationSteps):"):
            writer.WritePrintStr('=============================================')
            writer.WritePrintStrVar('SimulationStep: ', "SimulationStep + 1")
            writer.WriteBlankLine()

            # Transcript Elongation (TE)
            writer.WriteStatement("# Transcript Elongation (TE)")
            writer.WriteBlankLine()

            # TE - Allocate RNAP to transcript
            writer.WriteStatement("# TE - Allocate RNAP to transcript")
            writer.WriteStatement("RNAPPerTranscriptTF = tf.zeros(NumberOfUniqueTranscripts)")
            with writer.WriteStatement("for RNAPPosition in range(ActiveRNAPCount):"):
                writer.WriteStatement("RNAPPosition = tf.random.uniform(shape=[1,1], minval=1, maxval=NumberOfUniqueTranscripts, dtype='int32')")
                writer.WriteStatement("RNAPPerTranscriptTF = tf.tensor_scatter_nd_add(RNAPPerTranscriptTF, RNAPPosition, OneTF)")
            writer.WriteDebugAssert("tf.math.reduce_sum(RNAPPerTranscriptTF) == ActiveRNAPCount", 'Active RNAP is not properly allocated')
            writer.WriteStatement("RNAPPerTranscriptTF = tf.reshape(RNAPPerTranscriptTF, [-1, 1])")
            writer.WriteDebugPrintVar("RNAPPerTranscriptTF")
            writer.WriteBlankLine()

            # TE - Determine NT consumption
            writer.WriteStatement("# TE - Determine NT consumption")
            writer.WriteStatement("DeltaNTCountsTF = tf.linalg.matmul(TranscriptNTFreqsTF, RNAPPerTranscriptTF) * ElongationRate")
            writer.WriteStatement("DeltaNTCountsTF = tf.reshape(DeltaNTCountsTF, -1)")
            writer.WriteStatement("DeltaNTConcsTF = DeltaNTCountsTF / (CellVol * AvogadroNum)") # final unit: mol/L
            writer.WriteDebugVariable("NTConcsAvailTF", "tf.gather(MetaboliteConcsTF, NTConcsIndexTF)")
            writer.WriteDebugPrintVar("DeltaNTCountsTF")
            writer.WriteDebugPrintVar("DeltaNTConcsTF")
            writer.WriteDebugPrintVar("NTConcsAvailTF")
            writer.WriteDebugStatement("tf.debugging.assert_positive(NTConcsAvailTF - DeltaNTConcsTF), 'The cell is running out of NTs'")
            writer.WriteBlankLine()

            # TE - Update NT concs
            writer.WriteStatement("# TE - Update NT counts")
            writer.WriteStatement("MetaboliteConcsTF = tf.tensor_scatter_nd_sub(MetaboliteConcsTF, NTConcsIndexTF, DeltaNTConcsTF)")
            writer.WriteDebugVariable("NTConcsNewTF", "tf.gather(MetaboliteConcsTF, NTConcsIndexTF)")
            writer.WriteDebugPrintVar("NTConcsNewTF")
            writer.WriteDebugStatement("tf.debugging.assert_none_equal(NTConcsNewTF, NTConcsAvailTF), 'NT consumption is not properly applied'")
            writer.WritePrintVar("NTConcsNewTF")
            writer.WriteBlankLine()

            # TE - Update Transcript counts - TO BE IMPLEMENTED

            # Two Component Systems code (TCS)
            writer.WriteStatement("# Two Component Systems code (TCS)")

            # TCS - Run machine learned model
            writer.WriteStatement("# TCS - Run machine learned model")
            writer.WriteStatement("TCSMolCountsTF = tf.gather(MolCountsTF, TCSMolIndexTF)")
            writer.WriteStatement("TCSMolConcsTF = TCSMolCountsTF / (CellVol * AvogadroNum)") # TCSMolConcsTF == y_init
            writer.WriteStatement("TCSModelInput = tf.concat([TCSMolConcsTF, TCSODETimeStepTF], axis=0)")
            writer.WriteStatement("TCSMolConcsNewTF = TCSModel.predict(tf.reshape(TCSModelInput, (1, -1)))[0, :]")
            writer.WriteBlankLine()

            # TCS - Replace values < 0 to 0
            writer.WriteStatement("# TCS - Replace values < 0 to 0")
            writer.WriteStatement("TCSMolConcsNewZeroIndexTF = tf.where(tf.less(TCSMolConcsNewTF, 0))")
            writer.WriteStatement("TCSMolConcsNewReplaceTF = tf.zeros(TCSMolConcsNewZeroIndexTF.shape[0])")
            writer.WriteStatement("TCSMolConcsNewTF = tf.tensor_scatter_nd_update(TCSMolConcsNewTF, TCSMolConcsNewZeroIndexTF, TCSMolConcsNewReplaceTF)")
            writer.WriteDebugPrintVar("TCSMolConcsNewZeroIndexTF")
            writer.WriteDebugPrintVar("TCSMolConcsNewReplaceTF")
            writer.WriteDebugPrintVar("TCSMolConcsNewTF")
            writer.WriteDebugStatement("tf.debugging.assert_non_negative(TCSMolConcsNewTF, 'TCSMolConcsNewTF contains a negative value(s)')")

            writer.WriteDebugPrintVar("TCSMolCountsTF")
            writer.WriteDebugPrintVar("TCSMolConcsTF")
            writer.WriteDebugPrintVar("TCSMolConcsNewTF")
            writer.WriteDebugPrintVar("TCSMolCountsNewTF")
            writer.WriteBlankLine()

            # TCS - Update TCS molecule counts
            writer.WriteStatement("# Update two component systems molecule counts")
            writer.WriteStatement("TCSMolCountsNewTF = TCSMolConcsNewTF * (CellVol * AvogadroNum)")
            writer.WriteStatement("MolCountsTF = tf.tensor_scatter_nd_update(MolCountsTF, TCSMolIndexTF, TCSMolCountsNewTF)")
            writer.WriteDebugStatement("TCSMolCountsUpdated = tf.gather(MolCountsTF, TCSMolIndexTF)")
            # writer.WriteDebugStatement("tf.assert_equal(TCSMolCountsUpdated, TCSMolCountsNewTF, 'TCSMolCounts is not updated')")
            writer.WritePrintVar("TCSMolCountsNewTF")
            writer.WriteBlankLine()

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
    SetUpMatrix(Dataset, CompilerData)
    OutputFile = open(OutputCodeName, 'w')

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

    Compile(args.infiles,
            args.OutputFileName,
            args.data_dir,
            args.verbose)
