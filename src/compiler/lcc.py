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
import codegen
from codegen import Target
import inspect
from lccmodule import CellMX
from lccmodule import TCS
from lccmodule import TE
from lccmodule import RNADeg

LCC_VERSION = "0.1"

def lcc_dummy():
    pass
LCC_PATH = os.path.dirname(os.path.realpath(inspect.getsourcefile(lcc_dummy)))

class FCompilerData:
    def __init__(self):
        self.MetaboliteNames = []
        self.MetaboliteName2Index = {}
        self.MetaboliteConcs = None

        self.TranscriptNTFreqs = []

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

        self.RNAName2Index = {}
        self.RNANames = []
        self.RNAMonomerID2Index = {}
        self.RNAMonomerIDs = []
        self.RNAID2Index = {}
        self.RNAIDs = []
        self.RNAHalfLife = []
        self.RNASeq = []
        self.RNAType = []
        self.RNAMolWeight= []

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

    dump_dataset()
    return dataset


def SetUpCompilerData(Dataset, CompilerData):

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
            CompilerData.RNAHalfLife.append(HalfLife)
            CompilerData.RNASeq.append(Seq)
            CompilerData.RNAType.append(Type)
            CompilerData.RNAMolWeight.append(MolWeight)

        np.save('RNAIDs.npy', CompilerData.RNAIDs)
        return

    def GetMXRNANTFreq():
        NTFreqTable = list()
        # num_rnas = len(dataset['rnas.tsv'])
        RNAs = Dataset['rnas.tsv'] # Columns: [0]halfLife, [1]name, [2]seq, [3]type, [4]modifiedForms, [5]monomerId
        for RNA in RNAs:
            # HalfLife = RNA[0]
            # Name = RNA[1]
            Seq = RNA[2]
            # Type = RNA[3]
            # ModifiedForms = RNA[4]
            # MonomerID = RNA[5]
            NTTotalCount = len(Seq)
            NTCounts = np.array([Seq.count("A"), Seq.count("C"), Seq.count("G"), Seq.count("U")])
            if np.sum(NTCounts) != NTTotalCount:
                print("WARNING: RNA seq may contain non-ACGT text.", file=sys.stderr)
            NTFreq = NTCounts / NTTotalCount
            CompilerData.TranscriptNTFreqs.append(NTFreq)
        return

    SetUpDataIndexes()
    GetMXRNANTFreq()
    # CompilerData.TCSMolNames = TwoComponentSystems()


def WriteBody(writer, CompilerData):

    writer.DebugPrintSwitch = False # ON/OFF switch to print data being processed
    writer.DebugAssertSwitch = True


    writer.Variable_('LCCDataPath', "\"" + CompilerData.GetDataPath() + "\"")
    Lines = [
        "def main(GenomeFileName, verbose):",
    ]
    for Line in Lines:
        writer.Statement(Line)

    with writer:
        # Define simulation parameters
        writer.Statement("# Define simulation parameters")
        writer.Variable_("CellCycles", 1)
        writer.Variable_("SimulationSteps", 100)

        writer.BlankLine()

        writer.Statement("# This is a numpy code", TargetCode=Target.Numpy)
        writer.Statement("# This is a tensorflow code", TargetCode=Target.TensorFlow)

        # Define key constants
        writer.Statement("# Define key constants")
        writer.Variable_("AvogadroNum", 6.022141527E23)
        # Average E coli cell volume: 0.7 um3, which is 7e-16 liters
        writer.Variable_("CellVol", 7e-16)  # TO BE REPLACED AND MOVED INTO SIMULATION

        # Define accessory variables - TF version only
        writer.InitArrayWithOne('OneTF', [1])
        writer.InitArrayWithZero('ZeroTF', [1])
        writer.BlankLine()

        # GLOBAL INIT?

        # Load CellMX
        CellMX.Write_CellMX_Init(writer)
        writer.Statement("CellMX = FCellMX()")
        writer.BlankLine()

        # Save and load all Cell State matrices
        np.save("MetaboliteConcs", CompilerData.MetaboliteConcs)
        writer.Statement("# Load metabolite concentration table")
        writer.Statement("CellMX.MetaboliteConcs = np.load(\"MetaboliteConcs.npy\").astype('float32')")
        writer.Statement("CellMX.MetaboliteConcsTF = tf.convert_to_tensor(CellMX.MetaboliteConcs)")
        writer.DebugPVar("CellMX.MetaboliteConcs")
        writer.DebugPVar("CellMX.MetaboliteConcsTF")
        writer.BlankLine()

        # Save and load all molecule count
        np.save("MolCounts.npy", CompilerData.MolCounts) # Maybe provided in another folder later
        writer.Statement("# Load molecule count table")
        writer.Statement("MolCounts = np.load(\"MolCounts.npy\").astype('float32')")
        writer.Statement("CellMX.MolCountsTF = tf.convert_to_tensor(MolCounts)")
        writer.DebugPVar("MolCounts")
        writer.DebugPVar("CellMX.MolCountsTF")
        writer.BlankLine()

        # Save and load all RNA count (placeholder)
        writer.Statement("# RNADeg - Transcript Counts")
        writer.Variable_("DefaultCount", 100)
        writer.Statement("RNACounts = np.ones(" + str(len(CompilerData.RNAIDs)) + ").astype('int32') * DefaultCount")
        writer.Statement("CellMX.RNACountsTF = tf.convert_to_tensor(RNACounts)")
        writer.DebugPVar("CellMX.RNACountsTF")
        writer.BlankLine()


        # Load initialization functions for each process
        TE.Write_TE_Init(writer, CompilerData)
        TCS.Write_TCS_Init(writer, CompilerData)
        RNADeg.Write_RNADeg_Init(writer, CompilerData)

        # Load loop functions for each process
        TE.Write_TE_Loop(writer)
        TCS.Write_TCS_Loop(writer)
        RNADeg.Write_RNADeg_Loop(writer)

        # Run initialization functions for each process
        writer.Statement("TE_Init()")
        writer.Statement("TCS_Init()")
        writer.Statement("RNADeg_Init()")

        # Run simulation
        writer.Statement("# Run simulation")
        writer.PrintStrg("Simulation begins...")
        with writer.Statement("for SimulationStep in range(SimulationSteps):"):
            writer.PrintStrg('=============================================')
            writer.PrintStVa('SimulationStep: ', "SimulationStep + 1")
            writer.Statement("CellMX.SimStep.append(SimulationStep + 1)")
            writer.BlankLine()

            # Run loop functions for each process
            writer.Statement("TE_Loop()")
            writer.Statement("TCS_Loop()")
            writer.Statement("RNADeg_Loop()")
            writer.BlankLine()

        # Temporary TE visualization code
        writer.Statement("fig, ax = plt.subplots()")
        # writer.Statement("ax.plot(CellMX.SimStep, CellMX.TE_ACGU)")
        writer.Statement("lines = ax.plot(CellMX.SimStep, CellMX.TE_ACGU)")
        writer.Statement("labels = ['ATP', 'CTP', 'GTP', 'UTP']")
        writer.Statement("ax.legend(lines, labels)")
        writer.Statement("ax.set(xlabel='SimStep', ylabel='Concentration (M)', title='NTP level')")
        writer.Statement("ax.grid()")
        writer.Statement("plt.show()")

        writer.BlankLine()
        # End of simulation

        # Print input genome
        with writer.Statement("if GenomeFileName != \"\":"):
            writer.Statement("GenomeFile = open(GenomeFileName, 'w')")
            writer.Statement("InputGenomeFile = open('cell.fa')")
            with writer.Statement("for Line in InputGenomeFile:"):
                writer.Statement("Line = Line.strip()")
                writer.Statement("print(Line, file=GenomeFile)")
            writer.Statement("GenomeFile.close()")

def WriteMain(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)

    Writer.BlankLine()
    with Writer.Statement("if __name__ == '__main__':"):
        Writer.Statement("parser = ArgumentParser(description='')")
        Writer.Statement("parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='also print some statistics to stderr')")
        Writer.Statement("parser.add_argument('-g', '--genome', dest='GenomeFileName', type=str, default='', help='')")
        Writer.BlankLine()
        Writer.Statement("args = parser.parse_args()")
        Writer.Statement("main(args.GenomeFileName, args.verbose)")

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
        Writer = codegen.TFCodeWriter(OutputFile, 0)
    else:
        Writer = codegen.NumpyCodeWriter(OutputFile, 0)

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
