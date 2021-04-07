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
from lccmodule import TCS
from lccmodule import TE
from lccmodule import CellMX

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

        # Temporary visualization variables
        writer.Variable_("SimStep", [])
        writer.Variable_("TE_ACGU", [])

        # Load CellMX
        CellMX.Write_CellMX_Init(writer)

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

        TE.Write_TE_Init(writer, CompilerData)
        writer.Statement("TE_Init()")

        TCS.Write_TCS_Init(writer, CompilerData)
        writer.Statement("TCS_Init()")

        # Run simulation
        writer.Statement("# Run simulation")
        writer.PrintStrg("Simulation begins...")
        with writer.Statement("for SimulationStep in range(SimulationSteps):"):
            writer.PrintStrg('=============================================')
            writer.PrintStVa('SimulationStep: ', "SimulationStep + 1")
            writer.Statement("SimStep.append(SimulationStep + 1)")
            writer.BlankLine()

            TE.Write_TE_Loop(writer)
            writer.Statement("TE_Loop()")

            TCS.Write_TCS_Loop(writer)
            writer.Statement("TCS_Loop()")

        # Temporary TE visualization code
        writer.Statement("fig, ax = plt.subplots()")
        # writer.Statement("ax.plot(SimStep, TE_ACGU)")
        writer.Statement("lines = ax.plot(SimStep, TE_ACGU)")
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

    Writer = codegen.TFCodeWriter(OutputFile, 0)
    # Writer = codegen.NumpyCodeWriter(OutputFile, 0)

    WriteLicense(Writer)
    WriteImport(Writer); Writer.BlankLine()
    WriteBody(Writer, CompilerData); Writer.BlankLine()
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
                        default='cell'
                                ''
                                '',
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
