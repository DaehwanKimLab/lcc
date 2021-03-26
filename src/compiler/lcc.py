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


class FCompilerData:
    def __init__(self):
        self.MetaboliteNames = []
        self.MetaboliteName2Index = {}
        self.MetaboliteConcs = None

        self.TranscriptNTFreqs = None


class CodeWriter:
    def __init__(self, CodeFile, IndentLevel=1):
        self.IndentLevel = IndentLevel
        self.fp = CodeFile
        self.BuildIndentationPrefix()

    def BuildIndentationPrefix(self):
        self.IndentationPrefix = '\t' * self.IndentLevel

    def IncreaseIndent(self):
        self.IndentLevel += 1
        self.BuildIndentationPrefix()

    def DecreaseIndent(self):
        self.IndentLevel -= 1
        if self.IndentLevel < 0:
            self.IndentLevel = 0

        self.BuildIndentationPrefix()

    def GetIndentLevel(self):
        return self.IndentLevel

    def SetIndentLevel(self, IndentLevel):
        self.IndentLevel = IndentLevel
        self.BuildIndentationPrefix()

    def WriteStatement(self, Line):
        print("{Indent}{Line}".format(Indent=self.IndentationPrefix, Line=Line), file=self.fp)

    def WriteVariable(self, VariableName, Value):
        Line = '%s = %s' % (VariableName, Value)
        self.WriteStatement(Line)


def WriteLicense(code_file):
    for line in open("LICENSE.input"):
        line = line.strip()
        print(line, file=code_file)


def WriteImport(code_file):
    print("import os, sys", file=code_file)
    print("import numpy as np", file=code_file)
    print("import tensorflow as tf", file=code_file)
    print("from argparse import ArgumentParser, FileType", file=code_file)


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
    
    def RearrangeLstToDict(lst, pos_key):
        dict_out = {}
        for row in lst:
            row_without_key = []
            for column, i in zip(row, range(len(row))):
                if i == pos_key:
                    continue
                row_without_key.append(column)
            dict_item = {row[pos_key]: row_without_key}
        return dict_out

    Metabolites = Dataset['metaboliteConcentrations.tsv']
    CompilerData.MetaboliteConcs = np.zeros(len(Metabolites))
    for i, Value in enumerate(Metabolites):
        Name, Conc = Value
        assert Name not in CompilerData.MetaboliteName2Index
        CompilerData.MetaboliteName2Index[Name] = len(CompilerData.MetaboliteNames) # = i
        CompilerData.MetaboliteNames.append(Name) 
        CompilerData.MetaboliteConcs[i] = Conc



    def TranscriptionalElongation():

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

        # def GetMatrixRNAPPerTranscript():
        #     ActiveRNAPCount = 829
        #     RNAPPerRNA = np.zeros(len(RNAs))
        #     for position in range(ActiveRNAPCount):
        #         position = np.random.randint(1, len(RNAs))
        #         RNAPPerRNA[position] += 1
        #     if np.count_nonzero(RNAPPerRNA) == 0:
        #         print("WARNING: There is no RNAP on RNA.", file=sys.stderr)
        #     return RNAPPerRNA

        #     MatrixGetMatrixGenLocRNAP()
        #     return rnap_active
        #
        # def GetMatrixNtFlux():
        #     metabolites = dataset("metabolites.tsv")
        #     RearrangeLstToDict(metabolites, 1)
        #
        #     return acgu_flux
        #

        TranscriptNTFreqs = GetMatrixRNANTFreq()
        return TranscriptNTFreqs

        # mtrx_active_RNAP = GetMatrixActiveRNAP()
        #
        # mtrx_Nt_flux = GetMatrixNTFlux()

    CompilerData.TranscriptNTFreqs = TranscriptionalElongation()

"""
        # RNApol_avail =

        base_index_mw_A = dataset["metabolites.tsv"][0].index('ADENOSINE')
        base_index_mw_C = dataset["metabolites.tsv"][0].index('"CYTOSINE"')
        base_index_mw_G = dataset["metabolites.tsv"][0].index('"GUANINE"')
        base_index_mw_T = dataset["metabolites.tsv"][0].index('"THYMINE"')

        base_index_conc_A = dataset["metaboliteConcentrations.tsv"][0].index('"ADENOSINE"')
        base_index_conc_C = dataset["metaboliteConcentrations.tsv"][0].index('"CYTOSINE"')
        base_index_conc_G = dataset["metaboliteConcentrations.tsv"][0].index('"GUANINE"')
        base_index_conc_T = dataset["metaboliteConcentrations.tsv"][0].index('"THYMINE"')

        n_base_avail_A = dataset["metaboliteConcentrations.tsv"][base_index_mw_A, 1]\
                         * dataset["metabolite.tsv"][base_index_conc_A, 1]
        n_base_avail_C = dataset["metaboliteConcentrations.tsv"][base_index_mw_C, 1] \
                         * dataset["metabolite.tsv"][base_index_conc_C, 1]
        n_base_avail_G = dataset["metaboliteConcentrations.tsv"][base_index_mw_G, 1] \
                         * dataset["metabolite.tsv"][base_index_conc_G, 1]
        n_base_avail_T = dataset["metaboliteConcentrations.tsv"][base_index_mw_T, 1] \
                         * dataset["metabolite.tsv"][base_index_conc_T, 1]

        # turn dataset into dictionary or
        metabolites = dataset['metabolites.tsv']



        for metabolite in metabolites:
            for molecule in molecules_of_interest:
                if metabolite[0] == "ADENOSINE":
                    mol_weight_adenosine = metabolite[1]
            if metabolite[0] == "GUANINE":
                mol_weight_guanine = metabolite[1]
            if metabolite[0] == "CYTOSINE":
                mol_weight_cytosine = metabolite[1]
            if metabolite[0] == "URACIL":
                mol_weight_uracil = metabolite[1]

        # OR

        for metabolite in metabolites:
            if metabolite[0] == "ADENOSINE":
                mol_weight_adenosine = metabolite[1]
            if metabolite[0] == "GUANINE":
                mol_weight_guanine = metabolite[1]
            if metabolite[0] == "CYTOSINE":
                mol_weight_cytosine = metabolite[1]
            if metabolite[0] == "URACIL":
                mol_weight_uracil = metabolite[1]

        metaboliteConcentrations = dataset["metabololiteConcentrations.tsv"]

        for concentration in metaboliteConcentrations:
            if metaboliteConcentration[0] == "ADENOSINE":
                mol_weight_adenosine = metabolite[1]
            if metabolite[0] == "GUANINE":
                mol_weight_guanine = metabolite[1]
            if metabolite[0] == "CYTOSINE":
                mol_weight_cytosine = metabolite[1]
            if metabolite[0] == "URACIL":
                mol_weight_uracil = metabolite[1]

        dataset["metabolites.tsv"][1] # "ADENOSINE"	267.245
        dataset["metabolites.tsv"][1] # "GUANINE"	151.129
        dataset["metabolites.tsv"][1] # "CYTOSINE"    111.104
        dataset["metabolites.tsv"][1] # "THYMINE"	126.115

        dataset["metaboliteConcentrations.tsv"][1]    # "ADENOSINE"	1.30e-7
        dataset["metaboliteConcentrations.tsv"][1]    # "GUANINE"	1.90e-4

        # self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
        # self.bulkRnas = self.bulkMoleculesView(self.rnaIds)
        # self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
        # self.ppi = self.bulkMoleculeView('PPI[c]')
        # self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
        # self.flat_elongation = not sim._variable_elongation_transcription

"""


def WriteBody(CodeFile, CompilerData):
    Writer = CodeWriter(CodeFile, 0)
    Lines = [
        "def main(GenomeFileName, verbose):",
    ]
    for Line in Lines:
        Writer.WriteStatement(Line)

    Writer.IncreaseIndent()

    Writer.WriteVariable("CellCycles", 100)

    Writer.WriteVariable("AvogadroNum", 6.022141527E23)
    Writer.WriteStatement("CellVol = 7e-16 # Average E coli cell volume: 0.7 um3, which is 7e-16 liters.")

    np.save("MetaboliteConcs", CompilerData.MetaboliteConcs)
    Writer.WriteStatement("MetaboliteConcs = np.load('MetaboliteConcs.npy').astype('float32')")

    Writer.WriteStatement("MetaboliteConcsTF = tf.convert_to_tensor(MetaboliteConcs)")

    Writer.WriteStatement("print(MetaboliteConcs)")
    Writer.WriteStatement("print(MetaboliteConcsTF)")

    np.save("TranscriptNTFreqs", CompilerData.TranscriptNTFreqs)
    Writer.WriteStatement("TranscriptNTFreqs = np.load('TranscriptNTFreqs.npy').astype('float32')")
    Writer.WriteStatement("TranscriptNTFreqsTF = tf.convert_to_tensor(TranscriptNTFreqs)")
    Writer.WriteStatement("print(TranscriptNTFreqs)")
    Writer.WriteStatement("print(TranscriptNTFreqsTF)")

    NTIndexList = list()
    Writer.WriteStatement("NTCounts = np.zeros(4)")
    for i, NTName in enumerate(["ATP", "CTP", "GTP", "UTP"]):
        NTIndex = CompilerData.MetaboliteName2Index[NTName]
        NTIndexList.append(int(NTIndex))
        Writer.WriteStatement("NTCounts[%d] = MetaboliteConcs[%d] # %s" % (i, NTIndex, NTName))

    Writer.WriteStatement("NTConcsIndexTF = tf.reshape(tf.constant(" + str(NTIndexList) + "), [4, -1])")
    Writer.WriteStatement("print('NTConcsIndexTF = ', NTConcsIndexTF)")
    Writer.WriteStatement("print('NTCounts =', NTCounts)")

    Writer.WriteVariable('ElongationRate', 10)
    Writer.WriteVariable("ActiveRNAPCount", 829)

    Writer.WriteStatement("NumberOfUniqueTranscripts = len(TranscriptNTFreqs)")
    Writer.WriteStatement("RNAPPerTranscript = np.zeros(NumberOfUniqueTranscripts)")

    Writer.WriteStatement("for i in range(CellCycles):")

    Writer.IncreaseIndent()
    Writer.WriteStatement(  "for position in range(ActiveRNAPCount):")
    Writer.IncreaseIndent()
    Writer.WriteStatement(    "position = np.random.randint(1, NumberOfUniqueTranscripts)")
    Writer.WriteStatement(    "RNAPPerTranscript[position] += 1")
    Writer.DecreaseIndent()
    Writer.WriteStatement(  "if np.count_nonzero(RNAPPerTranscript) == 0:")
    Writer.IncreaseIndent()
    Writer.WriteStatement(    "print('WARNING: There is no RNAP on RNA.', file=sys.stderr)")
    Writer.DecreaseIndent()

    Writer.WriteStatement(  "print(RNAPPerTranscript)")

    Writer.WriteStatement(  "DeltaNTCounts = np.matmul(np.transpose(TranscriptNTFreqs), RNAPPerTranscript) * ElongationRate")
    Writer.WriteStatement(  "# DeltaNTCounts = np.sum(DeltaNTCounts, axis=0)")
    Writer.WriteStatement(  "print(\"DeltaNTCounts:\", DeltaNTCounts)")

    Writer.WriteStatement(  "DeltaNTCounts /= AvogadroNum")
    Writer.WriteStatement(  "DeltaNTCounts /= CellVol")

    Writer.WriteStatement(  "print('Available ACGU', tf.gather(MetaboliteConcsTF, NTConcsIndexTF))")
    Writer.WriteStatement(  "print(\"DeltaNTCounts (mol):\", DeltaNTCounts)")

    Writer.WriteStatement(  "MetaboliteConcsTF = tf.tensor_scatter_nd_sub(MetaboliteConcsTF, NTConcsIndexTF, DeltaNTCounts)")
    Writer.WriteStatement(  "print('After ACGU', tf.gather(MetaboliteConcsTF, NTConcsIndexTF))")

    Writer.WriteStatement(  "NTCounts -= DeltaNTCounts")



    Writer.WriteStatement(  "print(\"After one simulation unit,\")")
    Writer.WriteStatement(  "print(\"\tNTCounts =\", NTCounts)")
    Writer.DecreaseIndent()

    Writer.WriteStatement("i += 1")

    Writer.WriteStatement("if GenomeFileName != \"\":")
    Writer.IncreaseIndent()
    Writer.WriteStatement(  "GenomeFile = open(GenomeFileName, 'w')")
    Writer.WriteStatement(  "InputGenomeFile = open('cell.fa')")
    Writer.WriteStatement(  "for Line in InputGenomeFile:")
    Writer.IncreaseIndent()
    Writer.WriteStatement(    "Line = Line.strip()")
    Writer.WriteStatement(    "print(Line, file=GenomeFile)")
    Writer.DecreaseIndent()
    Writer.WriteStatement(  "GenomeFile.close()")



def WriteMain(code_file):
    lines = [
        "if __name__ == '__main__':",
        "\tparser = ArgumentParser(description='')",
        "\tparser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='also print some statistics to stderr')",
        "\tparser.add_argument('-g', '--genome', dest='GenomeFileName', type=str, default='', help='')",
        "",
        "\targs = parser.parse_args()",
        "\tmain(args.GenomeFileName, args.verbose)"
        "",
    ]
    for line in lines:
        print(line, file=code_file)


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
    for Line in open(GeneFileName):
        Line = Line.strip()
        Fields = Line.split('\t')
        if len(Fields) == 3:
            continue

        TmpSeq = Fields[2][1:-1]
        Seq += TmpSeq

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


    WriteLicense(OutputFile)
    WriteImport(OutputFile); NewLine(OutputFile)
    WriteBody(OutputFile, CompilerData); NewLine(OutputFile)
    WriteMain(OutputFile)

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
