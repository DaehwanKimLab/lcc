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


def WriteLicense(code_file):
    for line in open("LICENSE.input"):
        line = line.strip()
        print(line, file=code_file)


def WriteImport(code_file):
    print("import os, sys", file=code_file)
    print("import numpy as np", file=code_file)
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


def SetUpMatrix(Dataset):

    def RearrangeLstToDict(lst, pos_key):
        dict_out = {}
        for row in lst:
            row_without_key = []
            for column, i in zip(row, range(len(row))):
                if i == pos_key:
                    continue
                row_without_key.append(column)
            dict_item = {row[pos_key]: row_without_key}
        # dict_{}
        return dict_out

    def TranscriptionalElongation():

        def GetMatrixRNANTFreq():
            NTFreqTable = list()
            RNAs = Dataset['rnas.tsv']
            # num_rnas = len(dataset['rnas.tsv'])
            for RNA in RNAs:
                Seq = RNA[2]
                NTTotalCount = len(Seq)
                NTCounts = np.array([Seq.count("A"), Seq.count("C"), Seq.count("G"), Seq.count("U")])
                if np.sum(NTCounts) != NTTotalCount:
                    print("WARNING: RNA seq may contain non-ACGT text", file=sys.stderr)
                NTFreq = NTCounts / NTTotalCount
                NTFreqTable.append(NTFreq)
            return NTFreqTable

        # def GetMatrixGenLocRNAP():
        #
        # def GetMatrixActiveRNAP():
        #     MatrixGetMatrixGenLocRNAP()
        #     return rnap_active
        #
        # def GetMatrixNtFlux():
        #     metabolites = dataset("metabolites.tsv")
        #     RearrangeLstToDict(metabolites, 1)
        #
        #     return acgu_flux
        #
        mtrx_RNA_NT_freq = GetMatrixRNANTFreq()
        # mtrx_active_RNAP = GetMatrixActiveRNAP()
        #
        # mtrx_Nt_flux = GetMatrixNTFlux()

    TranscriptionalElongation()

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


def WriteBody(code_file, dataset):
    lines = [
        "def main(verbose):",
    ]
    for line in lines:
        print(line, file=code_file)

    for key, value in dataset.items():
        name = key.split(".")[0]
        print("\t%s = " % name, value, file=code_file)
        print("\tprint(%s)" % name, file=code_file)


def WriteMain(code_file):
    lines = [
        "if __name__ == '__main__':",
        "\tparser = ArgumentParser(description='')",
        "\tparser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='also print some statistics to stderr')",
        "",
        "\targs = parser.parse_args()",
        "\tmain(args.verbose)"
        "",
    ]
    for line in lines:
        print(line, file=code_file)


def NewLine(code_file):
    print("\t", file=code_file)


def Compile(code_fname,
            data_dir,
            verbose):

    dataset = LoadData(data_dir)
    SetUpMatrix(dataset)
    code_file = open(code_fname, 'w')

    WriteLicense(code_file)
    WriteImport(code_file); NewLine(code_file)
    WriteBody(code_file, dataset); NewLine(code_file)
    WriteMain(code_file)

    code_file.close()

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
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    if not args.data_dir:
        parser.print_help()
        exit(1)

    Compile("code.py",
            args.data_dir,
            args.verbose)
