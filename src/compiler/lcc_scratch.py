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


def TranscriptionalElongation():
    # base_types = ["A", "C", "G", "T"]
    atgc_table = list()

    gene_len = dataset("genes.tsv")[0]  # len
    genes_seq = dataset("genes.tsv")[2]  # seq

    for seq, length in enumerate(genes_seq, gene_len):
        atgc_count = [seq.count("A")], [seq.count("T")], seq.count("G"), seq.count("C")
        atgc_ratio = atgc_count / length
        atgc_table.append(atgc_ratio)

    RNApol_avail =

    base_index_mw_A = dataset(metabolites.tsv)[0].index('"ADENOSINE"')
    base_index_mw_C = dataset(metabolites.tsv)[0].index('"CYTOSINE"')
    base_index_mw_G = dataset(metabolites.tsv)[0].index('"GUANINE"')
    base_index_mw_T = dataset(metabolites.tsv)[0].index('"THYMINE"')

    base_index_conc_A = dataset(metaboliteConcentrations.tsv)[0].index('"ADENOSINE"')
    base_index_conc_C = dataset(metaboliteConcentrations.tsv)[0].index('"CYTOSINE"')
    base_index_conc_G = dataset(metaboliteConcentrations.tsv)[0].index('"GUANINE"')
    base_index_conc_T = dataset(metaboliteConcentrations.tsv)[0].index('"THYMINE"')

    n_base_avail_A = dataset(metaboliteConcentrations.tsv)[base_index_mw_A, 1] \
                     * dataset(metabolite.tsv)[base_index_conc_A, 1]
    n_base_avail_C = dataset(metaboliteConcentrations.tsv)[base_index_mw_C, 1] \
                     * dataset(metabolite.tsv)[base_index_conc_C, 1]
    n_base_avail_G = dataset(metaboliteConcentrations.tsv)[base_index_mw_G, 1] \
                     * dataset(metabolite.tsv)[base_index_conc_G, 1]
    n_base_avail_T = dataset(metaboliteConcentrations.tsv)[base_index_mw_T, 1] \
                     * dataset(metabolite.tsv)[base_index_conc_T, 1]

    dataset(metabolites.tsv[1])  # "ADENOSINE"	267.245
    dataset(metabolites.tsv[1])  # "GUANINE"	151.129
    dataset(metabolites.tsv[1])  # "CYTOSINE"    111.104
    dataset(metabolites.tsv[1])  # "THYMINE"	126.115

    dataset(metaboliteConcentrations.tsv[1])  # "ADENOSINE"	1.30e-7
    dataset(metaboliteConcentrations.tsv[1])  # "GUANINE"	1.90e-4

TranscriptionalElongation()