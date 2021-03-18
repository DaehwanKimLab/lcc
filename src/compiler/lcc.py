#!/usr/bin/env python3
#
# Copyright 2021,
# Donghoon Lee <>,
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

            dataset[fname]=list_of_rows[1:]

    def dump_dataset():
        for key, value in dataset.items():
            print(key, len(value))

    for fname in os.listdir(data_dir):
        if fname.endswith('.tsv'):
            parse_tsv(data_dir, fname)

    dump_dataset()
    return dataset

    

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
    ]
    for line in lines:
        print(line, file=code_file)

        
def NewLine(code_file):
    print("\t", file=code_file)


def Compile(code_fname,
            data_dir,
            verbose):

    dataset = LoadData(data_dir)
    
    code_file = open(code_fname, 'w')
    
    WriteLicense(code_file)
    WriteImport(code_file); NewLine(code_file)
    WriteBody(code_file, dataset); NewLine(code_file)
    WriteMain(code_file)
    
    code_file.close()
    

"""
"""
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
