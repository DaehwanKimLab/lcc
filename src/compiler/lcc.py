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


def WriteBody(code_file):
    lines = [
        "def main(verbose):",
        "\tprint(\"Hello World\")",
    ]
    for line in lines:
        print(line, file=code_file)

    
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
            verbose):
    print("TEST")
    
    code_file = open(code_fname, 'w')
    
    WriteLicense(code_file)
    WriteImport(code_file); NewLine(code_file)
    WriteBody(code_file); NewLine(code_file)
    WriteMain(code_file)
    

    code_file.close()
    

"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='')
    parser.add_argument('data_dir',
                        nargs='?',
                        type=str,
                        default="",
                        help='directory name for training dataset')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    Compile("code.py", 
            args.verbose)



