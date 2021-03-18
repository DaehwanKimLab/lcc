from __future__ import absolute_import, division, print_function

import os, sys
import numpy as np
import csv
import tensorflow as tf

# DL - Desired code to be written by the compiler

print("TensorFlow Version: %s" % tf.__version__)


# The structure of data files requires the following format:
#   First Row for Headings  ex) Name    Property_1  Property_2 Property_N
#   Rest of the Rows for Data   Gene1   123         456         789

data_files = [
    "reconstruction/ecoli/flat/metaboliteConcentrations.tsv", # concentrations of known metabolites
    "reconstruction/ecoli/flat/metabolites.tsv" # molecular weight dictionary for all metabolites
    "CP_log.csv"
]

# CellState matrix contains a specific type of values for all molecules
class CellState:
    def __init__(self):
        self.data_files = []
        self.n_data_files = 0
        self.data_types_list = []
        self.n_data_types = 0
        self.n_items = 0
        self.cellStateMatx = tf.zeros(0)
        self.cellStateLabel = tf.zeros(0)

    def InitializeCellStateMatrix(self, data_files):
        self.data_files = data_files
        self.n_data_files = len(self.data_files)

        for i in range(self.n_data_files):
            with open(self.data_files[i]) as data_file:
                csv_reader = csv.reader(data_file, delimiter="\t")
                list_of_rows = list(csv_reader)

            data_table = np.transpose(np.array(list_of_rows))
            n_rows = len(_data_table)

            for row in range(_n_rows):
                if row == 0:
                    row_tag = data_table[row, 0]
                else:
                    entry_tag = str(data_table[0]) + " of " + str(_mol_names[i])
                    self.cellStateLabel = tf.concat([self.cellStateLabel, [entry_tag[i]]], 0)
                    self.cellStateMatx = tf.concat([self.cellStateMatx, [data_table[i]]], 0)


Ecoli = CellState()
Ecoli.InitializeCellStateMatrix(data_files)


#                 # np.append(self.cellStateMatx, _value_matx)
#                 _value_matx = tf.constant(_value_matx, dtype=float)
#                 self.cellStateMatx = tf.concat([self.cellStateMatx, _value_matx], 0)
#                 # self.cellStateMatx.append(_value_matx)
#                 self.cellStateLabel.append(_label_matx)
#
#                 print(i)
#                 print(_value_matx)
#                 print(_label_matx)
#                 print(self.cellStateMatx)
#                 current_pos = len(self.cellStateMatx)
#         print(self.cellStateMatx)
#
#
# Ecoli = CellState()
# Ecoli.InitializeCellStateMatrix(data_files)
