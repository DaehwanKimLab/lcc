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
from os import listdir
import numpy as np
import csv
# import tensorflow as tf
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import CodeGen
from CodeGen import Target
import inspect
from CompilerData import FCompilerData
from lccvariable import Simulation
from lccvariable import Constant
from lccvariable import Environment
from lccvariable import CellState
from lccvariable import CellProcess
from lccvariable import Visualization_2D
# from lccvariable.cellstate import GenomeState
# from lccvariable.cellstate import RNAState
# from lccvariable.cellstate import ProteinMonomerState
# from lccvariable.cellstate import ComplexState
# from lccvariable.cellstate import LipidState
# from lccvariable.cellstate.genomestate import GeneState
# from lccvariable.cellstate.genomestate import PromoterState
from lccmodule.cellularprocesses import *
from lccmodule.simulation import ReactionExecution
from lccmodule.simulation import ConservationOfMass
from lccmodule.cellularprocesses.synthesis import Replication
from lccmodule.cellularprocesses.synthesis import Transcription
from lccmodule.cellularprocesses.synthesis import Translation
from lccmodule.cellularprocesses.degradation import DNADegradation
from lccmodule.cellularprocesses.degradation import RNADegradation
from lccmodule.cellularprocesses.degradation import ProteinDegradation
from lccmodule.cellularprocesses.modification import DNAModifications
from lccmodule.cellularprocesses.modification import RNAModifications
from lccmodule.cellularprocesses.modification import PostTranslationalModifications
from lccmodule.cellularprocesses.conversion import Complexation
from lccmodule.cellularprocesses.conversion import Equilibrium
from lccmodule.cellularprocesses.metabolism import MetabolicNetwork
from lccmodule.cellularprocesses.signaling import OneComponentSystems
from lccmodule.cellularprocesses.signaling import TwoComponentSystems
from lccmodule.cellularprocesses.signaling import ECFSigmaFactor
from lccmodule.cellularprocesses.biophysics import Cytokinesis
from lccmodule.cellularprocesses.biophysics import BiophysicalProperties
from lccmodule.cellularprocesses.biophysics import Flagellum



LCC_VERSION = "0.1"

def lcc_dummy():
    pass
LCC_PATH = os.path.dirname(os.path.realpath(inspect.getsourcefile(lcc_dummy)))

def WriteLicense(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)
    for line in open(LCC_PATH + "/LICENSE.input"):
        line = line.strip()
        Writer.Statement(line)
    Writer.SetIndentLevel(tmpLevel)
    Writer.BlankLine()


def WriteImport(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)

    Writer.Statement("import os, sys")
    Writer.Statement("import numpy as np")
    Writer.Statement("import tensorflow as tf")
    Writer.Statement("import matplotlib.pyplot as plt")
    Writer.Statement("from datetime import datetime")
    Writer.Statement("from os import listdir")
    Writer.Statement("from argparse import ArgumentParser, FileType")
    Writer.Statement("import abc")

    Writer.SetIndentLevel(tmpLevel)
    Writer.BlankLine()


def WriteBody(Writer, CompilerData):

    Writer.Switch4DebugCompiler = False
    Writer.Switch4DebugSimulation = False
    Writer.Switch4Graph = False

    Writer.Variable_('LCCDataPath', "\"" + CompilerData.GetDataPath() + "\"")
    Writer.BlankLine()

    Lines = [
        "def main(GenomeFileName, verbose):",
    ]
    for Line in Lines:
        Writer.Statement(Line)

    with Writer:
        Writer.BlankLine()

        # Specify the mode of data processing
        Writer.Statement("# This is a numpy code", TargetCode=Target.Numpy)
        Writer.Statement("# This is a tensorflow code", TargetCode=Target.TensorFlow)

        # Load all object classes for variables
        Writer.Statement("# Load all object classes.")
        Simulation.Write_Simulation(Writer, CompilerData)
        Constant.Write_Constant(Writer, CompilerData)
        Environment.Write_Environment(Writer, CompilerData)
        CellState.Write_CellState(Writer, CompilerData)
        CellProcess.Write_CellProcess_Init(Writer, CompilerData)

        # GenomeState.Write_GenomeState_Init(Writer)
        # GeneState.Write_GeneState_Init(Writer)
        # PromoterState.Write_PromoterState_Init(Writer)
        # RNAState.Write_RNAState_Init(Writer)
        # ProteinMonomerState.Write_ProteinMonomerState_Init(Writer)
        # ComplexState.Write_ComplexState_Init(Writer)
        # LipidState.Write_LipidState_Init(Writer)
        Writer.BlankLine()

        # Declare all simulation components
        Writer.Statement("# Declare all simulation components.")
        Writer.Statement("Sim = FSimulation()")
        Writer.Statement("Cst = FConstant()")
        Writer.Statement("Env = FEnvironment()")
        Writer.Statement("Cel = FCellState()")

        # Writer.Statement("DNA = FGenomeState()") # Subclass of Cel
        # Writer.Statement("RNA = FRNAState()") # Subclass of Cel
        # Writer.Statement("PRT = FProteinMonomerState()") # Subclass of Cel
        # Writer.Statement("CPX = FComplexState()") # Subclass of Cel
        # Writer.Statement("LIP = FLipidState()") # Subclass of Cel
        # Writer.Statement("GEN = FGeneState()") # Subclass of DNA
        # Writer.Statement("PRM = FPromoterState()") # Subclass of DNA
        Writer.BlankLine()

        # Declare all processes
        Writer.Statement("# Declare all processes.")
        Writer.Statement("Dict_CellProcess = dict()")

        # Dict_ProcessTypes = dict()
        # Dict_ProcessTypes['Syn'] = 'Synthesis'
        # Dict_ProcessTypes['Deg'] = 'Degradation'
        # Dict_ProcessTypes['Mod'] = 'Modification'
        # Dict_ProcessTypes['Cnv'] = 'Conversion'
        # Dict_ProcessTypes['Met'] = 'Metabolism'

        # Initialize all processes
        Replication.Write_SynDNA(Writer, CompilerData)
        Transcription.Write_SynRNA(Writer, CompilerData)
        Translation.Write_SynPRT(Writer, CompilerData)
        DNADegradation.Write_DegDNA(Writer, CompilerData)
        RNADegradation.Write_DegRNA(Writer, CompilerData)
        ProteinDegradation.Write_DegPRT(Writer, CompilerData)


        # Central Dogma
        ProcessTypes = ['Syn', 'Deg'] # add 'Mod' later
        ProcessTargets = ['DNA', 'RNA', 'PRT']

        for ProcessType in ProcessTypes:
            for ProcessTarget in ProcessTargets:
                ProcessTypeTarget = ProcessType + ProcessTarget
                Writer.Statement("{0} = F{0}()".format(ProcessTypeTarget))
                Writer.Statement("Dict_CellProcess['{0}'] = {0}".format(ProcessTypeTarget))
                Writer.BlankLine()

        # Initialize all simulation components.
        Writer.Statement("# Initialize all simulation components.") # Only if Initialize not included by __init__()
        # Writer.Statement("Cst.Initialize()")
        # Writer.Statement("Env.Initialize()")

        # Initialize all processes.
        with Writer.Statement("for CellProcessKey, CellProcessValue in Dict_CellProcess.items():"):
            Writer.PrintVari("CellProcessKey")
            # Writer.Statement("CellProcessKey.Init()")

        # Writer.Statement("DNA = FGenomeState()") # Subclass of Cel
        # Writer.Statement("RNA = FRNAState()") # Subclass of Cel
        # Writer.Statement("PRT = FProteinMonomerState()") # Subclass of Cel
        # Writer.Statement("CPX = FComplexState()") # Subclass of Cel
        # Writer.Statement("LIP = FLipidState()") # Subclass of Cel
        # Writer.Statement("GEN = FGeneState()") # Subclass of DNA
        # Writer.Statement("PRM = FPromoterState()") # Subclass of DNA
        Writer.BlankLine()

        #

        # Load all simulation initialization functions

        # Temporary parameters
        #
        # E coli cell volume: 0.7 um3 (average), which is 7e-16 liters
        Writer.Variable_("Cel.Vol", 7e-16)  # TO BE REPLACED AND MOVED INTO SIMULATION
        Writer.BlankLine()


        # Run simulation
        Writer.Statement("# Run simulation")
        # Writer.Statement("Sim.RunSimulation(Cel, Cst, Env, )")

        Writer.BlankLine()
        # End of simulation

        # Print input genome
        with Writer.Statement("if GenomeFileName != \"\":"):
            Writer.Statement("GenomeFile = open(GenomeFileName, 'w')")
            Writer.Statement("InputGenomeFile = open('cell.fa')")
            with Writer.Statement("for Line in InputGenomeFile:"):
                Writer.Statement("Line = Line.strip()")
                Writer.Statement("print(Line, file=GenomeFile)")
            Writer.Statement("GenomeFile.close()")

def WriteMain(Writer):
    tmpLevel = Writer.GetIndentLevel()
    Writer.SetIndentLevel(0)

    Writer.BlankLine()
    with Writer.Statement("if __name__ == '__main__':"):
        Writer.Statement("parser = ArgumentParser(description='')")
        Writer.Statement("parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='also print some statistics to stderr')")
        Writer.Statement("parser.add_argument('-g', '--genome', dest='GenomeFileName', type=str, default='', help='')")
        Writer.Statement("parser.add_argument('-e', '--eager', dest='TFEagerMode', action='store_true', default=False, help='Running in the eager mode')")
        Writer.Statement("parser.add_argument('--trace-graph', dest='TFTraceGraph', action='store_true', default=False, help='Tracing the TensorFlow Graph')")
        Writer.BlankLine()
        Writer.Statement("args = parser.parse_args()")
        Writer.BlankLine()
        Writer.Statement("print('TF Eager Mode:', args.TFEagerMode)", TargetCode=Target.TensorFlow)
        Writer.Statement("tf.config.run_functions_eagerly(args.TFEagerMode)", TargetCode=Target.TensorFlow)
        Writer.BlankLine()
        with Writer.Statement("if args.TFTraceGraph:"):
            Writer.Statement("tf.summary.trace_on(graph=True)  # , profiler=True)")
        Writer.BlankLine()
        Writer.Statement("main(args.GenomeFileName, args.verbose)")
        Writer.BlankLine()

        with Writer.Statement("if args.TFTraceGraph:"):
            Writer.Statement("timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')")
            Writer.Statement("tf_logdir = 'logs/graph/%s' % timestamp")
            Writer.Statement("summary_writer = tf.summary.create_file_writer(tf_logdir)")
            with Writer.Statement("with summary_writer.as_default():"):
                Writer.Statement("tf.summary.trace_export(name='lcc_trace', step=0, profiler_outdir=tf_logdir)")
        Writer.BlankLine()

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
    return Seq


def Compile(CodeFileNames,
            PrefixName,
            DataDir,
            SaveDir,
            TargetCodeModel,
            Verbose):

    OutputCodeName = PrefixName + ".py"

    CompilerData = FCompilerData()
    CompilerData.SetDataPath(os.path.realpath(DataDir))
    CompilerData.SetSavePath(os.path.realpath(SaveDir))
    CodeInfo = Parse(CodeFileNames)

    GenomeFileName = PrefixName + ".ecoli.fa"
    if 'Organism' in CodeInfo:
        OrganismName = CodeInfo['Organism'].replace(' ', '_')
        GenomeFileName = PrefixName + "." + OrganismName + ".fa"

    CompilerData.InputGenomeSeq = CompileToGenome(GenomeFileName, DataDir)

    Dataset = CompilerData.LoadRawData(DataDir)
    CompilerData.InitializeCompilerData(Dataset)
    CompilerData.SetUpCompilerData(Dataset)
    CompilerData.SaveCompilerData()
    OutputFile = open(OutputCodeName, 'w')

    # TODO: use factory
    if TargetCodeModel == Target.TensorFlow:
        Writer = CodeGen.TFCodeWriter(OutputFile, 0)
    else:
        Writer = CodeGen.NumpyCodeWriter(OutputFile, 0)

    WriteLicense(Writer)
    WriteImport(Writer); Writer.BlankLine()
    WriteBody(Writer, CompilerData); Writer.BlankLine()
    WriteMain(Writer)

    OutputFile.close()
    print('"%s" has been successfully generated.' % OutputCodeName)

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
    parser.add_argument('-S',
                        dest='save_dir',
                        type=str,
                        default='',
                        help='Library/Save directory')
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
            args.save_dir,
            TargetCodeModel,
            args.verbose)
