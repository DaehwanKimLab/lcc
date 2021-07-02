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
from enum import Enum
import numpy as np
import os

class Target(Enum):
    """
    Target Codes
    """
    All = 0
    TensorFlow = 1
    Numpy = 2

    @staticmethod
    def from_str(label):
        tmp_label = label.upper()
        if tmp_label in ('TENSORFLOW', 'TF'):
            return Target.TensorFlow
        elif tmp_label in ('NUMPY', 'NP'):
            return Target.Numpy
        else:
            raise NotImplementedError


class CodeWriter():
    def __init__(self, CodeFile, IndentLevel=0):
        self.IndentLevel = IndentLevel
        self.fp = CodeFile
        self.BuildIndentationPrefix()

        self.Switch4Comment = False
        self.Switch4DebugSimulationPrint = False
        self.Switch4DebugSimulationAssert = False
        self.Switch4Graph = False

    def __enter__(self):
        self.IncreaseIndent()

    def __exit__(self, type, value, traceback):
        self.DecreaseIndent()

    def BuildIndentationPrefix(self):
        # 4 spaces
        self.IndentationPrefix = '    ' * self.IndentLevel

    def IncreaseIndent(self):
        self.IndentLevel += 1
        self.BuildIndentationPrefix()

    def DecreaseIndent(self):
        if self.IndentLevel > 0:
            self.IndentLevel -= 1

        self.BuildIndentationPrefix()

    def GetIndentLevel(self):
        return self.IndentLevel

    def SetIndentLevel(self, IndentLevel):
        self.IndentLevel = IndentLevel
        self.BuildIndentationPrefix()

    def Statement(self, Line, **kwargs):
        """
            keyword args
                TargetCode=Target.All, Target.Numpy, Target.TensorFlow
        """

        print("{Indent}{Line}".format(Indent=self.IndentationPrefix, Line=Line), file=self.fp)
        return self

    def Comment__(self, Comment):
        if self.Switch4Comment:
            Line = '# %s' % Comment
            self.Statement(Line)

    def BlankLine(self):
        self.Statement("")

    def Pass_____(self):
        self.Statement("pass")

    def AbsMethod(self):
        self.Statement("@abc.abstractmethod")

    def ReturnVar(self, VariableName):
        Line = 'return %s' % VariableName
        self.Statement(Line)

    def Overwrite(self, DestinationVariable, SourceVariable):
        Line = '%s = %s' % (DestinationVariable, SourceVariable)
        self.Statement(Line)

    def PrintStrg(self, Str):
        Line = 'print("%s")' % Str
        self.Statement(Line)

    def PrintVari(self, VariableName):
        Line = 'print("%s = ", %s)' % (VariableName, VariableName)
        self.Statement(Line)

    def PrintStVa(self, Str, VariableName):
        Line = 'print("%s: ", %s)' % (Str, VariableName)
        self.Statement(Line)

    def PrintDict(self, DictVariableName):
        with self.Statement("for Key, Value in %s.items():" % DictVariableName):
            self.PrintStVa("Keys and Values in ", DictVariableName)
            self.Statement('print("%s: %s" % (Key, Value))')

    def PrintKeys(self, DictVariableName):
        with self.Statement("for Key, Value in %s.items():" % DictVariableName):
            self.PrintStVa("Keys in ", DictVariableName)
            self.Statement('print("%s" % Key)')

    def PrintVals(self, DictVariableName):
        with self.Statement("for Key, Value in %s.items():" % DictVariableName):
            self.PrintStVa("Values in ", DictVariableName)
            self.Statement('print("%s" % Val)')

    def DebugSTMT(self, Line):
        if self.Switch4DebugSimulationPrint or self.Switch4DebugSimulationAssert:
            self.Statement(Line)

    def DebugVari(self, VariableName, Value):
        if self.Switch4DebugSimulationPrint or self.Switch4DebugSimulationAssert:
            Line = '%s = %s' % (VariableName, Value)
            self.Statement(Line)

    def DebugPStr(self, Str):
        if self.Switch4DebugSimulationPrint:
            Line = 'print("%s")' % Str
            self.Statement(Line)

    def DebugPVar(self, VariableName):
        if self.Switch4DebugSimulationPrint:
            Line = 'print("%s = ", %s)' % (VariableName, VariableName)
            self.Statement(Line)

    def DebugPStV(self, Str, VariableName):
        if self.Switch4DebugSimulationPrint:
            Line = 'print("%s: ", %s)' % (Str, VariableName)
            self.Statement(Line)

    def DebugAsrt(self, Condition, ErrMsg=None):
        if self.Switch4DebugSimulationAssert:
            Line = "assert %s, '%s'" % (Condition, ErrMsg)
            self.Statement(Line)


    """
    Array creation functions
    """

    def InitOnes_(self, VariableName, Shape, Type='float32'):
        # var = tf.ones([4,3])
        Line = "# Not implemented"
        self.Statement(Line)

    def InitZeros(self, VariableName, Shape, Type='float32'):
        # var = tf.zeros([4,3])
        Line = "# Not implemented"
        self.Statement(Line)

    def InitVals_(self, VariableName, Value, Type='float32'):
        # var = tf.constant([3, 4, 5, 6, 7])
        Line = "# Not implemented"
        self.Statement(Line)

    """
    Array manipulation functions
    """
    def Reshape(self, DestVar, SrcVar, Shape):
        # var = tf.reshape(var, [4, -1])
        Line = "# Not implemented"
        self.Statement(Line)


class TFCodeWriter(CodeWriter):
    """
    Tensorflow code generator
    """

    def __init__(self, CodeFile, IndentLevel=0):
        self.Target = Target.TensorFlow
        super(TFCodeWriter, self).__init__(CodeFile, IndentLevel)

    def Variable_(self, VariableName, Value):
        Line = '%s = tf.constant([%s])' % (VariableName, Value)
        self.Statement(Line)

    def VarCmnt__(self, VariableName, Value, Comment):
        Line = '%s = tf.constant([%s]) # %s' % (VariableName, Value, Comment)
        self.Statement(Line)


    def InitOnes_(self, VariableName, Shape, Type='float32'):
        Line = "{VariableName} = tf.ones({Shape}, dtype=tf.{Type})"\
            .format(VariableName=VariableName, Shape=str(Shape), Type=Type)
        self.Statement(Line)

    def InitZeros(self, VariableName, Shape, Type='float32'):
        Line = '%s = tf.zeros(%s, dtype=tf.%s)' % (VariableName, str(Shape), Type)
        self.Statement(Line)

    def InitVals_(self, VariableName, Value, Type='float32'):
        Line = '%s = tf.constant(%s, dtype=%s)' % (VariableName, str(Value), Type)
        self.Statement(Line)

    def AddElwise(self, VariableName, Value):
        Line = '%s = tf.add(%s, %s)' % (VariableName, VariableName, str(Value))
        self.Statement(Line)

    def LoadSaved(self, SaveFilePath, VariableName, DataType=None):
        FileType = SaveFilePath.split('.')[-1]
        if FileType == 'npy':
            VariableLoad = 'np.load(r"%s")' % (SaveFilePath)
            self.Convert__('self.%s' % VariableName, VariableLoad, DataType)
            self.DebugPStV(VariableName, 'self.%s' % VariableName)
        else:
            self.Statement('LoadSaved function cannot handle the filetype: %s' % FileType)

    def Convert__(self, VariableNameNew, VariableNamePrev=None, DataType='float32'):
        if VariableNamePrev:
            Line = '%s = tf.convert_to_tensor(%s, dtype=tf.%s)' % (VariableNameNew, VariableNamePrev, DataType)
        else:
            Line = '%s = tf.convert_to_tensor(%s)' % (VariableNameNew, VariableNameNew)
        self.Statement(Line)

    def Reshape__(self, DestVar, SrcVar, Shape):
        Line = '%s = tf.reshape(%s, %s)' % (DestVar, SrcVar, str(Shape))
        self.Statement(Line)
        self.DebugPVar(DestVar)

    def Transpose(self, Variable):
        Line = '%s = tf.transpose(%s)' % (Variable, Variable)
        self.Statement(Line)
        self.DebugPVar(Variable)

    def Statement(self, Line, **kwargs):

        if 'TargetCode' in kwargs:
            TargetCode = kwargs['TargetCode']
        else:
            TargetCode = Target.TensorFlow

        if TargetCode not in [Target.All, Target.TensorFlow]:
            return self

        super().Statement(Line)
        # print("{Indent}{Line}".format(Indent=self.IndentationPrefix, Line=Line), file=self.fp)
        return self

    def TF_Graph_(self):
        if self.Switch4Graph:
            self.Statement("@tf.function")

    def TFForLoop(self, Timer, MaxTime):
        pass

    def SelectIdx(self, VariableName, N_MoleculesToDistribute, Indexes, Weights='None'):
        self.InitZeros(VariableName, 0)
        if Weights:
            self.OperGathr("Weights", Weights, Indexes)
            self.Statement("Weights = Weights / len(Weights)")
        else:
            self.InitOnes_("Weights", len(Indexes), Type='int32')
        with self.Statement("for i in range(%s):" % N_MoleculesToDistribute):
            self.Statement("Values = tf.data.experimental.sample_from_datasets(%s, weights=Weights)" % Indexes)
            self.Statement("%s = tf.concat([%s, Value], 0)" % (VariableName, VariableName))
        self.DebugPStV("Values Generated:", VariableName)
        return VariableName

    def Increment(self, TargetMX, N_MoleculesToDistribute, Indexes, IncrementValue, WeightMX='None'):
        self.SelectIdx('Values', N_MoleculesToDistribute, Indexes, Weights=WeightMX)
        self.DebugPStV("Before Increment:", TargetMX)
        with self.Statement("for i in range(%s):" % N_MoleculesToDistribute):
            self.Statement("j = Values[i]")
            self.Statement("SelectedFromIndex = %s[j]" % Indexes)
            self.Statement("SelectedFromIndex = tf.reshape(SelectedFromIndex, [-1, 1])")
            self.Statement("TargetMXDataType = tf.shape(%s).dtype" % TargetMX)
            self.Statement("One = tf.ones(1, TargetMXDataType)")
            self.Statement("%s = tf.tensor_scatter_nd_add(%s, SelectedFromIndex, One * %s)" % (TargetMX, TargetMX, IncrementValue))
        self.DebugPStV("After Increment:", TargetMX)

    def OperCncat(self, DestVar, SrcVar1, SrcVar2, Axis=0):
        # self.PrepCncat(DestVar, SrcVar1, SrcVar2)
        self.Statement("%s = tf.concat([%s, %s], %s)" % (DestVar, SrcVar1, SrcVar2, Axis))

    def PrepGathr(self, Source, Index):
        self.Reshape__(Source, Source, -1)
        self.Reshape__(Index, Index, [-1, 1])

    def OperGathr(self, VariableName, Source, Index):
        self.PrepGathr(Source, Index)
        self.Statement("%s = tf.gather(%s, %s)" % (VariableName, Source, Index))

    def PrepScNds(self, Target, Index, Values):
        self.Reshape__(Target, Target, -1)
        self.Reshape__(Values, Values, -1)
        self.Reshape__(Index, Index, [-1, 1])

    def OperScAdd(self, Target, Index, Values):
        self.PrepScNds(Target, Index, Values)
        self.Statement("%s = tf.tensor_scatter_nd_add(%s, %s, %s)" % (Target, Target, Index, Values))

    def OperScSub(self, Target, Index, Values):
        self.PrepScNds(Target, Index, Values)
        self.Statement("%s = tf.tensor_scatter_nd_sub(%s, %s, %s)" % (Target, Target, Index, Values))

    def OperScUpd(self, Target, Index, Values):
        self.PrepScNds(Target, Index, Values)
        self.Statement("%s = tf.tensor_scatter_nd_update(%s, %s, %s)" % (Target, Target, Index, Values))

    def PrepMXMul(self, MX1, MX2):
        with self.Statement("if (tf.rank(%s) = 2) and (tf.rank(%s) = 2):" % (MX1, MX2)):
            self.Statement("if Row_%s = ")

        with self.Statement("if tf.rank(%s) != 2:" % MX1):
            self.Statement("%s = tf.reshape(%s, [-1, tf.shape(%s)[0]])" % (MX1, MX1, MX1))
        with self.Statement("if tf.rank(%s) != 2:" % MX2):
            self.Statement("%s = tf.reshape(%s, [tf.shape(%s)[0], -1])" % (MX2, MX2, MX1))

    def OperMXAdd(self, VariableName, MX1, MX2):
        self.Statement("%s = tf.math.add(%s, %s)" % (VariableName, MX1, MX2))

    def OperMXMul(self, VariableName, MX1, MX2):
        self.PrepMXMul(MX1, MX2)
        self.Statement("%s = tf.linalg.matmul(%s, %s)" % (VariableName, MX1, MX2))
        self.Reshape__(VariableName, VariableName, -1)

    def OperRound(self, VariableName, MX):
        self.Statement("%s = tf.math.round(%s)" % (VariableName, MX))

    def Cast_____(self, VariableName, MX, Type='int32'):
        self.Statement("%s = tf.cast(%s, dtype=tf.%s)" % (VariableName, MX, Type))


    # # To operate in CellState Class?
    # def ID2IdxSim(self, MolIDList, MolID2Index):
    #     self.InitArrayWithZero("MolIndexList", 0)
    #     with self.Statement("for Molecule in %s:" % MolIDList):
    #         self.Statement("MolIndex = %s(Molecule)" % MolID2Index)
    #         self.OperCncat("MolIndexList", 'MolIndex')
    #     self.DebugAsrt("tf.debugging.assert_shapes([MolIndexList, %s])" % (MolIDList))
    #     self.DebugPVar("MolIndexList")
    #     self.ReturnVar("MolIndexList")

    # To get indexes from molecule ID when using the compiler
    # def ID2Idx___(self, MolIndexList_Cel, MolIDList_Compiler, MolID2Idx_Compiler):
    #     self.InitArrayWithZero("%s" % MolIndexList_Cel, 0)
    #     for i, MoleculeID in enumerate(MolIDList_Compiler):
    #         # Item assignment only works in Graph
    #         self.Statement("%s[%s] = %s # %s" % (MolIndexList_Cel, i, MolID2Idx_Compiler[MoleculeID], MoleculeID))
    #     self.DebugAsrt("tf.debugging.assert_shapes([%s, %s])" % (MolIndexList_Cel, MolIDList_Compiler), ErrMsg='Incomplete Indexing for "%s"' % MolIDList_Compiler)
    #     self.DebugPVar("%s" % MolIndexList_Cel)

    def LoadNP2TF(self, VariableName, SavedFile, DataType='float32'):
        self.Statement("%s = np.load(\"%s\").astype(\'%s\')" % (VariableName, SavedFile, DataType))
        self.Convert__("%s" % VariableName)

    def Conc2Cont(self, VariableName, MolIndexList, MolConcs):
        self.PrepGathr(MolConcs, MolIndexList)
        self.OperGathr("MolConcs4MolIndexList", MolConcs, MolIndexList)
        self.Statement("%s = MolConcs4MolIndexList * (Cel.CellVol * Cst.NA)" % VariableName)

    def Cont2Conc(self, VariableName, MolIndexList, MolCounts):
        self.PrepGathr(MolCounts, MolIndexList)
        self.OperGathr("MolCounts4MolIndexList", MolCounts, MolIndexList)
        self.Statement("%s = MolCounts4MolIndexList / (Cel.CellVol * Cst.NA)" % VariableName)

    # WRONG EQUATION - to be fixed with %s and more
    # def Cont2Conc(self, MolIDList, Mol2Index, CountsMX, MWsMX):
    #     self.InitArrayWithZero("MolIndices", 0)
    #     with self.Statement("for Molecule in %s:" % MolIDList):
    #         self.Statement("MolIndex = MolNames2Index(Molecule)")
    #         self.Statement("MolIndices = tf.concat([MolIndices, MolIndex], 0)")
    #     self.ReshapeMX("CountsMX", "CountsMX", -1)
    #     self.ReshapeMX("MWsMX", "MWsMX", -1)
    #     self.ReshapeMX("MolIndices", "MolIndices", [-1, 1])
    #     self.Statement("CountsForMolIndices = tf.gather(CountsMX, MolIndices)")
    #     self.Statement("MWsForMolIndices = tf.gather(MWsMX, MolIndices)")
    #     self.PrepMXMul("CountsForMolIndices", "MWsForMolIndices")
    #     self.Statement("ConcsForMolIndices = tf.matmul(CountsForMolIndices, MWsForMolIndices) * ")
    #
    #     self.InitArrayWithZero("MWs", 0)
    #     with self.Statement("for Molecule in MolIDList:"):
    #         self.Statement("MolIndex = MolNames2Index(Molecule)")
    #         self.Statement("MWs = tf.concat([MolIndices, MolIndex], 0))")
    #         self.Statement("matmul(Counts, MWs")
    #
    #     return self


class NumpyType:
    """
    Numpy type as string
    """
    TypeToString = dict({'float32': '"float32"', 'int32': '"int32"'})


class NumpyCodeWriter(CodeWriter):
    """
    Numpy code generator
    """

    def __init__(self, CodeFile, IndentLevel=0):
        self.Target = Target.Numpy
        super(NumpyCodeWriter, self).__init__(CodeFile, IndentLevel)

    def Variable_(self, VariableName, Value):
        Line = '%s = %s' % (VariableName, Value)
        self.Statement(Line)

    def InitOnes_(self, VariableName, Shape, Type='float32'):
        Line = "{VariableName} = np.ones({Shape}).astype({Type})"\
            .format(VariableName=VariableName, Shape=str(Shape), Type=NumpyType.TypeToString[Type])
        self.Statement(Line)

    def InitZeros(self, VariableName, Shape, Type='float32'):
        Line = "%s = np.zeros(%s).astype(%s)" % (VariableName, str(Shape), NumpyType.TypeToString[Type])
        self.Statement(Line)

    def InitArrayWithValue(self, VariableName, Value, Type='float32'):
        Line = "%s = np.array(%s).astype(%s)" % (VariableName, str(Value), Type)
        self.Statement(Line)

    def Reshape(self, DestVar, SrcVar, Shape):
        Line = '%s = np.reshape(%s, %s)' % (DestVar, SrcVar, str(Shape))
        self.Statement(Line)

    def Statement(self, Line, **kwargs):
        if 'TargetCode' in kwargs:
            TargetCode = kwargs['TargetCode']
        else:
            TargetCode = Target.Numpy

        if TargetCode not in [Target.All, Target.Numpy]:
            return

        super().Statement(Line)
        return self
