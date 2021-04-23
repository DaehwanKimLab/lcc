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
        self.DebugPrintSwitch = False
        self.DebugAssertSwitch = False

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

    def Variable_(self, VariableName, Value):
        Line = '%s = %s' % (VariableName, Value)
        self.Statement(Line)

    # def Comment(self, Comment):
    #     Line = '# %s' % Comment
    #     self.Statement(Line)
    #     return self

    def BlankLine(self):
        self.Statement("")

    def PrintStrg(self, Str):
        Line = 'print("%s")' % Str
        self.Statement(Line)
        return self

    def PrintVari(self, VariableName):
        Line = 'print("%s = ", %s)' % (VariableName, VariableName)
        self.Statement(Line)
        return self

    def PrintStVa(self, Str, VariableName):
        Line = 'print("%s: ", %s)' % (Str, VariableName)
        self.Statement(Line)
        return self

    def DebugSTMT(self, Line):
        if self.DebugPrintSwitch or self.DebugAssertSwitch:
            self.Statement(Line)
        return self

    def DebugVari(self, VariableName, Value):
        if self.DebugPrintSwitch or self.DebugAssertSwitch:
            Line = '%s = %s' % (VariableName, Value)
            self.Statement(Line)
        return self

    def DebugPStr(self, Str):
        if self.DebugPrintSwitch:
            Line = 'print("%s")' % Str
            self.Statement(Line)
        return self

    def DebugPVar(self, VariableName):
        if self.DebugPrintSwitch:
            Line = 'print("%s = ", %s)' % (VariableName, VariableName)
            self.Statement(Line)
        return self

    def DebugPStV(self, Str, VariableName):
        if self.DebugPrintSwitch:
            Line = 'print("%s: ", %s)' % (Str, VariableName)
            self.Statement(Line)
        return self

    def DebugAsrt(self, Condition, ErrMsg):
        if self.DebugAssertSwitch:
            Line = "assert %s, '%s'" % (Condition, ErrMsg)
            self.Statement(Line)
        return self


    """
    Array creation functions
    """

    def InitArrayWithOne(self, VariableName, Shape, Type='float32'):
        # var = tf.ones([4,3])
        Line = "# Not implemented"
        self.Statement(Line)

    def InitArrayWithZero(self, VariableName, Shape, Type='float32'):
        # var = tf.zeros([4,3])
        Line = "# Not implemented"
        self.Statement(Line)

    def InitArrayWithValue(self, VariableName, Value, Type='float32'):
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

    def InitArrayWithOne(self, VariableName, Shape, Type='float32'):
        Line = "{VariableName} = tf.ones({Shape}, dtype=tf.{Type})"\
            .format(VariableName=VariableName, Shape=str(Shape), Type=Type)
        self.Statement(Line)

    def InitArrayWithZero(self, VariableName, Shape, Type='float32'):
        Line = '%s = tf.zeros(%s, dtype=tf.%s)' % (VariableName, str(Shape), Type)
        self.Statement(Line)

    def InitArrayWithValue(self, VariableName, Value, Type='float32'):
        Line = '%s = tf.constant(%s, dtype=%s)' % (VariableName, str(Value), Type)
        self.Statement(Line)

    def ReshapeMX(self, DestVar, SrcVar, Shape):
        Line = '%s = tf.reshape(%s, %s)' % (DestVar, SrcVar, str(Shape))
        self.Statement(Line)

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

    # Random Value Generator
    def RndValues(self, VariableName, Size, Max):
        Line = '%s = tf.random.uniform(shape=[%s], minval=0, maxval=%s, dtype="int32")' % (VariableName, Size, Max)
        self.Statement(Line)

    # Randomly chosen molecules for finding molecular interaction - maybe more relevant to numpy later
    def RndIncrmt(self, TargetMX, NumberOfMoleculesToDistribute, Index, IncrementValue):
        self.RndValues('RandValuePicks', NumberOfMoleculesToDistribute, 'len(' + Index + ')')
        with self.Statement("for i in range(%s):" % NumberOfMoleculesToDistribute):
            self.Statement("j = RandValuePicks[i]")
            self.Statement("SelectedFromIndex = %s[j]" % Index)
            self.Statement("SelectedFromIndex = tf.reshape(SelectedFromIndex, [-1, 1])")
            self.Statement("TargetMXDataType = tf.shape(%s).dtype" % TargetMX)
            self.Statement("One = tf.ones(1, TargetMXDataType)")
            self.Statement("%s = tf.tensor_scatter_nd_add(%s, SelectedFromIndex, One * %s)" % (TargetMX, TargetMX, IncrementValue))

    # Weighted Value generator
    def WgtValues(self, VariableName, Size, Index, WeightMX):
        self.InitArrayWithZero(VariableName, 0)
        self.Statement("Weights = tf.gather(%s, %s)" % (WeightMX, Index))
        self.Statement("Weights = Weights / len(Weights)")
        with self.Statement("for i in range(%s):" % Size):
            self.Statement("WeightedIndexSelected = tf.data.experimental.sample_from_datasets(%s, weights=Weights)" % Index)
            self.Statement("%s = tf.concat([%s, WeightedIndexSelected], 0)" % (VariableName, VariableName))

    def WgtIncrmt(self, TargetMX, NumberOfMoleculesToDistribute, Index, IncrementValue):
        self.WgtValues('WeightedValuePicks', NumberOfMoleculesToDistribute, Index, TargetMX)
        with self.Statement("for i in range(%s):" % NumberOfMoleculesToDistribute):
            self.Statement("j = WeightedValuePicks[i]")
            self.Statement("SelectedFromIndex = %s[j]" % Index)
            self.Statement("SelectedFromIndex = tf.reshape(SelectedFromIndex, [-1, 1])")
            self.Statement("TargetMXDataType = tf.shape(%s).dtype" % TargetMX)
            self.Statement("One = tf.ones(1, TargetMXDataType)")
            self.Statement("%s = tf.tensor_scatter_nd_add(%s, SelectedFromIndex, One * %s)" % (TargetMX, TargetMX, IncrementValue))


    def PrepGathr(self, Target, Index):
        self.Statement("%s = tf.reshape(%s, -1)" % (Target, Target))
        self.Statement("%s = tf.reshape(%s, [-1, 1])" % (Index, Index))

    def OperGathr(self, VariableName, Target, Index):
        self.PrepGathr(Target, Index)
        self.Statement("%s = tf.gather(%s, %s)" % (VariableName, Target, Index))

    def PrepScNds(self, Target, Index, Values):
        self.Statement("%s = tf.reshape(%s, -1)" % (Target, Target))
        self.Statement("%s = tf.reshape(%s, -1)" % (Values, Values))
        self.Statement("%s = tf.reshape(%s, [-1, 1])" % (Index, Index))

    def OperScAdd(self, Target, Index, Values):
        self.PrepScNds(Target, Index, Values)
        self.Statement("%s = tf.tensor_scatter_nd_add(%s, %s, %s)" % (Target, Target, Index, Values))

    def OperScSub(self, Target, Index, Values):
        self.PrepScNds(Target, Index, Values)
        self.Statement("%s = tf.tensor_scatter_nd_sub(%s, %s, %s)" % (Target, Target, Index, Values))

    def PrepMXMul(self, MX1, MX2):
        with self.Statement("if tf.rank(%s) != 2:" % MX1):
            self.Statement("%s = tf.reshape(%s, [-1, tf.shape(%s)[0]])" % (MX1, MX1, MX1))
        with self.Statement("if tf.rank(%s) != 2:" % MX2):
            self.Statement("%s = tf.reshape(%s, [tf.shape(%s)[0], -1])" % (MX2, MX2, MX1))

    def OperMXMul(self, VariableName, MX1, MX2):
        self.PrepMXMul(MX1, MX2)
        self.Statement("%s = tf.linalg.matmul(%s, %s)" % (VariableName, MX1, MX2))
        self.Statement("%s = tf.reshape(%s, -1)" % (VariableName, VariableName))

    def Ls2Index_(self, MolIndexList, MolList, Mol2Index):
        self.InitArrayWithZero("%s" % MolIndexList)
        with self.Statement("for Molecule in %s:" % MolList):
            self.Statement("MolIndex = %s(Molecule)" % Mol2Index)
            self.Statement("%s = tf.concat([%s, MolIndex], 0)" % (MolIndexList, MolIndexList))
        self.DebugAsrt("tf.debugging.assert_shapes([%s, %s])" % (MolIndexList, MolList))
        self.DebugPVar("%s" % MolIndexList)

    def Conc2Cont(self, VariableName, MolList, Mol2IndexDict, MolIndexList, MolConcs):
        self.Ls2Index_(MolIndexList, MolList, Mol2IndexDict)
        self.PrepGathr(MolConcs, MolIndexList)
        self.OperGathr("MolConcs4MolIndexList", MolConcs, MolIndexList)
        self.Statement("%s = MolConcs4MolIndexList * CellMX.CellVol * AvogadroNum" % VariableName)

    def Cont2Conc(self, VariableName, MolList, Mol2IndexDict, MolIndexList, MolCounts):
        self.Ls2Index_(MolIndexList, MolList, Mol2IndexDict)
        self.PrepGathr(MolCounts, MolIndexList)
        self.OperGathr("MolCounts4MolIndexList", MolCounts, MolIndexList)
        self.Statement("%s = MolCounts4MolIndexList / CellMX.CellVol / AvogadroNum" % VariableName)

    # WRONG EQUATION - to be fixed with %s and more
    # def Cont2Conc(self, MolList, Mol2Index, CountsMX, MWsMX):
    #     self.InitArrayWithZero("MolIndices", 0)
    #     with self.Statement("for Molecule in %s:" % MolList):
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
    #     with self.Statement("for Molecule in MolList:"):
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

    def InitArrayWithOne(self, VariableName, Shape, Type='float32'):
        Line = "{VariableName} = np.ones({Shape}).astype({Type})"\
            .format(VariableName=VariableName, Shape=str(Shape), Type=NumpyType.TypeToString[Type])
        self.Statement(Line)

    def InitArrayWithZero(self, VariableName, Shape, Type='float32'):
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
