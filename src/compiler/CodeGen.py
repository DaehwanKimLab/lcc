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
        self.Switch4PrintString = False
        self.Switch4TFGraph = False

        # Metabolism
        self.Switch4PostSimulationStepCorrection = False
        self.Switch4Kinetics = False

        # Cell Division Test
        self.Switch4TestCellDivision = False

        # Print switches
        self.Switch4SimStepsExecuted = False
        self.Switch4ProcessSummary = False
        self.Switch4CellStateSummary = False
        self.Switch4PostSimulationStepCorrectionMessage = False
        self.Switch4Visualization2D = False

        # Debugging
        self.Switch4DebugSimulationPrint = False
        self.Switch4DebugSimulationAssert = False
        self.Switch4ProcessDebuggingMessages = False
        self.Switch4TFAssert = False
        self.Switch4SoftCheckCounts = False
        self.Switch4HardCheckCounts = False
        self.Switch4CheckDeltaCountsNeg = False
        self.Switch4ShowDeltaCounts = False
        self.Switch4ProcessTimer = False

        # Save Data
        self.Switch4Save = False
        self.Switch4SaveAllCounts = False
        self.Switch4SaveSpecificCounts = False


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

    def PrintLine(self, Length=100):
        Unit = '-'
        Line = 'print("%s")' % (Unit * Length)
        self.Statement(Line)

    def BlankLine(self):
        self.Statement("")

    def Pass_____(self):
        self.Statement("pass")

    def AbsMethod(self):
        self.Statement("@abc.abstractmethod")

    def LinkClObj(self, VariableName):
        Line = 'self.%s = %s' % (VariableName, VariableName)
        self.Statement(Line)

    def ReturnVar(self, *args):
        line = 'return'
        for VariableName in args:
            line += (' ' + str(VariableName) + ',')
        line = line[:-1]
        self.Statement(line)

    def Overwrite(self, DestinationVariable, SourceVariable):
        Line = '%s = %s' % (DestinationVariable, SourceVariable)
        self.Statement(Line)

    def PrintStrg(self, Str):
        if self.Switch4PrintString:
            Line = 'print("%s")' % Str
            self.Statement(Line)

    def PrintVar_(self, VariableName):
        if self.Switch4PrintString:
            Line = 'print(\t%s)' % VariableName
            self.Statement(Line)

    def PrintVaVa(self, VariableName):
        if self.Switch4PrintString:
            Line = 'print("\t%s = ", %s)' % (VariableName, VariableName)
            self.Statement(Line)

    def PrintStVa(self, Str, VariableName):
        if self.Switch4PrintString:
            self.Statement('print("%s: ")' % Str)
            self.Statement('print("\t", %s)' % VariableName)

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

    def Function_(self, FunctionName, *Args):
        if self.Switch4TFGraph:
            self.Statement("@tf.function")
        Line = "def {}".format(FunctionName)
        Line += "("
        Line += "self"
        if not Args:
            pass
        else:
            for Arg in Args:
                Line += ", " + Arg

        Line += "):"
        return self.Statement(Line)

    DefaultTensor = "[tf.constant(0)]"

    def WhileLoop(self, Condition, Body, LoopVariable=DefaultTensor):
        Line = 'tf.while_loop(%s, %s, %s)' % (Condition, Body, LoopVariable)
        self.Statement(Line)

    def WhileRtrn(self, ReturnVariable=DefaultTensor):
        self.ReturnVar(ReturnVariable)

    def IfElse___(self, Predicate, True_Function, False_Function):
        Line = 'tf.cond(%s, true_fn=%s, false_fn=%s)' % (Predicate, True_Function, False_Function)
        self.Statement(Line)

    # tf.print
    def PrintLine(self, Length=100):
        Unit = '-'
        Line = 'tf.print("%s", output_stream=sys.stdout)' % (Unit * Length)
        self.Statement(Line)

    def PrintStrg(self, Str):
        if self.Switch4PrintString:
            Line = 'tf.print("%s", output_stream=sys.stdout)' % Str
            self.Statement(Line)
        else:
            self.Pass_____()
            self.BlankLine()

    def PrintVar_(self, VariableName):
        if self.Switch4PrintString:
            Line = 'tf.print(\t%s, output_stream=sys.stdout)' % VariableName
            self.Statement(Line)
        else:
            self.Pass_____()
            self.BlankLine()

    def PrintVaVa(self, VariableName):
        if self.Switch4PrintString:
            Line = 'tf.print("\t%s = ", %s, output_stream=sys.stdout)' % (VariableName, VariableName)
            self.Statement(Line)
        else:
            self.Pass_____()
            self.BlankLine()

    def PrintStVa(self, Str, VariableName):
        if self.Switch4PrintString:
            # For a double liner
            self.Statement('tf.print("%s: ", output_stream=sys.stdout)' % Str)
            self.Statement('tf.print("\t", %s, output_stream=sys.stdout)' % VariableName)

            # For a single liner
            # self.Statement('tf.print("%s:\t", %s, output_stream=sys.stdout)' % (Str, VariableName))
        else:
            self.Pass_____()
            self.BlankLine()

    def Variable_(self, VariableName, Value, Type=None, Shape=None):
        Line = '%s = tf.constant([%s], dtype=%s, shape=%s)' % (VariableName, Value, Type, Shape)
        self.Statement(Line)

    def VarCmnt__(self, VariableName, Value, Comment):
        Line = '%s = tf.constant([%s]) # %s' % (VariableName, Value, Comment)
        self.Statement(Line)

    def VarRange_(self, VariableName, StartValue, LimitValue, Delta=1):
        Line = '%s = tf.range(%s, %s, delta=%s)' % (VariableName, str(StartValue), str(LimitValue), Delta)
        self.Statement(Line)

    def InitZeros(self, VariableName, Shape, Type='float32'):
        Line = '%s = tf.zeros(%s, dtype=tf.%s)' % (VariableName, str(Shape), Type)
        self.Statement(Line)

    def ZerosLike(self, VariableName, SrcVar, Type='float32'):
        Line = '%s = tf.zeros_like(%s, dtype=tf.%s)' % (VariableName, SrcVar, Type)
        self.Statement(Line)

    def InitOnes_(self, VariableName, Shape, Type='float32'):
        Line = "{VariableName} = tf.ones({Shape}, dtype=tf.{Type})"\
            .format(VariableName=VariableName, Shape=str(Shape), Type=Type)
        self.Statement(Line)

    def OnesLike_(self, VariableName, SrcVar, Type='float32'):
        Line = '%s = tf.ones_like(%s, dtype=tf.%s)' % (VariableName, SrcVar, Type)
        self.Statement(Line)

    def InitNOnes(self, VariableName, Shape, Type='float32'):
        Line = "{VariableName} = tf.math.multiply(tf.ones({Shape}, dtype=tf.{Type})"\
            .format(VariableName=VariableName, Shape=str(Shape), Type=Type)
        self.Statement(Line)

    def InitVals_(self, VariableName, Value, Type='float32'):
        Line = '%s = tf.constant(%s, dtype=tf.%s)' % (VariableName, str(Value), Type)
        self.Statement(Line)

    def VarRepeat(self, VariableName, Input, Repeats, Axis=None):
        Line = '%s = tf.repeat(%s, %s, axis=%s)' % (VariableName, Input, Repeats, Axis)
        self.Statement(Line)

    def VarFill__(self, VariableName, Shape, ScalarValue):
        Line = '%s = tf.fill(%s, %s)' % (VariableName, Shape, ScalarValue)
        self.Statement(Line)

    def VarTile__(self, VariableName, Input, Multiples):
        Line = '%s = tf.tile(%s, %s)' % (VariableName, Input, Multiples)
        self.Statement(Line)

    def Shift____(self, VariableName, Input, Shift, Axis):
        Line = '%s = tf.roll(%s, %s, %s)' % (VariableName, Input, Shift, Axis)
        self.Statement(Line)

    def ClearVal_(self, VariableName, Type='float32'):
        Line = '%s = tf.constant(0, dtype=tf.%s)' % (VariableName, Type)
        self.Statement(Line)

    def Shape____(self, DestVar, Input):
        Line = '%s = tf.shape(%s)' % (DestVar, Input)
        self.Statement(Line)

    def ShapeAxis(self, DestVar, Input, Axis):
        Line = '%s = tf.shape(%s)[%s]' % (DestVar, Input, Axis)
        self.Statement(Line)

    def LoadSaved(self, SaveFilePath, VariableName, DataType=None):
        FileType = SaveFilePath.split('.')[-1]
        if FileType == 'npy':
            VariableLoad = 'np.load(r"%s")' % (SaveFilePath)
            self.Convert__('self.%s' % VariableName, VariableLoad, DataType)
            self.DebugPStV(VariableName, 'self.%s' % VariableName)
        else:
            self.Statement('LoadSaved function cannot handle the filetype: %s' % FileType)

    def Convert__(self, VariableNameNew, VariableNamePrev, DataType='float32'):
        Line = '%s = tf.convert_to_tensor(%s, dtype=tf.%s)' % (VariableNameNew, VariableNamePrev, DataType)
        self.Statement(Line)

    # TODO: Make functions with reshaping within itself as an option. (e.g. Output = tf.reshape(tf.transpose(Input), Shape=[]))
    def Reshape__(self, DestVar, SrcVar, Shape):
        Line = '%s = tf.reshape(%s, %s)' % (DestVar, SrcVar, str(Shape))
        self.Statement(Line)
        self.DebugPVar(DestVar)

    def ReshType_(self, DestVar, SrcVar, Shape, DataType):
        self.Statement('%s = tf.cast(tf.reshape(%s, %s), dtype=tf.%s)' % (DestVar, SrcVar, str(Shape), DataType))
        self.DebugPVar(DestVar)

    def Transpose(self, DestVar, SrcVar):
        Line = '%s = tf.transpose(%s)' % (DestVar, SrcVar)
        self.Statement(Line)

    def Negative_(self, DestVar, SrcVar):
        self.Statement('%s = tf.math.negative(%s)' % (DestVar, SrcVar))
        self.DebugPVar(DestVar)

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

    def IdxRange_(self, VariableName, IdxSequences):
        self.VarRange_(VariableName, IdxSequences[0], IdxSequences[-1])

    def RndNumUni(self, VariableName, Shape=0, MinVal=0, MaxVal=1, Type='int32'):
        self.Statement("%s = tf.random.uniform(shape=%s, minval=%s, maxval=%s, dtype=tf.%s)" % (VariableName, Shape, MinVal, MaxVal, Type))

    def GreaterEq(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.greater_equal(%s, %s)" % (DestVar, MX1, MX2))

    def Greater__(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.greater(%s, %s)" % (DestVar, MX1, MX2))

    def LessEq___(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.less_equal(%s, %s)" % (DestVar, MX1, MX2))

    def Less_____(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.less(%s, %s)" % (DestVar, MX1, MX2))

    def Equal____(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.equal(%s, %s)" % (DestVar, MX1, MX2))

    def NotEqual_(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.not_equal(%s, %s)" % (DestVar, MX1, MX2))

    def Concat___(self, DestVar, SrcVar1, SrcVar2, Axis=0):
        # self.PrepCncat(DestVar, SrcVar1, SrcVar2)
        self.Statement("%s = tf.concat([%s, %s], %s)" % (DestVar, SrcVar1, SrcVar2, Axis))

    def PrepGathr(self, Source, Index):
        self.Reshape__("Source", Source, -1)
        self.Reshape__("Index", Index, [-1, 1])

    def Gather___(self, VariableName, Source, Index):
        self.Comment__("Reshape matrices for a gather operation")
        self.PrepGathr(Source, Index)
        # self.Statement("Source = tf.reshape(%s, [-1])   # Reshape source matrix for gather function" % Source)
        # self.Statement("Index = tf.reshape(%s, [-1, 1])   # Reshape index matrix for gather function" % Index)
        self.Statement("%s = tf.gather(%s, %s)" % (VariableName, "Source", "Index"))

    def PrepScNds(self, Target, Index, Values, Type):
        self.Comment__("Reshape matrices for a scatter operation")
        self.ReshType_("Target", Target, -1, Type)
        self.ReshType_("Index", Index, [-1, 1], 'int32')
        self.ReshType_("Values", Values, -1, Type)

    def ScatNdAdd(self, DestVar, Target, Index, Values, Type='int32'):
        self.PrepScNds(Target, Index, Values, Type)
        self.Statement("%s = tf.tensor_scatter_nd_add(%s, %s, %s)" % (DestVar, "Target", "Index", "Values"))

    def ScatNdSub(self, DestVar, Target, Index, Values, Type='int32'):
        self.PrepScNds(Target, Index, Values, Type)
        self.Statement("%s = tf.tensor_scatter_nd_sub(%s, %s, %s)" % (DestVar, "Target", "Index", "Values"))

    def ScatNdUpd(self, DestVar, Target, Index, Values, Type='int32'):
        self.PrepScNds(Target, Index, Values, Type)
        self.Statement("%s = tf.tensor_scatter_nd_update(%s, %s, %s)" % (DestVar, "Target", "Index", "Values"))

    def LogicAnd_(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.logical_and(%s, %s)" % (DestVar, MX1, MX2))
        # self.Statement("%s = %s & %s)" % (DestVar, MX1, MX2))

    def LogicOr__(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.logical_or(%s, %s)" % (DestVar, MX1, MX2))
        # self.Statement("%s = %s | %s)" % (DestVar, MX1, MX2))

    def Add______(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.add(%s, %s)" % (DestVar, MX1, MX2))

    def Subtract_(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.subtract(%s, %s)" % (DestVar, MX1, MX2))

    def Multiply_(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.multiply(%s, %s)" % (DestVar, MX1, MX2))

    def Divide___(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.divide(%s, %s)" % (DestVar, MX1, MX2))

    def FloorDiv_(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.floordiv(%s, %s)" % (DestVar, MX1, MX2))

    def Remainder(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.floormod(%s, %s)" % (DestVar, MX1, MX2))

    def Maximum__(self, DestVar, MX1, MX2):
        self.Statement("%s = tf.math.maximum(%s, %s)" % (DestVar, MX1, MX2))
        
    def Minimum__(self, DestVar, MX1, MX2):
        self.Comment__("Element-wise minimum operation")
        self.Statement("%s = tf.math.minimum(%s, %s)" % (DestVar, MX1, MX2))
        
    # def PrepMXMul(self, MX1, MX2):
    #     self.Comment__("Reshape matrices for matrix multiplication")
    #     with self.Statement("if tf.rank(%s) != 2:" % MX1):
    #         self.Statement("%s = tf.reshape(%s, [-1, tf.shape(%s)[0]])" % (MX1, MX1, MX1))
    #     with self.Statement("if tf.rank(%s) != 2:" % MX2):
    #         self.Statement("%s = tf.reshape(%s, [tf.shape(%s)[0, -1])" % (MX2, MX2, MX1))

    def MatrixMul(self, DestVar, MX1, MX2):
        # self.PrepMXMul(MX1, MX2)
        self.Statement("%s = tf.linalg.matmul(%s, %s)" % (DestVar, MX1, MX2))
        self.Reshape__(DestVar, DestVar, -1)

    def ReduceSum(self, DestVar, MX, Axis=None):
        self.Statement("%s = tf.math.reduce_sum(%s, axis=%s)" % (DestVar, MX, Axis))

    def ReduceMul(self, DestVar, MX, Axis=None):
        self.Statement("%s = tf.math.reduce_prod(%s, axis=%s)" % (DestVar, MX, Axis))

    def ReduceMax(self, DestVar, MX, Axis=None):
        self.Statement("%s = tf.math.reduce_max(%s, axis=%s)" % (DestVar, MX, Axis))

    def ReduceMin(self, DestVar, MX, Axis=None):
        self.Statement("%s = tf.math.reduce_min(%s, axis=%s)" % (DestVar, MX, Axis))

    def ReduceAvg(self, DestVar, MX, Axis=None):
        self.Statement("%s = tf.math.reduce_mean(%s, axis=%s)" % (DestVar, MX, Axis))

    def ReduceAll(self, DestVar, MX, Axis=None):
        # For boolean only
        self.Statement("%s = tf.math.reduce_all(%s, axis=%s)" % (DestVar, MX, Axis))

    def ReduceAny(self, DestVar, MX, Axis=None):
        # For boolean only
        self.Statement("%s = tf.math.reduce_any(%s, axis=%s)" % (DestVar, MX, Axis))

    def CumSum___(self, DestVar, MX, Axis=0, Exclusive=False, Reverse=False):
        # Cumulative sum
        self.Statement("%s = tf.math.cumsum(%s, axis=%s, exclusive=%s, reverse=%s)" % (DestVar, MX, Axis, Exclusive, Reverse))

    def NonZeros_(self, DestVar, MX, Axis=None):
        self.Statement("%s = tf.math.count_nonzero(%s, axis=%s)" % (DestVar, MX, Axis))

    def Round____(self, DestVar, MX):
        # Round without changing data type
        self.Statement("%s = tf.math.round(%s)" % (DestVar, MX))

    def Cast_____(self, DestVar, MX, Type='int32'):
        # Change data type
        self.Statement("%s = tf.cast(%s, dtype=tf.%s)" % (DestVar, MX, Type))

    def BoolToBin(self, BinMX, BoolMX):
        self.Cast_____(BinMX, BoolMX)

    def BoolMask_(self, DestVar, TargetMX, BoolMX, Axis=None):
        self.Statement("%s = tf.boolean_mask(%s, %s, axis=%s)" % (DestVar, TargetMX, BoolMX, Axis))

    def ToSparse_(self, DestVar, MX):
        pass
        # self.Statement("%s = tf.sparse.from_dense(%s)" % (DestVar, MX))

    def CondVal__(self, DestVar, MX, Equality, Value, IfTrue=None, IfFalse=None):
        self.Statement("%s = tf.where(%s %s %s, x=%s, y=%s)" % (DestVar, MX, str(Equality), Value, IfTrue, IfFalse))

    def ConvToBin(self, DestVar, MX, Equality, Value):
        self.Statement("%s = tf.where(%s %s %s, 1, 0)" % (DestVar, MX, str(Equality), Value))

    def GenIdx___(self, DestVar, MX):
        self.Statement("%s = tf.where(%s)" % (DestVar, MX))

    def GenIdxCnd(self, DestVar, MX, Equality, Value):
        self.Statement("%s = tf.where(%s %s %s)" % (DestVar, MX, str(Equality), Value))

    def ArgMax___(self, DestVar, MX, Axis=None, Type='int32'):
        self.Statement("%s = tf.math.argmax(%s, %s, output_type=tf.%s)" % (DestVar, MX, Axis, Type))

    def ArgMin___(self, DestVar, MX, Axis=None, Type='int32'):
        self.Statement("%s = tf.math.argmin(%s, %s, output_type=tf.%s)" % (DestVar, MX, Axis, Type))

    def Sort_____(self, DestVar, MX, Axis=1, Direction='DESCENDING'):
        # Direction may be either 'ASCENDING' or 'DESCENDING'
        self.Statement("%s = tf.sort(%s, axis=%s, direction=%s)" % (DestVar, MX, Axis, Direction))

    def RoundInt_(self, DestVar, MX, Type='int32'):
        self.Statement("%s = tf.cast(tf.math.round(%s), dtype=tf.%s)" % (DestVar, MX, Type))

    def Sigmoid__(self, DestVar, MX):
        self.Statement("%s = tf.math.sigmoid(%s))" % (DestVar, MX))

    def AsrtGrEq_(self, MX1, MX2, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_greater_equal(%s, %s, message=%s)" % (MX1, MX2, Message))

    def AsrtGr___(self, MX1, MX2, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_greater(%s, %s, message=%s)" % (MX1, MX2, Message))

    def AsrtLeEq_(self, MX1, MX2, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_less_equal(%s, %s, message=%s)" % (MX1, MX2, Message))

    def AsrtLe___(self, MX1, MX2, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_less(%s, %s, message=%s)" % (MX1, MX2, Message))

    def AsrtEq___(self, MX1, MX2, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_equal(%s, %s, message=%s)" % (MX1, MX2, Message))

    def AsrtNeg__(self, MX, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_negative(%s, message=%s)" % (MX, Message))

    def AsrtNoNeg(self, MX, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_non_negative(%s, message=%s)" % (MX, Message))

    def AsrtPos__(self, MX, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_positive(%s, message=%s)" % (MX, Message))

    def AsrtNoPos(self, MX, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_non_positive(%s, message=%s)" % (MX, Message))

    def AsrtNoEq_(self, MX1, MX2, Message=None):
        if self.Switch4TFAssert:
            self.Statement("tf.debugging.assert_none_equal(%s, %s, message=%s)" % (MX1, MX2, Message))

    # Routines
    def RndIdxUni(self, VariableName, N_MoleculesToDistribute, Indices):
        self.Comment__("Select random values")
        self.RndNumUni("RndVals", Shape=N_MoleculesToDistribute, MaxVal="len(%s)" % Indices, Type='int32')
        self.Comment__("Select indices corresponding to the random values")
        self.Gather___(VariableName, Indices, "RndVals")

    def RndIdxWg_(self, DestVar, N_MoleculesToDistribute, Indices, Weights):
        self.Reshape__("Counts", Weights, -1)
        self.Reshape__("Indices", Indices, -1)
        self.VarRepeat("PoolOfIndices", "Indices", "Counts")
        self.Reshape__("N_IndicesToDraw", N_MoleculesToDistribute, -1)
        self.RndIdxUni(DestVar, "N_IndicesToDraw", "PoolOfIndices")

    def RndIdxWgN(self, DestVar, N_MoleculesToDistribute, Indices, Weights):
        self.Reshape__("Counts", Weights, -1)
        self.Reshape__("Indices", Indices, -1)

        # This code removes minimally represented indices from being selected
        self.ReduceMax("MaxCount", "Counts", 0)
        self.Variable_("NormalizationFactor", 1000)
        self.Divide___("NormalizationPoint", "MaxCount", "NormalizationFactor")
        self.Cast_____("NormalizationPoint", "NormalizationPoint", 'int32')
        self.FloorDiv_("Weights", "Counts", "NormalizationPoint")

        self.VarRepeat("PoolOfIndices", "Indices", "Weights")
        self.Reshape__("N_IndicesToDraw", N_MoleculesToDistribute, -1)
        self.RndIdxUni(DestVar, "N_IndicesToDraw", "PoolOfIndices")

    # Ultimate version
    def RndIdxWgh(self, DestVar, N_MoleculesToDistribute, Indices, Weights):
        self.Reshape__("Counts", Weights, -1)
        self.Reshape__("Indices", Indices, -1)

        # Shape[0] on Counts

        # Cumsum on Counts

        # Generate random number ranging between first and last elements of cumsum int values

        # Repeat random number array then reshape to generate matrix shaped with [Count shape, -1]

        # Boolean operation on (repeated random number array < Cumsum (2D matrix)) for which count range it falls
        
        # Boolean operation on (Counts > 0) for idxs with counts
        
        # Logical_and on (repeated random number array < Cumsum, Counts > 0) then convert to binary

        # Argmax on the binary to get the first indices of 1's in each random number generated.
        
        


        self.ReduceMax("MaxCount", "Counts", 0)
        self.Variable_("NormalizationFactor", 1000)
        self.Divide___("NormalizationPoint", "MaxCount", "NormalizationFactor")
        self.Cast_____("NormalizationPoint", "NormalizationPoint", 'int32')
        self.FloorDiv_("Weights", "Counts", "NormalizationPoint")

        self.VarRepeat("PoolOfIndices", "Indices", "Weights")
        self.Reshape__("N_IndicesToDraw", N_MoleculesToDistribute, -1)
        self.RndIdxUni(DestVar, "N_IndicesToDraw", "PoolOfIndices")

    def Put0InLen(self, LengthMatrix, CountMatrixToAdd):
        # Generate boolean and binary matrices for element availability in nascent proteins length matrix.
        self.Less_____("Bool_LenAvailable", LengthMatrix, 0)
        self.BoolToBin("Bin_LenAvailable", "Bool_LenAvailable")
    
        # Generate a cumulative sum matrix for labeling the available elements with count.
        self.CumSum___("LenCumsum", "Bin_LenAvailable", 1)
    
        # Generate a boolean matrix of the count label matrix > 0, < count.
        self.Greater__("Bool_LenCumsumGreaterThanZero", "LenCumsum", 0)
        self.Reshape__("CountMatrixToAdd", CountMatrixToAdd, [-1, 1])
        self.LessEq___("Bool_LenCumsumLessThanOrEqualToCountMatrixToAdd", "LenCumsum", "CountMatrixToAdd")
        self.LogicAnd_("Bool_LenCumsum", "Bool_LenCumsumGreaterThanZero", "Bool_LenCumsumLessThanOrEqualToCountMatrixToAdd")
    
        # Generate a binary matrix for elements satisfying both availability and count conditions.
        self.LogicAnd_("Bool_LenSelected", "Bool_LenAvailable", "Bool_LenCumsum")
        self.BoolToBin("Bin_LenSelected", "Bool_LenSelected")
    
        # Turn -1 to 0 where nascent proteins are to be made in the nascent protein length table.
        self.Add______(LengthMatrix, LengthMatrix, "Bin_LenSelected")
        self.BlankLine()

    def Count2Bin(self, DestVar, Count, N_Columns):
        # Generate binary pair vector
        self.Statement("Length_Count = tf.shape(%s)" % Count)
        self.Multiply_("Length_BinaryPairArray", "Length_Count", "2")
        self.VarTile__("BinaryPairArray", "[1, 0]", "Length_BinaryPairArray")

        # Generate binary count vector
        self.Reshape__("Count", Count, [1, -1])
        self.Subtract_("Count_Inv", N_Columns, "Count")
        self.Concat___("BinaryPairCount", "Count", "Count_Inv", 0)
        self.Transpose("BinaryPairCount", "BinaryPairCount")
        self.Reshape__("BinaryPairCount", "BinaryPairCount", -1)

        # Generate binary matrix
        self.VarRepeat("BinaryRepeatArray", "BinaryPairArray", "BinaryPairCount")
        self.Reshape__("BinaryMatrix", "BinaryRepeatArray", [-1, N_Columns])
        self.Statement("%s = BinaryMatrix" % DestVar)

    # # To operate in CellState Class?
    # def ID2IdxSim(self, MolIDList, MolID2Index):
    #     self.InitArrayWithZero("MolIndexList", 0)
    #     with self.Statement("for Molecule in %s:" % MolIDList):
    #         self.Statement("MolIndex = %s(Molecule)" % MolID2Index)
    #         self.Concat___("MolIndexList", 'MolIndex')
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

    def Conc2Cont(self, DestVar, MolIndexList, MolConcs):
        self.PrepGathr(MolConcs, MolIndexList)
        self.Gather___("MolConcs4MolIndexList", MolConcs, MolIndexList)
        self.Statement("%s = MolConcs4MolIndexList * (Cel.CellVol * Cst.NA)" % DestVar)

    def Cont2Conc(self, DestVar, MolIndexList, MolCounts):
        self.PrepGathr(MolCounts, MolIndexList)
        self.Gather___("MolCounts4MolIndexList", MolCounts, MolIndexList)
        self.Statement("%s = MolCounts4MolIndexList / (Cel.CellVol * Cst.NA)" % DestVar)

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
