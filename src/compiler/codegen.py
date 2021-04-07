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

    def Reshape(self, DestVar, SrcVar, Shape):
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
