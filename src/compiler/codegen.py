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

class CodeWriter():
    def __init__(self, CodeFile, IndentLevel=0):
        self.IndentLevel = IndentLevel
        self.fp = CodeFile
        self.BuildIndentationPrefix()

    def __enter__(self):
        self.IncreaseIndent()

    def __exit__(self, type, value, traceback):
        self.DecreaseIndent()

    def BuildIndentationPrefix(self):
        self.IndentationPrefix = '\t' * self.IndentLevel

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

    def WriteStatement(self, Line):
        print("{Indent}{Line}".format(Indent=self.IndentationPrefix, Line=Line), file=self.fp)
        return self

    def WriteVariable(self, VariableName, Value):
        Line = '%s = %s' % (VariableName, Value)
        self.WriteStatement(Line)

    def WriteBlankLine(self):
        self.WriteStatement("")


    """
    Array creation functions
    """

    def InitArrayWithOne(self, VariableName, Shape, Type='float32'):
        # var = tf.ones([4,3])
        Line = "# Not implemented"
        self.WriteStatement(Line)

    def InitArrayWithZero(self, VariableName, Shape, Type='float32'):
        # var = tf.zeros([4,3])
        Line = "# Not implemented"
        self.WriteStatement(Line)

    def InitArrayWithValue(self, VariableName, Value, Type='float32'):
        # var = tf.constant([3, 4, 5, 6, 7])
        Line = "# Not implemented"
        self.WriteStatement(Line)

    """
    Array manipulation functions
    """
    def Reshape(self, DestVar, SrcVar, Shape):
        # var = tf.reshape(var, [4, -1])
        Line = "# Not implemented"
        self.WriteStatement(Line)



"""
Tensorflow code generator
"""
class TFCodeWriter(CodeWriter):
    def InitArrayWithOne(self, VariableName, Shape, Type='float32'):
        Line = "{VariableName} = tf.ones({Shape}, dtype=tf.{Type})"\
            .format(VariableName=VariableName, Shape=str(Shape), Type=Type)
        self.WriteStatement(Line)

    def InitArrayWithZero(self, VariableName, Shape, Type='float32'):
        Line = '%s = tf.zeros(%s, dtype=tf.%s)' % (VariableName, str(Shape), Type)
        self.WriteStatement(Line)

    def InitArrayWithValue(self, VariableName, Value, Type='float32'):
        Line = '%s = tf.constant(%s, dtype=%s)' % (VariableName, str(Value), Type)
        self.WriteStatement(Line)

    def Reshape(self, DestVar, SrcVar, Shape):
        Line = '%s = tf.reshape(%s, %s)' % (DestVar, SrcVar, str(Shape))
        self.WriteStatement(Line)

"""
Numpy type as string
"""
class NumpyType:
    TypeToString = dict({'float32': '"float32"', 'int32': '"int32"'})

"""
Numpy code generator
"""
class NumpyCodeWriter(CodeWriter):
    def InitArrayWithOne(self, VariableName, Shape, Type='float32'):
        Line = "{VariableName} = np.ones({Shape}).astype({Type})"\
            .format(VariableName=VariableName, Shape=str(Shape), Type=NumpyType.TypeToString[Type])
        self.WriteStatement(Line)

    def InitArrayWithZero(self, VariableName, Shape, Type='float32'):
        Line = "%s = np.zeros(%s).astype(%s)" % (VariableName, str(Shape), NumpyType.TypeToString[Type])
        self.WriteStatement(Line)

    def InitArrayWithValue(self, VariableName, Value, Type='float32'):
        Line = "%s = np.array(%s).astype(%s)" % (VariableName, str(Value), Type)
        self.WriteStatement(Line)

    def Reshape(self, DestVar, SrcVar, Shape):
        Line = '%s = np.reshape(%s, %s)' % (DestVar, SrcVar, str(Shape))
        self.WriteStatement(Line)
