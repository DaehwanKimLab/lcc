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

class CodeWriter:
    def __init__(self, CodeFile, IndentLevel=1):
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
        self.IndentationPrefix = '\t' * self.IndentLevel

    def IncreaseIndent(self):
        self.IndentLevel += 1
        self.BuildIndentationPrefix()

    def DecreaseIndent(self):
        self.IndentLevel -= 1
        if self.IndentLevel < 0:
            self.IndentLevel = 0

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

    # def WriteComment(self, Comment):
    #     Line = '# %s' % Comment
    #     self.WriteStatement(Line)
    #     return self

    def WriteBlankLine(self):
        self.WriteStatement("")
        return self

    def WritePrintStr(self, Str):
        Line = 'print("%s")' % Str
        self.WriteStatement(Line)
        return self

    def WritePrintVar(self, VariableName):
        Line = 'print("%s = ", %s)' % (VariableName, VariableName)
        self.WriteStatement(Line)
        return self

    def WritePrintStrVar(self, Str, VariableName):
        Line = 'print("%s: ", %s)' % (Str, VariableName)
        self.WriteStatement(Line)
        return self

    def WriteDebugStatement(self, Line):
        if self.DebugPrintSwitch or self.DebugAssertSwitch:
            self.WriteStatement(Line)
            return self

    def WriteDebugVariable(self, VariableName, Value):
        if self.DebugPrintSwitch or self.DebugAssertSwitch:
            Line = '%s = %s' % (VariableName, Value)
            self.WriteStatement(Line)
            return self

    def WriteDebugPrintStr(self, Str):
        if self.DebugPrintSwitch:
            Line = 'print("%s")' % Str
            self.WriteStatement(Line)
            return self

    def WriteDebugPrintVar(self, VariableName):
        if self.DebugPrintSwitch:
            Line = 'print("%s = ", %s)' % (VariableName, VariableName)
            self.WriteStatement(Line)
            return self

    def WriteDebugPrintStrVar(self, Str, VariableName):
        if self.DebugPrintSwitch:
            Line = 'print("%s: ", %s)' % (Str, VariableName)
            self.WriteStatement(Line)
            return self

    def WriteDebugAssert(self, Condition, ErrMsg):
         if self.DebugAssertSwitch:
             Line = "assert %s, '%s'" % (Condition, ErrMsg)
             self.WriteStatement(Line)
             return self