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

    def WriteBlankLine(self):
        print("", file=self.fp)
        return self