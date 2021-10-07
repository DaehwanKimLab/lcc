
# Comp is a short hand for CompilerData
def Write_Environment(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FEnvironment():"):
        with Writer.Function_("__init__"):
            Writer.Statement("pass")