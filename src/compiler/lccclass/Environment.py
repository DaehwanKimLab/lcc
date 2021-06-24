
# Comp is a short hand for CompilerData
def Write_Environment(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FEnvironment():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Statement("pass")