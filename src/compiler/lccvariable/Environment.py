

def Write_Environment_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("class FEnvironment():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Statement("pass")