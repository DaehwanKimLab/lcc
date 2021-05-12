

def Write_Class_Environment_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FEnvironment():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Statement("")