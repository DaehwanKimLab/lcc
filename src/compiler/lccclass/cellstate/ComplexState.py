

def Write_Class_ComplexState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FComplexState(FCellState):"):
        with Writer.Statement("def __init__(self, Species, CellID):"):
            Writer.Statement("super().__init__(Species, CellID)")