

def Write_Class_RNAState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FRNAState(FCellState):"):
        with Writer.Statement("def __init__(self, Species, CellID):"):
            Writer.Statement("super().__init__(Species, CellID)")