

def Write_Class_ProteinMonomerState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FProteinMonomerState(FCellState):"):
        with Writer.Statement("def __init__(self, Species, CellID):"):
            Writer.Statement("super().__init__(Species, CellID)")