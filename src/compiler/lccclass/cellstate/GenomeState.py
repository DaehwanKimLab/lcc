

def Write_Class_GenomeState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FGenomeState(FCellState):"):
        with Writer.Statement("def __init__(self, Species, CellID):"):
            Writer.Statement("super().__init__(Species, CellID)")