

def Write_Class_GeneState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FGeneState(FGenomeState):"):
        with Writer.Statement("def __init__(self, Species, CellID):"):
            Writer.Statement("super().__init__(Species, CellID)")