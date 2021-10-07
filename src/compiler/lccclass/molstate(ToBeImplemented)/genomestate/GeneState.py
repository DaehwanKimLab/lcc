

def Write_Class_GeneState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FGeneState(FGenomeState):"):
        with Writer.Function_("__init__", "Species", "CellID"):
            Writer.Statement("super().__init__(Species, CellID)")