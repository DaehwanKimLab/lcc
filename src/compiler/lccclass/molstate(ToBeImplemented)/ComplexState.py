

def Write_Class_ComplexState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FComplexState(FCellState):"):
        with Writer.Function_("__init__", "Species", "CellID"):
            Writer.Statement("super().__init__(Species, CellID)")