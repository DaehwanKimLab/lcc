

def Write_Class_ProteinMonomerState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FProteinMonomerState(FCellState):"):
        with Writer.Function_("__init__", "Species", "CellID"):
            Writer.Statement("super().__init__(Species, CellID)")