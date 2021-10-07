

def Write_Class_LipidState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FLipidState(FCellState):"):
        with Writer.Function_("__init__", "Species", "CellID"):
            Writer.Statement("super().__init__(Species, CellID)")