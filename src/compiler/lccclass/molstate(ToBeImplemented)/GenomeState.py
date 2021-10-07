

def Write_Class_GenomeState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FGenomeState(FCellState):"):
        with Writer.Function_("__init__", "Species", "CellID"):
            Writer.Statement("super().__init__(Species, CellID)")