

def Write_Class_PromoterState_Init(Writer):
    Writer.BlankLine()
    with Writer.Statement("class FPromoterState(FGenomeState):"):
        with Writer.Function_("__init__", "Species", "CellID"):
            Writer.Statement("super().__init__(Species, CellID)")