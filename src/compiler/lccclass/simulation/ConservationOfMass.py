'''



'''

def Write_ConservationOfMass_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Function_("Init", ""):
        Writer.Pass_____()
        Writer.Statement()


def Write_ConservationOfMass_Loop(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Function_("Loop"):
        Writer.MatrixMul("Cel.Count_Master", "Cel.MW_Master")