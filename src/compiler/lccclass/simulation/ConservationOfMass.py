'''



'''

def Write_ConservationOfMass_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def Init():"):
        Writer.Pass_____()
        Writer.Statement()


def Write_ConservationOfMass_Loop(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def Loop(self):"):
        Writer.MatrixMul("Cel.Count_Master", "Cel.MW_Master")