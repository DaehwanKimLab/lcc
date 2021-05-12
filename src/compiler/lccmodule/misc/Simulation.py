

def Write_Simulation_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def Simulation_Init():"):
        Writer.Variable_("self.SimulationStep", 100) #
        Writer.Variable_("self.CellCycle", 1)