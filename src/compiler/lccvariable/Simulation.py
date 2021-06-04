

def Write_Simulation_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("class FSimulation():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Statement("# Define simulation parameters")
            Writer.Variable_("self.CellCycles", 1)
            Writer.Variable_("self.SimulationSteps", 1)
            Writer.Variable_("self.SimulationStepsPerSecond", 1) # Simulation time resolution

# def Write_Class_Simulation_