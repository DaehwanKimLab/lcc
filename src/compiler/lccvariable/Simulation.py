
# Comp is a short hand for CompilerData
def Write_Simulation(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FSimulation():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Statement("# Define simulation parameters")
            Writer.Variable_("self.CellCycles", 0)
            Writer.Variable_("self.SimTimeRequested", 0)
            Writer.Variable_("self.SimTimeExecuted", 0)
            Writer.Variable_("self.SimStepsPerSecond", 0)  # Simulation time resolution
            Writer.Variable_("self.SimStepsRequested", 0)
            Writer.Variable_("self.SimStepsExecuted", 0)
            Writer.Statement("self.Initialize()")
            Writer.BlankLine()

        with Writer.Statement("def Initialize(self):"):
            # Implement user input for delta t value
            Writer.Variable_("self.SimTimeRequested", 100)
            Writer.Variable_("self.SimStepsPerSecond", 1)  # Simulation time resolution per second
            Writer.Statement("self.SimStepsRequested = self.SimTimeRequested / self.SimStepsPerSecond")


        with Writer.Statement("def RunSimulation(self, Cel, Cst, Env, Pro):"): # Pro being a cellular process
            Writer.PrintStrg("Simulation begins...")


            #
            # with Writer.Statement("for SimulationStep in range(Sim.SimStepsRequested):"):
            #     Writer.PrintStrg('=============================================')
            #     Writer.PrintStVa('SimulationStep: ', "SimulationStep + 1")
            #     Writer.BlankLine()

                # Run loop functions for each process
                # Writer.Statement("TE_Loop()")
                # Writer.Statement("TCS_Loop()")
                # Writer.Statement("RNADeg_Loop()")
                # Writer.Statement("PE_Loop()")
                # Writer.Statement("Metab_Loop()")
                # Writer.BlankLine()