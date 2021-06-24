
# Comp is a short hand for CompilerData
def Write_Simulation(Writer, Comp):
    Writer.BlankLine()
    with Writer.Statement("class FSimulation():"):
        with Writer.Statement("def __init__(self):"):
            Writer.Statement("# Define simulation parameters")
            Writer.Variable_("self.CellCycleRequested", 0)
            Writer.Variable_("self.CellCycleExecuted", 0)
            Writer.Variable_("self.SimTimeRequested", 0)
            Writer.Variable_("self.SimTimeExecuted", 0)
            Writer.Variable_("self.SimStepsPerSecond", 0)  # Simulation time resolution
            Writer.Variable_("self.SimStepsRequested", 0)
            Writer.Variable_("self.SimStepsExecuted", 0)
            Writer.BlankLine()

        #     Writer.Statement("self.Cel = None")
        #     Writer.Statement("self.Cst = None")
        #     Writer.Statement("self.Env = None")
        #     Writer.Statement("self.Exe = None")
        #     Writer.BlankLine()
        #
        # with Writer.Statement("def SetUpLink(self, Variable, Object):"):
        #     Writer.Overwrite("Variable", "Object")
        #     Writer.BlankLine()

        with Writer.Statement("def Initialize(self):"):
            # Implement user input for delta t value
            Writer.Variable_("self.SimTimeRequested", 100)  # Simulation time request in seconds
            Writer.Variable_("self.SimStepsPerSecond", 1)  # Simulation time resolution per second
            Writer.Statement("self.ConvertSimTimeToSimStep(self.SimTimeRequested, self.SimStepsRequested)")
            Writer.Statement("self.RunCellStateSetUp()")
            Writer.BlankLine()

        with Writer.Statement("def ConvertSimTimeToSimStep(self, Time, Step):"):
            Writer.Statement("self.SimStepsRequested = tf.math.multiply(self.SimTimeRequested, self.SimStepsPerSecond)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateSimTime(self):"):
            Writer.Statement("self.SimTimeExecuted = tf.math.divide(self.SimStepsExecuted, self.SimStepsPerSecond)")
            Writer.BlankLine()

        with Writer.Statement("def RunCellStateSetUp(self):"):
            # Generate reaction stoichiometry and rate matrices
            Writer.Statement("self.Cel.SetUpMatrices()")

            # Initialize all cell processes to set up Cel.Stoichs and Cel.Rates
            # with Writer.Statement("for Key_CellProcessStr, Value_CellProcessRef in Dict_CellProcess.items():"):
            #     Writer.Statement("Value_CellProcessRef.Init()")
            #     Writer.Pass_____()

            # Example
            Writer.Statement("SynDNA.InitProcess()")
            Writer.BlankLine()

            # Finalize RXN stoichiometry and rate matrices
            Writer.Statement("self.Cel.FinalizeMatrices()")

            # Load updated stoichiometry to Exe
            Writer.Statement("self.Exe.LoadStoichMatrix(self.Cel.Stoichs)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateRate(self):"):
            # Run all rate update codes
            # For each cell process


            # Reshape the rate matrix for matrix multiplication
            Writer.Statement("self.Cel.FinalizeMatrices()")

            # Load updated rates to Exe
            Writer.Statement("self.Exe.LoadRateMatrix(self.Cel.Rates)")
            Writer.BlankLine()

        with Writer.Statement("def RunReactions(self):"):
            Writer.Statement("self.Exe.MultiplyRXNMatrices()")
            Writer.Statement("self.Exe.RoundDeltaCountMatrix()")
            Writer.BlankLine()

        with Writer.Statement("def UpdateCounts(self):"):
            Writer.Statement("self.Exe.LoadCountMatrix(self.Cel.Counts)")
            Writer.Statement("self.Exe.AddCountMatrices()")
            Writer.Overwrite("self.Cel.Counts", "self.Exe.Count")
            Writer.BlankLine()

        with Writer.Statement("def Run(self):"):
            Writer.PrintStrg("Simulation Begins...")
            with Writer.Statement("while self.SimStepsExecuted < self.SimStepsRequested:"):
                # Increment simulation steps.
                Writer.AddElwise("self.SimStepsExecuted", 1)
                Writer.Statement("self.UpdateSimTime()")
                Writer.PrintVari("self.SimStepsExecuted")
                Writer.DebugPVar("self.SimTimeExecuted")

                # Models:
                # Rate Gauge Model (the current version)
                # This will be a user input later


                # Update rate matrix.
                Writer.Statement("self.UpdateRate()")


                # Execute reactions by multiplying reaction and rate matrices.
                # This calculates delta value for all molecular counts
                Writer.Statement("self.RunReactions()")

                # Hook for plugging in cellular processes not based on matrix calculations


                # Overwrite delta value for the applicable molecular counts
                Writer.Statement("self.UpdateCounts()")


                Writer.BlankLine()


            Writer.PrintStrg("Simulation Ended")
            Writer.Statement('print("Total # of Simulation Steps: %s" % self.SimStepsExecuted)')
            Writer.Statement('print("Total # of Simulation Time: %s" % self.SimTimeExecuted)')
            Writer.BlankLine()