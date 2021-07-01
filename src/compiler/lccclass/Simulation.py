
# Comp is a short hand for CompilerData
def Write_Simulation(Writer, Comp, ProGen):
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

        with Writer.Statement("def ConvertSimTimeToSimStep(self, Time, Step):"):
            Writer.Statement("self.SimStepsRequested = tf.math.multiply(self.SimTimeRequested, self.SimStepsPerSecond)")
            Writer.BlankLine()

        with Writer.Statement("def SetUpSimClock(self):"):
            # Implement user input for delta t value
            Writer.Variable_("self.SimTimeRequested", 100)  # Simulation time request in seconds
            Writer.Variable_("self.SimStepsPerSecond", 1)  # Simulation time resolution per second
            Writer.Statement("self.ConvertSimTimeToSimStep(self.SimTimeRequested, self.SimStepsRequested)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateSimTime(self):"):
            Writer.Statement("self.SimTimeExecuted = tf.math.divide(self.SimStepsExecuted, self.SimStepsPerSecond)")
            Writer.BlankLine()

        with Writer.Statement("def IncrementSimClock(self):"):
            Writer.AddElwise("self.SimStepsExecuted", 1)
            Writer.Statement("self.UpdateSimTime()")
            Writer.PrintVari("self.SimStepsExecuted")
            Writer.DebugPVar("self.SimTimeExecuted")
            Writer.BlankLine()

        with Writer.Statement("def SetUpStoichiometryMatrix(self, ProcessID):"):
            # Set up Cel.MX_Stoichiometries and Cel.MX_Rates
            Writer.Statement("ProcessID.SetUpStoichiometryMatrix()")
            Writer.BlankLine()

        with Writer.Statement("def SetUpCellStateMatrices(self):"):
            # Initialize reaction stoichiometry and rate matrices
            Writer.Statement("self.Cel.InitializeMatrices()")

            # Set up reaction stoichiometry matrix
            for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                Writer.Statement("self.SetUpStoichiometryMatrix(%s)" % ProcessID)

            # Finalize Reaction stoichiometry matrix
            Writer.Statement("self.Cel.FinalizeReactionMatrices()")

            # Load stoichiometry matrix to Exe
            Writer.Statement("self.LoadStoichiometryMatrix()")
            Writer.BlankLine()

        with Writer.Statement("def LoadStoichiometryMatrix(self):"):
            # Load updated stoichiometry to Exe
            Writer.Statement("self.Exe.LoadStoichMatrix(self.Cel.MX_Stoichiometries)")
            Writer.BlankLine()

        with Writer.Statement("def LoadRateMatrix(self):"):
            # Load updated stoichiometry to Exe
            Writer.Statement("self.Exe.LoadRateMatrix(self.Cel.MX_Rates)")
            Writer.BlankLine()

        with Writer.Statement("def LoadCountMatrix(self):"):
            Writer.Statement("self.Exe.LoadCountMatrix(self.Cel.MX_Counts)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateRates(self):"):  # Make the analogous method for stoich later
            # Reinitialize the rate matrix.
            Writer.Statement("self.Cel.ClearRateMatrix()")

            # Run all rate update codes
            for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                Writer.Statement("%s.UpdateRates()" % ProcessID)

            # Reshape the rate matrix for matrix multiplication
            Writer.Statement("self.Cel.FinalizeRateMatrix()")

            # Load updated rates to Exe
            Writer.Statement("self.LoadRateMatrix()")
            Writer.BlankLine()

        with Writer.Statement("def ExecuteReactions(self):"):
            Writer.Statement("self.Exe.MultiplyRXNMatrices()")
            Writer.Statement("self.Exe.RoundDeltaCountMatrix()")
            Writer.Statement("self.Exe.Solver()")
            Writer.Statement("self.Exe.OrdinaryDifferentialEquations()")
            Writer.BlankLine()

        with Writer.Statement("def UpdateCounts(self):"):
            Writer.Statement("self.LoadCountMatrix()")
            Writer.Overwrite("self.Cel.MX_Counts", "self.Exe.AddCountMatrices()")
            Writer.BlankLine()

        with Writer.Statement("def Initialize(self):"):
            Writer.Statement("self.SetUpSimClock()")
            Writer.Statement("self.SetUpCellStateMatrices()")
            Writer.BlankLine()

        with Writer.Statement("def Run(self):"):
            Writer.PrintStrg("Simulation Begins...")
            with Writer.Statement("while self.SimStepsExecuted < self.SimStepsRequested:"):
                # Increment simulation steps.
                Writer.Statement("self.IncrementSimClock()")

                # Models:
                # Rate Gauge Model (the current version)
                # This will be a user input later


                # Update rate matrix.
                Writer.Statement("self.UpdateRates()")


                # Execute reactions by multiplying reaction and rate matrices.
                # This calculates delta value for all molecular counts
                Writer.Statement("self.ExecuteReactions()")

                # Hook for plugging in cellular processes not based on matrix calculations


                # Overwrite delta value for the applicable molecular counts
                Writer.Statement("self.UpdateCounts()")


                Writer.BlankLine()


            Writer.PrintStrg("Simulation Ended")
            Writer.Statement('print("Total # of Simulation Steps: %s" % self.SimStepsExecuted)')
            Writer.Statement('print("Total # of Simulation Time: %s" % self.SimTimeExecuted)')
            Writer.BlankLine()
