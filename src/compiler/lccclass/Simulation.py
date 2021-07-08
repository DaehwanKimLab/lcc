
# Comp is a short hand for CompilerData
def Write_Simulation(Writer, Comp, ProGen):
    Writer.BlankLine()
    with Writer.Statement("class FSimulation():"):
        with Writer.Statement("def __init__(self, Bch, Cel, Cst, Env, Exe, Pol, Dict_CellProcesses):"):
            Writer.Statement("# Define simulation parameters")
            Writer.Variable_("self.CellCycleRequested", 0)
            Writer.Variable_("self.CellCycleExecuted", 0)
            
            Writer.Variable_("self.SimWallTimeRequested", 0)
            Writer.Variable_("self.SimWallTimeExecuted", 0)

            Writer.Variable_("self.SimStepsPerSecond", 0)  # Simulation time resolution
            Writer.Variable_("self.SimStepsRequested", 0)
            Writer.Variable_("self.SimStepsExecuted", 0)
            
            Writer.Variable_("self.SimTimerStart", 0)
            Writer.Variable_("self.SimTimerEnd", 0)
            Writer.Variable_("self.SimRunTimeExecuted", 0)

            Writer.Variable_("self.SimTimePercentReduction", 0)
            Writer.Variable_("self.SimPerformanceImprovement", 0)
            Writer.Variable_("self.SimTimesSpeed", 0)

            Writer.Variable_("self.SimStepsPrintResolution", 0)

            Writer.BlankLine()

            Writer.LinkClObj('Cel')
            Writer.LinkClObj('Cst')
            Writer.LinkClObj('Env')

            Writer.LinkClObj('Bch')
            Writer.LinkClObj('Pol')
            Writer.LinkClObj('Exe')

            Writer.BlankLine()

            for ProcessID in ProGen.Dict_CellProcesses.keys():
                Writer.Statement("self.%s = Dict_CellProcesses['%s']" % (ProcessID, ProcessID))

            Writer.BlankLine()

        with Writer.Statement("def ConvertSimTimeToSimStep(self, Time, Step):"):
            Writer.Statement("self.SimStepsRequested = tf.math.multiply(self.SimWallTimeRequested, self.SimStepsPerSecond)")
            Writer.BlankLine()

        with Writer.Statement("def StartSimTimer(self):"):
            Writer.Statement("self.SimTimerStart = tf.timestamp()")
            Writer.BlankLine()

        with Writer.Statement("def EndSimTimer(self):"):
            Writer.Statement("self.SimTimerEnd = tf.timestamp()")
            Writer.Statement("self.SimRunTimeExecuted = self.SimTimerEnd - self.SimTimerStart")
            Writer.BlankLine()

        with Writer.Statement("def CalculateSimTimeImprovement(self):"):
            Writer.Statement("DeltaSimTime = self.SimWallTimeExecuted - self.SimRunTimeExecuted")
            Writer.Statement("self.SimTimePercentReduction = (DeltaSimTime / self.SimWallTimeExecuted) * 100")
            Writer.Statement("self.SimPerformanceImprovement = (DeltaSimTime / self.SimRunTimeExecuted) * 100")
            Writer.Statement("self.SimTimesSpeed = self.SimWallTimeExecuted / self.SimRunTimeExecuted")
            Writer.BlankLine()

        with Writer.Statement("def ShowSimulationTime(self):"):
            Writer.Statement("self.EndSimTimer()")
            Writer.Statement("self.CalculateSimTimeImprovement()")
            Writer.PrintLine()
            Writer.PrintStrg("Simulation Time Summary (seconds):")
            Writer.Statement('print("Simulation Steps Executed: %s" % self.SimStepsExecuted)')
            Writer.Statement('print("     Simulation Wall Time: %s" % self.SimWallTimeExecuted)')
            Writer.Statement('print("      Simulation Run Time: %s" % self.SimRunTimeExecuted)')
            Writer.Statement('print("           X times faster: %s" % self.SimTimesSpeed)')
            Writer.BlankLine()

        with Writer.Statement("def SetUpSimClock(self):"):
            # Implement user input for delta t value
            Writer.Statement("self.StartSimTimer()")
            Writer.Variable_("self.SimWallTimeRequested", 2400)  # Simulation time request in seconds
            Writer.Variable_("self.SimStepsPerSecond", 1)  # Simulation time resolution per second
            Writer.Statement("self.ConvertSimTimeToSimStep(self.SimWallTimeRequested, self.SimStepsRequested)")
            Writer.Variable_("self.SimStepsPrintResolution", 100)
            Writer.BlankLine()

        with Writer.Statement("def UpdateSimTime(self):"):
            Writer.Statement("self.SimWallTimeExecuted = tf.math.divide(self.SimStepsExecuted, self.SimStepsPerSecond)")
            Writer.BlankLine()

        with Writer.Statement("def PrintSimStepsExecuted(self):"):
            with Writer.Statement("if tf.math.floormod(self.SimStepsExecuted, self.SimStepsPrintResolution) == 0:"):
                Writer.PrintVari("self.SimStepsExecuted")
                Writer.DebugPVar("self.SimWallTimeExecuted")
            Writer.BlankLine()

        with Writer.Statement("def IncrementSimClock(self):"):
            Writer.AddElwise("self.SimStepsExecuted", 1)
            Writer.Statement("self.UpdateSimTime()")
            if Writer.Switch4SimStepsExecuted:
                Writer.Statement("self.PrintSimStepsExecuted()")
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
                Writer.Statement("self.SetUpStoichiometryMatrix(self.%s)" % ProcessID)

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
                Writer.Statement("self.%s.UpdateRates()" % ProcessID)

            # Reshape the rate matrix for matrix multiplication
            Writer.Statement("self.Cel.FinalizeRateMatrix()")

            # Load updated rates to Exe
            Writer.Statement("self.LoadRateMatrix()")
            Writer.BlankLine()

        with Writer.Statement("def ExecuteReactions(self):"):
            Writer.Statement("self.Exe.RunReactions()")
            Writer.BlankLine()

        with Writer.Statement("def UpdateCounts(self):"):
            Writer.Statement("self.LoadCountMatrix()")
            Writer.Overwrite("self.Cel.MX_Counts", "self.Exe.AddCountMatrices()")
            Writer.BlankLine()

        with Writer.Statement("def DisplayCounts(self):"):
            for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                Writer.Statement("self.%s.DisplayCounts()" % ProcessID)
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

                # Display the molecular counts of molecules involved in reactions
                if Writer.Switch4ProcessSummary:
                    Writer.Statement("self.DisplayCounts()")


                Writer.BlankLine()


            Writer.PrintStrg("Simulation Ended")

            Writer.BlankLine()

            Writer.Statement("self.ShowSimulationTime()")

            Writer.BlankLine()
