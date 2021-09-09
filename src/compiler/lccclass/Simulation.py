# Simulation time request in seconds
SimWallTimeRequested = 45 * 60
SimStepsPrintResolution = 50

# Comp is a short hand for CompilerData
def Write_Simulation(Writer, Comp, ProGen):
    Writer.BlankLine()
    with Writer.Statement("class FSimulation():"):
        with Writer.Statement("def __init__(self, Cel, Cst, Env, Exe, Dict_CellProcesses):"):
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

            Writer.LinkClObj('Exe')

            Writer.BlankLine()

            for ProcessID in ProGen.Dict_CellProcesses.keys():
                Writer.Statement("self.%s = Dict_CellProcesses['%s']" % (ProcessID, ProcessID))

            Writer.BlankLine()

        # Simulation time control code

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
            Writer.Variable_("self.SimWallTimeRequested", SimWallTimeRequested)  # Simulation time request in seconds
            Writer.Variable_("self.SimStepsPerSecond", 1)  # Simulation time resolution per second
            Writer.Statement("self.ConvertSimTimeToSimStep(self.SimWallTimeRequested, self.SimStepsRequested)")
            Writer.Variable_("self.SimStepsPrintResolution", SimStepsPrintResolution)
            Writer.BlankLine()

        with Writer.Statement("def UpdateSimTime(self):"):
            Writer.Statement("self.SimWallTimeExecuted = tf.math.divide(self.SimStepsExecuted, self.SimStepsPerSecond)")
            Writer.BlankLine()

        with Writer.Statement("def PrintSimStepsExecuted(self):"):
            with Writer.Statement("if tf.math.floormod(self.SimStepsExecuted, self.SimStepsPrintResolution) == 0:"):
                Writer.PrintLine()
                Writer.PrintStrg("Simulation Steps Executed:")
                Writer.PrintVar_("self.SimStepsExecuted")
                Writer.DebugPVar("self.SimWallTimeExecuted")
            Writer.BlankLine()

        with Writer.Statement("def IncrementSimClock(self):"):
            Writer.Add______("self.SimStepsExecuted", "self.SimStepsExecuted", 1)
            Writer.Statement("self.UpdateSimTime()")
            if Writer.Switch4SimStepsExecuted:
                Writer.Statement("self.PrintSimStepsExecuted()")
            Writer.BlankLine()

        # with Writer.Statement("def SetUpStoichiometryMatrix(self, ProcessID):"):
        #     # Set up Cel.MX_Stoichiometries and Cel.MX_Rates
        #     Writer.Statement("ProcessID.SetUpStoichiometryMatrix()")
        #     Writer.BlankLine()

        # Simulation module control code

        with Writer.Statement("def SIM_ClearDeltaCounts(self):"):
            Writer.Statement("self.Cel.ClearDeltaCounts()")
            Writer.BlankLine()

        with Writer.Statement("def SIM_SetUpProcessSpecificVariables(self, ProcessID):"):
            # Set up Cel.MX_Stoichiometries and Cel.MX_Rates
            Writer.Statement("ProcessID.Init_ProcessSpecificVariables()")
            Writer.Statement("ProcessID.SetUp_ProcessSpecificVariables()")
            Writer.BlankLine()

        with Writer.Statement("def SIM_SetUpCellStateMatrices(self):"):
            # # Set up reaction stoichiometry matrix
            # for ProcessID, Module in ProGen.Dict_CellProcesses.items():
            #     Writer.Statement("self.SetUpStoichiometryMatrix(self.%s)" % ProcessID)
            #
            # # Finalize Reaction stoichiometry matrix
            # Writer.Statement("self.Cel.FinalizeReactionMatrices()")
            #
            # # Load stoichiometry matrix to Exe
            # Writer.Statement("self.LoadStoichiometryMatrix()")
            # Writer.BlankLine()

            # Set up process-specific variables
            for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                Writer.Statement("self.SIM_SetUpProcessSpecificVariables(self.%s)" % ProcessID)
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
            Writer.Statement("self.Exe.LoadCountMatrix(self.Cel.Counts)")
            Writer.BlankLine()

        with Writer.Statement("def ExecuteReactions(self):"):
            # Update rate matrix.
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

            # Multiply reaction stoichiometry and rate matrix
            Writer.Statement("self.Exe.RunReactions()")
            Writer.BlankLine()

        with Writer.Statement("def SIM_ExecuteProcesses(self):"):
            # Run all matrix operation codes
            for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                Writer.Statement("self.%s.ExecuteProcess()" % ProcessID)
            Writer.BlankLine()

        with Writer.Statement("def SIM_RunCellProcesses(self):"):
            # Execute reactions by multiplying reaction and rate matrices.
            # This calculates delta value for all molecular counts
            # Writer.Statement("self.ExecuteReactions()")

            # Run all matrix operations: This level of mask may be removed after getting rid of "reaction" style classes
            Writer.Statement("self.SIM_ExecuteProcesses()")

            # Update delta Counts after executing cell processes.
            Writer.Statement("")

            Writer.BlankLine()

        with Writer.Statement("def SIM_UpdateCounts(self):"):
            if Writer.Switch4ShowDeltaCounts:
                Writer.Statement("self.Cel.ShowDeltaCounts()")

            if Writer.Switch4CheckDeltaCountsNeg:
                Writer.Statement("self.Cel.CheckDeltaCountsNeg()")


            # Update Counts

            Writer.Add______("self.Cel.Counts", "self.Cel.Counts", "self.Cel.DeltaCounts")

            if Writer.Switch4SoftCheckCounts or Writer.Switch4HardCheckCounts:
                Writer.Statement("self.SIM_CheckCountsPos()")

            # Writer.Statement("self.LoadCountMatrix()")
            # Writer.Overwrite("self.Cel.Counts", "self.Exe.AddCountMatrices()")
            Writer.BlankLine()

        with Writer.Statement("def SIM_CheckCountsPos(self):"):
            Writer.Statement("self.Cel.CheckCountsPos()")
            Writer.BlankLine()

        with Writer.Statement("def SIM_ViewProcessDebuggingMessages(self):"):
            with Writer.Statement("if tf.math.floormod(self.SimStepsExecuted, self.SimStepsPrintResolution) == 0:"):
                Writer.PrintStrg("# Simulation Step Process Debugging Message:")
                for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                        Writer.PrintStrg("[%s]:" % ProcessID)
                        Writer.Statement("self.%s.ViewProcessDebuggingMessages()" % ProcessID)
                Writer.BlankLine()

        with Writer.Statement("def SIM_ViewProcessSummaries(self):"):
            with Writer.Statement("if tf.math.floormod(self.SimStepsExecuted, self.SimStepsPrintResolution) == 0:"):
                Writer.PrintStrg("### Simulation Step Process Summary ###")
                for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                        Writer.Statement("self.%s.ViewProcessSummary()" % ProcessID)
                Writer.Statement("self.CellDivision.PrintMessage()")
                Writer.BlankLine()

        with Writer.Statement("def SIM_CellDivision(self):"):
            Writer.Statement("self.CellDivision.ExecuteCellDivision()")
            if Writer.Switch4TestCellDivision:
                Writer.Statement("self.CellDivision.TestCellDivision()")
            Writer.BlankLine()

        with Writer.Statement("def ReplenishMetabolites(self, FinalCount):"):
            Writer.ScatNdUpd("FinalCount", "FinalCount", "self.MetaboliteIdxs", "self.MetaboliteCountsInitial")
            Writer.Reshape__("FinalCount_Replenished", "FinalCount", [-1, 1])
            Writer.ReturnVar("FinalCount_Replenished")
            Writer.BlankLine()

        with Writer.Statement("def PostSimulationStepCorrection(self):"):
            Writer.Statement("self.Metabolism.ReplenishMetabolites()")
            with Writer.Statement("if tf.math.floormod(self.SimStepsExecuted, self.SimStepsPrintResolution) == 0:"):
                Writer.PrintStrg("===== Post simulation step correction =====")
                Writer.Statement("self.Metabolism.Message_ReplenishMetabolites()")
            Writer.BlankLine()

        with Writer.Statement("def Initialize(self):"):
            Writer.PrintStrg("Simulation Initialization Begins...")
            Writer.Statement("self.SetUpSimClock()")
            Writer.Statement("self.SIM_SetUpCellStateMatrices()")
            Writer.PrintStrg("Simulation Initialization Completed.")
            Writer.BlankLine()

            # Writer.Statement("self.SIM_ViewProcessSummaries()")
            # Writer.BlankLine()

        with Writer.Statement("def Run(self):"):
            Writer.PrintStrg("Simulation Run Begins...")
            with Writer.Statement("while self.SimStepsExecuted < self.SimStepsRequested:"):
                # Increment simulation steps.
                Writer.Statement("self.IncrementSimClock()")

                Writer.Statement("self.SIM_ClearDeltaCounts()")

                if Writer.Switch4ProcessDebuggingMessages:
                    Writer.PrintStrg(">>> Before running the processes <<<")
                    Writer.Statement("self.SIM_ViewProcessDebuggingMessages()")

                # Models:
                # Rate Gauge Model (the current version)
                # This will be a user input later

                # Hook for plugging in cellular processes not based on matrix calculations
                Writer.Statement("self.SIM_RunCellProcesses()")

                if Writer.Switch4ProcessDebuggingMessages:
                    Writer.PrintStrg(">>> After running the processes <<<")
                    Writer.Statement("self.SIM_ViewProcessDebuggingMessages()")

                # Display the molecular counts of molecules involved in reactions
                if Writer.Switch4ProcessSummary:
                    Writer.Statement("self.SIM_ViewProcessSummaries()")

                # Overwrite delta value for the applicable molecular counts
                Writer.Statement("self.SIM_UpdateCounts()")

                if 'CellDivision' in ProGen.Dict_CellProcesses:
                    Writer.Statement("self.SIM_CellDivision()")
                    Writer.BlankLine()

                # Writer.Statement("self.SIM_CellFusion()")
                # Writer.BlankLine()

                # Apply post-simulation step corrections for selected reaction models:
                if Writer.Switch4PostSimulationStepCorrection:
                    Writer.Statement("self.PostSimulationStepCorrection()")

                Writer.BlankLine()

            Writer.PrintLine()

            Writer.PrintStrg("Simulation Run Completed.")
            Writer.BlankLine()

            Writer.Statement("self.ShowSimulationTime()")
            Writer.BlankLine()
