import os
import datetime

# Comp is a short hand for CompilerData
def Write_Simulation(Writer, Comp, ProGen):

    # Save File Configuration
    if Writer.Switch4Save:
        SavePath = os.path.realpath(Comp.SavePath)
        ID = ''
        Label = ''
        Idx = ''
        Type = ''

        if Writer.Switch4SaveAllCounts:
            Label = 'All'
            ID = 'ID_Master.npy'
            Type = 'Type_Master.npy'

        elif Writer.Switch4SaveSpecificCounts:
            Label = 'Specific'  # Change according to the object to save
            ID = ''
            Idx = ''
            Type = ''

        Path_ID = os.path.join(SavePath, ID)
        Path_Idx = os.path.join(SavePath, Idx)
        Path_Type = os.path.join(SavePath, Type)

        DateTime = datetime.datetime.now().strftime('%Y%m%d-%H%M')
        SaveFileName = SavePath + '/DL_EcoliSimulation_%s_%sCounts' % (DateTime, Label)
        SaveFileName_Count = SaveFileName + '_Data.csv'
        SaveFileName_Process = SaveFileName + '_Process.csv'

        SaveFileName_Supplement = None
        # SaveFileName_Supplement = SaveFileName + '_Supplement.csv'

    Writer.BlankLine()

    with Writer.Statement("class FSimulation():"):
        with Writer.Function_("__init__", "Cel", "Cst", "Env", "Exe", "Dict_CellProcesses"):
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

        with Writer.Function_("ConvertSimTimeToSimStep", "Time", "Step"):
            Writer.Statement("self.SimStepsRequested = tf.math.multiply(self.SimWallTimeRequested, self.SimStepsPerSecond)")
            Writer.BlankLine()

        with Writer.Function_("StartSimTimer"):
            Writer.Statement("self.SimTimerStart = tf.timestamp()")
            Writer.BlankLine()

        with Writer.Function_("EndSimTimer"):
            Writer.Statement("self.SimTimerEnd = tf.timestamp()")
            Writer.Statement("self.SimRunTimeExecuted = self.SimTimerEnd - self.SimTimerStart")
            Writer.BlankLine()

        with Writer.Function_("CalculateSimTimeImprovement"):
            Writer.Statement("DeltaSimTime = self.SimWallTimeExecuted - self.SimRunTimeExecuted")
            Writer.Statement("self.SimTimePercentReduction = (DeltaSimTime / self.SimWallTimeExecuted) * 100")
            Writer.Statement("self.SimPerformanceImprovement = (DeltaSimTime / self.SimRunTimeExecuted) * 100")
            Writer.Statement("self.SimTimesSpeed = self.SimWallTimeExecuted / self.SimRunTimeExecuted")
            Writer.BlankLine()

        with Writer.Function_("ShowSimulationTime"):
            Writer.Statement("self.EndSimTimer()")
            Writer.Statement("self.CalculateSimTimeImprovement()")
            Writer.PrintLine()
            Writer.PrintStrg("Simulation Time Summary:")
            Writer.Statement('print("Simulation Steps Executed (s): %s" % self.SimStepsExecuted[0])')
            Writer.Statement('print("     Simulation Wall Time (s): %s" % self.SimWallTimeExecuted[0])')
            Writer.Statement('print("      Simulation Run Time (s): %s" % self.SimRunTimeExecuted)')
            Writer.Statement('print("           X times faster    : %s" % self.SimTimesSpeed[0])')
            Writer.BlankLine()

        with Writer.Function_("SetUpSimClock", "SimWallTimeRequested", "SimStepsPrintResolution"):
            # Implement user input for delta t value
            Writer.Statement("self.StartSimTimer()")
            Writer.Variable_("self.SimWallTimeRequested", "SimWallTimeRequested")  # Simulation time request in seconds
            Writer.Variable_("self.SimStepsPerSecond", 1)  # Simulation time resolution per second
            Writer.Statement("self.ConvertSimTimeToSimStep(self.SimWallTimeRequested, self.SimStepsRequested)")
            Writer.Variable_("self.SimStepsPrintResolution", "SimStepsPrintResolution")
            Writer.BlankLine()

        with Writer.Function_("UpdateSimTime"):
            Writer.Statement("self.SimWallTimeExecuted = tf.math.divide(self.SimStepsExecuted, self.SimStepsPerSecond)")
            Writer.BlankLine()

        with Writer.Function_("PrintSimStepsExecuted"):
            Writer.IfElse___("self.SIM_DisplayTrigger", True_Function="self.PrintSimStepsExecuted_Display", False_Function="self.Pass")
            Writer.BlankLine()

        with Writer.Function_("PrintSimStepsExecuted_Display"):
            Writer.PrintLine()
            Writer.PrintStrg("Simulation Steps Executed:")
            Writer.PrintVar_("self.SimStepsExecuted")
            Writer.DebugPVar("self.SimWallTimeExecuted")
            Writer.BlankLine()

        with Writer.Function_("IncrementSimClock"):
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

        with Writer.Function_("SIM_ClearDeltaCounts"):
            Writer.Statement("self.Cel.ClearDeltaCounts()")
            Writer.BlankLine()

        with Writer.Function_("SIM_DisplayTrigger"):
            Writer.ReturnVar("tf.math.floormod(self.SimStepsExecuted, self.SimStepsPrintResolution) == 0")
            Writer.BlankLine()

        with Writer.Function_("SIM_SetUpProcessSpecificVariables", "ProcessID"):
            # Set up Cel.MX_Stoichiometries and Cel.MX_Rates
            Writer.Statement("ProcessID.Init_ProcessSpecificVariables()")
            Writer.Statement("ProcessID.SetUp_ProcessSpecificVariables()")
            Writer.BlankLine()

        with Writer.Function_("SIM_SetUpCellStateMatrices"):
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

        with Writer.Function_("LoadStoichiometryMatrix"):
            # Load updated stoichiometry to Exe
            Writer.Statement("self.Exe.LoadStoichMatrix(self.Cel.MX_Stoichiometries)")
            Writer.BlankLine()

        with Writer.Function_("LoadRateMatrix"):
            # Load updated stoichiometry to Exe
            Writer.Statement("self.Exe.LoadRateMatrix(self.Cel.MX_Rates)")
            Writer.BlankLine()

        with Writer.Function_("LoadCountMatrix"):
            Writer.Statement("self.Exe.LoadCountMatrix(self.Cel.Counts)")
            Writer.BlankLine()

        with Writer.Function_("ExecuteReactions"):
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

        with Writer.Function_("SIM_ExecuteProcesses"):
            # Run all matrix operation codes
            for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                if Writer.Switch4ProcessTimer:
                    Writer.Statement("Timer = tf.timestamp()")
                Writer.Statement("self.%s.ExecuteProcess()" % ProcessID)
                if Writer.Switch4ProcessTimer:
                    Writer.PrintStVa("@@ Run Time for %s" % ProcessID, "tf.timestamp() - Timer")
                    Writer.BlankLine()

            Writer.BlankLine()

        with Writer.Function_("SIM_RunCellProcesses"):
            # Execute reactions by multiplying reaction and rate matrices.
            # This calculates delta value for all molecular counts
            # Writer.Statement("self.ExecuteReactions()")

            # Run all matrix operations: This level of mask may be removed after getting rid of "reaction" style classes
            Writer.Statement("self.SIM_ExecuteProcesses()")

            # Update delta Counts after executing cell processes.
            Writer.Statement("")

            Writer.BlankLine()

        with Writer.Function_("SIM_UpdateCounts"):
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

        with Writer.Function_("SIM_CheckCountsPos"):
            Writer.Statement("self.Cel.CheckCountsPos()")
            Writer.BlankLine()

        with Writer.Function_("SIM_ViewProcessDebuggingMessages"):
            Writer.IfElse___("self.SIM_DisplayTrigger", True_Function="self.SIM_ViewProcessDebuggingMessages_Display", False_Function="self.Pass")
            Writer.BlankLine()

        with Writer.Function_("SIM_ViewProcessDebuggingMessages_Display"):
            Writer.PrintStrg("# Simulation Step Process Debugging Message:")
            Writer.BlankLine()

            for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                Writer.PrintStrg("[%s]:" % ProcessID)
                Writer.Statement("self.%s.ViewProcessDebuggingMessages()" % ProcessID)
                Writer.BlankLine()

        with Writer.Function_("SIM_ViewProcessSummaries"):
            Writer.IfElse___("self.SIM_DisplayTrigger", True_Function="self.SIM_ViewProcessSummaries_Display", False_Function="self.Pass")
            Writer.BlankLine()

        with Writer.Function_("SIM_ViewProcessSummaries_Display"):
            Writer.PrintStrg("### Simulation Step Process Summary ###")
            for ProcessID, Module in ProGen.Dict_CellProcesses.items():
                    Writer.Statement("self.%s.ViewProcessSummary()" % ProcessID)
            if "CellDivision" in ProGen.Dict_CellProcesses:
                Writer.Statement("self.CellDivision.PrintMessage()")
            Writer.BlankLine()

        with Writer.Function_("ViewReplicationCompletion"):
            if "Replication" in ProGen.Dict_CellProcesses:
                Writer.PrintStVa("% Replication completion",
                                 "self.Replication.DeterminePercentReplicationCompletion()")
            else:
                Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Function_("ViewMacromoleculeCounts"):
            TotalCountQueriesUsingCelMasterIdx = ['Genes', 'RNAs', 'Proteins', 'Complexes']
            for Query in TotalCountQueriesUsingCelMasterIdx:
                Writer.Statement("Counts = self.Cel.GetCounts(self.Cel.Idx_Master_%s)" % Query)
                Writer.ReduceSum("TotalCount", "Counts")
                Writer.PrintStVa("Total # of %s" % Query, "TotalCount")
                Writer.BlankLine()

        with Writer.Function_("ViewBuildingBlockCounts"):
            TotalCountQueriesUsingCelIdx = ['dNTPs', 'NTPs', 'AAs']
            for Query in TotalCountQueriesUsingCelIdx:
                BuildingBlock_Str = ""
                for BuildingBlock in Comp.BuildingBlock.Name2Key_BuildingBlocks[Query]:
                    BuildingBlock_Str += BuildingBlock
                Writer.PrintStVa("Total # of %s (%s)" % (Query, BuildingBlock_Str),
                                 "self.Cel.GetCounts(self.Cel.Idx_%s)" % Query)
                Writer.BlankLine()

        with Writer.Function_("ViewEnergyMoleculeCounts"):
            SingleCountQueries = ['ATP', 'NADH', 'NADPH', 'FADH2']
            for Query in SingleCountQueries:
                Writer.PrintStVa("Total # of %s" % Query, "self.Cel.GetCounts(self.Cel.Idx_%s)[0]" % Query)
                Writer.BlankLine()

        # TODO: TotalMass
        with Writer.Function_("ViewCellMass"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Function_("SIM_ViewCellStateSummary"):
            Writer.IfElse___("self.SIM_DisplayTrigger", True_Function="self.SIM_ViewProcessDebuggingMessages_Display", False_Function="self.Pass")


        with Writer.Function_("SIM_ViewCellStateSummary"):
            Writer.PrintStrg("### Simulation Step Current Status ###")
            Writer.Statement("self.ViewReplicationCompletion()")
            Writer.Statement("self.ViewMacromoleculeCounts()")
            # Writer.Statement("self.ViewBuildingBlockCounts()")
            # Writer.Statement("self.ViewEnergyMoleculeCounts()")
            # Writer.Statement("self.ViewCellMass()")
            Writer.BlankLine()

        with Writer.Function_("SIM_CellDivision"):
            Writer.Statement("self.CellDivision.ExecuteCellDivision()")
            if Writer.Switch4TestCellDivision:
                Writer.Statement("self.CellDivision.TestCellDivision()")
            Writer.BlankLine()

        if Writer.Switch4Save:
            with Writer.Function_("GenerateSupplementaryInfoSaveFile"):
                if SaveFileName_Supplement:
                    with Writer.Statement("with open('%s', 'w', newline='') as SaveFile:" % SaveFileName_Supplement):
                        Writer.Statement("ID = np.load(r'%s')" % Path_ID)
                        Writer.Statement("Type = np.load(r'%s')" % Path_Type)
                        Writer.Statement("Idx = np.arange(len(ID))")
                        Writer.Statement("Writer_Save = csv.writer(SaveFile)")
                        Writer.Statement("Writer_Save.writerow(list(np.array(ID)))")
                        Writer.Statement("Writer_Save.writerow(list(np.array(Type)))")
                        Writer.Statement("Writer_Save.writerow(list(np.array(Idx)))")
                else:
                    Writer.Pass_____()
                Writer.BlankLine()

            with Writer.Function_("GenerateHeaderRowInSaveFile"):
                with Writer.Statement("with open('%s', 'w', newline='') as SaveFile:" % SaveFileName_Count):
                    Writer.Statement("Header = np.load(r'%s')" % Path_ID)
                    Writer.Statement("Writer_Save = csv.writer(SaveFile)")
                    Writer.Statement("Writer_Save.writerow(list(np.array(Header)))")
                Writer.BlankLine()

            with Writer.Function_("AppendListAsRowInSaveFile", "Matrix"):
                with Writer.Statement("with open('%s', 'a+', newline='') as SaveFile:" % SaveFileName_Count):
                    Writer.Statement("Writer_Save = csv.writer(SaveFile)")
                    Writer.Statement("Writer_Save.writerow(list(np.array(Matrix)))")
                Writer.BlankLine()

            with Writer.Function_("SetUpIdxSave"):
                Writer.Pass_____()
                Writer.BlankLine()

            with Writer.Function_("SIM_SaveCounts"):
                if Writer.Switch4SaveAllCounts:
                    Writer.Statement("Count_Save = self.Cel.Counts")
                    Writer.BlankLine()

                elif Writer.Switch4SaveSpecificCounts:
                    Writer.Statement("Idx_Save = self.SetUpIdxSave()")
                    Writer.Statement("Count_Save = self.Cel.GetCounts(Idx_Save)")
                    Writer.BlankLine()

                # Save code
                Writer.Statement("self.AppendListAsRowInSaveFile(Count_Save)")
                # Writer.PrintStrg("### Simulation Step Current Count Saved ###")
                Writer.BlankLine()

        # with Writer.Statement("def ReplenishMetabolites(self, FinalCount):"):
        #     Writer.ScatNdUpd("FinalCount", "FinalCount", "self.MetaboliteIdxs", "self.MetaboliteCountsInitial")
        #     Writer.Reshape__("FinalCount_Replenished", "FinalCount", [-1, 1])
        #     Writer.ReturnVar("FinalCount_Replenished")
        #     Writer.BlankLine()

        with Writer.Function_("PostSimulationStepCorrection"):
            if "Metabolism" in ProGen.Dict_CellProcesses:
                Writer.Statement("self.Metabolism.ReplenishMetabolites()")
            if Writer.Switch4ProcessSummary and Writer.Switch4PostSimulationStepCorrectionMessage:
                Writer.IfElse___("self.SIM_DisplayTrigger", True_Function="self.PostSimulationStepCorrection_Display", False_Function="self.Pass")
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Function_("PostSimulationStepCorrection_Display"):
            Writer.PrintStrg("===== Post simulation step correction =====")
            Writer.Statement("self.Metabolism.Message_ReplenishMetabolites()")
            Writer.BlankLine()

        with Writer.Function_("Initialize", "SimWallTimeRequested", "SimStepsPrintResolution"):
            Writer.PrintStrg("Simulation Initialization Begins...")
            Writer.Statement("self.SetUpSimClock(SimWallTimeRequested, SimStepsPrintResolution)")
            Writer.Statement("self.SIM_SetUpCellStateMatrices()")
            if Writer.Switch4Save:
                Writer.Statement("self.GenerateHeaderRowInSaveFile()")
                Writer.Statement("self.GenerateSupplementaryInfoSaveFile()")
            Writer.PrintStrg("Simulation Initialization Completed.")
            Writer.BlankLine()

            # Writer.Statement("self.SIM_ViewProcessSummaries()")
            # Writer.BlankLine()

        # TODO: Replace the python while loop to tf.while loop
        with Writer.Function_("Run"):
            Writer.PrintStrg("Simulation Run Begins...")
            Writer.WhileLoop("self.Run_Condition", "self.Run_Body")
            Writer.BlankLine()

            Writer.PrintLine()
            Writer.PrintStrg("Simulation Run Completed.")
            Writer.BlankLine()

            Writer.Statement("self.ShowSimulationTime()")
            Writer.BlankLine()

        with Writer.Function_("Run_Condition", "Placeholder"):
            Writer.Less_____("SimulationCondition", "self.SimStepsExecuted", "self.SimStepsRequested")
            Writer.ReturnVar("SimulationCondition")
            Writer.BlankLine()

        with Writer.Function_("Run_Body", "Placeholder"):
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

            if Writer.Switch4CellStateSummary:
                Writer.Statement("self.SIM_ViewCellStateSummary()")

            # Overwrite delta value for the applicable molecular counts
            Writer.Statement("self.SIM_UpdateCounts()")

            if 'CellDivision' in ProGen.Dict_CellProcesses:
                Writer.Statement("self.SIM_CellDivision()")
                Writer.BlankLine()

            # Writer.Statement("self.SIM_CellFusion()")
            # Writer.BlankLine()

            if Writer.Switch4Save:
                Writer.Statement("self.SIM_SaveCounts()")
            if Writer.Switch4Save and Writer.Switch4Save:
                Writer.IfElse___("self.SIM_DisplayTrigger", True_Function="self.SIM_SaveCounts_Display", False_Function="self.Pass")
                Writer.BlankLine()

            # Apply post-simulation step corrections for selected reaction models:
            if Writer.Switch4PostSimulationStepCorrection:
                Writer.Statement("self.PostSimulationStepCorrection()")

            Writer.BlankLine()

            Writer.ReturnVar("[tf.constant(0)]")
            Writer.BlankLine()

        with Writer.Function_("SIM_SaveCounts_Display"):
            Writer.PrintLine()
            Writer.PrintStrg("Simulation Data Saved.")
            Writer.BlankLine()

        with Writer.Function_("Pass"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Function_("Visualize", "RequestedData"):
            if Writer.Switch4Visualization2D:
                Writer.PrintLine()
                Writer.PrintStrg("Displaying Requested Simulation Data in 2D plots.")
                Writer.Statement("Visualization2D.VisualizeData('%s', SimDataToVisualize['Molecules'])" % SaveFileName_Count)
                Writer.Statement("Visualization2D.VisualizeData('%s', SimDataToVisualize['ProcessAttributes'])" % SaveFileName_Process)
                Writer.PrintLine()
            else:
                Writer.Pass_____()
            Writer.BlankLine()