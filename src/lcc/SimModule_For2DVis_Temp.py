import os, sys
import numpy as np
from datetime import datetime
import csv
import SimFunctions as sim
import plot

N_SimSteps = 30000
SimStepTimeResolution = 100
DisplayCount = 0   # 0: Display Counts only, 1: Display both Counts & dCounts
Unit = 'nM'   # nM or uM supported for now

class FState:
    def __init__(self):
        self.Vol = 0

        # State Arrays
        # Temporary legend attribute
        self.Mol_Names = list()
        self.Legends = list()

        self.Count_All = np.zeros([1, 9])
        self.Count_All = np.zeros([1, 9])
        self.dCount_All = np.zeros([1, 9])

        # Spatial Simulation
        self.Dist_Names = list()
        self.Dist_All = list()
        self.Idx_Dist_L = None

        self.Pos_Names = list()
        self.Idx_Pos_E = None
        self.Pos_X = None
        self.Pos_Y = None
        self.Pos_Angle = None
        self.Pos_Threshold = None

        # Standard_Unregulated
        self.Const_k_Reactant_Standard_Unregulated = None
        self.Const_k_Product_Standard_Unregulated = None
        self.Idx_Reactant_0_Standard_Unregulated = None
        self.Idx_Reactant_1_Standard_Unregulated = None
        self.Idx_Product_0_Standard_Unregulated = None
        self.Idx_Product_1_Standard_Unregulated = None


        self.Idx_Mol_InStoichMatrix_Standard_Unregulated = None
        self.Const_StoichMatrix_Standard_Unregulated = None

        # Enz_Standard_Unregulated
        self.Idx_Enz_Enz_Standard_Unregulated = None
        self.Const_k_Reactant_Enz_Standard_Unregulated = None
        self.Const_k_Product_Enz_Standard_Unregulated = None
        self.Idx_Reactant_0_Enz_Standard_Unregulated = None
        self.Idx_Reactant_1_Enz_Standard_Unregulated = None
        self.Idx_Product_0_Enz_Standard_Unregulated = None
        self.Idx_Product_1_Enz_Standard_Unregulated = None

        self.Idx_Mol_InStoichMatrix_Enz_Standard_Unregulated = None
        self.Const_StoichMatrix_Enz_Standard_Unregulated = None

        # Enz_MichaelisMenten_Unregulated
        self.Idx_Enz_Enz_MichaelisMenten_Unregulated = None
        self.Const_kcat_Enz_MichaelisMenten_Unregulated = None
        self.Const_KM_Enz_MichaelisMenten_Unregulated = None
        self.Idx_EnzSub_Enz_MichaelisMenten_Unregulated = None

        self.Idx_Mol_InStoichMatrix_Enz_MichaelisMenten_Unregulated = None
        self.Const_StoichMatrix_Enz_MichaelisMenten_Unregulated = None

    def Initialize(self):
        self.Vol = 1

        self.Idx_Dist_L = None

        self.Dist_Names = ['L', ]
        self.Idx_Dist_L = np.asmatrix([0])
        Dist = sim.InitializeDistribution(1200, 800, 800, 600, 6.02214e+16)
        self.Dist_All.append(Dist)

        self.Pos_Names = ['E', ]
        self.Idx_Pos_E = np.asmatrix([0])
        # Currently support X, Y, Angle, Threshold
        self.Pos_X = np.array([400, ])
        self.Pos_Y = np.array([400, ])
        self.Pos_Angle = np.array([90.0, ])
        self.Pos_Threshold = np.array([5.905186e+14, ])

        # Standard_Unregulated
        self.Const_k_Reactant_Standard_Unregulated = np.array([1.000000e+09, 1.000000e+09, ])
        self.Const_k_Product_Standard_Unregulated = np.array([1.000000e+00, 1.000000e+00, ])
        self.Idx_Reactant_0_Standard_Unregulated = np.asmatrix([2, 1, ])
        self.Idx_Reactant_1_Standard_Unregulated = np.asmatrix([8, 8, ])
        self.Idx_Product_0_Standard_Unregulated = np.asmatrix([4, 3, ])
        self.Idx_Product_1_Standard_Unregulated = np.asmatrix([0, 0, ])

        self.Idx_Mol_InStoichMatrix_Standard_Unregulated = np.asmatrix([2, 4, 8, 1, 3, ])
        self.Const_StoichMatrix_Standard_Unregulated = np.asmatrix([[-1, 1, -1, 0, 0, ], [0, 0, -1, -1, 1, ], ])

        # Enz_Standard_Unregulated
        self.Idx_Enz_Enz_Standard_Unregulated = np.asmatrix([2, ])
        self.Const_k_Reactant_Enz_Standard_Unregulated = np.array([5.000000e+07, ])
        self.Const_k_Product_Enz_Standard_Unregulated = np.array([5.000000e-03, ])
        self.Idx_Reactant_0_Enz_Standard_Unregulated = np.asmatrix([7, ])
        self.Idx_Reactant_1_Enz_Standard_Unregulated = np.asmatrix([0, ])
        self.Idx_Product_0_Enz_Standard_Unregulated = np.asmatrix([6, ])
        self.Idx_Product_1_Enz_Standard_Unregulated = np.asmatrix([0, ])

        self.Idx_Mol_InStoichMatrix_Enz_Standard_Unregulated = np.asmatrix([7, 6, ])
        self.Const_StoichMatrix_Enz_Standard_Unregulated = np.asmatrix([[-1, 1, ], ])

        # Enz_MichaelisMenten_Unregulated
        self.Idx_Enz_Enz_MichaelisMenten_Unregulated = np.asmatrix([5, 5, 6, 6, ])
        self.Const_kcat_Enz_MichaelisMenten_Unregulated = np.array([1.000000e+00, 1.000000e+00, 2.000000e+02, 1.000000e+00, ])
        self.Const_KM_Enz_MichaelisMenten_Unregulated = np.array([1.000000e-19, 1.000000e-19, 1.000000e-09, 1.000000e-09, ])
        self.Idx_EnzSub_Enz_MichaelisMenten_Unregulated = np.asmatrix([1, 3, 2, 4, ])

        self.Idx_Mol_InStoichMatrix_Enz_MichaelisMenten_Unregulated = np.asmatrix([1, 2, 3, 4, ])
        self.Const_StoichMatrix_Enz_MichaelisMenten_Unregulated = np.asmatrix([[-1, 1, 0, 0, ], [0, 0, -1, 1, ], [1, -1, 0, 0, ], [0, 0, 1, -1, ], ])

        self.Mol_Names = ['Pseudo', 'A', 'Am', 'AL', 'AmL', 'R', 'BP', 'B', 'L', ]
        self.Legends = ['SimStep', 'Vol', 'Pseudo', 'A', 'Am', 'AL', 'AmL', 'R', 'BP', 'B', 'L', ]

        Idx_Mol = np.asmatrix([0, 1, 2, 3, 4, 5, 6, 7, 8, ])
        Count_Mol = np.array([6.022141e+23, 3.011070e+17, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.011070e+15, 0.000000e+00, 6.022141e+13, 6.022141e+16, ])
        MolarityFactor_Mol = np.array([1.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00, ])
        MolarityFactor_Mol = np.where(MolarityFactor_Mol == 1, self.Vol, 1)
        Count_Mol *= MolarityFactor_Mol
        np.put_along_axis(self.Count_All, Idx_Mol, Count_Mol, axis=1)

    def GetMolNames(self):
        return self.Mol_Names

    def GetDistributionNames(self):
        return self.Dist_Names

    def GetPosNames(self):
        return self.Pos_Names

    def ExportLegend(self):
        return self.Legends

    def ExportData(self, Time):
        Data = np.asmatrix(np.zeros(2 + 9))
        Data[0, 0:1] = Time
        Data[0, 1:2] = self.Vol
        Data[0, 2:11] = self.Count_All
        return Data

class FDataset:
    def __init__(self):
        self.Legend = list()
        self.Data = None

    def PrintLegend(self):
        print(self.Legend)

    def PrintData(self):
        print(self.Data)

class FSimulation:
    def __init__(self, InState, InDataset, InDataManager):
        self.N_SimSteps = 0
        self.SimStep = 0
        self.SimTimeResolutionPerSecond = 0

        self.State = InState
        self.Dataset = InDataset
        self.DataManager = InDataManager

        self.Idx_DistToCoord_L = None
        self.Idx_Restore_Pseudo = None
        self.Idx_Count_Homeostasis_Am = None
        self.Idx_Pos_Homeostasis_Am = None
        self.Homeostasis_Prev_Am = None

        # Debugging
        self.Debug_Idx_Molecules = list()
        self.UnitTxt = ''
        self.Unit = 0

    def Initialize(self, InN_SimSteps=1000, InTimeResolution=100):
        print('Simulation Initialized...')
        self.N_SimSteps = np.asmatrix([InN_SimSteps])
        self.SimTimeResolutionPerSecond = InTimeResolution

        self.Idx_DistToCoord_L = np.asmatrix([8])
        self.Idx_Restore_Pseudo = np.asmatrix([0])

        self.Idx_Count_Homeostasis_Am = np.asmatrix([2])
        self.Idx_Pos_Homeostasis_Am = np.asmatrix([0])
        self.Homeostasis_Prev_Am = 0

        self.State.Initialize()

        # Legend Export
        self.Dataset.Legend = self.State.ExportLegend()
        self.DataManager.SetLegend(self.Dataset.Legend)

        # Data Export
        self.ExportData()

        # Debugging
        self.Debug_SetIdxMoleculesToTrack()
        self.Debug_SetUnit(Unit)

    def SimLoop_WithSpatialSimulation(self):
        self.IncrementSimStep()
        # Run Spatial Simulation
        self.SpatialSimulation()

        # Run Reactions
        self.NonSpatialSimulation()
        # Update Substrate Count
        self.UpdateCounts()

        # Update Spatially Distributed Molecules On Count
        self.DistributionToCount()

        # Restore Substrate Count for Sustained Substrate Influx
        self.RestoreMoleculeCount()

    def SimLoop_WithoutSpatialSimulation(self):
        self.IncrementSimStep()
        self.NonSpatialSimulation()
        self.UpdateCounts()
        self.RestoreMoleculeCount()

    def SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution(self):
        self.IncrementSimStep()
        self.NonSpatialSimulation()
        self.UpdateCounts()
        self.RestoreMoleculeCount()
        self.DistributionToCount()

    def Run(self, Spatial=0):
        if Spatial:
            self.Run_WithSpatialSimulation()
        else:
            self.Run_WithoutSpatialSimulation()

    def Run_WithSpatialSimulation(self):
        print('Simulation Run Begins...')

        while self.SimStep < self.N_SimSteps:
            self.SimLoop_WithSpatialSimulation()
            # Trigger Event on Substrate Count
            self.TriggerEventMoleculeCount()

            # Save and Export Data
            self.ExportData()

    def Run_WithoutSpatialSimulation(self):
        print('Simulation Run_WithoutSpatialSimulation Begins...')

        while self.SimStep < self.N_SimSteps:
            self.SimLoop_WithoutSpatialSimulation()
            # Trigger Event on Substrate Count
            self.TriggerEventMoleculeCount()

            # Save and Export Data
            self.ExportData()

        print('Simulation Run_WithoutSpatialSimulation Completed')

    def ExportData(self):
        self.Dataset.Data = self.State.ExportData(self.SimStep/self.SimTimeResolutionPerSecond)
        self.DataManager.Add(self.Dataset.Data)

    def ApplySimTimeResolution(self, Rate):
        return Rate / self.SimTimeResolutionPerSecond

    def DistributionToCount(self):
        Count = self.GetCountFromDistributionByNameAndPos('L', 'E')
        np.put_along_axis(self.State.Count_All, self.Idx_DistToCoord_L, Count, axis=1)

    def RestoreMoleculeCount(self):
        np.put_along_axis(self.State.Count_All, self.Idx_Restore_Pseudo, 6.022141e+23 * self.State.Vol, axis=1)

    def TriggerEventMoleculeCount(self):
        Time = self.SimStep / self.SimTimeResolutionPerSecond

    def Homeostasis(self):
        print('Running simulation to achieve homeostasis for : Am, ')

        bNotHomeostasis_Am = True

        while (bNotHomeostasis_Am):
            self.SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution()

            Homeostasis_Now_Am = self.GetCountByName('Am')
            if Homeostasis_Now_Am > 0 and abs(Homeostasis_Now_Am - self.Homeostasis_Prev_Am) / Homeostasis_Now_Am < 1.000000e-07:
                bNotHomeostasis_Am = False
                self.State.Pos_Threshold[self.Idx_Pos_Homeostasis_Am] = Homeostasis_Now_Am * 9.999900e-01
                print('Homeostasis achieved for : Am @ {:.010f}'.format(self.Debug_ApplyUnit(Homeostasis_Now_Am)), self.UnitTxt)
                print('Homeostasis threshold set for : Am @ {:.010f}'.format(self.Debug_ApplyUnit(Homeostasis_Now_Am)), self.UnitTxt)
            self.Homeostasis_Prev_Am = Homeostasis_Now_Am

            # Trigger Event on Substrate Count
            self.TriggerEventMoleculeCount()

            # Save and Export Data
            self.ExportData()

    # Spatial Simulation related routines
    def SpatialSimulation(self):
        self.SpatialDiffusion()
        self.SpatialLocation()

    def SpatialDiffusion(self):
        self.State.Dist_All[0] = sim.DiffuseDistribution(self.State.Dist_All[0])

    def SpatialLocation(self):
        HomeostasisMolecule = self.GetCount(self.Idx_Count_Homeostasis_Am)
        self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle = sim.BacterialChemotaxis(np.array(HomeostasisMolecule), self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle, self.State.Pos_Threshold)

    def NonSpatialSimulation(self):
        self.StandardReactions()

        self.EnzymaticReactions()

    # Biochemical Reaction related routines
    def StandardReaction_Standard_Unregulated(self):
        Conc_Reactant_0 = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Reactant_0_Standard_Unregulated), self.State.Vol)
        Conc_Reactant_1 = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Reactant_1_Standard_Unregulated), self.State.Vol)
        Conc_Product_0 = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Product_0_Standard_Unregulated), self.State.Vol)
        Conc_Product_1 = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Product_1_Standard_Unregulated), self.State.Vol)
        Rate_Reactant = sim.Eqn_Standard_Unregulated_2(Conc_Reactant_0, Conc_Reactant_1, self.State.Const_k_Reactant_Standard_Unregulated)
        Rate_Reactant = self.ApplySimTimeResolution(Rate_Reactant)
        Rate_Reactant = sim.CheckRateAndConc_2(Rate_Reactant, Conc_Reactant_0, Conc_Reactant_1)
        Rate_Product = sim.Eqn_Standard_Unregulated_2(Conc_Product_0, Conc_Product_1, self.State.Const_k_Product_Standard_Unregulated)
        Rate_Product = self.ApplySimTimeResolution(Rate_Product)
        Rate_Product = sim.CheckRateAndConc_2(Rate_Product, Conc_Product_0, Conc_Product_1)
        Rate = Rate_Reactant - Rate_Product
        dConc_Mol_InStoichMatrix = sim.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_Standard_Unregulated, Rate)
        dCount_Mol_InStoichMatrix = sim.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)
        self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_Standard_Unregulated, dCount_Mol_InStoichMatrix)

    def EnzymaticReaction_Enz_Standard_Unregulated(self):
        Conc_Enz = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Enz_Enz_Standard_Unregulated), self.State.Vol)
        Conc_Reactant_0 = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Reactant_0_Enz_Standard_Unregulated), self.State.Vol)
        Conc_Reactant_1 = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Reactant_1_Enz_Standard_Unregulated), self.State.Vol)
        Conc_Product_0 = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Product_0_Enz_Standard_Unregulated), self.State.Vol)
        Conc_Product_1 = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Product_1_Enz_Standard_Unregulated), self.State.Vol)
        Rate_Reactant = sim.Eqn_Enz_Standard_Unregulated_2(Conc_Enz, Conc_Reactant_0, Conc_Reactant_1, self.State.Const_k_Reactant_Enz_Standard_Unregulated)
        Rate_Reactant = self.ApplySimTimeResolution(Rate_Reactant)
        Rate_Reactant = sim.CheckRateAndConc_2(Rate_Reactant, Conc_Reactant_0, Conc_Reactant_1)
        Rate_Product = sim.Eqn_Enz_Standard_Unregulated_2(Conc_Enz, Conc_Product_0, Conc_Product_1, self.State.Const_k_Product_Enz_Standard_Unregulated)
        Rate_Product = self.ApplySimTimeResolution(Rate_Product)
        Rate_Product = sim.CheckRateAndConc_2(Rate_Product, Conc_Product_0, Conc_Product_1)
        Rate = Rate_Reactant - Rate_Product
        dConc_Mol_InStoichMatrix = sim.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_Enz_Standard_Unregulated, Rate)
        dCount_Mol_InStoichMatrix = sim.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)
        self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_Enz_Standard_Unregulated, dCount_Mol_InStoichMatrix)

    def EnzymaticReaction_Enz_MichaelisMenten_Unregulated(self):
        Conc_Enz = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_Enz_Enz_MichaelisMenten_Unregulated), self.State.Vol)
        Conc_EnzSub = sim.CountToConc(np.take(self.State.Count_All, self.State.Idx_EnzSub_Enz_MichaelisMenten_Unregulated), self.State.Vol)
        Rate = sim.Eqn_Enz_MichaelisMenten_Unregulated(Conc_Enz, Conc_EnzSub, self.State.Const_kcat_Enz_MichaelisMenten_Unregulated, self.State.Const_KM_Enz_MichaelisMenten_Unregulated)
        Rate = self.ApplySimTimeResolution(Rate)
        dConc_Mol_InStoichMatrix = sim.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_Enz_MichaelisMenten_Unregulated, Rate)
        dCount_Mol_InStoichMatrix = sim.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)
        self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_Enz_MichaelisMenten_Unregulated, dCount_Mol_InStoichMatrix)

    def StandardReactions(self):
        self.StandardReaction_Standard_Unregulated()

    def EnzymaticReactions(self):
        self.EnzymaticReaction_Enz_Standard_Unregulated()
        self.EnzymaticReaction_Enz_MichaelisMenten_Unregulated()

    # Useful routines
    def GetSimStep(self):
        return self.SimStep

    def GetSimTime(self):
        return self.SimStep / self.SimTimeResolutionPerSecond

    def IncrementSimStep(self):
        self.SimStep += 1

    def GetCount(self, Idx):
        return np.take(self.State.Count_All, Idx)

    def GetConcentration(self, Idx):
        return sim.CountToConc(np.take(self.State.Count_All, Idx))

    def GetDistribution(self, Idx):
        return self.State.Dist_All[Idx]

    def AddTodCount(self, Idx, Values):
        dCountToAdd = np.zeros_like(self.State.dCount_All)
        np.put_along_axis(dCountToAdd, Idx, Values, axis=1)
        dCount_All_New = self.State.dCount_All + dCountToAdd
        ZeroTest = dCount_All_New + self.State.Count_All
        self.State.dCount_All =  np.where(ZeroTest < 0, dCount_All_New - ZeroTest, dCount_All_New)

    def UpdateCounts(self):
        self.State.Count_All += self.State.dCount_All
        self.CleardCounts()
    def CleardCounts(self):
        self.State.dCount_All = np.zeros_like(self.State.dCount_All)
    # Temporary routines
    def GetMolIdx(self, Name):
        return self.State.GetMolNames().index(Name)

    def GetCountByName(self, Name):
        return self.GetCount(self.GetMolIdx(Name))

    def GetConcentrationByName(self, Name):
        return self.GetConcentration(self.GetMolIdx(Name))

    def GetPosIdx(self, Name):
        return self.State.GetPosNames().index(Name)

    def GetPositionXY(self, Idx):
        return np.take(self.State.Pos_X, Idx), np.take(self.State.Pos_Y, Idx)

    def GetPositionXYByName(self, Name):
        Idx = self.GetPosIdx(Name)
        return self.GetPositionXY(Idx)

    def GetPositionXYAngle(self, Idx):
        return np.take(self.State.Pos_X, Idx), np.take(self.State.Pos_Y, Idx), np.take(self.State.Pos_Angle, Idx)

    def GetPositionXYAngleByName(self, Name):
        Idx = self.GetPosIdx(Name)
        return self.GetPositionXYAngle(Idx)

    def GetDistIdx(self, Name):
        return self.State.GetDistributionNames().index(Name)

    def GetDistributionByName(self, Name):
        Idx = self.GetDistIdx(Name)
        return self.GetDistribution(Idx)

    def GetCountFromDistribution(self, Dist, X, Y):
        return Dist[int(X)][int(Y)]

    def GetCountFromDistributionByNameAndPos(self, NameOfDist, NameOfPos):
        X, Y = self.GetPositionXYByName(NameOfPos)
        Dist = self.GetDistributionByName(NameOfDist)
        return self.GetCountFromDistribution(Dist, X, Y)

    # Temporary routines
    def OverElongationCorrection(self, Len_Elongated, Max):   # Some polymerization process may not have max
        Len_Over = np.where(Len_Elongated > Max, Len_Elongated - Max, 0)
        return Len_Elongated - Len_Over

    def BuildingBlockConsumption(self, Freq, N_Elongated_PerSpecies):
        Raw = sim.DetermineAmountOfBuildingBlocks(Freq, N_Elongated_PerSpecies)
        Rounded = np.around(Raw)

        # Discrepancy handling
        N_Elongated = np.sum(N_Elongated_PerSpecies)
        Discrepancy = np.sum(Rounded) - N_Elongated
        NUniq_BuildingBlocks = Freq.shape[1]
        Sets, Remainder = np.divmod(Discrepancy, NUniq_BuildingBlocks)
        return Rounded + np.ones(NUniq_BuildingBlocks) * np.int32(Sets) + np.concatenate((np.ones(np.int32(np.round(Remainder))), np.zeros(np.int32(np.around(NUniq_BuildingBlocks - Remainder)))))

    # Polymerase Reaction related
    def Initiation(self, Len_Template, Len_Target, Idx_Pol, Idx_Template, Idx_TemplateSubset, Weight, PolThreshold):
        # Get available, active polymerase count - TO BE UPDATED with more regulatory algorithms
        Count_Pol = self.GetCount(Idx_Pol)
        Count_Pol_Active = np.floor_divide(Count_Pol, 2).astype(int)
        Count_Pol_Occupied = np.count_nonzero(np.where(Len_Target != -1, 1, 0)) * PolThreshold
        Count_Pol_Avail = Count_Pol_Active - Count_Pol_Occupied
        Count_Pol_FunctionalUnit = np.floor_divide(Count_Pol_Avail, PolThreshold)
        Count_Pol_Avail = np.where(Count_Pol_FunctionalUnit > 0, Count_Pol_FunctionalUnit, 0)[0, 0]

        # Get final initiation weight by applying initiation site count
        Count_Template_Complete = self.GetCount(Idx_Template)
        Count_Template_Nascent = np.count_nonzero(np.where(Len_Template != -1, 1, 0), axis=0)   # Assumption: each nascent template has one highly efficient initiation site
        Count_TemplateSubset_Nascent = np.take(Count_Template_Nascent, Idx_TemplateSubset)
        Count_InitiationSite = Count_Template_Complete + Count_TemplateSubset_Nascent
        Weight_Initiation = Count_InitiationSite * Weight

        # Get randomly selected target indices
        Idx_Selected = sim.PickRandomIdx(Count_Pol_Avail, Idx_Template, Weight_Initiation)
        Len_Target_Initiated = sim.InsertZeroIntoNegOneElementInLenMatrix(Len_Target, Idx_Selected)
        # Export Data
        # N_Initiated

        return Len_Target_Initiated

    def Elongation(self, Len, Max, Rate, Weight, Freq, Idx_PolSub, Idx_BB):
        NUniq_BuildingBlocks = Freq.shape[1]
        NUniq_Species = Freq.shape[0]

        dLength = self.ApplySimTimeResolution(Rate)   # this is not necessarily true based on the reaction input
        Len_Elongated = np.where(Len >= 0, Len + dLength, Len)

        Len_Trimmed = self.OverElongationCorrection(Len_Elongated, Max)
        N_Elongated_PerSpecies = np.asmatrix(np.sum(Len_Trimmed - Len, axis=0))   # This step loses shape for some reason, hence apply matrix again
        N_Elongated = np.sum(N_Elongated_PerSpecies)

        Consumed_BB = self.BuildingBlockConsumption(Freq, N_Elongated_PerSpecies)
        # Update dCount for BuildingBlocks
        self.AddTodCount(Idx_BB, -Consumed_BB)

        # Update dCount for Polymerase Reaction Substrates (To be updated by the reaction matrix form
        self.AddTodCount(Idx_PolSub, N_Elongated)

        # Export Data
        # N_Elongated
        return Len_Trimmed

    def Termination(self, Len, Max, Idx_Target):   # Some polymerization process may not have max
        Bool_Completed = (Len == Max)
        N_Completed_PerSpecies = np.sum(Bool_Completed, axis=0)
        N_Completed = np.sum(N_Completed_PerSpecies)
        Len_Completed = np.where(Bool_Completed, -1, Len)
        # Update dCount for BuildingBlocks
        self.AddTodCount(Idx_Target, N_Completed_PerSpecies)

        # Export Data
        # N_Completed
        return Len_Completed

    def Debug_SetIdxMoleculesToTrack(self):
        # Add a list of molecules to track for debugging every simulation step
        Debug_Names_Molecules = []

        if Debug_Names_Molecules:
            for Name in Debug_Names_Molecules:
                self.Debug_Idx_Molecules.append(self.State.GetMolNames().index(Name))
            self.Debug_AllCountsSwitch = False
        else:
            self.Debug_Idx_Molecules = list(range(len(self.State.GetMolNames())))
            self.Debug_AllCountsSwitch = True

    def Debug_PrintCounts(self, Switch):
        self.Debug_PrintSimStepTime()
        for Idx in self.Debug_Idx_Molecules:
            self.Debug_PrintCount(Idx)
        print()
        if Switch:
            self.Debug_PrintSimStepTime()
            for Idx in self.Debug_Idx_Molecules:
                self.Debug_PrintdCount(Idx)
            print()

    def Debug_PrintSimStepTime(self):
        Time = self.GetSimTime()
        print(self.SimStep, '(', round(Time,3), 's)', end='	| ')

    def Debug_PrintCount(self, Idx):
        print(' ' + self.State.GetMolNames()[Idx], end=': ')
        print('{:010e}'.format(self.Debug_ApplyUnit(self.State.Count_All[0][Idx])), self.UnitTxt, end=' | ')

    def Debug_PrintdCount(self, Idx):
        print('d' + self.State.GetMolNames()[Idx], end=': ')
        print('{:010e}'.format(self.Debug_ApplyUnit(self.State.dCount_All[0][Idx])), self.UnitTxt, end=' | ')

    def Debug_SetUnit(self, Input):
        self.UnitTxt = Input
        if Unit == 'nM':
            self.Unit = 1e-9 * 6.022141e+23
        elif Unit == 'uM':
            self.Unit = 1e-6 * 6.022141e+23
        else:
            self.Unit = 1

    def Debug_ApplyUnit(self, Value):
        return Value / self.Unit

class FDataManager:
    def __init__(self):
        self.Legend = list()
        self.DataBuffer = list()

    def SetLegend(self, InLegend):
        self.Legend = InLegend

    def Add(self, InData):
        self.DataBuffer.append(InData)

    def SaveToFile(self, InFileName):
        with open(InFileName, 'w', newline='', encoding='utf-8') as OutFile:
            TsvWriter = csv.writer(OutFile, delimiter='\t')
            if self.Legend:
                TsvWriter.writerow(self.Legend)
            for Row in self.DataBuffer:
                TsvWriter.writerow(np.array(Row).flatten().tolist())

def main():   # add verbose

    State = FState()
    Data = FDataset()
    DataManager = FDataManager()
    Simulation = FSimulation(State, Data, DataManager)

    Simulation.Initialize(N_SimSteps, SimStepTimeResolution)
    Simulation.Run(Spatial=0) # 0: WithoutSpatialSimulation, 1: WithSpatialSimulation

    DataManager.SaveToFile('SimOut.tsv')

if __name__ == '__main__':
    main()
    plot.main()

