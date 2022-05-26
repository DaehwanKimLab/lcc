# Copyright 2015 gRPC authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""The Python implementation of the GRPC helloworld.Greeter server."""

# To build the .proto files:
# python -m grpc_tools.protoc -I./SimulationServer/protos --python_out=. --grpc_python_out=. ./SimulationServer/protos/<protofilename>.proto
# 

## lcc
import random
import SimModule
import numpy as np
from datetime import datetime
import time

## gRPC
import asyncio
from concurrent import futures
import logging

import grpc

import sys
sys.path.insert(0, './protos')
import lccsimulation_pb2
import lccsimulation_pb2_grpc


"""
    todo: move Sim related class definitions out somewhere else
"""

Screen_Size = W_S, H_S = 1200, 800

# Colors
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
GRAY1 = (230, 230, 230)
GRAY2 = (210, 210, 210)
GRAY3 = (150, 150, 150)
GRAY4 = (100, 100, 100)
YELLOW_FAINT = (200, 200, 150)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
YELLOW = (255, 255, 0)
BLUE = (0, 0, 255)
MAGENTA = (255, 0, 255)
CYAN = (0, 255, 255)

# Global variables
NA = 6.0221409e+23
pi = np.pi
uM = 1e-6
nM = 1e-9
UnitTxt = ''

# Determine global unit
Unit = nM
if Unit == nM:
    UnitTxt = 'nM'
elif Unit == uM:
    UnitTxt = 'uM'

GlucoseName = 'L'
HomeostasisMolName = ['Am', ]

class FEnvironment:
    def __init__(self, InX=W_S/2, InY=H_S/2, InShape='circle', InRadius=W_S*0.5, InThickness=10):
        self.X = InX
        self.Y = InY
        self.Shape = InShape
        self.Radius = int(InRadius)
        self.Thickness = InThickness
        self.TransparentCircleArea = None

    def Draw(self, shape=None):
        if shape == 'circle':
            pygame.draw.circle(Screen, YELLOW_FAINT, (self.X, self.Y), self.Radius)
        elif shape == 'lining':
            # pygame.draw.circle(Screen, YELLOW_FAINT, (self.X, self.Y), self.Radius)
            pygame.draw.circle(Screen, GRAY4, (self.X, self.Y), self.Radius, self.Thickness)
        else:
            pygame.draw.rect(Screen, YELLOW_FAINT, ((0, 0), (W_S, H_S)))

    def DrawTransparentArea(self):
        self.TransparentCircleArea = pygame.Surface((self.Radius * 2, self.Radius * 2), pygame.SRCALPHA)
        self.TransparentCircleArea.fill((255, 255, 255, 255))
        pygame.draw.circle(self.TransparentCircleArea, (0, 0, 0, 0), (self.Radius, self.Radius), self.Radius)

class FOrganism:
    def __init__(self, SimM, InName, InSpecies):
        self.SimM = SimM
        self.Name = InName
        self.Species = InSpecies
        self.X = None
        self.Y = None
        self.Angle = None

        self.BodyLength = 20
        self.BodyThickness = 10
        self.FlagellaLength_Fold = 2
        self.FlagellaThickness = 3

        # Memory for display & debugging
        self.Ligand_Prev = 0
        self.Am = 0
        self.SimCount = 0

        # Trajectory
        self.TrajectorySwitch = True
        self.Trajectory = dict()
        self.TrajectoryColor = list()

    def Draw(self):
        if self.Species == 'Ecoli':
            Color = ''
            # pygame.draw.circle(Screen, Color, (self.X, self.Y), 5)
            dX =  np.cos(self.Angle) * -self.BodyLength
            dY = -np.sin(self.Angle) * -self.BodyLength
            X_BodyEnd = self.X + dX
            Y_BodyEnd = self.Y + dY
            X_TailEnd = self.X + self.FlagellaLength_Fold * dX
            Y_TailEnd = self.Y + self.FlagellaLength_Fold * dY
            for i in range(self.X.size):
                if i == 0:
                    Color = RED
                else:
                    Color = YELLOW
                pygame.draw.line(Screen, Color, (self.X[i], self.Y[i]), (X_BodyEnd[i], Y_BodyEnd[i]), self.BodyThickness)
                pygame.draw.line(Screen, Color, (self.X[i], self.Y[i]), (X_TailEnd[i], Y_TailEnd[i]), self.FlagellaThickness)

    def SetPosition(self, Position):
        self.X = Position[0]
        self.Y = Position[1]
        self.Angle = Position[2]
        # Threshold value (Position[3]) is not used

    def Receptivity(self, N_SimulationsToPass=50):
        for _ in range(N_SimulationsToPass):
            self.SimM.SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution()
        self.SimM.SimLoop_WithSpatialSimulation()

    def IncrementSimCount(self):
        self.SimCount += 1

    def Initialize(self):
        self.SetPosition(self.SimM.GetPositionXYAngleByName(self.Name))
        self.InitializeTrajectory()
        self.IncrementSimCount()

    def ReportStatus(self):
        # for debugging
        self.Am = self.SimM.GetCountByName(HomeostasisMolName[0])
        self.SimCount += 1
        Ligand_Now = self.SimM.GetCountFromDistributionByNameAndPos(GlucoseName, self.Name)[0]
        SimStep = self.SimM.GetSimStep()
        Delta = 0
        if Ligand_Now != 0:
            Delta = (Ligand_Now - self.Ligand_Prev[0]) / Ligand_Now * 100
        # print('SimStep {:06d} [Chemotaxis  {:06d}] Ligand:{:.6f} {} ({}{:.4f}%) Am:{:.6f} {} (X:{:.2f}, Y:{:.2f}, {:3.1f} degree)'.format
        #in+ '   (SimStep, self.SimCount, Ligand_Now / Unit / NA, UnitTxt, ('+' if Delta >= 0 else ''), Delta, self.Am / Unit / NA, UnitTxt , self.X, self.Y, self.Angle / pi * 180))

    def Homeostasis(self, MolName=[]):
        self.SimM.Homeostasis(MolName)   # Input 'Am', 'qAm' here
        self.Ligand_Prev = self.SimM.GetCountFromDistributionByNameAndPos(GlucoseName, self.Name)
        self.SetPosition(self.SimM.GetPositionXYAngleByName(self.Name))
        self.InitializeTrajectory()
        self.ReportStatus()

    def InitializeTrajectory(self):
        for i in range(self.X.size):
            self.Trajectory[i] = [(self.X[i], self.Y[i])]
            self.TrajectoryColor.append(tuple(np.random.randint(0, 255, 3)))
        self.TrajectoryColor[0] = MAGENTA

    def AddToTrajectory(self):
        for i in range(self.X.size):
            self.Trajectory[i].append((self.X[i], self.Y[i]))

    def DrawTrajectory(self):
        for i in range(len(self.Trajectory)):
            pygame.draw.aalines(Screen, self.TrajectoryColor[i], False, self.Trajectory[i])

class FMolecule:
    def __init__(self, SimM, InName, InX, InY): # change: need to pass SimM to create these objects, 
        self.SimM = SimM
        self.Name = InName
        self.X = InX
        self.Y = InY
        self.DiffusionFactor = 2
        self.SpaceFactor = 50
        self.Color = ''
        self.Pattern = ''
        self.MaxBrightness = None
        self.NormalizationType = ''

        # Heatmap Drawing
        self.ReductionFactor = 5
        self.Max = 0
        self.InitializeHeatmapMax()

        # Particle Drawing
        self.Particle_N = 300
        self.Particle_PerLayer = 3
        self.Particle_Radius = 2
        self.Particle_SpreadFactor = 1.2
        self.Particle_XY_Static = []
        self.InitializeStaticParticles()

    def InitializeHeatmapMax(self):
        self.Max = np.max(self.SimM.GetDistributionByName(self.Name))

    def InitializeStaticParticles(self):
        for i in range(int(self.Particle_N / self.Particle_PerLayer)):
            for j in range(self.Particle_PerLayer):
                X = self.X + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))
                Y = self.Y + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))
                self.Particle_XY_Static.append((X, Y))

    def Reposition(self):
        self.X = random.randint(W_S * 2 / 5, W_S * 3 / 5)
        self.Y = random.randint(H_S * 2 / 5, H_S * 3 / 5)
        self.Particle_XY_Static = []
        self.InitializeStaticParticles()

    def Move(self, dX, dY):
        New_Particle_XY_Static = []
        for (X, Y) in self.Particle_XY_Static:
            X += dX
            Y += dY
            New_Particle_XY_Static.append((X, Y))
        self.Particle_XY_Static = New_Particle_XY_Static

    def Draw(self, Data, threshold=0):
        if self.Pattern == 'spots':
            Max = np.max(Data)
            if Max == 0:
                return
            else:
                Data = SimF.Normalize_Linear(Data)
                CoordsToDraw = np.where(Data > threshold)
                for X, Y, Value in zip(CoordsToDraw[0], CoordsToDraw[1], Data[CoordsToDraw]):
                    intensity = math.floor(Value * self.MaxBrightness)
                    if intensity > self.MaxBrightness:
                        intensity = self.MaxBrightness
                    color = self.GetColor(intensity)
                    pygame.draw.circle(Screen, color, (X, Y), self.Particle_Radius)

        elif self.Pattern == 'heatmap':
            Max = np.max(Data)
            if Max == 0:
                return
            else:
                if self.NormalizationType == 'linear':
                    Data = SimF.Normalize_Linear(Data)
                if self.NormalizationType == 'log':
                    Data = SimF.Normalize_P1Log(Data)
                for x in range(0, Data.shape[0], self.ReductionFactor):
                    for y in range(0, Data.shape[1], self.ReductionFactor):
                        PercentMolLevel = Data[x][y]
                        if PercentMolLevel < threshold or PercentMolLevel == 0:
                            continue
                        intensity = math.floor(PercentMolLevel * self.MaxBrightness)
                        if intensity > self.MaxBrightness:
                            intensity = self.MaxBrightness
                        color = self.GetColor(intensity)
                        pygame.draw.rect(Screen, color, ((x, y), (self.ReductionFactor, self.ReductionFactor)))

        elif self.Pattern == 'particle':
            for XY in self.Particle_XY_Static:
                pygame.draw.circle(Screen, BLUE, XY, self.Particle_Radius)

        else:
            assert True, 'Unsupported molecule distribution pattern for drawing: %s' % self.Pattern

    def GetColor(self, Intensity):
        if self.Color == 'Yellow':
            return (self.MaxBrightness, self.MaxBrightness, self.MaxBrightness - Intensity)
        if self.Color == 'Blue':
            return (self.MaxBrightness - Intensity, self.MaxBrightness - Intensity, self.MaxBrightness)

    def SetColor(self, Color, MaxBrightness):
        self.Color = Color
        self.MaxBrightness = MaxBrightness

    def SetPattern(self, Pattern, NormalizationType, ReductionFactor):
        self.Pattern = Pattern
        self.NormalizationType = NormalizationType
        self.ReductionFactor = ReductionFactor

def RandomInt(maxInt):
    return np.random.randint(0, high=maxInt)

def GenRandomColor():
    return lccsimulation_pb2.Vector3(X=RandomInt(256),Y=RandomInt(256),Z=RandomInt(256))    

class LCCSimulation(lccsimulation_pb2_grpc.LCCSimulationServicer):
    def __init__(self):
        self.IsRunning = False
    
    def Initialize(self, request, context):
        print('initializing')
        # Load model
        self.State = SimModule.FState()
        self.Data = SimModule.FDataset()
        self.DataManager = SimModule.FDataManager()
        self.SimM = SimModule.FSimulation(self.State, self.Data, self.DataManager)

        # Initialize model
        self.SimM.Initialize()

        L = FMolecule(self.SimM, 'L', 800, 500)
        L.SetColor('Yellow', 200)
        L.SetPattern('heatmap', 'linear', 5,)

        self.E = FOrganism(self.SimM, 'E', 'Ecoli')
        self.E.Initialize()
        self.E.Receptivity(100)

        # Setup Static Objects
        InitVisObjects = []
        ZeroVec = lccsimulation_pb2.Vector3(X=0,Y=0,Z=0)
        UnitVec = lccsimulation_pb2.Vector3(X=1,Y=1,Z=1)
        White = lccsimulation_pb2.Vector3(X=255,Y=255,Z=255)
        Blue = lccsimulation_pb2.Vector3(X=143, Y=186, Z=255)
        Yellow = lccsimulation_pb2.Vector3(X=255, Y=255, Z=102)

        # Static Glucose
        for Pos in L.Particle_XY_Static:
            NewObj = lccsimulation_pb2.VisObjectData(
                ID=0, 
                ObjType=lccsimulation_pb2.VisObjectType.M_GLUCOSE,
                Position=lccsimulation_pb2.Vector3(X=Pos[0], Y=Pos[1], Z=0),
                Rotation=ZeroVec,
                Scale=lccsimulation_pb2.Vector3(X=10,Y=10,Z=10),
                Color=White)
            
            InitVisObjects.append(NewObj)
        
        # Static Petri Dish
        InitVisObjects.append(lccsimulation_pb2.VisObjectData(
                ID=0, 
                ObjType=lccsimulation_pb2.VisObjectType.M_PETRI_DISH,
                Position=lccsimulation_pb2.Vector3(X=100, Y=100, Z=0),
                Rotation=ZeroVec,
                Scale=lccsimulation_pb2.Vector3(X=500,Y=500,Z=1),
                Color=Blue))

        # Dynamic EColi
        for i in range(len(self.E.X)):
            PosVec = lccsimulation_pb2.Vector3(X=self.E.X[i], Y=self.E.Y[i], Z=0)
            RotVec = lccsimulation_pb2.Vector3(X=0, Y=self.E.Angle[i] * (180/pi), Z=0)
            
            InitVisObjects.append(lccsimulation_pb2.VisObjectData(
                ID=i + 1,
                ObjType=lccsimulation_pb2.VisObjectType.M_ECOLI,
                Position=PosVec,
                Rotation=RotVec,
                Scale = lccsimulation_pb2.Vector3(X=0.2,Y=0.2,Z=0.2),
                Color = GenRandomColor(),
            ))
        
        # TODO: Create DNA Init Data Here

        return lccsimulation_pb2.InitData(InitObjects=InitVisObjects)

    # TODO: stream run
    def Run(self, request, context):
        print('running')
        self.IsRunning = True

        def response_messages():

            SimUnitTime = 0.2
            MinElapsedTime = 0.2
            ElapsedTime = 0
            PrevTime = datetime.now()

            while(True):
                if not self.IsRunning:
                    break
                start = time.time()

                CurrTime = datetime.now()
                ElapsedTime += (CurrTime - PrevTime).total_seconds()
                PrevTime = CurrTime

                while ElapsedTime >= SimUnitTime:
                    self.E.Receptivity(20)
                    self.E.SetPosition(self.SimM.GetPositionXYAngleByName('E'))
                    self.E.IncrementSimCount()

                    ElapsedTime -= SimUnitTime

                # Record data SimUnitState
                VisObjects = {} # map from id --> VisObjectData

                for i in range(len(self.E.X)):
                    PosVec = lccsimulation_pb2.Vector3(X=self.E.X[i], Y=self.E.Y[i], Z=0)
                    RotVec = lccsimulation_pb2.Vector3(X=0, Y=self.E.Angle[i] * (180/pi), Z=0)
                    
                    ObjID = i + 1; # ID 0 is used for static objects
                    VisObjects[ObjID] = lccsimulation_pb2.VisObjectData(
                        ID=ObjID,
                        ObjType=lccsimulation_pb2.VisObjectType.M_ECOLI,
                        Position=PosVec,
                        Rotation=RotVec,
                        # Scale =, # Scale doesn't change, leave it out
                        # Color =, # Color doesn't change, leave it out
                    )

                CurState = lccsimulation_pb2.SimUnitState(Time=self.E.SimCount, Objects=VisObjects)

                response = lccsimulation_pb2.RunData(State=CurState, Info="Any arbitrary message")

                time.sleep(max(1./14 - (time.time() - start), 0))
                yield response

        return response_messages()

    def Pause(self, request, context):
        print('pausing')
        self.IsRunning = False
        return lccsimulation_pb2.ControlSimulationResponse()

    def Stop(self, request, context):
        print('stopping')
        return lccsimulation_pb2.ControlSimulationResponse()


def serve():
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
    lccsimulation_pb2_grpc.add_LCCSimulationServicer_to_server(LCCSimulation(), server)

    server.add_insecure_port('[::]:50051')
    server.start()
    server.wait_for_termination()


if __name__ == '__main__':
    logging.basicConfig()
    serve()
