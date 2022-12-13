#!/usr/bin/env python3
import sys
from datetime import datetime
import numpy as np
import math
import pygame

np.random.seed(1)

# Colors
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
GRAY1 = (230, 230, 230)
GRAY2 = (210, 210, 210)
GRAY3 = (150, 150, 150)
GRAY4 = (100, 100, 100)
RED = (255, 0, 0)
RED_DARK = (139, 0, 0)
GREEN = (0, 255, 0)
GREEN_DARK = (0, 100, 0)
YELLOW = (255, 255, 0)
YELLOW_FAINT = (200, 200, 150)
YELLOW_DARK = (150, 150, 50)
BLUE = (0, 0, 255)
BLUE_DARK = (0, 0, 139)
MAGENTA = (255, 0, 255)
CYAN = (0, 255, 255)
ORANGE = (252, 102, 0)

pi = np.pi

def DisplayText(Text, XY, position='center', rotate=0, color=BLACK, type='mass'):
    if type == 'mass':
        Text = Font_MassObject.render(Text, True, color)
    elif type == 'individual':
        Text = Font_IndividualObject.render(Text, True, color)
    elif type == 'bold':
        Text = Font_Sans.render(Text, True, color)
    else:
        print('[ERROR] DisplayTextFunction: unsupported type parameter %s' % type)
    Text = pygame.transform.rotate(Text, rotate)
    Text_Rect = Text.get_rect()
    X, Y = XY
    if position == 'center':
        Text_Rect.center = (X, Y)
    if position == 'midleft':
        Text_Rect.midleft = (X, Y)
    elif position == 'midright':
        Text_Rect.midright = (X, Y)
    elif position == 'midtop':
        Text_Rect.midtop = (X, Y)
    else:
        print('[ERROR] DisplayTextFunction: unsupported position parameter %s' % position)
    Screen.blit(Text, Text_Rect)

def DisplayIteration(Iter):
    Text = 'Simulation Iteration: ' + str(round(Iter))
    DisplayText(Text, Screen.get_rect().midtop, position='midtop', type='bold')

def GetXYFromEllipse(Shape, Angle, ratio_factor=1.0):
    Center_X, Center_Y, Width, Height = Shape
    return Center_X + Width / 2 * ratio_factor * np.cos(Angle), Center_Y + Height / 2 * ratio_factor * np.sin(Angle)

def GetXYFromCenter(Center, Angle, Radius, anisotropy=1.0):
    Center_X, Center_Y = Center
    return Center_X + Radius * anisotropy * np.cos(Angle), Center_Y + Radius * np.sin(Angle)

def CheckIfWithinEllipse(X, Y, XYWH, membranebias_lower=1.0, membranebias_upper=1.0):
    assert membranebias_lower <= membranebias_upper
    Center_X, Center_Y, Width, Height = XYWH
    Evaluate = (X - Center_X) ** 2 / (Width / 2) ** 2 + (Y - Center_Y) ** 2 / (Height / 2) ** 2
    if Evaluate < 1:
        if membranebias_lower != 1.0 or membranebias_upper != 1.0:
            if Evaluate >= membranebias_lower and Evaluate <= membranebias_upper:
                return True
            else:
                return False
        else:
            return True
    else:
        return False


# pygame
pygame.init()

Title = 'MinCDE and FtsZ Ring'
pygame.display.set_caption(Title)

# Fonts
Font_Sans = pygame.font.Font('freesansbold.ttf', 20)
Font_Init = pygame.font.SysFont('arial', 50)
Font_MassObject = pygame.font.SysFont('arial', 20)
Font_IndividualObject = pygame.font.SysFont('arial', 15)

# Screen
Screen_Size = W_S, H_S = 1500, 800
Screen = pygame.display.set_mode(Screen_Size)
MID_X = W_S / 2
MID_Y = H_S / 2

# E coli Membrane
ScaleFactor_Dimension = 2.2
W_Ecoli = int(2000 / ScaleFactor_Dimension)
H_Ecoli = int(1000 / ScaleFactor_Dimension)
NumAngles = 100

# X-section along the axis of division of E coli
Offset_DivAxis_X = - W_S * 7 / 40
Offset_DivAxis_Y = 0

# X-section along the plane of division
Offset_DivPlane_X = W_S * 13 / 40
Offset_DivPlane_Y = 0

# Molecules
NumMinD = 5000
NumCytoMinD = NumMinD

NumMinE = 4000
NumCytoMinE = NumMinE

NumFtsA = 802

# Cytosolic compartment (currently support odd numbers only to have one absolute center compartment)
CytoCompartment_Rows = 3
CytoCompartment_Columns = 3
CytoCompartment_Width = W_Ecoli / CytoCompartment_Rows
CytoCompartment_Height = H_Ecoli / CytoCompartment_Columns

# Molecules distributed in cytosolic compartments
NumMinC = 200
NumCytoMinC = np.zeros([CytoCompartment_Rows, CytoCompartment_Columns])
NumCytoMinC[CytoCompartment_Rows // 2, CytoCompartment_Columns // 2] = NumMinC

NumFtsZ = 7898
NumCytoFtsZ = np.zeros([CytoCompartment_Rows, CytoCompartment_Columns])
NumCytoFtsZ[CytoCompartment_Rows // 2, CytoCompartment_Columns // 2] = NumFtsZ

# Drawing
ScaleFactor_Quantity = 0.4
ScaleFactor_Refresh = 50

Radius_Anchor = 4
Radius_MinD = 4
Radius_MinE = 4
Radius_MinC = 4
Radius_FtsZ = 4
Radius_FtsA = 4

Color_Anchor = YELLOW
Color_MinD = RED
Color_MinE = BLUE
Color_MinC = GREEN
Color_FtsZ = MAGENTA
Color_FtsA = BLACK

# Determine the center XY for each cytosolic compartment
CytoCompartment_X = np.zeros([CytoCompartment_Rows, CytoCompartment_Columns])
CytoCompartment_Y = np.zeros([CytoCompartment_Rows, CytoCompartment_Columns])
for row in range(CytoCompartment_Rows):
    Center_X = np.round(np.floor(CytoCompartment_Rows / 2))
    Offset_X = W_Ecoli / CytoCompartment_Rows * (row - Center_X)
    for column in range(CytoCompartment_Columns):
        Center_Y = np.round(np.floor(CytoCompartment_Columns / 2))
        Offset_Y = H_Ecoli / CytoCompartment_Columns * (column - Center_Y)

        CytoCompartment_X[row, column] = MID_X + Offset_DivAxis_X + Offset_X
        CytoCompartment_Y[row, column] = MID_Y + Offset_DivAxis_Y + Offset_Y

class FMembrane:
    def __init__(self, X=MID_X, Y=MID_Y, thickness=5, linecolor=GRAY4, bodycolor=YELLOW_FAINT, Label=False, width=0, height=0):
        self.X = X
        self.Y = Y
        self.W = width
        self.H = height

        self.X_DivPlane = X + Offset_DivPlane_X
        self.Y_DivPlane = Y + Offset_DivPlane_Y
        self.XY_DivPlane = self.X_DivPlane, self.Y_DivPlane
        self.XYR_DivPlane = self.X_DivPlane, self.Y_DivPlane, height / 2

        self.X_DivAxis = X + Offset_DivAxis_X
        self.Y_DivAxis = Y + Offset_DivAxis_Y
        self.XY_DivAxis = self.X_DivAxis, self.Y_DivAxis
        self.XYWH_DivAxis = self.X_DivAxis, self.Y_DivAxis, width, height


        self.Label = Label
        self.Thickness = thickness
        self.LineColor = linecolor
        self.BodyColor = bodycolor

    def Draw(self, axis=True, plane=True):
        if axis:
            pygame.draw.ellipse(Screen, self.BodyColor, (self.X_DivAxis - self.W / 2, self.Y_DivAxis - self.H / 2, self.W, self.H))
            pygame.draw.ellipse(Screen, self.LineColor, (self.X_DivAxis - self.W / 2, self.Y_DivAxis - self.H / 2, self.W, self.H), self.Thickness)

        if plane:
            pygame.draw.circle(Screen, self.BodyColor, (self.X_DivPlane, self.Y_DivPlane), self.H / 2)
            pygame.draw.circle(Screen, self.LineColor, (self.X_DivPlane, self.Y_DivPlane), self.H / 2, self.Thickness)


class FCluster:
    def __init__(self, Index, Angle, Reference, XYWH_DivAxis_Membrane, color):
        '''
        memb is X, Y, W, H of the membrane
        '''

        self.Index = Index
        self.Angle = Angle
        self.Reference = Reference
        self.Membrane = XYWH_DivAxis_Membrane

        # for Drawing
        self.X_Anchor, self.Y_Anchor = GetXYFromEllipse(XYWH_DivAxis_Membrane, Angle, ratio_factor=1.0)
        self.X_Cluster, self.Y_Cluster = GetXYFromEllipse(XYWH_DivAxis_Membrane, Angle, ratio_factor=0.75)
        self.Color_Anchor = color

        self.SpatialIndex = self.GetSpatialIdx()

    def GetSpatialIdx(self):
        Idx_X = np.logical_and(self.X_Cluster < CytoCompartment_X + CytoCompartment_Width / 2, self.X_Cluster > CytoCompartment_X - CytoCompartment_Width / 2)
        Idx_Y = np.logical_and(self.Y_Cluster < CytoCompartment_Y + CytoCompartment_Height / 2, self.Y_Cluster > CytoCompartment_Y - CytoCompartment_Height / 2)
        Idx_XY = np.logical_and(Idx_X, Idx_Y)
        return np.where(Idx_XY)

    def Draw(self):
        pass

class FMinCluster(FCluster):
    def __init__(self, Index, Angle, Reference, SizeD, SizeE, SizeC, XYWH_DivAxis_Membrane, color=GREEN):
        '''
        memb is X, Y, W, H of the membrane
        '''

        self.SizeD = SizeD
        self.SizeE = SizeE
        self.SizeC = SizeC

        FCluster.__init__(self, Index, Angle, Reference, XYWH_DivAxis_Membrane, color)

    def MembraneEffectFunc(self, Rad):
        return Rad * Rad

    def KonD(self):
        MaxMembraneEffect = self.MembraneEffectFunc(0.5 * pi)
        NewAngle = self.Angle if self.Angle < pi else self.Angle - pi
        MembraneEffect = self.MembraneEffectFunc(NewAngle - 0.5 * pi) / MaxMembraneEffect * 1.1

        # DK - debugging purposes
        # MembraneEffect *= abs(np.random.normal(2, 0.5))

        AggregationEffect = math.log(max(2, self.SizeD), 2)
        return MembraneEffect * AggregationEffect

    def KoffD(self):
        return self.KonD() * (self.SizeE / max(1, self.SizeD)) * np.random.normal(1.5, 0.2)

    def KonE(self):
        return self.KonD() * 0.2

    def KoffE(self):
        return max(0, self.SizeE - self.SizeD)

    def KonC(self):
        return self.KonD() * 0.2

    def KoffC(self):
        return max(0, self.SizeC - self.SizeD)

    def PrintInfo(self):
        print("rad: {:.2f}, MinD: {:5d}, Kon: {:.2f}".format(self.Angle, self.SizeD, self.Kon()))

    def Draw(self, show_anchor=False):
        X_Membrane, Y_Membrane, W_Membrane, H_Membrane = self.Membrane

        # Membrane anchor
        if show_anchor:
            pygame.draw.circle(surface=Screen, color=self.Color_Anchor, center=(self.X_Anchor, self.Y_Anchor), radius=Radius_Anchor)

        def DisplayMins(Size, additional_scale_factor=1.0, color=None, radius=1):
            Size = int(Size * ScaleFactor_Quantity * additional_scale_factor)
            X_Min = np.random.uniform(-1, 1, Size) * W_Membrane / 4 + X_Membrane + W_Membrane / 4 * (-1 if self.X_Anchor < X_Membrane else 1)
            Y_Min = np.random.normal(Y_Membrane, H_Membrane / 6, Size)
            for i in range(Size):
                if CheckIfWithinEllipse(X_Min[i], Y_Min[i], self.Membrane, membranebias_lower=0.8):
                    pygame.draw.circle(surface=Screen, color=color, center=(X_Min[i], Y_Min[i]), radius=radius)

        if self.SizeD > 0:
            DisplayMins(self.SizeD, additional_scale_factor=2, color=Color_MinD, radius=Radius_MinD)
        if self.SizeE > 0:
            DisplayMins(self.SizeE, additional_scale_factor=2, color=Color_MinE, radius=Radius_MinE)
        if self.SizeC > 0:
            DisplayMins(self.SizeC, additional_scale_factor=5, color=Color_MinC, radius=Radius_MinC)

class FFtsCluster(FCluster):
    def __init__(self, Index, Angle, Reference, SizeA, SizeZ, XYWH_DivAxis_Membrane, XYR_DivPlane_Membrane, color=GREEN, treadmilling=False):
        '''
        memb is X, Y, W, H of the membrane
        '''
        FCluster.__init__(self, Index, Angle, Reference, XYWH_DivAxis_Membrane, color)

        self.SizeA = SizeA
        self.SizeZ = SizeZ

        # DivPlane-related
        self.Color_FtsZ_DivPlane = np.random.randint(0, 255, 3)
        self.XYR_Membrane_DivPlane = self.X_DivPlane_Center, self.Y_DivPlane_Center, self.Radius = XYR_DivPlane_Membrane
        self.X_Anchor_DivPlane, self.Y_Anchor_DivPlane = GetXYFromCenter(
            (self.X_DivPlane_Center, self.Y_DivPlane_Center),
            self.Angle, self.Radius * 0.98)

        if treadmilling:
            self.ArcMode = True
            self.Angle_Stop = self.Angle

        else:
            self.ArcMode = False
            self.Angle_FtsZs_DivPlane = list()
            self.X_FtsZs_DivPlane = list()
            self.Y_FtsZs_DivPlane = list()



    def KonA(self):
        return 0

    def KoffA(self):
        return 0

    def KonZ(self, LocalCytoMinC):
        return self.SizeA / max(1, LocalCytoMinC)

    def KoffZ(self, LocalMembMinC):
        return LocalMembMinC

    def GetLocalMembMinC(self, MinClusters):
        if MinClusters[self.Index].SizeD > 0:
            return MinClusters[self.Index].SizeC
        else:
            DistanceFactorVsPole = np.cos(self.Angle) ** 2
            MinRef = 0
            if self.Angle > pi / 2 and self.Angle <= pi * 3 / 2:
                MinRef = int(NumAngles / 2)
            return MinClusters[MinRef].SizeC * DistanceFactorVsPole

    def Draw(self, show_anchor=False):
        X_Membrane, Y_Membrane, W_Membrane, H_Membrane = self.Membrane

        if show_anchor:
            pygame.draw.circle(surface=Screen, color=Color_Anchor, center=(self.X_Anchor, self.Y_Anchor), radius=Radius_Anchor)

        def DisplayFtss(Size, additional_scale_factor=1.0, color=None, radius=1, membranebias_lower=0.0, membranebias_upper=1.0):
            Size = int(Size * ScaleFactor_Quantity * additional_scale_factor)
            X_Fts = np.random.uniform(-1, 1, Size) * W_Membrane / (NumAngles / 2) + self.X_Anchor
            Y_Fts = np.random.uniform(-1, 1, Size) * H_Membrane / (NumAngles / 2) + self.Y_Anchor
            for i in range(Size):
                if CheckIfWithinEllipse(X_Fts[i], Y_Fts[i], self.Membrane, membranebias_lower=membranebias_lower, membranebias_upper=membranebias_upper):
                    pygame.draw.circle(surface=Screen, color=color, center=(X_Fts[i], Y_Fts[i]), radius=radius)

        def DisplayFtsZ(Size, additional_scale_factor=1.0, color=None, radius=1, membranebias_lower=0.0, membranebias_upper=1.0):
            Size = int(Size * ScaleFactor_Quantity * additional_scale_factor)
            X_Fts = np.random.uniform(-1, 1, Size) * W_Membrane / (NumAngles / 2) + self.X_Anchor
            Y_Fts = np.random.uniform(-1, 1, Size) * H_Membrane / 2 + Y_Membrane
            for i in range(Size):
                if CheckIfWithinEllipse(X_Fts[i], Y_Fts[i], self.Membrane, membranebias_lower=0, membranebias_upper=membranebias_upper):
                    pygame.draw.circle(surface=Screen, color=color, center=(X_Fts[i], Y_Fts[i]), radius=radius)

        if self.SizeA > 0:
            DisplayFtss(self.SizeA, additional_scale_factor=2, color=Color_FtsA, radius=Radius_FtsA, membranebias_lower=0.95)

        if self.SizeZ > 0:
            DisplayFtsZ(self.SizeZ, additional_scale_factor=0.5, color=Color_FtsZ, radius=Radius_FtsZ, membranebias_lower=0.9, membranebias_upper=0.95)

    def Draw_DivPlane(self, show_anchor=False):
        def DisplayFtsA_DivPlane():
            pygame.draw.circle(surface=Screen, color=Color_FtsA, center=(self.X_Anchor_DivPlane, self.Y_Anchor_DivPlane),
                               radius=Radius_Anchor)
            # Label = "FtsA" + "#" + str(self.Index) + "~FtsZ_{" + str(self.N_Children) + "}"
            # DisplayString(Label, self.X_DivPlane + self.Radius * 2.5 * np.cos(self.Angle),
            #               self.Y_DivPlane + self.Radius * 2.5 * np.sin(self.Angle), color=BLACK, Type='individual')

        def DisplayFtsZ_DivPlane():
            for i in range(len(self.Angle_FtsZs_DivPlane)):
                pygame.draw.circle(surface=Screen, color=Color_FtsZ, center=(
                self.X_FtsZs_DivPlane[i], self.Y_FtsZs_DivPlane[i]), radius=Radius_FtsZ, width=1)
            # DisplayString(self.Name + str(i), self.X_Children[i], self.Y_Children[i], color=BLACK)

        def DisplayFtsZ_DivPlane_Arc():
            def GetRectFromCircle(XY, Radius, scale_factor=1.0):
                X, Y = XY
                ScaledRadius = Radius * scale_factor
                return X - ScaledRadius, Y - ScaledRadius, ScaledRadius * 2, ScaledRadius * 2

            Rect = GetRectFromCircle((self.X_DivPlane_Center, self.Y_DivPlane_Center), self.Radius, scale_factor=0.95)
            pygame.draw.arc(surface=Screen, color=self.Color_FtsZ_DivPlane, rect=Rect, start_angle=self.Angle,
                            stop_angle=self.Angle_Stop, width=20)


        if show_anchor:
            DisplayFtsA_DivPlane()

        if self.ArcMode:
            DisplayFtsZ_DivPlane_Arc()
        else:
            if len(self.Angle_FtsZs_DivPlane) > 0:
                DisplayFtsZ_DivPlane()


class FLegends():
    def __init__(self):
        self.Names = ['MinClusterPivot', 'MinD', 'MinE', 'MinC', 'FtsZ', 'FtsA']
        self.Sizes = [Radius_Anchor, Radius_MinD, Radius_MinE, Radius_MinC, Radius_FtsZ, Radius_FtsA]
        self.Colors = [Color_Anchor, Color_MinD, Color_MinE, Color_MinC, Color_FtsZ, Color_FtsA]

    def Draw(self, show_anchor=False):
        AlignmentForText = 30
        LineSpacing = 25

        # DivAxis
        X, Y = int(W_S * 0.55), int(H_S * 0.8)
        for i in range(len(self.Names)):
            if not show_anchor and i == 0:
                continue
            pygame.draw.circle(Screen, self.Colors[i], (X, Y), self.Sizes[i])
            DisplayText(self.Names[i], (X + AlignmentForText, Y), position='midleft', color=self.Colors[i], type='mass')
            Y += LineSpacing

        # DivPlane
        X, Y = int(W_S * 0.75), int(H_S * 0.85)
        DisplayText('FtsZ polymers in different color', (X + AlignmentForText, Y), position='midleft', color=BLACK, type='mass')
        Y += LineSpacing

def Simulate():
    # Initialize Membrane
    Membrane = FMembrane(X=MID_X, Y=MID_Y, width=W_Ecoli, height=H_Ecoli, Label=True)

    # Initialize Molecules
    DeltaAngle = 2 * pi / NumAngles
    MinClusters = [FMinCluster(i, DeltaAngle * i, -1, 0, 0, 0, Membrane.XYWH_DivAxis, color=Color_Anchor) for i in range(NumAngles)]
    FtsClusters = [FFtsCluster(i, DeltaAngle * i, -1, int(NumFtsA / NumAngles), 0, Membrane.XYWH_DivAxis, Membrane.XYR_DivPlane, color=Color_FtsA, treadmilling=True) for i in range(NumAngles)]
    
    # Legends
    Legends = FLegends()

    def AddMinDs():
        global NumCytoMinD

        KonSum = 0
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Kon = MinCluster.KonD()
            KonSum += Kon

        if KonSum <= NumCytoMinD:
            ScaleFactor = 1
        else:
            ScaleFactor = NumCytoMinD / KonSum

        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Kon = MinCluster.KonD()
            Kon *= ScaleFactor
            ToAdd = min(int(Kon), NumCytoMinD)
            if ToAdd <= 0:
                continue

            MinCluster.SizeD += ToAdd
            NumCytoMinD -= ToAdd

    def RemoveMinDs():
        global NumCytoMinD
        global NumCytoMinE

        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Koff = MinCluster.KoffD()
            ToRemove = min(int(Koff), MinCluster.SizeD)
            if ToRemove <= 0:
                continue
            
            MinCluster.SizeD -= ToRemove
            NumCytoMinD += ToRemove

            if MinCluster.SizeD <= 0:
                assert MinCluster.SizeD == 0
                MinCluster.Reference = -1
                NumCytoMinE += MinCluster.SizeE
                MinCluster.SizeE = 0

        for i, MinCluster in enumerate(MinClusters):
            Reference = MinCluster.Reference
            while Reference >= 0:
                Cluster2 = MinClusters[Reference]
                assert Reference != Cluster2.Reference
                Reference = Cluster2.Reference
                if Reference < 0 and Cluster2.SizeD == 0:
                    MinCluster.Reference = -1

    def AddMinEs():
        global NumCytoMinE

        KonSum = 0
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Kon = MinCluster.KonD()
            KonSum += Kon

        if KonSum <= NumCytoMinE:
            ScaleFactor = 1
        else:
            ScaleFactor = NumCytoMinE / KonSum

        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Kon = MinCluster.KonE()
            Kon *= ScaleFactor
            ToAdd = min(int(Kon), NumCytoMinE)
            if ToAdd <= 0:
                continue

            MinCluster.SizeE += ToAdd
            NumCytoMinE -= ToAdd

    def RemoveMinEs():
        global NumCytoMinE
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Koff = MinCluster.KoffE()
            ToRemove = min(int(Koff), MinCluster.SizeE)
            if ToRemove <= 0:
                continue

            MinCluster.SizeE -= ToRemove
            NumCytoMinE += ToRemove

    def AddMinCs():
        global NumCytoMinC

        KonArray = np.zeros(len(MinClusters))

        KonSum = 0
        for i, MinCluster in enumerate(MinClusters):
            if MinCluster.Reference >= 0:
                continue

            Kon = MinCluster.KonC()
            KonArray[i] = Kon
            KonSum += Kon

        if KonSum <= np.sum(NumCytoMinC):
            KonArray *= 1
        else:
            KonArray *= np.sum(NumCytoMinC) / KonSum

        for i, MinCluster in enumerate(MinClusters):
            if MinCluster.Reference >= 0:
                continue

            ToAdd = min(int(KonArray[i]), int(NumCytoMinC[MinCluster.SpatialIndex][0]))
            if ToAdd <= 0:
                continue

            MinCluster.SizeC += ToAdd
            NumCytoMinC[MinCluster.SpatialIndex] -= ToAdd

    def RemoveMinCs():
        global NumCytoMinC
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Koff = MinCluster.KoffC()
            ToRemove = min(int(Koff), MinCluster.SizeC)
            if ToRemove <= 0:
                continue

            MinCluster.SizeC -= ToRemove
            NumCytoMinC[MinCluster.SpatialIndex] += ToRemove

    def AddFtsZs():
        global NumCytoFtsZ

        KonArray = np.zeros(len(FtsClusters))

        # Scaling
        KonSum = 0
        for i, FtsCluster in enumerate(FtsClusters):
            if FtsCluster.Reference >= 0:
                continue

            Kon = FtsCluster.KonZ(NumCytoMinC[FtsCluster.SpatialIndex][0])
            KonArray[i] = Kon
            KonSum += Kon

        if KonSum <= np.sum(NumCytoFtsZ):
            KonArray *= 1
        else:
            KonArray *= (np.sum(NumCytoFtsZ) / KonSum)

        # Add
        TotalAdd = 0
        for i, FtsCluster in enumerate(FtsClusters):
            if FtsCluster.Reference >= 0:
                continue

            ToAdd = min(int(KonArray[i]), int(NumCytoFtsZ[FtsCluster.SpatialIndex][0]))
            if ToAdd <= 0:
                continue

            FtsCluster.SizeZ += ToAdd
            TotalAdd += ToAdd
            NumCytoFtsZ[FtsCluster.SpatialIndex] -= ToAdd

        # For DivPlane
        Random = np.random.randint(0, NumAngles, TotalAdd)

        if FtsClusters[0].ArcMode:
            for j in range(TotalAdd):
                Idx = Random[j]
                FtsClusters[Idx].Angle_Stop += pi / 10 / FtsClusters[Idx].Radius

        else:
            for j in range(TotalAdd):
                Idx = Random[j]
                Angle = FtsClusters[Idx].Angle + len(FtsClusters[Idx].Angle_FtsZs_DivPlane) * pi / 180
                FtsClusters[Idx].Angle_FtsZs_DivPlane.append(Angle)
                X, Y = GetXYFromCenter(Membrane.XY_DivPlane, Angle, FtsClusters[Idx].Radius, RatioFactor=0.9)
                FtsClusters[Idx].X_FtsZs_DivPlane.append(X)
                FtsClusters[Idx].Y_FtsZs_DivPlane.append(Y)

            # DL: Debugging
            # if TotalAdd > 0:
            #     for FtsCluster in FtsClusters:
            #         if len(FtsCluster.X_FtsZs_DivPlane):
            #             print(FtsCluster.Index, len(FtsCluster.Angle_FtsZs_DivPlane), FtsCluster.Angle_FtsZs_DivPlane)

    def RemoveFtsZs():
        global NumCytoFtsZ

        TotalRemove = 0
        for i, FtsCluster in enumerate(FtsClusters):
            if FtsCluster.Reference >= 0:
                continue

            LocalMembMinC = FtsCluster.GetLocalMembMinC(MinClusters)
            Koff = FtsCluster.KoffZ(LocalMembMinC)
            ToRemove = min(int(Koff), FtsCluster.SizeZ)
            if ToRemove <= 0:
                continue

            FtsCluster.SizeZ -= ToRemove
            TotalRemove += ToRemove
            NumCytoFtsZ[FtsCluster.SpatialIndex] += ToRemove

        # For DivPlane
        def RemoveFromDivPlane_Arc(RemainingRemove):
            Random = np.random.randint(0, NumAngles, RemainingRemove)
            for j in range(RemainingRemove):
                Idx = Random[j]
                if FtsClusters[Idx].Angle_Stop > FtsClusters[Idx].Angle:
                    FtsClusters[Idx].Angle_Stop -= pi / 10 / FtsClusters[Idx].Radius
                    RemainingRemove -= 1
            if RemainingRemove > 0:
                RemoveFromDivPlane_Arc(RemainingRemove)

        def RemoveFromDivPlane(RemainingRemove):
            Random = np.random.randint(0, NumAngles, RemainingRemove)
            for j in range(RemainingRemove):
                Idx = Random[j]
                if len(FtsClusters[Idx].Angle_FtsZs_DivPlane):
                    FtsClusters[Idx].Angle_FtsZs_DivPlane.pop()
                    FtsClusters[Idx].X_FtsZs_DivPlane.pop()
                    FtsClusters[Idx].Y_FtsZs_DivPlane.pop()
                    RemainingRemove -= 1
            if RemainingRemove > 0:
                RemoveFromDivPlane(RemainingRemove)

        if FtsClusters[0].ArcMode:
            RemoveFromDivPlane_Arc(TotalRemove)

        else:
            RemoveFromDivPlane(TotalRemove)

        # DL: Debugging
        # if TotalRemove > 0:
        #     for FtsCluster in FtsClusters:
        #         if len(FtsCluster.X_FtsZs_DivPlane):
        #             print(FtsCluster.Index, len(FtsCluster.Angle_FtsZs_DivPlane), FtsCluster.Angle_FtsZs_DivPlane)

    def TreadmillFtsZs():
        pass

    def DiffuseFtsAs():
        pass

    def DiffuseCytoMolecules():
        global NumCytoMinC
        global NumCytoFtsZ

        def DiffuseDistribution_4Cell(Distribution, D=0.01, dTime=1):  # D must be less than 1/6
            # von Neumann neighborhood
            Diffusion_Quantity = Distribution * D
            Diffusion_Padded = np.pad(Diffusion_Quantity, (1, 1))

            # roll diffusion quantity
            Upward = np.roll(Diffusion_Padded, -1, axis=0)
            Downward = np.roll(Diffusion_Padded, 1, axis=0)
            Leftward = np.roll(Diffusion_Padded, -1, axis=1)
            Rightward = np.roll(Diffusion_Padded, 1, axis=1)

            # Flow back from the edges
            Upward[1, :] += Upward[0, :]
            Downward[-2, :] += Downward[-1, :]
            Leftward[:, 1] += Leftward[:, 0]
            Rightward[:, -2] += Rightward[:, -1]

            Distribution_Diffused = (Upward + Downward + Leftward + Rightward)[1:-1, 1:-1] - 4 * Diffusion_Quantity

            return Distribution + Distribution_Diffused * dTime

        NumCytoMinC = DiffuseDistribution_4Cell(NumCytoMinC)
        NumCytoFtsZ = DiffuseDistribution_4Cell(NumCytoFtsZ)

    def MergeMinClusters():
        MinimumSize = 50
        ci = 0
        while ci < len(MinClusters):
            C1 = MinClusters[ci]
            if C1.Reference >= 0:
                ci += 1
                continue

            ni = ci + 1
            while True:
                # niw as ni_wrapped
                niw = ni % len(MinClusters)
                C2 = MinClusters[niw]
                if C2.Reference >= 0:
                    ni += 1
                    continue
                else:
                    break

            if ci == niw:
                break

            if C1.SizeD >= MinimumSize and C2.SizeD >= MinimumSize:
                if C1.SizeD > C2.SizeD:
                    C2.Reference = ci
                    C1.SizeD += C2.SizeD
                    C2.SizeD = 0
                    C1.SizeE += C2.SizeE
                    C2.SizeE = 0
                    C1.SizeC += C2.SizeC
                    C2.SizeC = 0
                    ci = ci
                else:
                    C1.Reference = niw
                    C2.SizeD += C1.SizeD
                    C1.SizeD = 0
                    C2.SizeE += C1.SizeE
                    C1.SizeE = 0
                    C2.SizeC += C1.SizeC
                    C1.SizeC = 0
                    ci = niw
                
            else:
                ci = niw

            if ni >= len(MinClusters):
                break

    def FragmentMinClusters():
        MinimumSize = 50
        ci = 0
        while ci < len(MinClusters):
            C1 = MinClusters[ci]
            if C1.Reference >= 0:
                ci += 1
                continue

            if C1.SizeD >= MinimumSize * 1.5 or C1.SizeD <= MinimumSize:
                ci += 1
                continue

            ni = ci + (1 if np.random.randint(2) == 0 else -1) + len(MinClusters)
            niw = ni % len(MinClusters)
            C2 = MinClusters[niw]
            if C2.Reference < 0:
                ci = ci + 1
                continue

            C2.SizeD = int(C1.SizeD / 2.5)
            C2.SizeE = int(C1.SizeE / 2.5)
            C2.SizeC = int(C1.SizeC / 2.5)
            C1.SizeD = C1.SizeD - C2.SizeD
            C1.SizeE = C1.SizeE - C2.SizeE
            C1.SizeC = C1.SizeC - C2.SizeC

            C2.Reference = -1
            
            ci = ci + 1

    def DrawCytoMolecules():
        global NumCytoMinD
        global NumCytoMinE
        global NumCytoMinC
        global NumCytoFtsZ

        def DrawMolecules_SingleCytoCompartment(NumCytoMolecule, color=BLACK, radius=5, additional_scale_factor=0.5):
            CytoMolecules = int(NumCytoMolecule * ScaleFactor_Quantity * additional_scale_factor)
            X = np.random.uniform(-1, 1, CytoMolecules) * MID_X + MID_X
            Y = np.random.uniform(-1, 1, CytoMolecules) * MID_Y + MID_Y
            for i in range(CytoMolecules):
                if CheckIfWithinEllipse(X[i], Y[i], Membrane.XYWH_DivAxis):
                    pygame.draw.circle(surface=Screen, color=color, center=(X[i], Y[i]), radius=radius)

        def DrawMolecules_MultiCytoCompartments(NumCytoMolecule, color=BLACK, radius=5, additional_scale_factor=0.5):
            '''
            NumCytoMolecule is a 2D numpy array.
            '''
            for row in range(NumCytoMolecule.shape[0]):
                for column in range(NumCytoMolecule.shape[1]):
                    CytoMolecules = int(NumCytoMolecule[row, column] * ScaleFactor_Quantity * additional_scale_factor)
                    X = np.random.uniform(-1, 1, CytoMolecules) * W_Ecoli / NumCytoMolecule.shape[0] / 2 + CytoCompartment_X[row, column]
                    Y = np.random.uniform(-1, 1, CytoMolecules) * H_Ecoli / NumCytoMolecule.shape[1] / 2 + CytoCompartment_Y[row, column]
                    for i in range(CytoMolecules):
                        if CheckIfWithinEllipse(X[i], Y[i], Membrane.XYWH_DivAxis):
                            pygame.draw.circle(surface=Screen, color=color, center=(X[i], Y[i]), radius=radius)

        # CytoMinD and CytoMinE
        DrawMolecules_SingleCytoCompartment(NumCytoMinD, color=Color_MinD, radius=Radius_MinD)
        DrawMolecules_SingleCytoCompartment(NumCytoMinE, color=Color_MinE, radius=Radius_MinE)

        # CytoMinC and CytoFtsZ
        DrawMolecules_MultiCytoCompartments(NumCytoMinC, color=Color_MinC, radius=Radius_MinC, additional_scale_factor=5)
        DrawMolecules_MultiCytoCompartments(NumCytoFtsZ, color=Color_FtsZ, radius=Radius_FtsZ, additional_scale_factor=0.1)

    Iter = 0
    IterMax = 20000
    while Iter < IterMax:
        # MinD
        AddMinDs()
        MergeMinClusters()
        RemoveMinDs()
        FragmentMinClusters()

        # MinE
        AddMinEs()
        RemoveMinEs()

        # MinC
        AddMinCs()
        RemoveMinCs()

        # FtsZ
        AddFtsZs()
        RemoveFtsZs()
        TreadmillFtsZs()

        # FtsA
        DiffuseFtsAs()

        # Diffusion: MinC, FtsZ
        DiffuseCytoMolecules()

        if Iter % 100 == 0:
        # if True:
            NumMemMinD = 0
            NumMemMinE = 0
            NumMemMinC = 0
            NumMemFtsA = 0
            NumMemFtsZ = 0
            print("Iter {}".format(Iter))
            for i, MinCluster in enumerate(MinClusters):
                Angle = MinCluster.Angle
                SizeD = MinCluster.SizeD
                SizeE = MinCluster.SizeE
                SizeC = MinCluster.SizeC
                Reference = MinCluster.Reference
                if SizeD <= 0 and Reference < 0:
                    continue

                if Reference < 0:
                    KonD = MinCluster.KonD()
                    KoffD = MinCluster.KoffD()
                    KonE = MinCluster.KonE()
                    KoffE = MinCluster.KoffE()
                    KonC = MinCluster.KonC()
                    KoffC = MinCluster.KoffC()
                    Log = "MinCluster #{:<5d} rad: {:.2f}, MinD: {:5d}, KonD: {:.2f}, KoffD: {:.2f}, MinE: {:5d}, KonE: {:.2f}, KoffE: {:.2f}, MinC: {:5d}, KonC: {:.2f}, KoffC: {:.2f}".format(
                            i, Angle, SizeD, KonD, KoffD, SizeE, KonE, KoffE, SizeC, KonC, KoffC)
                    print(Log)
                    NumMemMinD += SizeD
                    NumMemMinE += SizeE
                    NumMemMinC += SizeC
                else:
                    print("MinCluster #{:<5d} Reference: #{:<5d}".format(i, Reference))

            for i, FtsCluster in enumerate(FtsClusters):
                Angle = FtsCluster.Angle
                SizeA = FtsCluster.SizeA
                SizeZ = FtsCluster.SizeZ
                Reference = FtsCluster.Reference
                if SizeZ <= 0 and SizeA <= 0:
                    continue

                if Reference < 0:
                    KonA = FtsCluster.KonA()
                    KoffA = FtsCluster.KoffA()
                    KonZ = FtsCluster.KonZ(NumCytoMinC[FtsCluster.SpatialIndex][0])
                    KoffZ = FtsCluster.KoffZ(FtsCluster.GetLocalMembMinC(MinClusters))
                    Log = "FtsCluster #{:<5d} rad: {:.2f}, FtsA: {:5d}, KonA: {:.2f}, KoffA: {:.2f}, FtsZ: {:5d}, KonZ: {:.2f}, KoffZ: {:.2f}".format(
                            i, Angle, SizeA, KonA, KoffA, SizeZ, KonZ, KoffZ)
                    print(Log + "\t| Min Ref: #{:<5d}, Fts Ref: #{:<5d}".format(MinClusters[FtsCluster.Index].Reference, Reference))
                    NumMemFtsZ += SizeZ
                    NumMemFtsA += SizeA
                else:
                    print("FtsCluster #{:<5d} Reference: #{:<5d}".format(i, Reference))

            print("cMinD: {:5d}, mMinD: {:5d}".format(NumCytoMinD, NumMemMinD))
            assert NumCytoMinD + NumMemMinD == NumMinD
            print("cMinE: {:5d}, mMinE: {:5d}".format(NumCytoMinE, NumMemMinE))
            assert NumCytoMinE + NumMemMinE == NumMinE
            Sum_NumCytoMinC = round(np.sum(NumCytoMinC))
            print("cMinC: {:5d}, mMinC: {:5d}".format(Sum_NumCytoMinC, NumMemMinC))
            assert Sum_NumCytoMinC + NumMemMinC == NumMinC
            Sum_NumCytoFtsZ = round(np.sum(NumCytoFtsZ))
            print("cFtsZ: {:5d}, mFtsZ: {:5d}, mFtsZ/FtsZ: {:2f}".format(Sum_NumCytoFtsZ, NumMemFtsZ, NumMemFtsZ/NumFtsZ))
            assert Sum_NumCytoFtsZ + NumMemFtsZ == NumFtsZ
            print("mFtsA: {:5d}".format(NumMemFtsA))
            assert NumMemFtsA == NumFtsA - (NumFtsA % NumAngles)
            print("\n")

        # Draw
        if Iter % ScaleFactor_Refresh == 0:
            Screen.fill(WHITE)
            Membrane.Draw()
            for MinCluster in MinClusters:
                if MinCluster.Reference < 0:
                    MinCluster.Draw()
                # MinCluster.Draw(show_anchor=True)
            for FtsCluster in FtsClusters:
                if FtsCluster.SizeA > 0 or FtsCluster.SizeZ > 0:
                    # FtsCluster.Draw()
                    FtsCluster.Draw(show_anchor=True)
                    FtsCluster.Draw_DivPlane(show_anchor=True)

            DrawCytoMolecules()

            DisplayIteration(Iter)
            Legends.Draw()
            # Legends.Draw(show_anchor=True)
            pygame.display.update()

        Iter += 1

    pygame.quit()
    sys.exit()


if __name__ == '__main__':
    Simulate()
