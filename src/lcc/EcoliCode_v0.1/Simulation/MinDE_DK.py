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

pi = np.pi

def DisplayTime(Time):
    Text = Font_Sans.render('Simulation Time: ' + str(round(Time)), True, BLACK)
    Text_Rect = Text.get_rect()
    Text_Rect.midtop = Screen.get_rect().midtop
    Screen.blit(Text, Text_Rect)

def GetXYFromEllipse(Shape, Angle, ratio_factor=1.0):
    Center_X, Center_Y, Width, Height = Shape
    return Center_X + Width / 2 * ratio_factor * np.cos(Angle), Center_Y + Height / 2 * ratio_factor * np.sin(Angle)

def GetXYFromCenter(Center, Angle, Radius, anisotropy=1.0):
    Center_X, Center_Y = Center
    return Center_X + Radius * anisotropy * np.cos(Angle), Center_Y + Radius * np.sin(Angle)

def CheckIfWithinEllipse(X, Y, XYWH, membranebias=0.8):
    Center_X, Center_Y, Width, Height = XYWH
    Evaluate = (X - Center_X) ** 2 / (Width / 2) ** 2 + (Y - Center_Y) ** 2 / (Height / 2) ** 2
    if Evaluate < 1:
        if membranebias:
            if Evaluate > membranebias:
                return True
            else:
                return False
        else:
            return True
    else:
        return False


# pygame
pygame.init()

Title = 'MinDE Visualization'
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
ScaleFactor_Dimension = 1.5
W_Ecoli = int(2000 / ScaleFactor_Dimension)
H_Ecoli = int(1000 / ScaleFactor_Dimension)

# Molecules
NumMinD = 5000
NumCytoMinD = NumMinD

NumMinE = 4000
NumCytoMinE = NumMinE

# Cytosolic compartment (currently support odd numbers only to have one absolute center compartment)
CytoCompartment_Rows = 3
CytoCompartment_Columns = 3
CytoCompartment_Width = W_Ecoli / CytoCompartment_Rows
CytoCompartment_Height = H_Ecoli / CytoCompartment_Columns

# Molecules distributed in cytosolic compartments
NumMinC = 200
NumCytoMinC = np.zeros([CytoCompartment_Rows, CytoCompartment_Columns])
NumCytoMinC[CytoCompartment_Rows // 2, CytoCompartment_Columns // 2] = NumMinC

NumFtsZ = 10000
NumCytoFtsZ = np.zeros([CytoCompartment_Rows, CytoCompartment_Columns])
NumCytoFtsZ[CytoCompartment_Rows // 2, CytoCompartment_Columns // 2] = NumFtsZ

# Drawing
ScaleFactor_Quantity = 0.4
ScaleFactor_Refresh = 30

Size_Anchor = 5 
Size_MinD = 5
Size_MinE = 5
Size_MinC = 5
Size_FtsZ = 5

Color_MinD = RED
Color_MinE = BLUE
Color_MinC = YELLOW
Color_FtsZ = MAGENTA

CytoCompartment_X = np.zeros([CytoCompartment_Rows, CytoCompartment_Columns])
CytoCompartment_Y = np.zeros([CytoCompartment_Rows, CytoCompartment_Columns])
for row in range(CytoCompartment_Rows):
    Center_X = np.round(np.floor(CytoCompartment_Rows / 2))
    Offset_X = W_Ecoli / CytoCompartment_Rows * (row - Center_X)
    for column in range(CytoCompartment_Columns):
        Center_Y = np.round(np.floor(CytoCompartment_Columns / 2))
        Offset_Y = H_Ecoli / CytoCompartment_Columns * (column - Center_Y)

        CytoCompartment_X[row, column] = MID_X + Offset_X
        CytoCompartment_Y[row, column] = MID_Y + Offset_Y

class FMembrane:
    def __init__(self, X=MID_X, Y=MID_Y, thickness=5, linecolor=GRAY4, bodycolor=YELLOW_FAINT, Label=False, width=0, height=0):
        self.X = X
        self.Y = Y
        self.W = width
        self.H = height
        self.XY = X, Y
        self.XYWH = X, Y, width, height

        self.Label = Label
        self.Thickness = thickness
        self.LineColor = linecolor
        self.BodyColor = bodycolor

    def Draw(self):
        pygame.draw.ellipse(Screen, self.BodyColor, (self.X - self.W / 2, self.Y - self.H / 2, self.W, self.H))
        pygame.draw.ellipse(Screen, self.LineColor, (self.X - self.W / 2, self.Y - self.H / 2, self.W, self.H), self.Thickness)


class FMinCluster:
    def __init__(self, Index, Angle, Reference, SizeD, SizeE, SizeC, memb=None, color=GREEN):
        '''
        memb is X, Y, W, H of the membrane
        '''

        self.Index = Index
        self.Angle = Angle
        self.Reference = Reference
        self.SizeD = SizeD
        self.SizeE = SizeE
        self.SizeC = SizeC

        # for Drawing
        self.X_Lipid, self.Y_Lipid = GetXYFromEllipse(memb, Angle, ratio_factor=1.0)
        self.X_Min, self.Y_Min = GetXYFromEllipse(memb, Angle, ratio_factor=0.75)
        self.Color = color

    def GetSpatialIdx(self):
        Idx_X = np.logical_and(self.X_Min < CytoCompartment_X + CytoCompartment_Width / 2, self.X_Min > CytoCompartment_X - CytoCompartment_Width / 2)
        Idx_Y = np.logical_and(self.Y_Min < CytoCompartment_Y + CytoCompartment_Height / 2, self.Y_Min > CytoCompartment_Y - CytoCompartment_Height / 2)
        Idx_XY = np.logical_and(Idx_X, Idx_Y)
        return np.where(Idx_XY)

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

    def Draw(self, XYWH):
        X_Membrane, Y_Membrane, W_Membrane, H_Membrane = XYWH
        pygame.draw.circle(surface=Screen, color=self.Color, center=(self.X_Lipid, self.Y_Lipid), radius=Size_Anchor)

        # MinD
        SizeD = int(self.SizeD * ScaleFactor_Quantity * 2)
        X_MinD = np.random.uniform(-1, 1, SizeD) * W_Membrane / 4 + MID_X + W_Membrane / 4 * (-1 if self.X_Lipid < MID_X else 1)
        Y_MinD = np.random.normal(MID_Y, H_Membrane / 6, SizeD)
        for i in range(SizeD):
            if CheckIfWithinEllipse(X_MinD[i], Y_MinD[i], XYWH, membranebias=0.8):
                pygame.draw.circle(surface=Screen, color=Color_MinD, center=(X_MinD[i], Y_MinD[i]), radius=Size_MinD)

        # MinE
        SizeE = int(self.SizeE * ScaleFactor_Quantity * 2)
        X_MinE = np.random.uniform(-1, 1, SizeE) * W_Membrane / 4 + MID_X + W_Membrane / 4 * (-1 if self.X_Lipid < MID_X else 1)
        Y_MinE = np.random.normal(MID_Y, H_Membrane / 6, SizeE)
        for i in range(SizeE):
            if CheckIfWithinEllipse(X_MinE[i], Y_MinE[i], XYWH, membranebias=0.8):
                pygame.draw.circle(surface=Screen, color=Color_MinE, center=(X_MinE[i], Y_MinE[i]), radius=Size_MinE)

        # MinC
        SizeC = int(self.SizeC * ScaleFactor_Quantity * 5)
        X_MinC = np.random.uniform(-1, 1, SizeC) * W_Membrane / 4 + MID_X + W_Membrane / 4 * (-1 if self.X_Lipid < MID_X else 1)
        Y_MinC = np.random.normal(MID_Y, H_Membrane / 6, SizeC)
        for i in range(SizeC):
            if CheckIfWithinEllipse(X_MinC[i], Y_MinC[i], XYWH, membranebias=0.8):
                pygame.draw.circle(surface=Screen, color=Color_MinC, center=(X_MinC[i], Y_MinC[i]), radius=Size_MinC)


def SimulateMinCDE():
    # Initialize Membrane
    Membrane = FMembrane(X=MID_X, Y=MID_Y, width=W_Ecoli, height=H_Ecoli, Label=True)

    # Initialize Molecules
    NumAngles = 100
    DeltaAngle = 2 * pi / NumAngles
    MinClusters = [FMinCluster(i, DeltaAngle * i, -1, 0, 0, 0, memb=Membrane.XYWH) for i in range(NumAngles)]

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

        KonSum = 0
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Kon = MinCluster.KonC()
            KonSum += Kon

        if KonSum <= np.sum(NumCytoMinC):
            ScaleFactor = 1
        else:
            ScaleFactor = np.sum(NumCytoMinC) / KonSum

        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Kon = MinCluster.KonC()
            Kon *= ScaleFactor
            LocalCytoIdx = MinCluster.GetSpatialIdx()
            ToAdd = min(int(Kon), int(NumCytoMinC[LocalCytoIdx][0]))
            if ToAdd <= 0:
                continue

            MinCluster.SizeC += ToAdd
            NumCytoMinC[LocalCytoIdx] -= ToAdd

    def RemoveMinCs():
        global NumCytoMinC
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Koff = MinCluster.KoffC()
            LocalCytoIdx = MinCluster.GetSpatialIdx()
            ToRemove = min(int(Koff), MinCluster.SizeC)
            if ToRemove <= 0:
                continue

            MinCluster.SizeC -= ToRemove
            NumCytoMinC[LocalCytoIdx] += ToRemove

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
                if CheckIfWithinEllipse(X[i], Y[i], Membrane.XYWH, membranebias=0):
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
                        if CheckIfWithinEllipse(X[i], Y[i], Membrane.XYWH, membranebias=0):
                            pygame.draw.circle(surface=Screen, color=color, center=(X[i], Y[i]), radius=radius)

        # CytoMinD and CytoMinE
        DrawMolecules_SingleCytoCompartment(NumCytoMinD, color=Color_MinD, radius=Size_MinD)
        DrawMolecules_SingleCytoCompartment(NumCytoMinE, color=Color_MinE, radius=Size_MinE)

        # CytoMinC and CytoFtsZ
        DrawMolecules_MultiCytoCompartments(NumCytoMinC, color=Color_MinC, radius=Size_MinC, additional_scale_factor=5)
        DrawMolecules_MultiCytoCompartments(NumCytoFtsZ, color=Color_FtsZ, radius=Size_FtsZ, additional_scale_factor=0.05)

    Iter = 0
    IterMax = 20000
    while Iter < IterMax:
        AddMinDs()

        MergeMinClusters()        

        RemoveMinDs()

        FragmentMinClusters()

        AddMinEs()

        RemoveMinEs()

        AddMinCs()

        RemoveMinCs()

        DiffuseCytoMolecules()

        if Iter % 100 == 0:
        # if Iter > 0:
            NumMemMinD = 0
            NumMemMinE = 0
            NumMemMinC = 0
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
                    Log = "#{:<5d} rad: {:.2f}, MinD: {:5d}, KonD: {:.2f}, KoffD: {:.2f}, MinE: {:5d}, KonE: {:.2f}, KoffE: {:.2f}, MinC: {:5d}, KonC: {:.2f}, KoffC: {:.2f}".format(
                            i, Angle, SizeD, KonD, KoffD, SizeE, KonE, KoffE, SizeC, KonC, KoffC, Reference)
                    print(Log)
                    NumMemMinD += SizeD
                    NumMemMinE += SizeE
                    NumMemMinC += SizeC
                else:
                    print("#{:<5d} Reference: #{:<5d}".format(i, Reference))

            print("cMinD: {:5d}, mMinD: {:5d}".format(NumCytoMinD, NumMemMinD))
            assert NumCytoMinD + NumMemMinD == NumMinD
            print("cMinE: {:5d}, mMinE: {:5d}".format(NumCytoMinE, NumMemMinE))
            assert NumCytoMinE + NumMemMinE == NumMinE
            Sum_NumCytoMinC = round(np.sum(NumCytoMinC))
            print("cMinC: {:5d}, mMinC: {:5d}".format(Sum_NumCytoMinC, NumMemMinC))
            assert Sum_NumCytoMinC + NumMemMinC == NumMinC
            Sum_NumCytoFtsZ = round(np.sum(NumCytoFtsZ))
            print("cFtsZ: {:5d}, mFtsZ: {:5d}".format(Sum_NumCytoFtsZ, NumMemFtsZ))
            assert Sum_NumCytoFtsZ + NumMemFtsZ == NumFtsZ
            print("\n")

        # Draw
        if Iter % ScaleFactor_Refresh == 0:
            Screen.fill(WHITE)
            Membrane.Draw()
            for MinCluster in MinClusters:
                MinCluster.Draw(Membrane.XYWH)
            DrawCytoMolecules()

            DisplayTime(Iter)
            pygame.display.update()

        Iter += 1

    pygame.quit()
    sys.exit()


if __name__ == '__main__':
    SimulateMinCDE()
