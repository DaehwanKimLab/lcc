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

# Molecule
NumMinD = 5000
NumCytoMinD = NumMinD
NumMinE = 4000
NumCytoMinE = NumMinE

# Drawing Molecules
ScaleFactor_Quantity = 0.4
ScaleFactor_Refresh = 30
Size_Anchor = 5 
Size_MinD = 5
Size_MinE = 5

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
    def __init__(self, Index, Angle, Reference, SizeD, SizeE, memb=None, color=GREEN):
        '''
        memb is X, Y, W, H of the membrane
        '''

        self.Index = Index
        self.Angle = Angle
        self.Reference = Reference
        self.SizeD = SizeD
        self.SizeE = SizeE

        # for Drawing
        self.X_Lipid, self.Y_Lipid = GetXYFromEllipse(memb, Angle, ratio_factor=1.0)
        self.X_Min, self.Y_Min = GetXYFromEllipse(memb, Angle, ratio_factor=0.75)
        self.Color = color

    def MembraneEffectFunc(self, Rad):
        # return math.cosh(Rad) - 0.9
        return Rad * Rad
        
    def KonD(self):
        MaxMembraneEffect = self.MembraneEffectFunc(0.5 * pi)
        NewAngle = self.Angle if self.Angle < pi else self.Angle - pi
        MembraneEffect = self.MembraneEffectFunc(NewAngle - 0.5 * pi) / MaxMembraneEffect * 1.1
        AggregationEffect = math.log(max(2, self.SizeD), 2)
        return MembraneEffect * AggregationEffect

    def KoffD(self):
        return self.KonD() * (self.SizeE / max(1, self.SizeD)) * np.random.normal(1.5, 0.2)

    def KonE(self):
        return self.KonD() * 0.2

    def KoffE(self):
        return max(0, self.SizeE - self.SizeD)

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
                pygame.draw.circle(surface=Screen, color=RED, center=(X_MinD[i], Y_MinD[i]), radius=Size_MinD)

        # MinE
        SizeE = int(self.SizeE * ScaleFactor_Quantity * 2)
        X_MinE = np.random.uniform(-1, 1, SizeE) * W_Membrane / 4 + MID_X + W_Membrane / 4 * (-1 if self.X_Lipid < MID_X else 1)
        Y_MinE = np.random.normal(MID_Y, H_Membrane / 6, SizeE)
        for i in range(SizeE):
            if CheckIfWithinEllipse(X_MinE[i], Y_MinE[i], XYWH, membranebias=0.8):
                pygame.draw.circle(surface=Screen, color=BLUE, center=(X_MinE[i], Y_MinE[i]), radius=Size_MinE)


def SimulateMinCDE():
    # Initialize Membrane
    Membrane = FMembrane(X=MID_X, Y=MID_Y, width=W_Ecoli, height=H_Ecoli, Label=True)

    # Initialize Molecules
    NumAngles = 100
    DeltaAngle = 2 * pi / NumAngles
    MinClusters = [FMinCluster(i, DeltaAngle * i, -1, 0, 0, memb=Membrane.XYWH) for i in range(NumAngles)]

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

            if C1.SizeD > 0 and C2.SizeD > 0 and False:
                print("C1 {}".format(ci))
                C1.PrintInfo()
                print("C2 {}".format(ni))
                C2.PrintInfo()

            if C1.SizeD >= MinimumSize and C2.SizeD >= MinimumSize:
                if C1.SizeD > C2.SizeD:
                    C2.Reference = ci
                    C1.SizeD += C2.SizeD
                    C1.SizeE += C2.SizeE
                    ci = ci
                else:
                    C1.Reference = niw
                    C2.SizeD += C1.SizeD
                    C2.SizeE += C1.SizeE
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

            if C1.SizeD >= MinimumSize * 2 or C1.SizeD <= MinimumSize:
                ci += 1
                continue

            ni = ci + 1
            niw = ni % len(MinClusters)
            C2 = MinClusters[niw]
            if C2.Reference < 0:
                ci = ni
                continue

            C1SizeD = int(C1.SizeD / 2)
            C1SizeE = int(C1.SizeE / 2)

            C2.SizeD = C1.SizeD - C1SizeD
            C2.SizeE = C1.SizeE - C1SizeE
            C1.SizeD = C1SizeD
            C1.SizeE = C1SizeE

            C2.Reference = -1
            
            ci = ni + 1

    def DrawCytoMin():
        global NumCytoMinD
        global NumCytoMinE

        # CytoMinD
        CytoMinD = int(NumCytoMinD * ScaleFactor_Quantity * 0.5)
        X_MinD = np.random.uniform(-1, 1, CytoMinD) * MID_X + MID_X
        Y_MinD = np.random.uniform(-1, 1, CytoMinD) * MID_Y + MID_Y
        for i in range(CytoMinD):
            if CheckIfWithinEllipse(X_MinD[i], Y_MinD[i], Membrane.XYWH, membranebias=0):
                pygame.draw.circle(surface=Screen, color=RED, center=(X_MinD[i], Y_MinD[i]), radius=Size_MinD)

        # CytoMinE
        CytoMinE = int(NumCytoMinE * ScaleFactor_Quantity * 0.5)
        X_MinE = np.random.uniform(-1, 1, CytoMinE) * MID_X + MID_X
        Y_MinE = np.random.uniform(-1, 1, CytoMinE) * MID_Y + MID_Y
        for i in range(CytoMinE):
            if CheckIfWithinEllipse(X_MinE[i], Y_MinE[i], Membrane.XYWH, membranebias=0):
                pygame.draw.circle(surface=Screen, color=BLUE, center=(X_MinE[i], Y_MinE[i]), radius=Size_MinE)

    Iter = 0
    while Iter < 20000:
        AddMinDs()

        MergeMinClusters()        

        RemoveMinDs()

        # FragmentMinClusters()

        AddMinEs()

        RemoveMinEs()

        if Iter % 100 == 0:
            NumMemMinD = 0
            NumMemMinE = 0
            print("Iter {}".format(Iter))
            for i, MinCluster in enumerate(MinClusters):
                Angle = MinCluster.Angle
                SizeD = MinCluster.SizeD
                SizeE = MinCluster.SizeE
                Reference = MinCluster.Reference
                if SizeD <= 0 and Reference < 0:
                    continue

                if Reference < 0:
                    KonD = MinCluster.KonD()
                    KoffD = MinCluster.KoffD()
                    KonE = MinCluster.KonE()
                    KoffE = MinCluster.KoffE()
                    print("#{:<5d} rad: {:.2f}, MinD: {:5d}, KonD: {:.2f}, KoffD: {:.2f}, MinE: {:5d}, KonE: {:.2f}, KoffE: {:.2f}".format(i, Angle, SizeD, KonD, KoffD, SizeE, KonE, KoffE))
                    NumMemMinD += SizeD
                    NumMemMinE += SizeE
                else:
                    print("#{:<5d} Reference: #{:<5d}".format(i, Reference))

            print("cMinD: {:5d}, mMinD: {:5d}".format(NumCytoMinD, NumMemMinD))
            assert NumCytoMinD + NumMemMinD == NumMinD
            print("cMinE: {:5d}, mMinE: {:5d}".format(NumCytoMinE, NumMemMinE))
            assert NumCytoMinE + NumMemMinE == NumMinE
            print("\n")

        # Draw
        if Iter % ScaleFactor_Refresh == 0:
            Screen.fill(WHITE)
            Membrane.Draw()
            for MinCluster in MinClusters:
                MinCluster.Draw(Membrane.XYWH)
            DrawCytoMin()

            DisplayTime(Iter)
            pygame.display.update()

        Iter += 1

    pygame.quit()
    sys.exit()


if __name__ == '__main__':
    SimulateMinCDE()
