# BSD 3-Clause License
# © 2023, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee, Chanhee Park and Daehwan Kim

#!/usr/bin/python3
import os, time
import math
from enum import Enum

import numpy as np
from matplotlib import image
from matplotlib import pyplot

import pygame

class MassObjectPattern(Enum):
    UNIFORM = 1
    DIFFUSION = 2

class MassObjectDimension(Enum):
    TWO = 1
    THREE = 2

class MassObject:
    def __init__(self, Type = MassObjectPattern.DIFFUSION):
        self.InoculationResolution = 1e-6 # 1um
        self.ScaleFactor = 1000

        self.Time = 0.0
        self.Type = Type
        self.Amount = 1 # 1 mol
        self.DiffusionSpeed = 1e-2 # 1cm
        self.MaxIntercept = 0
        self.UpdateInterceptSlope()

    def UpdateInterceptSlope(self):
        self.Radius = self.InoculationResolution + self.DiffusionSpeed * self.Time
        self.Intercept = 3.0 * self.Amount / (math.pi * self.Radius * self.Radius)
        self.Slope = - self.Intercept / self.Radius
        if self.Time <= 10.0:
            self.MaxIntercept = self.Intercept

    def Simulate(self, ElapsedTime):
        self.Time += ElapsedTime
        self.UpdateInterceptSlope()

    def Draw(self, Screen, Left, Right, Top, Bottom):
        CenterX = (Left + Right) / 2
        self.DrawMass(Screen, Left, CenterX, Top, Bottom)
        self.DrawPlot(Screen, CenterX, Right, Top, Bottom)

    def DrawMass(self, Screen, Left, Right, Top, Bottom):
        CenterX = (Left + Right) / 2
        CenterY = (Top + Bottom) / 2

        Radius = 30
        while CenterX + Radius < Right and CenterY + Radius < Bottom:
            RealRadius = Radius / self.ScaleFactor
            Conc = self.Intercept + self.Slope * RealRadius
            Conc = max(0, Conc)
            Factor = Conc / self.Intercept
            if Factor > 0:
                Color = (0, 0, 255 * Factor)
                pygame.draw.circle(Screen, Color, (CenterX, CenterY), Radius, 2)

            if Radius % 30 == 0:
                TextCenter = (CenterX + Radius, CenterY + Radius)
                self.DrawText(Screen, "R: {:.2f}, C: {:.2f}".format(RealRadius, Conc), TextCenter)

            Radius += 10

    def DrawPlot(self, Screen, Left, Right, Top, Bottom):
        Height = Bottom - Top
        Top2 = Top + Height * 0.1
        Bottom2 = Top2 + Height * 0.8
        Height2 = Bottom2 - Top2
        LeftY = Top2 + Height2 * (1.0 - self.Intercept / self.MaxIntercept)
        self.DrawText(Screen, "Y: {:.2f}".format(self.MaxIntercept), (Left, Top + 10))
        self.DrawText(Screen, "Y: {:.2f}".format(self.Intercept), (Left, LeftY + 10))
        self.DrawText(Screen, "R: {:.2f}".format(self.Radius), (Left + self.Radius * self.ScaleFactor, Bottom2 + 10))

        # X axis
        pygame.draw.line(Screen, (0, 0, 0), (Left, Bottom2), (Right, Bottom2))

        # Y axis
        pygame.draw.line(Screen, (0, 0, 0), (Left, Top), (Left, Bottom2))
       
        # Line
        pygame.draw.line(Screen, (255, 0, 0), (Left, LeftY), (Left + self.Radius * self.ScaleFactor, Bottom2))

    def DrawText(self, Screen, Text, Center):
        Font = pygame.font.Font('freesansbold.ttf', 12)
        Text = Font.render(Text, True, (0, 0, 0))
        TextRect = Text.get_rect()
        TextRect.center = Center
        Screen.blit(Text, TextRect)
            
    
class Octree:
    def __init__(self):
        self.Root = MassObject()

    def Simulate(self, ElapsedTime):
        self.Root.Simulate(ElapsedTime)

    def Draw(self, Screen, Left, Right, Top, Bottom):
        Screen.fill((255, 255, 255))
        self.Root.Draw(Screen, Left, Right, Top, Bottom)


if __name__ == "__main__":
    OctreeInst = Octree()
    
    pygame.init()
    ScreenSize = (1200, 600)
    Screen = pygame.display.set_mode(ScreenSize)

    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        OctreeInst.Draw(Screen, 0, ScreenSize[0], 0, ScreenSize[1])
        pygame.display.flip()

        OctreeInst.Simulate(0.5)
        time.sleep(0.5)

    pygame.quit()
