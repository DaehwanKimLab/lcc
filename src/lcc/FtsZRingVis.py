import sys
import pygame
import random
from datetime import datetime
import numpy as np
import SimModule
import SimState as SimS
import SimData as SimD
import SimFunctions as SimF

random.seed(1)
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

# Global variables
NA = 6.0221409e+23
pi = np.pi

# Utilities
def GetColorGradient(Fade, baseColor=None):
    assert Fade <= 255, 'ERROR: Fade factor is higher than 255!'
    if baseColor == 'Blue':
        return (Fade, Fade, 255)
    else:  # Assumes White
        return (Fade, Fade, Fade)

# pygame
pygame.init()

Screen_Size = W_S, H_S = 1000, 1000
Screen = pygame.display.set_mode(Screen_Size)

LEFT = 0
MID_X = W_S / 2
RIGHT = W_S

TOP = 0
MID_Y = H_S / 2
BOTTOM = H_S

CenterTop = (MID_X, TOP)
Center = (MID_X, MID_Y)

def AddPos(A, B):
    X_A, Y_A = A
    X_B, Y_B = B
    return (X_A + X_B, Y_A + Y_B)

def GetMidPoint(A, B):
    X_A, Y_A = A
    X_B, Y_B = B
    return ((X_A + X_B) / 2, (Y_A + Y_B) / 2)

def GetXYFromCenter(Center, Angle, Radius, RatioFactor=1.0):
    Center_X, Center_Y = Center
    return Center_X + Radius * RatioFactor * np.cos(Angle), Center_Y + Radius * RatioFactor * np.sin(Angle)

Title = 'FtsZ Ring Visualization'
pygame.display.set_caption(Title)

Font_Init = pygame.font.SysFont('arial', 50)
Font_MassObject = pygame.font.SysFont('arial', 20)
Font_IndividualObject = pygame.font.SysFont('arial', 15)

Font_Sans = pygame.font.Font('freesansbold.ttf', 20)
Font_Monospace = pygame.font.SysFont('monospace', 18, True)
Font_Radar = pygame.font.SysFont('arial', 11)

class FCompartment:
    def __init__(self, InX=W_S/2, InY=H_S/2, InShape='circle', InRadius=H_S*8/10/2, InThickness=5, linecolor=GRAY4, bodycolor=YELLOW_FAINT):
        self.X = InX
        self.Y = InY
        self.Shape = InShape
        self.Radius = int(InRadius)
        self.Thickness = InThickness
        self.LineColor = linecolor
        self.BodyColor = bodycolor

    def Draw(self):
        if self.Shape == 'circle':
            pygame.draw.circle(Screen, self.BodyColor, (self.X, self.Y), self.Radius)
            pygame.draw.circle(Screen, self.LineColor, (self.X, self.Y), self.Radius, self.Thickness)

class FObject:
    def __init__(self, Name, type='individual', id=0, x=W_S/2, y=H_S/2, angle=0, color=BLACK, radius=10, quantity=100, width=20, height=10, thickness=2):
        self.Name = Name
        self.Type = type
        self.ID = id
        self.X = x
        self.Y = y
        self.Angle = angle
        self.TotalQuantity = quantity
        self.CytosolicQuantity = quantity

        # ChildrenObject
        self.N_Children = 0
        self.X_Children = list()
        self.Y_Children = list()
        self.Angle_Children = list()

        # Draw
        self.Radius = radius
        self.Width = width
        self.Height = height
        self.Thickness = thickness
        self.Color = color

    def Draw(self):
        if self.Type == 'individual':
            pygame.draw.circle(surface=Screen, color=GREEN, center=(self.X, self.Y), radius=self.Radius)
            Label = self.Name + "#"+ str(self.ID) + "~FtsZ_{" + str(self.N_Children) + "}"
            self.DisplayString(Label, self.X+self.Radius*2.5*np.cos(self.Angle), self.Y+self.Radius*2.5*np.sin(self.Angle), color=BLACK)

            for i in range(len(self.X_Children)):
                pygame.draw.circle(surface=Screen, color=RED, center=(self.X_Children[i], self.Y_Children[i]), radius=self.Radius*0.5)
                # self.DisplayString(self.Name + str(i), self.X_Children[i], self.Y_Children[i], color=BLACK)

        elif self.Type == 'mass':
            Text = "Cytosolic " + self.Name + ": " + str(self.CytosolicQuantity)
            self.DisplayString(Text, W_S/2, H_S/2, color=BLACK)

        else:
            self.DisplayString("DRAW ERROR", W_S/2, H_S/4, color=BLACK)


    # def DrawLegend(self, Color, X, Y):
    #     for i in range(len(self.Radar_MolList)):
    #         dY = 0
    #         if len(self.Radar_MolList) > 1:
    #             dY = 16 * i - 8   # Hardcoded pattern for two molecules
    #         self.DisplayRadarText(self.Radar_MolList[i], X, Y + dY, color=self.Radar_MolColor[i])
    #         for j in range(self.Radar_Sampling):
    #             Spacing = self.Radar_Spacing * (j + 1)
    #             pygame.draw.circle(Screen, Color, (X, Y), Spacing, 1)
    #             RotationAngle = int(i * 90 / len(self.Radar_MolList))
    #             self.DisplayRadarTexts(self.Radar_MolList[i], X, Y, Spacing, rotate=RotationAngle, color=self.Radar_MolColor[i])
    #
    # def DisplayRingCompletionRate(self, X, Y, Rate):
    #     if Y >= 30 and Y < H_S:
    #         Value_Str = self.FormatValueToStr(Rate * 100) + '%'
    #         self.DisplayString(Value_Str, X, Y - 20)

    def FormatValueToStr(self, Value):
        if Value < 0.01:
            return SimF.SciFloat(Value, InPrecision=2, InExp_digits=2)
        else:
            return '{:.3f}'.format(Value)

    def DisplayString(self, String, X, Y, position='center', rotate=0, color=BLACK):
        Text=''
        if self.Type == 'mass':
            Text = Font_MassObject.render(String, True, color)
        elif self.Type == 'individual':
            Text = Font_IndividualObject.render(String, True, color)
        Text = pygame.transform.rotate(Text, rotate)
        Text_Rect = Text.get_rect()
        if position == 'center':
            Text_Rect.center = (X, Y)
        if position == 'midleft':
            Text_Rect.midleft = (X, Y)
        elif position == 'midright':
            Text_Rect.midright = (X, Y)
        Screen.blit(Text, Text_Rect)

class FControl:
    def __init__(self):
        self.MovementResolution = 2
        self.MessageWelcome = ''
        self.Pos_Welcome = GetMidPoint(Center, CenterTop)
        self.Message = ''
        self.MessageTimer = 3000

        self.InstructionSwitch = False
        self.Instructions = {
            'X' : 'Exit Visualization',
            'P' : 'Pause Visualization',
        }
        self.InitText = ''

        self.Time = 0

        # Pause
        self.MessagePause = 'PAUSE'
        self.PauseSwitch = False

    def DisplayInit(self):
        Screen.fill(GRAY1)
        Text = Font_Init.render(self.InitText, True, BLACK)
        Screen.blit(Text, Center)
        pygame.display.update()

    def DisplayTime(self):
        Text = Font_Sans.render('Simulation Time: ' + str(round(self.Time)), True, BLACK)
        Text_Rect = Text.get_rect()
        Text_Rect.midtop = Screen.get_rect().midtop
        Screen.blit(Text, Text_Rect)
    #
    # def DisplayNumberOfOrganisms(self, NumberOfOrganisms):
    #     Text = Font_Sans.render('# of Organisms: ' + str(round(NumberOfOrganisms)), True, BLUE)
    #     Text_Rect = Text.get_rect()
    #     Text_Rect.midtop = tuple(map(lambda i, j: i + j, Screen.get_rect().midtop, (0, Text.get_height() * 1.5)))
    #     Screen.blit(Text, Text_Rect)
    #

    def SetInstructionText(self):
        self.InstructionText = 'Instructions \n'
        for Key, Value in self.Instructions.items():
            Space = ' ' * (6 - len(Key))
            self.InstructionText = self.InstructionText + '  ' + Key + Space + ': ' + Value + '\n'

    def DisplayInstructions(self):
        TextLines = self.InstructionText.splitlines()
        Height = Font_Monospace.get_linesize() + 2
        X, Y = Screen.get_rect().topleft
        Color = None
        for i, TextLine in enumerate(TextLines):
            Color = BLACK
            Text = Font_Monospace.render(TextLine, True, Color)
            Text_Rect = Text.get_rect()
            Text_Rect.topleft = (X, Y + Height * i)
            Screen.blit(Text, Text_Rect)

    def DisplayPause(self):
        Text = Font_Sans.render(self.MessagePause, True, BLACK, WHITE)
        Text_Rect = Text.get_rect()
        Text_Rect.center = self.Pos_Welcome
        Screen.blit(Text, Text_Rect)


def main():
    print('\n%%%%%% Start FtsZ ring mock-visualization. %%%%%%')

    Control = FControl()
    SimUnitTime = 0.2

    # Switches
    Switch_Polymerization = True
    Switch_Depolymerization = True
    Switch_Constriction = True

    # Initialization
    EcoliBody = FCompartment(linecolor=GRAY2, bodycolor=GRAY2)
    Membrane = FCompartment()
    FtsZ = FObject('FtsZ', type='mass', quantity=1000)
    Dict_FtsA = dict()
    for i in range(8):
        Angle = np.deg2rad(45*i)
        X, Y = GetXYFromCenter((Membrane.X, Membrane.Y), Angle, Membrane.Radius)
        Dict_FtsA['FtsA' + str(i)] = FObject('FtsA', type='individual', id=i, x=X, y=Y, angle=Angle)
        # print(Angle, Dict_FtsA['FtsA' + str(i)].ID)

    ElapsedTime = 0
    PrevTime = datetime.now()

    SimState = True
    while SimState:

        if not Control.PauseSwitch:
            CurrTime = datetime.now()
            ElapsedTime += (CurrTime - PrevTime).total_seconds()
            PrevTime = CurrTime

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                SimState = False
            elif event.type == pygame.KEYDOWN:
                Control.MessageTimer = 5000
                if event.key == pygame.K_x:
                    SimState = False
                elif event.key == pygame.K_p:
                    Control.PauseSwitch = not Control.PauseSwitch

        if Control.PauseSwitch:
            Control.DisplayPause()

        while ElapsedTime >= SimUnitTime:

            ElapsedTime -= SimUnitTime
            Control.Time += 1

        # Children control
        KineticConstant = 0.05

        # Add Children
        if Switch_Polymerization:
            Raffle = np.random.uniform(0, 1, len(Dict_FtsA))
            for i, FtsA in enumerate(Dict_FtsA.values()):
                if FtsZ.CytosolicQuantity == 0:
                    break
                if Raffle[i] < KineticConstant * (FtsZ.CytosolicQuantity / FtsZ.TotalQuantity):
                    FtsA.N_Children += 1
                    if True:
                        FtsZAngle = FtsA.Angle + (-1) ** (len(FtsA.X_Children)) * np.deg2rad(0.7 * len(FtsA.X_Children))
                        X, Y = GetXYFromCenter((Membrane.X, Membrane.Y), FtsZAngle, Membrane.Radius, RatioFactor=0.97)
                        FtsA.X_Children.append(X)
                        FtsA.Y_Children.append(Y)
                        FtsA.Angle_Children.append(FtsZAngle)
                    FtsZ.CytosolicQuantity -= 1

        # Remove Children
        if Switch_Depolymerization:
            Raffle = np.random.uniform(0, 1, len(Dict_FtsA))
            for i, FtsA in enumerate(Dict_FtsA.values()):
                if Raffle[i] < KineticConstant * 0.4 * (FtsZ.TotalQuantity / FtsZ.CytosolicQuantity):
                    if FtsA.N_Children >= 1:
                        FtsA.N_Children -= 1
                        if True:
                            FtsA.X_Children.pop()
                            FtsA.Y_Children.pop()
                            FtsA.Angle_Children.pop()
                        FtsZ.CytosolicQuantity += 1

        # FtsZ Constriction
        if Switch_Constriction:
            ConstrictionRate = 0.9999 - 0.0001 * (FtsZ.CytosolicQuantity / FtsZ.TotalQuantity)
            if FtsZ.TotalQuantity - FtsZ.CytosolicQuantity > 200:
                # Membrane
                Membrane.Radius *= ConstrictionRate
                for i, FtsA in enumerate(Dict_FtsA.values()):
                    X, Y = GetXYFromCenter((Membrane.X, Membrane.Y), FtsA.Angle, Membrane.Radius)
                    FtsA.X = X
                    FtsA.Y = Y
                    # Children
                    X, Y = GetXYFromCenter((Membrane.X, Membrane.Y), FtsA.Angle_Children, Membrane.Radius, RatioFactor=0.97)
                    FtsA.X_Children = X.tolist()
                    FtsA.Y_Children = Y.tolist()


        # All the Drawings
        Screen.fill(WHITE)

        EcoliBody.Draw()
        Membrane.Draw()
        for FtsA in Dict_FtsA.values():
            FtsA.Draw()
        FtsZ.Draw()

        if Control.InstructionSwitch:
            Control.DisplayInstructions()

        Control.DisplayTime()

        pygame.display.update()

    pygame.quit()
    print('\n%%%%%% End FtsZ ring mock-visualization. %%%%%%')
    sys.exit()


if __name__ == '__main__':
    main()
