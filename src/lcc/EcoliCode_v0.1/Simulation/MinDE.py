import sys
import pygame
from datetime import datetime
import numpy as np
import math

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

def SciFloat(Float, InPrecision=4, InExp_digits=2):
    return np.format_float_scientific(Float, precision=InPrecision, unique=False, exp_digits=InExp_digits)

def DisplayString(String, X, Y, position='center', rotate=0, color=BLACK, Type='mass'):
    Text=''
    if Type == 'mass':
        Text = Font_MassObject.render(String, True, color)
    elif Type == 'individual':
        Text = Font_IndividualObject.render(String, True, color)
    elif Type == 'bold':
        Text = Font_Sans.render(String, True, BLACK)
    Text = pygame.transform.rotate(Text, rotate)
    Text_Rect = Text.get_rect()
    if position == 'center':
        Text_Rect.center = (X, Y)
    if position == 'midleft':
        Text_Rect.midleft = (X, Y)
    elif position == 'midright':
        Text_Rect.midright = (X, Y)
    Screen.blit(Text, Text_Rect)

# pygame
pygame.init()

Screen_Size = W_S, H_S = 1500, 800
Screen = pygame.display.set_mode(Screen_Size)

LEFT = 0
MID_X = W_S / 2
RIGHT = W_S

TOP = 0
MID_Y = H_S / 2
BOTTOM = H_S

CenterTop = (MID_X, TOP)
Center = (MID_X, MID_Y)

# X-section along the plane of division
Offset_X = W_S * 2 / 7
Offset_Y = 0

# X-section along the axis of division of E coli
Offset_X = - W_S * 9 / 40
Offset_Y = 0

def AddPos(A, B):
    X_A, Y_A = A
    X_B, Y_B = B
    return (X_A + X_B, Y_A + Y_B)

def GetMidPoint(A, B):
    X_A, Y_A = A
    X_B, Y_B = B
    return ((X_A + X_B) / 2, (Y_A + Y_B) / 2)

def GetDistanceBTWTwoPoints(A, B):
    X_A, Y_A = A
    X_B, Y_B = B
    return math.sqrt((X_A - X_B) ** 2 + (Y_A - Y_B) ** 2)

def GetXYFromCenter(Center, Angle, Radius, RatioFactor=1.0):
    Center_X, Center_Y = Center
    return Center_X + Radius * RatioFactor * np.cos(Angle), Center_Y + Radius * RatioFactor * np.sin(Angle)

def GetXYFromEllipse(Center, Angle, Width, Height, RatioFactor=1.0):
    Center_X, Center_Y = Center
    return Center_X + Width / 2 * RatioFactor * np.cos(Angle), Center_Y + Height / 2 * RatioFactor * np.sin(Angle)

Title = 'FtsZ Ring Visualization'
pygame.display.set_caption(Title)

Font_Init = pygame.font.SysFont('arial', 50)
Font_MassObject = pygame.font.SysFont('arial', 20)
Font_IndividualObject = pygame.font.SysFont('arial', 15)

Font_Sans = pygame.font.Font('freesansbold.ttf', 20)
Font_Monospace = pygame.font.SysFont('monospace', 18, True)
Font_Radar = pygame.font.SysFont('arial', 11)

class FCompartment:
    def __init__(self, X=MID_X, Y=MID_Y, thickness=5, linecolor=GRAY4, bodycolor=YELLOW_FAINT, Label=False, width=0, height=0):
        self.X = X
        self.Y = Y
        self.W = width
        self.H = height

        self.Label = Label
        self.Thickness = thickness
        self.LineColor = linecolor
        self.BodyColor = bodycolor

    def Draw(self):
        pass

class FMembrane(FCompartment):
    def __init__(self, X=MID_X, Y=MID_Y, thickness=5, linecolor=GRAY4, bodycolor=YELLOW_FAINT, Label=False, width=0, height=0):
        FCompartment.__init__(self, X=X, Y=Y, thickness=thickness, linecolor=linecolor, bodycolor=bodycolor, Label=Label, width=width, height=height)

    def Draw(self):
        pygame.draw.ellipse(Screen, self.BodyColor, (self.X - self.W / 2, self.Y - self.H / 2, self.W, self.H))
        pygame.draw.ellipse(Screen, self.LineColor, (self.X - self.W / 2, self.Y - self.H / 2, self.W, self.H), self.Thickness)
        # self.DisplayPercentCompletionOfDivision()

    # def DisplayPercentCompletionOfDivision(self):
    #     PercentCompletion = (1 - self.Radius * 2 / self.H)  * 100
    #     Text = "% Completion: " + "{:.1f}".format(PercentCompletion)
    #     DisplayString(Text, self.X, self.Y, color=BLACK, Type='mass')

class FObject:
    def __init__(self, Name, id=0, X=MID_X, Y=MID_Y):
        self.Name = Name
        self.ID = id
        self.X = X
        self.Y = Y

    def Draw(self):
        pass

class FMassObject(FObject):
    def __init__(self, Name, id=0, X=MID_X, Y=MID_Y, quantity=100):
        FObject.__init__(self, Name, id=id, X=X, Y=Y)
        self.InitialQuantity = quantity
        self.CurrentQuantity = quantity
        self.EngagedActivity = ""
        self.Engaged = 0

    def SetEngaged(self, EngagedActivity, Engaged):
        self.EngagedActivity = EngagedActivity
        self.Engaged = Engaged

    def Draw(self, membradius=0):
        Text = "Cytosolic " + self.Name[:-3] + ": " + str(self.CurrentQuantity)
        DisplayString(Text, self.X, self.Y + 30 * self.ID, color=BLACK, Type='mass')
    #     self.DisplayEngagedActivity()
    #
    # def DisplayEngagedActivity(self):
    #     Text = '%s engaged in %s: ' % (self.Name, self.EngagedActivity) + str(self.Engaged)
    #     DisplayString(Text, self.X, H_S * 12 / 13, color=BLACK)

class FIndividualObject(FObject):
    def __init__(self, Name, id=0, X=0, Y=0, angle=0, distance=0, color=GREEN, radius=10, thickness=2):
        FObject.__init__(self, Name, id=id, X=X, Y=Y)

        self.Angle = angle
        self.DistanceFromCenter = distance
        self.Spacing = 0.5

        # ChildrenObject
        self.X_MinD = list()
        self.Y_MinD = list()
        self.Group_MinD = list()
        self.Angle_MinD = list()
        self.Color_MinD = BLUE
        self.RatioFactor_MinD = 0.97

        self.X_MinE = list()
        self.Y_MinE = list()
        self.Angle_MinE = list()
        self.Group_MinE = list()
        self.Color_MinE = RED
        self.RatioFactor_MinE = 0.94

        # Draw
        self.Radius = radius
        self.Thickness = thickness
        self.Color = color

    def Draw(self):
        # Division Plane
        pygame.draw.circle(surface=Screen, color=self.Color, center=(self.X, self.Y), radius=self.Radius)
        Label = self.Name + "#" + str(self.ID) + "~MinD_{" + str(len(self.X_MinD)) + "}" + "~MinE_{" + str(len(self.X_MinE)) + "}"
        # DisplayString(Label, self.X+self.Radius*2.5*np.cos(self.Angle), self.Y+self.Radius*2.5*np.sin(self.Angle), color=BLACK, Type='individual')

    def DrawChildren(self):
        for i in range(len(self.X_MinD)):
            pygame.draw.circle(surface=Screen, color=self.Color_MinD, center=(self.X_MinD[i] - 5 * np.cos(self.Angle), self.Y_MinD[i] - 5 * np.sin(self.Angle)), radius=self.Radius*0.5)
            # DisplayString(self.Name + str(i), self.X_Children[i], self.Y_Children[i], color=BLACK)
        for i in range(len(self.X_MinE)):
            pygame.draw.circle(surface=Screen, color=self.Color_MinE, center=(self.X_MinE[i] - 5 * np.cos(self.Angle), self.Y_MinE[i] - 5 * np.sin(self.Angle)), radius=self.Radius*0.5)
            # DisplayString(self.Name + str(i), self.X_Children[i], self.Y_Children[i], color=BLACK)

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

    def DisplayPause(self):
        Text = Font_Sans.render(self.MessagePause, True, BLACK, WHITE)
        Text_Rect = Text.get_rect()
        Text_Rect.center = self.Pos_Welcome
        Screen.blit(Text, Text_Rect)

def main():
    print('\n%%%%%% Initialize MinDE mock-visualization. %%%%%%')
    print('Screen Size:\t', W_S, '\tby\t', H_S)

    Control = FControl()
    SimUnitTime = 0.2

    # Switches
    Switch_MinD_Recruitment = True
    Switch_MinE_Recruitment = True
    Switch_MinD_Removal = True
    Switch_MinE_Removal = True
    Switch_MinD_Spread = True
    Switch_MinE_Spread = True
    
    # Switch_MinD_Recruitment = False
    # Switch_MinE_Recruitment = False
    # Switch_MinD_Removal = False
    # Switch_MinE_Removal = False
    # Switch_MinD_Spread = False
    # Switch_MinE_Spread = False

    ScaleFactor_Dimension = 1.5
    W_Ecoli = int(2000 / ScaleFactor_Dimension)
    H_Ecoli = int(1000 / ScaleFactor_Dimension)

    # Initialization
    # EcoliBody = FCompartment(X=MID_X, Y=MID_Y, width=W_Ecoli, height=H_Ecoli, linecolor=GRAY2, bodycolor=GRAY2, Label=True)
    Membrane = FMembrane(X=MID_X, Y=MID_Y, width=W_Ecoli, height=H_Ecoli, Label=True)
    print('Ecoli Size :\t', W_Ecoli, '\tby\t', H_Ecoli)

    # Molecules
    ScaleFactor_Quantity = 50
    Quantity_MinD = int(5000 / ScaleFactor_Quantity)
    Quantity_MinE = int(4000 / ScaleFactor_Quantity)

    MinD = FMassObject('MinD_MO', id=0, quantity=Quantity_MinD)
    MinE = FMassObject('MinE_MO', id=1, quantity=Quantity_MinE)
    print('MinD Quantity: ', MinD.InitialQuantity)
    print('MinE Quantity: ', MinE.InitialQuantity)

    List_MembSite = list()
    N_MembraneSites = 100
    for i in range(N_MembraneSites):
        Angle = np.deg2rad(360 / N_MembraneSites * i)
        X, Y = GetXYFromEllipse((Membrane.X, Membrane.Y), Angle, Membrane.W, Membrane.H)
        DistanceFromCenter = GetDistanceBTWTwoPoints((Membrane.X, Membrane.Y), (X, Y))
        List_MembSite.append(FIndividualObject('Phospholipid', id=i, X=X, Y=Y, angle=Angle, distance=DistanceFromCenter))
        print('Phopholipid #', List_MembSite[i].ID, '\tX:', List_MembSite[i].X, '\tY:', List_MembSite[i].Y, '\tAngle:', List_MembSite[i].Angle, '\tDistance:', List_MembSite[i].DistanceFromCenter)

    ElapsedTime = 0
    PrevTime = datetime.now()

    Debug = False

    print('\n%%%%%% Start MinDE mock-visualization. %%%%%%')
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

        while ElapsedTime >= SimUnitTime:

            ElapsedTime -= SimUnitTime
            Control.Time += 1


        # Children control
        Factor_Global = 0.1
        Kon_MinD = 0.001 * Factor_Global
        Koff_MinD = 0.2 * Factor_Global
        Kon_MinE = 0.001 * Factor_Global
        Koff_MinE = 0.01 * Factor_Global
        Kdiff_MinD = 0.1 * Factor_Global
        Kdiff_MinE = 0.7 * Factor_Global

        Raffle_MinD_Recruitment = np.random.uniform(0, 1, len(List_MembSite))
        Raffle_MinE_Recruitment = np.random.uniform(0, 1, len(List_MembSite))
        Raffle_MinD_Removal = np.random.uniform(0, 1, len(List_MembSite))
        Raffle_MinE_Removal = np.random.uniform(0, 1, len(List_MembSite))
        Raffle_MinD_Spread = np.random.uniform(0, 1, len(List_MembSite))
        Raffle_MinE_Spread = np.random.uniform(0, 1, len(List_MembSite))
        Raffle_MinD_SpreadDirection = np.random.uniform(0, 1, len(List_MembSite))
        Raffle_MinE_SpreadDirection = np.random.uniform(0, 1, len(List_MembSite))
        for i, Memb in enumerate(List_MembSite):
            # Add MinD
            if Switch_MinD_Recruitment:
                if MinD.CurrentQuantity > 0:
                # Consider cytoplasmic MinD, local MinD, Distance from the center
                    if Raffle_MinD_Recruitment[i] < Kon_MinD * (MinD.CurrentQuantity / MinD.InitialQuantity) * (len(Memb.X_MinD) + 1) ** 2 * Memb.DistanceFromCenter ** 0.6:
                        MinDGroup = (-1) ** (len(Memb.X_MinD))
                        MinDAngle = Memb.Angle + MinDGroup * np.deg2rad(Memb.Spacing * len(Memb.X_MinD))
                        X, Y = GetXYFromEllipse((Membrane.X, Membrane.Y), MinDAngle, Membrane.W, Membrane.H, RatioFactor=Memb.RatioFactor_MinD)
                        Memb.X_MinD.append(X)
                        Memb.Y_MinD.append(Y)
                        Memb.Group_MinD.append(MinDGroup)
                        Memb.Angle_MinD.append(MinDAngle)
                        MinD.CurrentQuantity -= 1

            # Add MinE
            if Switch_MinE_Recruitment:
                if MinE.CurrentQuantity > 0 and len(Memb.X_MinD) > len(Memb.X_MinE):
                    # Consider cytoplasmic MinE, local MinD, Distance from the center
                    if Raffle_MinE_Recruitment[i] < Kon_MinE * (MinE.CurrentQuantity / MinE.InitialQuantity) * len(Memb.X_MinD) * Memb.DistanceFromCenter ** 0.2:
                        MinEGroup = (-1) ** (len(Memb.X_MinE))
                        MinEAngle = Memb.Angle + MinEGroup * np.deg2rad(Memb.Spacing * len(Memb.X_MinE))
                        X, Y = GetXYFromEllipse((Membrane.X, Membrane.Y), MinEAngle, Membrane.W, Membrane.H, RatioFactor=Memb.RatioFactor_MinE)
                        Memb.X_MinE.append(X)
                        Memb.Y_MinE.append(Y)
                        Memb.Group_MinE.append(MinEGroup)
                        Memb.Angle_MinE.append(MinEAngle)
                        MinE.CurrentQuantity -= 1

            # Remove MinD
            if Switch_MinD_Removal:
                if len(Memb.X_MinD) > 0:
                    if Raffle_MinD_Removal[i] < Koff_MinD * len(Memb.X_MinE) ** 4:
                        Memb.X_MinD.pop()
                        Memb.Y_MinD.pop()
                        Memb.Group_MinD.pop()
                        Memb.Angle_MinD.pop()
                        MinD.CurrentQuantity += 1

            # Remove MinE
            if Switch_MinE_Removal:
                if len(Memb.X_MinE) > 0:
                    if Raffle_MinE_Removal[i] < Koff_MinE:
                        Memb.X_MinE.pop()
                        Memb.Y_MinE.pop()
                        Memb.Group_MinE.pop()
                        Memb.Angle_MinE.pop()
                        MinE.CurrentQuantity += 1

            # Spread MinD
            if Switch_MinD_Spread:
                if len(Memb.X_MinD) > 0:
                    if Raffle_MinD_Spread[i] < Kdiff_MinD * len(Memb.X_MinD) ** 2:
                        # Current Location
                        Memb.X_MinD.pop()
                        Memb.Y_MinD.pop()
                        Memb.Group_MinD.pop()
                        Memb.Angle_MinD.pop()

                        # New Location
                        Direction = 1 if Raffle_MinD_SpreadDirection[i] < 0.5 else -1
                        MembToMove = List_MembSite[i + Direction] if i + Direction < len(List_MembSite) else List_MembSite[0]
                        MinDGroup = (-1) ** (len(MembToMove.X_MinD))
                        MinDAngle = MembToMove.Angle + MinDGroup * np.deg2rad(Memb.Spacing * len(MembToMove.X_MinD))
                        X, Y = GetXYFromEllipse((Membrane.X, Membrane.Y), MinDAngle, Membrane.W, Membrane.H, RatioFactor=Memb.RatioFactor_MinD)
                        MembToMove.X_MinD.append(X)
                        MembToMove.Y_MinD.append(Y)
                        MembToMove.Group_MinD.append(MinDGroup)
                        MembToMove.Angle_MinD.append(MinDAngle)

            # Spread MinE
            if Switch_MinE_Spread:
                if len(Memb.X_MinE) > 0:
                    if Raffle_MinE_Spread[i] < Kdiff_MinE * len(Memb.X_MinE) ** 2:
                        # Current Location
                        Memb.X_MinE.pop()
                        Memb.Y_MinE.pop()
                        Memb.Group_MinE.pop()
                        Memb.Angle_MinE.pop()

                        # New Location
                        Direction = 1 if Raffle_MinE_SpreadDirection[i] < 0.5 else -1
                        NewPosition = i + Direction if i + Direction < len(List_MembSite) else 0
                        NewPosition = i - Direction if len(List_MembSite[NewPosition].X_MinD) == 0 else i + Direction
                        NewPosition = NewPosition if NewPosition < len(List_MembSite) else 0
                        MembToMove = List_MembSite[NewPosition]
                        MinEGroup = (-1) ** (len(MembToMove.X_MinE))
                        MinEAngle = MembToMove.Angle + MinEGroup * np.deg2rad(Memb.Spacing * len(MembToMove.X_MinE))
                        X, Y = GetXYFromEllipse((Membrane.X, Membrane.Y), MinEAngle, Membrane.W, Membrane.H, RatioFactor=Memb.RatioFactor_MinE)
                        MembToMove.X_MinE.append(X)
                        MembToMove.Y_MinE.append(Y)
                        MembToMove.Group_MinE.append(MinEGroup)
                        MembToMove.Angle_MinE.append(MinEAngle)

        # All the Drawings
        Screen.fill(WHITE)

        # EcoliBody.Draw()
        Membrane.Draw()
        for Memb in List_MembSite:
            Memb.Draw()
            Memb.DrawChildren()
        MinD.Draw()
        MinE.Draw()

        Control.DisplayTime()

        pygame.display.update()

    pygame.quit()
    print('\n%%%%%% End MinDE mock-visualization. %%%%%%')
    sys.exit()


if __name__ == '__main__':
    main()