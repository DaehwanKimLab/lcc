import sys
import pygame
from datetime import datetime
import numpy as np

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

Screen_Size = W_S, H_S = 1700, 800
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
Offset_DivPlane_X = W_S * 2 / 7
Offset_DivPlane_Y = 0

# X-section along the axis of division of E coli
Offset_DivAxis_X = - W_S * 1 / 5
Offset_DivAxis_Y = 0

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
    def __init__(self, X=MID_X, Y=MID_Y, Radius=H_S*7/10/2, Thickness=5, linecolor=GRAY4, bodycolor=YELLOW_FAINT, Label=False):
        self.X = X
        self.Y = Y
        self.X_DivPlane = X + Offset_DivPlane_X
        self.Y_DivPlane = Y + Offset_DivPlane_Y
        self.X_DivAxis = X + Offset_DivAxis_X
        self.Y_DivAxis = Y + Offset_DivAxis_Y
        self.W_DivAxis = Radius * 3
        self.H_DivAxis = Radius * 2

        self.Label = Label
        self.Radius = int(Radius)
        self.Thickness = Thickness
        self.LineColor = linecolor
        self.BodyColor = bodycolor

    def Draw(self):
        # Division Plane
        pygame.draw.circle(Screen, self.BodyColor, (self.X_DivPlane, self.Y_DivPlane), self.Radius)

        # Labelking
        if self.Label:
            DisplayString("X-section along the Axis Of Division", self.X_DivAxis, H_S / 15, color=BLACK, Type='bold')
            DisplayString("X-section along the Plane Of Division", self.X_DivPlane, H_S / 15, color=BLACK, Type='bold')

class FMembrane(FCompartment):
    def __init__(self, X=MID_X, Y=MID_Y, Radius=H_S*7/10/2, Thickness=5, linecolor=GRAY4, bodycolor=YELLOW_FAINT, Label=False):
        FCompartment.__init__(self, X=X, Y=Y, Radius=Radius, Thickness=Thickness, linecolor=linecolor, bodycolor=bodycolor, Label=Label)

    def Draw(self):
        # Division Plane
        pygame.draw.circle(Screen, self.BodyColor, (self.X_DivPlane, self.Y_DivPlane), self.Radius)
        pygame.draw.circle(Screen, self.LineColor, (self.X_DivPlane, self.Y_DivPlane), self.Radius, self.Thickness)

        # Division Axis: Cell Shape
        pygame.draw.rect(Screen, self.BodyColor, (self.X_DivAxis - self.W_DivAxis / 2, self.Y_DivAxis - self.H_DivAxis / 2, self.W_DivAxis, self.H_DivAxis))
        pygame.draw.rect(Screen, self.LineColor, (self.X_DivAxis - self.W_DivAxis / 2, self.Y_DivAxis - self.H_DivAxis / 2, self.W_DivAxis, self.H_DivAxis), self.Thickness + 1)

        # Division Axis: Septum
        EdgeHandling = 1
        pygame.draw.line(Screen, self.LineColor, (self.X_DivAxis, self.Y_DivAxis - self.H_DivAxis / 2 + EdgeHandling), (self.X_DivAxis, self.Y_DivAxis + self.H_DivAxis / 2 - EdgeHandling), self.Thickness * 2)
        pygame.draw.line(Screen, self.BodyColor, (self.X_DivAxis, self.Y_DivAxis - self.H_DivAxis / 2), (self.X_DivAxis, self.Y_DivAxis + self.H_DivAxis / 2), 1)
        pygame.draw.line(Screen, self.BodyColor, (self.X_DivAxis, self.Y_DivAxis - self.Radius + self.Thickness + EdgeHandling), (self.X_DivAxis, self.Y_DivAxis + self.Radius - self.Thickness - EdgeHandling), self.Thickness * 2)

        self.DisplayPercentCompletionOfDivision()

    def DisplayPercentCompletionOfDivision(self):
        PercentCompletion = (1 - self.Radius * 2 / self.H_DivAxis)  * 100
        Text = "% Completion: " + "{:.1f}".format(PercentCompletion)
        DisplayString(Text, self.X_DivAxis, self.Y_DivAxis, color=BLACK, Type='mass')

class FObject:
    def __init__(self, Name, id=0, X=MID_X, Y=MID_Y):
        self.Name = Name
        self.ID = id
        self.X_DivPlane = X + Offset_DivPlane_X
        self.Y_DivPlane = Y + Offset_DivPlane_Y
        self.X_DivAxis = X + Offset_DivAxis_X
        self.Y_DivAxis = Y + Offset_DivAxis_Y

    def Draw(self):
        pass

class FMassObject(FObject):
    def __init__(self, Name, id=0, X=MID_X, Y=MID_Y, quantity=100):
        FObject.__init__(self, Name, id=id, X=X, Y=Y)
        self.TotalQuantity = quantity
        self.CytosolicQuantity = quantity
        self.EngagedActivity = ""
        self.Engaged = 0

    def SetEngaged(self, EngagedActivity, Engaged):
        self.EngagedActivity = EngagedActivity
        self.Engaged = Engaged

    def Draw(self, membradius=0):
        Text = "Cytosolic " + self.Name + ": " + str(self.CytosolicQuantity)
        DisplayString(Text, self.X_DivPlane, self.Y_DivPlane, color=BLACK, Type='mass')
        self.DisplayEngagedActivity()

    def DisplayEngagedActivity(self):
        Text = '%s engaged in %s: ' % (self.Name, self.EngagedActivity) + str(self.Engaged)
        DisplayString(Text, self.X_DivPlane, H_S * 12 / 13, color=BLACK)

class FIndividualObject(FObject):
    def __init__(self, Name, id=0, X=MID_X, Y=MID_Y, angle=0, color=BLACK, radius=10, quantity=100, thickness=2):
        FObject.__init__(self, Name, id=id, X=X, Y=Y)

        self.Angle = angle
        self.TotalQuantity = quantity
        self.CytosolicQuantity = quantity

        # ChildrenObject
        self.N_Children = 0
        self.X_Children_DivPlane = list()
        self.Y_Children_DivPlane = list()
        self.X_Children_DivAxis = list()
        self.Y_Children_DivAxis = list()
        self.Group_Children = list()
        self.Angle_Children = list()

        # Draw
        self.Radius = radius
        self.Thickness = thickness
        self.Color = color

    def Draw(self):
        # Division Plane
        pygame.draw.circle(surface=Screen, color=GREEN, center=(self.X_DivPlane, self.Y_DivPlane), radius=self.Radius)
        Label = self.Name + "#" + str(self.ID) + "~FtsZ_{" + str(self.N_Children) + "}"
        DisplayString(Label, self.X_DivPlane+self.Radius*2.5*np.cos(self.Angle), self.Y_DivPlane+self.Radius*2.5*np.sin(self.Angle), color=BLACK, Type='individual')

        # Division Axis
        Angle = np.rad2deg(self.Angle)
        if Angle == 90 or Angle == 270:
            pygame.draw.circle(surface=Screen, color=GREEN, center=(self.X_DivAxis - 5 * np.cos(self.Angle), self.Y_DivAxis - 5 * np.sin(self.Angle)), radius=self.Radius)

    def DrawChildren(self):
        # Division Plane
        for i in range(len(self.X_Children_DivPlane)):
            pygame.draw.circle(surface=Screen, color=RED, center=(self.X_Children_DivPlane[i] - 5 * np.cos(self.Angle), self.Y_Children_DivPlane[i] - 5 * np.sin(self.Angle)), radius=self.Radius*0.5, width=1)
            # DisplayString(self.Name + str(i), self.X_Children[i], self.Y_Children[i], color=BLACK)

        # Division Axis
        Angle = np.rad2deg(self.Angle)
        if Angle == 90 or Angle == 270:
            for i in range(len(self.X_Children_DivPlane)):
                pygame.draw.circle(surface=Screen, color=RED, center=(self.X_Children_DivAxis[i], self.Y_Children_DivAxis[i] - 5 * np.sin(self.Angle)), radius=self.Radius*0.5, width=1)
                if i > 4:
                    break
        else:
            pass

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
    print('\n%%%%%% Start FtsZ ring mock-visualization. %%%%%%')

    Control = FControl()
    SimUnitTime = 0.2

    # Switches
    Switch_Polymerization = True
    Switch_Depolymerization = True
    Switch_Constriction = True

    # Initialization
    EcoliBody = FCompartment(X=MID_X, Y=MID_Y, linecolor=GRAY2, bodycolor=GRAY2, Label=True)
    Membrane = FMembrane(X=MID_X, Y=MID_Y)
    FtsZ = FMassObject('FtsZ', quantity=1000)
    List_FtsA = list()
    N_FtsA = 20
    for i in range(N_FtsA):
        Angle = np.deg2rad(360 / N_FtsA * i)
        X, Y = GetXYFromCenter((Membrane.X, Membrane.Y), Angle, Membrane.Radius)
        List_FtsA.append(FIndividualObject('FtsA', id=i, X=X, Y=Y, angle=Angle))
        # print(Angle, List_FtsA['FtsA' + str(i)].ID)
    Treadmilling = 0

    ElapsedTime = 0
    PrevTime = datetime.now()

    Debug = False

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
        KineticConstant = 0.001 * N_FtsA

        # Add Children
        if Switch_Polymerization:
            Raffle = np.random.uniform(0, 1, len(List_FtsA))
            for i, FtsA in enumerate(List_FtsA):
                if FtsZ.CytosolicQuantity == 0:
                    break
                if Raffle[i] < KineticConstant * (FtsZ.CytosolicQuantity / FtsZ.TotalQuantity):
                    FtsA.N_Children += 1
                    if True:
                        FtsZGroup = (-1) ** (len(FtsA.X_Children_DivPlane))
                        FtsZAngle = FtsA.Angle + FtsZGroup * np.deg2rad(0.75 * len(FtsA.X_Children_DivPlane))
                        # if FtsZAngle < 0:
                        #     FtsZAngle += np.deg2rad(360)
                        FtsZRatioFactor = 0.92 + 0.02 * FtsZGroup
                        X, Y = GetXYFromCenter((Membrane.X, Membrane.Y), FtsZAngle, Membrane.Radius, RatioFactor=FtsZRatioFactor)
                        FtsA.X_Children_DivPlane.append(X + Offset_DivPlane_X)
                        FtsA.Y_Children_DivPlane.append(Y + Offset_DivPlane_Y)
                        FtsA.X_Children_DivAxis.append(X + Offset_DivAxis_X)
                        FtsA.Y_Children_DivAxis.append(Y + Offset_DivAxis_Y)
                        FtsA.Group_Children.append(FtsZGroup)
                        FtsA.Angle_Children.append(FtsZAngle)
                    FtsZ.CytosolicQuantity -= 1

        # Remove Children
        if Switch_Depolymerization:
            Raffle = np.random.uniform(0, 1, len(List_FtsA))
            for i, FtsA in enumerate(List_FtsA):
                if Raffle[i] < KineticConstant * 0.2 * (FtsZ.TotalQuantity / FtsZ.CytosolicQuantity):
                    if FtsA.N_Children >= 1:
                        FtsA.N_Children -= 1
                        if True:
                            FtsA.X_Children_DivPlane.pop()
                            FtsA.Y_Children_DivPlane.pop()
                            FtsA.X_Children_DivAxis.pop()
                            FtsA.Y_Children_DivAxis.pop()
                            FtsA.Group_Children.pop()
                            FtsA.Angle_Children.pop()
                        FtsZ.CytosolicQuantity += 1


        # FtsZ Constriction
        if Switch_Constriction:
            N_GlobalTreadmilling = 0
            N_LocalTreadmilling = 0
            for i, FtsA in enumerate(List_FtsA):
                TreadmillingArray_Prev = np.rad2deg(List_FtsA[i - 1].Angle_Children)
                TreadmillingArray = np.rad2deg(FtsA.Angle_Children)
                if TreadmillingArray_Prev.shape[0] > 0 and TreadmillingArray.shape[0] > 0:
                    Max_Now = np.max(TreadmillingArray)
                    Min_Prev = np.min(TreadmillingArray_Prev)
                    if np.all(Max_Now < Min_Prev):
                        if Debug:
                            print("** CORRECTION ** FtsA #:", FtsA.ID, "-------\nPrev:", TreadmillingArray_Prev, "\nNow:", TreadmillingArray)
                        TreadmillingArray_Prev -= 360
                    Max_Prev = np.max(TreadmillingArray_Prev)
                    N_LocalTreadmilling = np.count_nonzero(TreadmillingArray < Max_Prev)
                    N_GlobalTreadmilling += N_LocalTreadmilling

                if Debug:
                    print("------- FtsA #:", FtsA.ID, "-------\nPrev:", TreadmillingArray_Prev, "\nNow:", TreadmillingArray, "\nLocalTreadmilling:", N_LocalTreadmilling)
                TreadmillingArray_Prev = TreadmillingArray


            if Debug:
                print("================ GlobalTreadmilling:", N_GlobalTreadmilling, " ================")
            FtsZ.SetEngaged("Treadmilling", N_GlobalTreadmilling)

            ConstrictionRate = 0.000001 * N_GlobalTreadmilling
            if ConstrictionRate > 0:
                # Membrane
                Membrane.Radius *= (1 - ConstrictionRate)
                for i, FtsA in enumerate(List_FtsA):
                    # FtsA
                    X, Y = GetXYFromCenter((Membrane.X, Membrane.Y), FtsA.Angle, Membrane.Radius)
                    FtsA.X_DivPlane = X + Offset_DivPlane_X
                    FtsA.Y_DivPlane = Y + Offset_DivPlane_Y
                    FtsA.X_DivAxis = X + Offset_DivAxis_X
                    FtsA.Y_DivAxis = Y + Offset_DivAxis_Y

                    # FtsA Children
                    FtsZRatioFactor = 0.92 + 0.02 * np.array(FtsA.Group_Children)
                    X, Y = GetXYFromCenter((Membrane.X, Membrane.Y), np.array(FtsA.Angle_Children), Membrane.Radius, RatioFactor=FtsZRatioFactor)
                    FtsA.X_Children_DivPlane = (X + Offset_DivPlane_X).tolist()
                    FtsA.Y_Children_DivPlane = (Y + Offset_DivPlane_Y).tolist()
                    FtsA.X_Children_DivAxis = (X + Offset_DivAxis_X).tolist()
                    FtsA.Y_Children_DivAxis = (Y + Offset_DivAxis_Y).tolist()


        # All the Drawings
        Screen.fill(WHITE)

        EcoliBody.Draw()
        Membrane.Draw()
        for FtsA in List_FtsA:
            FtsA.Draw()
            FtsA.DrawChildren()
        FtsZ.Draw(Membrane.Radius)

        Control.DisplayTime()

        pygame.display.update()

    pygame.quit()
    print('\n%%%%%% End FtsZ ring mock-visualization. %%%%%%')
    sys.exit()


if __name__ == '__main__':
    main()
