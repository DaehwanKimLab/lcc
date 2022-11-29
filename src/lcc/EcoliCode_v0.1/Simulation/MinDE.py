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

# Utilities
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

def GetDistanceBTWTwoPoints_2D(A, B):
    X_A, Y_A = A
    X_B, Y_B = B
    return np.sqrt((X_A - X_B.T) ** 2 + (Y_A - Y_B.T) ** 2)

def GetXYFromCenter(Center, Angle, Radius, RatioFactor=1.0):
    Center_X, Center_Y = Center
    return Center_X + Radius * RatioFactor * np.cos(Angle), Center_Y + Radius * RatioFactor * np.sin(Angle)

def GetXYFromEllipse(Center, Angle, Width, Height, RatioFactor=1.0):
    Center_X, Center_Y = Center
    return Center_X + Width / 2 * RatioFactor * np.cos(Angle), Center_Y + Height / 2 * RatioFactor * np.sin(Angle)

def GetRandomXYWithinEllipse(Center, Width, Height):
    Center_X, Center_Y = Center
    X = np.random.uniform(Center_X - Width / 2, Center_X + Width / 2)
    Y = np.random.uniform(Center_Y - Height / 2, Center_Y + Height / 2)
    if CheckIfWithinEllipse(X, Y, Center_X, Center_Y, Width, Height):
        pass
    else:
        X, Y = GetRandomXYWithinEllipse(Center, Width, Height)
    return X, Y

def CheckIfWithinEllipse(X, Y, Center_X, Center_Y, Width, Height):
    if (X - Center_X) ** 2 / (Width / 2) ** 2 + (Y - Center_Y) ** 2 / (Height / 2) ** 2 < 1:
        return True
    else:
        return False

Title = 'MinDE Visualization'
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
        self.Boundary = X, Y, width, height

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
        self.Step = 0

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

    def PrintStep(self):
        Text = 'Simulation Step: ' + str(round(self.Step))
        print(Text)

    def PrintTime(self):
        Text = 'Simulation Time: ' + str(round(self.Time))
        print(Text)

    def DisplayPause(self):
        Text = Font_Sans.render(self.MessagePause, True, BLACK, WHITE)
        Text_Rect = Text.get_rect()
        Text_Rect.center = self.Pos_Welcome
        Screen.blit(Text, Text_Rect)

def main():
    print('\n%%%%%% Initialize MinDE mock-visualization. %%%%%%')
    print('Screen Size:\t', W_S, '\tby\t', H_S)

    Control = FControl()
    SimUnitTime = 0.1
    Debug = True

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

    # Diffusion Parameters
    Distance_Interaction = 10
    Distance_Diffusion = 5

    # Dissociation Parameters
    KD_MinD_Phospholipid            = 0.001
    KD_MinD_Phospholipid_And_MinD   = 0.00001
    KD_MinD_Phospholipid_And_MinE   = 0.9
    KD_MinE_Phospholipid_And_MinD   = 0.0001

    ScaleFactor_Dimension = 1.5
    W_Ecoli = int(2000 / ScaleFactor_Dimension)
    H_Ecoli = int(1000 / ScaleFactor_Dimension)

    # Initialization
    # EcoliBody = FCompartment(X=MID_X, Y=MID_Y, width=W_Ecoli, height=H_Ecoli, linecolor=GRAY2, bodycolor=GRAY2, Label=True)
    Membrane = FMembrane(X=MID_X, Y=MID_Y, width=W_Ecoli, height=H_Ecoli, Label=True)
    print('Ecoli Size :\t', W_Ecoli, '\tby\t', H_Ecoli)

    # Phospholipids
    N_MembraneSites = 200
    Phospholipid_X = np.zeros([1, N_MembraneSites])
    Phospholipid_Y = np.zeros([1, N_MembraneSites])
    for i in range(Phospholipid_X.shape[1]):
        Angle = np.deg2rad(360 / N_MembraneSites * i)
        Phospholipid_X[0, i], Phospholipid_Y[0, i] = GetXYFromEllipse((Membrane.X, Membrane.Y), Angle, Membrane.W, Membrane.H)
        if Debug:
            print('Phopholipid #', i, '\tX:', Phospholipid_X[0, i], '\tY:', Phospholipid_Y[0, i])
    Phospholipid = (Phospholipid_X, Phospholipid_Y)

    # Molecules
    ScaleFactor_Quantity = 10
    Quantity_MinD = int(5000 / ScaleFactor_Quantity)
    Quantity_MinE = int(4000 / ScaleFactor_Quantity)

    MinD_X = np.zeros([1, Quantity_MinD])
    MinD_Y = np.zeros([1, Quantity_MinD])
    for i in range(MinD_X.shape[1]):
        MinD_X[0, i], MinD_Y[0, i] = GetRandomXYWithinEllipse((Membrane.X, Membrane.Y), Membrane.W, Membrane.H)
        if Debug:
            print('MinD #', i, '\tX:', MinD_X[0, i], '\tY:', MinD_Y[0, i])
    MinD = (MinD_X, MinD_Y)

    MinE_X = np.zeros([1, Quantity_MinE])
    MinE_Y = np.zeros([1, Quantity_MinE])
    for i in range(MinE_X.shape[1]):
        MinE_X[0, i], MinE_Y[0, i] = GetRandomXYWithinEllipse((Membrane.X, Membrane.Y), Membrane.W, Membrane.H)
        if Debug:
            print('MinE #', i, '\tX:', MinE_X[0, i], '\tY:', MinE_Y[0, i])
    MinE = (MinE_X, MinE_Y)

    ElapsedTime = 0
    PrevTime = datetime.now()

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

        if Debug:
            Control.PrintStep()

        # KD_MinD_Phospholipid = 0.001
        # KD_MinD_Phospholipid_And_MinD = 0.00001
        # KD_MinD_Phospholipid_And_MinE = 0.9
        # KD_MinE_Phospholipid_And_MinD = 0.0001

        def OutOfBoundaryAdjustment(X, Y, Center_X, Center_Y, Width, Height, adjustmentfactor=5):
            OutOfBoundaryCheck = (X - Center_X) ** 2 / (Width / 2) ** 2 + (Y - Center_Y) ** 2 / (Height / 2) ** 2
            bOutOfBoundaryCheck = OutOfBoundaryCheck > 1
            Idx_OutOfBoundary = np.where(bOutOfBoundaryCheck)
            if np.any(bOutOfBoundaryCheck):
                X[Idx_OutOfBoundary] = np.where(X[Idx_OutOfBoundary] < Center_X,
                                                X[Idx_OutOfBoundary] + adjustmentfactor,
                                                X[Idx_OutOfBoundary] - adjustmentfactor)
                Y[Idx_OutOfBoundary] = np.where(Y[Idx_OutOfBoundary] < Center_Y,
                                                Y[Idx_OutOfBoundary] + adjustmentfactor,
                                                Y[Idx_OutOfBoundary] - adjustmentfactor)
            return X, Y

        def RandomMovement(X, movementfactor=5):
            X += np.random.uniform(-1, 1, X.shape[0]) * movementfactor
            return X

        def Diffusion(XY, Idx, movementfactor, Boundary):
            X, Y = XY
            Center_X, Center_Y, Width, Height = Boundary
            X[:, Idx] += np.random.uniform(-1, 1, Idx.shape[0]) * movementfactor
            Y[:, Idx] += np.random.uniform(-1, 1, Idx.shape[0]) * movementfactor
            X[:, Idx], Y[:, Idx] = OutOfBoundaryAdjustment(X[:, Idx], Y[:, Idx], Center_X, Center_Y, Width, Height)
            return (X, Y)

        def GetDistance(xy_mol1, xy_mol2):
            return GetDistanceBTWTwoPoints_2D(xy_mol1, xy_mol2)

        def GetInteractionMap(xy_mol1=None, xy_mol2=None, distance=None, interaction_distance=0, select_mol=0):
            '''
            Possible set of inputs:
                - xy_mol1, xy_mol2, interaction_distance, select_mol
                - distance, interaction_distance, select_mol
            '''
            if np.any(xy_mol1) and np.any(xy_mol2) and interaction_distance:
                distance = GetDistance(xy_mol1, xy_mol2)
            elif np.any(distance) and interaction_distance:
                pass
            else:
                print('[ERROR] GetInteractionMap: insufficient inputs')
                sys.exit(1)
            return np.any(distance < interaction_distance, axis=select_mol)

        def GetInteractionIdx(xy_mol1=None, xy_mol2=None, distance=None, interaction_distance=0, interaction_map=None, select_mol=0):
            '''
            Possible set of inputs:
                - xy_mol1, xy_mol2, interaction_distance, select_mol=0
                - distance, interaction_distance, select_mol=0
                - interaction_map, select_mol=0
            select_mol: 0 for mol1, 1 for mol2
            '''
            if np.any(xy_mol1) and np.any(xy_mol2) and interaction_distance:
                interaction_map = GetInteractionMap(xy_mol1=xy_mol1, xy_mol2=xy_mol2, interaction_distance=interaction_distance)
            elif np.any(distance) and interaction_distance:
                interaction_map = GetInteractionMap(distance=distance, interaction_distance=interaction_distance)
            elif interaction_map.size:
                pass
            else:
                print('[ERROR] GetInteractionIdx: insufficient inputs')
                sys.exit(1)
            return np.where(interaction_map)[select_mol]

        # MinD-Phospholipid interaction
        Distance_MinD_To_Phospholipid = GetDistance(MinD, Phospholipid)
        InteractionMap_MinD_With_Phospholipid = GetInteractionMap(distance=Distance_MinD_To_Phospholipid, interaction_distance=Distance_Interaction, select_mol=0)
        InteractionMap_MinD_Not_With_Phospholipid = np.invert(InteractionMap_MinD_With_Phospholipid)
        Idx_MinD_Not_With_Phospholipid = GetInteractionIdx(interaction_map=InteractionMap_MinD_Not_With_Phospholipid, select_mol=0)

        MinD = Diffusion(MinD, Idx_MinD_Not_With_Phospholipid, Distance_Diffusion, Membrane.Boundary)

        # MinE-MinD-Phospholipid interaction
        Idx_MinD_With_Phospholipid = GetInteractionIdx(interaction_map=InteractionMap_MinD_With_Phospholipid, select_mol=0)
        MinD_With_Phospholipid = MinD[0][:, Idx_MinD_With_Phospholipid], MinD[1][:, Idx_MinD_With_Phospholipid]
        InteractionMap_MinE_With_MinDWithPhospholipid = GetInteractionMap(xy_mol1=MinE, xy_mol2=MinD_With_Phospholipid, interaction_distance=Distance_Interaction, select_mol=0)
        InteractionMap_MinE_With_Phospholipid = GetInteractionMap(xy_mol1=MinE, xy_mol2=Phospholipid, interaction_distance=Distance_Interaction, select_mol=0)
        InteractionMap_MinE_With_MinDWithPhospholipid_And_Phospholipid = np.logical_and(InteractionMap_MinE_With_MinDWithPhospholipid, InteractionMap_MinE_With_Phospholipid)
        InteractionMap_MinE_NotWith_MinDWithPhospholipid_And_Phospholipid = np.invert(InteractionMap_MinE_With_MinDWithPhospholipid_And_Phospholipid, InteractionMap_MinE_With_Phospholipid)
        Idx_MinE_With_MinDWithPhospholipid_And_Phospholipid = GetInteractionIdx(interaction_map=InteractionMap_MinE_NotWith_MinDWithPhospholipid_And_Phospholipid, select_mol=0)

        MinE = Diffusion(MinE, Idx_MinE_With_MinDWithPhospholipid_And_Phospholipid, Distance_Diffusion, Membrane.Boundary)

        # Phospholipid-MinD dissociation


        # All the Drawings
        Screen.fill(WHITE)
        Membrane.Draw()
        for i in range(Phospholipid_X.shape[1]):
            pygame.draw.circle(surface=Screen, color=GREEN, center=(Phospholipid_X[0, i], Phospholipid_Y[0, i]), radius=5)
        for i in range(MinD_X.shape[1]):
            pygame.draw.circle(surface=Screen, color=BLUE, center=(MinD_X[0, i], MinD_Y[0, i]), radius=5)
        for i in range(MinE_X.shape[1]):
            pygame.draw.circle(surface=Screen, color=RED, center=(MinE_X[0, i], MinE_Y[0, i]), radius=5)

        Control.DisplayTime()
        pygame.display.update()

        Control.Step += 1

    pygame.quit()
    print('\n%%%%%% End MinDE mock-visualization. %%%%%%')
    sys.exit()


if __name__ == '__main__':
    main()
