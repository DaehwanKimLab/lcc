import sys
import pygame
import math
import random
from datetime import datetime
import Ingalls2013_Model6_13_BacterialChemotaxis_SimForVis

# Good for random walk demonstrations

# Set 1
ThresholdFactor = 0.9999999999
# random.seed(6)
# random.seed(10)
random.seed(12)
# random.seed(14)

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
pi = math.pi
uM = 1e-6
nM = 1e-9
UnitTxt = ''

# Determine global unit
Unit = nM
if Unit == nM:
    UnitTxt = 'nM'
elif Unit == uM:
    UnitTxt = 'uM'

# Utilities
def GetColorGradient(Fade, baseColor=None):
    assert Fade <= 255, 'ERROR: Fade factor is higher than 255!'
    if baseColor == 'Blue':
        return (Fade, Fade, 255)
    else:  # Assumes White
        return (Fade, Fade, Fade)

def Displacement(Distance, Angle):
    dX = Distance * math.sin(Angle)
    dY = Distance * math.cos(Angle)
    return dX, dY

# pygame
pygame.init()

Screen_Size = W_S, H_S = 1200, 800
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

# # Transparent control board
# ControlBoard = pygame.Surface(Screen_Size, pygame.SRCALPHA)
# ControlBoard.fill((0, 0, 0, 255))

Title = "Vis2D"
pygame.display.set_caption(Title)

Font_Sans = pygame.font.Font('freesansbold.ttf', 20)
Font_Monospace = pygame.font.SysFont('monospace', 15, True)

# Initialize model
Model = Ingalls2013_Model6_13_BacterialChemotaxis_SimForVis.FModel()


class FEnvironment:
    def __init__(self, InX=W_S/2, InY=H_S/2, InShape='circle', InRadius=W_S*0.3, InThickness=5):
        self.X = InX
        self.Y = InY
        self.Shape = InShape
        self.Radius = int(InRadius)
        self.Thickness = InThickness
        self.TransparentCircleArea = None

    def Draw(self, shape='circle'):
        if shape == 'circle':
            pygame.draw.circle(Screen, YELLOW_FAINT, (self.X, self.Y), self.Radius)
            pygame.draw.circle(Screen, GRAY4, (self.X, self.Y), self.Radius, self.Thickness)

    def DrawTransparentArea(self):
        self.TransparentCircleArea = pygame.Surface((self.Radius * 2, self.Radius * 2), pygame.SRCALPHA)
        self.TransparentCircleArea.fill((255, 255, 255, 255))
        pygame.draw.circle(self.TransparentCircleArea, (0, 0, 0, 0), (self.Radius, self.Radius), self.Radius)

    def CheckOutOfBound(self, X, Y):
        pass

class FOrganism:
    def __init__(self, InSpecies, InX, InY, InAngle=pi/2, InSpeedMax=5):
        self.Species = InSpecies
        self.X_Ori = InX
        self.Y_Ori = InY
        self.Angle = InAngle
        self.X = self.X_Ori
        self.Y = self.Y_Ori
        self.X_Prev = None
        self.Y_Prev = None
        self.Speed = InSpeedMax
        self.Speed_Max = InSpeedMax
        self.TumbleAngle = -pi/4
        self.FlagellaLength = 10

        # Memory
        self.Glucose_Prev = 0
        self.SimCount = 0

        # Trajectory
        self.TrajectorySwitch = True
        self.Trajectory = [(self.X_Ori, self.Y_Ori)]

        # Image
        self.Image = None
        self.Image_ScaleFactor = 50
        self.Image_Size_X = 500 / 550 * self.Image_ScaleFactor
        self.Image_Size_Y = 556 / 556 * self.Image_ScaleFactor
        self.Image_Rect = None
        self.LoadImage()

        # Mechanistic Mode
        self.Am = 0
        self.AmThreshold = 0
        # self.AmThreshold = 1.165

        # DK - debugging purposes
        self.MechanisticModeSwitch = True
        # self.MechanisticModeSwitch = False

    def LoadImage(self):
        if self.Species == 'Ecoli':
            self.Image = pygame.image.load('../../pygamelib/ecoli.png')
            self.Image = pygame.transform.scale(self.Image, (self.Image_Size_X, self.Image_Size_Y))
            self.Image_Rect = self.Image.get_rect()
            self.RotateImage(220)

    def RotateImage(self, AngleInDegree):
        self.Image = pygame.transform.rotate(self.Image, AngleInDegree)

    def CenterImage(self):
        # Rect = self.Image.get_rect()
        # self.Image.center = Rect.center
        pass

    def Reinitialize(self):
        self.X = self.X_Ori
        self.Y = self.Y_Ori
        self.Trajectory = [(self.X_Ori, self.Y_Ori)]

    def Draw(self):
        if self.Species == 'Ecoli':
            # BlitRotate2(Screen, self.Image, self.Image.topleft, math.degrees(self.TumbleAngle))
            # Image = pygame.surface.Surface((100, 100))
            Screen.blit(self.Image, (self.X - self.Image_Size_X / 2, self.Y - self.Image_Size_Y / 2))
        else:
            Color = YELLOW
            # pygame.draw.circle(Screen, Color, (self.X, self.Y), 5)
            dX =  math.cos(self.Angle) * -20
            dY = -math.sin(self.Angle) * -20
            pygame.draw.line(Screen, Color, (self.X, self.Y), (self.X + dX, self.Y + dY), 10)
            pygame.draw.line(Screen, Color, (self.X, self.Y), (self.X + 2 * dX, self.Y + 2 * dY), 3)

    def Chemotaxis(self, GlucoseLvl):
        if self.MechanisticModeSwitch:
            # At homeostasis, Am: 1.1653948157952327e-09

            # Perform 100 simulations
            for _ in range(100):
                self.SimCount += 1
                self.Am = Model.Simulate(GlucoseLvl)

            Delta = (GlucoseLvl - self.Glucose_Prev) / GlucoseLvl * 100
            print("[Chemotaxis  {:06d}] Glucose:{:.6f} {} ({}{:.4f}%) Am:{:.6f} {} (X:{:.2f} Y:{:.2f} {:3.1f} degree)".format
                  (self.SimCount, GlucoseLvl / Unit, UnitTxt, ("+" if Delta >= 0 else ""), Delta, self.Am / Unit, UnitTxt, self.X, self.Y, self.Angle / pi * 180))
            # print(self.Am)
            if self.Am < self.AmThreshold:
            # if self.Am < 1.165 * uM:
                self.Move(self.Speed)
            else:
                self.Tumble()

        else:
            if GlucoseLvl > self.Glucose_Prev:
                self.Move(self.Speed)
            else:
                self.Tumble()
                
        # Update
        self.Glucose_Prev = GlucoseLvl

    def Homeostasis(self, GlucoseLvl):
        while True:
            self.SimCount += 1
            PrevAm = Model.Simulate(GlucoseLvl)
            if PrevAm > 0 and abs(self.Am - PrevAm) / PrevAm < 1e-7:
                break
            self.Am = PrevAm
            self.AmThreshold = PrevAm * ThresholdFactor
            # self.AmThreshold = PrevAm * 0.999 # good for physiological condition
            print("[Homeostasis {:06d}] Glucose:{:.6f}{} Am:{:.6f}{}".format(self.SimCount, GlucoseLvl / Unit, UnitTxt, self.Am / Unit, UnitTxt))

        self.Glucose_Prev = GlucoseLvl


    def Move(self, Distance, Angle=None):
        if not Angle:
            Angle = self.Angle
        dX, dY = self.Displacement(Distance, Angle)
        self.X += dX
        self.Y += dY

    def Tumble(self):
        self.Angle += random.random() * pi
        # self.Angle += self.TumbleAngle
        while self.Angle < 0:
            self.Angle += 2 * pi

        # self.Move(0)
        self.Move(self.Speed / 3)
        # self.RandomMovement()
        self.Trajectory.append((self.X, self.Y))
        if self.Species == 'Ecoli':
            self.RotateImage(math.degrees(self.TumbleAngle))
            self.CenterImage()

    def Displacement(self, Distance, Angle):
        dX =  Distance * math.cos(Angle)
        dY = -Distance * math.sin(Angle)
        return dX, dY

    def RandomMovement(self):
        Distance = random.random() * (self.Speed / 10)
        Angle = random.randint(0, 100) * self.Angle
        self.Move(Distance, Angle)

    def DrawTrajectory(self):
        TrajectoryPoints = self.Trajectory[:]
        TrajectoryPoints.append((self.X, self.Y))

        pygame.draw.aalines(Screen, RED, False, TrajectoryPoints)

class FMolecule:
    def __init__(self, InX, InY, Max):
        self.X_Ori = InX
        self.Y_Ori = InY
        self.Max = Max
        self.DiffusionFactor = 2
        self.SpaceFactor = 50

        # Gradient Drawing
        self.GradLvl = 20
        self.BaseLvl = 50
        self.Skew = 0   # Choose 0 or any number < or > 0
        self.GradStep = self.GetGradientStep()
        self.GradStepList = self.GetGradientStepList()
        self.GradDensityList = self.GetGradientDensityList()
        self.DensityLimit = 5 * Unit
        self.GradBaseColor = 'Blue'
        self.GradColorList = self.GetGradientColorList(baseColor=self.GradBaseColor)

        # Particle Drawing
        self.Particle_N = 200
        self.Particle_PerLayer = 2
        self.Particle_Radius = 2
        self.Particle_SpreadFactor = 1.11
        self.Particle_XY_Static = []
        self.InitializeStaticParticles()

    def InitializeStaticParticles(self):
        for i in range(int(self.Particle_N / self.Particle_PerLayer)):
            for j in range(self.Particle_PerLayer):
                X = self.X_Ori + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))
                Y = self.Y_Ori + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))
                self.Particle_XY_Static.append((X, Y))

    def Reposition(self):
        self.X_Ori = random.randint(W_S * 2 / 5, W_S * 3 / 5)
        self.Y_Ori = random.randint(H_S * 2 / 5, H_S * 3 / 5)
        self.Particle_XY_Static = []
        self.InitializeStaticParticles()

    def Move(self, dX, dY):
        New_Particle_XY_Static = []
        for (X, Y) in self.Particle_XY_Static:
            X += dX
            Y += dY
            New_Particle_XY_Static.append((X, Y))
        self.Particle_XY_Static = New_Particle_XY_Static

    def GetGradientStep(self):
        return math.floor(255 / self.GradLvl)

    def GetGradientStepList(self):
        StepList = list()
        if self.Skew:
            for i in range(self.GradLvl):
                SkewedStep = None
                if self.Skew > 0:   # positive Skew
                    SkewedStep = self.Max / (i + 1)
                else:   # negative Skew
                    SkewedStep = 255 - (self.Max / (i + 1))
                StepList.append(math.floor(SkewedStep))
        else:
            StepList = list(range(0, 255, self.GradStep))
        return StepList

    def GetGradientDensityList(self):
        GradDensityList = list()
        for gradStep in self.GradStepList:
            GradDensityList.append(self.Max * gradStep / self.GradStepList[-1])
        return GradDensityList

    def GetGradientColorList(self, baseColor=None):
        ColorList = list()
        for gradient in self.GradStepList:
            gradient_Scaled = math.floor((gradient / self.Max) * 255)
            ColorList.append(GetColorGradient(255 - gradient, baseColor))
            # ColorList.append(GetColorGradient(255 - gradient_Scaled, baseColor))
        return ColorList

    def Draw(self, pattern='particle', dynamics='static'):

        if pattern == 'gradient':
            for i in range(self.GradLvl):
                # Skip the first circle
                if i == 0:
                    continue
                else:
                    Density = self.GradDensityList[i]
                    Radius = self.GetRadius(Density)

                    if Density < self.DensityLimit:
                        continue

                    Color = self.GradColorList[i]
                    pygame.draw.circle(Screen, Color, (self.X_Ori, self.Y_Ori), Radius)

        # # Old gradient
        # elif pattern == 'gradient':
        #     for i in range(self.GradLvl):
        #         Color = (self.BaseLvl, self.BaseLvl, self.BaseLvl + (255 - self.BaseLvl) * ((i + 1) /self.GradLvl))
        #         pygame.draw.circle(Screen, Color, (self.X_Ori, self.Y_Ori), 100 / (i + 1))
        #     pygame.draw.circle(Screen, BLUE, (self.X_Ori, self.Y_Ori), 5)

        elif ((pattern == 'particle') & (dynamics == 'dynamic')):
            for i in range(self.Particle_N):
                X = self.X_Ori + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))
                Y = self.Y_Ori + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))
                pygame.draw.circle(Screen, BLUE, (X, Y), self.Particle_Radius)

        elif ((pattern == 'particle') & (dynamics == 'static')):
            for XY in self.Particle_XY_Static:
                pygame.draw.circle(Screen, BLUE, XY, self.Particle_Radius)

        else:
            assert True, 'Unsupported molecule distribution pattern for drawing: %s' % pattern

    def GetAmount(self, X, Y):
        return self.Diffusion(X, Y)

    # TODO:Diffusion will be updated to a dynamic version
    def Diffusion(self, X, Y):
        # DK - debugging purposes
        Dist = math.sqrt(((X - self.X_Ori) / W_S) ** 2 + ((Y - self.Y_Ori) / H_S) ** 2)
        return self.Max / max(1, Dist * 30)
        # return self.Max / ((X - self.X_Ori) ** 2 + (Y - self.Y_Ori) ** 2)

    def GetRadius(self, Amount):   # Need to be updated according to the diffusion equation
        return ((self.Max / Amount) ** (1. / self.DiffusionFactor)) * self.SpaceFactor
        # return (self.Max / Amount) ** (1. / self.DiffusionFactor)

class FControl:
    def __init__(self):
        # self.FPS = 30
        self.MovementResolution = 30
        self.MessageWelcome = 'Welcome to Bacterial Chemotaxis!'
        self.Pos_Welcome = GetMidPoint(Center, CenterTop)
        self.Message = ''
        self.MessageTimer = 3000

        self.InstructionSwitch = False
        self.Instructions = {
            'LEFT'  : 'Move Glucose LEFTWARD',
            'RIGHT' : 'Move Glucose RIGHTWARD',
            'UP'    : 'Move Glucose UPWARD',
            'DOWN'  : 'Move Glucose DOWNWARD',
            '-'     : 'Glucose Level --',
            '+'     : 'Glucose Level ++',
            '['     : 'Glucose Movement Resolution --',
            ']'     : 'Glucose Movement Resolution ++',
            'G'     : 'Reposition Glucose Position',
            'O'     : 'Reinitialize Ecoli Position',
            'R'     : 'Reverse Ecoli Tumbling Direction',
            ';'     : 'Ecoli Tumbling Angle --',
            '"'     : 'Ecoli Tumbling Angle ++',
            '<'     : 'Ecoli Speed --',
            '>'     : 'Ecoli Speed ++',
            'M'     : 'Ecoli Mechanistic Mode Switch',
            'T'     : 'Display Trajectory Switch',
            'I'     : 'Display Instruction Switch',
            'S'     : 'Display Score Switch',
            'A'     : 'Display Status Switch',

            # 'D'     : 'Transparency Display Switch',
        }
        self.InstructionText = ''
        self.SetInstructionText()

        self.MoleculeGradientText = ''
        self.MoleculeGradientColor = list()

        self.Score = 0
        self.ScoreSwitch = False

        self.Time = 0
        self.ScoreSwitch = False

        # DK - debugging purposes
        self.StatusSwitch = True
        # self.StatusSwitch = False

        self.TransparencySwitch = False


    #
    #     self.ControlBoard = None
    #     self.InitializeControlBoard()
    #
    # def InitializeControlBoard(self):
    #     self.ControlBoard = pygame.Surface(Screen_Size, pygame.SRCALPHA)
    #     self.ControlBoard.fill((0, 0, 0, 255))
    #     # pygame.draw.rect(self.ControlBoard, (0, 0, 0, 255), (self.Radius, self.Radius), self.Radius)

    def SetInstructionText(self):
        self.InstructionText = 'Instructions\n'
        for Key, Value in self.Instructions.items():
            Space = " " * (6 - len(Key))
            self.InstructionText = self.InstructionText + "  " + Key + Space + ': ' + Value + '\n'

    def DisplayInstructions(self):
        TextLines = self.InstructionText.splitlines()
        Height = Font_Monospace.get_linesize() + 2
        X, Y = Screen.get_rect().topleft
        Color = None
        for i, TextLine in enumerate(TextLines):
            if 'Glucose' in TextLine:
                Color = BLUE
            elif 'Ecoli' in TextLine:
                Color = RED
            else:
                Color = BLACK
            Text = Font_Monospace.render(TextLine, True, Color)
            Text_Rect = Text.get_rect()
            Text_Rect.topleft = (X, Y + Height * i)
            Screen.blit(Text, Text_Rect)

    def SetMoleculeGradient(self, Molecule, MolName='Molecule'):
        self.MoleculeGradientText = '%s\n' % MolName
        self.MoleculeGradientColor.append(WHITE)

        for Density, Color in zip(Molecule.GradDensityList, Molecule.GradColorList):
            if Density < Molecule.DensityLimit:
                continue
            self.MoleculeGradientText = self.MoleculeGradientText + "{:.2f}".format(Density/Unit) + UnitTxt + '\n'
            self.MoleculeGradientColor.append(Color)

        self.MoleculeGradientColor.reverse()

    def DisplayMoleculeGradient(self):
        TextLines = self.MoleculeGradientText.splitlines()
        TextLines.reverse()
        Height = Font_Sans.get_linesize() + 2
        X, Y = Screen.get_rect().bottomright
        for i, (TextLine, Color) in enumerate(zip(TextLines, self.MoleculeGradientColor)):
            Text = Font_Sans.render(TextLine, True, BLACK, Color)
            Text_Rect = Text.get_rect()
            Text_Rect.bottomright = (X, Y - Height * i)
            Screen.blit(Text, Text_Rect)

    def DisplayWelcome(self):
        Text = Font_Sans.render(self.MessageWelcome, True, BLACK, WHITE)
        Text_Rect = Text.get_rect()
        Text_Rect.center = self.Pos_Welcome
        Screen.blit(Text, Text_Rect)

    def DisplayInput(self):
        Text = Font_Sans.render(self.Message, True, BLACK)
        Text_Rect = Text.get_rect()
        Text_Rect.bottomleft = Screen.get_rect().bottomleft
        Screen.blit(Text, Text_Rect)
        self.MessageTimer -= 1

    def SetMessage(self, Key):
        assert Key in self.Instructions
        return "Input '" + Key + "' : " + self.Instructions[Key] + "   "

    def DisplayScore(self):
        Text = Font_Sans.render('Score: ' + str(round(self.Score)), True, RED)
        Text_Rect = Text.get_rect()
        Text_Rect.midtop = tuple(map(lambda i, j: i + j, Screen.get_rect().midtop, (0, Text.get_height())))
        Screen.blit(Text, Text_Rect)

    def DisplayTime(self):
        Text = Font_Sans.render('Time: ' + str(round(self.Time / 1000)), True, BLACK)
        Text_Rect = Text.get_rect()
        Text_Rect.midtop = Screen.get_rect().midtop
        Screen.blit(Text, Text_Rect)

    def DisplayStatus(self, Glucose_Total, Glucose_Ecoli, Glucose_Prev_Ecoli, Am, switch=False):
        dGlucose = (Glucose_Ecoli - Glucose_Prev_Ecoli) / Glucose_Total * 100
        StatusText = "   Total Glucose : " + "{:.2f} ".format(Glucose_Total / Unit) + UnitTxt + "\n" \
                     + " Glucose @ Ecoli :" + "{:.2f} ".format(Glucose_Ecoli/ Unit) + UnitTxt + "\n" \
                     + " dGlucose @ Ecoli : " + ["", "+"][dGlucose > 0] + "{:.5f}".format(dGlucose) + " %"
        # StatusText = "   Total Glucose : " + str(Glucose_Total) + "\n  Glucose @ Ecoli :" + "{:.5f}".format(Glucose_Ecoli)
        if switch:
            StatusText += '\n[MECHANISTIC MODE]'
            StatusText = StatusText + "\nAm level of Ecoli : " + "{:.5f} ".format(Am / Unit) + UnitTxt
            # StatusText = StatusText + "\nAm level of Ecoli : " + "{:.5f}".format(Am)
        else:
            StatusText += '\n[ALGORITHMIC MODE]'
        TextLines = StatusText.splitlines()
        Height = Font_Monospace.get_linesize() + 2
        X, Y = Screen.get_rect().topright
        Color = BLACK
        for i, TextLine in enumerate(TextLines):
            if 'MODE' in TextLine:
                Color = RED
            Text = Font_Monospace.render(TextLine, True, Color)
            Text_Rect = Text.get_rect()
            Text_Rect.topright = (X, Y + Height * i)
            Screen.blit(Text, Text_Rect)

def main():
    global TransparencySwitch
    Control = FControl()

    SimUnitTime = 0.1

    PetriDish = FEnvironment()
    Glucose = FMolecule(W_S * 3 / 5 , H_S * 2 / 5, 100 * Unit)
    Ecoli = FOrganism('A', W_S / 3, H_S / 2)
    Glucose_Now = Glucose.GetAmount(Ecoli.X, Ecoli.Y)
    Ecoli.Homeostasis(Glucose_Now)

    # Control.SetMoleculeGradient(Glucose, 'Glucose')

    # if Control.TransparencySwitch:
    #     PetriDish.DrawTransparentArea()

    ElapsedTime = 0
    PrevTime = datetime.now()

    SimState = True
    while SimState:
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

                # Glucose Control
                elif event.key == pygame.K_LEFT:
                    Glucose.X_Ori -= Control.MovementResolution
                    Glucose.Move(-Control.MovementResolution, 0)
                    Control.Message = Control.SetMessage('LEFT')
                elif event.key == pygame.K_RIGHT:
                    Glucose.X_Ori += Control.MovementResolution
                    Glucose.Move(Control.MovementResolution, 0)
                    Control.Message = Control.SetMessage('RIGHT')
                elif event.key == pygame.K_UP:
                    Glucose.Y_Ori -= Control.MovementResolution
                    Glucose.Move(0, -Control.MovementResolution)
                    Control.Message = Control.SetMessage('UP')
                elif event.key == pygame.K_DOWN:
                    Glucose.Y_Ori += Control.MovementResolution
                    Glucose.Move(0, +Control.MovementResolution)
                    Control.Message = Control.SetMessage('DOWN')
                elif event.key == pygame.K_KP_PLUS or event.key == pygame.K_EQUALS:
                    Glucose.Max *= 2
                    Control.Message = Control.SetMessage('+')
                elif event.key == pygame.K_KP_MINUS or event.key == pygame.K_MINUS:
                    Glucose.Max /= 2
                    Control.Message = Control.SetMessage('-')
                elif event.key == pygame.K_LEFTBRACKET:
                    Control.MovementResolution -= 1
                    Control.Message = Control.SetMessage('[')
                elif event.key == pygame.K_RIGHTBRACKET:
                    Control.MovementResolution += 1
                    Control.Message = Control.SetMessage(']')
                elif event.key == pygame.K_g:
                    Glucose.Reposition()
                    Control.Message = Control.SetMessage('G')

                # Ecoli Control
                elif event.key == pygame.K_o:
                    Ecoli.Reinitialize()
                    Control.Message = Control.SetMessage('O')
                elif event.key == pygame.K_r:
                    Ecoli.TumbleAngle = -Ecoli.TumbleAngle
                    Control.Message = Control.SetMessage('R')
                elif event.key == pygame.K_SEMICOLON:
                    Ecoli.TumbleAngle -= pi/50
                    Control.Message = Control.SetMessage(';')
                elif event.key == pygame.K_QUOTE:
                    Ecoli.TumbleAngle += pi/50
                    Control.Message = Control.SetMessage('"')
                elif event.key == pygame.K_COMMA:
                    Ecoli.Speed -= Ecoli.Speed_Max / 10
                    Control.Message = Control.SetMessage('<')
                elif event.key == pygame.K_PERIOD:
                    Ecoli.Speed += Ecoli.Speed_Max / 10
                    Control.Message = Control.SetMessage('>')
                elif event.key == pygame.K_t:
                    Ecoli.TrajectorySwitch = not Ecoli.TrajectorySwitch
                    Control.Message = Control.SetMessage('T')
                elif event.key == pygame.K_m:
                    Ecoli.MechanisticModeSwitch = not Ecoli.MechanisticModeSwitch
                    Control.Message = Control.SetMessage('M')

                # Control panel
                elif event.key == pygame.K_i:
                    Control.InstructionSwitch = not Control.InstructionSwitch
                    Control.Message = Control.SetMessage('I')
                elif event.key == pygame.K_s:
                    Control.ScoreSwitch = not Control.ScoreSwitch
                    Control.Message = Control.SetMessage('S')
                elif event.key == pygame.K_a:
                    Control.StatusSwitch = not Control.StatusSwitch
                    Control.Message = Control.SetMessage('A')

                # # Irreversible transparency option disabled
                # elif event.key == pygame.K_d:
                #     Control.TransparencySwitch = not Control.TransparencySwitch
                #     Control.Message = Control.SetMessage('D')


        # if Control.TransparencySwitch:
        #     Screen.set_clip(None)

        Screen.fill(GRAY1)

        # if Control.TransparencySwitch:
        #     Topleft = (PetriDish.X - PetriDish.Radius, PetriDish.Y - PetriDish.Radius)
        #     ClipRect = pygame.Rect(Topleft, (PetriDish.Radius * 2, PetriDish.Radius * 2))
        #     # ClipRect = pygame.Rect((0, 0), Screen_Size)
        #     Screen.set_clip(ClipRect)

        PetriDish.Draw()
        Glucose.Draw()

        while ElapsedTime >= SimUnitTime:
            Glucose_Now = Glucose.GetAmount(Ecoli.X, Ecoli.Y)
            Ecoli.Chemotaxis(Glucose_Now)
            ElapsedTime -= SimUnitTime

        # if PetriDish.CheckOutOfBound(Ecoli.X, Ecoli.Y):
        #     Ecoli.X = Ecoli.X_Prev
        #     Ecoli.Y = Ecoli.Y_Prev
        #     Ecoli.Tumble()

        if Ecoli.TrajectorySwitch:
            Ecoli.DrawTrajectory()

        Ecoli.Draw()

        # if Control.TransparencySwitch:
        #     Screen.blit(PetriDish.TransparentCircleArea, Topleft)

        if Control.Time < 500:
            Control.DisplayWelcome()
        if Control.MessageTimer > 0:
            Control.DisplayInput()
        if Control.InstructionSwitch:
            Control.DisplayInstructions()
        if Control.ScoreSwitch:
            Control.Score += Glucose_Now / 10
            Control.DisplayScore()
            if Glucose_Now > (Glucose.Max * 0.999):
                Glucose.Reposition()

        Glucose_Now = Glucose.GetAmount(Ecoli.X, Ecoli.Y)
        if Control.StatusSwitch:
            Control.DisplayStatus(Glucose.Max, Glucose_Now, Ecoli.Glucose_Prev, Ecoli.Am, switch=Ecoli.MechanisticModeSwitch)

        Control.Time += 1
        Control.DisplayTime()
        # Control.DisplayMoleculeGradient()

        pygame.display.update()

    pygame.quit()
    sys.exit()

if __name__ == '__main__':
    main()

"""
def BlitRotate1(Surface, Image, X, Y, X_Ori, Y_Ori, AngleInDegree):
    # offset from pivot to center
    Image_Rect = Image.get_rect(topleft = (X - X_Ori, Y - Y_Ori))
    OffCenterToPivot = pygame.math.Vector2([X, Y]) - Image_Rect.center

    # rotated offset from pivot to center
    RotatedOffSet = OffCenterToPivot.rotate(-AngleInDegree)

    # rotated image center
    RotatedImage_Center = (X - RotatedOffSet.y, Y - RotatedOffSet.y)

    # get a rotated image
    RotatedImage = pygame.transform.rotate(Image,AngleInDegree)
    RotatedImage_Rect = RotatedImage.get_rect()

    # rotate and blit the image
    Surface.blit(RotatedImage, RotatedImage_Rect)
    pygame.draw.rect(Surface, (255, 255, 0), (*RotatedImage_Rect.topleft, *RotatedImage.get_size()), 2)

def BlitRotate2(Surface, Image, Topleft, AngleInDegree):
    # get a rotated image
    RotatedImage = pygame.transform.rotate(Image,AngleInDegree)
    RotatedImage_Rect = RotatedImage.get_rect(center = Image.get_rect(topleft=Topleft).center)

    # rotate and blit the image
    Surface.blit(RotatedImage, RotatedImage_Rect.topleft)
    pygame.draw.rect(Surface, (255, 255, 0), RotatedImage_Rect, 2)

"""
