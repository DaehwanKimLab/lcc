import sys
import pygame
import random
from datetime import datetime
import SimModule
import numpy as np

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

GlucoseName = "L"
EcoliName = "E"
HomeostasisMolName = "Am"

# Utilities
def GetColorGradient(Fade, baseColor=None):
    assert Fade <= 255, 'ERROR: Fade factor is higher than 255!'
    if baseColor == 'Blue':
        return (Fade, Fade, 255)
    else:  # Assumes White
        return (Fade, Fade, Fade)

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
Font_Monospace = pygame.font.SysFont('monospace', 18, True)

# Load model
State = SimModule.FState()
Data = SimModule.FDataset()
DataManager = SimModule.FDataManager()
Sim = SimModule.FSimulation(State, Data, DataManager)

# Initialize model
Sim.Initialize()

class FEnvironment:
    def __init__(self, InX=W_S/2, InY=H_S/2, InShape='circle', InRadius=W_S*0.3, InThickness=5):
        self.X = InX
        self.Y = InY
        self.Shape = InShape
        self.Radius = int(InRadius)
        self.Thickness = InThickness
        self.TransparentCircleArea = None

    def Draw(self, shape=None):
        if shape == 'circle':
            pygame.draw.circle(Screen, YELLOW_FAINT, (self.X, self.Y), self.Radius)
            pygame.draw.circle(Screen, GRAY4, (self.X, self.Y), self.Radius, self.Thickness)
        else:
            pygame.draw.rect(Screen, YELLOW_FAINT, ((0, 0), (W_S, H_S)))

    def DrawTransparentArea(self):
        self.TransparentCircleArea = pygame.Surface((self.Radius * 2, self.Radius * 2), pygame.SRCALPHA)
        self.TransparentCircleArea.fill((255, 255, 255, 255))
        pygame.draw.circle(self.TransparentCircleArea, (0, 0, 0, 0), (self.Radius, self.Radius), self.Radius)

    def CheckOutOfBound(self, X, Y):
        pass

class FOrganism:
    def __init__(self, InName, InSpecies):
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
        self.Glucose_Prev = 0
        self.Am = 0
        self.SimCount = 0

        # Trajectory
        self.TrajectorySwitch = True
        self.Trajectory = dict()
        self.TrajectoryColor = list()

    def Draw(self):
        if self.Species == 'Ecoli':
            Color = YELLOW
            # pygame.draw.circle(Screen, Color, (self.X, self.Y), 5)
            dX =  np.cos(self.Angle) * -self.BodyLength
            dY = -np.sin(self.Angle) * -self.BodyLength
            X_BodyEnd = self.X + dX
            Y_BodyEnd = self.Y + dY
            X_TailEnd = self.X + self.FlagellaLength_Fold * dX
            Y_TailEnd = self.Y + self.FlagellaLength_Fold * dY
            for i in range(self.X.size):
                if i == self.X.size - 1:
                    Color = RED
                pygame.draw.line(Screen, Color, (self.X[0, i], self.Y[0, i]), (X_BodyEnd[0, i], Y_BodyEnd[0, i]), self.BodyThickness)
                pygame.draw.line(Screen, Color, (self.X[0, i], self.Y[0, i]), (X_TailEnd[0, i], Y_TailEnd[0, i]), self.FlagellaThickness)


    def SetPosition(self, Position):
        self.X = Position[0]
        self.Y = Position[1]
        self.Angle = Position[2]
        # Threshold value (Position[3]) is not used

    def Receptivity(self, N_SimulationsToPass=100):
        # Perform 100 simulations
        # 99 simulations without spatial simulation
        for _ in range(N_SimulationsToPass):
            Sim.SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution()
        # 1 with spatial simulation
        Sim.SimLoop_WithSpatialSimulation()

    def ReportStatus(self):
        # for debugging
        self.Am = Sim.GetCountByName(HomeostasisMolName)
        self.SimCount += 1
        GlucoseLvl = Sim.GetCountFromDistributionByNameAndPos(GlucoseName, EcoliName)
        SimStep = Sim.GetSimStep()
        Delta = (GlucoseLvl - self.Glucose_Prev) / GlucoseLvl * 100
        # print("SimStep {:06d} [Chemotaxis  {:06d}] Glucose:{:.6f} {} ({}{:.4f}%) Am:{:.6f} {} (X:{:.2f}, Y:{:.2f}, {:3.1f} degree)".format
        #       (SimStep, self.SimCount, GlucoseLvl / Unit / NA, UnitTxt, ("+" if Delta >= 0 else ""), Delta, self.Am / Unit / NA, UnitTxt , self.X, self.Y, self.Angle / pi * 180))

    # def HomeostasisMessage(self):
    #     print("[Homeostasis {:06d}] Glucose:{:.6f}{} Am:{:.6f}{}".format(self.SimCount, GlucoseLvl / Unit / NA, UnitTxt, self.Am / Unit / NA, UnitTxt))

    def Homeostasis(self, MolNme):
        Sim.Homeostasis()   # Input 'Am' here
        self.Glucose_Prev = Sim.GetCountFromDistributionByNameAndPos(GlucoseName, EcoliName)
        self.SetPosition(Sim.GetPositionXYAngleByName(EcoliName))
        self.InitializeTrajectory()
        self.ReportStatus()

    def InitializeTrajectory(self):
        for i in range(self.X.size):
            self.Trajectory[i] = [(self.X[0, i], self.Y[0, i])]
            self.TrajectoryColor.append(tuple(np.random.randint(0, 255, 3)))

    def AddToTrajectory(self):
        for i in range(self.X.size):
            self.Trajectory[i].append((self.X[0, i], self.Y[0, i]))

    def DrawTrajectory(self):
        for i in range(len(self.Trajectory)):
            pygame.draw.aalines(Screen, self.TrajectoryColor[i], False, self.Trajectory[i])

class FMolecule:
    def __init__(self, InName, InX, InY):
        self.Name = InName
        self.X = InX
        self.Y = InY
        self.DiffusionFactor = 2
        self.SpaceFactor = 50

        # Heatmap Drawing
        self.ReductionFactor = 5

        # Particle Drawing
        self.Particle_N = 300
        self.Particle_PerLayer = 3
        self.Particle_Radius = 2
        self.Particle_SpreadFactor = 1.2
        self.Particle_XY_Static = []
        self.InitializeStaticParticles()

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

    def Draw(self, Data, pattern='heatmap'):
    # def Draw(self, pattern='particle', dynamics='static'):

        if pattern == 'heatmap':
            Max = np.max(Data)
            for x in range(0, Data.shape[0], self.ReductionFactor):
                for y in range(0, Data.shape[1], self.ReductionFactor):
                    intensity = (Data[x][y] / Max) * 255
                    # print(x, y, intensity)
                    color = (255 - intensity, 255 - intensity, 255)
                    pygame.draw.rect(Screen, color, ((x, y), (self.ReductionFactor, self.ReductionFactor)))

        elif pattern == 'particle':
            for XY in self.Particle_XY_Static:
                pygame.draw.circle(Screen, BLUE, XY, self.Particle_Radius)

        else:
            assert True, 'Unsupported molecule distribution pattern for drawing: %s' % pattern

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
            'P'     : 'Pause Visualization',

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

        # Pause
        self.MessagePause = 'PAUSE'
        self.PauseSwitch = False

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
        Text = Font_Sans.render('Simulation Time: ' + str(round(self.Time)), True, BLACK)
        Text_Rect = Text.get_rect()
        Text_Rect.midtop = Screen.get_rect().midtop
        Screen.blit(Text, Text_Rect)

    # TODO: Update display status
    # def DisplayStatus(self, Glucose_Total, Glucose_Ecoli, Glucose_Prev_Ecoli, Am):
    def DisplayStatus(self, Glucose_Ecoli, Glucose_Prev_Ecoli, Am_Ecoli):
        Glucose_Now = Glucose_Ecoli[0, -1]
        Glucose_Prev = Glucose_Prev_Ecoli[0, -1]
        Am = Am_Ecoli[-1, 0, 0]
        dGlucose = (Glucose_Now - Glucose_Prev) / Glucose_Now * 100

        # StatusText = "   Total Glucose : " + "{:.2f} ".format(Glucose_Total / Unit/ NA) + UnitTxt + "\n" \
        StatusText = " Glucose @ Ecoli :" + "{:.2f} ".format(Glucose_Now/ Unit / NA) + UnitTxt + "\n" \
                     + " dGlucose @ Ecoli : " + ("+" if dGlucose >= 0 else "") + "{:.5f}".format(dGlucose) + " %" \
                     + "\nAm level of Ecoli : " + "{:.5f} ".format(Am / Unit / NA) + UnitTxt   # Get the last E coli's info

        TextLines = StatusText.splitlines()
        Height = Font_Monospace.get_linesize() + 2
        X, Y = Screen.get_rect().topright
        Color = BLACK
        for i, TextLine in enumerate(TextLines):
            if 'Ecoli' in TextLine:
                Color = RED
            Text = Font_Monospace.render(TextLine, True, Color)
            Text_Rect = Text.get_rect()
            Text_Rect.topright = (X, Y + Height * i)
            Screen.blit(Text, Text_Rect)

    def DisplayPause(self):
        Text = Font_Sans.render(self.MessagePause, True, BLACK, WHITE)
        Text_Rect = Text.get_rect()
        Text_Rect.center = self.Pos_Welcome
        Screen.blit(Text, Text_Rect)


def main():
    Control = FControl()

    SimUnitTime = 0.1

    PetriDish = FEnvironment()

    # TODO: Communicate to initialize in Sim
    Glucose = FMolecule(GlucoseName, 800, 500)
    Ecoli = FOrganism(EcoliName, 'Ecoli')

    Ecoli.Homeostasis(HomeostasisMolName)

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

                # Glucose Control
                elif event.key == pygame.K_LEFT:
                    Glucose.X -= Control.MovementResolution
                    Glucose.Move(-Control.MovementResolution, 0)
                    Control.Message = Control.SetMessage('LEFT')
                elif event.key == pygame.K_RIGHT:
                    Glucose.X += Control.MovementResolution
                    Glucose.Move(Control.MovementResolution, 0)
                    Control.Message = Control.SetMessage('RIGHT')
                elif event.key == pygame.K_UP:
                    Glucose.Y -= Control.MovementResolution
                    Glucose.Move(0, -Control.MovementResolution)
                    Control.Message = Control.SetMessage('UP')
                elif event.key == pygame.K_DOWN:
                    Glucose.Y += Control.MovementResolution
                    Glucose.Move(0, +Control.MovementResolution)
                    Control.Message = Control.SetMessage('DOWN')
                elif event.key == pygame.K_KP_PLUS or event.key == pygame.K_EQUALS:
                    # Glucose.Max *= 2
                    Control.Message = Control.SetMessage('+')
                elif event.key == pygame.K_KP_MINUS or event.key == pygame.K_MINUS:
                    # Glucose.Max /= 2
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
                    # Ecoli.Reinitialize()
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
                    # Ecoli.Speed -= Ecoli.Speed_Max / 10
                    Control.Message = Control.SetMessage('<')
                elif event.key == pygame.K_PERIOD:
                    # Ecoli.Speed += Ecoli.Speed_Max / 10
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
                elif event.key == pygame.K_p:
                    Control.PauseSwitch = not Control.PauseSwitch
                    Control.Message = Control.SetMessage('P')

        Screen.fill(GRAY1)

        PetriDish.Draw()
        Glucose.Draw(Sim.GetDistributionByName(GlucoseName), pattern='heatmap')
        # Glucose.Draw(Sim.GetDistributionByName(GlucoseName), pattern='particle')


        if Control.PauseSwitch:
            Control.DisplayPause()

        else:
            while ElapsedTime >= SimUnitTime:

                # Previously Ecoli.Chemotaxis, now without decision making
                Ecoli.Receptivity()
                Ecoli.SetPosition(Sim.GetPositionXYAngleByName(EcoliName))
                Ecoli.ReportStatus()

                ElapsedTime -= SimUnitTime
                Control.Time += 1

        if Ecoli.TrajectorySwitch:
            Ecoli.AddToTrajectory()
            Ecoli.DrawTrajectory()

        Ecoli.Draw()

        Glucose_Now = Sim.GetCountFromDistributionByNameAndPos(GlucoseName, EcoliName)

        if Control.Time < 50:
            Control.DisplayWelcome()
        if Control.MessageTimer > 0:
            Control.DisplayInput()
        if Control.InstructionSwitch:
            Control.DisplayInstructions()
        if Control.ScoreSwitch:
            Control.Score += Glucose_Now / 10
            Control.DisplayScore()
            # if Glucose_Now > (Glucose.Max * 0.999):
            #     Glucose.Reposition()

        if Control.StatusSwitch:
            Control.DisplayStatus(Glucose_Now, Ecoli.Glucose_Prev, Ecoli.Am)

        # Update Glucose Prev
        Ecoli.Glucose_Prev = Glucose_Now

        Control.DisplayTime()
        # Control.DisplayMoleculeGradient()


        pygame.display.update()

    pygame.quit()
    sys.exit()

if __name__ == '__main__':
    main()
