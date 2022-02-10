import sys
import pygame
import math
import random

# Colors
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
GRAY1 = (150, 150, 150)
GRAY2 = (100, 100, 100)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
YELLOW = (255, 255, 0)
BLUE = (0, 0, 255)
MAGENTA = (255, 0, 255)
CYAN = (0, 255, 255)

# Global variables
pi = math.pi

# pygame
pygame.init()

Screen_Size = W_S, H_S = 1200, 800
Screen = pygame.display.set_mode(Screen_Size)
Center = W_S / 2, H_S /2

# # Transparent control board
# ControlBoard = pygame.Surface(Screen_Size, pygame.SRCALPHA)
# ControlBoard.fill((0, 0, 0, 255))


Title = "Vis2D"
pygame.display.set_caption(Title)

Font_Input = pygame.font.Font('freesansbold.ttf', 20)
Font_Instruction = pygame.font.SysFont('monospace', 15, True)

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
            pygame.draw.circle(Screen, GRAY1, (self.X, self.Y), self.Radius)
            pygame.draw.circle(Screen, GRAY2, (self.X, self.Y), self.Radius, self.Thickness)

    def DrawTransparentArea(self):
        self.TransparentCircleArea = pygame.Surface((self.Radius * 2, self.Radius * 2), pygame.SRCALPHA)
        self.TransparentCircleArea.fill((255, 255, 255, 255))
        pygame.draw.circle(self.TransparentCircleArea, (0, 0, 0, 0), (self.Radius, self.Radius), self.Radius)

    def CheckOutOfBound(self, X, Y):
        pass

class FOrganism:
    def __init__(self, InSpecies, InX, InY, InA=0, InSpeedMax=0.05):
        self.Species = InSpecies
        self.X_Ori = InX
        self.Y_Ori = InY
        self.A = InA   # Angle
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

    def LoadImage(self):
        if self.Species == 'Ecoli':
            self.Image = pygame.image.load('pygamelib/ecoli.png')
            self.Image = pygame.transform.scale(self.Image, (self.Image_Size_X, self.Image_Size_Y))
            self.Image_Rect = self.Image.get_rect()
            self.RotateImage(220)

    def RotateImage(self, AngleInDegree):
        self.Image = pygame.transform.rotate(self.Image, AngleInDegree)

    def CenterImage(self):
        # Rect = self.Image.get_rect()
        # self.Image.center = Rect.center
        pass

    def Draw(self):
        if self.Species == 'Ecoli':
            # BlitRotate2(Screen, self.Image, self.Image.topleft, math.degrees(self.TumbleAngle))
            # Image = pygame.surface.Surface((100, 100))
            Screen.blit(self.Image, (self.X - self.Image_Size_X / 2, self.Y - self.Image_Size_Y / 2))
        else:
            Color = YELLOW
            # pygame.draw.circle(Screen, Color, (self.X, self.Y), 5)
            dX = math.sin(self.A) * -13
            dY = math.cos(self.A) * -13
            pygame.draw.line(Screen, Color, (self.X, self.Y), (self.X + dX, self.Y + dY), 7)
            pygame.draw.line(Screen, Color, (self.X, self.Y), (self.X + 2 * dX, self.Y + 2 * dY), 3)

    def Chemotaxis(self, GlucoseLvl):
        if GlucoseLvl > self.Glucose_Prev:
            self.Move(self.Speed)
        else:
            self.Tumble()

        # Update
        self.Glucose_Prev = GlucoseLvl

    def Move(self, Distance):
        dX, dY = Displacement(Distance, self.A)
        self.X += dX
        self.Y += dY

    def Tumble(self):
        self.A += self.TumbleAngle
        self.Move(self.Speed / 10)
        self.Trajectory.append((self.X, self.Y))
        if self.Species == 'Ecoli':
            self.RotateImage(math.degrees(self.TumbleAngle))
            self.CenterImage()

    def DrawTrajectory(self):
        TrajectoryPoints = self.Trajectory[:]
        TrajectoryPoints.append((self.X, self.Y))

        pygame.draw.aalines(Screen, RED, False, TrajectoryPoints)

class FMolecule:
    def __init__(self, InX, InY, Max):
        self.X_Ori = InX
        self.Y_Ori = InY
        self.Max = Max
        self.Amount_now = 0
        self.Amount_prev = 0

        # Gradient Drawing
        self.GradLvl = 5
        self.BaseLvl = 50

        # Particle Drawing
        self.Particle_N = 70
        self.Particle_Radius = 2
        self.Particle_SpreadFactor = 1.2

    def Draw(self, pattern='gradient'):
        if pattern == 'gradient':
            for i in range(self.GradLvl):
                Color = (self.BaseLvl, self.BaseLvl, self.BaseLvl + (255 - self.BaseLvl) * ((i + 1) /self.GradLvl))
                pygame.draw.circle(Screen, Color, (self.X_Ori, self.Y_Ori), 100 / (i + 1))
            pygame.draw.circle(Screen, BLUE, (self.X_Ori, self.Y_Ori), 5)

        elif pattern == 'particle':
            for i in range(self.Particle_N):
                X = self.X_Ori + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))
                Y = self.Y_Ori + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))
                pygame.draw.circle(Screen, BLUE, (X, Y), self.Particle_Radius)

        else:
            assert True, 'Unsupported Species for drawing: %s' % pattern

    def GetAmount(self, X, Y):
        return self.Diffusion(X, Y)

    # Diffusion will be updated to a dynamic version
    def Diffusion(self, X, Y):
        return self.Max / ((X - self.X_Ori) ** 2 + (Y - self.Y_Ori) ** 2)

class FControl:
    def __init__(self):
        # self.FPS = 30
        self.MovementResolution = 30
        self.Message = 'Welcome to Bacterial Chemotaxis!'
        self.MessageTimer = 3000

        self.InstructionSwitch = False
        self.Instructions = {
            'LEFT'  : 'Move Glucose LEFTWARD',
            'RIGHT' : 'Move Glucose RIGHTWARD',
            'UP'    : 'Move Glucose UPWARD',
            'DOWN'  : 'Move Glucose DOWNWARD',
            '['     : 'Reduce Glucose Movement Resolution',
            ']'     : 'Increment Glucose Movement Resolution',
            'O'     : 'Reinitialize Ecoli Position',
            'R'     : 'Reverse Ecoli Tumbling Direction',
            ';'     : 'Reduce Ecoli Tumbling Angle',
            '"'     : 'Increment Ecoli Tumbling Angle',
            '<'     : 'Reduce Ecoli Speed',
            '>'     : 'Increment Ecoli Speed',
            'T'     : 'Trajectory Mode Switch',
            'I'     : 'Instruction Display Switch',
            'C'     : 'Transparency Display Switch',
        }
        self.InstructionText = ''
        self.SetInstructionText()

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
        if self.MessageTimer > 0:
            TextLines = self.InstructionText.splitlines()
            Height = Font_Instruction.get_linesize() + 2
            X, Y = Screen.get_rect().topleft
            for i, TestLine in enumerate(TextLines):
                Text = Font_Instruction.render(TestLine, True, BLACK)
                Text_Rect = Text.get_rect()
                Text_Rect.topleft = (X, Y + Height * i)
                Screen.blit(Text, Text_Rect)

    def DisplayInput(self):
        if self.MessageTimer > 0:
            Text = Font_Input.render(self.Message, True, BLACK)
            Text_Rect = Text.get_rect()
            Text_Rect.topright = Screen.get_rect().topright
            Screen.blit(Text, Text_Rect)
            self.MessageTimer -= 1
        if self.MessageTimer == 0:
            self.Message = ''

    def SetMessage(self, Key):
        assert Key in self.Instructions
        return "Input '" + Key + "' : " + self.Instructions[Key] + "   "

def Displacement(Distance, Angle):
    dX = Distance * math.sin(Angle)
    dY = Distance * math.cos(Angle)
    return dX, dY

def main():
    global TransparencySwitch
    Control = FControl()

    PetriDish = FEnvironment()
    Glucose = FMolecule(W_S * 3 / 5 , H_S * 2 / 5, 1000)
    Ecoli = FOrganism('A', W_S / 3, H_S / 3)

    PetriDish.DrawTransparentArea()

    SimState = True
    while SimState:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                SimState = False
            elif event.type == pygame.KEYDOWN:
                Control.MessageTimer = 5000
                if event.key == pygame.K_x:
                    SimState = False
                elif event.key == pygame.K_LEFT:
                    Glucose.X_Ori -= Control.MovementResolution
                    Control.Message = Control.SetMessage('LEFT')
                elif event.key == pygame.K_RIGHT:
                    Glucose.X_Ori += Control.MovementResolution
                    Control.Message = Control.SetMessage('RIGHT')
                elif event.key == pygame.K_UP:
                    Glucose.Y_Ori -= Control.MovementResolution
                    Control.Message = Control.SetMessage('UP')
                elif event.key == pygame.K_DOWN:
                    Glucose.Y_Ori += Control.MovementResolution
                    Control.Message = Control.SetMessage('DOWN')
                elif event.key == pygame.K_LEFTBRACKET:
                    Control.MovementResolution -= 1
                    Control.Message = Control.SetMessage('[')
                elif event.key == pygame.K_RIGHTBRACKET:
                    Control.MovementResolution += 1
                    Control.Message = Control.SetMessage(']')
                elif event.key == pygame.K_o:
                    Ecoli.X = Ecoli.X_Ori
                    Ecoli.Y = Ecoli.Y_Ori
                    Ecoli.Trajectory = [(Ecoli.X_Ori, Ecoli.Y_Ori)]
                    Control.Message = Control.SetMessage('O')
                elif event.key == pygame.K_r:
                    Ecoli.TumbleAngle = -Ecoli.TumbleAngle
                    Control.Message = Control.SetMessage('R')
                elif event.key == pygame.K_SEMICOLON:
                    Ecoli.TumbleAngle -= pi/10
                    Control.Message = Control.SetMessage(';')
                elif event.key == pygame.K_QUOTE:
                    Ecoli.TumbleAngle += pi/10
                    Control.Message = Control.SetMessage('"')
                elif event.key == pygame.K_LESS:
                    Ecoli.Speed -= Ecoli.Speed_Max / 10
                    Control.Message = Control.SetMessage('<')
                elif event.key == pygame.K_GREATER:
                    Ecoli.Speed += Ecoli.Speed_Max / 10
                    Control.Message = Control.SetMessage('>')
                elif event.key == pygame.K_t:
                    Ecoli.TrajectorySwitch = not Ecoli.TrajectorySwitch
                    Control.Message = Control.SetMessage('T')
                elif event.key == pygame.K_i:
                    Control.InstructionSwitch = not Control.InstructionSwitch
                    Control.Message = Control.SetMessage('I')

                # Irreversible
                elif event.key == pygame.K_c:
                    Control.TransparencySwitch = not Control.TransparencySwitch
                    Control.Message = Control.SetMessage('C')


        if Control.TransparencySwitch:
            Screen.set_clip(None)

        Screen.fill(WHITE)

        if Control.TransparencySwitch:
            Topleft = (PetriDish.X - PetriDish.Radius, PetriDish.Y - PetriDish.Radius)
            ClipRect = pygame.Rect(Topleft, (PetriDish.Radius * 2, PetriDish.Radius * 2))
            # ClipRect = pygame.Rect((0, 0), Screen_Size)
            Screen.set_clip(ClipRect)

        PetriDish.Draw()

        # For Gradient
        # Glucose.Draw()

        # For Particle
        Glucose.Draw(pattern='particle')

        Glucose_Now = Glucose.GetAmount(Ecoli.X, Ecoli.Y)
        Ecoli.Chemotaxis(Glucose_Now)

        # if PetriDish.CheckOutOfBound(Ecoli.X, Ecoli.Y):
        #     Ecoli.X = Ecoli.X_Prev
        #     Ecoli.Y = Ecoli.Y_Prev
        #     Ecoli.Tumble()

        if Ecoli.TrajectorySwitch:
            Ecoli.DrawTrajectory()

        Ecoli.Draw()

        if Control.TransparencySwitch:
            Screen.blit(PetriDish.TransparentCircleArea, Topleft)

        if Control.MessageTimer > 0:
            Control.DisplayInput()
        if Control.InstructionSwitch:
            Control.DisplayInstructions()



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