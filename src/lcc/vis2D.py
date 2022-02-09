import sys
import pygame
import math

def main():
    pygame.init()

    WHITE = (255, 255, 255)
    GRAY = (50, 50, 50)
    YELLOW = (255, 255, 0)
    BLUE = (0, 0, 255)

    # Screen
    Screen_Size = W_S, H_S = 800, 600
    Center = W_S /2, H_S / 2
    Screen = pygame.display.set_mode(Screen_Size)

    Title = "Vis2D"
    pygame.display.set_caption(Title)

    # Petridish
    Petridish_Radius = 250
    def Draw_Petridish():
        pygame.draw.circle(Screen, GRAY, Center, Petridish_Radius)
        pygame.draw.circle(Screen, WHITE, Center, Petridish_Radius, 5)
        return

    # Glucose
    X_G = W_S / 2
    Y_G = H_S / 2
    Glucose_Max = 1000
    Glucose_Prev = 0
    Glucose_Now = 0
    Glucose_GradLvl = 5
    BaseLvl = 50

    def Draw_Glucose(Pos):
        for i in range(Glucose_GradLvl):
            Color = (BaseLvl, BaseLvl, BaseLvl + int((255 - BaseLvl) / (Glucose_GradLvl - i)))
            pygame.draw.circle(Screen, Color, Pos, 100 / (i + 1))

    def GetGlucose(X, Y):
        return Glucose_Max / ((X - X_G) ** 2 + (Y - Y_G) ** 2)

    # E coli
    Angle_E = 0
    Length_E = 5
    X_E = W_S / 3
    Y_E = H_S / 3
    X_E_prev = X_E
    Y_E_prev = Y_E
    def Draw_Ecoli():
        pygame.draw.circle(Screen, YELLOW, (X_E, Y_E), 5)
        # pygame.draw.rect(Screen, WHITE, (X_E, Y_E, X_E + 1, Y_E + 1), 1)

    def Move_Ecoli(Now, Prev):
        if Now > Prev:
            return 0.05, 0
        else:
            return 0.005, 30

    SimState = True
    while SimState:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                SimState = False

        Draw_Petridish()
        Draw_Glucose(Center)

        Glucose_Now = GetGlucose(X_E, Y_E)
        dDistance, dAngle = Move_Ecoli(Glucose_Now, Glucose_Prev)
        X_E += dDistance
        Y_E += dDistance
        Angle_E += dAngle
        # X_E, Y_E = CheckOutOfBound(X_E, Y_E)

        Draw_Ecoli()

        # Update Previous Glucose level
        Glucose_Prev = Glucose_Now

        pygame.display.update()

    pygame.quit()
    sys.exit()

if __name__ == '__main__':
    main()




#
# def CheckOutOfBound(X, Y):
#     if X ** 2 + Y ** 2 > Petridish_Radius:
#         return X_E_prev, Y_E_prev
#     else:
#         return X, Y
#
# def UpdatePrevGlucose():
#     Glucose_Prev = Glucose_Now