#include "writer.h"
#include <algorithm>

using namespace std;

void FWriter::SimVis2D() {
    std::cout << "Generating SimVis2D..." << std::endl;

    // write SimVis2D.py
    std::ofstream ofs(Option.SimVis2DFile.c_str());
    std::string endl = "\n";

    ofs << "import sys" << endl;
    ofs << "import pygame" << endl;
    ofs << "import random" << endl;
    ofs << "from datetime import datetime" << endl;
    ofs << "import numpy as np" << endl;
    ofs << "import SimModule" << endl;
    ofs << "import SimFunctions as SimF" << endl;
    ofs << endl;
    ofs << "random.seed(1)" << endl;
    ofs << "np.random.seed(1)" << endl;
    ofs << endl;
    ofs << "# Colors" << endl;
    ofs << "BLACK = (0, 0, 0)" << endl;
    ofs << "WHITE = (255, 255, 255)" << endl;
    ofs << "GRAY1 = (230, 230, 230)" << endl;
    ofs << "GRAY2 = (210, 210, 210)" << endl;
    ofs << "GRAY3 = (150, 150, 150)" << endl;
    ofs << "GRAY4 = (100, 100, 100)" << endl;
    ofs << "RED = (255, 0, 0)" << endl;
    ofs << "RED_DARK = (139, 0, 0)" << endl;
    ofs << "GREEN = (0, 255, 0)" << endl;
    ofs << "GREEN_DARK = (0, 100, 0)" << endl;
    ofs << "YELLOW = (255, 255, 0)" << endl;
    ofs << "YELLOW_FAINT = (200, 200, 150)" << endl;
    ofs << "BLUE = (0, 0, 255)" << endl;
    ofs << "BLUE_DARK = (0, 0, 139)" << endl;
    ofs << "MAGENTA = (255, 0, 255)" << endl;
    ofs << "CYAN = (0, 255, 255)" << endl;
    ofs << endl;
    ofs << "# Global variables" << endl;
    ofs << "NA = 6.0221409e+23" << endl;
    ofs << "pi = np.pi" << endl;
    ofs << endl;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");

    // Define indices for spatially distributed molecules and containers
    ofs << "GlucoseName = 'L'" << endl;
    ofs << "AutoinducerName = 'qL'" << endl;

    // Pathway dependent key molecules to monitor
    ofs << "HomeostasisMolName = [";
    for (auto& location: MolLoc) {
        if (location->Name == "L") { ofs << "'" << "Am" << "', "; }
        else if (location->Name == "qL") { ofs << "'" << "qAm" << "', "; }
    }
    ofs << "]" << endl;
    ofs << endl;

    ofs << "# Utilities" << endl;
    ofs << "def GetColorGradient(Fade, baseColor=None):" << endl;
    ofs << in+ "assert Fade <= 255, 'ERROR: Fade factor is higher than 255!'" << endl;
    ofs << in+ "if baseColor == 'Blue':" << endl;
    ofs << in+ in+ "return (Fade, Fade, 255)" << endl;
    ofs << in+ "else:  # Assumes White" << endl;
    ofs << in+ in+ "return (Fade, Fade, Fade)" << endl;
    ofs << endl;
    ofs << "# pygame" << endl;
    ofs << "pygame.init()" << endl;
    ofs << endl;
    ofs << "# Load model" << endl;
    ofs << "State = SimModule.FState()" << endl;
    ofs << "Data = SimModule.FDataset()" << endl;
    ofs << "DataManager = SimModule.FDataManager()" << endl;
    ofs << "SimM = SimModule.FSimulation(State, Data, DataManager)" << endl;
    ofs << endl;
    ofs << "# Initialize model" << endl;
    ofs << "SimM.Initialize()" << endl;
    ofs << endl;
    ofs << "# Determine global unit" << endl;
    ofs << "Unit = SimM.Unit" << endl;
    ofs << "UnitTxt = SimM.UnitTxt" << endl;
    ofs << endl;
    ofs << "Screen_Size = W_S, H_S = " << "SimM.GetDistWidth()" << "," << "SimM.GetDistHeight()" << endl;
    ofs << "Screen = pygame.display.set_mode(Screen_Size)" << endl;
    ofs << endl;
    ofs << "LEFT = 0" << endl;
    ofs << "MID_X = W_S / 2" << endl;
    ofs << "RIGHT = W_S" << endl;
    ofs << endl;
    ofs << "TOP = 0" << endl;
    ofs << "MID_Y = H_S / 2" << endl;
    ofs << "BOTTOM = H_S" << endl;
    ofs << endl;
    ofs << "CenterTop = (MID_X, TOP)" << endl;
    ofs << "Center = (MID_X, MID_Y)" << endl;
    ofs << endl;
    ofs << "def AddPos(A, B):" << endl;
    ofs << in+ "X_A, Y_A = A" << endl;
    ofs << in+ "X_B, Y_B = B" << endl;
    ofs << in+ "return (X_A + X_B, Y_A + Y_B)" << endl;
    ofs << endl;
    ofs << "def GetMidPoint(A, B):" << endl;
    ofs << in+ "X_A, Y_A = A" << endl;
    ofs << in+ "X_B, Y_B = B" << endl;
    ofs << in+ "return ((X_A + X_B) / 2, (Y_A + Y_B) / 2)" << endl;
    ofs << endl;
    ofs << "# # Transparent control board" << endl;
    ofs << "# ControlBoard = pygame.Surface(Screen_Size, pygame.SRCALPHA)" << endl;
    ofs << "# ControlBoard.fill((0, 0, 0, 255))" << endl;
    ofs << endl;
    ofs << "Title = 'Vis2D'" << endl;
    ofs << "pygame.display.set_caption(Title)" << endl;
    ofs << endl;
    ofs << "Font_Sans = pygame.font.Font('freesansbold.ttf', 20)" << endl;
    ofs << "Font_Monospace = pygame.font.SysFont('monospace', 18, True)" << endl;
    ofs << "Font_Radar = pygame.font.SysFont('arial', 11)" << endl;
    ofs << endl;

    ofs << "class FEnvironment:" << endl;
    ofs << in+ "def __init__(self, InX=W_S/2, InY=H_S/2, InShape='circle', InRadius=W_S*0.5, InThickness=10):" << endl;
    ofs << in+ in+ "self.X = InX" << endl;
    ofs << in+ in+ "self.Y = InY" << endl;
    ofs << in+ in+ "self.Shape = InShape" << endl;
    ofs << in+ in+ "self.Radius = int(InRadius)" << endl;
    ofs << in+ in+ "self.Thickness = InThickness" << endl;
    ofs << in+ in+ "self.TransparentCircleArea = None" << endl;
    ofs << endl;
    ofs << in+ "def Draw(self, shape=None):" << endl;
    ofs << in+ in+ "if shape == 'circle':" << endl;
    ofs << in+ in+ in+ "pygame.draw.circle(Screen, YELLOW_FAINT, (self.X, self.Y), self.Radius)" << endl;
    ofs << in+ in+ "elif shape == 'lining':" << endl;
    ofs << in+ in+ in+ "# pygame.draw.circle(Screen, YELLOW_FAINT, (self.X, self.Y), self.Radius)" << endl;
    ofs << in+ in+ in+ "pygame.draw.circle(Screen, GRAY4, (self.X, self.Y), self.Radius, self.Thickness)" << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "pygame.draw.rect(Screen, YELLOW_FAINT, ((0, 0), (W_S, H_S)))" << endl;
    ofs << endl;
    ofs << in+ "def DrawTransparentArea(self):" << endl;
    ofs << in+ in+ "self.TransparentCircleArea = pygame.Surface((self.Radius * 2, self.Radius * 2), pygame.SRCALPHA)"
        << endl;
    ofs << in+ in+ "self.TransparentCircleArea.fill((255, 255, 255, 255))" << endl;
    ofs << in+ in+
           "pygame.draw.circle(self.TransparentCircleArea, (0, 0, 0, 0), (self.Radius, self.Radius), self.Radius)"
        << endl;
    ofs << endl;
    ofs << "class FOrganism:" << endl;
    ofs << in+ "def __init__(self, InName, InSpecies):" << endl;
    ofs << in+ in+ "self.Name = InName" << endl;
    ofs << in+ in+ "self.Species = InSpecies" << endl;
    ofs << in+ in+ "self.X = None" << endl;
    ofs << in+ in+ "self.Y = None" << endl;
    ofs << in+ in+ "self.Angle = None" << endl;
    ofs << endl;
    ofs << in+ in+ "# Draw" << endl;
    ofs << in+ in+ "self.BodyLength = 20" << endl;
    ofs << in+ in+ "self.BodyThickness = 10" << endl;
    ofs << in+ in+ "self.FlagellaLength_Fold = 2" << endl;
    ofs << in+ in+ "self.FlagellaThickness = 3" << endl;
    ofs << endl;
    ofs << in+ in+ "# Homeostasis Molecule Display" << endl;
    ofs << in+ in+ "self.HomeostasisMolecule_Switch = False" << endl;
    ofs << in+ in+ "self.HomeostasisMolecule_Colors = list()" << endl;
    ofs << endl;
    ofs << in+ in+ "# Radar" << endl;
    ofs << in+ in+ "self.Radar_Switch = False" << endl;
    ofs << in+ in+ "self.Radar_Sampling = 0" << endl;
    ofs << in+ in+ "self.Radar_Spacing = 0" << endl;
    ofs << endl;
    ofs << in+ in+ "# Memory for display & debugging" << endl;
    ofs << in+ in+ "self.Ligand_Prev = 0" << endl;
    ofs << in+ in+ "self.Am = 0" << endl; // TODO: hardcoded
    ofs << in+ in+ "self.SimCount = 0" << endl;
    ofs << endl;
    ofs << in+ in+ "# Trajectory" << endl;
    ofs << in+ in+ "self.TrajectorySwitch = True" << endl;
    ofs << in+ in+ "self.Trajectory = dict()" << endl;
    ofs << in+ in+ "self.TrajectoryColor = list()" << endl;
    ofs << endl;
    ofs << in+ "def Draw(self):" << endl;
    ofs << in+ in+ "if self.Species == 'Ecoli':" << endl;
    ofs << in+ in+ in+ "dX =  np.cos(self.Angle) * -self.BodyLength" << endl;
    ofs << in+ in+ in+ "dY = -np.sin(self.Angle) * -self.BodyLength" << endl;
    ofs << in+ in+ in+ "X_BodyEnd = self.X + dX" << endl;
    ofs << in+ in+ in+ "Y_BodyEnd = self.Y + dY" << endl;
    ofs << in+ in+ in+ "X_TailEnd = self.X + self.FlagellaLength_Fold * dX" << endl;
    ofs << in+ in+ in+ "Y_TailEnd = self.Y + self.FlagellaLength_Fold * dY" << endl;
    ofs << in+ in+ in+ "Threshold = None" << endl;
    ofs << in+ in+ in+ "Current = None" << endl;
    ofs << in+ in+ in+ "bChemotaxis = None" << endl;
    ofs << in+ in+ in+ "if self.HomeostasisMolecule_Switch:" << endl;
    ofs << in+ in+ in+ in+ "Threshold = SimM.Debug_ApplyUnit(SimM.State.Pos_Threshold)" << endl;
    ofs << in+ in+ in+ in+ "Current = SimM.Debug_ApplyUnit(SimM.GetCount(SimM.Idx_Count_Homeostasis).transpose())" << endl;
    ofs << in+ in+ in+ in+ "bChemotaxis = SimM.EvaluateChemotaxisThreshold()" << endl;
    ofs << endl;
    ofs << in+ in+ in+ "for i in range(self.X.size):" << endl;
    ofs << in+ in+ in+ in+ "if i == 0:" << endl;
    ofs << in+ in+ in+ in+ in+ "Color = RED" << endl;
    ofs << in+ in+ in+ in+ "else:" << endl;
    ofs << in+ in+ in+ in+ in+ "Color = YELLOW" << endl;
    ofs << in+ in+ in+ in+ "pygame.draw.line(Screen, Color, (self.X[i], self.Y[i]), (X_BodyEnd[i], Y_BodyEnd[i]), self.BodyThickness)" << endl;
    ofs << in+ in+ in+ in+ "pygame.draw.line(Screen, Color, (self.X[i], self.Y[i]), (X_TailEnd[i], Y_TailEnd[i]), self.FlagellaThickness)" << endl;
    ofs << in+ in+ in+ in+ "if self.Radar_Switch:" << endl;
    ofs << in+ in+ in+ in+ in+ "self.DrawRadar(Color, self.X[i], self.Y[i])" << endl;
    ofs << in+ in+ in+ in+ "if self.HomeostasisMolecule_Switch:" << endl;
    ofs << in+ in+ in+ in+ in+ "self.DisplayHomeostasisMolecules(self.X[i], self.Y[i], Threshold[:, i], Current[:, i], bChemotaxis[:, i])" << endl;
    ofs << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "pygame.draw.circle(Screen, BLACK, (self.X, self.Y), 5)" << endl;
    ofs << endl;
    ofs << in+ "def DrawRadar(self, Color, X, Y):" << endl;
    ofs << in+ in+ "for i in range(len(self.Radar_MolList)):" << endl;
    ofs << in+ in+ in+ "dY = 0" << endl;
    ofs << in+ in+ in+ "if len(self.Radar_MolList) > 1:" << endl;
    ofs << in+ in+ in+ in+ "dY = 16 * i - 8   # Hardcoded pattern for two molecules" << endl;
    ofs << in+ in+ in+ "self.DisplayRadarText(self.Radar_MolList[i], X, Y + dY, color=self.Radar_MolColor[i])" << endl;
    ofs << in+ in+ in+ "for j in range(self.Radar_Sampling):" << endl;
    ofs << in+ in+ in+ in+ "Spacing = self.Radar_Spacing * (j + 1)" << endl;
    ofs << in+ in+ in+ in+ "pygame.draw.circle(Screen, Color, (X, Y), Spacing, 1)" << endl;
    ofs << in+ in+ in+ in+ "RotationAngle = int(i * 90 / len(self.Radar_MolList))" << endl;
    ofs << in+ in+ in+ in+ "self.DisplayRadarTexts(self.Radar_MolList[i], X, Y, Spacing, rotate=RotationAngle, color=self.Radar_MolColor[i])" << endl;
    ofs << endl;
    ofs << in+ "def DisplayHomeostasisMolecules(self, X, Y, Threshold, Current, bChemotaxis):" << endl;
    ofs << in+ in+ "if len(HomeostasisMolName) == 0:" << endl;
    ofs << in+ in+ in+ "return" << endl;
    ofs << endl;
    ofs << in+ in+ "# Determine coordinates" << endl;
    ofs << in+ in+ "Spacing = 20" << endl;
    ofs << in+ in+ "if self.Radar_Switch:" << endl;
    ofs << in+ in+ in+ "Spacing += self.Radar_Spacing * self.Radar_Sampling" << endl;
    ofs << in+ in+ "X_Category = X - 30" << endl;
    ofs << in+ in+ "X_Mols = [X, X + 50]" << endl;
    ofs << in+ in+ "Y_Legends = Y - Spacing - 60" << endl;
    ofs << in+ in+ "Y_Threshold = Y - Spacing - 40" << endl;
    ofs << in+ in+ "Y_Current = Y - Spacing - 20" << endl;
    ofs << in+ in+ "Y_Chemotaxis = Y - Spacing" << endl;
    ofs << in+ in+ "if Y_Legends < 0:" << endl;
    ofs << in+ in+ in+ "Y_Legends = Y + Spacing" << endl;
    ofs << in+ in+ in+ "Y_Threshold = Y + Spacing + 20" << endl;
    ofs << in+ in+ in+ "Y_Current = Y + Spacing + 40" << endl;
    ofs << in+ in+ in+ "Y_Chemotaxis = Y + Spacing + 60" << endl;
    ofs << endl;
    ofs << in+ in+ "# Hardcoded for one and two threshold cases" << endl;
    ofs << in+ in+ "# Categories" << endl;
    ofs << in+ in+ "self.DisplayString('Threshold:', X_Category, Y_Threshold, position='midright', color=BLACK)" << endl;
    ofs << in+ in+ "self.DisplayString('Current:', X_Category, Y_Current, position='midright', color=BLACK)" << endl;
    ofs << in+ in+ "for i in range(len(HomeostasisMolName)):" << endl;
    // name
    ofs << in+ in+ in+ "Text = HomeostasisMolName[i] + ' (' + UnitTxt + ')'" << endl;
    ofs << in+ in+ in+ "self.DisplayString(Text, X_Mols[i], Y_Legends, color=self.Radar_MolColor[i])" << endl;
    // threshold
    ofs << in+ in+ in+ "Text = self.FormatValueToStr(Threshold[i])" << endl;
    ofs << in+ in+ in+ "self.DisplayString(Text, X_Mols[i], Y_Threshold, color=self.Radar_MolColor[i])" << endl;
    // current
    ofs << in+ in+ in+ "Text = self.FormatValueToStr(Current[i])" << endl;
    ofs << in+ in+ in+ "self.DisplayString(Text, X_Mols[i], Y_Current, color=self.Radar_MolColor[i])" << endl;
    // chemotaxis result
    ofs << in+ in+ in+ "Text = 'RUN' if bChemotaxis[i] else 'TUMBLE'" << endl;
    ofs << in+ in+ in+ "Color = RED if bChemotaxis[i] else BLACK" << endl;
    ofs << in+ in+ in+ "self.DisplayString(Text, X_Mols[i], Y_Chemotaxis, color=Color)" << endl;

    ofs << endl;
    ofs << in+ "def DisplayRadarTexts(self, MolName, X, Y, Spacing, rotate=0, color=BLACK):" << endl;
    ofs << in+ in+ "if rotate == 0 or rotate == 90:" << endl;
    ofs << in+ in+ in+ "Up = Y - Spacing" << endl;
    ofs << in+ in+ in+ "Down = Y + Spacing" << endl;
    ofs << in+ in+ in+ "Left = X - Spacing" << endl;
    ofs << in+ in+ in+ "Right = X + Spacing" << endl;
    ofs << in+ in+ in+ "self.DisplayRadarText(MolName, X, Up, color=color)" << endl;
    ofs << in+ in+ in+ "self.DisplayRadarText(MolName, X, Down, color=color)" << endl;
    ofs << in+ in+ in+ "self.DisplayRadarText(MolName, Left, Y, color=color)" << endl;
    ofs << in+ in+ in+ "self.DisplayRadarText(MolName, Right, Y, color=color)" << endl;
    ofs << endl;
    ofs << in+ in+ "elif rotate > 0 and rotate < 90:" << endl;
    ofs << in+ in+ in+ "Spacing_Y = Spacing * np.cos(np.deg2rad(rotate))" << endl;
    ofs << in+ in+ in+ "Spacing_X = Spacing * np.sin(np.deg2rad(rotate))" << endl;
    ofs << in+ in+ in+ "Up = int(Y - Spacing_Y)" << endl;
    ofs << in+ in+ in+ "Down = int(Y + Spacing_Y)" << endl;
    ofs << in+ in+ in+ "Left = int(X - Spacing_X)" << endl;
    ofs << in+ in+ in+ "Right = int(X + Spacing_X)" << endl;
    ofs << in+ in+ in+ "self.DisplayRadarText(MolName, Left, Up, color=color)" << endl;
    ofs << in+ in+ in+ "self.DisplayRadarText(MolName, Right, Up, color=color)" << endl;
    ofs << in+ in+ in+ "self.DisplayRadarText(MolName, Left, Down, color=color)" << endl;
    ofs << in+ in+ in+ "self.DisplayRadarText(MolName, Right, Down, color=color)" << endl;
    ofs << endl;
    ofs << in+ "def DisplayRadarText(self, MolName, X, Y, unit=False, color=BLACK, position='center'):" << endl;
    ofs << in+ in+ "if X >= 0 and Y >= 0 and X < W_S and Y < H_S:" << endl;
    ofs << in+ in+ in+ "Value = SimM.Debug_ApplyUnit(SimM.GetCountFromDistributionByNameOfDistAndXY(MolName, np.array([X]), np.array(Y)))" << endl;
    ofs << in+ in+ in+ "Value_Str = self.FormatValueToStr(Value[0])" << endl;
    ofs << in+ in+ in+ "if unit:" << endl;
    ofs << in+ in+ in+ in+ "Value_Str += (' ' + UnitTxt)" << endl;
    ofs << in+ in+ in+ "self.DisplayString(Value_Str, X, Y, color=color, position=position)" << endl;
    ofs << endl;
    ofs << in+ "def FormatValueToStr(self, Value):" << endl;
    ofs << in+ in+ "if Value < 0.01:" << endl;
    ofs << in+ in+ in+ "return SimF.SciFloat(Value, InPrecision=2, InExp_digits=2)" << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "return '{:.3f}'.format(Value)" << endl;
    ofs << endl;
    ofs << in+ "def DisplayString(self, String, X, Y, position='center', rotate=0, color=BLACK):" << endl;
    ofs << in+ in+ "Text = Font_Radar.render(String, True, color)" << endl;
    ofs << in+ in+ "Text = pygame.transform.rotate(Text, rotate)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "if position == 'center':" << endl;
    ofs << in+ in+ in+ "Text_Rect.center = (X, Y)" << endl;
    ofs << in+ in+ "if position == 'midleft':" << endl;
    ofs << in+ in+ in+ "Text_Rect.midleft = (X, Y)" << endl;
    ofs << in+ in+ "elif position == 'midright':" << endl;
    ofs << in+ in+ in+ "Text_Rect.midright = (X, Y)" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "def SetPosition(self, Position):" << endl;
    ofs << in+ in+ "self.X = Position[0]" << endl;
    ofs << in+ in+ "self.Y = Position[1]" << endl;
    ofs << in+ in+ "self.Angle = Position[2]" << endl;
    ofs << in+ in+ "# Threshold value (Position[3]) is not used" << endl;
    ofs << endl;
    ofs << in+ "def SetRadar(self, switch=True, sampling=6, spacing=30):" << endl;
    ofs << in+ in+ "self.HomeostasisMolecule_Switch = switch" << endl;
    ofs << in+ in+ "self.HomeostasisMolecule_Color = [RED_DARK, GREEN_DARK]" << endl;
    ofs << in+ in+ "self.Radar_Switch = switch" << endl;
    ofs << in+ in+ "self.Radar_Sampling = sampling" << endl;
    ofs << in+ in+ "self.Radar_Spacing = spacing" << endl;

    // temporary hardcoding
    ofs << in+ in+ "self.Radar_MolList = [";
    for (int i = 0; i < MolLoc.size(); i++) { ofs << "'" << MolLoc[i]->Name << "', "; }
    ofs << "]" << endl;
    
    std::vector<std::string> Radar_MolColor = {"BLACK", "BLUE"};
    ofs << in+ in+ "self.Radar_MolColor = [";
    for (int i = 0; i < MolLoc.size(); i++) { ofs << Radar_MolColor[i] << ", "; }
    ofs << "]" << endl;
    
    ofs << endl;
    ofs << in+ "def Receptivity(self, N_SimulationsToPass=50):" << endl;
    ofs << in+ in+ "for _ in range(N_SimulationsToPass):" << endl;
    ofs << in+ in+ in+ "SimM.SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution()" << endl;
    ofs << in+ in+ "SimM.SimLoop_WithSpatialSimulation()" << endl;
    ofs << endl;
    ofs << in+ "def ReportStatus(self):" << endl;
    ofs << in+ in+ "# for debugging" << endl;
    ofs << in+ in+ "self.Am = SimM.GetCountByName(HomeostasisMolName[0])" << endl;
//    ofs << in+ in+ "self.qAm = SimM.GetCountByName(HomeostasisMolName[1])" << endl;
    ofs << in+ in+ "self.SimCount += 1" << endl;
    ofs << in+ in+ "Ligand_Now = SimM.GetCountFromDistributionByNamesOfDistAndPos(GlucoseName, self.Name)[0]" << endl; // TODO: HARDCODED
    ofs << in+ in+ "SimStep = SimM.GetSimStep()" << endl;
    ofs << in+ in+ "Delta = 0" << endl;
    ofs << in+ in+ "if Ligand_Now != 0:" << endl;
    ofs << in+ in+ in+ "Delta = (Ligand_Now - self.Ligand_Prev[0]) / Ligand_Now * 100" << endl;
    ofs << in+ in+ "# print('SimStep {:06d} [Chemotaxis  {:06d}] Ligand:{:.6f} {} ({}{:.4f}%) Am:{:.6f} {} (X:{:.2f}, Y:{:.2f}, {:3.1f} degree)'.format" << endl;
    ofs << in+ in+ "#in+ '   (SimStep, self.SimCount, Ligand_Now / Unit / NA, UnitTxt, ('+' if Delta >= 0 else ''), Delta, self.Am / Unit / NA, UnitTxt , self.X, self.Y, self.Angle / pi * 180))" << endl;
    ofs << endl;
    ofs << in+ "def Homeostasis(self, MolName=[]):" << endl;
    ofs << in+ in+ "SimM.Homeostasis(MolName)   # Input 'Am', 'qAm' here" << endl;
    ofs << in+ in+ "self.Ligand_Prev = SimM.GetCountFromDistributionByNamesOfDistAndPos(GlucoseName, self.Name)" << endl; // TODO: HARDCODED
    ofs << in+ in+ "self.SetPosition(SimM.GetPositionXYAngleByName(self.Name))" << endl;
    ofs << in+ in+ "self.InitializeTrajectory()" << endl;
    ofs << in+ in+ "self.ReportStatus()" << endl;
    ofs << endl;
    ofs << in+ "def InitializeTrajectory(self):" << endl;
    ofs << in+ in+ "for i in range(self.X.size):" << endl;
    ofs << in+ in+ in+ "self.Trajectory[i] = [(self.X[i], self.Y[i])]" << endl;
    ofs << in+ in+ in+ "self.TrajectoryColor.append(tuple(np.random.randint(0, 255, 3)))" << endl;
    ofs << in+ in+ "self.TrajectoryColor[0] = MAGENTA" << endl;
    ofs << endl;
    ofs << in+ "def AddToTrajectory(self):" << endl;
    ofs << in+ in+ "for i in range(self.X.size):" << endl;
    ofs << in+ in+ in+ "self.Trajectory[i].append((self.X[i], self.Y[i]))" << endl;
    ofs << endl;
    ofs << in+ "def DrawTrajectory(self):" << endl;
    ofs << in+ in+ "for i in range(len(self.Trajectory)):" << endl;
    ofs << in+ in+ in+ "pygame.draw.aalines(Screen, self.TrajectoryColor[i], False, self.Trajectory[i])" << endl;
    ofs << endl;
    ofs << "class FMolecule:" << endl;
    ofs << in+ "def __init__(self, InName, InX, InY):" << endl;
    ofs << in+ in+ "self.Name = InName" << endl;
    ofs << in+ in+ "self.X = InX" << endl;
    ofs << in+ in+ "self.Y = InY" << endl;
    ofs << in+ in+ "self.DiffusionFactor = 2" << endl;
    ofs << in+ in+ "self.SpaceFactor = 50" << endl;
    ofs << in+ in+ "self.Color = ''" << endl;
    ofs << in+ in+ "self.Pattern = ''" << endl;
    ofs << in+ in+ "self.MaxBrightness = None" << endl;
    ofs << in+ in+ "self.NormalizationType = ''" << endl;
    ofs << endl;
    ofs << in+ in+ "# Heatmap Drawing" << endl;
    ofs << in+ in+ "self.ReductionFactor = 5" << endl;
    ofs << in+ in+ "self.MaxAmount = 0" << endl;
    ofs << in+ in+ "self.bMaxStatic = False" << endl;
    ofs << in+ in+ "self.bContourLine = False" << endl;
    ofs << in+ in+ "self.ContourLinePoints = list()" << endl;
    ofs << in+ in+ "self.InitializeHeatmapMax()" << endl;
    ofs << endl;
    ofs << in+ in+ "# Particle Drawing" << endl;
    ofs << in+ in+ "self.Particle_N = 300" << endl;
    ofs << in+ in+ "self.Particle_PerLayer = 3" << endl;
    ofs << in+ in+ "self.Particle_Radius = 2" << endl;
    ofs << in+ in+ "self.Particle_SpreadFactor = 1.2" << endl;
    ofs << in+ in+ "self.Particle_XY_Static = []" << endl;
    ofs << in+ in+ "self.InitializeStaticParticles()" << endl;
    ofs << endl;
    ofs << in+ "def InitializeHeatmapMax(self):" << endl;
    ofs << in+ in+ "self.MaxAmount = np.max(SimM.GetDistributionByName(self.Name))" << endl;
    ofs << endl;
    ofs << in+ "def InitializeStaticParticles(self):" << endl;
    ofs << in+ in+ "for i in range(int(self.Particle_N / self.Particle_PerLayer)):" << endl;
    ofs << in+ in+ in+ "for j in range(self.Particle_PerLayer):" << endl;
    ofs << in+ in+ in+ in+ "X = self.X + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))" << endl;
    ofs << in+ in+ in+ in+ "Y = self.Y + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))" << endl;
    ofs << in+ in+ in+ in+ "self.Particle_XY_Static.append((X, Y))" << endl;
    ofs << endl;
    ofs << in+ "def Reposition(self):" << endl;
    ofs << in+ in+ "self.X = np.random.randint(W_S * 2 / 5, W_S * 3 / 5)" << endl;
    ofs << in+ in+ "self.Y = np.random.randint(H_S * 2 / 5, H_S * 3 / 5)" << endl;
    ofs << in+ in+ "self.Particle_XY_Static = []" << endl;
    ofs << in+ in+ "self.InitializeStaticParticles()" << endl;
    ofs << endl;
    ofs << in+ "def Move(self, dX, dY):" << endl;
    ofs << in+ in+ "New_Particle_XY_Static = []" << endl;
    ofs << in+ in+ "for (X, Y) in self.Particle_XY_Static:" << endl;
    ofs << in+ in+ in+ "X += dX" << endl;
    ofs << in+ in+ in+ "Y += dY" << endl;
    ofs << in+ in+ in+ "New_Particle_XY_Static.append((X, Y))" << endl;
    ofs << in+ in+ "self.Particle_XY_Static = New_Particle_XY_Static" << endl;
    ofs << endl;
    ofs << in+ "def Draw(self, Data, threshold=0):" << endl;
    ofs << in+ in+ "if np.max(Data) == 0:" << endl;
    ofs << in+ in+ in+ "return"<< endl;
    ofs << endl;
    ofs << in+ in+ "if self.Pattern == 'spots':" << endl;
    ofs << in+ in+ in+ "Data_Normalized = SimF.Normalize_Linear(Data)" << endl;
    ofs << in+ in+ in+ "CoordsToDraw = np.where(Data_Normalized > threshold)" << endl;
    ofs << in+ in+ in+ "for X, Y, Value in zip(CoordsToDraw[0], CoordsToDraw[1], Data_Normalized[CoordsToDraw]):" << endl;
    ofs << in+ in+ in+ in+ "intensity = np.floor(Value * self.MaxBrightness)" << endl;
    ofs << in+ in+ in+ in+ "if intensity > self.MaxBrightness:" << endl;
    ofs << in+ in+ in+ in+ in+ "intensity = self.MaxBrightness" << endl;
    ofs << in+ in+ in+ in+ "color = self.GetColor(intensity)"<< endl;
    ofs << in+ in+ in+ in+ "pygame.draw.circle(Screen, color, (X, Y), self.Particle_Radius)" << endl;
    ofs << endl;
    ofs << in+ in+ "elif self.Pattern == 'heatmap':" << endl;
    ofs << in+ in+ in+ "Data_Normalized = None" << endl;
    ofs << in+ in+ in+ "ContourLine = list()" << endl;
    ofs << in+ in+ in+ "if self.NormalizationType == 'linear':" << endl;
    ofs << in+ in+ in+ in+ "Data_Normalized = SimF.Normalize_Linear(Data)" << endl;
    ofs << in+ in+ in+ "elif self.NormalizationType == 'log':" << endl;
    ofs << in+ in+ in+ in+ "Data_Normalized = SimF.Normalize_P1Log(Data)" << endl;
    ofs << endl;
    ofs << in+ in+ in+ "for x in range(0, Data_Normalized.shape[0], self.ReductionFactor):" << endl;
    ofs << in+ in+ in+ in+ "for y in range(0, Data_Normalized.shape[1], self.ReductionFactor):" << endl;
    ofs << in+ in+ in+ in+ in+ "PercentMolLevel = Data_Normalized[x][y]" << endl;
    ofs << in+ in+ in+ in+ in+ "if PercentMolLevel < threshold or PercentMolLevel == 0:" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "continue" << endl;
    ofs << in+ in+ in+ in+ in+ "intensity = np.floor(PercentMolLevel * self.MaxBrightness)" << endl;
    ofs << in+ in+ in+ in+ in+ "if intensity > self.MaxBrightness:" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "intensity = self.MaxBrightness" << endl;
    ofs << in+ in+ in+ in+ in+ "color = self.GetColor(intensity)"<< endl;
    ofs << in+ in+ in+ in+ in+ "pygame.draw.rect(Screen, color, ((x, y), (self.ReductionFactor, self.ReductionFactor)))" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ in+ "if self.bContourLine:" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "if not self.bMaxStatic:" << endl;
    ofs << in+ in+ in+ in+ in+ in+ in+ "self.SetContourLinePoints(Max=np.Max(Data))" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "if self.CheckContourLine(Data[x][y]):" << endl;
    ofs << in+ in+ in+ in+ in+ in+ in+ "pygame.draw.rect(Screen, BLACK, ((x, y), (self.ReductionFactor, self.ReductionFactor)))" << endl;
    ofs << endl;
    ofs << in+ in+ "elif self.Pattern == 'particle':" << endl;
    ofs << in+ in+ in+ "for XY in self.Particle_XY_Static:" << endl;
    ofs << in+ in+ in+ in+ "pygame.draw.circle(Screen, BLUE, XY, self.Particle_Radius)" << endl;
    ofs << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "assert True, 'Unsupported molecule distribution pattern for drawing: %s' % self.Pattern" << endl;
    ofs << endl;

    ofs << in+ "def CheckContourLine(self, Value):" << endl;
    ofs << in+ in+ "for (Min, Max) in self.ContourLinePoints:" << endl;
    ofs << in+ in+ in+ "if Value > Min and Value < Max:" << endl;
    ofs << in+ in+ in+ in+ "return True" << endl;
    ofs << in+ in+ "return False" << endl;
    ofs << endl;

    ofs << in+ "def GetColor(self, Intensity):" << endl;
    ofs << in+ in+ "if self.Color == 'Red':" << endl;
    ofs << in+ in+ in+ "return (self.MaxBrightness, self.MaxBrightness - Intensity, self.MaxBrightness - Intensity)" << endl;
    ofs << in+ in+ "elif self.Color == 'Green':" << endl;
    ofs << in+ in+ in+ "return (self.MaxBrightness - Intensity, self.MaxBrightness, self.MaxBrightness - Intensity)" << endl;
    ofs << in+ in+ "elif self.Color == 'Blue':" << endl;
    ofs << in+ in+ in+ "return (self.MaxBrightness - Intensity, self.MaxBrightness - Intensity, self.MaxBrightness)" << endl;
    ofs << in+ in+ "elif self.Color == 'Yellow':" << endl;
    ofs << in+ in+ in+ "return (self.MaxBrightness, self.MaxBrightness, self.MaxBrightness - Intensity)" << endl;
    ofs << in+ in+ "elif self.Color == 'Cyan':" << endl;
    ofs << in+ in+ in+ "return (self.MaxBrightness - Intensity, self.MaxBrightness, self.MaxBrightness)" << endl;
    ofs << in+ in+ "elif self.Color == 'Magenta':" << endl;
    ofs << in+ in+ in+ "return (self.MaxBrightness, self.MaxBrightness - Intensity, self.MaxBrightness)" << endl;
    ofs << endl;

    ofs << in+ "def SetColor(self, Color, MaxBrightness):" << endl;
    ofs << in+ in+ "self.Color = Color" << endl;
    ofs << in+ in+ "self.MaxBrightness = MaxBrightness" << endl;
    ofs << endl;

    ofs << in+ "def SetPattern(self, Pattern, NormalizationType, ReductionFactor, maxstatic=False, contourline=False):" << endl;
    ofs << in+ in+ "self.Pattern = Pattern" << endl;
    ofs << in+ in+ "self.NormalizationType = NormalizationType" << endl;
    ofs << in+ in+ "self.ReductionFactor = ReductionFactor" << endl;
    ofs << in+ in+ "self.bMaxStatic = maxstatic" << endl;
    ofs << in+ in+ "self.bContourLine = contourline" << endl;
    ofs << in+ in+ "if self.bContourLine:" << endl;
    ofs << in+ in+ in+ "self.SetContourLinePoints()" << endl;
    ofs << endl;

    ofs << in+ "def SetContourLinePoints(self, NumberOfPoints=7, Max=None):" << endl;
    ofs << in+ in+ "if not Max:" << endl;
    ofs << in+ in+ in+ "Max = self.MaxAmount" << endl;
    ofs << in+ in+ "for i in range(NumberOfPoints):" << endl;
    ofs << in+ in+ in+ "Max /= (1.005 * (i + 1))" << endl;
    ofs << in+ in+ in+ "self.ContourLinePoints.append((Max * 0.99, Max * 1.001))" << endl;
    ofs << endl;

    ofs << "class FControl:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "# self.FPS = 30" << endl;
    ofs << in+ in+ "self.MovementResolution = 2" << endl;
    ofs << in+ in+ "self.MessageWelcome = 'Welcome to Bacterial Chemotaxis!'" << endl;
    ofs << in+ in+ "self.Pos_Welcome = GetMidPoint(Center, CenterTop)" << endl;
    ofs << in+ in+ "self.Message = ''" << endl;
    ofs << in+ in+ "self.MessageTimer = 3000" << endl;
    ofs << endl;
    ofs << in+ in+ "self.InstructionSwitch = False" << endl;
    ofs << in+ in+ "self.Instructions = {" << endl;
    ofs << in+ in+ in+ "'T' : 'Display Trajectory Switch'," << endl;
    ofs << in+ in+ in+ "'I' : 'Display Instruction Switch'," << endl;
//    ofs << in+ in+ in+ "'S' : 'Display Score Switch'," << endl;
    ofs << in+ in+ in+ "'A' : 'Display Status Switch'," << endl;
    ofs << in+ in+ in+ "'R' : 'Display Radar Switch'," << endl;
    ofs << in+ in+ in+ "'H' : 'Display Homeostasis Switch'," << endl;
    ofs << in+ in+ in+ "'P' : 'Pause Visualization'," << endl;
    ofs << in+ in+ "}" << endl;
    ofs << in+ in+ "self.InstructionText = ''" << endl;
    ofs << in+ in+ "self.SetInstructionText()" << endl;
    ofs << endl;
    ofs << in+ in+ "self.MoleculeGradientText = ''" << endl;
    ofs << in+ in+ "self.MoleculeGradientColor = list()" << endl;
    ofs << endl;
//    ofs << in+ in+ "self.Score = 0" << endl;
//    ofs << in+ in+ "self.ScoreSwitch = False" << endl;
//    ofs << endl;
    ofs << in+ in+ "self.Time = 0" << endl;
    ofs << endl;
    ofs << in+ in+ "# DK - debugging purposes" << endl;
    ofs << in+ in+ "self.StatusSwitch = True" << endl;
    ofs << in+ in+ "# self.StatusSwitch = False" << endl;
    ofs << endl;
    ofs << in+ in+ "# Pause" << endl;
    ofs << in+ in+ "self.MessagePause = 'PAUSE'" << endl;
    ofs << in+ in+ "self.PauseSwitch = False" << endl;
    ofs << endl;
    ofs << in+ "def SetInstructionText(self):" << endl;
    ofs << in+ in+ "self.InstructionText = 'Instructions \\" << "n'" << endl;
    ofs << in+ in+ "for Key, Value in self.Instructions.items():" << endl;
    ofs << in+ in+ in+ "Space = ' ' * (6 - len(Key))" << endl;
    ofs << in+ in+ in+ "self.InstructionText = self.InstructionText + '  ' + Key + Space + ': ' + Value + '\\" << "n'" << endl;
    ofs << endl;
    ofs << in+ "def DisplayInstructions(self):" << endl;
    ofs << in+ in+ "TextLines = self.InstructionText.splitlines()" << endl;
    ofs << in+ in+ "Height = Font_Monospace.get_linesize() + 2" << endl;
    ofs << in+ in+ "X, Y = Screen.get_rect().topleft" << endl;
    ofs << in+ in+ "Color = None" << endl;
    ofs << in+ in+ "for i, TextLine in enumerate(TextLines):" << endl;
    ofs << in+ in+ in+ "Color = BLACK" << endl;
    ofs << in+ in+ in+ "Text = Font_Monospace.render(TextLine, True, Color)" << endl;
    ofs << in+ in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ in+ "Text_Rect.topleft = (X, Y + Height * i)" << endl;
    ofs << in+ in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "def DisplayWelcome(self):" << endl;
    ofs << in+ in+ "Text = Font_Sans.render(self.MessageWelcome, True, BLACK, WHITE)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "Text_Rect.center = self.Pos_Welcome" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "def DisplayInput(self):" << endl;
    ofs << in+ in+ "Text = Font_Sans.render(self.Message, True, BLACK)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "Text_Rect.bottomleft = Screen.get_rect().bottomleft" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << in+ in+ "self.MessageTimer -= 1" << endl;
    ofs << endl;
    ofs << in+ "def SetMessage(self, Key):" << endl;
    ofs << in+ in+ "assert Key in self.Instructions" << endl;
    ofs << in+ in+ "return 'Input [' + Key + '] : ' + self.Instructions[Key] + '   '" << endl;
    ofs << endl;
//    ofs << in+ "def DisplayScore(self):" << endl;
//    ofs << in+ in+ "Text = Font_Sans.render('Score: ' + str(round(self.Score)), True, RED)" << endl;
//    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
//    ofs << in+ in+ "Text_Rect.midtop = tuple(map(lambda i, j: i + j, Screen.get_rect().midtop, (0, Text.get_height())))" << endl;
//    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
//    ofs << endl;
    ofs << in+ "def DisplayTime(self):" << endl;
    ofs << in+ in+ "Text = Font_Sans.render('Simulation Time: ' + str(round(self.Time)), True, BLACK)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "Text_Rect.midtop = Screen.get_rect().midtop" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "# TODO: Update display status" << endl;
    ofs << in+ "# def DisplayStatus(self, Ligand_Total, Ligand_Local, Ligand_Prev_Local, Am):" << endl;
    ofs << in+ "def DisplayStatus(self, Ligand_Local, Ligand_Prev_Local, Am_Local):" << endl;
    ofs << in+ in+ "Ligand_Now = Ligand_Local[0]" << endl;
    ofs << in+ in+ "Ligand_Prev = Ligand_Prev_Local[0]" << endl;
    ofs << in+ in+ "Am = Am_Local[-1, 0]" << endl;
    ofs << in+ in+ "dLigand = 0" << endl;
    ofs << in+ in+ "if Ligand_Now != 0:" << endl;
    ofs << in+ in+ in+ "dLigand = (Ligand_Now - Ligand_Prev) / Ligand_Now * 100" << endl;
    ofs << endl;
    ofs << in+ in+ "StatusText = 'Ligand @ RED :' + '{:.5f} '.format(Ligand_Now/ Unit) + UnitTxt + '\\n' \\" << endl;
    ofs << in+ in+ in+ in+ in+ " + 'dLigand @ RED : ' + ('+' if dLigand >= 0 else '') + '{:.5f}'.format(dLigand) + ' % \\n' \\" << endl;
    ofs << in+ in+ in+ in+ in+ " + 'Am level of RED : ' + '{:.5f} '.format(Am / Unit) + UnitTxt   # Get the last E coli's info " << endl;
    ofs << endl;
    ofs << in+ in+ "TextLines = StatusText.splitlines()" << endl;
    ofs << in+ in+ "Height = Font_Monospace.get_linesize() + 2" << endl;
    ofs << in+ in+ "X, Y = Screen.get_rect().topright" << endl;
    ofs << in+ in+ "Color = BLACK" << endl;
    ofs << in+ in+ "for i, TextLine in enumerate(TextLines):" << endl;
    ofs << in+ in+ in+ "if 'RED' in TextLine:" << endl;
    ofs << in+ in+ in+ in+ "Color = RED" << endl;
    ofs << in+ in+ in+ "Text = Font_Monospace.render(TextLine, True, Color)" << endl;
    ofs << in+ in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ in+ "Text_Rect.topright = (X, Y + Height * i)" << endl;
    ofs << in+ in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "def DisplayPause(self):" << endl;
    ofs << in+ in+ "Text = Font_Sans.render(self.MessagePause, True, BLACK, WHITE)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "Text_Rect.center = self.Pos_Welcome" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << endl;
    ofs << "def main():" << endl;
    ofs << in+ "Control = FControl()" << endl;
    ofs << endl;
    ofs << in+ "SimUnitTime = 0.25" << endl;
    ofs << endl;
    ofs << in+ "PetriDish = FEnvironment()" << endl;
    ofs << endl;
    ofs << in+ "# TODO: Communicate to initialize in Sim" << endl;

    // Distribution settings (hardcoding)
    //                                              L,              qL
    std::vector<std::string>    Color               {"Yellow",      "Blue",     "Green"     };
    std::vector<int>            MaxBrightness       {200,           170,        50          };

    std::vector<std::string>    Pattern             {"heatmap",     "heatmap",  "particles" };
    std::vector<std::string>    NormalizationType   {"linear",      "log",      "particles" };
    std::vector<int>            ReductionFactor     {5,             2,          1           };

    std::vector<std::string>    MaxStatic;
    std::vector<std::string>    ContourLine;

//    if (Option.bDebug)
//    {
                                Color =             {"Green",       "Blue",     "Green"     };
//                                MaxStatic =         {"True",        "False",    "particles" };
//                                ContourLine =       {"True",        "False",    "False" };
//    }

    // Instantiate Molecules for Distribution
    // auto& MolLoc = Context.GetSubList_LocationList("Molecule");
    for (int i = 0; i < MolLoc.size(); i++) {
        // instantiate
        ofs << in+ MolLoc[i]->Name << " = FMolecule('" << MolLoc[i]->Name << "', ";
        ofs << MolLoc[i]->Coord[0] << ", " << MolLoc[i]->Coord[1] << ")" << endl;

        // set color
        ofs << in+ MolLoc[i]->Name << ".SetColor('" << Color[i] << "', " << MaxBrightness[i] << ")" << endl;

        // set pattern
        ofs << in+ MolLoc[i]->Name << ".SetPattern('" << Pattern[i] << "', '" << NormalizationType[i] << "', " << ReductionFactor[i];
        if (!MaxStatic.empty())     { ofs << ", maxstatic=" << MaxStatic[i]; }
        if (!ContourLine.empty())   { ofs << ", contourline=" << ContourLine[i]; }
        ofs << ",)" << endl;

        ofs << endl;
    }

    // Distribution settings (hardcoding)
    std::vector<std::string>    Radar_Switch;

//    if (Option.bDebug)
//    {
                                Radar_Switch =      {"True",};
//    }

    // Instantiate Organisms
    auto OrgNames = Context.GetUniqueNames_LocationList("Organism");
    for (int i = 0; i < OrgNames.size(); i++) {
        ofs << in+ OrgNames[i] << " = FOrganism('" << OrgNames[i] << "', " << "'Ecoli'" << ")" << endl; // TODO: Get Species later

        // set radar
        if (!Radar_Switch.empty())     { ofs << in+ OrgNames[i] << ".SetRadar(switch=" << Radar_Switch[i] << ")" << endl; }
        ofs << endl;
    }

    // Run Homeostasis for "HomeostasisMolName"
    for (auto& OrganismName : OrgNames) {
        ofs << in+ OrganismName << ".Homeostasis(HomeostasisMolName)" << endl; // TODO: HARDCODED
    }
    ofs << endl;
    ofs << in+ "ElapsedTime = 0" << endl;
    ofs << in+ "PrevTime = datetime.now()" << endl;
    ofs << endl;
    ofs << in+ "SimState = True" << endl;
    ofs << in+ "while SimState:" << endl;
    ofs << endl;
    ofs << in+ in+ "if not Control.PauseSwitch:" << endl;
    ofs << in+ in+ in+ "CurrTime = datetime.now()" << endl;
    ofs << in+ in+ in+ "ElapsedTime += (CurrTime - PrevTime).total_seconds()" << endl;
    ofs << in+ in+ in+ "PrevTime = CurrTime" << endl;
    ofs << endl;
    ofs << in+ in+ "for event in pygame.event.get():" << endl;
    ofs << in+ in+ in+ "if event.type == pygame.QUIT:" << endl;
    ofs << in+ in+ in+ in+ "SimState = False" << endl;
    ofs << in+ in+ in+ "elif event.type == pygame.KEYDOWN:" << endl;
    ofs << in+ in+ in+ in+ "Control.MessageTimer = 5000" << endl;
    ofs << in+ in+ in+ in+ "if event.key == pygame.K_x:" << endl;
    ofs << in+ in+ in+ in+ in+ "SimState = False" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "# Organism Control" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_t:" << endl;
    for (auto& OrganismName : OrgNames) {
        ofs << in+ in+ in+ in+ in+ OrganismName << ".TrajectorySwitch = not " << OrganismName << ".TrajectorySwitch" << endl;
    }
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('T')" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ "# Control panel" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_i:" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.InstructionSwitch = not Control.InstructionSwitch" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('I')" << endl;
//    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_s:" << endl;
//    ofs << in+ in+ in+ in+ in+ "Control.ScoreSwitch = not Control.ScoreSwitch" << endl;
//    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('S')" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_a:" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.StatusSwitch = not Control.StatusSwitch" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('A')" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_r:" << endl;
    for (auto& OrganismName : OrgNames) {
        ofs << in+ in+ in+ in+ in+ OrganismName << ".Radar_Switch = not " << OrganismName << ".Radar_Switch" << endl;
    }
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('R')" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_h:" << endl;
    for (auto& OrganismName : OrgNames) {
        ofs << in+ in+ in+ in+ in+ OrganismName << ".HomeostasisMolecule_Switch = not " << OrganismName << ".HomeostasisMolecule_Switch" << endl;
    }
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('H')" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_p:" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.PauseSwitch = not Control.PauseSwitch" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('P')" << endl;
    ofs << endl;
    ofs << in+ in+ "Screen.fill(GRAY1)" << endl;
    ofs << endl;
    ofs << in+ in+ "# PetriDish.Draw()" << endl;
    ofs << in+ in+ "PetriDish.Draw(shape='circle')" << endl;
    ofs << endl;

    for (auto& Mol : MolLoc) {
        ofs << in+ in+ Mol->Name << ".Draw(SimM.GetDistributionByName('" << Mol->Name << "'))" << endl;
    }
    ofs << endl;
    ofs << in+ in+ "if Control.PauseSwitch:" << endl;
    ofs << in+ in+ in+ "Control.DisplayPause()" << endl;
    ofs << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "while ElapsedTime >= SimUnitTime:" << endl;
    ofs << endl;

    for (auto& OrganismName : OrgNames) {
        ofs << in+ in+ in+ in+ OrganismName << ".Receptivity(50)" << endl;
        ofs << in+ in+ in+ in+ OrganismName << ".SetPosition(SimM.GetPositionXYAngleByName('" << OrganismName << "'))" << endl;
        ofs << in+ in+ in+ in+ OrganismName << ".ReportStatus()" << endl;
    }
    ofs << endl;

    ofs << in+ in+ in+ in+ "ElapsedTime -= SimUnitTime" << endl;
    ofs << in+ in+ in+ in+ "Control.Time += 1" << endl;
    ofs << endl;

    for (auto& OrganismName : OrgNames) {
        ofs << in+ in+ OrganismName << ".AddToTrajectory()" << endl;
        ofs << in+ in+ "if " << OrganismName << ".TrajectorySwitch:" << endl;
        ofs << in+ in+ in+ OrganismName << ".DrawTrajectory()" << endl;
    }
    ofs << endl;

    ofs << endl;

    for (auto& OrganismName : OrgNames) {
        ofs << in+ in+ OrganismName << ".Draw()" << endl;
    }
    ofs << in+ in+ "PetriDish.Draw(shape='lining')" << endl;
    ofs << endl;

    // Show first molecule only
    std::string MolName = MolLoc[0]->Name;
    std::string OrgName = OrgNames[0];

    ofs << in+ in+ MolName << "_Now = SimM.GetCountFromDistributionByNamesOfDistAndPos('" << MolName << "', '" << OrgName << "')" << endl;

    ofs << endl;

    if (!Option.bDebug) {
        ofs << in+ in+ "if Control.Time < 50:" << endl;
        ofs << in+ in+ in+ "Control.DisplayWelcome()" << endl;
    }
    ofs << in+ in+ "if Control.MessageTimer > 0:" << endl;
    ofs << in+ in+ in+ "Control.DisplayInput()" << endl;
    ofs << in+ in+ "if Control.InstructionSwitch:" << endl;
    ofs << in+ in+ in+ "Control.DisplayInstructions()" << endl;
//    ofs << in+ in+ "if Control.ScoreSwitch:" << endl;
//    ofs << in+ in+ in+ "Control.Score += " << MolName << "_Now / 10" << endl;
//    ofs << in+ in+ in+ "Control.DisplayScore()" << endl;
    ofs << endl;
    ofs << in+ in+ "if Control.StatusSwitch:" << endl;
    ofs << in+ in+ in+ "Control.DisplayStatus(" << MolName << "_Now, " << OrgName << ".Ligand_Prev, " << OrgName << ".Am)" << endl;
    ofs << endl;
    ofs << in+ in+ "# Update Ligand Prev" << endl;
    ofs << in+ in+ OrgName << ".Ligand_Prev = " << MolName << "_Now" << endl;
    ofs << endl;
    ofs << in+ in+ "Control.DisplayTime()" << endl;
    ofs << in+ in+ "# Control.DisplayMoleculeGradient()" << endl;
    ofs << endl;
    ofs << in+ in+ "pygame.display.update()" << endl;
    ofs << endl;
    ofs << in+ "pygame.quit()" << endl;
    ofs << in+ "sys.exit()" << endl;
    ofs << endl;
    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "main()" << endl;
    ofs << endl;

    std::cout << "  Visualization_2D program has been generated: ";
}
