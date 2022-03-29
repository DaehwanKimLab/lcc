import numpy as np
import random

NA = 6.0221409e+23
pi = np.pi
random.seed(1)

def ConcToCount(Conc_Molecule, Volume):
    return (Conc_Molecule * NA) * Volume

def CountToConc(Count_Molecule, Volume):
    return (Count_Molecule / NA) / Volume

### if N_MOLECULEALLOWED = 1 ###########################################################################################
def CheckRateAndConc_1(Rate, Conc_1):
    return np.min((Rate, Conc_1), axis=0)

# Standard reactions
def Eqn_Standard_Unregulated_1(Conc_1, k):
    return k * Conc_1  

def Eqn_Standard_Inhibition_Allosteric_1(Conc_1, Conc_Inhibitor, k, Ki, n):
    return k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1  

def Eqn_Standard_Activation_Allosteric_1(Conc_1, Conc_Activator, k, Ka, n):
    return k * (1 + (Conc_Activator / Ka) ** n) * Conc_1  

# Enzymatic, standard reactions
def Eqn_Enz_Standard_Unregulated_1(Conc_Enzyme, Conc_1, k):
    return Conc_Enzyme * k * Conc_1  
    # return Conc_Enzyme * k (for saturation)

def Eqn_Enz_Standard_Inhibition_Allosteric_1(Conc_Enzyme, Conc_1, Conc_Inhibitor, k, Ki, n):
    return Conc_Enzyme * k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1  

def Eqn_Enz_Standard_Activation_Allosteric_1(Conc_Enzyme, Conc_1, Conc_Activator, k, Ka, n):
    return Conc_Enzyme * k * (1 + (Conc_Activator / Ka) ** n) * Conc_1  
########################################################################################################################

### if N_MOLECULEALLOWED = 2 ###########################################################################################
def CheckRateAndConc_2(Rate, Conc_1, Conc_2):
    return np.min((Rate, Conc_1, Conc_2), axis=0)

# Standard reactions
def Eqn_Standard_Unregulated_2(Conc_1, Conc_2, k):
    return k * Conc_1 * Conc_2 

def Eqn_Standard_Inhibition_Allosteric_2(Conc_1, Conc_2, Conc_Inhibitor, k, Ki, n):
    return k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1 * Conc_2 

def Eqn_Standard_Activation_Allosteric_2(Conc_1, Conc_2, Conc_Activator, k, Ka, n):
    return k * (1 + (Conc_Activator / Ka) ** n) * Conc_1 * Conc_2 

# Enzymatic, standard reactions
def Eqn_Enz_Standard_Unregulated_2(Conc_Enzyme, Conc_1, Conc_2, k):
    return Conc_Enzyme * k * Conc_1 * Conc_2 
    # return Conc_Enzyme * k (for saturation)

def Eqn_Enz_Standard_Inhibition_Allosteric_2(Conc_Enzyme, Conc_1, Conc_2, Conc_Inhibitor, k, Ki, n):
    return Conc_Enzyme * k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1 * Conc_2 

def Eqn_Enz_Standard_Activation_Allosteric_2(Conc_Enzyme, Conc_1, Conc_2, Conc_Activator, k, Ka, n):
    return Conc_Enzyme * k * (1 + (Conc_Activator / Ka) ** n) * Conc_1 * Conc_2 
########################################################################################################################

### if N_MOLECULEALLOWED = 3 ###########################################################################################
def CheckRateAndConc_3(Rate, Conc_1, Conc_2, Conc_3):
    return np.min((Rate, Conc_1, Conc_2, Conc_3), axis=0)

# Standard reactions
def Eqn_Standard_Unregulated_3(Conc_1, Conc_2, Conc_3, k):
    return k * Conc_1 * Conc_2 * Conc_3

def Eqn_Standard_Inhibition_Allosteric_3(Conc_1, Conc_2, Conc_3, Conc_Inhibitor, k, Ki, n):
    return k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1 * Conc_2 * Conc_3

def Eqn_Standard_Activation_Allosteric_3(Conc_1, Conc_2, Conc_3, Conc_Activator, k, Ka, n):
    return k * (1 + (Conc_Activator / Ka) ** n) * Conc_1 * Conc_2 * Conc_3

# Enzymatic, standard reactions
def Eqn_Enz_Standard_Unregulated_3(Conc_Enzyme, Conc_1, Conc_2, Conc_3, k):
    return Conc_Enzyme * k * Conc_1 * Conc_2 * Conc_3
    # return Conc_Enzyme * k (for saturation)

def Eqn_Enz_Standard_Inhibition_Allosteric_3(Conc_Enzyme, Conc_1, Conc_2, Conc_3, Conc_Inhibitor, k, Ki, n):
    return Conc_Enzyme * k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1 * Conc_2 * Conc_3

def Eqn_Enz_Standard_Activation_Allosteric_3(Conc_Enzyme, Conc_1, Conc_2, Conc_3, Conc_Activator, k, Ka, n):
    return Conc_Enzyme * k * (1 + (Conc_Activator / Ka) ** n) * Conc_1 * Conc_2 * Conc_3
########################################################################################################################

### if N_MOLECULEALLOWED = 4 ###########################################################################################
def CheckRateAndConc_4(Rate, Conc_1, Conc_2, Conc_3, Conc_4):
    return np.min((Rate, Conc_1, Conc_2, Conc_3, Conc_4), axis=0)

# Standard reactions
def Eqn_Standard_Unregulated_4(Conc_1, Conc_2, Conc_3, Conc_4, k):
    return k * Conc_1 * Conc_2 * Conc_3 * Conc_4

def Eqn_Standard_Inhibition_Allosteric_4(Conc_1, Conc_2, Conc_3, Conc_4, Conc_Inhibitor, k, Ki, n):
    return k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1 * Conc_2 * Conc_3 * Conc_4

def Eqn_Standard_Activation_Allosteric_4(Conc_1, Conc_2, Conc_3, Conc_4, Conc_Activator, k, Ka, n):
    return k * (1 + (Conc_Activator / Ka) ** n) * Conc_1 * Conc_2 * Conc_3 * Conc_4

# Enzymatic, standard reactions
def Eqn_Enz_Standard_Unregulated_4(Conc_Enzyme, Conc_1, Conc_2, Conc_3, Conc_4, k):
    return Conc_Enzyme * k * Conc_1 * Conc_2 * Conc_3 * Conc_4
    # return Conc_Enzyme * k (for saturation)

def Eqn_Enz_Standard_Inhibition_Allosteric_4(Conc_Enzyme, Conc_1, Conc_2, Conc_3, Conc_4, Conc_Inhibitor, k, Ki, n):
    return Conc_Enzyme * k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1 * Conc_2 * Conc_3 * Conc_4

def Eqn_Enz_Standard_Activation_Allosteric_4(Conc_Enzyme, Conc_1, Conc_2, Conc_3, Conc_4, Conc_Activator, k, Ka, n):
    return Conc_Enzyme * k * (1 + (Conc_Activator / Ka) ** n) * Conc_1 * Conc_2 * Conc_3 * Conc_4
########################################################################################################################

# Enzymatic, Michaelis Menten reactions
def Eqn_Enz_MichaelisMenten_Unregulated(Conc_Enzyme, Conc_Substrate, kcat, KM):
    return (kcat * Conc_Enzyme * Conc_Substrate) / (KM + Conc_Substrate)

def Eqn_Enz_MichaelisMenten_CompetitiveInhibition(Conc_Enzyme, Conc_Substrate, Conc_Inhibitor, kcat, KM, Ki):
    return (kcat * Conc_Enzyme * Conc_Substrate) / (KM * (1 + (Conc_Inhibitor / Ki)) + Conc_Substrate)

def Eqn_Enz_MichaelisMenten_Inhibition_Allosteric(Conc_Enzyme, Conc_Substrate, Conc_Inhibitor, kcat, KM, Ki, n):
    return (kcat * Conc_Enzyme / (1 + (Conc_Inhibitor / Ki) ** n)) * (Conc_Substrate / (KM + Conc_Substrate))

def Eqn_Enz_MichaelisMenten_Activation_Allosteric(Conc_Enzyme, Conc_Substrate, Conc_Activator, kcat, KM, Ka, n):
    return (kcat * Conc_Enzyme * (1 + (Conc_Activator / Ka) ** n)) * (Conc_Substrate / (KM + Conc_Substrate))


def MatrixMultiplication_Rev(Freq, Rate):
    return np.matmul(Rate, Freq)

def GetDerivativeFromStoichiometryMatrix(Freq, Rate):
    return MatrixMultiplication_Rev(Freq, Rate)

def DetermineAmountOfBuildingBlocks(Freq, Rate):
    return MatrixMultiplication_Rev(Freq, Rate)

def PickRandomIdx(Quantity, Indices, Weight=1):
    # Adjust Quantity and Weight if Weight is completely zero
    Sum_Weight = np.sum(Weight)
    Weight = Weight + np.where(Sum_Weight == 0, 1, 0)
    Quantity = Quantity * np.where(Sum_Weight == 0, 0, 1)

    # Generate cumulative sum on weight and pick a random number in its range
    Weight_Cumsum = np.cumsum(Weight)
    Weight_Cumsum_Min = Weight_Cumsum[0]
    Weight_Cumsum_Max = Weight_Cumsum[-1]
    Weight_Cumsum_Min = np.where(Weight_Cumsum_Min == Weight_Cumsum_Max, Weight_Cumsum_Min - 1, Weight_Cumsum_Min)
    RanNums = np.asmatrix(np.random.randint(Weight_Cumsum_Min, high=Weight_Cumsum_Max, size=Quantity)).transpose()
    # Generate a matrix of the random numbers for comparison to indices
    RanNums_Matrix = np.reshape(np.repeat(RanNums, Indices.shape[1]), [-1, Indices.shape[1]])
    Bin_RanNumLessThanWeightCumsum = np.where(RanNums_Matrix < Weight_Cumsum, 1, 0)
    Idx_Rnd = np.argmax(Bin_RanNumLessThanWeightCumsum, 1)
    return Idx_Rnd

def InsertZeroIntoNegOneElementInLenMatrix(Len, Indices):
    # Generate an array of counts for each index
    Count_Indices = np.zeros(Len.shape[1])
    np.put_along_axis(Count_Indices, Indices, 1, axis=0)
    # Generate a cumulative sum matrix of available position in the Len Matrix
    Bool_LenAvailable = np.less(Len, 0)  # used again later
    Bin_LenAvailable = Bool_LenAvailable.astype(int)
    LenCumsum = np.cumsum(Bin_LenAvailable, axis=0)
    # Get an overlap between availability and new count positions
    Bool_LenCumsumGreaterThanZero = np.greater(LenCumsum, 0)
    Bool_LenCumsumLessThanOrEqualToCountOfIndices = np.less_equal(LenCumsum, Count_Indices)
    Bool_LenCumsum = np.logical_and(Bool_LenCumsumGreaterThanZero, Bool_LenCumsumLessThanOrEqualToCountOfIndices)
    Bin_Len_Selected = np.logical_and(Bool_LenAvailable, Bool_LenCumsum).astype(int)
    return Len + Bin_Len_Selected


# Spatial simulation functions

def InitializeDistribution(Width, Height, X_Ori, Y_Ori, MaxAmount):
    Dist_Init = np.zeros((Width, Height))
    for X in range(Width):
        for Y in range(Height):
            Dist_Init[X][Y] = InitialDiffusionPattern(X, Y, Width, Height, X_Ori, Y_Ori, MaxAmount)
    return Dist_Init

def InitialDiffusionPattern(X, Y, Width, Height, X_Ori, Y_Ori, Max):
    Dist = np.sqrt(((X - X_Ori) / Width) ** 2 + ((Y - Y_Ori) / Height) ** 2)
    return Max / max(1, Dist * 30)

def Eqn_Transporter_Unregulated(Conc_Inside, Conc_Outside, ki, ko):
    return Conc_Outside * ko - Conc_Inside * ki

def Eqn_Diffusion_Spatial_4Cell(Distribution, D):
    # von Neumann neighborhood
    Diffusion_Quantity = Distribution * D
    Diffusion_Padded = np.pad(Diffusion_Quantity, (1, 1))

    # roll diffusion quantity
    Upward = np.roll(Diffusion_Padded, -1, axis=0)
    Downward = np.roll(Diffusion_Padded, 1, axis=0)
    Leftward = np.roll(Diffusion_Padded, -1, axis=1)
    Rightward = np.roll(Diffusion_Padded, 1, axis=1)

    # Flow back from the edges
    Upward[1, :] += Upward[0, :]
    Downward[-2, :] += Downward[-1, :]
    Leftward[:, 1] += Leftward[:, 0]
    Rightward[:, -2] += Rightward[:, -1]

    Distribution_Diffused = (Upward + Downward + Leftward + Rightward)[1:-1, 1:-1] - 4 * Diffusion_Quantity

    return Distribution_Diffused

def Eqn_Diffusion_Spatial_8Cell(Distribution, D):
    # Moore neighborhood
    Diffusion_Quantity = Distribution * D
    Diffusion_Padded = np.pad(Diffusion_Quantity, (1, 1))

    # roll diffusion quantity
    Upward = np.roll(Diffusion_Padded, -1, axis=0)
    Downward = np.roll(Diffusion_Padded, 1, axis=0)
    Leftward = np.roll(Diffusion_Padded, -1, axis=1)
    Rightward = np.roll(Diffusion_Padded, 1, axis=1)
    UpLeftward = np.roll(Upward, -1, axis=1)
    UpRightward = np.roll(Upward, 1, axis=1)
    DownLeftward = np.roll(Downward, -1, axis=1)
    DownRightward = np.roll(Downward, 1, axis=1)

    Distribution_Diffused = (Upward + Downward + Leftward + Rightward + UpLeftward + UpRightward + DownLeftward + DownRightward)[1:-1, 1:-1] - 8 * Diffusion_Quantity

    return Distribution_Diffused

def DiffuseDistribution_4Cell(Distribution, D=0.2, dTime=1):
    Distribution_Updated = Distribution + Eqn_Diffusion_Spatial_4Cell(Distribution, D) * dTime
    # return Distribution
    return Distribution_Updated

def DiffuseDistribution_8Cell(Distribution, D=0.2, dTime=1):
    Distribution_Updated = Distribution + Eqn_Diffusion_Spatial_8Cell(Distribution, D) * dTime
    # return Distribution
    return Distribution_Updated

def Displacement_2D(Distance, Angle):
    dX =  Distance * np.cos(Angle)
    dY = -Distance * np.sin(Angle)
    return dX, dY

# Temporary chemotaxis simulation functions. 
def BacterialChemotaxis(ThresholdSubject, X, Y, Angle, Threshold):
    # optimize the following code into additions and subtractions
    Decision = ThresholdSubject < Threshold
    X_Run, Y_Run, Angle_Run  = Chemotaxis_Run(X, Y, Angle)
    X_Tumble, Y_Tumble, Angle_Tumble = Chemotaxis_Tumble(X, Y, Angle)
    X_New = np.where(Decision, X_Run, X_Tumble)
    Y_New = np.where(Decision, Y_Run, Y_Tumble)
    Angle_New = np.where(Decision, Angle_Run, Angle_Tumble)
    return X_New, Y_New, Angle_New

def Chemotaxis_Run(X, Y, Angle, Distance = 15):
    dX, dY = Displacement_2D(Distance, Angle)
    return X + dX, Y + dY, Angle

def Chemotaxis_Tumble(X, Y, Angle, Distance = 5):
    dX, dY = Displacement_2D(Distance, Angle)
    NewAngle = np.random.random_sample(Angle.shape) * 2 * pi
    return X + dX, Y + dY, NewAngle

def CorrectOutOfBounds(X, Y, Width, Height):
    X_Corrected = np.where(X < np.array([0]), 0, X)
    X_Corrected = np.where(X_Corrected >= np.array([Width]), Width - 1, X_Corrected)
    Y_Corrected = np.where(Y < np.array([0]), 0, Y)
    Y_Corrected = np.where(Y_Corrected >= np.array([Height]), Height - 1, Y_Corrected)
    return X_Corrected, Y_Corrected

