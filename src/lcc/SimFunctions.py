import numpy as np
import math

NA = 6.0221409e+23

def ConcToCount(Conc_Molecule, Volume):
    return (Conc_Molecule * NA) * Volume

def CountToConc(Count_Molecule, Volume):
    return (Count_Molecule / NA) / Volume

def CheckRateAndConc(Rate, Conc_1, Conc_2, Conc_3):
    return np.min((Rate, Conc_1, Conc_2, Conc_3), axis=0)

# Standard reactions
def Eqn_Standard_Unregulated(Conc_1, Conc_2, Conc_3, k):
    return k * Conc_1 * Conc_2 * Conc_3

def Eqn_Standard_Inhibition_Allosteric(Conc_1, Conc_2, Conc_3, Conc_Inhibitor, k, Ki, n):
    return k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1 * Conc_2 * Conc_3

def Eqn_Standard_Activation_Allosteric(Conc_1, Conc_2, Conc_3, Conc_Activator, k, Ka, n):
    return k * (1 + (Conc_Activator / Ka) ** n) * Conc_1 * Conc_2 * Conc_3

# Enzymatic, standard reactions
def Eqn_Enz_Standard_Unregulated(Conc_Enzyme, Conc_1, Conc_2, Conc_3, k):
    return Conc_Enzyme * k * Conc_1 * Conc_2 * Conc_3
    # return Conc_Enzyme * k (for saturation)

def Eqn_Enz_Standard_Inhibition_Allosteric(Conc_Enzyme, Conc_1, Conc_2, Conc_3, Conc_Inhibitor, k, Ki, n):
    return Conc_Enzyme * k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_1 * Conc_2 * Conc_3

def Eqn_Enz_Standard_Activation_Allosteric(Conc_Enzyme, Conc_1, Conc_2, Conc_3, Conc_Activator, k, Ka, n):
    return Conc_Enzyme * k * (1 + (Conc_Activator / Ka) ** n) * Conc_1 * Conc_2 * Conc_3

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
    Dist = math.sqrt(((X - X_Ori) / Width) ** 2 + ((Y - Y_Ori) / Height) ** 2)
    return Max / max(1, Dist * 30)

def Eqn_Diffusion(D, Amount_Source, Amount_Neighbor):
    return D * (Amount_Neighbor - Amount_Source)

def Eqn_Diffusion_Spatial(Distribution, D):
    Row_Zero = np.zeros((1, Distribution.shape[1]))
    Column_Zero = np.zeros((Distribution.shape[0], 1))

    # expanded version for readability
    # Roll
    Upward = np.roll(Distribution, -1, axis=0)
    Downward = np.roll(Distribution, 1, axis=0)
    Leftward = np.roll(Distribution, -1, axis=1)
    Rightward = np.roll(Distribution, 1, axis=1)

    # Diffuse
    Diffusion_Upward = Eqn_Diffusion(D, Distribution, Upward)
    Diffusion_Downward = Eqn_Diffusion(D, Distribution, Downward)
    Diffusion_Leftward = Eqn_Diffusion(D, Distribution, Leftward)
    Diffusion_Rightward = Eqn_Diffusion(D, Distribution, Rightward)

    # Correction by Slicing and Adding back zeros for flows from null direction
    Upward_Corrected = np.concatenate((Diffusion_Upward[:-1, :], Row_Zero), axis=0)
    Downward_Corrected = np.concatenate((Row_Zero, Diffusion_Downward[1:, :]), axis=0)
    Leftward_Corrected = np.concatenate((Diffusion_Leftward[:, :-1], Column_Zero), axis=1)
    Rightward_Corrected = np.concatenate((Column_Zero, Diffusion_Rightward[:, 1:]), axis=1)

    return Upward_Corrected + Downward_Corrected + Leftward_Corrected + Rightward_Corrected

def DiffuseDistribution(Distribution, D=0.01, dTime=1):
    Distribution_Updated = Distribution + Eqn_Diffusion_Spatial(Distribution, D) * dTime
    return Distribution_Updated
