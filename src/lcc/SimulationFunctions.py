import numpy as np


def ConcToCount(Conc_Molecule, Volume):
    return Conc_Molecule * Volume

def CountToConc(Count_Molecule, Volume):
    return Count_Molecule / Volume

def StandardEqn(Conc_Enzyme, Conc_Reactant_1, Conc_Reactant_2, Conc_Product_1, Conc_Product_2, k, krev):
    return Conc_Enzyme * ((k * Conc_Reactant_1 * Conc_Reactant_2) - (krev * Conc_Product_1 * Conc_Product_2))

# def StandardEqn(Conc_Enzyme, Conc_Reactant, Conc_Product, k, krev):
#     return Conc_Enzyme * ((k * Conc_Reactant) - (krev * Conc_Product))

def StandardEqn_Inhibition_Allosteric(Conc_Enzyme, Conc_Reactant, Conc_Product, Conc_Inhibitor, k, krev, Ki, n):
    return Conc_Enzyme * ((k / (1 + (Conc_Inhibitor / Ki) ** n) * Conc_Reactant) - (krev * Conc_Product))

def StandardEqn_Activation_Allosteric(Conc_Enzyme, Conc_Reactant, Conc_Product, Conc_Activator, k, krev, Ka, n):
    return Conc_Enzyme * ((k * (1 + (Conc_Activator / Ka) ** n) * Conc_Reactant) - (krev * Conc_Product))

def MichaelisMentenEqn(Conc_Enzyme, Conc_Substrate, kcat, KM):
    return (kcat * Conc_Enzyme * Conc_Substrate) / (KM + Conc_Substrate)

# def MichaelisMentenEqn_CompetitiveInhibition(Conc_Enzyme, Conc_Substrate, Conc_Inhibitor, kcat, KM, Ki):
#     return (kcat * Conc_Enzyme * Conc_Substrate) / (KM * (1 + (Conc_Inhibitor / Ki) + Conc_Substrate))

def MichaelisMentenEqn_Inhibition_Allosteric(Conc_Enzyme, Conc_Substrate, Conc_Inhibitor, kcat, KM, Ki, n):
    return (kcat * Conc_Enzyme / (1 + (Conc_Inhibitor / Ki) ** n)) * (Conc_Substrate / (KM + Conc_Substrate))

def MichaelisMentenEqn_Activation_Allosteric(Conc_Enzyme, Conc_Substrate, Conc_Inhibitor, kcat, KM, Ki, n):
    return (kcat * Conc_Enzyme * (1 + (Conc_Inhibitor / Ki) ** n)) * (Conc_Substrate / (KM + Conc_Substrate))

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
