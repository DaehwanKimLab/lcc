import numpy as np

NA = 6.0221409e23
pi = np.pi
np.random.seed(1)


def ConcToCount(Conc_Molecule, Volume):
    return (Conc_Molecule * NA) * Volume


def CountToConc(Count_Molecule, Volume):
    return (Count_Molecule / NA) / Volume

""" Old Implementation: Attempting to phase out -- 20220603 ATM"""
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


### Regulatory/Allostery ###########################################################################################
class RegulatoryEquations:
    """ Regulatory equations
    Equations taken from: http://be150.caltech.edu/2020/content/lessons/05_ffls.html#Dynamical-equations-for-FFLs


    Methods
    -------
    andLogic2Activators
        Calculates the regulatory effect of two activators with AND logic
    orLogic2Activators
        Calculates the regulatory effect of two activators with OR logic
    andLogic2Repressors
        Calculates the regulatory effect of two repressors with AND logic
    orLogic2Repressors
        Calculates the regulatory effect of two repressors with OR logic
    andLogic1Activator1Repressor
        Calculates the regulatory effect of one activator and one repressor with AND logic
    orLogic1Activator1Repressor
        Calculates the regulatory effect of one activator and one repressor with OR logic


    Notes
    -----
    Kd = 1 is provided in each equation as an input although it is not utilized in the calcluation.
    This is until I figure out whether or not it should replace the 1 in the denom.

    """
    def __init__(self) -> None:
        pass

    def andLogic2Activators(self, conc_1, conc_2, Kd1, Kd2, n1 = 1, n2 = 1, Kd = 1):
        return ((conc_1 / Kd1) ** n1 * (conc_2 / Kd2) ** n2) / ( 1 +(conc_1 / Kd1) ** n1 * (conc_2 / Kd2) ** n2 ) 

    def orLogic2Activators(self, conc_1, conc_2, Kd1, Kd2, n1 = 1, n2 = 1, Kd = 1):
        return ((conc_1 / Kd1) ** n1 + (conc_2 / Kd2) ** n2) / ( 1 +(conc_1 / Kd1) ** n1 + (conc_2 / Kd2) ** n2 ) 

    def andLogic2Repressors(self, conc_1, conc_2, Kd1, Kd2, n1 = 1, n2 = 1, Kd = 1):
        return  1 / ( 1 +(1 + conc_1 / Kd1) ** n1 + (1 + conc_2 / Kd2) ** n2 ) 

    def orLogic2Repressors(self, conc_1, conc_2, Kd1, Kd2, n1 = 1, n2 = 1, Kd = 1):
        return  (1 + (conc_1 / Kd1) ** n1 + (conc_2 / Kd2) ** n2) / ( 1 +(1 + conc_1 / Kd1) ** n1 + (1 + conc_2 / Kd2) ** n2 ) 

    def andLogic1Activator1Repressor(self, conc_act, conc_rep, k_act, k_rep, n_act = 1, n_rep = 1, Kd = 1):
        return ((conc_act / k_act) ** n_act) / ( 1 + (conc_act / k_act) ** n_act + (conc_rep / k_rep) ** n_rep ) 
        
    def orLogic1Activator1Repressor(self, conc_act, conc_rep, k_act, k_rep, n_act = 1, n_rep = 1, Kd = 1):
        return (1 + (conc_act / k_act) ** n_act) / ( 1 + (conc_act / k_act) ** n_act + (conc_rep / k_rep) ** n_rep ) 



### if N_MOLECULEALLOWED = N ###########################################################################################
class ReactionEquations:
    """ ReactionEquations

        Reaction equations describing the mechanics of biochemical reactions.

        Attributes:
        ----------
            arr1DConc : `np.ndarray` 
                Input array (1,n) containing the concentrations of all reactants. Note, concentrations of 
                non-consumed reactants such as enzymes should be included.
            
            Rate,
                ????

            dReactionRateCoeff : `float`,              
                Mass action reaction rate k(T). Describes the direction and rate of a reaction at a given temperature.
                Units depend on order of reaction but are of the form: mol^1−(m+n)·L(m+n)^−1·s^−1. Where m, n are:  
                `rate = k(T)*[A]^m*[B]^n`.  Note m,n are not necessarily equal to the stoichiometric coefficients
                of the reaction.

            arr1DdHillCoeff : `float`, Default = 1.0
                Quantifies the degree of interaction (cooperativity) between ligand binding sites. 
                Noncooperative binding (dHillCoeff = 1.0) indicates all binding events are independent.
                Positive cooperativity (dHillCoeff > 1.0) affinity increases as ligands bind.
                Negative cooperativity (dHillCoeff < 1.0) affinity decreases as ligands bind.
                Note, The hill coefficient should always be positive (i.e. > 0.0).  

            arr4DAllostery : `np.ndarray`
                arr4DAllostery[0] : `int` , Default = np.zeros() 
                    The `substrateTargetIndex`. Holds index in the substrate array which holds the target substrate of the current regulatory element.
                arr4DAllostery[1] : `float` , Default = np.zeros()
                    The `regulatoryElementConcentration`.  Holds the concentration of the current regulatory element, if none present this should be 0.0
                arr4DAllostery[2] : `float` , Default = np.ones()
                    The `regulatoryElementDissociationConstant`.  The rate of dissociation of the regulatory element. Kd = k-unbinding/kbinding (I think)
                    Note, this element must be non-zero, o/w allostery will divide by zero (and throw error)
                arr4DAllostery[3] : `int` , Default = np.ones()
                    The `regulatoryType`.  The direction of the regulation (activation/inhibition).  Activation = 1, inhibition = -1. This value should always be either 1 or -1.
            
            arr1DdKMichaelis : `np.ndarray(float)`, Default = np.zeros              
                The concentration at which the reaction proceeds (i.e. Vmax / 2 ).  Default setting (i.e. KM = 0) suggests that the reaction velocity is independent of substrate concentration.

            arr1DdRelativeConc : `np.ndarray(float)`
                The relative concentration of each substrate as determined by: [relative_X] = [total_X] * [allosteric_X].
                Fundamentally this accounts for situations where the rate limiting concentration is not the absolute minimum concentration due to allosteric effects.               

            rateLimitingConcentration : `float`
                The minimum of arr1DRelativeConc

        Methods
        -------
        CheckRateAndConc
            ????

        MassActionReactionRate
            Determines the rate of the current reaction using mass action kinetics.
            Under mass action, a reaction rate is the product concentration of all reactants (A, B, ..., N) and the
            corresponding rate coefficient at temperature T.
            i.e. rate = dReactionRateCoefficent * [A][B][...][N]
        
        Reaction 
            The Cumulative reaction described by both mass action reaction rate and allosteric effects.
            Note the default values of allosteric elements equate to 1 such that the reaction can be
            summarized as an unregulated mass action scheme.

        Allostery 
            Calculates the allosteric effects on the current reaction. 
            Currently, if the substrate of a reaction has an allosteric effect on one of its regulators, this is ignored.

        Saturation
            Michaelis-Menten like hyperbolic saturation.  [S]^n / KM^n + [S]^n. Note cooperativity is also taken into consideration (i.e. n) by means
            of the Hill coefficient. Note: The assumptions for the Michaelis Meneten model are not being considered (i.e. assumed to hold).
        
        RelativeConcentration
            Determines the relative concentration of reaction substrates.  Note the relative concentration can be larger than the absolute concentration. 
            This is intended in order to allosterically modulate KM (michaelis constant) without actually changing its value. These relative concentrations
            are not propogated out, but instead are used internally to calculate the change in product for the current reaction.

        ReactionMichaelisMenten
            calculates the change in product for a reaction base on the rate limiting concentration, the reaction rate (k-value), and the saturation.
            Note: this is truely a Pseudo-MichaelisMenten and instead assuming the enzyme is the rate limiting concentration, the minimum
            relative concentration is used. Theoretically, however, the rate limiting concentration should almost always be the enzyme concentration.

        Notes
        -----
        Current implementations do not consider saturation or diffusion limited reaction rates...
        Fixed (although with limitations) >> Currently at most one allosteric modulation can be considered..

        Current version uses list comprehension.  This should be vectorized.
        Currently there is no check if the deltaProduct is greater than the stoichiometry of substrates...
        
        TODO
        ----
        Implement AND/OR operations i.e.
            AND -- allosteric effect only occurs if both [A] AND [B] are present
            OR -- allosteric effect occurs if either [A] OR [B] are present

        Vectorize list comprehensions
        Initialize default values (np.zeros/np.ones where needed)
        (?) Check deltaProduct greater than stoichiometry allows (?)

    """

    def __init__(
        self,
        arr1DConc: np.ndarray,  # Conc. of all reaction components (now incl. enz)
        Rate,  # How is rate different from k?
        dReactionRateCoeff: np.ndarray,  # float,  # k
        arr1DdHillCoeff: np.ndarray,  # float = 1.0,  # Default hill coeff
        arr4DAllostery: np.ndarray,
        arr1DdKMichaelis: np.ndarray,
    ):
        #
        self.Rate = Rate
        self.dReactionRateCoeff = dReactionRateCoeff

        self.arr1DConc = arr1DConc
        self.arr1DdKMichaelis = arr1DdKMichaelis
        self.arr1DdHillCoeff = arr1DdHillCoeff  

        # Need (n, 4) array [0] holds the index of the substrate in arr1Conc, [1] conc allosteric regulator
        # [2] dissociation constant of allosteric regulator [3] type of allosteric regulation, 
        # [4] holds index of in arr5DAllostery of substrate should it exist -- regulator of regulators
        # Removed the regulator of regulators... It gets into a circular logic problem.
        self.arr4DAllostery = arr4DAllostery

        # Currently relative concentraitons of allosteric components override absolute concentrations... 
        # This may be a problem if a substrate is also a regulator...

        # Get Relative concentrations
        # Relative conc. regulators
        #self.arr5DAllostery[1] = self.RelativeConcentration(self.arr5DAllostery[1], 4)
        
        # Relatice conc. Reaction components
        self.arr1DdRelativeConc = self.RelativeConcentration(self.arr1DConc, 0)
        # Get Rate Limiting concentration
        self.rateLimitingConcentration = np.min(self.arr1DdRelativeConc)

    # I'm not sure why we are taking the minimum of concentration and rate...
    def CheckRateAndConc(self):
        """ ???? """
        return np.min(
            (self.Rate, self.arr1DConc), axis=0
        )  # may need to get min of array 1st...

    def MassActionReactionRate(self):
        """Mass action reaction rate:
        Rate of reaction is the product of all reactants and the reaction coefficient.
        """
        return self.dReactionRateCoeff * np.prod(self.arr1DConc)       

    def Reaction(self):
        """Total Reaction: 
        The MassActionReactionRate modulated by allosteric elements.
        """
        # Need to update this b/c I changed Allostery
        return self.dReactionRateCoeff * np.prod(self.arr1DdRelativeConc)

    # [Relative] and Rate limiting step
    def RelativeConcentration(self, arr1DConc, allosteryMatchIndex  = 0):
        """Combine allosteric effects with concentration to get relative concentration"""
        return np.array(
            [
                arr1DConc[i]
                * self.Allostery(
                    arr4DAllostery = np.array([np.where(self.arr4DAllostery[allosteryMatchIndex] == i)])  # Subset of arr4D with allosteric elements for current substrate
                )
                for i in range(arr1DConc.shape[0]) 
                if self.arr4DAllostery[allosteryMatchIndex].all() == i
            ]
        )

    def Allostery(self, arr4DAllostery):
        """ Calculate the impact of allosteric elements on reaction
        Pass in subset of Arr4D for one regulated component.

        """
        # note what was formerly "n" was removed as I don't think this cooperative effect applies here... 
        ## Actually the above is false and I should add a hill coeff in.  Its totally possible that a protein can bind multiple
        ## of the same inhibitor/activator and they may have cooperative effects.

        # Loop through, get individual allosteric elements,
        # overall allostery is the product of those elements
        # Note: I'm not sure if this works for regulator of regulators is going to work... 
        # We need to allostery(regulators), get allosteric interactions between regulators...
        return np.prod(
            np.array(
                [
                    (1 + (arr4DAllostery[i, 1] / arr4DAllostery[i, 2]))
                    ** arr4DAllostery[i, 3]
                    for i in range(arr4DAllostery[0].shape[0]-1) ########################### Not sure if this is right but it un-broke it...
                ]
            )
        )

    def Saturation(self):
        """ Calculate Saturation (i.e. relative availability) + cooperativity (via hill coeff)"""

        return np.prod(
            np.array(
                [
                    self.arr1DConc[i] ** self.arr1DdHillCoeff[i]
                    / (
                        self.arr1DdKMichaelis[i] ** self.arr1DdHillCoeff[i]
                        + self.arr1DConc[i] ** self.arr1DdHillCoeff[i]
                    )
                    for i in range(self.arr1DConc.shape[0])
                ]
            )
        )

    def MichaelisMentenReaction(self):
        """Total reaction using MichaelisMenten-like kinetics"""
        return (
            self.dReactionRateCoeff * self.rateLimitingConcentration * self.Saturation()
        )

###### End of code update

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
    Weight_Cumsum_Min = np.where(
        Weight_Cumsum_Min == Weight_Cumsum_Max, Weight_Cumsum_Min - 1, Weight_Cumsum_Min
    )
    RanNums = np.asmatrix(
        np.random.randint(Weight_Cumsum_Min, high=Weight_Cumsum_Max, size=Quantity)
    ).transpose()
    # Generate a matrix of the random numbers for comparison to indices
    RanNums_Matrix = np.reshape(
        np.repeat(RanNums, Indices.shape[1]), [-1, Indices.shape[1]]
    )
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
    Bool_LenCumsumLessThanOrEqualToCountOfIndices = np.less_equal(
        LenCumsum, Count_Indices
    )
    Bool_LenCumsum = np.logical_and(
        Bool_LenCumsumGreaterThanZero, Bool_LenCumsumLessThanOrEqualToCountOfIndices
    )
    Bin_Len_Selected = np.logical_and(Bool_LenAvailable, Bool_LenCumsum).astype(int)
    return Len + Bin_Len_Selected


# Spatial simulation functions
def InitializeDistribution(
    Width,
    Height,
    X_Ori,
    Y_Ori,
    MaxAmount,
    BasalAmount=0,
    shape="",
    size=0,
    pattern="diffuse",
):
    """
    available shape options: circle, square
    temporary size options are interpreted as its diameter or length of each side, respectively.
    """

    if shape != "":
        assert (
            size > 0 and pattern != ""
        ), "Distribution Error: When shape is defined, size and pattern must be provided"

    Dist_Init = np.full((Width, Height), BasalAmount, dtype="float32")

    # Early termination when there is nothing to distribute
    if MaxAmount == 0:
        return Dist_Init

    # Determine Indices to fill based on the shape
    X_Min = int(X_Ori - size / 2)
    X_Max = int(X_Ori + size / 2)
    Y_Min = int(Y_Ori - size / 2)
    Y_Max = int(Y_Ori + size / 2)

    X_Array = []
    Y_Array = []

    if shape == "circle":
        for X in range(X_Min, X_Max):
            for Y in range(Y_Min, Y_Max):
                if np.sqrt((X_Ori - X) ** 2 + (Y_Ori - Y) ** 2) < size / 2:
                    X_Array.append(X)
                    Y_Array.append(Y)
    elif shape == "square":
        for X in range(X_Min, X_Max):
            for Y in range(Y_Min, Y_Max):
                X_Array.append(X)
                Y_Array.append(Y)
    else:
        for X in range(Width):
            for Y in range(Height):
                X_Array.append(X)
                Y_Array.append(Y)

    Coords_ToFill = (np.array(X_Array), np.array(Y_Array))

    # Fill the shape with distribution
    if pattern == "solid":
        Dist_Init[Coords_ToFill] = MaxAmount
    elif shape != "" and pattern == "diffuse":
        Dist_Init[Coords_ToFill] = [
            InitialPattern(X, Y, size, size, X_Ori, Y_Ori, MaxAmount)
            for X, Y in zip(Coords_ToFill[0], Coords_ToFill[1])
        ]
    else:
        Dim_Shorter = min(Width, Height)
        Dist_Init[Coords_ToFill] = [
            InitialPattern(X, Y, Dim_Shorter, Dim_Shorter, X_Ori, Y_Ori, MaxAmount)
            for X, Y in zip(Coords_ToFill[0], Coords_ToFill[1])
        ]
    return Dist_Init


def InitialPattern(
    X, Y, Width, Height, X_Ori, Y_Ori, Max,
):
    Dist = np.sqrt(((X - X_Ori) / Width) ** 2 + ((Y - Y_Ori) / Height) ** 2)
    return Max / max(1, Dist * 30)


def Eqn_Transporter_Unregulated(Conc_Transporter, Conc_Inside, Conc_Outside, ki, ko):
    return Conc_Transporter * (Conc_Inside * ko - Conc_Outside * ki)


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

    Distribution_Diffused = (Upward + Downward + Leftward + Rightward)[
        1:-1, 1:-1
    ] - 4 * Diffusion_Quantity

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

    Distribution_Diffused = (
        Upward
        + Downward
        + Leftward
        + Rightward
        + UpLeftward
        + UpRightward
        + DownLeftward
        + DownRightward
    )[1:-1, 1:-1] - 8 * Diffusion_Quantity

    return Distribution_Diffused


def DiffuseDistribution_4Cell(Distribution, D=0.1, dTime=1):  # D must be less than 1/6
    Distribution_Updated = (
        Distribution + Eqn_Diffusion_Spatial_4Cell(Distribution, D) * dTime
    )
    # return Distribution
    return Distribution_Updated


def DiffuseDistribution_8Cell(Distribution, D=0.1, dTime=1):  # D less than 1/10
    Distribution_Updated = (
        Distribution + Eqn_Diffusion_Spatial_8Cell(Distribution, D) * dTime
    )
    # return Distribution
    return Distribution_Updated


def Displacement_2D(Distance, Angle):
    dX = Distance * np.cos(Angle)
    dY = -Distance * np.sin(Angle)
    return dX, dY


# Temporary chemotaxis simulation functions.
def BacterialChemotaxis(ThresholdSubject, X, Y, Angle, Threshold):
    # optimize the following code into additions and subtractions
    Decision = np.any(ThresholdSubject < Threshold, axis=0)
    X_Run, Y_Run, Angle_Run = Chemotaxis_Run(X, Y, Angle)
    X_Tumble, Y_Tumble, Angle_Tumble = Chemotaxis_Tumble(X, Y, Angle)
    X_New = np.where(Decision, X_Run, X_Tumble)
    Y_New = np.where(Decision, Y_Run, Y_Tumble)
    Angle_New = np.where(Decision, Angle_Run, Angle_Tumble)
    return X_New, Y_New, Angle_New


def Chemotaxis_Run(X, Y, Angle, Distance=15):
    dX, dY = Displacement_2D(Distance, Angle)
    return X + dX, Y + dY, Angle


def Chemotaxis_Tumble(X, Y, Angle, Distance=5):
    dX, dY = Displacement_2D(Distance, Angle)
    NewAngle = GetRandomAngle(Angle.shape)
    return X + dX, Y + dY, NewAngle


def GetRandomAngle(AngleShape):
    return np.random.random_sample(AngleShape) * 2 * pi


def CorrectOutOfBounds(X, Y, Angle, Width, Height):
    X_LessThanZero = X < np.array([1])
    X_MoreThanWidth = X >= np.array([Width - 1])
    Y_LessThanZero = Y < np.array([1])
    Y_MoreThanHeight = Y >= np.array([Height - 1])
    Angle_ToChange = np.any(
        np.stack([X_LessThanZero, X_MoreThanWidth, Y_LessThanZero, Y_MoreThanHeight]),
        axis=0,
    )

    X_Corrected = np.where(X_LessThanZero, 1, X)
    X_Corrected = np.where(X_MoreThanWidth, Width - 2, X_Corrected)
    Y_Corrected = np.where(Y_LessThanZero, 1, Y)
    Y_Corrected = np.where(Y_MoreThanHeight, Height - 2, Y_Corrected)
    Angle_Corrected = np.where(Angle_ToChange, GetRandomAngle(Angle.shape), Angle)
    return X_Corrected, Y_Corrected, Angle_Corrected


def Normalize_Linear(Data):
    return Data / np.max(Data)


def Normalize_P1Log(Data):
    Data_Modified = np.log(Data + 1)
    return Data_Modified / np.max(Data_Modified)


# Debugging tools
def SciFloat(Float, InPrecision=4, InExp_digits=2):
    return np.format_float_scientific(
        Float, precision=InPrecision, unique=False, exp_digits=InExp_digits
    )
