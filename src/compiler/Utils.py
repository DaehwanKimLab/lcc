"""
Reaction Matrix Building functions
"""



def ParseBuildingBlocks_(MolGroup, Stoichiometry, Comp):
    MolGroupParsed = []
    StoichiometryParsed = []
    if MolGroup == 'dNTP':
        MolGroupParsed = Comp.BuildingBlock.Name_dNTPs
        StoichiometryParsed = FlatList(
            Comp.Chromosome.Freq_NTsInChromosomesInGenome)  # This is a temporary solution for a single chromosome
    elif MolGroup == 'NTP':
        MolGroupParsed = Comp.BuildingBlock.Name_NTPs
        StoichiometryParsed = Comp.RNA.Freq_NTsInRNAs
    elif MolGroup == 'AA':
        MolGroupParsed = Comp.BuildingBlock.Name_AAs
        StoichiometryParsed = Comp.Protein.Freq_AAsInProteins
    StoichiometryParsed = [i * Stoichiometry for i in StoichiometryParsed]
    return MolGroupParsed, StoichiometryParsed

def RefineBuildingBlocks(IDs, Stoichiometries, Comp, ReactionOrder=None):
    assert len(IDs) == len(Stoichiometries), "#'s of Molecule IDs and Stoichiometries do not match"
    BuildingBlocks = ['dNTP', 'NTP', 'AA']
    # MolTypes = ['Chromosome', 'Gene', 'Promoter', 'RNA', 'Protein', 'Complex', 'Metabolite']
    IDs_Refined = []
    Stoichiometries_Refined = []
    ReactionOrder_Refined = []
    for ID, Stoichiometry in zip(IDs, Stoichiometries):
        if ID in Comp.Master.ID_Master:
            IDs_Refined.append(ID)
            Stoichiometries_Refined.append(Stoichiometry)
        elif ID in BuildingBlocks:
            ID_BuildingBlocks, Stoich_BuildingBlocks = ParseBuildingBlocks_(ID, Stoichiometry, Comp)
            for ID_BuildingBlock in ID_BuildingBlocks:
                assert ID_BuildingBlock in Comp.Master.ID_Master, '%s not found in Master ID list' % ID_BuildingBlock
            IDs_Refined.append(ID_BuildingBlocks)
            Stoichiometries_Refined.append(Stoich_BuildingBlocks)
        else:
            print('Molecule ID not defined in the organism: %s' % ID)
    IDs_Refined = FlatList(IDs_Refined)
    Stoichiometries_Refined = FlatList(Stoichiometries_Refined)
    return IDs_Refined, Stoichiometries_Refined


def PrepareReactionForMatrix(Stoichiometry, Comp):
    Idx, Stoich = None, None
    for MolIDs, Stoichiometries in Stoi:
        Idx, Stoich = RefineBuildingBlocks(MolIDs, Stoichiometries, Comp)
    return Idx, Stoich

def ParseStoichiometry(Stoichiometry, Comp):
    [MolIDs, Coeffs] = Stoichiometry
    Idx, Stoich = None, None
    Stoichiometry_Parsed = RefineBuildingBlocks(MolIDs, Coeffs, Comp)
    return Stoichiometry_Parsed

def ParseStoichiometries(Stoichiometries, Comp):
    Stoichiometries_Parsed = []
    for Stoichiometry in Stoichiometries:
        Stoichiometry_Parsed = ParseStoichiometry(Stoichiometry)
        Stoichiometries_Parsed.append(Stoichiometry_Parsed)
    return Stoichiometries_Parsed

def ParseRXNRate(Rate):
    # Parse Rxn rate (to be expanded)
    RXNRateParsed = Rate
    return RXNRateParsed

def FlatList(List):
    ListFlattened = []
    for i in List:
        if isinstance(i, str):
            ListFlattened.append(i)
        elif isinstance(i, int):
            ListFlattened.append(i)
        elif isinstance(i, float):
            ListFlattened.append(i)
        else:
            for j in i:
                ListFlattened.append(j)
    return ListFlattened


def GetMolIdx(Molecules, MolIdxRef):
    MolIdxList = []

    # For ID2Idx cases
    if isinstance(MolIdxRef, dict):
        for Molecule in Molecules:
            MolIdx = MolIdxRef[Molecule]
            MolIdxList.append(MolIdx)
        return MolIdxList
    # For Type2Idx cases
    elif isinstance(MolIdxRef, list):
        for MolIdx, Type in enumerate(MolIdxRef):
            if Type == Molecules:
                MolIdxList.append(MolIdx)
        return MolIdxList
    else:
        print("Inappropriate reference type used in GetMolIdx function parameter: %s" % MolIdxRef)

def ReverseNumberSign(ListOfNumbers):
    NewListOfNumbers = [-Number for Number in ListOfNumbers]
    return NewListOfNumbers

