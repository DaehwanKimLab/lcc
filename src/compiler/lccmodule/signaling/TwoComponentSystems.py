
def SetUpData_TwoComponentSystems(self, Dataset): # TCS for short
    # TwoComponentSystems = Dataset['twoComponentSystems.tsv']
    TwoComponentSystems = Dataset['TwoComponentSystemsTemporary_DL.tsv'] # temporary data table
    # Get the ID and Index list of molecule involved in TCS
    self.N_TCSMolecules = len(TwoComponentSystems)
    self.TCSMolIndexes = np.zeros(self.N_TCSMolecules)
    for i, Value in enumerate(TwoComponentSystems):
        Name, Count = Value
        self.TCSMolNames.append(Name)
        # self.TCSMolIndexes[i] = self.ProteinID2Index[Name]
    # Get the index list of TCS molecules
    return
