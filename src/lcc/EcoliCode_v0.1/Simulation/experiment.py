"""

"""
import numpy as np
from copy import deepcopy

class Titration:
    """ Perform titration on a model.
    """
    def __init__(self, model, inputName: str, outputName: list(str), inputGridSearchRange = np.logspace(-9, 2, 12), maxSteps = 1000) -> None:
        self.model = model
        self.inputName = inputName
        self.outputName = outputName
        self.inputGridSearchRange = inputGridSearchRange
        self.maxSteps = maxSteps
        self.titrationResults = {}

    def run(self):
        # This is sloppy, but may work for now.
        # Basically I need to modify the inital condition of the concentration
        for conc in self.inputGridSearchRange:
            # Set up, copy base model
            currModel = deepcopy(self.model)
            # Input the experimental params
            currModel.time_totalSteps = self.maxSteps
            currModel.perturb = False
            currModel.dictCountArrays[self.inputName][0] = conc
            currModel.molecules[self.inputName] = conc
            # Run sim with updated params
            outT, outC, _, _ = currModel.run()

            # Save sim results
            self.titrationResults[f'{self.inputName}_{conc}'] = {'time':outT, 'conc':outC}
            ## REPEAT
        
        # Return the desired output results.
        return self.titrationResults