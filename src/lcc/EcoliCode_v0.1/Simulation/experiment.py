"""

"""
import numpy as np
from copy import deepcopy
from util import printBlockMessage

class Titration:
    """ Perform titration on a model.
    """
    def __init__(self, model, inputName: str, outputName: list, inputGridSearchRange = np.logspace(-9, 2, 12), maxSteps = 2000) -> None:
        self.model = model
        self.inputName = inputName
        self.outputName = outputName
        self.inputGridSearchRange = inputGridSearchRange
        self.maxSteps = maxSteps
        self.titrationResults = {}

    def _findMaxTimeKey(self):
        assert self.titrationResults == True, "Titration results appears to be empty! Did the `run` method get called?" 
        # Gets a tuple from each exp with (expName, maxTimepoint)
        titrationTimes = [(key, self.titrationResults[key]['time'][-1]) for key in self.titrationResults.keys()]
        # return the key from the titration with the largest 'time' value
        return max(titrationTimes, key=lambda tup: tup[1])[0]

    def run(self):
        assert self.outputName is not None, "The experiment has no output parameter! Please select an output and run again."
        # This is sloppy, but may work for now.
        # Basically I need to modify the inital condition of the concentration
        printBlockMessage("Beginning Titration experiment...", ts = True)
        for conc in self.inputGridSearchRange:
            print(f"Beginning simulation for {self.inputName} with concentration {conc}")
            # Set up, copy base model
            currModel = deepcopy(self.model)
            # Input the experimental params
            currModel.time_totalSteps = self.maxSteps
            currModel.perturb = False
            currModel.silent = True
            currModel.dictCountArrays[self.inputName][0] = conc
            currModel.molecules[self.inputName] = conc
            # Run sim with updated params
            outT, outC, _, _ = currModel.run()

            # Save sim results
            self.titrationResults[f'{self.inputName}_{conc}'] = {'time':outT, 'conc':{key:outC[key] for key in outC.keys() if key in self.outputName}}
            ## REPEAT
        
        printBlockMessage("Titration experiment complete!", ts = True)
        # Return the desired output results.
        return self.titrationResults