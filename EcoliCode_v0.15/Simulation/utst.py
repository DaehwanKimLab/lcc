import unittest
import numpy as np

# Unit test for ReactionEquations class
class TestReactionEquations(unittest.TestCase):
    def setUp(self) -> None:
        super().setUp()
        self.time, self.arr2DCounts, self.arr2DMetaboliteRate, self.arr2DDeltaCounts, self.arr1DStepAdjustment, self.metaboliteLegend = run()
        self.dCounts = abs(self.arr2DDeltaCounts[0]) - (abs(self.arr2DDeltaCounts[1]) + abs(self.arr2DDeltaCounts[2]))

    def testMassConservationLower(self):
        """Ensure within a single step that no value loses more than 2 zeptomole (i.e. * 10 ** -21). 1 zmol == -602 molecules"""
        self.assertEqual( sum(counttoZeptomoles(self.dCounts) < -2), 0, msg = "Warning!! Mass differential (Loss) is greater than 2 zeptomole!")

    def testMassConservationUpper(self):
        """Ensure within a single step that no value gains more than 2 zeptomole (i.e. * 10 ** -21). 1 zmol == -602 molecules"""
        self.assertEqual( sum(counttoZeptomoles(self.dCounts) > 2), 0, msg = "Warning!! Mass differential (Gain) is greater than 2 zeptomole!")


if __name__ == '__main__':
    unittest.main()

