# BSD 3-Clause License
# © 2023, The University of Texas Southwestern Medical Center. All rights reserved.
# Donghoon M. Lee, Chanhee Park, and Daehwan Kim

class FSimulator():
    def __init__(self):
        self.Iter = 0
        self.Debug_Info = 100

    def Initialize(self):
        None

    def Simulate(self, TotalTime = 24 * 60.0 * 60.0, DeltaTime = 1.0):
        print()
        print("-- Initial Conditions --")
        self.Info()
        print()

        while self.Iter < TotalTime / DeltaTime:

            # if self.Iter == 105:
            #     print("Debugging: Iter", self.Iter)

            self.SimulateDelta(DeltaTime)
            self.Iter += 1
            if self.Iter % self.Debug_Info == 0:
                print()
                print("-- Iteration {} ({:.3f}s) --".format(self.Iter, self.Iter * DeltaTime))
                self.Info()
                print()

        print("\n")
        print("-- Summary --")
        self.Summary()

    def SimulateDelta(self, DeltaTime = 1.0):
        None

    def AddToDataset(self):
        None

    def GetDataset(self):
        return None

    def Info(self):
        None

    def Summary(self):
        None

