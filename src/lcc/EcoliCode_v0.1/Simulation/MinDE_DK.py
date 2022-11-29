#!/usr/bin/env python3
import sys
from datetime import datetime
import numpy as np
import math

np.random.seed(1)

pi = np.pi

# Molecules
NumMinD = 5000
NumCytoMinD = NumMinD
NumMinE = 4000
NumCytoMinE = NumMinE

class FMinCluster:
    def __init__(self, Index, Angle, Reference, SizeD, SizeE):
        self.Index = Index
        self.Angle = Angle
        self.Reference = Reference
        self.SizeD = SizeD
        self.SizeE = SizeE

    def CustomFunc(self, Rad):
        return math.cosh(Rad) - 0.9
        
    def KonD(self):
        Max = self.CustomFunc(0.5 * pi)
        NewAngle = self.Angle if self.Angle < pi else self.Angle - pi
        Num = self.CustomFunc(NewAngle - 0.5 * pi)
        Num = Num / Max * 1.1 * math.log(max(2, self.SizeD), 2)
        return Num

    def KoffD(self):
        return self.KonE()

    def KonE(self):
        Max = self.CustomFunc(0.5 * pi)
        if self.Angle < 0.5 * pi:
            NewAngle = self.Angle
        elif self.Angle < pi:
            NewAngle = pi - self.Angle
        elif self.Angle < 1.5 * pi:
            NewAngle = self.Angle - pi
        else:
            NewAngle = 2 * pi - self.Angle

        Num = self.CustomFunc(NewAngle)
        Num = Num / Max * 1.6 * math.log(max(1, self.SizeD), 2) * math.log(max(2, self.SizeE), 2)
        return Num

    def KoffE(self):
        return self.KonE() / 2

    def PrintInfo(self):
        print("rad: {:.2f}, MinD: {:5d}, Kon: {:.2f}".format(self.Angle, self.SizeD, self.Kon()))


def SimulateMinCDE():
    NumAngles = 100
    DeltaAngle = 2 * pi / NumAngles
    MinClusters = [FMinCluster(i, DeltaAngle * i, -1, 0, 0) for i in range(NumAngles)]

    def AddMinDs():
        global NumCytoMinD
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Kon = MinCluster.KonD()
            ToAdd = min(int(Kon), NumCytoMinD)
            if ToAdd <= 0:
                continue

            MinCluster.SizeD += ToAdd
            NumCytoMinD -= ToAdd

    def RemoveMinDs():
        global NumCytoMinD
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Koff = MinCluster.KoffD()
            ToRemove = min(int(Koff), MinCluster.SizeD)
            if ToRemove <= 0:
                continue
            
            MinCluster.SizeD -= ToRemove
            NumCytoMinD += ToRemove

    def AddMinEs():
        global NumCytoMinE
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Kon = MinCluster.KonE()
            ToAdd = min(int(Kon), NumCytoMinE)
            if ToAdd <= 0:
                continue

            MinCluster.SizeE += ToAdd
            NumCytoMinE -= ToAdd

    def RemoveMinEs():
        global NumCytoMinE
        for MinCluster in MinClusters:
            if MinCluster.Reference >= 0:
                continue

            Koff = MinCluster.KoffE()
            ToRemove = min(int(Koff), MinCluster.SizeE)
            if ToRemove <= 0:
                continue

            MinCluster.SizeE -= ToRemove
            NumCytoMinE += ToRemove

    def MergeMinClusters():
        MinimumSize = 50
        ci = 0
        while ci < len(MinClusters):
            C1 = MinClusters[ci]
            if C1.Reference >= 0:
                ci += 1
                continue

            ni = ci + 1
            while True:
                # niw as ni_wrapped
                niw = ni % len(MinClusters)
                C2 = MinClusters[niw]
                if C2.Reference >= 0:
                    ni += 1
                    continue
                else:
                    break

            if C1.SizeD > 0 and C2.SizeD > 0 and False:
                print("C1 {}".format(ci))
                C1.PrintInfo()
                print("C2 {}".format(ni))
                C2.PrintInfo()

            if C1.SizeD >= MinimumSize and C2.SizeD >= MinimumSize:
                if C1.SizeD > C2.SizeD:
                    C2.Reference = ci
                    C1.SizeD += C2.SizeD
                    C1.SizeE += C2.SizeE
                    ci = ci
                else:
                    C1.Reference = niw
                    C2.SizeD += C1.SizeD
                    C2.SizeE += C1.SizeE
                    ci = niw
            else:
                ci = niw

            if ni >= len(MinClusters):
                break

    def FragmentMinClusters():
        MinimumSize = 50
        ci = 0
        while ci < len(MinClusters):
            C1 = MinClusters[ci]
            if C1.Reference >= 0:
                ci += 1
                continue

            if C1.SizeD >= MinimumSize * 2 or C1.SizeD <= MinimumSize:
                ci += 1
                continue

            ni = ci + 1
            niw = ni % len(MinClusters)
            C2 = MinClusters[niw]
            if C2.Reference < 0:
                ci = ni
                continue

            C1SizeD = int(C1.SizeD / 2)
            C1SizeE = int(C1.SizeE / 2)

            C2.SizeD = C1.SizeD - C1SizeD
            C2.SizeE = C1.SizeE - C1SizeE
            C1.SizeD = C1SizeD
            C1.SizeE = C1SizeE

            C2.Reference = -1
            
            ci = ni + 1


    Iter = 0
    while Iter < 2000:
        AddMinDs()

        MergeMinClusters()        

        RemoveMinDs()

        FragmentMinClusters()

        AddMinEs()

        # RemoveMinEs()

        if Iter % 10 == 0:
            NumMemMinD = 0
            NumMemMinE = 0
            print("Iter {}".format(Iter))
            for i, MinCluster in enumerate(MinClusters):
                Angle = MinCluster.Angle
                SizeD = MinCluster.SizeD
                SizeE = MinCluster.SizeE
                Reference = MinCluster.Reference
                if SizeD <= 0 and Reference < 0:
                    continue

                if Reference < 0:
                    KonD = MinCluster.KonD()
                    KonE = MinCluster.KonE()
                    print("#{:<5d} rad: {:.2f}, MinD: {:5d}, KonD: {:.2f}, MinE: {:5d}, KonE: {:.2f}".format(i, Angle, SizeD, KonD, SizeE, KonE))
                    NumMemMinD += SizeD
                    NumMemMinE += SizeE
                else:
                    print("#{:<5d} Reference: #{:<5d}".format(i, Reference))

            print("cMinD: {:5d}, mMinD: {:5d}".format(NumCytoMinD, NumMemMinD))
            assert NumCytoMinD + NumMemMinD == NumMinD
            print("cMinE: {:5d}, mMinE: {:5d}".format(NumCytoMinE, NumMemMinE))
            assert NumCytoMinE + NumMemMinE == NumMinE
            print("\n")

        Iter += 1



if __name__ == '__main__':
    SimulateMinCDE()
