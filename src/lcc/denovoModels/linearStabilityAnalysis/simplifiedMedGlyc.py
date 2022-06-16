import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Avogadro
from scipy import sqrt
from scipy.optimize import fixed_point
import sympy as sm

from fun import Plotly_Dynamics, generateRGBList

#plotting
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from random import randint
import os
import re
from plotly.subplots import make_subplots

# Micromolar conversion 
def uMtoCount(inVal:float):
    return inVal*10**-6*Avogadro
def counttouM(inVal:float):
    return inVal / Avogadro / 10**-6

# constants
SIM_DURATION_SECONDS = 10
STEPS_PER_SECOND = 100

TOTAL_STEPS = SIM_DURATION_SECONDS * STEPS_PER_SECOND

# kcats
KCAT_HEXOKINASE = 4e3
KCAT_PFK = 1.7e4
KCAT_PFK_ACTIVE = 1.7e4
KCAT_PFK_INACTIVE = 1.7e4
KCAT_PK = 170
##
KCAT_SUMMARY_MAKE_PEPASE = 1e4 # Arbitrary
# ks
K_SUMMARY_MAKE_PEPASE = 1e4

# KMs
KM_HEXOKINASE_GLUCOSE = uMtoCount(5)
KM_PFK_ACTIVE = uMtoCount(0.21)
KM_PFK_INACTIVE = uMtoCount(3.1e7)
KM_PK_PEP = uMtoCount(700)

#####
KM_SUMMARY_MAKE_PEPASE = uMtoCount(30) # just guessing.

# Concentrations
## Metabolites
CONC_GLUCOSE_T0 = uMtoCount(500)
CONC_ADP_T0 = uMtoCount(1700)
CONC_ATP_T0 = uMtoCount(5000)
## Intermediates
CONC_ATP = CONC_ATP_T0
CONC_ADP = CONC_ADP_T0
CONC_F6P = uMtoCount(15)
CONC_F16P = 0
CONC_PEP = 0
CONC_PYRUVATE = 0
## Enzymes
CONC_PFK = 4000 / Avogadro
CONC_PFK_ACTIVE = 0
CONC_PFK_INACTIVE = 0
CONC_HEXOKINASE = uMtoCount(0.01)
CONC_SUMMARY_MAKE_PEPASE = uMtoCount(1)
CONC_PYRUVATE_KINASE = uMtoCount(1)
CONC_PHOSPHOGLUCOSEISOMERASE = uMtoCount(1)

## Current exploratory values
KI_ATP_PYRUVATEKINASE = uMtoCount(5000)

useOneStatePFK = True

def saturation(conc, KM):
    return (conc / (KM + conc))

# Attempting to base
def hexokinase(Conc_Substrate = CONC_GLUCOSE_T0, Conc_Enzyme = CONC_HEXOKINASE, kcat = KCAT_HEXOKINASE, KM = KM_HEXOKINASE_GLUCOSE):
    # Note: glucose is not decreasing...
    return (kcat * Conc_Enzyme * saturation(Conc_Substrate, KM))

if useOneStatePFK:
    def pfk(Conc_Substrate = CONC_F6P, Conc_Enzyme = CONC_PFK, kcat = KCAT_PFK_ACTIVE, KM = KM_PFK_ACTIVE):
        return (kcat * Conc_Enzyme * saturation(Conc_Substrate, KM))
else:
    def pfk_active(Conc_Substrate = CONC_F6P, Conc_Enzyme = CONC_PFK_ACTIVE, kcat = KCAT_PFK_ACTIVE, KM = KM_PFK_ACTIVE):
        return (kcat * Conc_Enzyme *  saturation(Conc_Substrate, KM))
    def pfk_inactive(Conc_Substrate = CONC_F6P, Conc_Enzyme = CONC_PFK_INACTIVE, kcat = KCAT_PFK_INACTIVE, KM = KM_PFK_INACTIVE):
        return (kcat * Conc_Enzyme * saturation(Conc_Substrate, KM))

def summaryMakePEPase(Conc_Substrate = CONC_F16P,Conc_Enzyme = CONC_SUMMARY_MAKE_PEPASE, kcat = KCAT_SUMMARY_MAKE_PEPASE, KM = KM_SUMMARY_MAKE_PEPASE):
    return(kcat * Conc_Enzyme * saturation(Conc_Substrate, KM))

def pyruvateKinase(Conc_Substrate = CONC_PEP, Conc_Inhibitor = CONC_ATP, Ki = KI_ATP_PYRUVATEKINASE, Conc_Enzyme = CONC_PYRUVATE_KINASE, kcat = KCAT_PK, KM = KM_PK_PEP):
    return (kcat * Conc_Enzyme / (1 + (Conc_Inhibitor / Ki)) * saturation(Conc_Substrate, KM))

def run():
    # reimplement simulation: 
    metaboliteLegend = ['F6P', 'F16BP', 'PEP', 'pyruvate']
    enzymeLegend = ['HK', 'PFK', 'PS', 'PK']
    time = np.linspace(0, SIM_DURATION_SECONDS, TOTAL_STEPS)

    # For now only worry about F6P, F16BP, PEP, pyruvate
    arr2DCounts= np.zeros([TOTAL_STEPS, 4]) 
    arr2DEnzymePercentSubstrateSaturation = np.zeros([TOTAL_STEPS, 4])
    arr2DEnzymePercentMaxActivty = np.zeros([TOTAL_STEPS, 4])
    # 4 enzymes atm
    arr2DEnzymeStepTurnover = np.zeros([TOTAL_STEPS, 4])

    for step in range(1,TOTAL_STEPS):
        # Get turnover for each enzyme
        arr2DEnzymeStepTurnover[step][0] = hexokinase() / STEPS_PER_SECOND                               # In: Gluc  ; Out: F6P
        arr2DEnzymeStepTurnover[step][1] = pfk(arr2DCounts[step-1][0]) / STEPS_PER_SECOND                 # In: F6P   ; Out: F16BP
        arr2DEnzymeStepTurnover[step][2] = summaryMakePEPase(arr2DCounts[step-1][1]) / STEPS_PER_SECOND   # In: F16BP ; Out: PEP
        arr2DEnzymeStepTurnover[step][3] = pyruvateKinase(arr2DCounts[step-1][2]) / STEPS_PER_SECOND      # In: PEP   ; Out: Pyruvate

        # Get percent substrate Saturation:
        arr2DEnzymePercentSubstrateSaturation[step] = saturation(CONC_GLUCOSE_T0, KM_HEXOKINASE_GLUCOSE)
        arr2DEnzymePercentSubstrateSaturation[step] = saturation(arr2DCounts[step-1][0], KM_PFK_ACTIVE)
        arr2DEnzymePercentSubstrateSaturation[step] = saturation(arr2DCounts[step-1][1], KM_SUMMARY_MAKE_PEPASE)
        arr2DEnzymePercentSubstrateSaturation[step] = saturation(arr2DCounts[step-1][2], KM_PK_PEP)

        # Get Percent Max Activity
        arr2DEnzymePercentMaxActivty[step][0] = round(arr2DEnzymeStepTurnover[step][0] / (KCAT_HEXOKINASE * CONC_HEXOKINASE) * 100 * STEPS_PER_SECOND, 2)
        arr2DEnzymePercentMaxActivty[step][1] = round(arr2DEnzymeStepTurnover[step][1] / (KCAT_PFK * CONC_PFK) * 100 * STEPS_PER_SECOND, 2)
        arr2DEnzymePercentMaxActivty[step][2] = round(arr2DEnzymeStepTurnover[step][2] / (KCAT_SUMMARY_MAKE_PEPASE * CONC_SUMMARY_MAKE_PEPASE) * 100 * STEPS_PER_SECOND, 2)
        arr2DEnzymePercentMaxActivty[step][3] = round(arr2DEnzymeStepTurnover[step][3] / (KCAT_PK * CONC_PYRUVATE_KINASE) * 100 * STEPS_PER_SECOND, 2)

        # Update counts
        arr2DCounts[step][0] = max(0, arr2DCounts[step-1][0] + arr2DEnzymeStepTurnover[step][0] - arr2DEnzymeStepTurnover[step][1])
        arr2DCounts[step][1] = max(0, arr2DCounts[step-1][1] + arr2DEnzymeStepTurnover[step][1] - arr2DEnzymeStepTurnover[step][2])
        arr2DCounts[step][2] = max(0, arr2DCounts[step-1][2] + arr2DEnzymeStepTurnover[step][2] - arr2DEnzymeStepTurnover[step][3])
        arr2DCounts[step][3] = max(0, arr2DCounts[step-1][3] + arr2DEnzymeStepTurnover[step][3])

    return time, arr2DCounts.T, arr2DEnzymeStepTurnover.T, arr2DEnzymePercentSubstrateSaturation.T, arr2DEnzymePercentMaxActivty.T, enzymeLegend, metaboliteLegend







