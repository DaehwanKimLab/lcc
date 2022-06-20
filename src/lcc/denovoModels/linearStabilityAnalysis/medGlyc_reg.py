import numpy as np
from scipy.constants import Avogadro

# Micromolar conversion 
def uMtoCount(inVal:float):
    return inVal*10**-6*Avogadro
def counttouM(inVal:float):
    return inVal / Avogadro / 10**-6

# constants
SIM_DURATION_SECONDS = 0.6
#SIM_DURATION_SECONDS = 600
#SIM_DURATION_SECONDS = 2000
STEPS_PER_SECOND = 100

TOTAL_STEPS = SIM_DURATION_SECONDS * STEPS_PER_SECOND

# kcats
KCAT_HEXOKINASE =  4000
KCAT_PFK = 1.7e4
KCAT_PFK_ACTIVE = 1.7e4
KCAT_PFK_INACTIVE = 1.7e4
KCAT_PK = 170
##
KCAT_SUMMARY_MAKE_PEPASE = 5e3 # Arbitrary
# ks
K_SUMMARY_MAKE_PEPASE = 1e4

# KMs
KM_HEXOKINASE_GLUCOSE = uMtoCount(5)  # yeast
KM_HK_ATP = uMtoCount(100)            # close guesses from lit (Yeast)
KM_PFK_ACTIVE_F6P = uMtoCount(0.21)
KM_PFK_ACTIVE_ATP = uMtoCount(60)
KM_PFK_INACTIVE = uMtoCount(3.1e7)
KM_PK_PEP = uMtoCount(700)
KM_PK_ADP = uMtoCount(50) # Guesstimate from lit, most bact seem around this

#####
KM_SUMMARY_MAKE_PEPASE = uMtoCount(3) # just guessing.

# Concentrations
## Metabolites
CONC_GLUCOSE_T0 = uMtoCount(2200)
CONC_ADP_T0 = uMtoCount(560)        # https://www.nature.com/articles/nchembio.186.pdf
CONC_ATP_T0 = uMtoCount(9600)       # https://www.nature.com/articles/nchembio.186.pdf
## Intermediates
CONC_ATP = CONC_ATP_T0              # Goal: 9.6mM - https://www.nature.com/articles/nchembio.186.pdf
CONC_ADP = CONC_ADP_T0              # Goal: 0.56mM - https://www.nature.com/articles/nchembio.186.pdf
CONC_F6P = 0#uMtoCount(15)          # Goal: 8.8mM - https://www.nature.com/articles/nchembio.186.pdf
CONC_F16P = 0                       # Goal: 15.0mM - https://www.nature.com/articles/nchembio.186.pdf
CONC_PEP = 0                        # Goal: 0.18mM - https://www.nature.com/articles/nchembio.186.pdf
CONC_PYRUVATE = 0                   
## Enzymes
CONC_PFK = uMtoCount(33e-3) # 33nM https://masspy.readthedocs.io/en/v0.1.6/education/sb2/chapters/sb2_chapter14.html
CONC_PFK_ACTIVE = 0
CONC_PFK_INACTIVE = 0
CONC_HEXOKINASE = uMtoCount(50e-3)
CONC_SUMMARY_MAKE_PEPASE = uMtoCount(6e-2)
CONC_PYRUVATE_KINASE = uMtoCount(1)

## Current exploratory values
KA_ADP_PFK = uMtoCount(25)
#KI_PEP_PFK = uMtoCount(260)             # use Keft_pep from file:///D:/Documents/Gradschool/KimLab/lit/someEcoliODEs.pdf
KI_PEP_PFK = uMtoCount(10)
#KA_F16BP_PK = uMtoCount(390)            # file:///D:/Documents/Gradschool/KimLab/lit/someEcoliODEs.pdf
KA_F16BP_PK = uMtoCount(3)            # file:///D:/Documents/Gradschool/KimLab/lit/someEcoliODEs.pdf
KI_ATP_PK = uMtoCount(84000)            # Ki_r of ATP ~ 84mM : file:///D:/Documents/Gradschool/KimLab/lit/someEcoliODEs.pdf

KI_F6P_HK = uMtoCount(8000)


PFK_STATE_TRANSITION_CONSTANT = 14.4   #https://www.researchgate.net/publication/23230008_Kinetic_Model_of_Phosphofructokinase-1_from_Escherichia_coli?enrichId=rgreq-40a860a066e1fe2f49a81d0a30616e86-XXX&enrichSource=Y292ZXJQYWdlOzIzMjMwMDA4O0FTOjIyMTQ0NDY4MTI3NzQ0MEAxNDI5ODA3OTM0NDgw&el=1_x_2&_esc=publicationCoverPdf

useOneStatePFK = True

def saturation(conc, KM):
    return (conc / (KM + conc))
def activate(conc, Ka):
    return(1 + (conc/Ka))
def competAlloOneActOneInhib(actConc, Ka, inhibConc, Ki):
    return((1 + actConc/Ka)/(1 + inhibConc/Ki))


# Attempting to base
def hexokinase(Conc_Substrate1 = CONC_GLUCOSE_T0, Conc_Substrate2 = CONC_ATP,Conc_Inhibitor = CONC_F6P, Ki =KI_F6P_HK , Conc_Enzyme = CONC_HEXOKINASE, kcat = KCAT_HEXOKINASE, KM1 = KM_HEXOKINASE_GLUCOSE,KM2 =  KM_HK_ATP):
    # Note: glucose is not decreasing...
    return (kcat * Conc_Enzyme / activate(Conc_Inhibitor, Ki) * saturation(Conc_Substrate1, KM1) * saturation(Conc_Substrate2, KM2))

if useOneStatePFK:
    #def pfk(Conc_Substrate1 = CONC_F6P, Conc_Substrate2 = CONC_ATP, Conc_Activator = CONC_ADP,Conc_Inhibitor = CONC_PEP, Ka = KA_ADP_PFK, Ki = KI_PEP_PFK, Conc_Enzyme = CONC_PFK, kcat = KCAT_PFK_ACTIVE, KM1 = KM_PFK_ACTIVE_F6P, KM2 = KM_PFK_ACTIVE_ATP):
    #    return (kcat * Conc_Enzyme * activate(Conc_Activator, Ka) * saturation(Conc_Substrate1, KM1) * saturation(Conc_Substrate2, KM = KM2))
    def pfk(Conc_Substrate1 = CONC_F6P, Conc_Substrate2 = CONC_ATP, Conc_Activator = CONC_ADP,Conc_Inhibitor = CONC_PEP, Ka = KA_ADP_PFK, Ki = KI_PEP_PFK, Conc_Enzyme = CONC_PFK, kcat = KCAT_PFK_ACTIVE, KM1 = KM_PFK_ACTIVE_F6P, KM2 = KM_PFK_ACTIVE_ATP):
        return (kcat * Conc_Enzyme / activate(Conc_Inhibitor, Ki) * saturation(Conc_Substrate1, KM1) * saturation(Conc_Substrate2, KM = KM2))
    #def pfk(Conc_Substrate1 = CONC_F6P, Conc_Substrate2 = CONC_ATP, Conc_Activator = CONC_ADP,Conc_Inhibitor = CONC_PEP, Ka = KA_ADP_PFK, Ki = KI_PEP_PFK, Conc_Enzyme = CONC_PFK, kcat = KCAT_PFK_ACTIVE, KM1 = KM_PFK_ACTIVE_F6P, KM2 = KM_PFK_ACTIVE_ATP):
    #    return (kcat * Conc_Enzyme * competAlloOneActOneInhib(Conc_Activator, Ka, Conc_Inhibitor, Ki) * saturation(Conc_Substrate1, KM1) * saturation(Conc_Substrate2, KM = KM2))
else:
    def pfk_active(Conc_Substrate = CONC_F6P, Conc_Enzyme = CONC_PFK_ACTIVE, kcat = KCAT_PFK_ACTIVE, KM = KM_PFK_ACTIVE_F6P):
        return (kcat * Conc_Enzyme *  saturation(Conc_Substrate, KM))
    def pfk_inactive(Conc_Substrate = CONC_F6P, Conc_Enzyme = CONC_PFK_INACTIVE, kcat = KCAT_PFK_INACTIVE, KM = KM_PFK_INACTIVE):
        return (kcat * Conc_Enzyme * saturation(Conc_Substrate, KM))

def summaryMakePEPase(Conc_Substrate = CONC_F16P,Conc_Enzyme = CONC_SUMMARY_MAKE_PEPASE, kcat = KCAT_SUMMARY_MAKE_PEPASE, KM = KM_SUMMARY_MAKE_PEPASE):
    return(kcat * Conc_Enzyme * saturation(Conc_Substrate, KM))

def pyruvateKinase(Conc_Substrate1 = CONC_PEP, Conc_Substrate2 = CONC_ADP, Conc_Inhibitor = CONC_ATP,Conc_Activator = CONC_F16P, Ka= KA_F16BP_PK, Ki = KI_ATP_PK, Conc_Enzyme = CONC_PYRUVATE_KINASE, kcat = KCAT_PK, KM1 = KM_PK_PEP, KM2 = KM_PK_ADP):
    return (kcat * Conc_Enzyme / (1 + (Conc_Inhibitor / Ki)) * activate(Conc_Activator, Ka) * saturation(Conc_Substrate1, KM1)* saturation(Conc_Substrate2, KM2))

def run():
    # reimplement simulation: 
    metaboliteLegend = ['F6P', 'F16BP', 'PEP', 'pyruvate', 'ATP', 'ADP']
    enzymeLegend = ['HK', 'PFK', 'PEPSYN', 'PK']
    enzymeMetaboliteLegend = ['HK-Gluc', 'PFK-F6P', 'PFK-ATP', 'PEPSYN-F16BP', 'PK-PEP', 'PK-ADP', 'HK-ATP']
    time = np.linspace(0, SIM_DURATION_SECONDS, TOTAL_STEPS)

    # For now only worry about F6P, F16BP, PEP, pyruvate
    arr2DCounts= np.zeros([TOTAL_STEPS, len(metaboliteLegend)]) 
    # Init ATP/ADP Counts
    arr2DCounts[0][4] = CONC_ATP_T0
    arr2DCounts[0][5] = CONC_ADP_T0

    arr2DEnzymePercentSubstrateSaturation = np.zeros([TOTAL_STEPS, len(enzymeMetaboliteLegend)])
    arr2DEnzymePercentMaxActivty = np.zeros([TOTAL_STEPS, len(enzymeLegend)])
    arr2DEnzymeStepTurnover = np.zeros([TOTAL_STEPS, len(enzymeLegend)])

    for step in range(1,TOTAL_STEPS):
        # Get turnover for each enzyme
        arr2DEnzymeStepTurnover[step][0] = hexokinase(Conc_Substrate2=arr2DCounts[step-1][4], Conc_Inhibitor=arr2DCounts[step-1][0]) / STEPS_PER_SECOND                               # In: Gluc + ATP  ; Out: F6P
        arr2DEnzymeStepTurnover[step][1] = pfk(arr2DCounts[step-1][0], arr2DCounts[step-1][4], arr2DCounts[step-1][5],arr2DCounts[step-1][2]) / STEPS_PER_SECOND                               # In: F6P + ATP (ADP reg); Out: F16BP
        arr2DEnzymeStepTurnover[step][2] = summaryMakePEPase(arr2DCounts[step-1][1]) / STEPS_PER_SECOND                                         # In: F16BP ; Out: PEP
        arr2DEnzymeStepTurnover[step][3] = pyruvateKinase(arr2DCounts[step-1][2], arr2DCounts[step-1][5], arr2DCounts[step-1][4]) / STEPS_PER_SECOND      # In: PEP  (ATP reg) ; Out: Pyruvate

        # Get percent substrate Saturation:
        arr2DEnzymePercentSubstrateSaturation[step][0] = saturation(CONC_GLUCOSE_T0, KM_HEXOKINASE_GLUCOSE) * 100           # HK-Gluc
        arr2DEnzymePercentSubstrateSaturation[step][6] = saturation(arr2DCounts[step-1][5], KM_HK_ATP) * 100                # HK-ATP
        arr2DEnzymePercentSubstrateSaturation[step][1] = saturation(arr2DCounts[step-1][0], KM_PFK_ACTIVE_F6P) * 100        # PFK-F6P
        arr2DEnzymePercentSubstrateSaturation[step][2] = saturation(arr2DCounts[step-1][4], KM_PFK_ACTIVE_ATP) * 100        # PFK-ATP

        arr2DEnzymePercentSubstrateSaturation[step][3] = saturation(arr2DCounts[step-1][1], KM_SUMMARY_MAKE_PEPASE) * 100
        arr2DEnzymePercentSubstrateSaturation[step][4] = saturation(arr2DCounts[step-1][2], KM_PK_PEP) * 100
        arr2DEnzymePercentSubstrateSaturation[step][5] = saturation(arr2DCounts[step-1][5], KM_PK_ADP) * 100

        #arr2DEnzymePercentSubstrateSaturation[step][7] = competAlloOneActOneInhib(arr2DCounts[step-1][5], KA_ADP_PFK, arr2DCounts[step-1][2], KI_PEP_PFK) * 100

        # Get Percent Max Activity
        arr2DEnzymePercentMaxActivty[step][0] = round(arr2DEnzymeStepTurnover[step][0] / (KCAT_HEXOKINASE * CONC_HEXOKINASE) * 100 * STEPS_PER_SECOND, 2)
        arr2DEnzymePercentMaxActivty[step][1] = round(arr2DEnzymeStepTurnover[step][1] / (KCAT_PFK * CONC_PFK) * 100 * STEPS_PER_SECOND, 2)
        arr2DEnzymePercentMaxActivty[step][2] = round(arr2DEnzymeStepTurnover[step][2] / (KCAT_SUMMARY_MAKE_PEPASE * CONC_SUMMARY_MAKE_PEPASE) * 100 * STEPS_PER_SECOND, 2)
        arr2DEnzymePercentMaxActivty[step][3] = round(arr2DEnzymeStepTurnover[step][3] / (KCAT_PK * CONC_PYRUVATE_KINASE) * 100 * STEPS_PER_SECOND, 2)

        # Update counts
        arr2DCounts[step][0] = max(0, arr2DCounts[step-1][0] + arr2DEnzymeStepTurnover[step][0] - arr2DEnzymeStepTurnover[step][1])
        arr2DCounts[step][1] = max(0, arr2DCounts[step-1][1] + arr2DEnzymeStepTurnover[step][1] - arr2DEnzymeStepTurnover[step][2])
        arr2DCounts[step][2] = max(0, arr2DCounts[step-1][2] + 2*arr2DEnzymeStepTurnover[step][2] - arr2DEnzymeStepTurnover[step][3])
        
        #arr2DCounts[step][3] = max(0, arr2DCounts[step-1][3] + arr2DEnzymeStepTurnover[step][3])
        # Pyruvate -- ensure it doesn't grow beyond 10mM
        arr2DCounts[step][3] = min(uMtoCount(10000), max(0, arr2DCounts[step-1][3] + arr2DEnzymeStepTurnover[step][3]))
        
        # ATP
        arr2DCounts[step][4] = max(0, arr2DCounts[step-1][4] + 2*arr2DEnzymeStepTurnover[step][3] - arr2DEnzymeStepTurnover[step][1] - arr2DEnzymeStepTurnover[step][0])
        # ADP
        arr2DCounts[step][5] = max(0, arr2DCounts[step-1][5] - 2*arr2DEnzymeStepTurnover[step][3] + arr2DEnzymeStepTurnover[step][1] + arr2DEnzymeStepTurnover[step][0])

    return time, arr2DCounts.T, arr2DEnzymeStepTurnover.T, arr2DEnzymePercentSubstrateSaturation.T, arr2DEnzymePercentMaxActivty.T, enzymeLegend, metaboliteLegend, enzymeMetaboliteLegend







