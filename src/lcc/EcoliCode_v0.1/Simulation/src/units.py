from scipy.constants import Avogadro

dictSIPrefix = {
    'none' : 0,
    'milli' : -3,
    'micro': -6,
    'nano': -9,
    'pico': -12,
    'femto': -15,
    'atto': -18,
    'zepto': -21
}

def mol2cnt(inVal:float, siPrefix='none'):
    """Convert moles to counts. If an si-prefix is provided, scale counts based on prefix"""
    return inVal*Avogadro*10**dictSIPrefix[siPrefix]
def cnt2mol(inVal:float, siPrefix='none'):
    """Convert counts to moles. If an si-prefix is provided, scale moles based on prefix"""
    return inVal/Avogadro/10**dictSIPrefix[siPrefix]



def umoltoCountPerEcoli(inVal, Vol = 0.6e-6):
    """Convert input in umol to molecule counts per ecoli (i.e. scaled to ecoli volume)"""
    return inVal*10**-6*Avogadro/Vol

def counttouMPerEcoli(inVal, Vol = 0.6e-6):
    """Convert input molecule counts to uM"""
    return inVal / Avogadro / 1e-6 / Vol
