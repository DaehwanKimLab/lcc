
from fun import Plotly_Simulation_Overview
# from simplifiedMedGlyc import run
#from medGlyc_atpadp_reg import run
from medGlyc_reg import run

if __name__ == '__main__':
    
    time, arr2DCounts, arr2DEnzymeStepTurnover, arr2DEnzymePercentSubstrateSaturation, arr2DEnzymePercentMaxActivty, enzymeLegend, metaboliteLegend, enzymeMetaboliteLegend = run()

    fig=Plotly_Simulation_Overview(
            Time = time,
            metaboliteCountData = arr2DCounts,
            enzymeTurnoverData = arr2DEnzymeStepTurnover,
            enzymeSubstrateSaturationData= arr2DEnzymePercentSubstrateSaturation,
            enzymePercentMaxActivityData= arr2DEnzymePercentMaxActivty,
            enzymeLegend = enzymeLegend,
            metaboliteLegend = metaboliteLegend,
            enzymeMetaboliteLegend = enzymeMetaboliteLegend)
    
    fig.show()
