
from fun import Plotly_Simulation_Overview
from simplifiedMedGlyc import run

if __name__ == '__main__':
    
    time, arr2DCounts, arr2DEnzymeStepTurnover, arr2DEnzymePercentSubstrateSaturation, arr2DEnzymePercentMaxActivty, enzymeLegend, metaboliteLegend = run()

    fig=Plotly_Simulation_Overview(
            Time = time,
            metaboliteCountData = arr2DCounts,
            enzymeTurnoverData = arr2DEnzymeStepTurnover,
            enzymeSubstrateSaturationData= arr2DEnzymePercentSubstrateSaturation,
            enzymePercentMaxActivityData= arr2DEnzymePercentMaxActivty,
            enzymeLegend = enzymeLegend,
            metaboliteLegend = metaboliteLegend)
    
    fig.show()
