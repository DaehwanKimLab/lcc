# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

# Tweaking the dash tutorials has been the basis for much of this app (I've listed two but all 5 parts are decent)
#https://dash.plotly.com/layout 
#https://dash.plotly.com/basic-callbacks


# from operator import index
from turtle import ht
from dash import Dash, html, dcc, Input, Output, State
# from matplotlib.pyplot import legend
# import plotly.express as px
# import plotly.graph_objects as go
import pandas as pd
import numpy as np
# from random import randint
# import os
from fun import preProcessSimOut, Plotly_Dynamics, CompileLPP, extractKVals, Plotly_Simulation_Overview
from simplifiedMedGlyc import run
from functools import reduce


app = Dash(__name__)

## Apparently this is bad practice
#df = pd.read_csv('src\lcc\SimOut.tsv', sep = '\t', index_col = False)
#df = preProcessSimOut(df)

time, arr2DCounts, arr2DEnzymeStepTurnover, arr2DEnzymePercentSubstrateSaturation, arr2DEnzymePercentMaxActivty, enzymeLegend, metaboliteLegend = run()

app.config.suppress_callback_exceptions = True

app.layout = html.Div([
    # The plot and time slider
    dcc.Graph(id='graph-with-slider'),
    dcc.RangeSlider(
        time.min(),
        time.max(),
        step=None,
        value=[time.min(), time.max()],
        #marks={str(step): str(step) for step in range(0, int(time.max()), 100)},
        id='step-slider'
    ),
])

@app.callback(
    Output('graph-with-slider', 'figure'),
    Input('step-slider', 'value'))
def update_figure(selected_step):
    
    # convert data to the Plotly_Dynamics format
    X = time[selected_step[0]:selected_step[1]]
    Y1 = arr2DCounts[selected_step[0]:selected_step[1], :]
    Y2 = arr2DEnzymeStepTurnover[selected_step[0]:selected_step[1], :]
    Y3 = arr2DEnzymePercentSubstrateSaturation[selected_step[0]:selected_step[1], :]
    Y4 = arr2DEnzymePercentMaxActivty[selected_step[0]:selected_step[1], :]
    L1 = enzymeLegend
    L2 = metaboliteLegend

    fig=Plotly_Simulation_Overview(
        Time = X,
        metaboliteCountData = Y1,
        enzymeTurnoverData = Y2,
        enzymeSubstrateSaturationData= Y3,
        enzymePercentMaxActivityData= Y4,
        enzymeLegend = L1,
        metaboliteLegend = L2)

    fig.update_layout(transition_duration=500)

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
