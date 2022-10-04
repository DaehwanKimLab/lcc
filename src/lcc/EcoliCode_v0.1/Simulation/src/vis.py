import math
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from heapq import heapify
from random import randint
from plotly.subplots import make_subplots
# Custom
from modelEq import saturation, alloAct, alloInhib

## Austin 
def generateRGBList(n:int) -> list:
    ''' generate n rgb strings and return the list for plotly colors'''
    outputRGBs = [f'rgb({randint(0,255)},{randint(0,255)},{randint(0,255)})' for val in range(n)]
    return outputRGBs

def pplot(
    Time,
    metaboliteCountData:dict, 
    qssa:dict,
    yAxisScaleMetabolite:str = 'log',
    showQSSA:bool = True, 
    deltaT0Scale:bool = False,
    excludedMolecules: list = [],
    ):
    """
    
    """
    FORMAT_TITLEFONT = dict(family = 'Arial',size = 16, color = 'rgb(0,0,0)')

    X = Time

    # Check for excluded molecules:
    for mol in excludedMolecules:
        if mol in metaboliteCountData.keys():
            del metaboliteCountData[mol]
    
    ### TEMPORARY QUICKFIX:
    # Ensure no "0" values (for the sake of log scale):
    for key in metaboliteCountData.keys():
        if metaboliteCountData[key][0] == 0:
            metaboliteCountData[key][0] = 1e-18

    if deltaT0Scale:
        metaboliteCountData = {key:metaboliteCountData[key] - metaboliteCountData[key][0]}

    Y = [metaboliteCountData]

    metaboliteColors = {key:generateRGBList(1) for key in metaboliteCountData.keys()}
    
    fig = make_subplots(
        rows = 1, cols = 6, 
        specs = [
            [{"colspan":6}, None, None, None, None, None],
 #           [{"colspan":3}, None, None,{"colspan":3}, None, None],
 #           [{"colspan":3}, None, None,{"colspan":3}, None, None],
        ],
        shared_xaxes= True,
        subplot_titles=("Molecule Concentrations",)
    )


    # [metabolites]
    for key in metaboliteCountData.keys():
        fig.add_trace(
            # The line
            go.Scatter(
                x = X, y = Y[0][key],
                mode = 'lines', name = key, text = key,
                line = dict(color = metaboliteColors[key][0]),
                connectgaps = True,
            ), row = 1, col = 1)
    
        fig.add_trace(go.Scatter(
            x=[( X[0]+X[-1] + np.random.randint(-100,100)) // 2 ], # Add variable name (with a bit of random jitter so they don't overlap)
            y=[Y[0][key][len(X) // 2]],
            mode='text', name = key, text = key,
            marker = dict(color = metaboliteColors[key][0]),
            showlegend = False
        ), row = 1, col = 1)

        # Points at t0 and tn (start and end)
        fig.add_trace(go.Scatter(
            x=[X[0], X[-1]],
            y=[Y[0][key][0], Y[0][key][-1]],
            mode='markers', name = key, text = key,
            marker = dict(color = metaboliteColors[key][0]),
            showlegend = False
        ), row = 1, col = 1)

    if showQSSA:
        for key in qssa:
            # If violation:
            if sum(qssa[key]):
                # for each qssa, get the boundaries of each violation
                qssaViolations = [i for i in range(1,len(qssa[key])) if qssa[key][i] != qssa[key][i-1]]
                # Ensure the violations have start and end:
                if len(qssaViolations) % 2 != 0:
                    qssaViolations.append(len(qssa[key]))
                heapify(qssaViolations)

                # Then plot the violations with rects
                while qssaViolations:
                    l = qssaViolations.pop()
                    r = qssaViolations.pop()
                    fig.add_vrect(x0=l, x1=r, line_width=0, fillcolor=metaboliteColors[key.split("_")[0]][0], opacity=0.2, annotation_text=key, annotation_font_color="red")


    fig.update_xaxes(title_text = "steps", titlefont = FORMAT_TITLEFONT,
            showline=True,
            showgrid=True,
            showticklabels=True,
            linecolor='rgb(0, 0, 0)',
            linewidth=2, 
            ticks='outside',
            tickfont=dict(
                family='Arial',
                size=16,
                color='rgb(0, 0, 0)',
            ))

    fig.update_yaxes(title_text = "Concentration (uM)", 
            titlefont = FORMAT_TITLEFONT,
            type = yAxisScaleMetabolite,
            # TODO: Implement dynamic range
            range = [-20, 5],
            showline=True,
            showgrid=True,
            showticklabels=True,
            linecolor='rgb(0, 0, 0)',
            linewidth=2, 
            ticks='outside',
            tickfont=dict(
                family='Arial',
                size=16,
                color='rgb(0, 0, 0)',
            ), row = 1, col = 1)  

    fig.update_layout(
        autosize=True,   
        legend = dict(
            orientation = "h",
        ),
        plot_bgcolor='white'
        )

    if yAxisScaleMetabolite == 'log':
        fig.add_hline(y = 1.66e-18, line_dash = "dot",line_color = "red")
        fig.add_annotation(x = max(X), y = math.log10(0.5e-18), text = "Zero Concentration Threshold", showarrow = False)

    return fig

def pplot_sat_sim(
    enzName, XName, YName, YType,
    Time,
    counts,
    km,
    kreg,
    start = -12,
    stop = 8,
    stepsPerLog10 = 10
    ):
    # Set up the scale
    conc = [10**(i/stepsPerLog10) for i in range(start*stepsPerLog10,stop*stepsPerLog10)]

    if YType == 'sat':
        Z = [[saturation(j, km[f'{enzName}_{XName}'])*saturation(i, km[f'{enzName}_{YName}'])*100 for j in conc] for i in conc]
    elif YType == 'AlloI':
        Z = [[saturation(j, km[f'{enzName}_{XName}'])*alloInhib(i, kreg[f'Ki_{YName}_{enzName}'])*100 for j in conc] for i in conc]
    elif YType == 'AlloA':
        Z = [[saturation(j, km[f'{enzName}_{XName}'])*alloAct(i, kreg[f'Ka_{YName}_{enzName}'])*100 for j in conc] for i in conc]
    else:
        print(f"ERROR: the YType Value provided was not valid. \nCurrent Y-Type:{YType}\nAccepted Y-Types: sat, alloI, AlloA")

    XCounts = counts[XName]
    YCounts = counts[YName]

    fig = go.Figure(data = go.Contour( 
        z = Z, 
        x = conc,
        y = conc,
        colorscale = 'viridis',
        colorbar = dict(nticks = 20, title = f'{enzName} Relative Activity (%)', titleside='right'),
        line_width = 1,
        contours = dict(start = 0, end = 100, size = 5, showlabels=True, labelfont= dict(size = 12, color = 'white'))
        ))
    
    # Add simulation line
    fig.add_trace(
        go.Scatter(
                x = XCounts, y = YCounts,
                mode = 'lines',name = "Time: ", text = Time,
                line = dict(color = 'red'),
                connectgaps = True,
            ), 
        )


    fig.update_xaxes(title_text = f'{XName} (uM)',type = 'log', range= [start, stop])
    fig.update_yaxes(title_text = f'{YName} (uM)', type = 'log', range= [start, stop])
    
    return fig
