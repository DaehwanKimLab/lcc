import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from random import randint
import os
import re

NA = 6.0221409e+23

# Temporary preprocessing function.
def ppAll():
    return preProcessSimOut(pd.read_csv('src\lcc\SimOut.tsv', sep = '\t', index_col = False))

## 
def preProcessSimOut(pdfSimOut):
    ''' preprocessing of simout.tsv
    1) drop the 'vol' and 'Pseudo' columns for now as they won't be used/needed
    
    SHORTCUT: it will change in the future.
    select only the columns I actually want 
    
    2) melt everything else together for now
    '''
    #pdf = pdfSimOut.drop(labels = ['Vol', 'Pseudo'])  
    pdf = pdfSimOut.loc[:,['SimStep', 'ADP', 'ATP', 'F16BP', 'fructose6P', 'PEP']]
    pdf = pdf.melt(id_vars = ['SimStep'], value_vars = pdf.drop(labels = 'SimStep', axis = 1).columns, var_name = 'molecule', value_name = 'concentration')
    return pdf

## Austin 
def generateRGBList(n:int) -> list:
    ''' generate n rgb strings and return the list for plotly colors'''
    outputRGBs = [f'rgb({randint(0,255)},{randint(0,255)},{randint(0,255)})' for val in range(n)]
    return outputRGBs

### Austin-- convert matplotlib to plotly for easier labelling
## https://plotly.com/python/line-charts/#label-lines-with-annotations

def Plotly_Dynamics(Title, Time, Data, Legend):
    ''' Generate a plotly plot wiht time on x axis and colored line for each data'''
    X = Time
    Y = Data
    # Temp conversion to mM
    Y = np.multiply(Y, 1000)
    Y = np.divide(Y, NA)

    colors = generateRGBList(len(Legend))

    #assert len(X) == Y.shape[-1]

    # make figure
    fig = go.Figure()

    for i in range(len(Legend)):
        if Legend[i] == "Pseudo":
            continue
        fig.add_trace(
            # The line
            go.Scatter(
                x = X, y = Y[i],
                mode = 'lines', name = Legend[i],
                line = dict(color = colors[i]),
                connectgaps = True,
                ))
        # Points at t0 and tn (start and end)
        fig.add_trace(go.Scatter(
            x=[X[0], X[-1]],
            y=[Y[i][0], Y[i][-1]],
            mode='markers', name = Legend[i],
            marker = dict(color = colors[i]),
            showlegend = False
        ))

    # Formatting 
    fig.update_layout(
        # xaxis
        xaxis=dict(
            title = "Time (s)",
            titlefont = dict(
                family = 'Arial',
                size = 16,
                color = 'rgb(0,0,0)'
            ),
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
            ),
        ),
        #  yaxis
        yaxis=dict(
            title = "Concentration (mM)",
            titlefont = dict(
                family = 'Arial',
                size = 16,
                color = 'rgb(0,0,0)'
            ),
            zeroline = True,
            showgrid=True,
            showline=True,
            showticklabels=True,
            linecolor='rgb(0, 0, 0)',
            linewidth=2,
            ticks='outside',
            tickfont=dict(
                family='Arial',
                size=16,
                color='rgb(0, 0, 0)',
            )
        ),
        autosize=True,   
        legend = dict(
            orientation = "h",
        ),
        plot_bgcolor='white'
        )

    annotations = []
    # Adding labels
    for y_trace, label in zip(Y, Legend):
        # labeling the left_side of the plot
        annotations.append(dict(xref='paper', axref = 'x', x=0.1, y=y_trace[0]+.1,
                                    xanchor='right', yanchor='middle',
                                    text=' {:2g}mM'.format(round(y_trace[0], 2)),
                                    font=dict(family='Arial',
                                                size=16),
                                    showarrow=False))
        # labeling the right_side of the plot
        annotations.append(dict(xref='paper', axref = "x", x=0.95, y=y_trace[-1],
                                    xanchor='left', yanchor='middle',
                                    text=label +' {:2g}mM'.format(round(y_trace[-1], 2)),
                                    font=dict(family='Arial',
                                                size=16),
                                    showarrow=False))
    # Title
    annotations.append(dict(xref='paper', yref='paper', x=0.0, y=1.05,
                                xanchor='left', yanchor='bottom',
                                text='[Metabolites] over time',
                                font=dict(family='Arial',
                                            size=30,
                                            color='rgb(37,37,37)'),
                                showarrow=False))

    fig.update_layout(annotations=annotations)
    return fig

def CompileLPP():
    '''Compile the lpp file.
    Eventually will want the file I wish to compile to be selectable.    
    '''
    os.system("build\Debug\lcc.exe denovoModels\glucose_metabolism_lowRes_abridged.lpp")


def extractKVals(lstOfStrings:list):
    '''Extract kvalues and protein names from lpp file.
    
    Currently this is only designed/tested on the sample file
    (denovoModels\glucose_metabolism_lowRes_abridged.lpp)
    Essentially I just find lines with 'protein ', then assuming protein begins at index 0:
    read the protein name up to the '\('.
    kcat in the first comma, so do some splitting and extract the number.

    Returns a dictionary of format {'proteinName':{'kcat':kcatValue}}
    '''
    
    if not lstOfStrings: 
        return None
    outputDict = {'':{}}
    for l in lstOfStrings:
        if re.match('protein ', l):
            outputDict[l[8:re.search('\(', l).span()[0]]] = {'kcat':l.split(',')[1].split('=')[1].strip()}
    del outputDict[''] 
    return outputDict