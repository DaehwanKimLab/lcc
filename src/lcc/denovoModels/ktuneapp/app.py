# Just getting the dash tutorial up and running.
#https://dash.plotly.com/layout 
#https://dash.plotly.com/basic-callbacks

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from operator import index
from dash import Dash, html, dcc, Input, Output, State
from matplotlib.pyplot import legend
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from random import randint
import os

## 
def preProcessSimOut(pdfSimOut):
    ''' preprocessing of simout.tsv
    1) drop the 'vol' and 'Pseudo' columns for now as they won't be used/needed
    
    SHORTCUT: it will change in the future.
    select only the columns I actually want 
    
    2) melt everything else together for now
    '''
    #pdf = pdfSimOut.drop(labels = ['Vol', 'Pseudo'])  
    pdf = pdfSimOut.loc[:,['SimStep', 'ADP', 'ATP', 'NAD', 'NADH', 'P', 'pyruvate', 'CoA', 'acetylCoA', 'FAD', 'FADH2']]
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
            go.Scatter(
                x = X, y = Y[i],
                mode = 'lines', name = Legend[i],
                line = dict(color = colors[i]),
                connectgaps = True,
                ))
        # endpoints
        fig.add_trace(go.Scatter(
            x=[X[0], X[-1]],
            y=[Y[i][0], Y[i][-1]],
            mode='markers', name = Legend[i],
            marker = dict(color = colors[i]),
            showlegend = False
        ))

    fig.update_layout(
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

NA = 6.0221409e+23
COMPILECLICKS = 0
K1 = 0

app = Dash(__name__)

df = pd.read_csv('src\lcc\SimOut.tsv', sep = '\t', index_col = False)
df = preProcessSimOut(df)

app.layout = html.Div([
    dcc.Graph(id='graph-with-slider'),
    dcc.RangeSlider(
        df['SimStep'].min(),
        df['SimStep'].max(),
        step=None,
        value=[df['SimStep'].min(), df['SimStep'].max()],
        marks={str(step): str(step) for step in df['SimStep'].unique()},
        id='step-slider'
    ),
    dcc.Input(
        id = 'k1-input',
        type = 'number',
        children = 'compile',
        value = 0
    ),
    html.Button('Compile', id='compile-button', n_clicks=0),
    html.Br(),
    html.Div(id = 'k1-output', children='Compile')
])


@app.callback(
    Output('k1-output', 'children'),
    Input('compile-button', 'n_clicks'),
    State('k1-input', 'value'))
def update_k1(n_clicks, k1_value):
    return f'k1 = {k1_value}'

@app.callback(
    Output('graph-with-slider', 'figure'),
    Input('step-slider', 'value'))
def update_figure(selected_step):
    filtered_df = df[df['SimStep'] >= selected_step[0]]
    filtered_df = filtered_df[filtered_df['SimStep'] <= selected_step[1]]
    
    # convert data to the Plotly_Dynamics format
    Legend = filtered_df['molecule'].unique()
    X = filtered_df['SimStep'].to_numpy()
    Y = [filtered_df['concentration'].where(filtered_df['molecule'] == molecule).dropna().rename(molecule).to_numpy() for molecule in Legend]
    
    fig = Plotly_Dynamics(Title='', Time = X, Data = Y, Legend = Legend)

    fig.update_layout(transition_duration=500)

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
