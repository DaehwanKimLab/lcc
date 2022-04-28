# Just getting the dash tutorial up and running.
#https://dash.plotly.com/layout 
#https://dash.plotly.com/basic-callbacks

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

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
from fun import preProcessSimOut, Plotly_Dynamics, CompileLPP, extractKVals
from functools import reduce

app = Dash(__name__)

## Apparently this is bad practice
df = pd.read_csv('src\lcc\SimOut.tsv', sep = '\t', index_col = False)
df = preProcessSimOut(df)

app.config.suppress_callback_exceptions = True

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
    html.Button('Compile L++', id='compile-button', n_clicks=0),
    html.Br(),
    html.Div(id = 'k1-output', children='Compile'),
    html.Br(),
    html.Button('get kcat', id='kcat-button', n_clicks=0),
    html.Div(id = 'k-val-container'),
    html.Div(dcc.Input(id = 'empty')),
    html.Br(),
    dcc.Input(
        id = 'lcc-input',
        type = 'text',
        value = 'src\lcc\denovoModels\glucose_metabolism_lowRes_abridged.lpp'
    ),
    html.Button('Load L++', id='load-lpp-button', n_clicks=0),
    html.Br(),
    html.Div([
        dcc.Markdown(id = 'lcc-text', children = 'lccFile')
    ]),
    dcc.Store(id = 'protein-k-values')
])

@app.callback(
    [Output('lcc-text', 'children'),
    Output('protein-k-values', 'data')],
    Input('load-lpp-button', 'n_clicks'),
    State('lcc-input', 'value'))
def update_lcctext(n_clicks,lcc_filepath):
    lpp = open(lcc_filepath, 'r')
    lines = lpp.readlines()
    lpp.close()
    return lines, extractKVals(lines)


@app.callback(
    Output('k-val-container', 'children'),
    Input('kcat-button', 'n_clicks'),
    State('protein-k-values', 'data'))
def update_kvals(button,data):
    if not data: 
        return None
    
    out = []
    for protein in data.keys():
        out.append(html.H4(protein, style={'display':'inline-block','margin-right':20}))
        out.append(dcc.Input(id = f'Protein_{protein}', type='number', value = float(data[protein]['kcat'])))
    return html.Div(out)

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
