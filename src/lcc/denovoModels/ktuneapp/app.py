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
from fun import preProcessSimOut, Plotly_Dynamics, CompileLPP

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
    html.Button('Compile L++', id='compile-button', n_clicks=0),
    html.Br(),
    html.Div(id = 'k1-output', children='Compile'),
    html.Br(),
    dcc.Input(
        id = 'lcc-input',
        type = 'text',
        value = 'src\lcc\denovoModels\glucose_metabolism_lowRes_abridged.lpp',
        children = 'lccFile'
    ),
    html.Br(),
    html.Div(id = 'lcc-output', children = 'lccFile')
    # html.Div([
    #     dcc.Markdown(id = 'lcc-text', children = 'lccFile')
    # ])
    
])

@app.callback(
    Output('lcc-output', 'lccFile'),
    Input('lcc-input', 'value'))
def update_lcctext(lcc_filepath):
    lpp = open(lcc_filepath, 'r')
    lines = lpp.readlines()
    lpp.close()
    return lines

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
