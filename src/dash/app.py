
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import plotly.express as px

from data import get_metabolomics_data

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# load metabolomics data matrix
metabolomics_df = get_metabolomics_data()

available_metabolobites = metabolomics_df.columns.tolist()

app.layout = html.Div([
    html.Div([

        html.Div([
            dcc.Dropdown(
                id='biomolecule_id',
                options=[{'label': i, 'value': i} for i in available_metabolobites],
                value=available_metabolobites[0]
            ),
        ],

        style={'width': '48%', 'display': 'inline-block'}),

    ]),

    dcc.Graph(id='biomolecule-barplot'),

])

@app.callback(
    Output('biomolecule-barplot', 'figure'),
    [Input('biomolecule_id', 'value')])

def update_graph(biomolecule_id):

    fig = px.bar(metabolomics_df, x=metabolomics_df.index, y=biomolecule_id)
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)

    print(available_metabolobites)
