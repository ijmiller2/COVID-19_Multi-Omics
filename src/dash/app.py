
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from data import get_metabolomics_data
from plot import biomolecule_bar, pca_plot

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# load metabolomics data matrix
metabolomics_df = get_metabolomics_data()

# load combined data matrix (with clinical metadata)
combined_df = get_metabolomics_data(with_metadata=True)

# define static pca plot
pca_figure = pca_plot(metabolomics_df, combined_df)

available_metabolobites = metabolomics_df.columns.tolist()

app.layout = html.Div([

    dcc.Graph(figure=pca_figure),

    html.Div([

        html.Div([
            html.Label('Select Biomolecule'),
            dcc.Dropdown(
                id='biomolecule_id',
                options=[{'label': i, 'value': i} for i in available_metabolobites],
                value=available_metabolobites[0]
            ),
        ],

        style={'width': '48%', 'display': 'inline-block'}),

    ]),

    dcc.Graph(id='biomolecule-barplot')

], style={'columnCount': 2})

@app.callback(
    Output('biomolecule-barplot', 'figure'),
    [Input('biomolecule_id', 'value')])

def update_biomolecule_barplot(biomolecule_id):

    x = metabolomics_df.index
    y = biomolecule_id
    fig = biomolecule_bar(metabolomics_df, x, y)

    return fig

#@app.callback(
#    Output('pca', 'figure'),
#    Input()
#)
#def update_pca(combined_df):
#
#    fig = pca(combined_df)
#
#    return fig

if __name__ == '__main__':
    app.run_server(debug=True)
