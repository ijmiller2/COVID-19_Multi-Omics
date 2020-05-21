
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from data import get_metabolomics_data
from plot import biomolecule_bar, pca_plot, boxplot

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
                # only passing in quant value columns
                value=available_metabolobites[0]
            ),
        ],

        style={'width': '49%', 'display': 'inline-block',
        'border': 'thin lightgrey solid',
        'backgroundColor': 'rgb(250, 250, 250)',
        'padding': '10px 5px'}),

    ]),

    html.Div(
        [html.Div(dcc.Graph(id='biomolecule-barplot', className="six columns")),
        html.Div(dcc.Graph(id='biomolecule-boxplot', className="six columns"))],
        className="row"),

])

@app.callback(
    Output('biomolecule-barplot', 'figure'),
    [Input('biomolecule_id', 'value')])

def update_biomolecule_barplot(biomolecule_id):

    biomolecule_name = biomolecule_id
    x = metabolomics_df.index
    y = biomolecule_id
    fig = biomolecule_bar(combined_df, x, y, biomolecule_name)

    return fig

@app.callback(
    Output('biomolecule-boxplot', 'figure'),
    [Input('biomolecule_id', 'value')])

def update_biomolecule_boxplot(biomolecule_id):

    biomolecule_name = biomolecule_id
    x = combined_df.index
    y = biomolecule_id
    fig = boxplot_figure = boxplot(combined_df, biomolecule_name)

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
