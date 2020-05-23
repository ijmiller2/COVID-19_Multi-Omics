
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from data import get_metabolomics_data
from plot import pca_scores_plot, pca_loadings_plot, biomolecule_bar, boxplot

#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_stylesheets=[dbc.themes.BOOTSTRAP]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# load metabolomics data matrix
metabolomics_df = get_metabolomics_data()

# load combined data matrix (with clinical metadata)
combined_df = get_metabolomics_data(with_metadata=True)

# define static pca plot
pca_scores_figure = pca_scores_plot(metabolomics_df, combined_df)
pca_loadings_figure = pca_loadings_plot(metabolomics_df)

available_metabolobites = metabolomics_df.columns.tolist()

first_card = dbc.Card(
    [
        dbc.CardHeader("PCA Scores Plot",
                            style={"background-color":"#5bc0de",
                                    "font-weight":"bold",
                                    "font-size":"large"}),
        dbc.CardBody(dcc.Graph(figure=pca_scores_figure))

        ])

second_card = dbc.Card(
    [
        dbc.CardHeader("PCA Loadings Plot",
                            style={"background-color":"#5bc0de",
                                    "font-weight":"bold",
                                    "font-size":"large"}),
        dbc.CardBody(dcc.Graph(figure=pca_loadings_figure))
    ])

biomolecule_dropdown = dbc.Card(
    [
        dbc.CardHeader("Select Biomolecule",
                            style={"background-color":"#5bc0de",
                                        "font-weight":"bold",
                                        "font-size":"large"}),
        dbc.CardBody(
            dcc.Dropdown(
                id='biomolecule_id',
                options=[{'label': i, 'value': i} for i in available_metabolobites],
                # only passing in quant value columns
                value=available_metabolobites[0])
                )
    ])

third_card = dbc.Card(
    [
        dbc.CardHeader("Biomolecule Barplot",
                            style={"background-color":"#5bc0de",
                                        "font-weight":"bold",
                                        "font-size":"large"}),
        dbc.CardBody(dcc.Graph(id='biomolecule-barplot'))
    ])

fourth_card = dbc.Card(
    [
        dbc.CardHeader("Biomolecule Boxplot",
                            style={"background-color":"#5bc0de",
                                    "font-weight":"bold",
                                    "font-size":"large"}),
        dbc.CardBody(dcc.Graph(id='biomolecule-boxplot'))
    ])

COONLAB_LOGO="https://coonlabs.com/wp-content/uploads/2016/07/coon-logo-white.png"
navbar = dbc.NavbarSimple(
    children=[

        dbc.NavItem(dbc.NavLink("About", href="#")),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("Meet the Team", header=True),
                dbc.DropdownMenuItem("UW–Madison", href="https://coonlabs.com/"),
                dbc.DropdownMenuItem("Albany Medical Center", href="https://www.amc.edu/Profiles/jaitova.cfm"),
                dbc.DropdownMenuItem("Morgridge Institute for Research", href="https://morgridge.org/research/regenerative-biology/bioinformatics/")
            ],
            nav=True,
            in_navbar=True,
            label="More",
        ),
        html.Div(style={"width":"200px"}),
        html.A(html.Img(src=COONLAB_LOGO, style={"height":"40px", "float":"right"}
                ), href="https://coonlabs.com/"),

    ],
    brand="NIH National Center for Quantitative Biology of Complex Systems",
    brand_style={"font-size":"xx-large"},
    brand_href="https://www.ncqbcs.com/",
    color="#5bc0de",
    dark=True,
)

app.layout = dbc.Container([

    navbar,

    html.Hr(),

    dbc.Row(dbc.Col(html.H1("COVID-19 Multi-Omics Data Dashboard"), width={"size": 6, "offset": 3})),

    html.Hr(),

    dbc.Row([dbc.Col(first_card, md=6, align="center"), dbc.Col(second_card, md=6, align="center")], className="mb-3"),

    dbc.Row([dbc.Col(biomolecule_dropdown, width=4)], className="mb-3"),

    dbc.Row([dbc.Col(third_card, md=6, align="center"), dbc.Col(fourth_card, md=6, align="center")], className="mb-3")

], fluid=True)

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

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')
