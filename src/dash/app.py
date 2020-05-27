
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from data import get_omics_data
from plot import pca_scores_plot, pca_loadings_plot, biomolecule_bar, boxplot

#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_stylesheets=[dbc.themes.BOOTSTRAP]

app = dash.Dash(
    __name__,
    external_stylesheets=external_stylesheets)
app.title = 'COVID-19 Multi-Omics'

# load metabolomics data matrix
#print("Loading metabolomics data...")
metabolomics_df, metabolomics_quant_range = get_omics_data(dataset='metabolomics', with_metadata=True)
#print("Metabolomics data shape: {}".format(metabolomics_df.shape))
#print("Loading lipidomics data...")
lipidomics_df, lipidomics_quant_range = get_omics_data(dataset='lipidomics', with_metadata=True)
#print("Lipidomics data shape: {}".format(lipidomics_df.shape))
#print("Loading proteomics data...")
proteomics_df, proteomics_quant_range = get_omics_data(dataset='proteomics', with_metadata=True)
#print("Proteomics data shape: {}".format(proteomics_df.shape))

available_datasets = ['Proteins', 'Lipids', 'Metabolites']
# start with proteomics data
available_biomolecules = proteomics_df.columns[:proteomics_quant_range].sort_values().tolist()

plotly_config = {"toImageButtonOptions":{'format':'svg',
                'filename': 'dash_plot'},
                "displaylogo": False,
                "displayModeBar":True}

first_card = dbc.Card(
    [
        dbc.CardHeader("PCA SCORES PLOT",
                            style={"background-color":"#5bc0de",
                                    "font-weight":"bold",
                                    "font-size":"large"}),
        dbc.CardBody(dcc.Graph(id='pca-scores-figure',
            config=plotly_config))

        ])

second_card = dbc.Card(
    [
        dbc.CardHeader("PCA LOADINGS PLOT",
                            style={"background-color":"#5bc0de",
                                    "font-weight":"bold",
                                    "font-size":"large"}),
        dbc.CardBody(dcc.Graph(id='pca-loadings-figure',
        config=plotly_config))
    ])

control_panel = dbc.Card(
    [
        dbc.CardHeader("CONTROL PANEL",
                            style={"background-color":"#5bc0de",
                                        "font-weight":"bold",
                                        "font-size":"large"}),
        dbc.CardBody(
            [html.P("Select Dataset", className="card-title", style={"font-weight":"bold"}),
            dcc.Dropdown(
                id='dataset_id',
                options=[{'label': i, 'value': i} for i in available_datasets],
                # only passing in quant value columns
                value=available_datasets[0]),
            html.Hr(),
            html.P("Select Biomolecule", className="card-title", style={"font-weight":"bold"}),

            # NOTE: This is dcc object not dbc
            dcc.Dropdown(
                id='biomolecule_id',
                options=[{'label': i, 'value': i} for i in available_biomolecules],
                # only passing in quant value columns
                value=available_biomolecules[0],
                className="dropdown-item p-0"),

            #dcc.Dropdown(id='biomolecule_id')

                ])
    ])

third_card = dbc.Card(
    [
        dbc.CardHeader("BIOMOLECULE BARPLOT",
                            style={"background-color":"#5bc0de",
                                        "font-weight":"bold",
                                        "font-size":"large"}),
        dbc.CardBody(dcc.Graph(id='biomolecule-barplot',
        config=plotly_config))
    ])

fourth_card = dbc.Card(
    [
        dbc.CardHeader("BIOMOLECULE BOXPLOT",
                            style={"background-color":"#5bc0de",
                                    "font-weight":"bold",
                                    "font-size":"large"}),
        dbc.CardBody(dcc.Graph(id='biomolecule-boxplot',
        config=plotly_config))
    ])

COONLAB_LOGO="https://coonlabs.com/wp-content/uploads/2016/07/coon-logo-white.png"
navbar = dbc.NavbarSimple(

        [

        dbc.NavItem(dbc.NavLink(html.Span(
                "About",
                id="tooltip-about",
                style={"cursor": "pointer"},
            ))),

        dbc.Tooltip(
        "Link to Preprint",
        target="tooltip-about"
        ),

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

        dbc.NavItem(html.Div(html.A(html.Img(src=COONLAB_LOGO,
                style={"height":"40px"}),
                href="https://coonlabs.com/"),
                className="d-none d-lg-block ml-4"))

        ],
    brand="NIH NATIONAL CENTER FOR QUANTITATIVE BIOLOGY OF COMPLEX SYSTEMS",
    brand_style={"font-size":"xx-large", "font-style":"italic"},
    brand_href="https://www.ncqbcs.com/",
    color="#5bc0de",
    dark=True,
    className="mt-1",
    fluid=True
)

app.layout = dbc.Container([

    navbar,

    html.Hr(),

    dbc.Row(dbc.Col(html.H1("COVID-19 Multi-Omics Data Dashboard"), width={"size": 6, "offset": 3})),

    html.Hr(),

    dbc.Row(
        [dbc.Col(
        dbc.Nav(
    [
        html.H3("TYPE OF ANALYSIS", style={"font-weight":"bold", "color":"black"}),

        dbc.NavItem(dbc.NavLink("PCA", active=True, href="#", style={"background-color":"grey"})),

        dbc.NavItem(dbc.NavLink(

            html.Span(
                    "Linear Regression",
                    id="tooltip-lr",
                    style={"cursor": "pointer", "color":"grey"},
                ),disabled=False, href="linear_regression")),

        dbc.NavItem(dbc.NavLink(
            html.Span(
                    "Differential Expression",
                    id="tooltip-de",
                    style={"cursor": "pointer", "color":"grey"},
                ),disabled=False, href="differential_expression")),

                dbc.NavItem(dbc.NavLink(
                    html.Span(
                            "Pathway Analysis",
                            id="tooltip-pa",
                            style={"cursor":"pointer", "color":"grey"},
                        ),disabled=False, href="pathway_analysis")),

        # tooltip for linear regression
        dbc.Tooltip(
        "Coming Soon!",
        target="tooltip-lr"
        ),

        # tooltip for differential expression
        dbc.Tooltip(
        "Coming Soon!",
        target="tooltip-de"
        ),

        # tooltip for pathway analysis
        dbc.Tooltip(
        "Coming Soon!",
        target="tooltip-pa"
        ),

        html.Hr(),
        control_panel
    ],
    vertical="md",
    pills=True
        ), md=2, className="mb-3"),

        #dbc.Col(control_panel, md=6)
        dbc.Col(first_card, md=6),
        dbc.Col(second_card, md=4)
        ],

        className="mb-3"),

    dbc.Row([dbc.Col(third_card, md=7, align="center"), dbc.Col(fourth_card, md=5, align="center")], className="mb-3")

], fluid=True)

dataset_dict = {
        "Proteins":"proteomics",
        "Lipids":"lipidomics",
        "Metabolites":"metabolomics",
        "Transcripts":"transcriptomics"
    }

df_dict = {
    "proteomics":proteomics_df,
    "lipidomics":lipidomics_df,
    "metabolomics":metabolomics_df,
}

quant_value_range_dict = {
    "proteomics":proteomics_quant_range,
    "lipidomics":lipidomics_quant_range,
    "metabolomics":metabolomics_quant_range,
}

@app.callback(
    dash.dependencies.Output('biomolecule_id', 'options'),
    [Input('dataset_id', 'value')])

def update_biomolecule_options(dataset_id):

    dataset = dataset_dict[dataset_id]
    df = df_dict[dataset]
    quant_value_range = quant_value_range_dict[dataset]

    # get list of columns for dataset
    available_biomolecules = df.columns[:quant_value_range].sort_values().tolist()

    options = [{'label': i, 'value': i} for i in available_biomolecules]

    return options

@app.callback(
    Output('biomolecule_id', 'value'),
    [Input('dataset_id', 'value')])

def update_default_biomolecule(dataset_id):

    dataset = dataset_dict[dataset_id]
    df = df_dict[dataset]
    quant_value_range = quant_value_range_dict[dataset]

    # get list of columns for dataset
    available_biomolecules = df.columns[:quant_value_range].sort_values().tolist()

    default_biomolecule = available_biomolecules[0]

    return default_biomolecule

@app.callback(
    Output('pca-scores-figure', 'figure'),
    [Input('dataset_id', 'value')])

def update_pca_scores_plot(dataset_id):

    dataset = dataset_dict[dataset_id]
    df = df_dict[dataset]
    quant_value_range = quant_value_range_dict[dataset]

    fig = pca_scores_plot(df, quant_value_range)

    return fig

@app.callback(
    Output('pca-loadings-figure', 'figure'),
    [Input('dataset_id', 'value')])

def update_pca_loadings_plot(dataset_id):

    dataset = dataset_dict[dataset_id]
    df = df_dict[dataset]
    quant_value_range = quant_value_range_dict[dataset]

    fig = pca_loadings_plot(df, quant_value_range, dataset_id)

    return fig

@app.callback(
    Output('biomolecule-barplot', 'figure'),
    [Input('biomolecule_id', 'value'),
    Input('dataset_id', 'value')])

def update_biomolecule_barplot(biomolecule_id, dataset_id):

    dataset = dataset_dict[dataset_id]
    df = df_dict[dataset]

    biomolecule_name = biomolecule_id
    x = df.index
    y = biomolecule_id
    fig = biomolecule_bar(df, x, y, biomolecule_name)

    return fig

@app.callback(
    Output('biomolecule-boxplot', 'figure'),
    [Input('biomolecule_id', 'value'),
    Input('dataset_id', 'value')])

def update_biomolecule_boxplot(biomolecule_id, dataset_id):

    dataset = dataset_dict[dataset_id]
    df = df_dict[dataset]

    biomolecule_name = biomolecule_id
    x = df.index
    y = biomolecule_id
    fig = boxplot(df, biomolecule_name)

    return fig

print("Starting server...")

if __name__ == '__main__':
    app.run_server(
        debug=True,
        host='0.0.0.0',
        #port='8080'
        )
