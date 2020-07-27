
import dash
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import datetime
import pandas as pd

from data import get_omics_data, get_biomolecule_names, get_combined_data, get_p_values, get_volcano_data
from plot import volcano_plot

# importing app through index page
from app import app

print()
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("Loading data for heatmap...")
print()

# load metabolomics data matrix
print("Loading metabolomics data...")
from app import metabolomics_df, metabolomics_quant_range
print("Metabolomics data shape: {}".format(metabolomics_df.shape))
print("Loading lipidomics data...")
from app import lipidomics_df, lipidomics_quant_range
print("Lipidomics data shape: {}".format(lipidomics_df.shape))
print("Loading proteomics data...")
from app import proteomics_df, proteomics_quant_range
print("Proteomics data shape: {}".format(proteomics_df.shape))
print("Loading transcriptomics data...")
from app import transcriptomics_df, transcriptomics_quant_range
print("Transcriptomics data shape: {}".format(transcriptomics_df.shape))

available_datasets = ['Proteins', 'Lipids', 'Metabolites', 'Transcripts']

# define dataset dictionaries
# define dataset dictionaries
from app import dataset_dict, df_dict, quant_value_range_dict
from app import metabolomics_biomolecule_names_dict
from app import lipidomics_biomolecule_names_dict
from app import proteomics_biomolecule_names_dict
from app import transcriptomics_biomolecule_names_dict

global_names_dict = {
    "proteomics":proteomics_biomolecule_names_dict,
    "lipidomics":lipidomics_biomolecule_names_dict,
    "metabolomics":metabolomics_biomolecule_names_dict,
    "transcriptomics":transcriptomics_biomolecule_names_dict,
    "combined":{**proteomics_biomolecule_names_dict,
                **lipidomics_biomolecule_names_dict,
                **metabolomics_biomolecule_names_dict,
                **transcriptomics_biomolecule_names_dict}
}



control_panel = dbc.Card(
    [
        dbc.CardHeader("CONTROL PANEL",
                            style={"background-color":"#5bc0de",
                                        "font-weight":"bold",
                                        "font-size":"large"}),
        dbc.CardBody(
            [html.P("Select Dataset", className="card-title", style={"font-weight":"bold"}),
            dcc.Dropdown(
                id='dataset-hm',
                options=[{'label': i, 'value': i} for i in available_datasets],
                # only passing in quant value columns
                value=available_datasets[0]),

                ])
    ])

first_card = html.Iframe(
            id="cgrammer-iframe",
            src="https://amp.pharm.mssm.edu/clustergrammer/viz/5eed203e8ec9bb2d622075e9/proteomics.txt",
            height='600',
            width='900',
            style={"border-color":"transparent"}
        )

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

#app.layout = dbc.Container([
layout = dbc.Container([

    navbar,

    html.Hr(),

    dbc.Row(dbc.Col(html.H1("COVID-19 Multi-Omics Data Dashboard"), width={"size": 6, "offset": 3})),

    html.Hr(),

    dbc.Row(
        [dbc.Col(
        dbc.Nav(
    [
        html.H3("TYPE OF ANALYSIS", style={"font-weight":"bold", "color":"black"}),

        dbc.NavItem(dbc.NavLink(html.Span("PCA"),
            disabled=False,
            href="pca",
            style={"color":"grey"})),

        dbc.NavItem(dbc.NavLink(

            html.Span(
                    "Linear Regression",
                    id="tooltip-lr",
                    style={"cursor": "pointer", "color":"grey"},
                ),disabled=False, href="linear_regression")),
        
        dbc.NavItem(dbc.NavLink(

            html.Span(
                    "Differential Expression",
                    id="tooltip-lr",
                    style={"cursor": "pointer", "color":"grey"},
                ),disabled=False, href="differential_expression")),

        dbc.NavItem(dbc.NavLink("Clustergrammer", active=True, href="clustergrammer", style={"background-color":"grey"})),

        html.Hr(),
        control_panel
    ],
    vertical="md",
    pills=True
        ), md=2, className="mb-3"),

        dbc.Col(first_card, md=7),
        ],

        className="mb-3"),

], fluid=True, style={"height": "100vh"})


@app.callback(
    Output('cgrammer-iframe', 'src'),
    [Input('dataset-hm', 'value')])
def update_heatmap(dataset):

    url_dict = {
    "Proteins": "https://amp.pharm.mssm.edu/clustergrammer/viz/5eed203e8ec9bb2d622075e9/proteomics.txt",
    "Lipids": "https://amp.pharm.mssm.edu/clustergrammer/viz/5f061a208ec9bb6fb2f14a1d/lipidomics.txt",
    "Metabolites": "https://amp.pharm.mssm.edu/clustergrammer/viz/5f061c7f8ec9bb6fb2f14a37/metabolomics.txt",
    "Transcripts": "https://amp.pharm.mssm.edu/clustergrammer/viz/5f061cfc8ec9bb6fb2f14a45/transcriptomics.txt"

    }

    url = url_dict[dataset]

    return url

if __name__ == '__main__':

    import dash_bootstrap_components as dbc
    external_stylesheets=[dbc.themes.BOOTSTRAP]

    app = dash.Dash(
        __name__,
        external_stylesheets=external_stylesheets)
    app.title = 'Clustergrammer'

    app.layout = layout

    @app.callback(
        Output('cgrammer-iframe', 'src'),
        [Input('dataset-hm', 'value')])
    def update_heatmap(dataset):

        url_dict = {
        "Proteins": "https://amp.pharm.mssm.edu/clustergrammer/viz/5eed203e8ec9bb2d622075e9/proteomics.txt",
        "Lipids": "https://amp.pharm.mssm.edu/clustergrammer/viz/5f061a208ec9bb6fb2f14a1d/lipidomics.txt",
        "Metabolites": "https://amp.pharm.mssm.edu/clustergrammer/viz/5f061c7f8ec9bb6fb2f14a37/metabolomics.txt",
        "Transcripts": "https://amp.pharm.mssm.edu/clustergrammer/viz/5f061cfc8ec9bb6fb2f14a45/transcriptomics.txt"

        }

        url = url_dict[dataset]

        return url

    app.run_server(
        debug=True,
        host='0.0.0.0',
        #port='8080'
        )
