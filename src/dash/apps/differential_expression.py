
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
print("Loading data for differential_expression...")
print()
# load metabolomics data matrix
print("Loading metabolomics data...")
metabolomics_df, metabolomics_quant_range = get_omics_data(dataset='metabolomics', with_metadata=True)
print("Metabolomics data shape: {}".format(metabolomics_df.shape))
print("Loading lipidomics data...")
lipidomics_df, lipidomics_quant_range = get_omics_data(dataset='lipidomics', with_metadata=True)
print("Lipidomics data shape: {}".format(lipidomics_df.shape))
print("Loading proteomics data...")
proteomics_df, proteomics_quant_range = get_omics_data(dataset='proteomics', with_metadata=True)
print("Proteomics data shape: {}".format(proteomics_df.shape))

# make biomolecule_name_dict
metabolomics_biomolecule_names_dict = get_biomolecule_names(dataset='metabolomics')
lipidomics_biomolecule_names_dict = get_biomolecule_names(dataset='lipidomics')
proteomics_biomolecule_names_dict = get_biomolecule_names(dataset='proteomics')

# drop unknown lipids (to test speed up)
"""lipidomics_drop_list = []
for biomolecule_id in lipidomics_df.columns[:lipidomics_quant_range]:
    if "Unknown" in lipidomics_biomolecule_names_dict[biomolecule_id]:
        lipidomics_drop_list.append(biomolecule_id)
lipidomics_df.drop(lipidomics_drop_list, axis=1, inplace=True)
lipidomics_quant_range = lipidomics_quant_range - len(lipidomics_drop_list)"""

# define dataset dictionaries
dataset_dict = {
        "Proteins":"proteomics",
        "Lipids":"lipidomics",
        "Metabolites":"metabolomics",
        "Transcripts":"transcriptomics",
        "Combined":"combined"
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

global_names_dict = {
    "proteomics":proteomics_biomolecule_names_dict,
    "lipidomics":lipidomics_biomolecule_names_dict,
    "metabolomics":metabolomics_biomolecule_names_dict,
    "combined":{**proteomics_biomolecule_names_dict,
                **lipidomics_biomolecule_names_dict,
                **metabolomics_biomolecule_names_dict}
}

# get combined omics df and quant value range
print("Creating combined omics df...")
df_dict, quant_value_range_dict = get_combined_data(df_dict, quant_value_range_dict)

# start with proteomics data
sorted_biomolecule_names_dict = {k: v for k, v in sorted(global_names_dict['combined'].items(), key=lambda item: item[1])}
#available_biomolecules = proteomics_biomolecule_names_dict.values()
#available_biomolecules = proteomics_df.columns[:proteomics_quant_range].sort_values().tolist()
default_biomolecule = list(sorted_biomolecule_names_dict.keys())[0]

plotly_config = {"toImageButtonOptions":{'format':'svg',
                'filename': 'dash_plot'},
                "displaylogo": False}

dataset = 'combined'
combined_omics_df = df_dict[dataset]
quant_value_range = quant_value_range_dict[dataset]

print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("Getting pvalue data..")
pvalues_df = get_p_values()
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("Building volcano plot..")

confounders='ICU_1;Gender;Age_less_than_90'
volcano_df = get_volcano_data(pvalues_df[(pvalues_df['confounders']==confounders) & (pvalues_df['comparison']=='COVID_vs_NONCOVID')],
    df_dict,
    quant_value_range, global_names_dict)
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("Rendering volcano plot..")
fig = volcano_plot(volcano_df)
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

biomolecule_options = [{'label': value, 'value': key} for key, value in sorted_biomolecule_names_dict.items() if int(key) in pvalues_df['biomolecule_id'].to_list()]

control_panel = dbc.Card(
    [
        dbc.CardHeader("CONTROL PANEL",
                            style={"background-color":"#5bc0de",
                                        "font-weight":"bold",
                                        "font-size":"large"}),
        dbc.CardBody(
            [

            html.P("Biomolecule Data", className="card-title", style={"font-weight":"bold"}),

            # NOTE: This is dcc object not dbc
            dcc.Dropdown(
                id='biomolecule_id-de',
                # label maps to biomolecule name, value to biomolecule_id
                options=biomolecule_options,
                # only passing in quant value columns
                value=default_biomolecule,
                className="dropdown-item p-0"),

                ])
    ])

first_card = dbc.Card(
    [
        dbc.CardHeader("VOLCANO PLOT",
                            style={"background-color":"#5bc0de",
                                        "font-weight":"bold",
                                        "font-size":"large"}),
        dbc.CardBody(dcc.Graph(figure=fig, id='volcano',
        config=plotly_config))
    ])

# drop p_value_id
pvalues_columns = pvalues_df.columns[1:]
table_example = dash_table.DataTable(
    id='table',
    style_table={'overflowX': 'auto'},
    columns=[{"name": i, "id": i} for i in pvalues_columns],
    data=pvalues_df[pvalues_df['biomolecule_id']==int(default_biomolecule)].to_dict('records'),
)
second_card = dbc.Card(
    [
        dbc.CardHeader("BIOMOLECULE DATA",
                            style={"background-color":"#5bc0de",
                                        "font-weight":"bold",
                                        "font-size":"large"}),
        dbc.CardBody(table_example)
    ], className="mb-3")

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

        dbc.NavItem(dbc.NavLink("Differential Expression", active=True, href="differential_expression", style={"background-color":"grey"})),

        dbc.NavItem(dbc.NavLink(
            html.Span(
                    "Pathway Analysis",
                    id="tooltip-pa",
                    style={"cursor":"pointer", "color":"grey"},
                ),disabled=False, href="#")),

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

        dbc.Col(first_card, md=7),
        ],

        className="mb-3"),

    dbc.Row([dbc.Col(html.Div(), md=2, align="center"), dbc.Col(second_card, md=7, align="center")], className="mb-3")


], fluid=True, style={"height": "100vh"})




@app.callback(
    Output('table', 'data'),
    [Input('biomolecule_id-de', 'value')])
def update_table(biomolecule_id):

    data_list = pvalues_df[pvalues_df['biomolecule_id']==int(biomolecule_id)].to_dict('records')

    for data in data_list:
        for key, value in data.items():

            if key in ['p_value', 'q_value']:
                data[key] = '%.3E' % value

            elif key in ['log2_FC', 'neg_log10_p_value', 'effect_size']:
                data[key] = '%.3f' % value

    return data_list # list of dicts
