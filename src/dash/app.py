import dash

import dash_bootstrap_components as dbc
external_stylesheets=[dbc.themes.BOOTSTRAP]

app = dash.Dash(__name__,
    suppress_callback_exceptions=True,
    external_stylesheets=external_stylesheets)

app.title = 'COVID-19 Multi-Omics'
server = app.server

### load shared data ###
from data import get_omics_data, get_biomolecule_names
import datetime

print()
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("Loading data for app...")
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
print("Loading transcriptomics data...")
transcriptomics_df, transcriptomics_quant_range = get_omics_data(dataset='transcriptomics', with_metadata=True)
print("Transcriptomics data shape: {}".format(transcriptomics_df.shape))

# make biomolecule_name_dict
metabolomics_biomolecule_names_dict = get_biomolecule_names(dataset='metabolomics')
lipidomics_biomolecule_names_dict = get_biomolecule_names(dataset='lipidomics')
proteomics_biomolecule_names_dict = get_biomolecule_names(dataset='proteomics')
transcriptomics_biomolecule_names_dict = get_biomolecule_names(dataset='transcriptomics')

# define dataset dictionaries
dataset_dict = {
        "Proteins":"proteomics",
        "Lipids":"lipidomics",
        "Metabolites":"metabolomics",
        "Transcripts":"transcriptomics",
        "Combined Biomolecules":"combined"
    }

df_dict = {
    "proteomics":proteomics_df,
    "lipidomics":lipidomics_df,
    "metabolomics":metabolomics_df,
    "transcriptomics":transcriptomics_df
}

quant_value_range_dict = {
    "proteomics":proteomics_quant_range,
    "lipidomics":lipidomics_quant_range,
    "metabolomics":metabolomics_quant_range,
    "transcriptomics":transcriptomics_quant_range
}

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
