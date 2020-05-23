
from sqlalchemy import create_engine, MetaData, Table, select, join
import pandas as pd

# SQLite path
db_path = 'sqlite:///../../data/SQLite Database/Covid-19 Study DB.sqlite'

def get_metabolomics_data(with_metadata=False):
    # Create an engine that connects to the Covid-19 Study DB.sqlite file: engine
    engine = create_engine(db_path)

    # Establish connection
    connection = engine.connect()

    # pull table into df
    metabolomics_measurements_df = pd.read_sql_query("SELECT * from metabolomics_measurements", connection)

    # pull table into df
    metabolomics_runs_df = pd.read_sql_query("SELECT * from metabolomics_runs", connection)

    # pull table into df
    rawfiles_df = pd.read_sql_query("SELECT * from rawfiles WHERE ome_id=3 AND sample_ID<>-1 and keep=1", connection)

    # pull table into df
    deidentified_patient_metadata_df = pd.read_sql_query("SELECT * from deidentified_patient_metadata", connection)

    # make sure the merge by columns are all the same type -> pandas seems sensitive to this
    metabolomics_measurements_df = metabolomics_measurements_df.astype({'replicate_id': 'int32'})
    metabolomics_runs_df = metabolomics_runs_df.astype({'replicate_id': 'int32', 'rawfile_id': 'int32'})
    rawfiles_df = rawfiles_df.astype({'rawfile_id': 'int32', 'sample_id': 'int32'})
    deidentified_patient_metadata_df = deidentified_patient_metadata_df.astype({'sample_id': 'int32'})

    joined_df = metabolomics_measurements_df\
                .join(metabolomics_runs_df.set_index('replicate_id'), on='replicate_id')\
                .join(rawfiles_df.set_index('rawfile_id'), on='rawfile_id')\
                .join(deidentified_patient_metadata_df.set_index('sample_id'), on='sample_id')

    # drop samples that are missing COVID or ICU status
    joined_df.dropna(subset=['ICU_1','COVID'], inplace=True)

    # pivot to wide format
    wide_df = joined_df.pivot_table(index='sample_id', columns='biomolecule_id', values='normalized_abundance')

    # get biomolecule names
    biomolecules_df = pd.read_sql_query("SELECT * from biomolecules", connection)

    # build biomolecule name dict and drop list
    biomolecule_name_dict = {}
    biomolecule_drop_list = []
    for index, row in biomolecules_df.iterrows():
        biomolecule_id = row['biomolecule_id']
        standardized_name = row['standardized_name']
        biomolecule_name_dict[biomolecule_id] = standardized_name

        keep = row['keep']
        if keep!="1":
            biomolecule_drop_list.append(biomolecule_id)

    # drop biomolecules
    wide_df.drop(biomolecule_drop_list, axis=1, inplace=True)

    # replace wide_df column names
    new_col_names = []
    for col in wide_df.columns:
        new_col_names.append(biomolecule_name_dict[col])
    wide_df.columns = new_col_names

    # optional return matrix with clinical metadata
    if with_metadata:

        combined_df = wide_df.join(deidentified_patient_metadata_df.set_index('sample_id'), on='sample_id')#.dropna()
        return combined_df

    return wide_df
