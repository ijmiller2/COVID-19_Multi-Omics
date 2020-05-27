
from sqlalchemy import create_engine, MetaData, Table, select, join
import pandas as pd

# SQLite path
db_path = 'sqlite:///../../data/SQLite Database/20200525/Covid-19 Study DB.sqlite'

omics_id_dict = {
        "proteomics":1,
        "lipidomics":2,
        "metabolomics":3,
        "transcriptomics":4
    }

def get_omics_data(with_metadata=False, dataset="proteomics"):

    omics_id = omics_id_dict[dataset]

    # Create an engine that connects to the Covid-19 Study DB.sqlite file: engine
    engine = create_engine(db_path)

    # Establish connection
    connection = engine.connect()

    # pull table into df
    query = "SELECT * from {}_measurements".format(dataset)
    omics_measurements_df = pd.read_sql_query(query, connection)

    # pull table into df
    query = "SELECT * from {}_runs".format(dataset)
    omics_runs_df = pd.read_sql_query(query, connection)

    # pull table into df
    query = "SELECT * from rawfiles WHERE ome_id={} AND sample_ID<>-1 and keep=1".format(omics_id)
    rawfiles_df = pd.read_sql_query(query, connection)

    # pull table into df
    deidentified_patient_metadata_df = pd.read_sql_query("SELECT * from deidentified_patient_metadata", connection)

    # make sure the merge by columns are all the same type -> pandas seems sensitive to this
    omics_measurements_df = omics_measurements_df.astype({'replicate_id': 'int32'})
    omics_runs_df = omics_runs_df.astype({'replicate_id': 'int32', 'rawfile_id': 'int32'})
    rawfiles_df = rawfiles_df.astype({'rawfile_id': 'int32', 'sample_id': 'int32'})
    deidentified_patient_metadata_df = deidentified_patient_metadata_df.astype({'sample_id': 'int32'})

    joined_df = omics_measurements_df\
                .join(omics_runs_df.set_index('replicate_id'), on='replicate_id')\
                .join(rawfiles_df.set_index('rawfile_id'), on='rawfile_id')\
                .join(deidentified_patient_metadata_df.set_index('sample_id'), on='sample_id')

    # drop samples that are missing COVID or ICU status
    joined_df.dropna(subset=['ICU_1','COVID'], inplace=True)

    # pivot to wide format
    wide_df = joined_df.pivot_table(index='sample_id', columns='biomolecule_id', values='normalized_abundance')
    wide_df.columns = [str(col) for col in wide_df.columns]

    query = "SELECT * from biomolecules WHERE omics_id={}".format(omics_id)
    # get biomolecule names
    biomolecules_df = pd.read_sql_query(query, connection)

    # close DB connection
    connection.close()

    # build biomolecule name dict and drop list
    biomolecule_name_dict = {}
    biomolecule_drop_list = []
    for index, row in biomolecules_df.iterrows():
        biomolecule_id = str(row['biomolecule_id'])
        standardized_name = row['standardized_name']
        biomolecule_name_dict[biomolecule_id] = standardized_name

        keep = row['keep']
        if keep!="1":
            biomolecule_drop_list.append(biomolecule_id)

    # drop biomolecules
    wide_df.drop(biomolecule_drop_list, axis=1, inplace=True)

    # replace wide_df column names
    #new_col_names = []
    #for col in wide_df.columns:
    #    new_col_names.append(biomolecule_name_dict[str(col)])
    #wide_df.columns = new_col_names

    # record quant value range
    quant_value_range = wide_df.shape[1]

    # optional return matrix with clinical metadata
    if with_metadata:

        combined_df = wide_df.join(deidentified_patient_metadata_df.set_index('sample_id'), on='sample_id')#.dropna()
        return combined_df, quant_value_range

    return wide_df, quant_value_range

def get_biomolecule_names(dataset='proteomics'):

    omics_id = omics_id_dict[dataset]

    # Create an engine that connects to the Covid-19 Study DB.sqlite file: engine
    engine = create_engine(db_path)

    # Establish connection
    connection = engine.connect()

    query = "SELECT * from biomolecules WHERE omics_id={} and KEEP=1".format(omics_id)
    # get biomolecule names
    biomolecules_df = pd.read_sql_query(query, connection)

    # build biomolecule name dict and drop list
    biomolecule_name_dict = {}
    for index, row in biomolecules_df.iterrows():
        biomolecule_id = str(row['biomolecule_id'])
        standardized_name = row['standardized_name']
        biomolecule_name_dict[biomolecule_id] = standardized_name

    # return dictionary with biomolecule ids and standard names

    if not dataset=="proteomics":

        # close DB connection
        connection.close()

        return biomolecule_name_dict

    # for proteomics data, return fasta headers instead

    query = "SELECT * from metadata"
    # get biomolecule names
    metadata_df = pd.read_sql_query(query, connection)

    ## NOTE: Could swap this with fasta headers once they're available
    gene_name_df = metadata_df[metadata_df['metadata_type'] == 'gene_name']
    gene_name_df = gene_name_df.astype({'biomolecule_id': 'str'})

    for biomolecule_id in biomolecule_name_dict:
        if biomolecule_id in gene_name_df['biomolecule_id'].tolist():
            # update to gene name
            gene_name = gene_name_df[gene_name_df['biomolecule_id']==biomolecule_id]['metadata_value'].values[0]
            if not gene_name == "":
                biomolecule_name_dict[biomolecule_id] = gene_name

    # close DB connection
    connection.close()

    return biomolecule_name_dict
