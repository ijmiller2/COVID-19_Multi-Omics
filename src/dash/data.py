
from sqlalchemy import create_engine, MetaData, Table, select, join
import pandas as pd
import re
import numpy as np

# SQLite path
db_path = 'sqlite:///../../data/SQLite Database/20200617/Covid-19 Study DB.sqlite'

omics_id_dict = {
        "proteomics":1,
        "lipidomics":2,
        "metabolomics":3,
        "transcriptomics":5
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
    if dataset == "metabolomics":
        query = "SELECT * from rawfiles WHERE ome_id=3 OR ome_id=4 AND sample_ID<>-1 AND keep=1".format(omics_id)

    else:
        query = "SELECT * from rawfiles WHERE ome_id={} AND sample_ID<>-1 and keep=1".format(omics_id)

    rawfiles_df = pd.read_sql_query(query, connection)
    ## NOTE: For some reason, SQL filter not working for keep on raw files
    rawfiles_df = rawfiles_df[rawfiles_df['keep']==1]

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

    if dataset == "metabolomics":
        query = "SELECT * from biomolecules WHERE omics_id=3 OR omics_id=4".format(omics_id)
    else:
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

    dataset_abr_prefix = "[{}] ".format(dataset[0].upper())

    print("Getting biomolecule names for dataset: {}".format(dataset))
    omics_id = omics_id_dict[dataset]

    # Create an engine that connects to the Covid-19 Study DB.sqlite file: engine
    engine = create_engine(db_path)

    # Establish connection
    connection = engine.connect()

    if dataset == "metabolomics":
        query = "SELECT * from biomolecules WHERE omics_id=3 OR omics_id=4 and KEEP=1".format(omics_id)
    else:
        query = "SELECT * from biomolecules WHERE omics_id={} and KEEP=1".format(omics_id)

    # get biomolecule names
    biomolecules_df = pd.read_sql_query(query, connection)

    # build biomolecule name dict and drop list
    biomolecule_name_dict = {}
    for index, row in biomolecules_df.iterrows():
        biomolecule_id = str(row['biomolecule_id'])
        standardized_name = row['standardized_name']
        biomolecule_name_dict[biomolecule_id] = dataset_abr_prefix + standardized_name

    # return dictionary with biomolecule ids and standard names

    if not dataset=="proteomics":

        # close DB connection
        connection.close()

        return biomolecule_name_dict

    # for proteomics data, return fasta headers instead

    query = "SELECT * from metadata"
    # get biomolecule names
    metadata_df = pd.read_sql_query(query, connection)

    fasta_header_df = metadata_df[metadata_df['metadata_type'] == 'fasta_header']
    fasta_header_dict = {}
    for index, row in fasta_header_df.iterrows():
        biomolecule_id = str(int(row['biomolecule_id'])) # string was getting truncated by rstrip(".0")
        fasta_header = str(row['metadata_value'])
        fasta_header_dict[biomolecule_id] = fasta_header

    for biomolecule_id in biomolecule_name_dict:

        fasta_header = fasta_header_dict[biomolecule_id]
        fasta_header = re.search("\s(.*?)\sO[SX]=", fasta_header).group(1)
        biomolecule_name_dict[biomolecule_id] = dataset_abr_prefix + fasta_header

    # close DB connection
    connection.close()

    return biomolecule_name_dict

def get_combined_data(df_dict, quant_range_dict, with_transcripts=False):

    # load metabolomics data matrix
    metabolomics_df, metabolomics_quant_range = df_dict['metabolomics'], quant_range_dict['metabolomics']
    lipidomics_df, lipidomics_quant_range = df_dict['lipidomics'], quant_range_dict['lipidomics']
    proteomics_df, proteomics_quant_range = df_dict['proteomics'], quant_range_dict['proteomics']
    if with_transcripts:
        transcriptomics_df, transcriptomics_quant_range = get_omics_data(dataset='transcriptomics', with_metadata=True)

    # get quant columns
    lipidomics_quant_columns = lipidomics_df.columns[:lipidomics_quant_range]
    lipidomics_quant_df = lipidomics_df[lipidomics_quant_columns]

    metabolomics_quant_columns = metabolomics_df.columns[:metabolomics_quant_range]
    metabolomics_quant_df = metabolomics_df[metabolomics_quant_columns]

    proteomics_quant_columns = proteomics_df.columns[:proteomics_quant_range]
    proteomics_quant_df = proteomics_df[proteomics_quant_columns]

    if with_transcripts:
        transcriptomics_quant_columns = transcriptomics_df.columns[:transcriptomics_quant_range]
        transcriptomics_quant_df = transcriptomics_df[transcriptomics_quant_columns]

    # get clinical_metadata_df
    clinical_metadata_columns = proteomics_df.columns[proteomics_quant_range:]
    clinical_metadata_df = proteomics_df[clinical_metadata_columns]
    clinical_metadata_df

    # join quant values together
    if with_transcripts:
        combined_df = proteomics_quant_df.join(lipidomics_quant_df).join(metabolomics_quant_df).join(transcriptomics_quant_df)
    else:
        combined_df = proteomics_quant_df.join(lipidomics_quant_df).join(metabolomics_quant_df)
    combined_quant_range = combined_df.shape[1]
    combined_quant_columns = combined_df.columns[:combined_quant_range]
    # now join with clinical metadata
    combined_df = combined_df.join(clinical_metadata_df)

    # drop any samples with missing values in quant columns
    combined_df.dropna(subset=combined_quant_columns,inplace=True)

    # also return df_dict with combined_df
    df_dict['combined'] = combined_df

    # update quant_range_dict
    quant_range_dict['combined'] = combined_quant_range

    return df_dict, quant_range_dict

def get_p_values():

    # Create an engine that connects to the Covid-19 Study DB.sqlite file: engine
    engine = create_engine(db_path)

    # Establish connection
    connection = engine.connect()

    query = "SELECT * from biomolecules WHERE KEEP=1"
    # get biomolecule names
    biomolecules_df = pd.read_sql_query(query, connection)

    # build biomolecule name dict and drop list
    biomolecule_name_dict = {}
    for index, row in biomolecules_df.iterrows():
        biomolecule_id = str(row['biomolecule_id'])
        standardized_name = row['standardized_name']
        biomolecule_name_dict[biomolecule_id] = standardized_name

    query = "SELECT * from pvalues"
    # get biomolecule names
    pvalues_df = pd.read_sql_query(query, connection)

    return pvalues_df

def get_volcano_data(pvalues_df, df_dict, quant_value_range,
    global_names_dict, comparison_column='COVID'):

    group_1_quant_value_dict = {}
    formula = 'normalized_abundance ~ COVID * ICU_1 + Gender + Age_less_than_90 vs. normalized_abundance ~ ICU_1 + Gender + Age_less_than_90'

    pvalues_df = pvalues_df[(pvalues_df['comparison']=='COVID_vs_NONCOVID') & (pvalues_df['formula']==formula)]
    print("Volcano data pvalues shape: {}".format(pvalues_df.shape))

    group_1 = 1
    group_2 = 0

    combined_df = df_dict['combined']
    quant_value_columns = combined_df.columns[:quant_value_range]

    group_1_quant_value_dict = {}

    for sample_id, row in combined_df[combined_df[comparison_column]==group_1].iterrows():

        for biomolecule_id in quant_value_columns:
            quant_value = row[biomolecule_id]

            if not biomolecule_id in group_1_quant_value_dict:
                group_1_quant_value_dict[biomolecule_id] = [quant_value]

            else:
                group_1_quant_value_dict[biomolecule_id].append(quant_value)

    group_2_quant_value_dict = {}

    for sample_id, row in combined_df[combined_df[comparison_column]==group_2].iterrows():

        for biomolecule_id in quant_value_columns:
            quant_value = row[biomolecule_id]

            if not biomolecule_id in group_2_quant_value_dict:
                group_2_quant_value_dict[biomolecule_id] = [quant_value]

            else:
                group_2_quant_value_dict[biomolecule_id].append(quant_value)


    FC_dict = {}
    for biomolecule_id in quant_value_columns:

        group_1_quant_values = group_1_quant_value_dict[biomolecule_id]
        group_2_quant_values = group_2_quant_value_dict[biomolecule_id]

        # in log2 space; subtract
        FC = np.mean(group_1_quant_values) - np.mean(group_2_quant_values)

        FC_dict[biomolecule_id] = FC

    FC_list = []
    ome_list = []
    standardized_name_list = []
    for index, row in pvalues_df.iterrows():
        biomolecule_id = str(row['biomolecule_id'])

        ## NOTE: should create biomolecule_ome_dict
        if biomolecule_id in df_dict['proteomics'].columns:
            ome_list.append("proteomics")
        elif biomolecule_id in df_dict['lipidomics'].columns:
            ome_list.append("lipidomics")
        elif biomolecule_id in df_dict['metabolomics'].columns:
            ome_list.append("metabolomics")
        else:
            #print("Biomolecule {} not mapped to ome!".format(biomolecule_id)) # don't currenly have targeted metabolomics included
            ome_list.append(np.nan)
            #break

        try:
            FC = FC_dict[biomolecule_id]
        except:
            # may have been a dropped biomolecule or from a different ome
            FC = np.nan

        FC_list.append(FC)

        try:
            standardized_name = global_names_dict['combined'][biomolecule_id]
        except:
            # may have been a dropped biomolecule or from a different ome
            standardized_name = np.nan

        standardized_name_list.append(standardized_name)

    pvalues_df['log2_FC'] = FC_list
    pvalues_df['ome_type'] = ome_list
    pvalues_df['standardized_name'] = standardized_name_list
    pvalues_df['neg_log10_p_value'] = pvalues_df['p_value'].apply(np.log10).apply(np.negative)

    return pvalues_df
