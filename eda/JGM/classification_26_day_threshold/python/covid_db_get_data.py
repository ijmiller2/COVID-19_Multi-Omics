from sqlalchemy import create_engine, MetaData, Table, select, join
import pandas as pd
import re


# SQLite path
db_path = 'sqlite:///C:\\covid_proteomics\\Covid-19 Study DB.sqlite'

def get_omics_data(with_metadata=False, dataset="proteomics"):

    omics_id_dict = {
        "proteomics":1,
        "lipidomics":2,
        "metabolomics":3,
        "transcriptomics":5
    }

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
    new_col_names = []
    for col in wide_df.columns:
        new_col_names.append(biomolecule_name_dict[str(col)])
    wide_df.columns = new_col_names

    # record quant value range
    quant_value_range = wide_df.shape[1]

    # optional return matrix with clinical metadata
    if with_metadata:

        combined_df = wide_df.join(deidentified_patient_metadata_df.set_index('sample_id'), on='sample_id')#.dropna()
        return combined_df, quant_value_range
    
    # close DB connection
    connection.close()
    
    return wide_df, quant_value_range

if __name__ == "__main__":

    dataset=args['dataset'].lower()
    if not dataset in ['proteomics', 'metabolomics', 'lipidomics', 'transcriptomics']:
        print("Dataset {} not recognized...".format(dataset))
        print("Please select from: {}".format(",".join(dataset)))

    omics_df, quant_value_range = get_omics_data(with_metadata=True, dataset=dataset)
    quant_columns = omics_df.columns[:quant_value_range]
    quant_df = omics_df[quant_columns]
    print("{} combined data set has shape: {}".format(dataset, omics_df.shape))
    print("{} quant data set has shape: {}".format(dataset, quant_df.shape))

    date = datetime.today().strftime('%Y-%m-%d')
    outpath = "{}_{}.tsv".format(date, dataset)
    print("Writing out data to: {}".format(outpath))
    omics_df.to_csv(outpath, sep="\t")
