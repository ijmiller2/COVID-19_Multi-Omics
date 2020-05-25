
import pandas as pd
from tqdm import tqdm

proteomics_measurements_path = "/Users/Ian/Desktop/UW_2020/" + \
    "COVID19_multi-omics/COVID-19_Multi-Omics/data/proteomics/" + \
    "proteomics_measurements_wide.xlsx"
proteomics_LFQ_measurements_df = pd.read_excel(
    proteomics_measurements_path,
    sheet_name="LFQ",
    index_col="biomolecule_id")

proteomics_LFQ_measurements_df.stack().head()

proteomics_LFQ_measurements_stacked_df = proteomics_LFQ_measurements_df.stack()
proteomics_LFQ_measurements_stacked_df = proteomics_LFQ_measurements_stacked_df.reset_index()
proteomics_LFQ_measurements_stacked_df.columns = ['biomolecule', 'value_type', 'value']

# populate a dict of dict
# now iterate through columns and break into target columns
# [measurement_id	replicate_id	biomolecule_id	raw_LFQ_abundance
#normalized_LFQ_abundance	raw_iBAQ_abundance	normalized_iBAQ_abundance]

test_dict = {
    "7593": {
            "77": {
                "measurement_id":1,
                "raw_LFQ":237240000000,
                "normalized_LFQ":165150000000
            }
    }
}

measurement_dict = {}
measurement_id = 1

for index, row in tqdm(proteomics_LFQ_measurements_stacked_df.iterrows(),
    total=proteomics_LFQ_measurements_stacked_df.shape[0]):

    biomolecule_id = row['biomolecule']
    value_type = row['value_type']
    quant_value = row['value']
    # raw or normalized LFQ
    quant_value_type = "_".join(value_type.split("_")[0:2]) # "raw_LFQ_77" -> "raw_LFQ"

    replicate_id = value_type.split("_")[-1]

    # new biomolecule
    if not biomolecule_id in measurement_dict:

        measurement_dict[biomolecule_id] = {}
        measurement_dict[biomolecule_id][replicate_id] = {}

        measurement_dict[biomolecule_id][replicate_id]['measurement_id'] = measurement_id
        measurement_dict[biomolecule_id][replicate_id][quant_value_type] = quant_value

        measurement_id += 1

    # new replicate/measurement
    elif not replicate_id in measurement_dict[biomolecule_id]:

        measurement_dict[biomolecule_id][replicate_id] = {}
        measurement_dict[biomolecule_id][replicate_id]['measurement_id'] = measurement_id
        measurement_dict[biomolecule_id][replicate_id][quant_value_type] = quant_value

        measurement_id += 1

    # new quant value type, not a new measurement
    else:

        measurement_dict[biomolecule_id][replicate_id][quant_value_type] = quant_value


### now repeat with iBAQ measurements ###
proteomics_measurements_path = "/Users/Ian/Desktop/UW_2020/" + \
    "COVID19_multi-omics/COVID-19_Multi-Omics/data/proteomics/" + \
    "proteomics_measurements_wide.xlsx"
proteomics_iBAQ_measurements_df = pd.read_excel(
    proteomics_measurements_path,
    sheet_name="iBAQ",
    index_col="biomolecule_id")

proteomics_iBAQ_measurements_df.stack().head()

proteomics_iBAQ_measurements_stacked_df = proteomics_iBAQ_measurements_df.stack()
proteomics_iBAQ_measurements_stacked_df = proteomics_iBAQ_measurements_stacked_df.reset_index()
proteomics_iBAQ_measurements_stacked_df.columns = ['biomolecule', 'value_type', 'value']

for index, row in tqdm(proteomics_iBAQ_measurements_stacked_df.iterrows(),
    total=proteomics_iBAQ_measurements_stacked_df.shape[0]):

    biomolecule_id = row['biomolecule']
    value_type = row['value_type']
    quant_value = row['value']
    # raw or normalized LFQ
    quant_value_type = "_".join(value_type.split("_")[0:2]) # "raw_LFQ_77" -> "raw_LFQ"

    replicate_id = value_type.split("_")[-1]

    # should just be adding new quant values from iBAQ data
    # should not be any new measurements
    measurement_dict[biomolecule_id][replicate_id][quant_value_type] = quant_value

# write out data
outpath = "/Users/Ian/Desktop/UW_2020/" + \
    "COVID19_multi-omics/COVID-19_Multi-Omics/data/proteomics/" + \
    "proteomics_measurements.csv"

# table column names
outfile_col_name_list = ["measurement_id", "replicate_id", "biomolecule_id",
            "raw_LFQ", "normalized_LFQ", "raw_iBAQ", "normalized_iBAQ"]

with open(outpath, "w") as outfile:
    # write headers
    outfile.write(",".join(outfile_col_name_list) + "\n")

    for biomolecule_id in tqdm(measurement_dict, total=len(measurement_dict)):
        for replicate_id in measurement_dict[biomolecule_id]:

            measurement_id = measurement_dict[biomolecule_id][replicate_id]['measurement_id']
            raw_LFQ = measurement_dict[biomolecule_id][replicate_id]['raw_LFQ']
            normalized_LFQ = measurement_dict[biomolecule_id][replicate_id]['normalized_LFQ']
            raw_iBAQ = measurement_dict[biomolecule_id][replicate_id]['raw_iBAQ']
            normalized_iBAQ = measurement_dict[biomolecule_id][replicate_id]['normalized_iBAQ']

            # create ordered list of values based on colnames/keys
            outlist = [measurement_id, replicate_id, biomolecule_id,
                        raw_LFQ, normalized_LFQ, raw_iBAQ, normalized_iBAQ]

            #print("\t".join([str(i) for i in outlist]))
            outfile.write(",".join([str(i) for i in outlist]) + "\n")
