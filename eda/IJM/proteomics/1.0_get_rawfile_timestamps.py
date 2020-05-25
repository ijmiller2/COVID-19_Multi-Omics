
from os.path import getmtime, getctime
import glob
import pandas as pd

# wd : /Volumes/projects/All_20200428_COVID_plasma_multiomics/Proteomics/COVID19 study single-shot DDA/
proteomics_raw_file_paths = glob.glob("*/*.raw")

rawfile_id = 629
['rawfile_id',  'timestamp', 'rawfile_name',  'sample_id', 'run_type',  'keep',  batch ome_id]
for path in proteomics_raw_file_paths:

    # skip human controls
    if not "_HC_" in path:

        timestamp = str(getctime(path)).rstrip(".0")
        datetime = pd.to_datetime(timestamp, unit='s')
        formatted_datetime = datetime.strftime('%Y%m%d%H%S')
        sample_id = path.split("_")[-1].rstrip(".raw")
        batch = path.split("/")[0].split()[1]
        unique_identifier = path.split("/")[-1].rstrip(".raw")
        keep = 1
        ome_id = 1
        run_type = ""
        outlist = [rawfile_id, formatted_datetime, unique_identifier, sample_id, run_type, keep, batch, ome_id]

        print("\t".join([str(i) for i in outlist]))

    rawfile_id += 1
