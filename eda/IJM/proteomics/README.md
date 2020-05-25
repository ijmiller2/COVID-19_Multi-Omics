# Notebooks and EDA for proteomics datasets

1.0_get_rawfile_timestamps.py

description: Generate rawfiles table for proteomics data
Relevant Issue(s): NA
date created: 5/23/20
date last modified: 5/23/20
input: "All_20200428_COVID_plasma_multiomics/Proteomics/COVID19 study single-shot DDA/*/*.raw"
output:
  - data/proteomics_raw_file_timestamps.tsv
  - [proteomics_raw_file_timestamps.xlsx](https://docs.google.com/spreadsheets/d/1GftO-cTfZURquhVuCGQBXe3QcwwABnYW/edit#gid=1832301628)

2.0_stack_proteomics_measurements.py

description: Process and reshape data for proteomics_measurements table upload.
Relevant Issue(s): NA
date created: 5/24/20
date last modified: 5/24/20
input: 
  - data/proteomics/proteomics_measurements_wide.xlsx
output:
  - data/proteomics/proteomics_raw_file_timestamps.tsv
  - data/proteomics/proteomics_measurements.csv