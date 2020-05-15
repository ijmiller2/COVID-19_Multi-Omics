##### README ###### 

## "01_KAO_Establishing_connection_to_db_extracting_timeStamp.R"
# description: Establishes DB connection using RSQLite package and 
#     fetches time stamp information for Raw files.
# Relevant Issue(s): 
# date created: 5/12/20
# date last modified: 5/12/20
# input:
#  - Covid-19 Study DB.sqlite
# output:
#  - 

## "X1_KAO_Updating_GC_keep_status_in_db.R"
# description: Sets 1:10 split GC files keep column to 0 (FALSE)
# date created: 5/13/20
# date last modified: 5/13/20
# input: 
# - Covid-19 Study DB.sqlite
# output: 
# - Covid-19 Study DB.sqlite (modified)

## "02_KAO_Runorder_correction_for_GC_metabolomics_data.R"
# description: Performs run order correction of the GC data and explores results.
#   this analysis is exploratory only and does not modify the db.
# date created: 5/13/20
# date last modified: 5/14/20
# input:
# - Covid-19 Study DB.sqlite

## "X2_KAO_GC_metabolomics_runtime_correction.R"
# description: Performs run order correction and modifies db. 
# date created: 5/14/20
# date last modified: 5/15/20
# input:
# - Covid-19 Study DB.sqlite
# ouptu:
# - Covid-19 Study DB.sqlite (modified metabolite_measurements table)


