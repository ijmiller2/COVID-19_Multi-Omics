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
# output:
# - Covid-19 Study DB.sqlite (modified metabolite_measurements table)

## "X3_KAO_Updating_GC_metabolite_tier_in_DB.R"
# description: extracts the mean tier information by molecule. This tier
#   information is useful for filtering out poor quality metabolites and 
#   is added to the sqlite db metadata table
# data created: 5/15/20
# date last modified: 5/15/20
# input:
# - Covid-19 Study DB.sqlite
# output:
# - Covid-19 Study DB.sqlite (modified metadata table)

## "03_KAO_Exploring_GC_feature_quality.R"
# description: Explores 4 metrics of GC-metabolomics feature quality -
#   1) duplicate moleucles, 2) mean tier quality, 3) RSDs of QC sampels 
#   within and between batches, 4) dynamic range. 
# date created: 5/16/20
# date last modified: 5/18/20
# input:
# - Covid-19 Study DB.sqlite

