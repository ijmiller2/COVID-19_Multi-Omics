##### README ###### 

## "0_pathway_toolkit.R"
# description: This script contains 2 function useful for pathway 
#   analysis in R. Priciple is categorical terms are used to make 
#   a master list of id-to-category relationships. Then 2nd function
#   uses the mater list (AKA reference set) when performing enrichment
#   analysis usign fisher exact test, outputs the enrichemnt score/pvalue
#   and adjusted p-value for the categorical terms. This code was originally
#   produced for the dental informatics project. 
# issue: #9
# date created: 11/07/2017
# date last modified: 05/30/2020


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
# date last modified: 5/27/20
# input:
# - Covid-19 Study DB.sqlite

## "X4_KAO_updating_biomolecules_keep_column_GC_metabolites.R"
# description: modifies DB to update metabolite keep column to 
#   denote features which should be excluded from downstream analysis.
#   5/27/20 added more filter - tier information. 
# CAUTION: script iterates and caution should be used when executing. 
# date created: 5/18/20
# date last modified: 5/27/20 
# input: 
# - Covid-19 Study DB.sqlite
# output:
# - Covid-19 Study DB.sqlite (biomolecules 'keep' colulmn updated)


## "04_KAO_Exploring_GC_data_after_feature_filtering.R"
# description: looks at GC data by PCA after features have been
# filtered.
# date created: 5/18/2020
# date last modified: 5/19/2020
# input: 
# - Covid-19 Study DB.sqlite

## "05_KAO_Batch_effects_in_lipidomics_data.R"
# description: looks at batch effect of lipidomics data 
#   and in doing so, catches an initial error in the db 
#   entries for lipidomics features due to the way features 
#   were named resulting in duplicate identifiers. I will 
#   work with Dain to update the lipidomics values. 
# date created 5/19/2020
# date last modified: 5/20/2020
# input:
# - Covid-19 Study DB.sqlite
# - Lipidomics/Lipidomics_quant_results/Final_Results.csv


## X5_KAO_creating_new_lipidomics_table_to_match_original.R
# description: In file 05_KAO_Batch_effects_in_lipidomics_data.R,
#   I found that the biomolecule ids did not match up across 
#   the tables in the data frame. This code creates a csv that 
#   looks like the lipidomics_measurements table, but with 
#   updated biomolecule ids (no duplicates) and batch correction
#   to the lipiomics data - run-time correction similar to the
#   GC metabolomics data. Lipid standardized names are also updated
#   in this document. 
# issue: #7
# date created: 5/20/2020
# date last modified: 5/26/2020
# input:
# - Covid-19 Study DB.sqlite
# - Lipidomics/Lipidomics_quant_results/Final_Results.csv
# output:
# -  "../../data/lipidomics_measurements_20200523.csv"
# - above csv file was used to modify db

## "X6_KAO_Creating_pvalues_table.R"
# description: this script contains a function for likelyhood ratio testing between
#   two linear regression models. This script runs this function on all metabolomics,
#   and proteomics measurements. creates a p_value and then a q_value. These data
#   were added to the databse as pvalues table.
# issue: #4
# date created: 5/26/2020
# date last modified: 5/27/2020
# input: 
# - Covid-19 Study DB.sqlite
# output:
# - Covid-19 Study DB.sqlite, added table pvalues 

## "06_KAO_exploring_pvalue_histograms.R
# description: this is an exploratory data analysis of pvalues 
# generated by LR test. This looks at the overall effect 
# of confounders. 
# date created: 5/27/2020
# date last modified: 5/27/2020
# input: 
# - Covid-19 Study DB.sqlite

## "X7_KAO_updating_metadata_biomolecule_id.R"
# description: updates metadata table with non-duplicate biomolecule ids 
#   for lipidomics features. Also see file:  "X5_KAO_creating_new_lipidomics
#   _table_to_match_original.R" 
# date created: 5/30/2020
# date last modified: 5/30/2020
# input: 
# - Covid-19 Study DB.sqlite
# output:
# - Covid-19 Study DB.sqlite, modified metadata biomolecule_id 

