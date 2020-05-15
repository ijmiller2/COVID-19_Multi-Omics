##### MODIFYING DB 2020-05-15 #####

## edit by KAO
## A common occurance in daily GC runs is signal loss or gain 
## throughout the day's analysis. We assume this is due in part
## to changing derivitization over time. In many cases it is minor
## but for a handful of metabolite features, run-time correction 
## can make a significant difference in data reproducibility.
##
## The file "eda/KAO/02_KAO_Runorder_correction_for_GC_metabolomics_data.R"
## explores the use of robust linear regression to perform run order 
## correction of the 'raw_abundance' GC data; in that file I assess metrics
## like QC control sample RSDs across batches before and after correction.
## This correction improves QC RSDs and minimally effects the 
## overall data distributions.
## 
## Below is code to perform run time correction and update the db with
## the normalized values 


library(DBI)
library(RSQLite)
library(MASS)

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

##### Create local data frame with run_time information and metabolite values #####

timeStamp_df<- dbGetQuery(con, "SELECT unique_identifier, measurement_id, timestamp, raw_abundance, biomolecule_id, batch
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           WHERE keep = 1
           ")

dbDisconnect(con)

#### Adding variables to data frame #####
# converting timestamp to POSIXct class
timeStamp_df$timestamp <- as.POSIXct(strptime(timeStamp_df$timestamp, "%Y%m%d%H%M%S"))

# calculating difftime for linear regression model
timeStamp_df$diff_time <- as.numeric(difftime(timeStamp_df$timestamp, timeStamp_df$timestamp[1]))

# creating columne where the log2 raw abudance values are in non-log 
timeStamp_df$raw_abundance_nonlog <- 2^timeStamp_df$raw_abundance

# creating an NA column that will eventually hold the run-time corrected values
timeStamp_df$runtime_corrected <- NA

# creating a weights column that more heavily weights control values in the rlm model
timeStamp_df$weights <- 0
timeStamp_df$weights[grepl("Control", timeStamp_df$unique_identifier)] <- 4

#### For loop for calculating runtime corrected values #####
set.seed(1234)

for(i in unique(timeStamp_df$biomolecule_id)){
  for(k in unique(timeStamp_df$batch)){
    index <- timeStamp_df$biomolecule_id == i & timeStamp_df$batch == k

    fit <- tryCatch(rlm(raw_abundance_nonlog~diff_time, data = timeStamp_df[index,], weights = weights),warning = function(w) fit =NULL, error = function(e) fit = NULL) 
    if(!is.null(fit) ){
      intercept <- summary(fit)$coef[1,1]
      slope <- summary(fit)$coef[2,1]
    } else {
      intercept <- mean(timeStamp_df$raw_abundance_nonlog[index], na.rm =T)
      slope <- 0
    }
    # normalize by run-order
    # sampleValue - ((slope * runOrder) + intercept) + biomolecule mean
    new_value <- (timeStamp_df$raw_abundance_nonlog[index] - ((slope * timeStamp_df$diff_time[index]) + intercept) )+ mean(timeStamp_df$raw_abundance_nonlog[timeStamp_df$biomolecule_id == i & grepl("Control", timeStamp_df$unique_identifier)], na.rm =T)
    # in cases where run-order corrected value is below 0, replacing with pre-corrected values
    # this step was decided after carefully looking at instances of negative values which were very low pre correction. 
    # in total only ~ 1000 values have this occurance 
    new_value[new_value<0] <- timeStamp_df$raw_abundance_nonlog[index][new_value<0]
    timeStamp_df$runtime_corrected[index] <-log2(new_value)
  
  }
}
###### Simplify df: measurement ID and normalized values ####
df <- data.frame(measurement_id = timeStamp_df$measurement_id, runtime_corrected = timeStamp_df$runtime_corrected)

###### Establish connection with DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

##### Test query metabolite values #####

dbGetQuery(con, "SELECT normalized_abundance
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           WHERE keep = 1
           ")

##### R Function to replace values ##### 
replaceValues <- function(x) dbExecute(con, paste("UPDATE metabolomics_measurements
          SET normalized_abundance = ",unlist(x[2])," WHERE measurement_id = ", unlist(x[1]), sep = ""))

##### Iterate over data frame to update normalized_abundance values in sqlite db ######

apply(df, 1, replaceValues)

##### Check that values match ##### 

normalized_abundance <- dbGetQuery(con, "SELECT normalized_abundance
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           WHERE keep = 1
           ")

table(normalized_abundance == timeStamp_df$runtime_corrected) ## not sure why false
summary(normalized_abundance[,1]) 
summary(timeStamp_df$runtime_corrected) 
# these numbers match

dbDisconnect(con)
