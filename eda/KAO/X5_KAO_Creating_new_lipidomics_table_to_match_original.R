#### X5_KAO_creating_new_lipidomics_table_to_match_original ###### 

## To fix the mismatching labels between the current 
## db lipidomics_measurements and the original 
## data frame. I am planning on creating a long 
## format of the original data frame that matches 
## the number and raw file indices from the current
## db. Then I will update the biomolecules_id in the 
## lipidomics_measurements table. 

library(DBI)
library(RSQLite)

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

##### Creating df of lipid raw abundance ###### 
lipids <- dbGetQuery(con, "SELECT unique_identifier, batch, timestamp, lipidomics_runs.replicate_id, biomolecule_id, raw_abundance, measurement_id
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           ")

lipid_ids <- dbGetQuery(con, "SELECT * 
                        FROM biomolecules
                        WHERE omics_id = 2")

dbDisconnect(con)

## load in original lipid file. 
lipids_original <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Lipidomics/Lipidomics_quant_results/Final_Results_with_unique_identifiers.csv", stringsAsFactors = F)

##### Checking file biomolecule order, ideally would not need to re-add biomolecules to table ### 

# confirming no duplicates in the naming this time
any(duplicated(lipids_original$Standardized_name)) #FALSE

# checking that the biomolecule order is the same (doing so by checking dimensions and RT)
lipids_original$Standardized_name [1:10] 
lipid_ids2 <- lipid_ids[-c(1:30), ] # as found in eda/KAO/05_KAO_Batch_effects_in_lipidiomics_data.R, removing first 30 duplicates

dim(lipids_original)
dim(lipid_ids2)

## check matches to RT 
biomolecule_id_rt <- sub(" -", "", sub(" +", "", sub("_Duplicate","", sub("Unknown Lipid RT", "", lipid_ids2$standardized_name), fixed = T), fixed = T), fixed = T)
as.numeric(biomolecule_id_rt) == lipids_original$Retention.Time..min.

## Adding current biomolecules.biomolecule_id information into what will become new input table
lipids_original <- cbind(lipid_ids2$biomolecule_id, lipids_original)

###### converting original to long format ######

lipids_original_long <- reshape(lipids_original, idvar = "standardized_name", 
                                ids = lipids_original$Standardized_name, 
                                times = names(lipids_original)[-c(1:9)], 
                                timevar = "unique_identifier",
                                varying = list(names(lipids_original)[-c(1:9)]),
                                direction = "long")

# Checking if the values match (mostly to) the current data table. 
lipids$raw_abundance[1:10] == lipids_original_long$X20200430_KAO_Control_1.raw..F1.[1:10]

# renaming raw_abundance measruement column (currently named 'X20200430_KAO_Control_1.raw..F1.')
names(lipids_original_long[grep("Control_1.raw", names(lipids_original_long))]) <- "raw_abundance"

# Checking that the file name identifiers match with current db table 
all(sapply(strsplit(sub("X", "", lipids_original_long$unique_identifier), ".raw"), function(x) x[1]) == lipids$unique_identifier)

# adding replicate_id to new table
lipids_original_long$replicate_id <- lipids$replicate_id 

dim(lipids) # 1324005 4 
dim(lipids_original_long) # 1324005 13

## lengths match, can use biomolecule_id in new table to replace current metabolomics_measurements.biomolecule_id in db

mismatch <- lipids_original_long$`lipid_ids2$biomolecule_id` != lipids$biomolecule_id
table(mismatch)

lipids_original_long$measurement_id <- lipids$measurement_id

#### Considering batch correction prior to re-uploading table #####
#Performing same steps as GC-metabolomics run-time correction

# converting timestamp to POSIXct class
lipids$timestamp <- as.POSIXct(strptime(lipids$timestamp, "%Y%m%d%H%M%S"))

# calculating difftime for linear regression model
lipids$diff_time <- as.numeric(difftime(lipids$timestamp, lipids$timestamp[1]))

# creating columne where the log2 raw abudance values are in non-log 
lipids$raw_abundance_nonlog <- lipids$raw_abundance

# creating an NA column that will eventually hold the run-time corrected values
lipids$runtime_corrected <- NA

# creating a weights column that more heavily weights control values in the rlm model
lipids$weights <- 0
lipids$weights[grepl("Control", lipids$unique_identifier)] <- 4
lipids$weights[grepl("half|quarter", lipids$unique_identifier)] <- 0


lipids$biomolecule_id <- lipids_original_long$`lipid_ids2$biomolecule_id`

#### For loop for calculating runtime corrected values #####
set.seed(1234)

for(i in unique(lipids$biomolecule_id)){
  #for(k in unique(lipids$batch)){
    k = 4 # redoing batch 4 only, realizing quarter/half samples were weighted. (seed was set 1234)
    index <- lipids$biomolecule_id == i & lipids$batch == k
    
    fit <- tryCatch(rlm(raw_abundance_nonlog~diff_time, data = lipids[index,], weights = weights),warning = function(w) fit =NULL, error = function(e) fit = NULL) 
    if(!is.null(fit) ){
      intercept <- summary(fit)$coef[1,1]
      slope <- summary(fit)$coef[2,1]
    } else {
      intercept <- mean(lipids$raw_abundance_nonlog[index], na.rm =T)
      slope <- 0
    }
    # normalize by run-order
    # sampleValue - ((slope * runOrder) + intercept) + biomolecule mean
    new_value <- (lipids$raw_abundance_nonlog[index] - ((slope * lipids$diff_time[index]) + intercept) )+ mean(lipids$raw_abundance_nonlog[lipids$biomolecule_id == i & grepl("Control", lipids$unique_identifier)], na.rm =T)
    # in cases where run-order corrected value is below 0, replacing with pre-corrected values
    # this step was decided after carefully looking at instances of negative values which were very low pre correction. 
    # in total only ~ 1000 values have this occurance 
    new_value[new_value<0] <- lipids$raw_abundance_nonlog[index][new_value<0]
    lipids$runtime_corrected[index] <-log2(new_value)
    
  }
#}

## save the table 
#write.csv(lipids, file = "../../data/lipids_run_time_corrected_long_format.csv", row.names = F)





#### Exploratory code to check run time corrections ##### 
# summary - looks reasonable. 
boxplot(log2(lipids$raw_abundance)~lipids$batch)
boxplot(lipids$runtime_corrected~lipids$batch)
plot(log2(lipids$raw_abundance), lipids$runtime_corrected)

hist(log2(lipids$raw_abundance))
hist(lipids$runtime_corrected)

## convert to wide 
post <- reshape(lipids, timevar = "biomolecule_id", v.names = "runtime_corrected",
               idvar = "unique_identifier", direction = "wide" )

# selecting only data which is either controls or samples
post <- post[!(grepl("resolvin|Method", post$unique_identifier)|grepl("half|quarter", post$unique_identifier)), ]

post[!grepl("Control", post$unique_identifier),]$unique_identifier

names(post)
# performing principle component analysis on normalized data 
pca_post <- prcomp(post[,-c(1:9)])
plot(pca_post$x, col = post$batch, pch = 19)

pre <- reshape(lipids, timevar = "biomolecule_id", v.names = "raw_abundance",
                idvar = "unique_identifier", direction = "wide" )

# selecting only data which is either controls or samples
pre <- pre[!(grepl("resolvin|Method", pre$unique_identifier)|grepl("half|quarter", pre$unique_identifier)), ]

pre[!grepl("Control", pre$unique_identifier),]$unique_identifier

names(pre)
# performing principle component analysis on normalized data 
pca_pre <- prcomp(log2(pre[,-c(1:9)]))
plot(pca_pre$x, col = pre$batch, pch = 19, main = "Raw_abundance")

par(mfrow = c(1,2))
plot(pca_pre$x, col = pre$batch, pch = 19, main = "Raw_abundance")
plot(pca_post$x, col = post$batch, pch = 19, main = "Run-time corrected")

####### Creating data frame for upload into db ###### 
## need measurement_id replicate_id biomolecule_id raw_abundance normalized_abundance 

df_upload <- data.frame(measurement_id = lipids$measurement_id, replicate_id = lipids$replicate_id, 
                        biomolecule_id = lipids_original_long$`lipid_ids2$biomolecule_id`, 
                        raw_abundance = lipids$raw_abundance, normalized_abundance = lipids$runtime_corrected, stringsAsFactors = F)


dim(df_upload) #1324005 5
dim(lipids) #1324005 11
dim(lipids_original_long) #1324005 14

names(df_upload)

write.csv(df_upload, file = "../../data/lipidomics_measurements_20200523.csv", row.names = F)

###### The above table was uploaded to the DB to directly replace the lipidomics measurements.

###### Creating a name update function for lipid biomolecule names ###### 
df <- data.frame(biomolecule_id = lipids_original$`lipid_ids2$biomolecule_id`, standardized_name = lipids_original$Standardized_name, stringsAsFactors = F)


##### R Function to replace values ##### 
replaceValues <- function(x) dbExecute(con, paste("UPDATE biomolecules
          SET standardized_name = '",unlist(x[2]),"' WHERE biomolecule_id = ", unlist(x[1]), sep = ""))


#replaceValues <- function(x) paste("UPDATE biomolecules
#          SET standardized_name = '",unlist(x[2]),"' WHERE biomolecule_id = ", unlist(x[1]), sep = "")

##### Iterate over data frame to update normalized_abundance values in sqlite db ######
dbReadTable(con, "biomolecules")

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

apply(df, 1, replaceValues)

dbDisconnect(con)

