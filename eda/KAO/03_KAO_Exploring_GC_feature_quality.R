##### 03 Consider features to keep in GC metabolomics data #####

## Some features in the current GC metabolomics data set should be
## removed due to poor reproducibility or duplicate identification
## (i.e., metabolites with multiple derivitization products). This 
## script aims to explore different aspects of feature quality. 
##
## The duplicates have been desiganated in the original data table
## as remove = T; this information is contained in the metadata table
## in the sqlite db. Tier information for each feature is another 
## important consideration. Tier 4 and 5 indicate the feature was 
## not found in a particular file and was imputed with a value 
## near the noise. Features with many tier 4 and 5 values should be 
## considered for removal. Lastly, QC samples within and between 
## batchs help determine reproducibility of a quantitative value. 
## RSDs of QC samples will also be used for filtering out a poor
## quality feature. Overall dynamic range might be a final consideration:
## unknown features which have limited dynamic range might be
## indicative of a feature which is an artifact and not biological
## in nature. 
## 

library(DBI)
library(RSQLite)
colors <- read.csv("../../reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

dbListTables(con)
dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "metadata")[1:10,]
dbReadTable(con, "metabolomics_measurements")[0,]
dbReadTable(con, "metabolomics_runs")[0,]
dbReadTable(con, "biomolecules")[0,]


##### 
df<- dbGetQuery(con, "SELECT unique_identifier, normalized_abundance, metabolomics_measurements.biomolecule_id, batch 
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN metadata on metadata.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE keep = 1 
           AND ome_id = 3 
           ")

biomolecule_df <- dbGetQuery(con, "SELECT *
                             FROM biomolecules
                             INNER JOIN metadata on metadata.biomolecule_id = biomolecules.biomolecule_id
                             WHERE omics_id = 3
                             AND (metadata_type = 'tier_mean' OR metadata_type = 'Remove_duplicate') ", stringsAsFactors = F)

dbDisconnect(con)

##### Converting to wide format ######

# being explicit about items in the data frame, in case db changes 
biomolecule_df_2[0,]
biomolecule_df_2 <- data.frame(biomolecule_id = biomolecule_df$biomolecule_id, 
                               standardized_name = as.character(biomolecule_df$standardized_name),
                               metadata_type = as.character(biomolecule_df$metadata_type), 
                               metadata_value = as.character(biomolecule_df$metadata_value),
                               stringsAsFactors = F)

bio_df_wide <- reshape(biomolecule_df_2, v.names = "metadata_value", timevar = "metadata_type", idvar = c("biomolecule_id"), direction = "wide")
row.names(bio_df_wide) <- bio_df_wide$biomolecule_id

str(bio_df_wide)
##### 1) Looking at duplicate features ##### 

bio_df_wide$metadata_value.Remove_duplicate <- bio_df_wide$metadata_value.Remove_duplicate == "T"

##### 2) Look at tier for features ######

# first converting data into numeric
bio_df_wide$metadata_value.tier_mean <- as.numeric(bio_df_wide$metadata_value.tier_mean)

## PLOT: Histogram of the mean teir 
hist(bio_df_wide$metadata_value.tier_mean, col = 1, breaks = 50,
     main = "Histogram of feature mean quant tier \n(1 = best, 5 = worst)",
     xlab = "Feature mean quant tier")

##### 3) Look at control RSD within batch and between batch #### 
## FUNCTION to calculate RSD 
rsd <- function(x) sd(2^x)/mean(2^x)

# select for controls
controls <- grepl("Control", df$unique_identifier)

# calculate rsd for controls by biomolecule and batch
rsd_controls_by_batch <- aggregate(df$normalized_abundance[controls], 
                                   by = list(batch = df$batch[controls], biomolecule_id = df$biomolecule_id[controls]), 
                                   rsd)
# calculate rsd for controls by biomolecule across all batches
rsd_controls <- aggregate(df$normalized_abundance[controls], 
                          by = list(biomolecule_id = df$biomolecule_id[controls]),
                          rsd)

# reshape by batch into wide format 
rsd_controls_by_batch_wide <- reshape(rsd_controls_by_batch, v.names = "x", timevar = "batch", idvar = c("biomolecule_id"), direction = "wide")
names(rsd_controls_by_batch_wide) <- sub("x.", "rsd_batch_", names(rsd_controls_by_batch_wide))

# append overall control rsds to df
rsd_controls_by_batch_wide$overall_rsd <- rsd_controls$x


## PLOT: boxplot of rsds by batch
boxplot(rsd_controls_by_batch_wide[,-1], col = c(1:7, "black"), las = 2, mar = c(8,5,5,2),
        main = "RSDs of QC samples by batch\nGC-metabolomics features",
        ylab = "RSD of QC samples")

## PLOT: plot of RSD by features
plot(rsd_controls_by_batch_wide[,1], rsd_controls_by_batch_wide[,2], col= 1, pch=19,
     main = "RSDs for QC by biomolecule_id",
     ylab = "RSD of QC samples", 
     xlab = "biomolecule_id")
for(i in 3:8){points(rsd_controls_by_batch_wide[,1], rsd_controls_by_batch_wide[,i], col =i-1, pch=19)}
points(rsd_controls_by_batch_wide[,1], rsd_controls_by_batch_wide$overall_rsd, col = "black", pch=19)
abline(h = 0.3, lty = 2)

# table for rsd > 0.3 for # of batches and whether or not they also have overall rsd >0.3 
table(rowSums(rsd_controls_by_batch_wide[,2:8]>0.3), rsd_controls_by_batch_wide$overall_rsd > 0.3)

# propose filtering out if overall RSD > 0.3 or more than 2 batches RSD is greater than 0.3 
filter <- rowSums(rsd_controls_by_batch_wide[,2:8]> 0.3)>2 | rsd_controls_by_batch_wide$overall_rsd > 0.3
table(filter)
bio_df_wide$standardized_name[filter]

##### 4) dynamic range of features ######

# extracting out whare are for patient samples
samples <- grepl("COVID", df$unique_identifier)

# function for calculating biomolecule dynamic range
dynamic_range <- function(x) max(x)-min(x)

# calculating dynamic range for samples 
dynamic_range_samples <- aggregate(df$normalized_abundance, by = list(biomolecule_id = df$biomolecule_id), dynamic_range)

## PLOT: Histogram of dynamic range for each biomolecule
hist(dynamic_range_samples$x, breaks = 50, main = "Dynamic range for GC-metabolites",
     xlab = "Dynamic range (log2(peak height))", col = 1)

# only one biomolecule had dynamic range less than 3
table(dynamic_range_samples$x < 3)
bio_df_wide$standardized_name[dynamic_range_samples$x < 3]
