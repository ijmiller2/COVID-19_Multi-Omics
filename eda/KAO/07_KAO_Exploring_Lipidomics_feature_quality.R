###### 07_KAO_Exploring_lipidomics_feature_quality ##### 

## To filter for quality lipidomics features will:
## 1) exlclude features with RT less than 1.15
## 2) explore if feature was quantified in control 
## 3) determine intra and inter-batch CVs
## 4) look at slopes of quarter and half samples (only if in control)
## 5) imputed in many samples? 

library(DBI)
library(RSQLite)
colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

df<- dbGetQuery(con, "SELECT unique_identifier, normalized_abundance, lipidomics_measurements.biomolecule_id, batch 
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           WHERE keep = 1
           ")

biomolecule_df <- dbGetQuery(con, "SELECT *
                             FROM biomolecules
                             INNER JOIN metadata on metadata.biomolecule_id = biomolecules.biomolecule_id
                             WHERE omics_id = 2", stringsAsFactors = F)

dbDisconnect(con)

##### Load in data w/ slope of control dilution and whether the feature was found in control ##### 

df_slope <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Lipidomics/Lipidomics_quant_results/Final_Results_Slopes_KAO.csv", stringsAsFactors = F)

##### Converting to wide format ######

# being explicit about items in the data frame, in case db changes 
biomolecule_df_2 <- data.frame(biomolecule_id = rep(unique(df$biomolecule_id), each = length(unique(biomolecule_df$metadata_type) )), 
                               standardized_name = as.character(biomolecule_df$standardized_name),
                               metadata_type = as.character(biomolecule_df$metadata_type), 
                               metadata_value = as.character(biomolecule_df$metadata_value),
                               stringsAsFactors = F)

bio_df_wide <- reshape(biomolecule_df_2, v.names = "metadata_value", timevar = "metadata_type", idvar = c("biomolecule_id"), direction = "wide")
row.names(bio_df_wide) <- bio_df_wide$biomolecule_id

# checking that these dimentions match with the df_slope dim 
dim(bio_df_wide) ## incorrect;  7235 8
dim(df_slope) # 7235 11
length(unique(df$biomolecule_id)) #7235
dim(metadata_df) # 7235

##### combining meta data into one data frame #### 

bio_df_wide$slope_of_control <- df_slope$slope

bio_df_wide$Present_in_control1_batch4 <- df_slope$Present_in_control1_batch4

# clean up environment
rm(df_slope)

##### Starting a keep df - will update throughout this script ####
df_keep <- data.frame(biomolecule_id = unique(biomolecule_df$biomolecule_id),
                      keep = 1)

##### 1) Exploring number of features below 1.15 ###### 
## looking at Compound discoverer data viewer, features less than 1.15 RT had 
## incomplete points across the chromatogram 

RT_too_low <- bio_df_wide$`metadata_value.Retention Time (min)` < 1.15
table(RT_too_low) # 32 true

# modify the keep string; if 1 turn to 0, then append "RT_too_low"
df_keep$keep[RT_too_low & df_keep$keep == "1"] <- sub("1", "0", df_keep$keep[RT_too_low & df_keep$keep == "1"])

df_keep$keep[RT_too_low] <- paste(df_keep$keep[RT_too_low], "RT_too_low", sep = ";")        

##### 2) Explore features quantified in control ##### 
bio_df_wide$Present_in_control1_batch4[is.na(bio_df_wide$Present_in_control1_batch4)] <- 0
in_control <-bio_df_wide$Present_in_control1_batch4 == 1
table(in_control) #1853

##### 3) Explore intra and inter batch cv for features found in control 

## FUNCTION to calculate RSD 
rsd <- function(x) sd(2^x)/mean(2^x)

# select for controls
controls <- grepl("Control", df$unique_identifier) & !grepl("quarter|half", df$unique_identifier)

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
boxplot(rsd_controls_by_batch_wide[in_control,-1], col = c(1:7, "black"), las = 2, mar = c(8,5,5,2),
        main = "RSDs of QC samples by batch\nLipidomics features",
        ylab = "RSD of QC samples")

## PLOT: plot of RSD by features
plot(rsd_controls_by_batch_wide[in_control,1], rsd_controls_by_batch_wide[in_control,2], col= 1, pch=19,
     main = "RSDs for QC by biomolecule_id",
     ylab = "RSD of QC samples", 
     xlab = "biomolecule_id")
for(i in 3:8){points(rsd_controls_by_batch_wide[in_control,1], rsd_controls_by_batch_wide[in_control,i], col =i-1, pch=19)}
points(rsd_controls_by_batch_wide[in_control,1], rsd_controls_by_batch_wide$overall_rsd[in_control], col = "black", pch=19)
abline(h = 0.3, lty = 2)

# table for rsd > 0.3 for # of batches and whether or not they also have overall rsd >0.3 
table(rowSums(rsd_controls_by_batch_wide[in_control,2:8]>0.3), rsd_controls_by_batch_wide$overall_rsd[in_control] > 0.3)

# propose filtering out if overall control RSD > 0.3 or more than 2 batches RSD is greater than 0.3 
# only for features found in control
filter <- in_control & (rowSums(rsd_controls_by_batch_wide[,2:8]> 0.3)>1 | rsd_controls_by_batch_wide$overall_rsd > 0.3)
table(filter)
bio_df_wide$standardized_name[filter]

QC_intrabatch_RSD_over30perc <- in_control & rowSums(rsd_controls_by_batch_wide[,2:8]> 0.3)>1
QC_interbatch_RSD_over30perc <- in_control & rsd_controls_by_batch_wide$overall_rsd > 0.3

table(QC_interbatch_RSD_over30perc, QC_intrabatch_RSD_over30perc)
## Update keep df

# modify the keep string; if 1 turn to 0, then append "QC_intrabatch_RSD_over30perc"
df_keep$keep[QC_intrabatch_RSD_over30perc & df_keep$keep == "1"] <- sub("1", "0", df_keep$keep[QC_intrabatch_RSD_over30perc & df_keep$keep == "1"])

df_keep$keep[QC_intrabatch_RSD_over30perc] <- paste(df_keep$keep[QC_intrabatch_RSD_over30perc], "QC_intrabatch_RSD_over30perc", sep = ";")        


# modify the keep string; if 1 turn to 0, then append "QC_interbatch_RSD_over30perc"
df_keep$keep[QC_interbatch_RSD_over30perc & df_keep$keep == "1"] <- sub("1", "0", df_keep$keep[QC_interbatch_RSD_over30perc & df_keep$keep == "1"])

df_keep$keep[QC_interbatch_RSD_over30perc] <- paste(df_keep$keep[QC_interbatch_RSD_over30perc], "QC_interbatch_RSD_over30perc", sep = ";")    


##### 4) Explore if slopes are negative in control samples at 1/2 dilution and 1/4 dilution ######

negative_QC_slope <- bio_df_wide$slope_of_control < 0 & in_control

table(negative_QC_slope, QC_interbatch_RSD_over30perc|QC_intrabatch_RSD_over30perc)

bio_df_wide$standardized_name[negative_QC_slope] ## lots of TGs 

## PLOT: lipids with negative slope 
barplot(table(negative_QC_slope, bio_df_wide$`metadata_value.Lipid Class`)[,-1],
        main = "Lipids with negative slope in QC dilution series", las = 2)


plot(1:length(which(bio_df_wide$`metadata_value.Lipid Class` == "TG")), 
     bio_df_wide$`metadata_value.Area (max)`[bio_df_wide$`metadata_value.Lipid Class` == "TG"], 
     pch = 19, col = as.numeric(negative_QC_slope[bio_df_wide$`metadata_value.Lipid Class` == "TG"])+1,
     ylab = "max area", 
     xlab = "TGs",
     main = "QC slope is negative (orange)")


plot(bio_df_wide$`metadata_value.Retention Time (min)`[bio_df_wide$`metadata_value.Lipid Class` == "TG"], 
     bio_df_wide$`metadata_value.Features Found`[bio_df_wide$`metadata_value.Lipid Class` == "TG"], 
     pch = 19, col = as.numeric(negative_QC_slope[bio_df_wide$`metadata_value.Lipid Class` == "TG"])+1,
     ylab = "features found", 
     xlab = "TGs - Rention time (min)",
     main = "QC slope is negative (orange)")

## Well, probably better filter out at this time... it could be due to ionization competion 

# modify the keep string; if 1 turn to 0, then append "negative_QC_slope"
df_keep$keep[negative_QC_slope & df_keep$keep == "1"] <- sub("1", "0", df_keep$keep[negative_QC_slope & df_keep$keep == "1"])

df_keep$keep[negative_QC_slope] <- paste(df_keep$keep[negative_QC_slope], "negative_QC_slope", sep = ";")    

##### 5) imputed in many samples? ######
# to explore imputation will use 2 lines of investigation, 
# 1) bimodal distributions, 2) low number of features found

library(gamlss.mx)

hist(as.numeric(bio_df_wide$`metadata_value.Features Found`), breaks = 100)

biomodal_test <- function(x) {
        fit <- tryCatch(gamlssMX(normalized_abundance ~ 1, data = df[df$biomolecule_id == x[1], ], family = NO, K = 2),  error = function(e) NULL)
        if(is.null(fit)) c(NA, NA, NA, NA) else c(fit$prob, unname(fit$models[[1]]$mu.coefficients), unname(fit$models[[2]]$mu.coefficients))
}

#iteratively b/c issues with n digits
bimodal_probabilities_1 <- apply(df_keep[1:754,], 1, biomodal_test)
bimodal_probabilities_1p5 <- apply(df_keep[755:774,], 1, biomodal_test)
bimodal_probabilities_2 <- apply(df_keep[775:nrow(df_keep),], 1, biomodal_test) ## something weird about number of digits 

# combine resuslt from bimodal probabilites, row 1 and 2 is probababy dist 1 and 2. row 3 and 4 is means for dist 1 and 2. 
bimodal_probabilities <- cbind(bimodal_probabilities_1,bimodal_probabilities_1p5,bimodal_probabilities_2)

plot(bimodal_probabilities[1, ]~ bio_df_wide$`metadata_value.Features Found`)
plot(bimodal_probabilities[2,], bimodal_probabilities[4,])

# filter for things which both means are below 10, 
# if means are distint (2 log2 apart) then filter if prbobability of means less than 15 are greater than .75. 
imputationHigh <- bimodal_probabilities[3,]-bimodal_probabilities[4,] >2 &
        !(bimodal_probabilities[3,] > 15 & bimodal_probabilities[4,] > 15) &
        ((bimodal_probabilities[3, ] < bimodal_probabilities[4, ] & bimodal_probabilities[3, ] < 15 & bimodal_probabilities[1,] > 0.75) |
        (bimodal_probabilities[3, ] > bimodal_probabilities[4, ] & bimodal_probabilities[4,] < 15 & bimodal_probabilities[2,] > 0.75 )) | 
        (bimodal_probabilities[3,] < 10 & bimodal_probabilities[4,] < 10) 

imputationHigh[is.na(imputationHigh)] <- FALSE # there are a few na fits 

# features found is less than 4 for unidentified features 
few_features_found <- (bio_df_wide$`metadata_value.Features Found` < 4 & bio_df_wide$`metadata_value.Lipid Class` == "")


# modify the keep string; if 1 turn to 0, then append "Imputation_high"
df_keep$keep[imputationHigh & df_keep$keep == "1"] <- sub("1", "0", df_keep$keep[imputationHigh & df_keep$keep == "1"])

df_keep$keep[imputationHigh] <- paste(df_keep$keep[imputationHigh], "Imputation_high", sep = ";")        


# modify the keep string; if 1 turn to 0, then append "unknowns_with_few_features_found"
df_keep$keep[few_features_found & df_keep$keep == "1"] <- sub("1", "0", df_keep$keep[few_features_found & df_keep$keep == "1"])

df_keep$keep[few_features_found] <- paste(df_keep$keep[few_features_found], "Unknown_with_fewer_than_4_features_found", sep = ";")    


write.csv(df_keep, file = "P:/All_20200428_COVID_plasma_multiomics/Lipidomics/Lipidomics_quant_results/Final_keep_table.csv")


###### updating keep table in database #######
df_keep <-  df_keep[df_keep$keep != "1",]


# funtion to replace biomolecules keep values
replaceValues <- function(x)  dbExecute( con, paste("UPDATE biomolecules
          SET keep = '",unlist(x[2]),"' WHERE biomolecule_id = ", x[1], sep = "")) 

##### Iterate over smaller data frame to update biomolecule keep values in sqlite db ######

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#dbExecute(con, "UPDATE biomolecules SET keep = '1'")                        

apply(df_keep, 1, replaceValues)

## confirm 
biomolecule_df <- dbGetQuery(con, "SELECT *
          FROM biomolecules
          ")

dbDisconnect(con)

