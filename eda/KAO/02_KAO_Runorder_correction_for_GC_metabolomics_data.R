##### 02 Run-time correction for GC metabolomics data #####

## A common occurance in daily GC runs is signal loss or gain 
## throughout the day's analysis. We assume this is due in part
## to changing derivitization over time. In many cases it is minor
## but for a handful of metabolite features, run-time correction 
## can make a significant difference in data reproducibility.
## 
## Below is code to perform run time correction, this uses 
## a robust linear regression to minimize effect of outliers.
## 


library(DBI)
library(RSQLite)
library(MASS)

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### List tables in db #####

dbListTables(con)
dbReadTable(con, "rawfiles")
dbReadTable(con, "metabolomics_runs")
dbReadTable(con, "metabolomics_measurements")

##### 
timeStamp_df<- dbGetQuery(con, "SELECT unique_identifier, timestamp, raw_abundance, biomolecule_id, batch
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           WHERE keep = 1
           ")

dbDisconnect(con)

#### Adding variables to data frame #####
# converting timestamp to POSIXct class
timeStamp_df[,2] <- as.POSIXct(strptime(timeStamp_df[,2], "%Y%m%d%H%M%S"))

# calculating difftime for linear regression model
timeStamp_df$difftime <- as.numeric(difftime(timeStamp_df$timestamp, timeStamp_df$timestamp[1]))

# creating columne where the log2 raw abudance values are in non-log 
timeStamp_df$raw_abundance_nonlog <- 2^timeStamp_df$raw_abundance

# creating an NA column that will eventually hold the run-time corrected values
timeStamp_df$runtime_corrected <- NA

# creating a weights column that more heavily weights control values in the rlm model
timeStamp_df$weights <- 0
timeStamp_df$weights[grepl("Control", timeStamp_df$unique_identifier)] <- 4

#### For assessment of the rlm approach calculating the rsd of QC samples pre and post runt time correction 
rsd_pre <- matrix(NA, nrow = length(unique(timeStamp_df$batch)), ncol = length(unique(timeStamp_df$biomolecule_id)))
rsd_post <-matrix(NA, nrow = length(unique(timeStamp_df$batch)), ncol = length(unique(timeStamp_df$biomolecule_id)))

#### test code for troubleshooting fits ##### 
# plot(timeStamp_df$raw_abundance_nonlog[index]~timeStamp_df$difftime[index])
# abline(fit)

#### For loop for calculating runtime corrected values #####
set.seed(1234)

for(i in unique(timeStamp_df$biomolecule_id)){
  for(k in unique(timeStamp_df$batch)){
    index <- timeStamp_df$biomolecule_id == i & timeStamp_df$batch == k
    rsd_pre[k,i] <- sd(timeStamp_df$raw_abundance_nonlog[index & grepl("Control", timeStamp_df$unique_identifier)])/mean(timeStamp_df$raw_abundance_nonlog[index & grepl("Control", timeStamp_df$unique_identifier)])
    fit <- tryCatch(rlm(raw_abundance_nonlog~difftime, data = timeStamp_df[index,], weights = weights),warning = function(w) fit =NULL, error = function(e) fit = NULL) 
    if(!is.null(fit) ){
      intercept <- summary(fit)$coef[1,1]
      slope <- summary(fit)$coef[2,1]
    } else {
      intercept <- mean(timeStamp_df$raw_abundance_nonlog[index], na.rm =T)
      slope <- 0
    }
    # normalize by run-order
    # sampleValue - ((slope * runOrder) + intercept) + biomolecule mean
    new_value <- (timeStamp_df$raw_abundance_nonlog[index] - ((slope * timeStamp_df$difftime[index]) + intercept) )+ mean(timeStamp_df$raw_abundance_nonlog[timeStamp_df$biomolecule_id == i & grepl("Control", timeStamp_df$unique_identifier)], na.rm =T)
    # in cases where run-order corrected value is below 0, replacing with pre-corrected values
    # this step was decided after carefully looking at instances of negative values which were very low pre correction. 
    # in total only ~ 1000 values have this occurance 
    new_value[new_value<0] <- timeStamp_df$raw_abundance_nonlog[index][new_value<0]
    timeStamp_df$runtime_corrected[index] <-log2(new_value)
    rsd_post[k,i] <- sd(2^timeStamp_df$runtime_corrected[index & grepl("Control", timeStamp_df$unique_identifier)])/mean(2^timeStamp_df$runtime_corrected[index & grepl("Control", timeStamp_df$unique_identifier)])
    
  }
}


##### Assessing whether the run time improved data quality ##### 

## PLOT: runtime_corrected value vs. time
par(las = 1)
plot(timeStamp_df$runtime_corrected~ timeStamp_df$difftime, 
     main = "Distribution of intensity values across runs",
     ylab = "log2(run-time corrected intensity)", 
     xlab = "run time (sec)")

## PLOT: runtime_corrected value vs. original value
plot(timeStamp_df$runtime_corrected, timeStamp_df$raw_abundance, col = as.numeric(timeStamp_df$batch),
     main = "Comparing run-time corrected vs. original values", 
     ylab = "log2(intensity) Original",
     xlab = "log2(intensity) Run-time corrected")
legend("topleft", col = 1:7, paste("Batch",unique(timeStamp_df$batch)), pch =1, ncol = 2, bty = "n")

## PLOT: RSD of control samples for each batch pre and post
par(mfrow=c(2,1), mar = c(4,5,4,4))
boxplot(t(rsd_pre), ylim = c(0,2.5), col = 1:7, 
        ylab = "RSD Control Samples",
        xlab = "Batch",
        main = "RSDs for QC samples original")
boxplot(t(rsd_post),  ylim = c(0,2.5), col = 1:7, 
        ylab = "RSD Control Samples",
        xlab = "Batch",
        main = "RSDs for QC samples after run-time correction")

## PLOT: change in RSD of control samples after run-time correction
par(mfrow = c(1,1))
hist(rsd_pre - rsd_post, breaks =100, col = 2, 
     main ="Change in metabolite RSD w/in controls \nwith run-time correction", 
     xlab = "Change in RSD (pre - post correction)") 


#### Creating a wide-format data frame to facilitate PCA #####
pre <- reshape(timeStamp_df, timevar = "biomolecule_id", v.names = "raw_abundance",
               idvar = "unique_identifier", direction = "wide" )
post <- reshape(timeStamp_df, timevar = "biomolecule_id", v.names = "runtime_corrected",
               idvar = "unique_identifier", direction = "wide" )
pre_exprs <- pre[,-c(1:7)]
post_exprs <- post[,-c(1:7)]

#### Pricipal component analysis pre and post correction #####
# PCA analysis of samples with original intensity values
pca_pre <- prcomp(pre_exprs)

# PCA anlysis of samples with run-time corrected values 
pca_post <- prcomp(post_exprs)

## PLOT: PCA pre and PCA post correction 
par(mfrow = c(1,2), mar = c(5,4,4,2))
plot(pca_pre$x, col = pre$batch, pch = 19, #xlim = c(-30,30), ylim = c(-30,30),
     main = "PCA original")
legend("topleft", col = 1:7, paste("Batch",unique(timeStamp_df$batch)), 
       pch =19, ncol = 1, bty = "n", cex = 0.7)

plot(pca_post$x, col =pre$batch, pch = 19, #xlim = c(-30,30), ylim = c(-30,30),
     main = "PCA run-time corrected")
legend("topleft", col = 1:7, paste("Batch",unique(timeStamp_df$batch)), 
       pch =19, ncol = 1, bty = "n", cex = 0.7)
