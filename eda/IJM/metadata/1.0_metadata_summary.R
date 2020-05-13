
library(ggplot2)
library(cowplot)
library(readxl)
library(ggstance) # https://stackoverflow.com/questions/40334056/horizontal-barchart-with-facet-grid-free-x-not-working

setwd("~/Desktop/UW_2020/COVID-19_multi-omics/COVID-19_Multi-Omics/")
metadata_df <- read_excel("data/processed/deidentified_covid_metadata.xlsx")
nrow(metadata_df) # 127

colnames(metadata_df)
#[1] "Sample ID"                    "Length of hospital stay"      "Hospital free days at day 28" "Gender"                      
#[5] "Age"                          "ICU=1,floor =2"               "floor to ICU=1"               "APACHEII"                    
#[9] "Charlson score"               "Mech Ventilation y=1 n=0"     "Days Vented"                  "Days of ICU"                 
#[13] "COVID, y=1, n=0"     

table(metadata_df$`ICU=1,floor =2`)
#1  2 
#66 61 
table(metadata_df$`COVID, y=1, n=0`)
#0   1 
#26 101 

# drop samples where stay is TBD
metadata_df <- subset(metadata_df, `Length of hospital stay` != "TBD" | `Hospital free days at day 28` != "TBD")
nrow(metadata_df) # 112

# convert to numeric 
metadata_df$`Length of hospital stay` <- as.numeric(metadata_df$`Length of hospital stay`)

# sort samples levels by covid diagnosis
#metadata_df$`Sample ID` <- factor(metadata_df$`Sample ID`, levels=metadata_df$`Sample ID`[order(metadata_df$`COVID, y=1, n=0`)])
# sort samples levels by covid diagnosis, then length of hospital stay
metadata_df$`Sample ID` <- factor(metadata_df$`Sample ID`, levels=metadata_df$`Sample ID`[with(metadata_df, order(`COVID, y=1, n=0`, `Length of hospital stay`, `Days Vented`))])

# Plot 1
# x = subjectid, y1 = days hospitalized, y2, days vented, facet by ICU, color by COVID +/-
#ggplot(metadata_df, aes(y=`Sample ID`, x=`Length of hospital stay`)) +
#  geom_barh(stat='identity', aes(fill=factor(`COVID, y=1, n=0`))) +
#  geom_point(data=metadata_df, x=`Days Vented`, y=`Sample ID`)
#  facet_grid(`ICU=1,floor =2` ~ ., scales = "free_y")

icu_df <- subset(metadata_df, `ICU=1,floor =2` == 1)
icu_days_vented_df <- icu_df[complete.cases(icu_df$`Days Vented`),]
icu_days_vented_df$`Days Vented` <- as.numeric(icu_days_vented_df$`Days Vented`)

a <- ggplot(icu_df, aes(x=`Sample ID`, y=`Length of hospital stay`)) +
    geom_bar(stat='identity', aes(fill=factor(`COVID, y=1, n=0`)), show.legend = F) +
    geom_point(data=icu_days_vented_df, aes(x=`Sample ID`, y=`Days Vented`)) +
    coord_flip() + 
  labs(title="ICU n=66", xlab="") +
  scale_fill_manual(name = "Diagnosis", labels = c("0", "1"), values=c("skyblue", "tomato"))
a

no_icu_df <- subset(metadata_df, `ICU=1,floor =2` == 2)
no_icu_days_vented_df <- no_icu_df[complete.cases(no_icu_df$`Days Vented`),]
no_icu_days_vented_df$`Days Vented` <- as.numeric(no_icu_days_vented_df$`Days Vented`)

b <- ggplot(no_icu_df, aes(x=`Sample ID`, y=`Length of hospital stay`)) +
  geom_bar(stat='identity', aes(fill=factor(`COVID, y=1, n=0`))) +
  geom_point(data=no_icu_days_vented_df, aes(x=`Sample ID`, y=`Days Vented`)) +
  coord_flip() + 
  labs(title="NO ICU n=61", xlab="") +
  scale_fill_manual(name = "Diagnosis", labels = c("NO-COVID\n(n=26)", "COVID\n(n=101)"), values=c("skyblue", "tomato"))
b

# NOTE: some of the sample id labels appear to be off (based on "Reason for admission" column)

# Plot 2, age and gender distributions
c <- ggplot(metadata_df[complete.cases(metadata_df$Gender),], aes(x=Age)) +
  geom_histogram(bins=20) +
  facet_grid(cols = vars(Gender))
c

# plot 3 - chronic disease burden

d <- ggplot(metadata_df, aes(x=factor(`ICU=1,floor =2`), y=`Charlson score`, fill=factor(`COVID, y=1, n=0`))) +
  geom_boxplot(show.legend = F) +
  labs(title="Charlson score", x="ICU (1) vs. NO-ICU (2)") +
  scale_fill_manual(name = "Diagnosis", labels = c("NO COVID", "COVID"), values=c("skyblue", "tomato"))
d

plot_grid(a, b, c, d)

e <-  ggplot(metadata_df, aes(x=factor(`ICU=1,floor =2`), y=`APACHEII`, fill=factor(`COVID, y=1, n=0`))) +
  geom_boxplot(show.legend = F) +
  labs(title="APACHEII Score")
e

