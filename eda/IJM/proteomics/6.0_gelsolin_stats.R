
# load proteomics data with clinical metadata
proteomics_data_path = "~/Desktop/UW_2020/COVID19_multi-omics/COVID-19_Multi-Omics/data/proteomics_measurements_w_clinical_metadata.csv"

# note the 'color_by' column represents the four patient groups:
#COVID_ICU    COVID_NONICU    NONCOVID_ICU NONCOVID_NONICU 
#51              51              15              10 

proteomics_df <- read.table(proteomics_data_path, sep=",", header=T)
covid_df <- subset(proteomics_df, COVID == 1)
no_covid_df <- subset(proteomics_df, COVID == 0)

#### 1.	COVID	pGSC ~ Hospital Free Days at Day 45 ####

test_1 <- lm(X7974 ~ Hospital_free_days_45, data=covid_df)
summary(test_1)

# 1.03e-06
# Multiple R-squared:  0.2133
# dof: 100

#### 2.	NO COVID	pGSC ~ Hospital Free Days at Day 45 ####

test_2 <- lm(X7974 ~ Hospital_free_days_45, data=no_covid_df)
summary(test_2)

# p-value: 0.0388
# Multiple R-squared:  0.1728
# dof: 23

#### 3.	COVID_NONICU-COVID_ICU ####

gelsolin.lm <- lm(X7974 ~ color_by, data = proteomics_df)
gelsolin.av <- aov(gelsolin.lm)
summary(gelsolin.av)

test_3 <- TukeyHSD(gelsolin.av)
test_3

# p-value: 4e-07

#### 4. COVID	pGSC ~ Vent Free Days at Day 28 ####

test_4 <- lm(X7974 ~ Vent_free_days, data=covid_df)
summary(test_4)

# p-value: 2.35e-05
# Multiple R-squared:  0.1645
# dof: 100

#### 5. NO COVID	pGSC ~ Vent Free Days at Day 28 ####

test_5 <- lm(X7974 ~ Vent_free_days, data=no_covid_df)
summary(test_5)

# p-value: 0.371
# Multiple R-squared:  0.03492
# dof: 23

#### 6. COVID	pGSC ~ APACHE II ####

test_6 <- lm(X7974 ~ APACHEII, data=covid_df)
summary(test_6)

# p-value: 0.000992
# Multiple R-squared:  0.1774
# dof: 56

#### 7. COVID	pGSC ~ SOFA ####

test_7 <- lm(X7974 ~ SOFA, data=covid_df)
summary(test_7)

# p-value: 3.04e-05
# Multiple R-squared:  0.2732
# dof: 55

#### 8. COVID	pGSC ~ SAPS II ####

test_8 <- lm(X7974 ~ SAPS2, data=covid_df)
summary(test_8)

# p-value: 3.04e-05
# Multiple R-squared:  0.2732
# dof: 55

#### 9. COVID	pGSC ~ P/F ratio ####

test_9 <- lm(X7974 ~ PF_ratio_numeric, data=covid_df)
summary(test_9)

# p-value: 0.0253
# Multiple R-squared:  0.102
# dof: 47