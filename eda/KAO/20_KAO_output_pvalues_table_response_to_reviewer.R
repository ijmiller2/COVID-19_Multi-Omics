### NA_comparing_pvalues_robustlm_vs_lm.R #### 


library(DBI)
library(RSQLite)


con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")
# pull pvalues 

pvalues_COVID_robust <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, p_value, q_value, omics_id
          FROM pvalues 
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.comparison = 'COVID' AND pvalues.test = 'LR_test_robust' AND biomolecules.keep = 1
                  ")

pvalues_COVID <- dbGetQuery(con, "SELECT  pvalues.biomolecule_id, p_value, q_value
          FROM pvalues 
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.comparison = 'COVID_vs_NONCOVID' AND biomolecules.keep = 1
                  ")

pvalues_HFD_robust <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, p_value, q_value
          FROM pvalues 
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.formula = 'Hospital_free_days_45 ~ normalized_abundance + Age_less_than_90 + Gender vs. Hospital_free_days_45 ~ Age_less_than_90 + Gender' AND pvalues.test = 'LR_test_robust' AND biomolecules.keep = 1
                  ")
pvalues_HFD <- dbGetQuery(con, "SELECT pvalues.biomolecule_id, p_value, q_value
          FROM pvalues  
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.formula = 'Hospital_free_days_45 ~ normalized_abundance + Age_less_than_90 + Gender vs. Hospital_free_days_45 ~ Age_less_than_90 + Gender' AND biomolecules.keep = 1 AND pvalues.test = 'LR_test'
                  ")
biomolecules <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name
                           FROM biomolecules
                           WHERE keep = 1")

# disconnect
dbDisconnect(con) 


###  Reading in Log2(fold change covid/non-covid)
FC <- read.csv("data/COVID_fc_by_biomolecule_ID.csv")

## load in transcript pvalues 

transcript_pvalues <- read.delim("P:/All_20200428_COVID_plasma_multiomics/Transcriptomics/AllGenes.txt", sep = " ")
transcript_pvalues <- data.frame(standardized_name = row.names(transcript_pvalues), PP = transcript_pvalues[,1], stringsAsFactors = F)

## Elastic net results ## 

EN <- do.call(rbind, lapply(list.files("./data/ElasticNet", "tsv", recursive = T, full.names=TRUE), read.delim, stringsAsFactors = F))
names(EN) <- c("standardized_name", "EN_coefficient")


#### merge into one data frame #### 

pvalues_merge <- merge(pvalues_COVID, pvalues_COVID_robust, by = "biomolecule_id", all = T)
names(pvalues_merge) <- c("biomolecule_id", "p_value_COVID", "q_value_COVID", "p_value_COVID_robust", "q_value_COVID_robust", "omics_id")

pvalues_merge2 <- merge(pvalues_merge, pvalues_HFD, by = "biomolecule_id", all = T)
names(pvalues_merge2)[7:8] <- c("p_value_HFD", "q_value_HFD")
head(pvalues_merge2)

pvalues_merge3 <- merge(pvalues_merge2, pvalues_HFD_robust, by = "biomolecule_id", all = T)
names(pvalues_merge3)[9:10] <- c("p_value_HFD_robust", "q_value_HFD_robust")

pvalues_merge4 <- merge(pvalues_merge3, biomolecules, by = "biomolecule_id", all = T)
names(pvalues_merge4)

pvalues_merge5 <- merge(pvalues_merge4, transcript_pvalues, by = "standardized_name", all = T)
names(pvalues_merge5)[12] <- "Postior_probabilty_EBseq_COVID"

pvalues_merge6 <- merge(pvalues_merge5, EN, by = "standardized_name", all = T)
names(pvalues_merge6)

names(FC)
FC_2 <- FC[,c(2,7)]

pvalues_merge7 <- merge(pvalues_merge6, FC_2, by = "biomolecule_id", all = T)
names(pvalues_merge7)

write.csv(pvalues_merge7, file ="data/pvalues_all_20200928_for_supplemental_table.csv")


pdf("plots/COVID_pvalues_robustlm_vs_lm.pdf", useDingbats = F)
plot(-log(pvalues_merge$q_value_COVID), -log(pvalues_merge$q_value_COVID_robust))
abline(0,1, col = "red")
#points(-log(pvalues_merge$q_value_COVID)[pvalues_merge$q_value_COVID < 0.05 & pvalues_merge$q_value_COVID_robust > 0.05], -log(pvalues_merge$q_value_COVID_robust)[pvalues_merge$q_value_COVID < 0.05 & pvalues_merge$q_value_COVID_robust > 0.05], pch = 19 ,col = "blue")
#legend("topleft", pch = 19, col = "blue", paste("n =", table(pvalues_merge$q_value_COVID < 0.05 & pvalues_merge$q_value_COVID_robust > 0.05)[2]))
dev.off()

pvalues_HFD_merge <- merge(pvalues_HFD, pvalues_HFD_robust, by = "biomolecule_id")
names(pvalues_HFD_merge) <- c("biomolecule_id", "p_value_HFD", "q_value_HFD", "p_value_HFD_robust", "q_value_HFD_robust")

plot(-log(pvalues_HFD_merge$q_value_HFD), -log(pvalues_HFD_merge$q_value_HFD_robust),
     ylab = "-log(q-values) calculated with robust linear regression",
     xlab = "-log(q-values) calculated with standard linear regression")
abline(0,1, col = "red")
dev.off()
