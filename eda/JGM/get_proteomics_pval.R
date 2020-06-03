library(RSQLite)
con <- dbConnect(SQLite(), dbname = "C:/covid_proteomics/Covid-19 Study DB.sqlite")

dbListTables(con)

# Extract tables from SQLite database. The square prakets provide the first row or column names of table. 
dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "metadata")[1:10,]
dbReadTable(con, "proteomics_measurements")[0,]
dbReadTable(con, "metabolomics_runs")[0,]
dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "deidentified_patient_metadata")     
dbReadTable(con, "pvalues")[0,]
dbReadTable(con, "omes")
dbReadTable(con, "rawfiles")[,3]
dbReadTable(con, "rawfiles")[700:770,]

df<- dbGetQuery(con, "SELECT q_value, pvalues.biomolecule_id, test, biomolecules.standardized_name
           FROM pvalues
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = pvalues.biomolecule_id
           INNER JOIN proteomics_measurements ON proteomics_measurements.biomolecule_id = biomolecules.biomolecule_id
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           WHERE rawfiles.keep = 1 
           AND ome_id = 1
           ")

df<- dbGetQuery(con, "SELECT q_value, comparison, pvalues.biomolecule_id, test, biomolecules.standardized_name
           FROM pvalues
           INNER JOIN biomolecules ON biomolecules.biomolecule_id = pvalues.biomolecule_id
           AND omics_id= 1")



head(df)
nrow(df)

proteinslist<-strsplit(unlist(df$standardized_name), ";")

getfirst=function(x){unlist(x)[1]}

df$firstprot<-unlist(lapply(FUN=getfirst, proteinslist))
head(df)


### read the peptide correlation analysis results
setwd("C:/covid_proteomics/pecora/no_control/")

o<-read.delim("out.txt", stringsAsFactors = F)

head(o)
?merge
head(reg.out)
merged_pecora<-merge(reg.out, df, by.x="protein", by.y='firstprot', all.x=T)
head(merged_pecora)
write.table(merged_pecora, "pecora_merged.txt", 
            quote=F, sep="\t",
            row.names = F)

