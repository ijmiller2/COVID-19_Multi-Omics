# .R code for GC Metabolomics and LC-MS/MS Proteomics

### GC_metabolomics_heatmap.R

**description**:   
    1. Load data from SQLlite data base (Metabolomics normalized measurements, deidentified patient meta data, and metabolite class meta)  
    2. Data formatting from long to wide format  
    3. Meta formattting to pad all single digit patient id's with a zero  
    4. Generate heatmap with annotations (177 metabolites) for Controls, NONCOVID, and COVID patients  
    5. Remove Controls  
    6. Exploratory figures for COVID and NONCOVID metabolomics  
    7. Fold Changes to median NONCOVID  
    8. Filter features witth fold changer greater than 4.2 in at least one patient  
    
**Relevant Issue(s)**: Unsupervised analysis #3  
**date created**: 5/26/20  
**date last modified**: 5/26/20  
**input**:  
  - sqlite_db: data/SQLite Database/Covid-19 Study DB.sqlite  

**output**:  
  - NA  
 
### Proteomics_heatmap.R

**description**:   
    1. Load data from SQLlite data base (Proteomics normalized measurements, deidentified patient meta data, and proteomics meta)  
    2. Data formatting from long to wide format  
    3. Meta formattting to pad all single digit patient id's with a zero  
    4. Fix raw file names for HC and pooled plasma.
    5. Generate heatmap with annotations (517 proteins) for Human Controls (HC), pooled plasma, NONCOVID, and COVID patients  
    6. Remove HC and pooled plasma samples  
    7. Exploratory figures for COVID and NONCOVID proteomics  
    8. Fold Changes to median NONCOVID  
    9. Fold Changes to median NONCOVID not in the ICU 
    
**Relevant Issue(s)**: Unsupervised analysis #3  
**date created**: 5/27/20  
**date last modified**: 5/27/20  
**input**:  
  - sqlite_db: data/SQLite Database/Covid-19 Study DB.sqlite  

**output**:  
  - NA 

### Protein_peptide_90minSingleShot_COVID_CoagulationCascade.R 

description: Generates csv's for each protein in MaxQuant proteingroups.csv, where each csv contains the  peptide measured in MS/MS analysis for that specific protein.
Relevant Issue(s): NA
date created: 5/26/20
date last modified: 5/26/20  
**input**:  
  - protein groups csv: P:/All_20200428_COVID_plasma_multiomics/Proteomics/90 min single-shot reps for Anji/combined_includingIsoforms_Deane_coagulationPathway/txt/peptides_Plasma_90SingleShot.csv  
  - peptides csv: P:/All_20200428_COVID_plasma_multiomics/Proteomics/90 min single-shot reps for Anji/combined_includingIsoforms_Deane_coagulationPathway/txt/peptides_Plasma_90SingleShot.csv
  - Uniprot gene csv: P:/All_20200428_COVID_plasma_multiomics/Proteomics/90 min single-shot reps for Anji/UniprotDatabase_CoagulationCascade_List.csv  

**output**:  
  - peptide csv for each gene: All_20200428_COVID_plasma_multiomics/Proteomics/90 min single-shot reps for Anji/CoagulationCascade_peptide_data/
