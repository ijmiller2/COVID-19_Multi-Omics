# COVID-19_Multi-Omics
A repository for collaborative multi-omics data analysis and the covid-omics.app.

Note that the relational SQLite database is available for download via ftp from [MassIVE - accession number MSV000085703](https://doi.org/10.25345/C5F74G) under the "other" directory. For more information on data processing, analysis, and availability, please see the [open access manuscript](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30371-9#secsectitle0105).

### Directory structure
```bash
├── README.md
├── eda
│   ├── BJA
│   ├── EAT
│   ├── IJM
│   ├── JGM
│   ├── KAO
│   ├── MNB
│   └── VL
├── reference
│   ├── color_palette.txt
│   ├── data_dictionary.txt
│   └── db_schema.png
└── src
    ├── README.md
    ├── dash
        ├── README.md
        ├── app.py
        ├── apps
        │   ├── clustergrammer.py
        │   ├── differential_expression.py
        │   ├── linear_regression.py
        │   └── pca.py
        ├── assets
        │   └── favicon.ico
        ├── change_log.md
        ├── data.py
        ├── index.py
        ├── nav.py
        ├── plot.py
        └── requirements.txt
    ├── db_queries
        ├── db_queries.R
        └── db_queries.py
```
