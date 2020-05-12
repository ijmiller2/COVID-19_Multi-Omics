# COVID-19_Multi-Omics
A repository for collaborative multi-omics data analysis

### Directory structure
```bash
├── .env <- Store password, API keys here
├── .gitignore <- List of files and folders for git to ignore
├── README.md
├── data <- Keep local data here; no need to commit this
│   ├── external
│   ├── interim
│   ├── processed
│   └── raw
├── eda <- Directory for exploratory data analysis and notebooks
│   └── IJM
│       ├── lipidomics
│       │   └── README.md
│       ├── metabolomics
│       │   └── README.md
│       ├── mutlti-omics
│       │   └── README.md
│       └── proteomics
│           └── README.md
├── figures <- Keep figures here; no need to commit this
├── reference <- Data dictionaries, manuals, and all other explanatory materials.
│   ├── color_palette.txt
│   └── data_dictionary.txt
└── src <- Keep reusable functions and source code here
    ├── README.md
    ├── analysis
    │   ├── analysis.R
    │   └── analysis.py
    ├── db_queries
    │   ├── db_queries.R
    │   └── db_queries.py
    ├── processing
    │   ├── processing.R
    │   └── processing.py
    └── visualization
        ├── visualization.R
        └── visualization.py
```
