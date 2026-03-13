# A Review of AMR of NTS in Africa

## Project Overview

This repository contains all R Markdown (`.Rmd`) scripts used for the review of *Salmonella* antimicrobial resistance (AMR) across Africa. Each `.Rmd` file is self-contained — it loads data, runs the analysis, and renders all figures inline. Knitting any `.Rmd` produces a fully reproducible HTML report with embedded plots.

**Important:** All `.Rmd` files must be placed in the same folder as `AMR_clean.csv`. The setup chunk in each file sets the working directory — update `root.dir` to match your local path before knitting.

---

## Key Findings

- **Temporal trend:** AMR prevalence increased significantly over time (quasibinomial GLM: coefficient = 0.075, *p* < 0.001)
- **Best-fit model:** Region, species, antimicrobial compound, and Antimicrobial susceptibility testing  method all independently explained resistance variation
- **Regional resistance:** Tetracyclines, sulfonamides, and penicillins showed the highest pooled resistance prevalence across regions
- **MDR patterns:** Multi-drug resistance varied substantially by host species and antibiotic class
- **Geographic distribution:** Studies concentrated in East and North Africa, with notable data gaps in Central Africa

---

## Repository Structure

```
📁 01_Data/
   └── AMR_sal.csv
    └── AMR_clean.csv

📁 02_RScripts/
   └── 01_map_of_studies.Rmd
    └── 02_temporal_trends_analysis.Rmd
     └── 03_regional_resistance_patterns.Rmd
      └── 04_antibiotic_class_resistance.Rmd
       └── 05_multidrug_resistance.Rmd
        └── 06_AMR_africa_drug_acronyms.Rmd

📁 03_Figures/
   └── Figure_1b
    └── Figure_2
     └── Figure_3
      └── Figure_4
       └── Figure_5
        └── Figure_6
         SupplementaryFigure_1
         SupplementaryTable_1

```

---

## Reproducing the Analysis

### 1. Prerequisites

```r
install.packages(c(
  "tidyverse", "janitor", "ggpubr", "patchwork", "ggthemes",
  "tmap", "leaflet", "sf", "rnaturalearth", "viridis", "mapview",
  "RColorBrewer", "wesanderson", "ggsci", "grid", "gridExtra",
  "flextable", "officer", "gt", "kableExtra", "tinytex",
  "AER", "lme4", "boot", "rmarkdown", "knitr"
))
```


### 2. Run Order

```r
rmarkdown::render("DataCleaning/00_data_cleaning.Rmd")
rmarkdown::render("MapOfStudies/01_map_of_studies.Rmd")
rmarkdown::render("TemporalTrends/02_temporal_trends_analysis.Rmd")
rmarkdown::render("RegionalResistance/03_regional_resistance_patterns.Rmd")
rmarkdown::render("AntibioticClassResistance/04_antibiotic_class_resistance.Rmd")
rmarkdown::render("MDRAnalysis/05_multidrug_resistance.Rmd")
rmarkdown::render("DrugAcronyms/06_AMR_africa_drug_acronyms.Rmd")
```

---

## Data Availability

All scripts read from `AMR_clean.csv`. Key variables:

| Variable | Description |
|----------|-------------|
| `doi` | Study identifier |
| `country` / `region` | Geographic location |
| `species` | Host species (Cattle, Chicken, Goats, Pigs, Sheep, Turkey, Environment) |
| `antimicrobial` / `antimicrobial_compound` | Drug name and abbreviation |
| `antibiotic_class` | Antibiotic class |
| `who_classification` | WHO importance classification |
| `no_isolate` | Total isolates tested |
| `no_isolates_resistant` | Resistant isolates |
| `mdr_percentage` | Multi-drug resistance percentage |
| `sampling_end_year` | Year of sampling |
| `x_coordinate` / `y_coordinate` | Geographic coordinates |

---

## Citation

To be added upon publication acceptance.
