# Dant et al. (2025) — Greenhouse and Landscape Analyses

This repository contains the data and R scripts used in **Dant et al. (2025)** to evaluate phenotypic and environmental variation in *Centaurea melitensis* (Malta starthistle) across California.  
The project combines greenhouse trait experiments with large-scale environmental clustering analyses of California census tracts.

---

##  Repository Contents

| File | Description |
|------|--------------|
| **Dantetal_2025_Greenhouse_Script.R** | Performs nested ANOVAs (Environmental Zone / Location) and post-hoc Tukey tests on greenhouse trait data. Generates violin + boxplots with Compact Letter Displays (CLDs) and compares Urban vs Non-Urban treatments. |
| **Dantetal_2025_Cluster_Script.R** | Conducts hierarchical clustering (Ward.D2) and principal-component analyses (PCA) on census-tract environmental, climatic, and socio-economic data to define nine Environmental Zones. Includes mapping, distance analyses, and shapefile export. |
| **Greenhouse_data_Soil_treatment.csv** | Raw morphological trait data from the soil treatment experiment. |
| **Greenhouse_data_Sand_treatment.csv** | Raw physiological and morphological trait data from the sand treatment experiment. |
| **final_shape.shp** (+ `.shx`, `.dbf`, `.prj`, `.cpg`, `.sbn`, `.sbx`, `.xml`) | Spatial shapefile representing 7,829 California census tracts with environmental, climatic, and socio-economic variables. Used for Environmental Zone clustering. |
| **README.md** | Documentation and overview of the repository contents. |

> **Note:** The shapefile must include all associated sidecar files (`.shp`, `.shx`, `.dbf`, `.prj`, `.cpg`, `.sbn`, `.sbx`, `.xml`) to function properly.

---

##  Methods Summary

Environmental and socio-economic data for 7,829 California census tracts were standardized and clustered using Ward.D2 hierarchical clustering to define **nine Environmental Zones (EZs)**.  
Principal Component Analysis (PCA) was used to visualize the multivariate structure of environmental variation.  

Greenhouse trait data from *C. melitensis* populations were analyzed using nested ANOVAs (EZ / Location) and Urban vs Non-Urban contrasts to test how environmental heterogeneity relates to morphological and physiological variation.  

Analyses were implemented in R (≥ 4.1) using the following key packages:
`tidyverse`, `sf`, `corrplot`, `RColorBrewer`, `dendextend`, `cluster`, `factoextra`, `forcats`, `ggplot2`, `reshape2`, `lmerTest`, `emmeans`, `multcomp`, and `multcompView`.

---

## 🧾 Abbreviations

| Term | Definition |
|------|-------------|
| **SLA** | Specific Leaf Area (cm²·g⁻¹) |
| **EZ** | Environmental Zone (1–9) |
| **ANOVA** | Analysis of Variance |
| **PCA** | Principal Component Analysis |
| **CLD** | Compact Letter Display (post-hoc grouping letters) |

---

##  Citation

If you use this repository, please cite:  
**Dant, A. et al. (2025).** *Applying a classification approach to categorizing urbanized landscapes in California and invasion by Maltese starthistle.* (Manuscript in preparation).

---

## 🧭 Contact

For questions or collaboration inquiries, please contact:  
**Anthony Dant** — [adant@arizona.edu](mailto:adant@arizona.edu)
