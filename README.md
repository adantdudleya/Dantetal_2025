# Dantetal_2026

Data and R code accompanying:

**"A classification approach to categorize urbanized landscapes in California and
trait variation in Maltese starthistle"**
Anthony E. Dant et al. *Ecological Applications* (accepted; Manuscript ID EAP25-0912.R1)

To cite this deposit, use the Author-Year, repository name, and DOI in your Open
Research Statement, and give the full citation in the manuscript References (add the
DOI once the deposit is minted).

## Files in this deposit

| File | Description |
|------|-------------|
| `README.md` | This file: file- and column-level metadata for every file, run order, and data-source citations. |
| `Dantetal_2025_Cluster_Script.R` | Landscape classification and mapping: clusters California census tracts into 9 Environmental Zones, runs the landscape PCA, and produces the EZ map, distance regression, city-composition figure, the statewide GBIF map, and the collection-site map. |
| `Dantetal_2025_Greenhouse_Script.R` | Greenhouse trait analysis: nested ANOVAs (trait ~ Environmental Zone / Location) with Tukey post-hoc, the three-way Binary vs impervious-Gradient vs EZ nested comparison, and the two side-by-side comparison figures. |
| `final_shape.shp` (+ `.shx`, `.dbf`, `.cpg`, `.sbn`, `.sbx`, `.shp.xml`) | Raw landscape input: one polygon per California census tract with the 19 environmental variables. **Derived** — see Data Access below. |
| `urban_typology2.shp` (+ `.shx`, `.dbf`, `.prj`) | Census-tract polygons with the Environmental Zone assignment (`Envrn_Z`). **Derived** from `final_shape.shp` by `Dantetal_2025_Cluster_Script.R`. |
| `Greenhouse_data_Soil_treatment.csv` | Common-garden trait measurements, soil (normal) treatment; one row per plant. |
| `Greenhouse_data_Sand_treatment.csv` | Common-garden trait measurements, sand treatment (SLA, root traits); one row per plant. |
| `Impervious_Surface_Cenmel_Population.csv` | The 17 collection sites: coordinates, region, and mean NLCD impervious surface. |
| `gbif_cenmel_cache.csv` | *Centaurea melitensis* GBIF occurrences used for the maps. Included in the deposit as a derived product; GBIF is cited as the source in Data Access / References. |

*File-naming note:* all deposit file names use only letters, numbers, and underscores
(no parentheses, accents, or asterisks), per ESA guidance.

---

## Run order

1. **`Dantetal_2025_Cluster_Script.R`** — run first: landscape classification, the EZ
   map, distance regression, city-composition figure, the statewide GBIF map, and the
   collection-site map. Writes `urban_typology2.shp`.
2. **`Dantetal_2025_Greenhouse_Script.R`** — the trait analysis and comparison figures.

Each script runs top to bottom. GBIF occurrences are read from the deposited
`gbif_cenmel_cache.csv` (or downloaded once and cached if that file is absent).

---

## Software

R, with: `tidyverse, sf, ggplot2, dplyr, RColorBrewer, dendextend, cluster,
factoextra, forcats, reshape2, corrplot, patchwork, maps, rgbif, lmerTest, emmeans,
multcomp, multcompView, tibble`.

---

## Spatial data formats

Spatial data are provided as ESRI Shapefiles. ESA's guidance notes that GeoPackage
(`.gpkg`) is preferred over Shapefiles (not a strict requirement); the shapefiles can
be converted with `sf::st_write(x, "file.gpkg")` if desired. CRS and EPSG codes for each shapefile are given in the column-level metadata section above.

---

---

## Column-level metadata

*Prepared 2026-08-20. Missing values are represented by an empty cell (read as `NA` in R)
unless noted otherwise.*

**Environmental data sources** — all accessed **April 2025**. A candidate set of 49
variables was filtered to the **19** used (variables with |r| > 0.7 removed as redundant;
see Appendix S1).
- **Climate:** WorldClim **v2.1** (worldclim.org) — bioclimatic variables BIO1–BIO19.
- **Pollution / socio-economic:** **CalEnviroScreen 4.0** (2021, CalEPA) indicator
  percentiles (0–100), stored in CalEPA's pre-existing census-tract polygons.
- **Land cover (via the MRLC Consortium):** impervious surface — USGS **2024** Annual
  NLCD Collection 1 Science Products; tree canopy — USGS / U.S. Forest Service **2024**
  NLCD Tree Canopy Cover Products; all other land-cover classes — USGS Annual Land Cover
  database. Per-tract amounts were computed with **ArcGIS Pro v3.03** (Tabulate Area) and
  converted to percent of the total pixels in each tract.

**File types.** All spatial data are **vector** (ESRI Shapefile); the deposit contains **no raster files**, so ESA's raster-specific metadata requirements do not apply. Each file is described at the file level below, followed by its complete attribute/column table.

---

### 1. `final_shape.shp` — raw landscape input (census-tract polygons)
Input to `Dantetal_2025_Cluster_Script.R`. One record per California census tract.

- **Format:** ESRI Shapefile (`.shp/.shx/.dbf/.prj/…`), polygon geometry
- **Records:** 7829 tracts
- **Coordinate reference system:** NAD 1983 California (Teale) Albers, **EPSG:3310**
  (units: meters). The shipped `.prj` is empty; the CRS is documented here.
- **Character encoding:** Latin-1 (DBF contains e.g. "Piñon Hills")
- **Source:** U.S. Census TIGER/Line census-tract polygons, with the 19 environmental attributes joined from the sources listed at the top of this document (CalEnviroScreen 4.0; NLCD via MRLC, 2024 products; WorldClim v2.1).

| Column | Definition | Units | Type |
|---|---|---|---|
| `Tract` | Census tract identifier (GEOID) | — | numeric |
| `County` | California county name | — | text |
| `ApproxLoc` | Approximate locality / place name for the tract | — | text |
| `PM2_5_P` | CalEnviroScreen indicator: PM2.5 concentration | CES 4.0 statewide percentile (0–100) | numeric |
| `DieselPM_P` | CalEnviroScreen indicator: diesel particulate matter | CES 4.0 statewide percentile (0–100) | numeric |
| `PesticideP` | CalEnviroScreen indicator: pesticide use | CES 4.0 statewide percentile (0–100) | numeric |
| `Tox_Rel_P` | CalEnviroScreen indicator: toxic releases from facilities | CES 4.0 statewide percentile (0–100) | numeric |
| `TrafficP` | CalEnviroScreen indicator: vehicular traffic density | CES 4.0 statewide percentile (0–100) | numeric |
| `CleanupP` | CalEnviroScreen indicator: cleanup sites | CES 4.0 statewide percentile (0–100) | numeric |
| `ImpWatBodP` | CalEnviroScreen indicator: impaired water bodies | CES 4.0 statewide percentile (0–100) | numeric |
| `Shape_Leng` | Polygon perimeter length | m | numeric |
| `Shape_Area` | Polygon area | m² | numeric |
| `Barrenland` | NLCD barren land within tract (raw count of 30 m cells) | count | numeric |
| `Shrubland` | NLCD shrubland within tract (raw count of 30 m cells) | count | numeric |
| `Grassland` | NLCD grassland within tract (raw count of 30 m cells) | count | numeric |
| `Total_Wetl` | NLCD woody + herbaceous wetland within tract (raw count of 30 m cells) | count | numeric |
| `Crop_Pastu` | NLCD cropland + pasture within tract (raw count of 30 m cells) | count | numeric |
| `Annual_Tem` | WorldClim BIO1 — Annual mean temperature | °C | numeric |
| `Warmest_Qu` | WorldClim BIO5 — Max temperature of warmest month | °C | numeric |
| `Precip` | WorldClim BIO12 — Annual precipitation | mm | numeric |
| `Impervious` | Impervious surface (percent of tract area, NLCD) | % | numeric |
| `Treecover` | Tree canopy cover (percent of tract area, NLCD) | % | numeric |
| `Mean_Diurn` | WorldClim BIO2 — Mean diurnal range | °C | numeric |
| `ISO` | WorldClim BIO3 — Isothermality (BIO2/BIO7 × 100) | index (%) | numeric |
| `Temp_Seaso` | WorldClim BIO4 — Temperature seasonality (SD × 100) | °C × 100 | numeric |
| `Min_Temp_C` | WorldClim BIO6 — Min temperature of coldest month | °C | numeric |
| `Temp_Annua` | WorldClim BIO7 — Temperature annual range | °C | numeric |
| `Mean_Temp_` | WorldClim BIO8 — Mean temperature of wettest quarter | °C | numeric |
| `Mean_Temp1` | WorldClim BIO9 — Mean temperature of driest quarter | °C | numeric |
| `Mean_Tem_1` | WorldClim BIO10 — Mean temperature of warmest quarter | °C | numeric |
| `Mean_Tem_2` | WorldClim BIO11 — Mean temperature of coldest quarter | °C | numeric |
| `Preciptiat` | WorldClim BIO13 — Precipitation of wettest month | mm | numeric |
| `Precipti_1` | WorldClim BIO14 — Precipitation of driest month | mm | numeric |
| `Precipti_2` | WorldClim BIO15 — Precipitation seasonality (CV) | index (%) | numeric |
| `Precipti_3` | WorldClim BIO16 — Precipitation of wettest quarter | mm | numeric |
| `Precipti_4` | WorldClim BIO17 — Precipitation of driest quarter | mm | numeric |
| `Precipti_5` | WorldClim BIO18 — Precipitation of warmest quarter | mm | numeric |
| `PREC_CQ` | WorldClim BIO19 — Precipitation of coldest quarter | mm | numeric |
| `AsthmaP` | CalEnviroScreen indicator: asthma rate | CES 4.0 statewide percentile (0–100) | numeric |
| `PovertyP` | CalEnviroScreen indicator: poverty | CES 4.0 statewide percentile (0–100) | numeric |
| `Urban_core` | Distance from tract centroid to nearest urban core | m | numeric |
| `GWThreatP` | CalEnviroScreen indicator: groundwater threats | CES 4.0 statewide percentile (0–100) | numeric |

---

### 2. `urban_typology2.shp` — classified landscape (census-tract polygons with EZ)
Written by `Dantetal_2025_Cluster_Script.R`; read by the map sections and by the
greenhouse collection-site map. One record per census tract.

- **Format:** ESRI Shapefile, polygon geometry
- **Records:** 7829 tracts
- **Coordinate reference system:** NAD 1983 California (Teale) Albers, **EPSG:3310** (units: meters)
- **Note:** the `.cpg` sidecar declares UTF-8 but the DBF bytes are Latin-1; read with Latin-1.
- **Character encoding:** Latin-1
- **Source:** Derived from `final_shape.shp` by `Dantetal_2025_Cluster_Script.R` (variable standardization, hierarchical clustering of 19 environmental variables, and Environmental Zone assignment).
- Land-cover appears in two forms: raw NLCD cell counts (`Brrnlnd`, `Shrblnd`,
  `Grsslnd`, `WdyanHW`, `CrplnaP`) and the percent-of-area versions used in analysis
  (`Brrn(%)`, `Grss(%)`, `Shrb(%)`, `WaHW(%)`, `CraP(%)`).

| Column | Definition | Units | Type |
|---|---|---|---|
| `Tract` | Census tract identifier (GEOID) | — | numeric |
| `County` | California county name | — | text |
| `ApprxLc` | Approximate locality / place name for the tract | — | text |
| `PrtM2_5` | CalEnviroScreen indicator: PM2.5 concentration | CES 4.0 statewide percentile (0–100) | numeric |
| `DslPrtM` | CalEnviroScreen indicator: diesel particulate matter | CES 4.0 statewide percentile (0–100) | numeric |
| `PstcdUs` | CalEnviroScreen indicator: pesticide use | CES 4.0 statewide percentile (0–100) | numeric |
| `ToxcRls` | CalEnviroScreen indicator: toxic releases from facilities | CES 4.0 statewide percentile (0–100) | numeric |
| `VhclrTr` | CalEnviroScreen indicator: vehicular traffic density | CES 4.0 statewide percentile (0–100) | numeric |
| `ClnpSts` | CalEnviroScreen indicator: cleanup sites | CES 4.0 statewide percentile (0–100) | numeric |
| `ImprdWB` | CalEnviroScreen indicator: impaired water bodies | CES 4.0 statewide percentile (0–100) | numeric |
| `Shp_Lng` | Polygon perimeter length | m | numeric |
| `Shap_Ar` | Polygon area | m² | numeric |
| `Brrnlnd` | NLCD barren land within tract (raw count of 30 m cells) | count | numeric |
| `Shrblnd` | NLCD shrubland within tract (raw count of 30 m cells) | count | numeric |
| `Grsslnd` | NLCD grassland within tract (raw count of 30 m cells) | count | numeric |
| `WdyanHW` | NLCD woody + herbaceous wetland within tract (raw count of 30 m cells) | count | numeric |
| `CrplnaP` | NLCD cropland + pasture within tract (raw count of 30 m cells) | count | numeric |
| `Bio01` | WorldClim BIO1 — Annual mean temperature | °C | numeric |
| `Bio05` | WorldClim BIO5 — Max temperature of warmest month | °C | numeric |
| `Bio12` | WorldClim BIO12 — Annual precipitation | mm | numeric |
| `ImpS(%)` | Impervious surface (percent of tract area, NLCD) | % | numeric |
| `Trcv(%)` | Tree canopy cover (percent of tract area, NLCD) | % | numeric |
| `Bio02` | WorldClim BIO2 — Mean diurnal range (mean of monthly max−min) | °C | numeric |
| `Bio03` | WorldClim BIO3 — Isothermality (BIO2/BIO7 × 100) | index (%) | numeric |
| `Bio04` | WorldClim BIO4 — Temperature seasonality (SD × 100) | °C × 100 | numeric |
| `Bio06` | WorldClim BIO6 — Min temperature of coldest month | °C | numeric |
| `Bio07` | WorldClim BIO7 — Temperature annual range (BIO5−BIO6) | °C | numeric |
| `Bio08` | WorldClim BIO8 — Mean temperature of wettest quarter | °C | numeric |
| `Bio09` | WorldClim BIO9 — Mean temperature of driest quarter | °C | numeric |
| `Bio10` | WorldClim BIO10 — Mean temperature of warmest quarter | °C | numeric |
| `Bio11` | WorldClim BIO11 — Mean temperature of coldest quarter | °C | numeric |
| `Bio13` | WorldClim BIO13 — Precipitation of wettest month | mm | numeric |
| `Bio14` | WorldClim BIO14 — Precipitation of driest month | mm | numeric |
| `Bio15` | WorldClim BIO15 — Precipitation seasonality (CV) | index (%) | numeric |
| `Bio16` | WorldClim BIO16 — Precipitation of wettest quarter | mm | numeric |
| `Bio17` | WorldClim BIO17 — Precipitation of driest quarter | mm | numeric |
| `Bio18` | WorldClim BIO18 — Precipitation of warmest quarter | mm | numeric |
| `Bio19` | WorldClim BIO19 — Precipitation of coldest quarter | mm | numeric |
| `AsthmaP` | CalEnviroScreen indicator: asthma rate | CES 4.0 statewide percentile (0–100) | numeric |
| `PvrtyIn` | CalEnviroScreen indicator: poverty | CES 4.0 statewide percentile (0–100) | numeric |
| `DstntUC` | Distance from tract centroid to nearest urban core | m | numeric |
| `GrndwtT` | CalEnviroScreen indicator: groundwater threats | CES 4.0 statewide percentile (0–100) | numeric |
| `Brrn(%)` | Barren land cover (percent of tract area) | % | numeric |
| `Grss(%)` | Grassland cover (percent of tract area) | % | numeric |
| `Shrb(%)` | Shrubland cover (percent of tract area) | % | numeric |
| `WaHW(%)` | Woody + herbaceous wetland cover (percent of tract area) | % | numeric |
| `CraP(%)` | Cropland + pasture cover (percent of tract area) | % | numeric |
| `Envrn_Z` | Environmental Zone assignment (1–9) from hierarchical clustering | categorical (1–9) | text |

---

### 3. `Greenhouse_data_Soil_treatment.csv` — common-garden traits (soil treatment)
- **Format:** CSV, 535 rows × 12 columns (one row per plant)
- **`Site`** (join key to the impervious-surface table) has been added, derived
  from `Individual` as the text before the first digit (reproduces the 17 site
  codes; "Tuna" = Santa Cruz Island, which the analysis excludes).
- **Source:** Common-garden greenhouse measurements of *C. melitensis* grown from field-collected seed.

| Column | Definition | Units | Type |
|---|---|---|---|
| `Individual` | Plant identifier; site code followed by a within-site number (e.g. "Dark 10s") | — | text |
| `Site` | Collection-site code; join key to the impervious-surface table | — | text |
| `Location` | Collection locality (full site name) | — | text |
| `Land Classification` | Binary urban class of the source population: Urban or NonUrban | — | text |
| `Soil` | Greenhouse soil treatment (Normal or Sand); Ch.1 analyses use Normal only | — | text |
| `Number of flowerheads` | Number of flowerheads produced by the plant | count | integer |
| `Biomass` | Total dry biomass of the plant | g | numeric |
| `# Number of leaves` | Number of leaves | count | integer |
| `Height of the largest leaf (mm)` | Length of the longest leaf | mm | numeric |
| `Width of largest leaf (mm)` | Width of the widest leaf | mm | numeric |
| `Environmental_Zone` | Environmental Zone (1–9) of the source population | categorical (1–9) | integer |
| `Month Picked` | Month the source seed/plant was collected in the field | — | text |

---

### 4. `Greenhouse_data_Sand_treatment.csv` — common-garden traits (sand treatment)
- **Format:** CSV, 222 rows × 15 columns (one row per plant)
- **`Site`** (join key) has been added, derived from `Individual Name` (same rule as the soil table).
- **Source:** Common-garden greenhouse measurements of *C. melitensis* grown from field-collected seed.

| Column | Definition | Units | Type |
|---|---|---|---|
| `Individual Name` | Plant identifier; site code followed by a within-site number | — | text |
| `Site` | Collection-site code; join key to the impervious-surface table | — | text |
| `Above Ground Biomass` | Above-ground dry biomass | g | numeric |
| `Below Ground Biomass` | Below-ground (root) dry biomass | g | numeric |
| `Root_to_shoot_Ratio` | Ratio of below-ground to above-ground biomass | ratio (unitless) | numeric |
| `Area` | Leaf area used for SLA | cm² | numeric |
| `Wet Mass (g)` | Fresh (wet) mass of the leaf sample | g | numeric |
| `Dry mass` | Dry mass of the leaf sample | g | numeric |
| `Treatment` | Binary urban class of the source population: Urban or Non-Urban | — | text |
| `Environmental_Zone` | Environmental Zone (1–9) of the source population | categorical (1–9) | integer |
| `Root Diameter` | Mean root diameter | cm | numeric |
| `SLA` | Specific leaf area (leaf area / dry mass) | cm²/g | numeric |
| `SLA Date` | Date SLA was measured | date (M/D/YYYY) | date |
| `Collection Month` | Month the source seed/plant was collected in the field | — | text |
| `Location` | Collection locality (full site name) | — | text |

---

### 5. `Impervious_Surface_Cenmel_Population.csv` — site impervious surface & coordinates
Replaces the former `Master_Data_6_30_2025.csv` for the revision analyses.
- **Format:** CSV, 17 rows × 5 columns (one row per collection site)
- **Source:** Site coordinates and mean NLCD impervious surface for the 17 collection sites.

| Column | Definition | Units | Type |
|---|---|---|---|
| `Site` | Collection-site code; JOIN KEY to the greenhouse trait tables | — | text |
| `Latitude` | Site latitude (WGS84) | decimal degrees | numeric |
| `Longitude` | Site longitude (WGS84) | decimal degrees | numeric |
| `Location` | Geographic region of California (Northern / Central / Southern); note trailing spaces in some values | — | text |
| `Impervious_Surface` | Mean impervious surface at the collection site (NLCD) | % | numeric |

### 6. `gbif_cenmel_cache.csv` — GBIF occurrences for the maps
*Centaurea melitensis* occurrence records from GBIF (United States), filtered to
coordinate uncertainty ≤ 1000 m; used for the statewide GBIF map and the
collection-site map inset. **Derived** — see Data Access. Columns are as written
by the scripts (`lon`, `lat`, `unc`).

- **Format:** CSV, one row per occurrence record
- **Source:** GBIF, retrieved via the query in Data Access below; cited in the References.

| Column | Definition | Units | Type |
|---|---|---|---|
| `lon` | Occurrence longitude (WGS84) | decimal degrees | numeric |
| `lat` | Occurrence latitude (WGS84) | decimal degrees | numeric |
| `unc` | Coordinate uncertainty reported by GBIF; blank/`NA` if not reported | m | numeric |

---

## Data Access & References (derived and external data sources)

The landscape data were assembled from four external sources (all accessed **April
2025**); a candidate set of 49 variables was reduced to the 19 used, removing
variables correlated at |r| > 0.7 (see Appendix S1). Transformations:

- **`final_shape.shp`** — census-tract polygons carrying: CalEnviroScreen 4.0
  pollution and socio-economic indicator percentiles (in CalEPA's pre-existing tract
  polygons); NLCD land-cover, impervious-surface, and tree-canopy layers summarized
  per tract with ArcGIS Pro v3.03 (Tabulate Area) and converted to percent of tract
  pixels; and WorldClim v2.1 bioclimatic variables (BIO1–BIO19) extracted per tract.
- **`urban_typology2.shp`** — produced from `final_shape.shp` by the cluster script:
  the 19 variables were standardized, hierarchically clustered, and each tract was
  assigned to one of 9 Environmental Zones (`Envrn_Z`).
- **`gbif_cenmel_cache.csv`** — GBIF occurrences of *Centaurea melitensis*, retrieved
  with `rgbif::occ_search` using the query: `scientificName = "Centaurea melitensis"`,
  `country = "US"`, `hasCoordinate = TRUE`, `hasGeospatialIssue = FALSE`,
  `year = "1990,2026"`, coordinate uncertainty ≤ 1000 m. These records are deposited here as `gbif_cenmel_cache.csv`, and GBIF is cited as
  the source in the References — a GBIF download DOI is optional (GBIF recommends one
  for attribution) and not required when the data are deposited.

**References**

- OEHHA (California Office of Environmental Health Hazard Assessment). 2021.
  CalEnviroScreen 4.0. California Environmental Protection Agency, Sacramento, CA.
  https://oehha.ca.gov/calenviroscreen/report/calenviroscreen-40
- U.S. Geological Survey / Multi-Resolution Land Characteristics (MRLC) Consortium.
  2024. Annual National Land Cover Database (NLCD) Collection 1 Science Products
  (land cover, impervious surface). https://www.mrlc.gov
- U.S. Geological Survey and U.S. Forest Service. 2024. National Land Cover Database
  Tree Canopy Cover. MRLC Consortium. https://www.mrlc.gov
- Fick, S. E., and R. J. Hijmans. 2017. WorldClim 2: new 1-km spatial resolution
  climate surfaces for global land areas. International Journal of Climatology
  37:4302–4315. https://doi.org/10.1002/joc.5086  (data v2.1, https://www.worldclim.org)
- GBIF.org. GBIF Occurrence data for *Centaurea melitensis* (United States).
  https://www.gbif.org  [accessed April 2025; query above]
