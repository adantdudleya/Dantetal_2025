############################################################
# LANDSCAPE CLUSTERING AND PCA ANALYSIS OF CALIFORNIA TRACTS
# ----------------------------------------------------------
# This script performs an integrated landscape analysis:
# - Loads environmental and socio-economic data from census tracts
# - Renames variables for interpretability
# - Scales and correlates variables
# - Performs hierarchical clustering and PCA
# - Visualizes multi-dimensional relationships among environmental factors
# - Assigns Environmental Zones and maps spatial distributions
# - Exports final shapefiles and summary tables for downstream use
############################################################


#############################################
## LOAD LIBRARIES
#############################################
# These packages handle data manipulation, visualization, spatial analysis,
# clustering, and correlation visualization.
library(tidyverse)   # Data wrangling and visualization
library(sf)          # Spatial data handling (simple features)
library(corrplot)    # Correlation matrix visualization
library(RColorBrewer)# Color palettes
library(dendextend)  # Customizing dendrograms
library(cluster)     # Clustering algorithms
library(factoextra)  # PCA and clustering visualization
library(forcats)     # Factor manipulation
library(ggplot2)     # Core plotting
library(reshape2)    # Data reshaping for long-format analysis


#############################################
## UNIFIED COLOR PALETTE & THEME
#############################################
# Define a consistent visual theme and color scheme for all plots.
# Colors are based on RColorBrewer’s "Paired" palette (9 groups = 9 clusters).
my_colors <- brewer.pal(9, "Paired")
names(my_colors) <- as.character(1:9)

# Establish a custom minimalist theme with uniform text and grid formatting.
my_theme <- theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 18, face = "bold", color = "black"),
    axis.text = element_text(size = 15, face = "bold", color = "black"),
    legend.title = element_text(size = 16, face = "bold", color = "black"),
    legend.text = element_text(size = 14, face = "bold", color = "black"),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank()
  )


#############################################################
## READ DATA: Shapefile of 7,829 Census tracts in California
#############################################################
# Import the spatial shapefile containing environmental and socio-economic attributes.
calshapes_2 <- read_sf("final_shape.shp")


#############################################
## RENAME ALL COLUMNS AT ONCE
#############################################
# Create a dictionary mapping abbreviated column names to descriptive labels.
# This step standardizes variable naming for readability and consistency across analyses.
rename_cols <- c(
  "PM2_5_P"     = "Particulate Matter 2.5",
  "DieselPM_P"  = "Diesel Particulate Matter",
  "PesticideP"  = "Pesticide Use",
  "Tox_Rel_P"   = "Toxic Releases from Facilities",
  "TrafficP"    = "Vehicular Traffic",
  "CleanupP"    = "Cleanup Sites",
  "GWThreatP"   = "Groundwater Threats",
  "ImpWatBodP"  = "Impaired Water Bodies",
  "Total_Wetl"  = "Woody and Herbaceous Wetland",
  "Crop_Pastu"  = "Cropland and Pastureland",
  "Annual_Tem"  = "Annual Mean Temperature",
  "Warmest_Qu"  = "Max Temp of the Warmest Month",
  "Precip"      = "Annual Precipitation",
  "Impervious"  = "Impervious Surface",
  "Mean_Diurn"  = "Mean Diurnal Range",
  "ISO"         = "Isothermality",
  "Temp_Seaso"  = "Temperature Seasonality",
  "Min_Temp_C"  = "Min Temp of the Coldest Month",
  "Temp_Annua"  = "Temperature Annual Range",
  "Mean_Temp_"  = "Mean Temperature of Wettest Quarter",
  "Mean_Temp1"  = "Mean Temperature of Driest Quarter",
  "Mean_Tem_1"  = "Mean Temperature of Warmest Quarter",
  "Mean_Tem_2"  = "Mean Temperature of Coldest Quarter",
  "Preciptiat"  = "Precipitation of Wettest Month",
  "Precipti_1"  = "Precipitation of Driest Month",
  "Precipti_2"  = "Precipitation Seasonality",
  "Precipti_3"  = "Precipitation of Wettest Quarter",
  "Precipti_4"  = "Precipitation of Driest Quarter",
  "Precipti_5"  = "Precipitation of Warmest Quarter",
  "PREC_CQ"     = "Precipitation of Coldest Quarter",
  "Treecover"   = "Tree Cover",
  "Urban_core"  = "Distance to Urban Core",
  "AsthmaP"     = "Asthma Rate",
  "Grassland_Percent" = "Grassland",
  "Shrubland"   = "Shrub",
  "Woody and Herbaceous Wetlands_Percent" = "Woody and Herbaceous Wetland",
  "Cropland and Pastureland_Percent"      = "Cropland and Pastureland",
  "PovertyP"    = "Poverty Index"
)

# Apply renaming and clean up whitespace in field names.
names(calshapes_2)[names(calshapes_2) %in% names(rename_cols)] <-
  rename_cols[names(calshapes_2)[names(calshapes_2) %in% names(rename_cols)]]
names(calshapes_2) <- trimws(names(calshapes_2), which = "both")


#############################################
## CALCULATE PERCENTAGES FOR LAND COVER
#############################################
# Convert absolute area values into percentages relative to each tract’s total area.
land_cover_cols <- c(
  "Barrenland", "Grassland", "Shrub", 
  "Woody and Herbaceous Wetland", "Cropland and Pastureland"
)

for (col in land_cover_cols) {
  if (col %in% names(calshapes_2)) {
    calshapes_2[[col]] <- (calshapes_2[[col]] / calshapes_2$Shape_Area) * 100
  }
}


#############################################
## SELECT FINAL COLUMNS FOR CLUSTERING
#############################################
# Select relevant variables spanning pollution, socio-economic,
# climatic, and land-cover categories for clustering and PCA.
landscape_vars <- c(
  "Tract",
  "Particulate Matter 2.5",
  "Pesticide Use",
  "Vehicular Traffic",
  "Diesel Particulate Matter",
  "Toxic Releases from Facilities",
  "Cleanup Sites",
  "Groundwater Threats",
  "Impaired Water Bodies",
  "Poverty Index",
  "Min Temp of the Coldest Month",
  "Annual Precipitation",
  "Max Temp of the Warmest Month",
  "Impervious Surface",
  "Tree Cover",
  "Grassland",
  "Shrub",
  "Woody and Herbaceous Wetland",
  "Cropland and Pastureland",
  "Asthma Rate"
)

# Create a clean subset containing only relevant fields and geometry.
calshapes_new <- calshapes_2 %>%
  select(all_of(landscape_vars), geometry)


#############################################
## VALIDATE STRUCTURE
#############################################
# Print summary structure and check for missing values before clustering.
message("✔ Column count: ", ncol(calshapes_new))  # Expected 21


if (anyNA(calshapes_new)) {
  message("⚠ There ARE missing values!")
  print(colSums(is.na(calshapes_new)))
} else {
  message("✅ No missing values in calshapes_new.")
}


#############################################
## CORRELATION PLOT
#############################################
# Compute and visualize variable correlations to detect redundancy or covariance.
scaled_dat <- calshapes_new %>%
  select(-Tract) %>%
  st_set_geometry(NULL)

M <- cor(scaled_dat, use = "pairwise.complete.obs")

corrplot(M,
         order = "hclust",
         type = "upper",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         tl.col = "black",
         tl.cex = 0.8)


#############################################
## FIND HIGH CORRELATIONS
#############################################
# Identify highly correlated variable pairs (|r| ≥ 0.5) for interpretation or pruning.
M_long <- melt(round(M, 3), varnames = c("Var1", "Var2"), value.name = "Correlation") %>%
  filter(Var1 != Var2) %>%
  rowwise() %>%
  mutate(Pair = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = "___")) %>%
  ungroup() %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-Pair) %>%
  filter(abs(Correlation) >= 0.5) %>%
  arrange(desc(abs(Correlation)))

View(M_long)


#############################################
## HIERARCHICAL CLUSTERING
#############################################
# Perform Ward.D2 hierarchical clustering based on scaled Euclidean distances
# to group census tracts with similar environmental profiles.
dend <- calshapes_new %>%
  select(-Tract) %>%
  st_set_geometry(NULL) %>%
  scale() %>%
  dist(method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram()

# Color dendrogram branches to visualize 9 clusters.
dend_colored <- dend %>%
  set("labels", NA) %>%
  color_branches(k = 9, col= my_colors) %>%
  set("branches_lwd", 0.5)

plot(dend_colored, main = "Colored Dendrogram")


#############################################
## CLUSTER ASSIGNMENT
#############################################
# Assign each census tract to one of 9 Environmental Zones based on dendrogram cut.
k <- 9
cluster <- cutree(dend, k = k)

calshapes_new <- calshapes_new %>%
  mutate(Environmental_Zone = as.factor(cluster))

calshapes_2 <- calshapes_2 %>%
  mutate(Environmental_Zone = as.factor(cluster))


#############################################
## PCA ANALYSIS
#############################################
# Conduct Principal Component Analysis (PCA) to visualize dominant
# environmental gradients and variable loadings across tracts.

# Define subset of relevant landscape variables (same as used in clustering)
landscape_vars <- c(
  "Particulate Matter 2.5", "Pesticide Use", "Vehicular Traffic",
  "Diesel Particulate Matter", "Toxic Releases from Facilities",
  "Cleanup Sites", "Groundwater Threats", "Impaired Water Bodies",
  "Poverty Index", "Min Temp of the Coldest Month", "Annual Precipitation",
  "Max Temp of the Warmest Month", "Impervious Surface", "Tree Cover",
  "Grassland", "Shrub", "Woody and Herbaceous Wetland",
  "Cropland and Pastureland", "Asthma Rate"
)

# Drop geometry, scale, and compute PCA
landscape_dat <- calshapes_new %>%
  st_drop_geometry() %>%
  select(all_of(landscape_vars)) %>%
  na.omit()

landscape_pca <- prcomp(landscape_dat, center = TRUE, scale. = TRUE)

# Define variable groupings by category for color-coding
var_groups <- c(
  "Particulate Matter 2.5" = "Pollution",
  "Pesticide Use" = "Pollution",
  "Vehicular Traffic" = "Pollution",
  "Cleanup Sites" = "Pollution",
  "Groundwater Threats" = "Pollution",
  "Impaired Water Bodies" = "Pollution",
  "Toxic Releases from Facilities" = "Pollution",
  "Diesel Particulate Matter" = "Pollution",
  "Poverty Index" = "Socio-Economic",
  "Asthma Rate" = "Socio-Economic",
  "Min Temp of the Coldest Month" = "Climatic",
  "Max Temp of the Warmest Month" = "Climatic",
  "Annual Precipitation" = "Climatic",
  "Impervious Surface" = "Landcover",
  "Tree Cover" = "Landcover",
  "Grassland" = "Landcover",
  "Shrub" = "Landcover",
  "Woody and Herbaceous Wetland" = "Landcover",
  "Cropland and Pastureland" = "Landcover"
)

# Map colors by variable group for PCA visualization
var_names <- rownames(landscape_pca$rotation)
col_vector <- var_groups[var_names]

# Plot PCA (PC1 vs PC2 and PC3 vs PC4)
fviz_pca_var(
  landscape_pca,
  col.var = col_vector,
  palette = "Set1",
  labelsize = 5,
  repel = TRUE,
  arrowsize = 1,
  title = "Landscape PCA: PC1 vs PC2"
) + my_theme

fviz_pca_var(
  landscape_pca,
  axes = c(3, 4),
  col.var = col_vector,
  palette = "Set1",
  labelsize = 5,
  repel = TRUE,
  arrowsize = 1,
  title = "Landscape PCA: PC3 vs PC4"
) + my_theme


#############################################
## PCA OF BIOCLIM VARIABLES
#############################################
# Conduct separate PCA using bioclimatic (temperature & precipitation) variables.
bioclim_vars <- c(
  "Annual Mean Temperature", "Mean Diurnal Range", "Isothermality",
  "Temperature Seasonality", "Max Temp of the Warmest Month",
  "Min Temp of the Coldest Month", "Temperature Annual Range",
  "Mean Temperature of Wettest Quarter", "Mean Temperature of Driest Quarter",
  "Mean Temperature of Warmest Quarter", "Mean Temperature of Coldest Quarter",
  "Annual Precipitation", "Precipitation of Wettest Month",
  "Precipitation of Driest Month", "Precipitation Seasonality",
  "Precipitation of Wettest Quarter", "Precipitation of Driest Quarter",
  "Precipitation of Warmest Quarter", "Precipitation of Coldest Quarter"
)

bioclim_data <- calshapes_2 %>%
  st_drop_geometry() %>%
  select(all_of(bioclim_vars)) %>%
  na.omit()

bioclim_pca <- prcomp(bioclim_data, center = TRUE, scale. = TRUE)

# Assign grouping: first 11 temperature, next 8 precipitation variables
bioclim_groups <- c(rep("Temperature", 11), rep("Precipitation", 8))

bioclim_var_names <- rownames(bioclim_pca$rotation)

fviz_pca_var(
  bioclim_pca,
  col.var = bioclim_groups,
  palette = "Set1",
  labelsize = 5,
  repel = TRUE,
  arrowsize = 1,
  title = "Bioclim PCA: PC1 vs PC2"
) + my_theme


#############################################
## CORRELATION PLOT - BIOCLIM VARIABLES
#############################################
# Visualize intercorrelations among bioclimatic variables.
bioclim_cor <- cor(bioclim_data, use = "pairwise.complete.obs")
bioclim_cor_rounded <- round(bioclim_cor, 3)

corrplot(
  bioclim_cor,
  order = "hclust",
  type = "upper",
  col = brewer.pal(n = 10, name = "RdYlBu"),
  tl.col = "black",
  tl.cex = 0.8,
  title = "Correlation Plot - Bioclim Variables",
  mar = c(0,0,2,0)
)


#############################################
## FOREST PLOT
#############################################

# Handle shifted negative temperature for plotting
min_val <- min(calshapes_new$`Min Temp of the Coldest Month`, na.rm = TRUE)
shift_amt <- abs(min_val) + 1
calshapes_new <- calshapes_new %>%
  mutate(`Min Temp of the Coldest Month` = 
           `Min Temp of the Coldest Month` + shift_amt)

cat("Shifted Min Temperature by ", shift_amt, "\n")

df_long <- calshapes_new %>%
  st_drop_geometry() %>%
  select(Environmental_Zone, all_of(landscape_vars)) %>%
  mutate(across(all_of(landscape_vars), log1p)) %>%
  pivot_longer(-Environmental_Zone, names_to = "Variable", values_to = "LogValue")

overall_means <- df_long %>%
  group_by(Variable) %>%
  summarise(OverallMean = mean(LogValue, na.rm = TRUE), .groups = "drop")

zone_means <- df_long %>%
  group_by(Environmental_Zone, Variable) %>%
  summarise(ZoneMean = mean(LogValue, na.rm = TRUE), .groups = "drop") %>%
  left_join(overall_means, by = "Variable") %>%
  mutate(Diff = ZoneMean - OverallMean)

ref_order <- zone_means %>%
  filter(Environmental_Zone == "1") %>%
  arrange(Diff) %>%
  pull(Variable)


# Ensure proper ordering of fill groups
zone_means <- zone_means %>%
  mutate(Environmental_Zone = factor(as.character(Environmental_Zone), levels = as.character(1:9)))

# Build plot once, after setting factor
forest_plot <- ggplot(
  zone_means,
  aes(
    x = Diff,
    y = fct_relevel(Variable, ref_order),
    fill = Environmental_Zone,
    group = Environmental_Zone
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = my_colors) +
  labs(
    title = " ",
    x = "Difference from Overall Mean (log scale)",
    y = NULL,
    fill = "Environmental Zone"
  ) +
  my_theme



print(forest_plot)




zone_means <- zone_means %>%
  mutate(
    Zone_Reversed = as.character(10 - as.numeric(as.character(Environmental_Zone)))  # 1→9, 2→8, ..., 9→1
  )


# Paired color palette reversed
my_colors_reversed <- rev(brewer.pal(9, "Paired"))
names(my_colors_reversed) <- as.character(1:9)  # Match flipped zone levels


forest_plot <- ggplot(
  zone_means,
  aes(
    x = Diff,
    y = fct_relevel(Variable, ref_order),
    fill = Zone_Reversed,
    group = Zone_Reversed
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),  # Your preferred dodge width
    width = 1  # Your preferred bar width
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = my_colors_reversed,
    labels = as.character(1:9),  # Keep legend correct
    name = "Environmental Zone"
  ) +
  labs(
    title = "Forest Plot: Differences from Overall Mean by Environmental Zone",
    x = "Difference from Overall Mean (log scale)",
    y = NULL
  ) +
  my_theme

print(forest_plot)



#############################################
## DISTANCE ANALYSIS
#############################################
# Purpose:
# This section quantifies how Environmental Zones are distributed 
# along a distance gradient from urban cores. It converts the 
# "Distance to Urban Core" variable into kilometers, bins it 
# into 30 km intervals, and calculates the relative proportion 
# of each Environmental Zone within those distance bins.

# Convert distance from meters to kilometers and bin into 30 km intervals
calshapes_new <- calshapes_2 %>%
  mutate(
    Distance_km = `Distance to Urban Core` / 1000,
    distance_bin = cut(Distance_km, breaks = seq(0, 210, by = 30),
                       include.lowest = TRUE, right = FALSE)
  )

# Prepare non-spatial data for summarization
cal_distance <- calshapes_new %>%
  st_drop_geometry() %>%
  mutate(Environmental_Zone = as.factor(cluster))

# Compute percentage of tracts per Environmental Zone within each distance bin
zone_dist_pct <- cal_distance %>%
  group_by(distance_bin, Environmental_Zone) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(distance_bin) %>%
  mutate(percent = 100 * count / sum(count))

# Plot the relative distribution of Environmental Zones by distance
ggplot(zone_dist_pct, aes(x = distance_bin, y = percent, fill = Environmental_Zone)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.9) +
  scale_fill_manual(values = my_colors) +
  labs(
    x = "Distance from Urban Center (km)",
    y = "Percentage of Tracts",
    fill = "Environmental Zone"
  ) +
  my_theme


#############################################
## EXPORT SHAPEFILE
#############################################
# Purpose:
# Export the final spatial dataset, including Environmental Zone labels 
# and standardized BIOCLIM variable names (Bio01–Bio19) for compatibility 
# with ecological modeling tools such as Maxent or ResistanceGA.

# Define renaming scheme for BIOCLIM variable codes
rename_cols_shp <- c(
  "Annual Mean Temperature"               = "Bio01",
  "Mean Diurnal Range"                    = "Bio02",
  "Isothermality"                         = "Bio03",
  "Temperature Seasonality"               = "Bio04",
  "Max Temperature of Warmest Month"      = "Bio05",
  "Min Temperature of Coldest Month"      = "Bio06",
  "Temperature Annual Range"              = "Bio07",
  "Mean Temperature of Wettest Quarter"   = "Bio08",
  "Mean Temperature of Driest Quarter"    = "Bio09",
  "Mean Temperature of Warmest Quarter"   = "Bio10",
  "Mean Temperature of Coldest Quarter"   = "Bio11",
  "Annual Precipitation"                  = "Bio12",
  "Precipitation of Wettest Month"        = "Bio13",
  "Precipitation of Driest Month"         = "Bio14",
  "Precipitation Seasonality"             = "Bio15",
  "Precipitation of Wettest Quarter"      = "Bio16",
  "Precipitation of Driest Quarter"       = "Bio17",
  "Precipitation of Warmest Quarter"      = "Bio18",
  "Precipitation of Coldest Quarter"      = "Bio19"
)

# Apply new names to shapefile fields
names(calshapes_2)[names(calshapes_2) %in% names(rename_cols_shp)] <-
  rename_cols_shp[names(calshapes_2)[names(calshapes_2) %in% names(rename_cols_shp)]]

# Export to ESRI Shapefile format
st_write(calshapes_2, "urban_typology.shp", driver = "ESRI Shapefile")


#############################################
## EXPORT DISTANCE CSV
#############################################
# Purpose:
# Save the Environmental Zone proportions by distance bin 
# as a CSV file for supplementary analysis or visualization.
write.csv(zone_dist_pct, "environmental_zone_distribution_by_distance.csv", row.names = FALSE)


#############################################
## VISUALIZATION BY REGION
#############################################
# Purpose:
# Generate separate regional maps of Environmental Zones 
# for California and its major metropolitan regions. 
# Each map uses consistent color palettes and legend formatting.

# Define regions of interest (entire state + 3 metro areas)
regions <- list(
  California = calshapes_2,
  Los_Angeles = calshapes_2 %>% filter(County == "Los Angeles"),
  San_Diego = calshapes_2 %>% filter(County == "San Diego"),
  SF_Bay_Area = calshapes_2 %>% filter(County %in% c(
    "San Francisco", "Alameda", "Contra Costa", "San Mateo", "Marin", "Santa Clara"
  ))
)

# Iterate through each region to visualize and summarize Environmental Zones
for (region_name in names(regions)) {
  region <- regions[[region_name]]
  
  if (nrow(region) > 0) {
    
    cat(paste0("\nTotal census tracts in ", region_name, ": ", nrow(region), "\n"))
    
    # Add black border for subregions (none for full state map)
    border_color <- ifelse(region_name == "California", NA, "black")
    
    # Map Environmental Zones
    p <- ggplot(region) +
      geom_sf(aes(fill = Environmental_Zone, alpha = 0.1), color = border_color) +
      scale_fill_manual(values = my_colors) +
      labs(
        title = paste("Environmental Zones -", gsub("_", " ", region_name)),
        fill = "Environmental Zone",
        x = "",
        y = ""
      ) +
      coord_sf() +
      theme_classic() +
      theme(
        axis.text.x = element_text(size = 18, color = "black", face = "bold"),
        axis.text.y = element_text(size = 18, color = "black", face = "bold"),
        axis.title.x = element_text(size = 25, color = "black", face = "bold"),
        axis.title.y = element_text(size = 25, color = "black", face = "bold"),
        plot.title = element_text(size = 25, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 25, color = "black"),
        legend.title = element_text(size = 25, face = "bold", color = "black"),
        legend.text = element_text(size = 25, face = "bold", color = "black")
      )
    
    # Display regional map
    print(p)
    
    # Print summary table of tracts per Environmental Zone
    print(
      region %>%
        st_drop_geometry() %>%
        count(Environmental_Zone, sort = TRUE) %>%
        rename(`Environmental Zone` = Environmental_Zone, Count = n)
    )
  }
}

  
               