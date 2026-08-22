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
library(patchwork)   # Figures
library(rgbif)       # GBIF
library(maps)        #For California inset map


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
# Requires patchwork (for the combined (a)/(b) figure).
library(patchwork)

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

# variance explained, works for prcomp, princomp, and FactoMineR::PCA objects
eig <- get_eigenvalue(landscape_pca)

# ---------------------------------------------------------------------------
# Combined biplot: PC1 vs PC2 = panel (a), PC3 vs PC4 = panel (b),
# with a boxed (a)/(b) tag in the top-right corner of each panel.
# ---------------------------------------------------------------------------
make_pca_panel <- function(axes, xlab, ylab, tag) {
  fviz_pca_var(
    landscape_pca,
    axes      = axes,
    col.var   = col_vector,
    palette   = "Set1",
    labelsize = 5,
    repel     = TRUE,
    arrowsize = 1
  ) +
    labs(title = NULL, x = xlab, y = ylab) +
    # boxed (a)/(b) tag, top-right corner (nudge with hjust/vjust)
    annotate(
      "label", x = Inf, y = Inf, label = tag,
      hjust = 1.3, vjust = 1.3, size = 6,
      label.size = 0.5, label.r = unit(0, "pt"),
      fill = "white", color = "black"
    ) +
    my_theme
}

pca_a <- make_pca_panel(
  axes = c(1, 2),
  xlab = sprintf("PC1 (%.1f%%)", eig[1, "variance.percent"]),
  ylab = sprintf("PC2 (%.1f%%)", eig[2, "variance.percent"]),
  tag  = "(a)"
)

pca_b <- make_pca_panel(
  axes = c(3, 4),
  xlab = sprintf("PC3 (%.1f%%)", eig[3, "variance.percent"]),
  ylab = sprintf("PC4 (%.1f%%)", eig[4, "variance.percent"]),
  tag  = "(b)"
)

# Side by side, sharing one legend. To drop the legend entirely (as in the
# reference figure), set legend.position = "none" below.
pca_combined <- (pca_a | pca_b) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(pca_combined)

ggsave(
  filename = "landscape_pca_biplot.png",
  plot     = pca_combined,
  width    = 15,
  height   = 7,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# Vector PDF (device inferred from the .pdf extension; best for publication)
ggsave(
  filename = "landscape_pca_biplot.pdf",
  plot     = pca_combined,
  width    = 15,
  height   = 7,
  units    = "in",
  bg       = "white"
)


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



# Paired color palette reversed

my_colors_reversed <- rev(brewer.pal(9, "Paired"))
names(my_colors_reversed) <- as.character(1:9)  # Match flipped zone levels


zone_means <- zone_means %>%
  mutate(
    Zone_Reversed = as.character(10 - as.numeric(as.character(Environmental_Zone)))  # 1→9, 2→8, ..., 9→1
  )


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
st_write(calshapes_2, "urban_typology2.shp", driver = "ESRI Shapefile")


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




# ##########################################################################
#FIGURE/ANALYSIS: Distance regression + city composition figure + chi-square
# ##########################################################################


my_colors <- brewer.pal(9, "Paired")
names(my_colors) <- as.character(1:9)

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

# =============================================================================
# CITY x ENVIRONMENTAL ZONE COUNTS  (folded in so this script is self-contained)
# -----------------------------------------------------------------------------
# Previously `city_ez_counts` was built only in the exploratory revision script
# (05_Exploratory_Revision_Analyses.R), which forced you to run that first.
# It is rebuilt here from calshapes_2 alone -- identical regions, counts, and
# within-city percentages -- so Sections 2 and 3 below no longer depend on it.
# Count feeds the chi-square; Percent feeds panel (e).
# =============================================================================
city_ez_counts <- calshapes_2 %>%
  st_drop_geometry() %>%
  mutate(City = dplyr::case_when(
    County == "Los Angeles" ~ "Los Angeles",
    County == "San Diego"   ~ "San Diego",
    County %in% c("San Francisco", "Alameda", "Contra Costa",
                  "San Mateo", "Marin", "Santa Clara") ~ "SF Bay Area",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(City)) %>%
  count(City, Environmental_Zone, name = "Count") %>%
  group_by(City) %>%
  mutate(Percent = 100 * Count / sum(Count)) %>%
  ungroup() %>%
  mutate(
    City = factor(City, levels = c("Los Angeles", "San Diego", "SF Bay Area")),
    Environmental_Zone = factor(Environmental_Zone, levels = as.character(1:9))
  )

# =============================================================================
# 1. DISTANCE REGRESSION  (final version; two earlier drafts removed)
# =============================================================================
cat("\n========================================\n")
cat("DISTANCE REGRESSION PLOT\n")
cat("========================================\n")

# Prepare data
cal_distance <- calshapes_2 %>%
  st_drop_geometry() %>%
  mutate(
    Distance_km = `Distance to Urban Core` / 1000,
    Environmental_Zone = as.factor(Environmental_Zone),
    distance_bin = cut(Distance_km, breaks = seq(0, 210, by = 30),
                       include.lowest = TRUE, right = FALSE)
  )

# Calculate bin midpoints
bin_midpoints <- cal_distance %>%
  group_by(distance_bin) %>%
  summarise(bin_midpoint = median(Distance_km, na.rm = TRUE), .groups = "drop")

# Percentage of each EZ within each bin
zone_dist_pct <- cal_distance %>%
  group_by(distance_bin, Environmental_Zone) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(distance_bin) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup() %>%
  left_join(bin_midpoints, by = "distance_bin")

# Print the median calculation breakdown
cat("\n--- Median Distance per Bin ---\n")
print(
  cal_distance %>%
    group_by(distance_bin) %>%
    summarise(
      n_tracts = n(),
      min_km = round(min(Distance_km), 1),
      median_km = round(median(Distance_km), 1),
      max_km = round(max(Distance_km), 1),
      .groups = "drop"
    ),
  n = 20
)

# --- Calculate significance per EZ ---
ez_sig <- data.frame()
for (ez in unique(as.character(zone_dist_pct$Environmental_Zone))) {
  ez_data <- zone_dist_pct %>% filter(as.character(Environmental_Zone) == ez)
  if (nrow(ez_data) >= 2) {
    mod <- lm(percent ~ bin_midpoint, data = ez_data)
    s <- summary(mod)
    p_val <- if (nrow(ez_data) >= 3) s$coefficients[2, 4] else 1
    ez_sig <- bind_rows(ez_sig, data.frame(
      Environmental_Zone = ez,
      slope = coef(mod)[2],
      r2 = s$r.squared,
      p_val = p_val,
      significant = ifelse(p_val < 0.05, "Significant", "Not Significant")
    ))
  }
}

cat("\n--- Regression Statistics ---\n")
print(ez_sig)

# Add significance to the data
zone_dist_pct <- zone_dist_pct %>%
  left_join(ez_sig %>% select(Environmental_Zone, significant),
            by = "Environmental_Zone")

# --- Generate fitted line data manually ---
fitted_lines <- data.frame()
for (ez in unique(as.character(zone_dist_pct$Environmental_Zone))) {
  ez_data <- zone_dist_pct %>%
    filter(as.character(Environmental_Zone) == ez)
  if (nrow(ez_data) >= 2) {
    mod <- lm(percent ~ bin_midpoint, data = ez_data)
    pred_x <- seq(min(ez_data$bin_midpoint), max(ez_data$bin_midpoint), length.out = 50)
    pred_y <- predict(mod, newdata = data.frame(bin_midpoint = pred_x))
    sig_val <- ez_sig$significant[as.character(ez_sig$Environmental_Zone) == ez]
    if (length(sig_val) == 0) sig_val <- "Not Significant"
    fitted_lines <- bind_rows(fitted_lines, data.frame(
      bin_midpoint = pred_x,
      percent = pred_y,
      Environmental_Zone = ez,
      significant = sig_val
    ))
  }
}

cat("Fitted lines rows:", nrow(fitted_lines), "\n")

# --- Plot ---
distance_regression <- ggplot() +
  geom_line(data = fitted_lines,
            aes(x = bin_midpoint, y = percent, color = Environmental_Zone,
                linetype = significant),
            linewidth = 2.5) +
  geom_point(data = zone_dist_pct,
             aes(x = bin_midpoint, y = percent, color = Environmental_Zone),
             size = 7) +
  scale_color_manual(values = my_colors, name = "Environmental\nZone") +
  scale_linetype_manual(values = c("Significant" = "solid",
                                   "Not Significant" = "dashed"),
                        name = "Regression") +
  labs(
    x = "Distance from Urban Center (km)",
    y = "Percentage of Tracts (%)"
  ) +
  my_theme

print(distance_regression)
ggsave("Fig5_Distance_Regression.png", distance_regression,
       width = 12, height = 8, dpi = 300, bg = "white")
cat("Saved: Fig5_Distance_Regression.png\n")

# =============================================================================
# 2. COMBINED FIGURE 4: CA MAP + OUTTAKES + STACKED BAR
# =============================================================================

cat("\n\n========================================\n")
cat("COMBINED FIGURE 4\n")
cat("========================================\n")

# --- Define city bounding boxes for outtake rectangles ---

# SF Bay Area
sf_counties <- c("San Francisco", "Alameda", "Contra Costa",
                  "San Mateo", "Marin", "Santa Clara")
sf_data <- calshapes_2 %>% filter(County %in% sf_counties)
sf_bbox <- st_bbox(sf_data)

# San Diego
sd_data <- calshapes_2 %>% filter(County == "San Diego")
sd_bbox <- st_bbox(sd_data)

# Los Angeles
la_data <- calshapes_2 %>% filter(County == "Los Angeles")
la_bbox <- st_bbox(la_data)

# --- Panel (a): California with outtake boxes ---

# Create rectangle data for the three cities
outtake_boxes <- data.frame(
  city = c("SF Bay Area", "San Diego", "Los Angeles"),
  xmin = c(sf_bbox["xmin"], sd_bbox["xmin"], la_bbox["xmin"]),
  xmax = c(sf_bbox["xmax"], sd_bbox["xmax"], la_bbox["xmax"]),
  ymin = c(sf_bbox["ymin"], sd_bbox["ymin"], la_bbox["ymin"]),
  ymax = c(sf_bbox["ymax"], sd_bbox["ymax"], la_bbox["ymax"]),
  label = c("(b)", "(c)", "(d)")
)

panel_a <- ggplot() +
  geom_sf(data = calshapes_2, aes(fill = Environmental_Zone),
          color = NA, alpha = 0.7) +
  # Outtake rectangles
  geom_rect(data = outtake_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", linewidth = 1.2, linetype = "solid") +
  # Labels for outtake boxes
  geom_text(data = outtake_boxes,
            aes(x = xmax + 0.3, y = (ymin + ymax) / 2, label = label),
            fontface = "bold", size = 5, hjust = 0) +
  scale_fill_manual(values = my_colors) +
  labs(x = "", y = "", title = "(a) California") +
  coord_sf() +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 9, face = "bold", color = "black"),
    axis.line = element_line(color = "black")
  )

# --- Panel (b): SF Bay Area ---
panel_b <- ggplot() +
  geom_sf(data = sf_data, aes(fill = Environmental_Zone),
          color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(x = "", y = "", title = "(b) San Francisco") +
  coord_sf() +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 9, face = "bold", color = "black")
  )

# --- Panel (c): San Diego ---
panel_c <- ggplot() +
  geom_sf(data = sd_data, aes(fill = Environmental_Zone),
          color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  labs(x = "", y = "", title = "(c) San Diego") +
  coord_sf() +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 9, face = "bold", color = "black")
  )

# --- Fix LA axis labels: use fewer breaks ---
panel_d <- ggplot() +
  geom_sf(data = la_data, aes(fill = Environmental_Zone),
          color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(n.breaks = 4) +
  labs(x = "", y = "", title = "(d) Los Angeles") +
  coord_sf() +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 9, face = "bold", color = "black")
  )

# --- Make stacked bar smaller with padding ---
panel_e <- ggplot(city_ez_counts,
                  aes(x = City, y = Percent, fill = Environmental_Zone)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = my_colors, name = "Environmental\nZone", drop = FALSE) +
  labs(x = NULL, y = "Proportion (%)", title = "(e)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, face = "bold", color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    plot.margin = margin(10, 60, 10, 60)  # add horizontal padding to center it
  )

right_side <- (panel_b | panel_d) / panel_c / panel_e +
  plot_layout(heights = c(1, 1, 0.5))

combined_fig4 <- panel_a | right_side +
  plot_layout(widths = c(1, 1.5))



print(combined_fig4)
ggsave("Fig4_Combined_Outtakes_StackedBar.png", combined_fig4,
       width = 14, height = 16, dpi = 300, bg = "white")

# =============================================================================
# 3. CHI-SQUARE TEST (overall only)
# =============================================================================

cat("\n\n========================================\n")
cat("CHI-SQUARE TEST\n")
cat("========================================\n")

# Build contingency table
contingency_wide <- city_ez_counts %>%
  select(City, Environmental_Zone, Count) %>%
  pivot_wider(names_from = Environmental_Zone, values_from = Count, values_fill = 0)

contingency_matrix <- as.matrix(contingency_wide[, -1])
rownames(contingency_matrix) <- contingency_wide$City

cat("\nContingency table:\n")
print(contingency_matrix)

# Overall chi-square
chi_result <- chisq.test(contingency_matrix)
cat("\n--- Chi-Square Test: EZ composition across three cities ---\n")
cat("X² =", round(chi_result$statistic, 2), "\n")
cat("df =", chi_result$parameter, "\n")
cat("p  =", format(chi_result$p.value, digits = 6, scientific = TRUE), "\n")

# Standardized residuals (shows which city-EZ combos deviate most)
cat("\n--- Standardized Residuals (>|2| = notable departure) ---\n")
print(round(chi_result$stdres, 2))

# =============================================================================
# DONE
# =============================================================================
cat("\n\n========================================\n")
cat("DISTANCE + CITY FIGURES COMPLETE\n")
cat("========================================\n")
cat("Requires in memory: calshapes_2 (from 01_Landscape_Classification_PCA_Map.R)\n")
cat("Figures saved:\n")
cat("  Fig5_Distance_Regression.png\n")
cat("  Fig4_Combined_Outtakes_StackedBar.png\n")
cat("Caption note: 'Environmental Zone delineation described in Figures 2 and 3'\n")
cat("========================================\n")


# ##########################################################################
# FIGURE/ANALYSIS: Statewide GBIF occurrence map
# ##########################################################################

# ============================================================================
# LAYOUT + TEXT SIZES  (points are literal at the width below)
# ============================================================================
fig_width_in  <- 7.10   # 180 mm full-page width (use 3.15 for single column)
fig_height_in <- 7.60

MAIN_AXIS_TITLE  <- 15
MAIN_AXIS_TEXT   <- 13
LEGEND_TITLE     <- 14
LEGEND_TEXT      <- 13
INSET_AXIS_TITLE <- 15   # <- MUCH larger inset axis titles (tune to taste)
INSET_AXIS_TEXT  <- 12   # <- MUCH larger inset tick labels

pad_lon <- 0.7          # degrees of padding around the sampling extent
pad_lat <- 0.5

# Inset position within the figure (npc units, 0-1; top-right corner)
inset_left <- 0.64; inset_bottom <- 0.62; inset_right <- 1.00; inset_top <- 1.00

# Legend position on the main panel (npc units)
legend_pos <- c(0.16, 0.30)
# ============================================================================


################################################################################
# 1. COLOR PALETTE
################################################################################
full_palette <- brewer.pal(9, "Paired")
zone_colors  <- setNames(full_palette, as.character(1:9))


################################################################################
# 2. SITE COORDINATES (impervious table) + EZ / Urban-Non-Urban (greenhouse)
#    EZ and Land Classification come from the included soil treatment ("Normal");
#    the excluded "Sand" rows carried inconsistent Environmental_Zone values.
################################################################################
imp_sites <- read.csv("Impervious_Surface_Cenmel_Population.csv")

plant <- read.csv("Greenhouse_data_Soil_treatment.csv")
site_meta <- plant %>%
  filter(!is.na(Site), Site != "", tolower(Soil) == "normal") %>%
  group_by(Site) %>%
  summarise(
    EZ             = dplyr::first(Environmental_Zone),
    Classification = dplyr::first(Land.Classification),
    .groups = "drop"
  ) %>%
  mutate(
    EZ             = factor(as.character(EZ), levels = as.character(1:9)),
    Classification = gsub("Non[ ._]?Urban", "Non-Urban", Classification),
    Classification = factor(Classification, levels = c("Non-Urban", "Urban"))
  )

site_coords <- imp_sites %>%
  dplyr::select(Site, Longitude, Latitude) %>%
  inner_join(site_meta, by = "Site") %>%
  filter(!is.na(Longitude), !is.na(Latitude))

cat("Chapter 1 collection sites:", nrow(site_coords), "\n")
print(table(site_coords$Classification, useNA = "ifany"))


################################################################################
# 3. CROP EXTENT (drives the main-panel window AND the inset rectangle)
################################################################################
sites_sf  <- st_as_sf(site_coords, coords = c("Longitude", "Latitude"), crs = 4326)
site_bbox <- st_bbox(sites_sf)

xmin_crop <- as.numeric(site_bbox["xmin"]) - pad_lon
xmax_crop <- as.numeric(site_bbox["xmax"]) + pad_lon
ymin_crop <- as.numeric(site_bbox["ymin"]) - pad_lat
ymax_crop <- as.numeric(site_bbox["ymax"]) + pad_lat


################################################################################
# 4. LOAD + CLEAN + CROP THE ENVIRONMENTAL ZONE SHAPEFILE
################################################################################
ez_shp <- read_sf("urban_typology2.shp")
cat("Shapefile loaded:", nrow(ez_shp), "census tracts\n")

ez_zone_col <- grep("Envrn_Z", names(ez_shp), value = TRUE)[1]
if (is.na(ez_zone_col)) stop("Could not find Envrn_Z column in urban_typology2.shp")
ez_shp$EZ <- factor(ez_shp[[ez_zone_col]], levels = as.character(1:9))

if (st_crs(ez_shp) != st_crs(4326)) ez_shp <- st_transform(ez_shp, 4326)

ez_shp <- ez_shp[!st_is_empty(ez_shp), ]
ez_shp <- st_make_valid(ez_shp)
ez_shp <- ez_shp[, c("EZ", "geometry")]
ez_shp <- ez_shp[!is.na(ez_shp$EZ), ]
ez_shp <- st_buffer(ez_shp, dist = 0)

crop_polygon <- st_as_sfc(st_bbox(c(
  xmin = xmin_crop, xmax = xmax_crop,
  ymin = ymin_crop, ymax = ymax_crop
), crs = 4326))

ez_cropped <- st_intersection(ez_shp, crop_polygon)
cat("Cropped to", nrow(ez_cropped), "tracts\n")


################################################################################
# 5. GBIF OCCURRENCES FOR THE INSET (cached), CLIPPED TO CALIFORNIA
################################################################################
gbif_cache <- "gbif_cenmel_cache.csv"

if (file.exists(gbif_cache)) {
  occs <- read.csv(gbif_cache)
  cat("Loaded", nrow(occs), "GBIF records from cache.\n")
} else {
  cat("Downloading C. melitensis occurrences from GBIF...\n")
  gbif_raw <- occ_search(
    scientificName     = "Centaurea melitensis",
    country            = "US",
    hasCoordinate      = TRUE,
    hasGeospatialIssue = FALSE,
    year               = "1990,2026",
    limit              = 100000
  )
  occs <- gbif_raw$data %>%
    transmute(lon = decimalLongitude, lat = decimalLatitude,
              unc = coordinateUncertaintyInMeters) %>%
    filter(!is.na(lon), !is.na(lat), is.na(unc) | unc <= 1000)
  write.csv(occs, gbif_cache, row.names = FALSE)
  cat("Saved", nrow(occs), "records to", gbif_cache, "\n")
}

# --- Keep ONLY occurrences that fall inside the California boundary -----------
ca_poly <- st_as_sf(maps::map("state", "california", fill = TRUE, plot = FALSE))
st_crs(ca_poly) <- 4326
ca_poly <- st_make_valid(st_union(ca_poly))
occs_sf <- st_as_sf(occs, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
occs <- occs[lengths(st_intersects(occs_sf, ca_poly)) > 0, ]
cat("GBIF occurrences inside California:", nrow(occs), "\n")


################################################################################
# 6. INSET — statewide GBIF map, clipped to CA, with the sampling-extent box
################################################################################
ca_outline <- map_data("state", region = "california")

fig_inset <- ggplot() +
  geom_polygon(data = ca_outline,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey40", linewidth = 0.5) +
  geom_point(data = occs, aes(x = lon, y = lat),
             shape = 21, size = 1.1, fill = "gold",
             color = "black", alpha = 0.55, stroke = 0.2) +
  annotate("rect",
           xmin = xmin_crop, xmax = xmax_crop,
           ymin = ymin_crop, ymax = ymax_crop,
           fill = NA, color = "black", linewidth = 1.0) +
  scale_x_continuous(breaks = c(-124, -120, -116)) +
  scale_y_continuous(breaks = c(34, 38, 42)) +
  coord_fixed(ratio = 1.25,
              xlim = c(-125, -114), ylim = c(32, 42.5), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    axis.title       = element_text(size = INSET_AXIS_TITLE, face = "bold", color = "black"),
    axis.text        = element_text(size = INSET_AXIS_TEXT,  color = "black"),
    axis.ticks       = element_line(linewidth = 0.4, color = "black"),
    panel.grid       = element_blank(),
    plot.background   = element_rect(fill = "white", color = "black", linewidth = 0.6),
    plot.margin      = margin(2, 3, 2, 2)
  )


################################################################################
# 7. MAIN PANEL — collection sites on the cropped EZ choropleth
################################################################################
fig_main <- ggplot() +
  geom_sf(data = ez_cropped, aes(fill = EZ), color = NA, alpha = 0.6) +
  geom_point(data = site_coords,
             aes(x = Longitude, y = Latitude, fill = EZ, shape = Classification),
             size = 3.2, color = "black", stroke = 0.9) +
  scale_fill_manual(values = zone_colors, name = "Environmental\nZone", drop = FALSE) +
  scale_shape_manual(values = c("Non-Urban" = 21, "Urban" = 24),
                     name = "Binary\nClassification", drop = FALSE) +
  coord_sf(xlim = c(xmin_crop, xmax_crop),
           ylim = c(ymin_crop, ymax_crop), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(
    axis.title        = element_text(size = MAIN_AXIS_TITLE, face = "bold", color = "black"),
    axis.text         = element_text(size = MAIN_AXIS_TEXT,  face = "bold", color = "black"),
    axis.line         = element_line(linewidth = 0.6, color = "black"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.title      = element_text(size = LEGEND_TITLE, face = "bold"),
    legend.text       = element_text(size = LEGEND_TEXT,  face = "bold"),
    legend.position   = legend_pos,
    legend.key        = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    legend.spacing.y  = unit(2, "pt"),
    legend.margin     = margin(3, 6, 3, 6)
  ) +
  guides(
    fill  = guide_legend(order = 1, override.aes = list(alpha = 0.9, shape = 22, size = 3.5)),
    shape = guide_legend(order = 2, override.aes = list(size = 3.2, fill = "grey50"))
  )


################################################################################
# 8. COMBINE (main + inset) AND SAVE (PNG + vector PDF)
################################################################################
fig <- fig_main +
  inset_element(fig_inset,
                left = inset_left, bottom = inset_bottom,
                right = inset_right, top = inset_top,
                align_to = "full")

print(fig)

ggsave("Ch1_Collection_Sites_EZ_Map.png", fig,
       width = fig_width_in, height = fig_height_in, dpi = 300, bg = "white")
ggsave("Ch1_Collection_Sites_EZ_Map.pdf", fig,
       width = fig_width_in, height = fig_height_in, bg = "white")

cat("\nSaved: Ch1_Collection_Sites_EZ_Map.png / .pdf at",
    fig_width_in, "x", fig_height_in, "in\n")
cat("Inset axis text:", INSET_AXIS_TEXT, "/", INSET_AXIS_TITLE, "pt\n")
