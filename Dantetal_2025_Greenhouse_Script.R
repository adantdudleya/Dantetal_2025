
###########################################################################
# This script performs statistical analyses and visualizations
# of greenhouse trait data for *Centaurea melitensis* (Cenmel),
# comparing Environmental Zones and Urban/Non-Urban treatments.
# Key analyses include:
#  - Nested ANOVAs by zone and location
#  - Post-hoc Tukey tests using emmeans
#  - Violin and boxplot visualizations of trait variation
#  - Treatment-level (Urban vs Non-Urban) and Land Classification analyses
###########################################################################

# -----------------------------
# Load Required Libraries
# -----------------------------
library(tidyverse)     
library(lmerTest)      
library(emmeans)       
library(multcomp)      
library(multcompView)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(sf)
library(patchwork)
library(maps)
library(rgbif)


# -----------------------------
# Define Environmental Zone Colors
# -----------------------------
# Zones 1, 2, and 7 are excluded due to no populations being sampled from any.
# This section defines a custom color palette for zones 3–6, 8, and 9.
available_zones <- c("3", "4", "5", "6", "8", "9")
full_palette <- brewer.pal(9, "Paired")
selected_indices <- c(3, 4, 5, 6, 8, 9)
my_colors <- full_palette[selected_indices]
names(my_colors) <- available_zones


# -----------------------------
# Data Loading and Cleaning
# -----------------------------

# ---- (Soil Treatment) ----

dat_plant <- read.csv("Greenhouse_data_Soil_treatment.csv", fileEncoding = "UTF-8-BOM")

# Inspect structure to verify column names and data types
str(dat_plant)

# Clean and filter dataset
dat_plant <- dat_plant %>%
  rename(
    # Standardize soil-treatment column names for the trait mapping
    Number_of_Leaves      = X..Number.of.leaves,
    Largest_Leaf_Height   = Height.of.the.largest.leaf..mm.,
    Largest_Leaf_Width    = Width.of.largest.leaf..mm.,
    Number_of_Flowerheads = Number.of.flowerheads
  ) %>%
  filter(Location != "Santa Cruz Island")  # Optional filtering step

# Convert Environmental Zone to factor
dat_plant$Environmental_Zone <- as.factor(dat_plant$Environmental_Zone)


# ---- (Sand Treatment) ----

dat_SLA <- read.csv("Greenhouse_data_Sand_treatment.csv", fileEncoding = "UTF-8-BOM")

# Inspect structure to verify columns
str(dat_SLA)

# Clean and filter dataset
dat_SLA <- dat_SLA %>%
  rename(
    # Standardize sand-treatment column names for the trait mapping
    Above_Ground_Biomass = Above.Ground.Biomass,
    Below_Ground_Biomass = Below.Ground.Biomass,
    Root_Diameter        = Root.Diameter,
    Specific_Leaf_Area   = SLA,
    Leaf_Area            = Area,
    Wet_Mass_g           = Wet.Mass..g.,
    Dry_Mass_g           = Dry.mass
  ) %>%
  filter(Location != "Santa Cruz Island")

# Convert Environmental Zone to factor
dat_SLA$Environmental_Zone <- as.factor(dat_SLA$Environmental_Zone)


# Display the number of unique collection sites
length(unique(dat_plant$Location))
length(unique(dat_SLA$Location))


# -----------------------------
# Trait Mapping Setup
# -----------------------------
# Define which dataset each trait belongs to and map

# Define dataset associations
trait_list <- list(
  "Specific Leaf Area (cm2/g)"  = dat_SLA,
  "Biomass (g)"                 = dat_plant,
  "Root to Shoot Ratio"         = dat_SLA,
  "Root Diameter (cm)"          = dat_SLA,
  "Number of Leaves"            = dat_plant,
  "Number of Flowerheads"       = dat_plant,
  "Width of Longest Leaf (mm)"  = dat_plant,
  "Length of Longest Leaf (mm)" = dat_plant
)

# Map trait labels to their respective column names
trait_column_map <- c(
  "Specific Leaf Area (cm2/g)"  = "Specific_Leaf_Area",
  "Biomass (g)"                 = "Biomass",
  "Root to Shoot Ratio"         = "Root_to_shoot_Ratio",
  "Root Diameter (cm)"          = "Root_Diameter",
  "Number of Leaves"            = "Number_of_Leaves",
  "Number of Flowerheads"       = "Number_of_Flowerheads",
  "Width of Longest Leaf (mm)"  = "Largest_Leaf_Width",
  "Length of Longest Leaf (mm)" = "Largest_Leaf_Height"
)

# ============================================================
# Nested ANOVA and Trait Visualization Script (original density version)
# ============================================================

# -----------------------------
# Nested ANOVA by Environmental Zone / Location
# -----------------------------
# For each trait, this section:
#  1. Fits a nested ANOVA (Environmental Zone / Location)
#  2. Performs post-hoc Tukey tests
#  3. Generates violin + boxplots with CLD letters

for (pretty_name in names(trait_list)) {
  df <- trait_list[[pretty_name]]
  raw_name <- trait_column_map[[pretty_name]]

  # Ensure proper factor levels
  df$Environmental_Zone <- factor(df$Environmental_Zone, levels = names(my_colors))
  df$Location <- as.factor(df$Location)

  # Nested ANOVA (Environmental Zone / Location)
  nested_model <- lm(as.formula(paste0("`", raw_name, "` ~ Environmental_Zone / Location")), data = df)
  cat("\n\n====================")
  cat(paste("\nNested ANOVA for", pretty_name))
  cat("\n====================\n")
  print(anova(nested_model))

  # Post-hoc Tukey comparison among Environmental Zones
  posthoc_model <- lm(as.formula(paste0("`", raw_name, "` ~ Environmental_Zone")), data = df)
  emm <- emmeans(posthoc_model, ~ Environmental_Zone)
  cld_out <- multcomp::cld(emm, Letters = letters, adjust = "tukey") %>% as.data.frame()
  cat("\nTukey groups (Environmental Zone):\n")
  print(cld_out)
}

# ##########################################################################
# REVISION ANALYSIS: Three-way classification comparison (Binary vs Gradient vs EZ)
# ##########################################################################

################################################################################
# CHAPTER 1 SUPPLEMENTARY ANALYSIS — COMPLETE SCRIPT
# 
# Input files: Greenhouse_data_Soil_treatment.csv, Greenhouse_data_Sand_treatment.csv, Impervious_Surface_Cenmel_Population.csv,
#              urban_typology2.shp (for map only)
#
# Produces:
#   1. Three-way trait comparison table (Binary vs Gradient vs EZ)
#   2. Violin plots for gradient classification (8 PNGs)
#   3. Collection site map with EZ overlay
#
# 
# 
################################################################################

library(tidyverse)
library(emmeans)
library(multcomp)
library(multcompView)
library(RColorBrewer)
library(sf)

# =============================================================================
# PART 0: SHARED SETUP
# =============================================================================

full_palette    <- brewer.pal(9, "Paired")
zone_colors     <- setNames(full_palette, as.character(1:9))
available_zones <- c("3", "4", "5", "6", "8", "9")

gradient_colors <- c(
  "Natural"       = "white",
  "Low"           = "gold",
  "Moderate-High" = "grey"
)

# =============================================================================
# PART 1: LOAD AND PREPARE DATA
# =============================================================================

cat("\n========================================\n")
cat("LOADING DATA\n")
cat("========================================\n")

# --- Plant trait data (Soil treatment) ---
dat_plant <- read.csv("Greenhouse_data_Soil_treatment.csv", fileEncoding = "UTF-8-BOM") %>%
  rename(
    Number_of_leaves      = X..Number.of.leaves,
    Largest_leaf_width    = Width.of.largest.leaf..mm.,
    Largest_leaf_height   = Height.of.the.largest.leaf..mm.,
    Number.of.flowerheads = Number.of.flowerheads
  ) %>%
  filter(Soil != "Sand", Location != "Santa Cruz Island")

dat_plant$Environmental_Zone <- as.factor(dat_plant$Environmental_Zone)
dat_plant$Location <- as.factor(dat_plant$Location)

cat("Plant data loaded:", nrow(dat_plant), "individuals\n")

# --- SLA trait data (Sand treatment) ---
dat_SLA <- read.csv("Greenhouse_data_Sand_treatment.csv", fileEncoding = "UTF-8-BOM") %>%
  filter(Location != "Santa Cruz Island")

dat_SLA$Environmental_Zone <- as.factor(dat_SLA$Environmental_Zone)
dat_SLA$Location <- as.factor(dat_SLA$Location)

cat("SLA data loaded:", nrow(dat_SLA), "individuals\n")

# --- Impervious surface lookup from Master Data ---
master <- read.csv("Impervious_Surface_Cenmel_Population.csv")

site_imperv <- master %>%
  group_by(Site) %>%
  summarise(Impervious_Surface = first(Impervious_Surface), .groups = "drop")


# --- Merge impervious and create gradient classification ---
dat_plant <- dat_plant %>%
  left_join(site_imperv, by = "Site") %>%
  mutate(
    Imperv_Gradient = cut(
      Impervious_Surface,
      breaks = c(-Inf, 5, 20, Inf),
      labels = c("Natural", "Low", "Moderate-High"),
      include.lowest = TRUE
    )
  )

dat_SLA <- dat_SLA %>%
  left_join(site_imperv, by = "Site") %>%
  mutate(
    Imperv_Gradient = cut(
      Impervious_Surface,
      breaks = c(-Inf, 5, 20, Inf),
      labels = c("Natural", "Low", "Moderate-High"),
      include.lowest = TRUE
    )
  )

# Verify
cat("\n--- Plant gradient distribution ---\n")
print(table(dat_plant$Imperv_Gradient, useNA = "always"))
cat("\n--- SLA gradient distribution ---\n")
print(table(dat_SLA$Imperv_Gradient, useNA = "always"))

cat("\nSites per gradient category:\n")
print(
  dat_plant %>%
    group_by(Imperv_Gradient) %>%
    summarise(
      n_individuals = n(),
      n_sites = n_distinct(Site),
      sites = paste(sort(unique(Site)), collapse = ", "),
      imperv_range = paste0(round(min(Impervious_Surface), 1), "–",
                            round(max(Impervious_Surface), 1), "%"),
      .groups = "drop"
    )
)

# =============================================================================
# PART 2: TRAIT DEFINITIONS
# =============================================================================

plant_traits <- list(
  "Biomass (g)"                  = "Biomass",
  "Number of Flowerheads"        = "Number.of.flowerheads",
  "Number of Leaves"             = "Number_of_leaves",
  "Width of Longest Leaf (mm)"   = "Largest_leaf_width",
  "Length of Longest Leaf (mm)"  = "Largest_leaf_height"
)

sla_traits <- list(
  "Specific Leaf Area (cm2/g)"   = "SLA",
  "Root to Shoot Ratio"          = "Root_to_shoot_Ratio",
  "Root Diameter (cm)"           = "Root.Diameter"
)

# =============================================================================
# PART 3: THREE-WAY NESTED ANOVA COMPARISON
# =============================================================================

cat("\n\n========================================\n")
cat("THREE-WAY COMPARISON: Binary vs Gradient vs EZ\n")
cat("========================================\n")

results <- data.frame()

# --- Plant traits (binary = Land.Classification) ---
for (pretty_name in names(plant_traits)) {
  col_name <- plant_traits[[pretty_name]]
  df <- dat_plant
  
  df$Land.Classification <- recode(df$Land.Classification, "NonUrban" = "Non-Urban")
  df$Land.Classification <- factor(df$Land.Classification, levels = c("Non-Urban", "Urban"))
  df$Environmental_Zone  <- factor(df$Environmental_Zone, levels = available_zones)
  
  # Model 1: Binary
  m_binary <- tryCatch({
    mod <- lm(as.formula(paste0("`", col_name, "` ~ Land.Classification / Location")), data = df)
    a <- anova(mod)
    list(
      class_F = a["Land.Classification", "F value"],
      class_p = a["Land.Classification", "Pr(>F)"],
      class_df = a["Land.Classification", "Df"],
      loc_F = a["Land.Classification:Location", "F value"],
      loc_p = a["Land.Classification:Location", "Pr(>F)"],
      loc_df = a["Land.Classification:Location", "Df"],
      resid_df = a["Residuals", "Df"]
    )
  }, error = function(e) list(class_F=NA, class_p=NA, class_df=NA, 
                              loc_F=NA, loc_p=NA, loc_df=NA, resid_df=NA))
  
  # Model 2: Gradient
  m_gradient <- tryCatch({
    mod <- lm(as.formula(paste0("`", col_name, "` ~ Imperv_Gradient / Location")), data = df)
    a <- anova(mod)
    list(
      class_F = a["Imperv_Gradient", "F value"],
      class_p = a["Imperv_Gradient", "Pr(>F)"],
      class_df = a["Imperv_Gradient", "Df"],
      loc_F = a["Imperv_Gradient:Location", "F value"],
      loc_p = a["Imperv_Gradient:Location", "Pr(>F)"],
      loc_df = a["Imperv_Gradient:Location", "Df"],
      resid_df = a["Residuals", "Df"]
    )
  }, error = function(e) list(class_F=NA, class_p=NA, class_df=NA,
                              loc_F=NA, loc_p=NA, loc_df=NA, resid_df=NA))
  
  # Model 3: Environmental Zones
  m_ez <- tryCatch({
    mod <- lm(as.formula(paste0("`", col_name, "` ~ Environmental_Zone / Location")), data = df)
    a <- anova(mod)
    list(
      class_F = a["Environmental_Zone", "F value"],
      class_p = a["Environmental_Zone", "Pr(>F)"],
      class_df = a["Environmental_Zone", "Df"],
      loc_F = a["Environmental_Zone:Location", "F value"],
      loc_p = a["Environmental_Zone:Location", "Pr(>F)"],
      loc_df = a["Environmental_Zone:Location", "Df"],
      resid_df = a["Residuals", "Df"]
    )
  }, error = function(e) list(class_F=NA, class_p=NA, class_df=NA,
                              loc_F=NA, loc_p=NA, loc_df=NA, resid_df=NA))
  
  results <- bind_rows(results, data.frame(
    Trait = pretty_name, Dataset = "Plant",
    # Binary
    Binary_F = m_binary$class_F, Binary_p = m_binary$class_p, Binary_df = m_binary$class_df,
    Binary_Loc_F = m_binary$loc_F, Binary_Loc_p = m_binary$loc_p, Binary_Loc_df = m_binary$loc_df,
    Binary_Resid_df = m_binary$resid_df,
    # Gradient
    Gradient_F = m_gradient$class_F, Gradient_p = m_gradient$class_p, Gradient_df = m_gradient$class_df,
    Gradient_Loc_F = m_gradient$loc_F, Gradient_Loc_p = m_gradient$loc_p, Gradient_Loc_df = m_gradient$loc_df,
    Gradient_Resid_df = m_gradient$resid_df,
    # EZ
    EZ_F = m_ez$class_F, EZ_p = m_ez$class_p, EZ_df = m_ez$class_df,
    EZ_Loc_F = m_ez$loc_F, EZ_Loc_p = m_ez$loc_p, EZ_Loc_df = m_ez$loc_df,
    EZ_Resid_df = m_ez$resid_df
  ))
  cat("Completed:", pretty_name, "\n")
}

# --- SLA traits (binary = Treatment) ---
for (pretty_name in names(sla_traits)) {
  col_name <- sla_traits[[pretty_name]]
  df <- dat_SLA
  
  df$Treatment <- recode(df$Treatment, "NonUrban" = "Non-Urban")
  df$Treatment <- factor(df$Treatment, levels = c("Non-Urban", "Urban"))
  df$Environmental_Zone <- factor(df$Environmental_Zone, levels = available_zones)
  
  # Model 1: Binary
  m_binary <- tryCatch({
    mod <- lm(as.formula(paste0("`", col_name, "` ~ Treatment / Location")), data = df)
    a <- anova(mod)
    list(
      class_F = a["Treatment", "F value"],
      class_p = a["Treatment", "Pr(>F)"],
      class_df = a["Treatment", "Df"],
      loc_F = a["Treatment:Location", "F value"],
      loc_p = a["Treatment:Location", "Pr(>F)"],
      loc_df = a["Treatment:Location", "Df"],
      resid_df = a["Residuals", "Df"]
    )
  }, error = function(e) list(class_F=NA, class_p=NA, class_df=NA,
                              loc_F=NA, loc_p=NA, loc_df=NA, resid_df=NA))
  
  # Model 2: Gradient
  m_gradient <- tryCatch({
    mod <- lm(as.formula(paste0("`", col_name, "` ~ Imperv_Gradient / Location")), data = df)
    a <- anova(mod)
    list(
      class_F = a["Imperv_Gradient", "F value"],
      class_p = a["Imperv_Gradient", "Pr(>F)"],
      class_df = a["Imperv_Gradient", "Df"],
      loc_F = a["Imperv_Gradient:Location", "F value"],
      loc_p = a["Imperv_Gradient:Location", "Pr(>F)"],
      loc_df = a["Imperv_Gradient:Location", "Df"],
      resid_df = a["Residuals", "Df"]
    )
  }, error = function(e) list(class_F=NA, class_p=NA, class_df=NA,
                              loc_F=NA, loc_p=NA, loc_df=NA, resid_df=NA))
  
  # Model 3: Environmental Zones
  m_ez <- tryCatch({
    mod <- lm(as.formula(paste0("`", col_name, "` ~ Environmental_Zone / Location")), data = df)
    a <- anova(mod)
    list(
      class_F = a["Environmental_Zone", "F value"],
      class_p = a["Environmental_Zone", "Pr(>F)"],
      class_df = a["Environmental_Zone", "Df"],
      loc_F = a["Environmental_Zone:Location", "F value"],
      loc_p = a["Environmental_Zone:Location", "Pr(>F)"],
      loc_df = a["Environmental_Zone:Location", "Df"],
      resid_df = a["Residuals", "Df"]
    )
  }, error = function(e) list(class_F=NA, class_p=NA, class_df=NA,
                              loc_F=NA, loc_p=NA, loc_df=NA, resid_df=NA))
  
  results <- bind_rows(results, data.frame(
    Trait = pretty_name, Dataset = "SLA",
    Binary_F = m_binary$class_F, Binary_p = m_binary$class_p, Binary_df = m_binary$class_df,
    Binary_Loc_F = m_binary$loc_F, Binary_Loc_p = m_binary$loc_p, Binary_Loc_df = m_binary$loc_df,
    Binary_Resid_df = m_binary$resid_df,
    Gradient_F = m_gradient$class_F, Gradient_p = m_gradient$class_p, Gradient_df = m_gradient$class_df,
    Gradient_Loc_F = m_gradient$loc_F, Gradient_Loc_p = m_gradient$loc_p, Gradient_Loc_df = m_gradient$loc_df,
    Gradient_Resid_df = m_gradient$resid_df,
    EZ_F = m_ez$class_F, EZ_p = m_ez$class_p, EZ_df = m_ez$class_df,
    EZ_Loc_F = m_ez$loc_F, EZ_Loc_p = m_ez$loc_p, EZ_Loc_df = m_ez$loc_df,
    EZ_Resid_df = m_ez$resid_df
  ))
  cat("Completed:", pretty_name, "\n")
}

# --- Format and display ---
results_display <- results %>%
  mutate(
    Binary_sig   = ifelse(Binary_p < 0.05, "*", ""),
    Gradient_sig = ifelse(Gradient_p < 0.05, "*", ""),
    EZ_sig       = ifelse(EZ_p < 0.05, "*", ""),
    across(ends_with("_F"), ~ round(.x, 5)),
    across(ends_with("_p"), ~ round(.x, 5))
  )

cat("\n================================================================\n")
cat("FULL RESULTS TABLE\n")
cat("================================================================\n\n")

# Classification-level results
cat("--- Classification Level ---\n")
print(
  results_display[, c("Trait",
                      "Binary_df", "Binary_F", "Binary_p", "Binary_sig",
                      "Gradient_df", "Gradient_F", "Gradient_p", "Gradient_sig",
                      "EZ_df", "EZ_F", "EZ_p", "EZ_sig")],
  row.names = FALSE
)

# Location-level results
cat("\n--- Location (nested) ---\n")
print(
  results_display[, c("Trait",
                      "Binary_Loc_df", "Binary_Loc_F", "Binary_Loc_p",
                      "Gradient_Loc_df", "Gradient_Loc_F", "Gradient_Loc_p",
                      "EZ_Loc_df", "EZ_Loc_F", "EZ_Loc_p")],
  row.names = FALSE
)

# Residual df
cat("\n--- Residual df ---\n")
print(
  results_display[, c("Trait", "Binary_Resid_df", "Gradient_Resid_df", "EZ_Resid_df")],
  row.names = FALSE
)

cat("\n=== Summary ===\n")
cat("Binary:  ", sum(results$Binary_p < 0.05, na.rm = TRUE), "/", nrow(results), "\n")
cat("Gradient:", sum(results$Gradient_p < 0.05, na.rm = TRUE), "/", nrow(results), "\n")
cat("EZ:      ", sum(results$EZ_p < 0.05, na.rm = TRUE), "/", nrow(results), "\n")

write.csv(results_display, "Supplementary_Gradient_Comparison_FULL.csv", row.names = FALSE)
cat("Saved: Supplementary_Gradient_Comparison_FULL.csv\n")

# =============================================================================
# NESTED CLASSIFICATION COMPARISON COMPLETE
# =============================================================================
cat("\n========================================\n")
cat("NESTED CLASSIFICATION COMPARISON COMPLETE\n")
cat("========================================\n")
cat("Nested models: trait ~ Classification / Location, under Binary, impervious\n")
cat("Gradient, and Environmental Zone schemes (Bonferroni-corrected across traits).\n")
cat("Output: Supplementary_Gradient_Comparison_FULL.csv\n")
cat("========================================\n")

# ##########################################################################
# REVISION ANALYSIS: Combined side-by-side violin figures (EZ vs Binary, EZ vs Gradient)
# ##########################################################################

################################################################################
# COMBINED VIOLIN PLOTS — TWO SEPARATE FIGURES
#
# Figure 1: EZ (left) | Gradient (right)
# Figure 2: EZ (left) | Binary (right)
#
# Input files: Greenhouse_data_Soil_treatment.csv, Greenhouse_data_Sand_treatment.csv, Impervious_Surface_Cenmel_Population.csv
# Greenhouse_data_Soil_treatment.csv = Soil treatment (Biomass, Leaves, Flowerheads, Leaf dims)
# Greenhouse_data_Sand_treatment.csv = Sand treatment (SLA, Root:Shoot, Root Diameter)
################################################################################

library(tidyverse)
library(emmeans)
library(multcomp)
library(multcompView)
library(RColorBrewer)
library(patchwork)

# =============================================================================
# COLORS & THEME
# =============================================================================

available_zones <- c("3", "4", "5", "6", "8", "9")
full_palette <- brewer.pal(9, "Paired")
ez_colors <- full_palette[c(3, 4, 5, 6, 8, 9)]
names(ez_colors) <- available_zones

gradient_colors <- c(
  "Natural"       = "white",
  "Low"           = "gold",
  "Moderate-High" = "grey"
)

binary_colors <- c("Non-Urban" = "#F4A582", "Urban" = "#92C5DE")

shared_theme <- theme_bw(base_size = 14) +
  theme(
    panel.border  = element_blank(),
    panel.grid    = element_blank(),
    axis.line     = element_line(color = "black"),
    legend.position = "none",
    axis.title.y  = element_text(size = 16, face = "bold"),
    axis.title.x  = element_text(size = 16, face = "bold"),
    axis.text.x   = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 15, face = "bold", color = "black")
  )

# =============================================================================
# LOAD DATA (loading + merging in one pipeline)
# =============================================================================

cat("Loading data...\n")

master <- read.csv("Impervious_Surface_Cenmel_Population.csv")
site_imperv <- master %>%
  group_by(Site) %>%
  summarise(Impervious_Surface = first(Impervious_Surface), .groups = "drop")

dat_soil <- read.csv("Greenhouse_data_Soil_treatment.csv", fileEncoding = "UTF-8-BOM") %>%
  filter(Soil != "Sand", Location != "Santa Cruz Island") %>%
  mutate(
    Environmental_Zone = factor(Environmental_Zone, levels = available_zones),
    Location = as.factor(Location)
  ) %>%
  left_join(site_imperv, by = "Site") %>%
  mutate(
    Imperv_Gradient = cut(Impervious_Surface,
                           breaks = c(-Inf, 5, 20, Inf),
                           labels = c("Natural", "Low", "Moderate-High"),
                           include.lowest = TRUE)
  )

dat_sand <- read.csv("Greenhouse_data_Sand_treatment.csv", fileEncoding = "UTF-8-BOM") %>%
  filter(Location != "Santa Cruz Island") %>%
  mutate(
    Environmental_Zone = factor(Environmental_Zone, levels = available_zones),
    Location = as.factor(Location)
  ) %>%
  left_join(site_imperv, by = "Site") %>%
  mutate(
    Imperv_Gradient = cut(Impervious_Surface,
                           breaks = c(-Inf, 5, 20, Inf),
                           labels = c("Natural", "Low", "Moderate-High"),
                           include.lowest = TRUE)
  )

cat("Soil treatment:", nrow(dat_soil), "rows\n")
cat("Sand treatment:", nrow(dat_sand), "rows\n")
cat("Soil Imperv_Gradient:\n"); print(table(dat_soil$Imperv_Gradient))
cat("Sand Imperv_Gradient:\n"); print(table(dat_sand$Imperv_Gradient))

# =============================================================================
# TRAIT DEFINITIONS
# =============================================================================

all_traits <- list(
  "Length of Longest Leaf (mm)"  = list(data = "soil", col = "Height.of.the.largest.leaf..mm.", binary = "Land.Classification"),
  "Specific Leaf Area (cm2/g)"  = list(data = "sand", col = "SLA", binary = "Treatment"),
  "Number of Flowerheads"       = list(data = "soil", col = "Number.of.flowerheads", binary = "Land.Classification"),
  "Root Diameter (cm)"          = list(data = "sand", col = "Root.Diameter", binary = "Treatment"),
  "Width of Longest Leaf (mm)"  = list(data = "soil", col = "Width.of.largest.leaf..mm.", binary = "Land.Classification"),
  "Biomass (g)"                 = list(data = "soil", col = "Biomass", binary = "Land.Classification"),
  "Number of Leaves"            = list(data = "soil", col = "X..Number.of.leaves", binary = "Land.Classification"),
  "Root to Shoot Ratio"         = list(data = "sand", col = "Root_to_shoot_Ratio", binary = "Treatment")
)

get_data <- function(dataset_name) {
  if (dataset_name == "sand") return(dat_sand)
  if (dataset_name == "soil") return(dat_soil)
}

# Only the traits that passed the Bonferroni correction (shown in the figures)
traits_to_show <- c(
  "Length of Longest Leaf (mm)",
  "Specific Leaf Area (cm2/g)",
  "Number of Flowerheads",
  "Root Diameter (cm)"
)

# =============================================================================
# FUNCTION: Make EZ violin plot
# =============================================================================

make_ez_violin <- function(df, col_name, pretty_name) {
  df$Environmental_Zone <- factor(df$Environmental_Zone, levels = available_zones)
  
  posthoc_model <- lm(as.formula(paste0("`", col_name, "` ~ Environmental_Zone")), data = df)
  emm <- emmeans(posthoc_model, ~ Environmental_Zone)
  cld_out <- tryCatch({
    multcomp::cld(emm, Letters = letters, adjust = "tukey") %>% as.data.frame()
  }, error = function(e) NULL)
  
  df_plot <- df %>% mutate(value = .data[[col_name]])
  
  if (!is.null(cld_out)) {
    violin_tips <- df_plot %>%
      group_by(Environmental_Zone) %>%
      summarise(
        y_pos = {
          vals <- value[!is.na(value)]
          if (length(vals) >= 2) max(density(vals)$x) else max(vals, na.rm = TRUE)
        },
        .groups = "drop"
      )
    cld_out <- left_join(cld_out, violin_tips, by = "Environmental_Zone")
    max_val <- max(df_plot$value, na.rm = TRUE)
    cld_out$label_y <- cld_out$y_pos + 0.02 * max_val
  }
  
  p <- ggplot(df_plot, aes(x = Environmental_Zone, y = value, fill = Environmental_Zone)) +
    geom_violin(width = 0.9, alpha = 0.4, trim = FALSE) +
    geom_boxplot(width = 0.10, color = "grey20", fill = "white", outlier.size = 3) +
    scale_fill_manual(values = ez_colors) +
    labs(x = "Environmental Zone", y = pretty_name) +
    shared_theme
  
  if (!is.null(cld_out)) {
    p <- p + geom_text(data = cld_out,
                        aes(x = Environmental_Zone, y = label_y, label = .group),
                        inherit.aes = FALSE, fontface = "bold", size = 5)
  }
  return(p)
}

# =============================================================================
# FUNCTION: Make binary violin plot
# =============================================================================

make_binary_violin <- function(df, col_name, pretty_name, binary_col) {
  df[[binary_col]] <- recode(df[[binary_col]], "NonUrban" = "Non-Urban")
  df[[binary_col]] <- factor(df[[binary_col]], levels = c("Non-Urban", "Urban"))
  
  model <- lm(as.formula(paste0("`", col_name, "` ~ `", binary_col, "`")), data = df)
  a <- anova(model)
  p_val <- a[binary_col, "Pr(>F)"]
  sig_label <- ifelse(p_val < 0.05,
                       paste0("p = ", format(round(p_val, 4), nsmall = 4), " *"),
                       paste0("p = ", format(round(p_val, 4), nsmall = 4), " (n.s.)"))
  
  p <- ggplot(df, aes(x = .data[[binary_col]], y = .data[[col_name]],
                        fill = .data[[binary_col]])) +
    geom_violin(width = 0.9, alpha = 0.4, trim = FALSE) +
    geom_boxplot(width = 0.05, fill = "white", color = "grey20", outlier.size = 3) +
    scale_fill_manual(values = binary_colors) +
    labs(title = NULL, x = NULL, y = pretty_name) +
    shared_theme +
    theme(plot.title = element_text(size = 10, face = "italic", hjust = 0.5))
  
  return(p)
}

# =============================================================================
# FUNCTION: Make gradient violin plot
# =============================================================================

make_gradient_violin <- function(df, col_name, pretty_name) {
  df$Imperv_Gradient <- factor(df$Imperv_Gradient,
                                levels = c("Natural", "Low", "Moderate-High"))
  
  model <- lm(as.formula(paste0("`", col_name, "` ~ Imperv_Gradient")), data = df)
  a <- anova(model)
  p_val <- a["Imperv_Gradient", "Pr(>F)"]
  sig_label <- ifelse(p_val < 0.05,
                       paste0("p = ", format(round(p_val, 4), nsmall = 4), " *"),
                       paste0("p = ", format(round(p_val, 4), nsmall = 4), " (n.s.)"))
  
  emm <- emmeans(model, ~ Imperv_Gradient)
  cld_out <- tryCatch({
    multcomp::cld(emm, Letters = letters, adjust = "tukey") %>% as.data.frame()
  }, error = function(e) NULL)
  
  if (!is.null(cld_out)) {
    violin_tips <- df %>%
      group_by(Imperv_Gradient) %>%
      summarise(
        y_pos = {
          vals <- .data[[col_name]][!is.na(.data[[col_name]])]
          if (length(vals) >= 2) max(density(vals)$x) else max(vals, na.rm = TRUE)
        },
        .groups = "drop"
      )
    cld_out <- left_join(cld_out, violin_tips, by = "Imperv_Gradient")
    max_val <- max(df[[col_name]], na.rm = TRUE)
    cld_out$label_y <- cld_out$y_pos + 0.02 * max_val
  }
  
  p <- ggplot(df, aes(x = Imperv_Gradient, y = .data[[col_name]],
                        fill = Imperv_Gradient)) +
    geom_violin(width = 0.9, alpha = 0.4, trim = FALSE) +
    geom_boxplot(width = 0.10, color = "grey20", fill = "white", outlier.size = 3) +
    scale_fill_manual(values = gradient_colors) +
    scale_x_discrete(labels = c("Natural" = "<5%",
                                 "Low" = "5-19%",
                                 "Moderate-High" = "20-49%")) +
    labs(title = NULL, x = "Impervious Surface Gradient", y = pretty_name) +
    shared_theme +
    theme(plot.title = element_text(size = 10, face = "italic", hjust = 0.5))
  
  if (!is.null(cld_out)) {
    p <- p + geom_text(data = cld_out,
                        aes(x = Imperv_Gradient, y = label_y, label = .group),
                        inherit.aes = FALSE, fontface = "bold", size = 5)
  }
  return(p)
}

# =============================================================================
# FIGURE 1: EZ (left) | Gradient (right)
# =============================================================================

cat("\n========================================\n")
cat("FIGURE 1: EZ + GRADIENT\n")
cat("========================================\n")

ez_grad_rows <- list()
panel_labels <- letters[1:length(traits_to_show)]

for (i in seq_along(traits_to_show)) {
  trait <- traits_to_show[i]
  info <- all_traits[[trait]]
  df <- get_data(info$data)
  col <- info$col
  
  cat("Processing:", trait, "\n")
  
  left <- make_ez_violin(df, col, trait)
  
  right <- make_gradient_violin(df, col, trait) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("(", panel_labels[i], ")"),
             hjust = 1.2, vjust = 1.5, fontface = "bold", size = 5)
  
  ez_grad_rows[[i]] <- left | right
}

combined_ez_grad <- wrap_plots(ez_grad_rows, ncol = 1)
print(combined_ez_grad)
ggsave("Figure_EZ_vs_Gradient.png", combined_ez_grad,
       width = 14, height = 4 * length(traits_to_show),
       dpi = 300, bg = "white")
cat("Saved: Figure_EZ_vs_Gradient.png\n")

# =============================================================================
# FIGURE 2: EZ (left) | Binary (right)
# =============================================================================

cat("\n========================================\n")
cat("FIGURE 2: EZ + BINARY\n")
cat("========================================\n")

ez_bin_rows <- list()

for (i in seq_along(traits_to_show)) {
  trait <- traits_to_show[i]
  info <- all_traits[[trait]]
  df <- get_data(info$data)
  col <- info$col
  binary_col <- info$binary
  
  cat("Processing:", trait, "\n")
  
  left <- make_ez_violin(df, col, trait)
  
  right <- make_binary_violin(df, col, trait, binary_col) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("(", panel_labels[i], ")"),
             hjust = 1.2, vjust = 1.5, fontface = "bold", size = 5)
  
  ez_bin_rows[[i]] <- left | right
}

combined_ez_bin <- wrap_plots(ez_bin_rows, ncol = 1)
print(combined_ez_bin)
ggsave("Figure_EZ_vs_Binary.png", combined_ez_bin,
       width = 14, height = 4 * length(traits_to_show),
       dpi = 300, bg = "white")
cat("Saved: Figure_EZ_vs_Binary.png\n")

cat("\n========================================\n")
cat("DONE — TWO FIGURES SAVED\n")
cat("========================================\n")
cat("Figure_EZ_vs_Gradient.png  (EZ | Impervious Surface Gradient)\n")
cat("Figure_EZ_vs_Binary.png    (EZ | Urban vs Non-Urban)\n")

