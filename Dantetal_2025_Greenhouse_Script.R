
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

dat_plant <- read.csv("Greenhouse_data_Soil_treatment.csv")

# Inspect structure to verify column names and data types
str(dat_plant)

# Clean and filter dataset
dat_plant <- dat_plant %>%
  rename(
    # Standardize column names for consistency with trait mapping
    Above_Ground_Biomass = Above.Ground.Biomass,
    Below_Ground_Biomass = Below.Ground.Biomass,
    Root_to_shoot_Ratio  = Root_to_shoot_Ratio,
    Root_Diameter        = Root.Diameter,
    Specific_Leaf_Area   = SLA,
    Leaf_Area            = Area,
    Wet_Mass_g           = Wet.Mass..g.,
    Dry_Mass_g           = Dry.mass
  ) %>%
  filter(Location != "Santa Cruz Island")  # Optional filtering step

# Convert Environmental Zone to factor
dat_plant$Environmental_Zone <- as.factor(dat_plant$Environmental_Zone)


# ---- (Sand Treatment) ----

dat_SLA <- read.csv("Greenhouse_data_Sand_treatment.csv")

# Inspect structure to verify columns
str(dat_SLA)

# Clean and filter dataset
dat_SLA <- dat_SLA %>%
  rename(
    Number_of_Leaves     = X..Number.of.leaves,
    Largest_Leaf_Height  = Height.of.the.largest.leaf..mm.,
    Largest_Leaf_Width   = Width.of.largest.leaf..mm.,
    Number_of_Flowerheads = Number.of.flowerheads,
    Biomass              = Biomass
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
  "Specific Leaf Area (cm2/g)"  = dat_plant,
  "Biomass (g)"                 = dat_SLA,
  "Root to Shoot Ratio"         = dat_plant,
  "Root Diameter (cm)"          = dat_plant,
  "Number of Leaves"            = dat_SLA,
  "Number of Flowerheads"       = dat_SLA,
  "Width of Longest Leaf (mm)"  = dat_SLA,
  "Length of Longest Leaf (mm)" = dat_SLA
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
  
  # Nested ANOVA
  nested_model <- lm(as.formula(paste0("`", raw_name, "` ~ Environmental_Zone / Location")), data = df)
  cat("\n\n====================")
  cat(paste("\nNested ANOVA for", pretty_name))
  cat("\n====================\n")
  print(anova(nested_model))
  
  # Post-hoc comparison by Environmental Zone
  posthoc_model <- lm(as.formula(paste0("`", raw_name, "` ~ Environmental_Zone")), data = df)
  emm <- emmeans(posthoc_model, ~ Environmental_Zone)
  
  cld_out <- multcomp::cld(emm, Letters = letters, adjust = "tukey") %>%
    as.data.frame() %>%
    rename(Environmental_Zone = Environmental_Zone)
  
  # Prepare data for plotting
  df_plot <- df %>%
    mutate(value = .data[[raw_name]])
  
  # Estimate violin "tip" height for CLD label placement
  # (Original behavior: uses raw density for label placement)
  violin_tips <- df_plot %>%
    group_by(Environmental_Zone) %>%
    summarise(
      y_pos = {
        vals <- value[!is.na(value)]
        if (length(vals) >= 2) {
          dens <- density(vals, na.rm = TRUE)
          max(dens$x)
        } else {
          max(vals, na.rm = TRUE) # fallback if single observation
        }
      },
      .groups = "drop"
    )
  
  cld_out <- left_join(cld_out, violin_tips, by = "Environmental_Zone")
  
  # Plot violin + boxplot per zone
  p <- ggplot(df_plot, aes(x = Environmental_Zone, y = value, fill = Environmental_Zone)) +
    geom_violin(width = 0.9, alpha = 0.4, trim = FALSE) +
    geom_boxplot(width = 0.10, color = "grey20", fill = "white", outlier.size = 3) +
    geom_text(
      data = cld_out,
      aes(x = Environmental_Zone, y = y_pos + 0.02 * max(df_plot$value, na.rm = TRUE), label = .group),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 12
    ) +
    scale_fill_manual(values = my_colors) +
    labs(
      x = "Environmental Zone",
      y = pretty_name
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none",
      axis.title.y = element_text(size = 40, face = "bold"),
      axis.title.x = element_text(size = 40, face = "bold"),
      axis.text.x = element_text(size = 35, face = "bold", color = "black"),
      axis.text.y = element_text(size = 35, face = "bold", color = "black")
    )
  
  print(p)
}


# -----------------------------
# Sand treatment trait dataset (Urban vs Non-Urban)
# -----------------------------
# Performs ANOVA for SLA-based traits between treatments
# and plots violin + boxplots colored by treatment class.

sla_traits <- c(
  "Specific Leaf Area (cm2/g)",
  "Root to Shoot Ratio",
  "Root Diameter (cm)"
)

dat_SLA$Treatment <- recode(dat_SLA$Treatment, "NonUrban" = "Non-Urban")
dat_SLA$Treatment <- factor(dat_SLA$Treatment, levels = c("Non-Urban", "Urban"))

for (pretty_name in sla_traits) {
  trait <- trait_column_map[[pretty_name]]
  df <- dat_SLA
  
  model <- lm(as.formula(paste(trait, "~ Treatment")), data = df)
  cat("\n\n====================")
  cat(paste("\nANOVA for", pretty_name, "by Treatment"))
  cat("\n====================\n")
  print(anova(model))
  
  p <- ggplot(df, aes(x = Treatment, y = .data[[trait]], fill = Treatment)) +
    geom_violin(width = 0.9, alpha = 0.4, trim = FALSE) +
    geom_boxplot(width = 0.05, fill = "white", color = "grey20", outlier.size = 3) +
    scale_fill_manual(values = c("Non-Urban" = "#F4A582", "Urban" = "#92C5DE")) +
    labs(x = NULL, y = pretty_name) +
    theme_classic() +
    theme(
      axis.title.y = element_text(size = 40, face = "bold"),
      axis.text.x = element_text(size = 40, face = "bold", color = "black"),
      axis.text.y = element_text(size = 35, face = "bold", color = "black"),
      legend.position = "none"
    )
  
  print(p)
}


# -----------------------------
# Soil treatment trait dataset (Urban vs Non-Urban)
# -----------------------------
# Repeats ANOVA + violin plots for morphological traits.

plant_traits <- c(
  "Biomass (g)", "Number of Flowerheads",
  "Number of Leaves", "Width of Longest Leaf (mm)",
  "Length of Longest Leaf (mm)"
)

for (pretty_name in plant_traits) {
  trait <- trait_column_map[[pretty_name]]
  df <- dat_plant
  
  df$Land.Classification <- recode(df$Land.Classification, "NonUrban" = "Non-Urban")
  df$Land.Classification <- factor(df$Land.Classification, levels = c("Non-Urban", "Urban"))
  
  model <- lm(as.formula(paste(trait, "~ Land.Classification")), data = df)
  cat("\n\n====================")
  cat(paste("\nANOVA for", pretty_name, "by Land.Classification"))
  cat("\n====================\n")
  print(anova(model))
  
  p <- ggplot(df, aes(x = Land.Classification, y = .data[[trait]], fill = Land.Classification)) +
    geom_violin(width = 0.9, alpha = 0.4, trim = FALSE) +
    geom_boxplot(width = 0.05, fill = "white", color = "grey20", outlier.size = 3) +
    scale_fill_manual(values = c("Non-Urban" = "#F4A582", "Urban" = "#92C5DE")) +
    labs(x = NULL, y = pretty_name) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 40, face = "bold", color = "black"),
      axis.text.y = element_text(size = 40, face = "bold", color = "black"),
      axis.title = element_text(size = 35, face = "bold", color = "black"),
      legend.position = "none"
    )
  
  print(p)
}












