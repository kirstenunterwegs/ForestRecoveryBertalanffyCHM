########################################################################################
#
# Comparison of photogammetric nDOM and lidar derived CHM in Berchtesgaden National Park
#
########################################################################################

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Data Sources
# -------------------------------------------------------------------------
# 
# nDOM: German abbreviation for "normalisiertes Digitales Oberhöhen-Modell" 
#       (equivalent to a canopy height model, CHM), derived from aerial imagery.
# CHM:  Canopy Height Model derived from LiDAR data.
#
# -------------------------------------------------------------------------
# Comparison of Canopy Height Estimates
# -------------------------------------------------------------------------
#
# We compared canopy height estimates obtained from:
#   1. A LiDAR-derived Canopy Height Model (CHM) collected in August 2021, and
#   2. A photogrammetric nDOM (normalized Digital Surface Model), 
#      derived from aerial orthophotos acquired in summer 2020.
#
# Both datasets have a spatial resolution of 1 m.
#
# Height models were extracted for all combinations of:
#   - forest type,
#   - elevation bin, and
#   - gap size (based on gap mapping from Krüger et al., 2024).
#
# For each unique class combination, we randomly selected 20 sample points,
# resulting in a total of 2,400 paired height samples.
#
# -------------------------------------------------------------------------


# --- libraries

library(terra)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(cowplot)


# --- set working directory 

setwd("~/data") # or how you name your directory


######################################
# --- get stratified sample points ---
######################################

# load data layers 

elevation <- rast("02_dataRaw/bdg_environment/elevation_below1800_200steps.tif")
forest_type <- rast("02_dataRaw/bdg_environment/forest_type2020_reclass_1m.tif")
gaps <- rast("02_dataRaw/gaps_bdg/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
closed_forest <- vect("02_dataRaw/bdg_environment/closed_forest.gpkg")


# get gap size for stratified sampling on raster level

gaps.df <- terra::as.data.frame(gaps, na.rm = TRUE)
names(gaps.df) <- ("gap_ID")
saveRDS(gaps.df, "02_dataRaw/gaps_bdg/gap.df.rds")

gaps.reclass <- gaps.df %>%
  group_by(gap_ID) %>%
  count() %>%
  mutate(category = cut(n,
                        breaks = c(0, 500, 1000, 5000, 20000, Inf),
                        labels = c(1, 2, 3, 4, 5),
                        right = FALSE))

  
# Create a reclassification matrix

reclass_matrix <- gaps.reclass %>% 
  select(gap_ID, category) %>%
  mutate(
    gap_ID = as.numeric(gap_ID),
    category = as.numeric(factor(category, levels = c(1, 2, 3, 4, 5)))
  ) %>%
  distinct() %>%
  as.matrix()
  

# reclassify gaps into size categories

gap_size <- classify(gaps, rcl = reclass_matrix, right = FALSE) 

writeRaster(gap_size, "03_work/data_processed/nDOM_CHM_comp/gap_size_2021.tif")


# Stack the layers

layers_stack <- c(elevation, forest_type, gap_size)
names(layers_stack) <- c("elevation", "ftype", "gap_size")


# mask stack to forest area

layers_stack_forest <- mask(layers_stack, closed_forest)

# Extract values and create a data frame

values <- terra::as.data.frame(layers_stack_forest, xy=TRUE)

# filter rows with only NAs

values <- values[rowSums(is.na(values)) != ncol(values), ]

saveRDS(values, "03_work/data_processed/nDOM_CHM_comp/stratification_df.rds")


# Define the number of samples per stratum

samples_per_stratum <- 20  # Adjust as needed

# Stratified sampling

set.seed(123)  # Set seed for reproducibility
sampled_points <- values %>%
  group_by(elevation, ftype, gap_size) %>%
  sample_n(samples_per_stratum, replace = TRUE) %>%
  ungroup()

# Convert the data frame to points

spatial_points <- vect(sampled_points, geom=c("x", "y"), crs = crs(elevation))
writeVector(spatial_points, "03_work/data_processed/nDOM_CHM_comp/stratification_points.gpkg")

# Plot the points

plot(elevation)
points(spatial_points, col = 'red', pch = 16)


################################################
# --- Extract nDOM and CHM for sampel points---
################################################

# nDOM: German shortcut for "normalisiertes Digitales Oberhöhen-Modell" (CHM in English) - from areal images
# CHM: Canopy height model derived from lidar data

# the CHM layer from Berchtesgaden needs to be requested from the Berchtesgaden National Park Research division
# the nDOM for Berchtesgaden can be obtained from the Agency for Digitisation, High-Speed Internet and Surveying via the Open data Portal (https://geodaten.bayern.de/opengeodata/)

nDOM <- rast("03_work/data_processed/CHM/bdg/20/nDOM_20.tif")
CHM <- rast("E:/Projects/GapDynamics/data/processed/CHM_data/lidar/2021/berchtesgaden_2021_chm.tif")

# correct artifacts in Height Models (values <0 and >50 become threshold value )

nDOM_clamped <- clamp(nDOM, lower = 0, upper = 50, values = TRUE)
CHM_clamped <- clamp(CHM, lower = 0, upper = 50, values = TRUE)


# Extract raster values for CHM

CHM_values <- as.data.frame(extract(CHM_clamped, spatial_points, bind=TRUE, xy=TRUE)) %>% 
  mutate(type = "CHM")

# Extract raster values for nDOM

nDOM_values <- as.data.frame(extract(nDOM_clamped, spatial_points, bind=TRUE, xy=TRUE)) %>% 
  mutate(type = "nDOM") 
nDOM_z <- nDOM_values[,"Band_1"]

# create the comparison dataframe:

comparison_df1 <- cbind(CHM_values, nDOM_z)

nDOM_values <- nDOM_values %>%
  rename(Z = Band_1)

comparison_df2 <- rbind(CHM_values, nDOM_values)

# change labels for plotting & interpretation

comparison_df1 <- comparison_df1 %>%
  mutate(
    elevation = case_when(
      elevation == 1 ~ "600-800",
      elevation == 2 ~ "800-1000",
      elevation == 3 ~ "1000-1200",
      elevation == 4 ~ "1200-1400",
      elevation == 5 ~ "1400-1600",
      elevation == 6 ~ "1600-1800",
      elevation == 7 ~ "> 1800"
    ),
    ftype = case_when(
      ftype == 1 ~ "Beech",
      ftype == 2 ~ "Spruce-fir-beech",
      ftype == 4 ~ "Spruce",
      ftype == 5 ~ "Larch-pine"
    ),
    gap_size = case_when(
      gap_size == 1 ~ "<0.05 ha",
      gap_size == 2 ~ "0.05 - 0.1 ha",
      gap_size == 3 ~ "0.01 - 0.5 ha",
      gap_size == 4 ~ "0.5 - 2 ha",
      gap_size == 5 ~ ">2 ha"
    )
  )

comparison_df2 <- comparison_df2 %>%
  mutate(
    elevation = case_when(
      elevation == 1 ~ "600-800",
      elevation == 2 ~ "800-1000",
      elevation == 3 ~ "1000-1200",
      elevation == 4 ~ "1200-1400",
      elevation == 5 ~ "1400-1600",
      elevation == 6 ~ "1600-1800",
      elevation == 7 ~ "> 1800"
    ),
    ftype = case_when(
      ftype == 1 ~ "Beech",
      ftype == 2 ~ "Spruce-fir-beech",
      ftype == 4 ~ "Spruce",
      ftype == 5 ~ "Larch-pine"
    ),
    gap_size = case_when(
      gap_size == 1 ~ "<0.05 ha",
      gap_size == 2 ~ "0.05 - 0.1 ha",
      gap_size == 3 ~ "0.01 - 0.5 ha",
      gap_size == 4 ~ "0.5 - 2 ha",
      gap_size == 5 ~ ">2 ha"
    )
  )

saveRDS(comparison_df1, "03_work/data_processed/nDOM_CHM_comp/comparison_df1.rds")
saveRDS(comparison_df2, "03_work/data_processed/nDOM_CHM_comp/comparison_df2.rds")


# --- calculate height difference 

comparison_df1$diff <- comparison_df1$Z - comparison_df1$nDOM_z
comparison_df1$diff_category <- ifelse(comparison_df1$diff < -4, 1,
                                       ifelse(comparison_df1$diff > 4, 2, 0))

points_comparison <- vect(comparison_df1, geom=c("x", "y"), crs="EPSG:25832")
writeVector(points_comparison, "03_work/data_processed/nDOM_CHM_comp/points_comparison.gpkg", overwrite=T)

saveRDS(comparison_df1, "03_work/data_processed/nDOM_CHM_comp/comparison_df1_forest.rds")
saveRDS(comparison_df2, "03_work/data_processed/nDOM_CHM_comp/comparison_df2_forest.rds")

# comparison_df1 <- readRDS("03_work/data_processed/nDOM_CHM_comp/comparison_df1.rds")
# comparison_df2 <- readRDS("03_work/data_processed/nDOM_CHM_comp/comparison_df2.rds")


##########################################
# -- Analyse & plot height differences --
#########################################


# -- sort levels

# Define the desired order for gap_size and elevation

gap_size_levels <- c("<0.05 ha", "0.05 - 0.1 ha", "0.01 - 0.5 ha", "0.5 - 2 ha", ">2 ha")
elevation_levels <- c("600-800", "800-1000", "1000-1200", "1200-1400", "1400-1600", "1600-1800", "> 1800")

# Set the levels for gap_size and elevation

comparison_df1$gap_size <- factor(comparison_df1$gap_size, levels = gap_size_levels)
comparison_df1$elevation <- factor(comparison_df1$elevation, levels = elevation_levels)

comparison_df2$gap_size <- factor(comparison_df2$gap_size, levels = gap_size_levels)
comparison_df2$elevation <- factor(comparison_df2$elevation, levels = elevation_levels)

# --- Error metrics

# Mean Absolute Error (MAE) - less sensitive to outliers

MAE <- mean(na.omit(abs(comparison_df1$Z - comparison_df1$nDOM_z))) # 1.811365 

MAE_gap.size <- comparison_df1 %>% group_by(gap_size) %>%
  summarize(MAE_mean = mean(na.omit(abs(Z - nDOM_z))),
            MAE_median = median(na.omit(abs(Z - nDOM_z))))

MAE_elevation <- comparison_df1 %>% group_by(elevation) %>%
  summarize(MAE_mean = mean(na.omit(abs(Z - nDOM_z))),
            MAE_median = median(na.omit(abs(Z - nDOM_z))))

MAE_ftype <- comparison_df1 %>% group_by(ftype) %>%
  summarize(MAE_mean = mean(na.omit(abs(Z - nDOM_z))),
            MAE_median = median(na.omit(abs(Z - nDOM_z))))

# Root Mean Squared Error (RMSE) - penalizing larger errors

RMSE <- sqrt(mean(na.omit((comparison_df1$Z - comparison_df1$nDOM_z))^2))# 3.987143

RMSE_gap.size <- comparison_df1 %>% group_by(gap_size) %>% 
  summarize(RMSE_mean = sqrt(mean(na.omit((Z - nDOM_z)^2))))

RMSE_gap.elevation <- comparison_df1 %>% group_by(elevation) %>% 
  summarize(RMSE_mean = sqrt(mean(na.omit((Z - nDOM_z)^2))))


# ---------- Build label data frames ----------

# Overall labels (for non-faceted scatter)
lab_overall <- sprintf("MAE = %.2f\nRMSE = %.2f", MAE, RMSE)

# Per gap size

lab_gap <- MAE_gap.size %>%
  left_join(RMSE_gap.size, by = "gap_size") %>%
  transmute(
    gap_size,
    label = sprintf("MAE = %.2f\nRMSE = %.2f", MAE_mean, RMSE_mean),
    x = 0.5,  
    y = 8.5
  )

# Per elevation

lab_elev <- MAE_elevation %>%
  left_join(RMSE_gap.elevation, by = "elevation") %>%
  transmute(
    elevation,
    label = sprintf("MAE = %.2f\nRMSE = %.2f", MAE_mean, RMSE_mean),
    x = 0.5,
    y = 8.5
  )

# RMSE 

RMSE_ftype <- comparison_df1 %>%
  group_by(ftype) %>%
  summarize(RMSE_mean = sqrt(mean(na.omit((Z - nDOM_z)^2))))

lab_ftype <- MAE_ftype %>%
  left_join(RMSE_ftype, by = "ftype") %>%
  transmute(
    ftype,
    label = sprintf("MAE = %.2f\nRMSE = %.2f", MAE_mean, RMSE_mean),
    x = 1,    
    y = 19   
  )

# For boxplots faceted by elevation/gap_size you can reuse lab_elev and lab_gap:

lab_elev_box <- lab_elev %>% mutate(x = 1, y = 19)
lab_gap_box  <- lab_gap  %>% mutate(x = 1, y = 19)



# --- plots: --- 

# reverse labels for intuitive display

comparison_df1 <- comparison_df1 %>%
  mutate(
    elevation_plot = factor(elevation, levels = rev(elevation_levels)),
    gap_size_plot  = factor(gap_size,  levels = rev(gap_size_levels))
  )

lab_elev <- lab_elev %>%
  mutate(elevation_plot = factor(elevation, levels = rev(elevation_levels)))

lab_gap <- lab_gap %>%
  mutate(gap_size_plot = factor(gap_size, levels = rev(gap_size_levels)))


# --- 1) top: full-width simple scatter with overall MAE/RMSE ---

scatter_top <- ggplot(comparison_df1, aes(x = nDOM_z , y = Z)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "nDOM vegetation height (m)",
       y = "CHM vegetation height (m)") +
  theme_minimal() +
  xlim(0, 9) +  
  ylim(0, 9) +
  annotate("text", x = 0.5, y = 9, label = lab_overall, hjust = 0, vjust = 1)

# --- 2) left column: elevation (vertical stack; highest on top) ---

scatter_elevation_col <- ggplot(comparison_df1, aes(x = nDOM_z , y = Z)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "nDOM vegetation height (m)",
       y = "CHM vegetation height (m)") +
  theme_minimal() +
  facet_wrap(~elevation_plot, ncol = 1, drop = FALSE) +
  xlim(0, 9) +  
  ylim(0, 9) +
  geom_text(data = lab_elev, aes(x = 6, y = 8, label = label, group = elevation_plot),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5)+ 
  plot_annotation(title = "Per elevation bin") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# --- 3) right column: gap size (vertical stack; largest on top) ---

scatter_gapsize_col <- ggplot(comparison_df1, aes(x = nDOM_z , y = Z)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "nDOM vegetation height (m)",
       y = "CHM vegetation height (m)") +
  theme_minimal() +
  facet_wrap(~gap_size_plot, ncol = 1, drop = FALSE) +
  xlim(0, 9) +  
  ylim(0, 9) +
  geom_text(data = lab_gap, aes(x = 6.5, y = 8, label = label, group = gap_size_plot),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5)+ 
  plot_annotation(title = "Per gap size") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ---- assemble full plot -----

left_header  <- ggdraw() + draw_label("Per elevation bin", fontface = "bold", x = 0.5, hjust = 0.5)
right_header <- ggdraw() + draw_label("Per gap size",     fontface = "bold", x = 0.5, hjust = 0.5)

scatter_plot_all <- scatter_top /
  (left_header | right_header) /
  (scatter_elevation_col | scatter_gapsize_col) +
  plot_layout(heights = c(1, 0.08, 2))  

scatter_plot_all


# --- Boxplots for direct data comparison ---

# per elevation

box_elev <- ggplot(comparison_df2, aes(x = type, y = Z, fill = type)) +
  geom_boxplot() +
  facet_grid(~elevation) +
  labs( x = "",
        y = "Canopy Height") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.x = element_blank())+
  ylim(0, 20) +  # 35 values > 20 m
  geom_text(data = lab_elev_box, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.2)+ 
  ggtitle("Per elevation bin")

# per forest type

box_ftype <- ggplot(comparison_df2, aes(x = type, y = Z, fill = type)) +
  geom_boxplot() +
  facet_grid(~ftype) +
  labs( x = "",
        y = "Canopy Height") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.x = element_blank())+
  ylim(0, 20) +# 35 values > 20 m
  geom_text(data = lab_ftype, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.5)+ 
  ggtitle("Per forest type")

# per gap size

box_gapsize <- ggplot(comparison_df2, aes(x = type, y = Z, fill = type)) +
  geom_boxplot() +
  facet_grid(~gap_size) +
  labs( x = "",
        y = "Canopy Height") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.x = element_blank())+
  ylim(0, 20) +# 35 values > 20 m
  geom_text(data = lab_gap_box, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.5)+ 
  ggtitle("Per gap size")


comparisons_heights_all <- (box_elev / box_ftype / box_gapsize) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.margin = margin(t = 4, r = 6, b = 4, l = 6))

comparisons_heights_all


# defineoutput directory
output_dir <- "03_work/analysis/comparison_photoDOM_lidarCHM/results/"

# Save plots
ggsave(filename = file.path(output_dir, "scatter_plot_all.png"), plot = scatter_plot_all, width = 8, height = 12, units = "in")
ggsave(filename = file.path(output_dir, "boxplot_comparisons_all.png"), plot = comparisons_heights_all, width = 7, height = 10, units = "in")

