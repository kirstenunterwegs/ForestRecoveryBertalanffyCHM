######################################################################
#
# Spatial Overview of Sample Points for the analysis
#
######################################################################

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
#
# This script creates overview maps of the study region, sample locations,
# and example CHM tiles used to illustrate post-disturbance canopy structure.
#
# Specifically, it:
#   1. Loads Bavarian administrative boundaries and the spatial sample
#      points coloured by management type.
#   2. Produces an inset map of mainland Europe highlighting Bavaria.
#   3. Loads three example CHM patches (~3, 15, 30 years post disturbance)
#      and converts them to comparable raster data frames.
#   4. Builds zoomed CHM panels (A–C) with consistent height scaling and
#      labels, linked to their locations in Bavaria.
#   5. Combines sample points, CHM panels, and inset map into a main
#      overview figure and a management-type specific overview.
#
# Output:
#   - 03_work/results/overview_map.tif
#   - 03_work/results/supporting_information/overview_map_per_mngt.tif
#
# Notes:
#   - Example CHM tiles are *not* included in the repository and must be
#     obtained from the Bavarian State Institute of Forestry or the Bavarian
#     Agency for Digitisation, High-Speed Internet and Surveying.
# -------------------------------------------------------------------------


# ---libraries

library(ggplot2)
library(terra)
library(MetBrewer)
library(data.table)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(patchwork)
library(cowplot)

# set directory

setwd("~/NAS/Projects/ForestRecovery/")



# -------------------------------
# load data layer for plotting
# -------------------------------

# Bavaria & surrondings
bavaria = st_read("02_dataRaw/admin/bayern_Grenze.gpkg")

# points coloured by mngt category
sample_points <- readRDS("03_work/data_processed/analysis.dt/spatial_sub250.rds")

sample_points[, mngt_type := factor(mngt_type,
                                    levels = 0:6,
                                    labels = c("Set aside", "Large private", "Small private",
                                               "State forest", "Municipal forest", "Federal forest", "Other"))]
obs_points <- st_as_sf(sample_points, coords = c("x", "y"), crs = 25832)



# --- note: these CHM tiles are not available in this data repository and need to be obtained from the The Bavarian State Institute of Forestry
#           or from the Bavarian Agency for Digitisation, High-Speed Internet and Surveying

chm_3 = rast("03_work/data_processed/map/CHM_patches/patch_2088676_2023_30m_innerbuffer.tif")
chm_15= rast("03_work/data_processed/map/CHM_patches/patch_844352_2018_30m_innerbuffer.tif")
chm_30 = rast("03_work/data_processed/map/CHM_patches/patch_331527_2020_30m_innerbuffer.tif")


# CHMS - viridis: 1* for ~ 3 years post dist; 1* 15 years post dist ; 1* 30 years post dist



# -- load background map --

world <- ne_countries(scale = "medium", returnclass = "sf")

europe <- world[
  world$continent == "Europe" & 
    world$sovereignt != "Russia" &
    world$sovereignt != "Iceland" &
    world$type != "Dependency", 
]

centroids <- st_centroid(europe)
coords <- st_coordinates(centroids)

# Keep only mainland Europe based on rough bounding box
europe_mainland <- europe[
  coords[, "X"] > -10 &  # Exclude far west (e.g., Azores, Madeira)
    coords[, "Y"] > 35,    # Exclude south (e.g., Canary Islands)
]

# Plot
overview_p= ggplot(data = europe_mainland) +
  geom_sf(fill = "white", color = "black") +
  geom_sf(data = bavaria, fill = "black", color = "black", alpha = 0.8, size = 0.5) +
  coord_sf(xlim = c(-10, 45), ylim = c(35, 72)) +
  theme_void() 

overview_p

# Overview sample location:

europe_mainland_utm <- st_transform(europe_mainland, crs = 25832)


bbox <- st_bbox(bavaria)
padding_x <- 0.7 # Add some padding (in degrees)
padding_y <- 0.7


# --- define management colours ---

mngt_colors <- c(
  "State forest" = "#A00E00",
  "Other"     = "#beaed4"  ,      
  "Federal forest" = "#DC7400" , 
  "Set aside"    = "#F6C200",
  "Municipal forest"     = "#529A6F",
  "Small private"    = "#132B69",
  "Large private" = "#e7298a"
)
# ----------------------------------


# Get centroids of CHM extents

chm_3_centroid <- st_centroid(st_as_sfc(st_bbox(chm_3)))
chm_15_centroid <- st_centroid(st_as_sfc(st_bbox(chm_15)))
chm_30_centroid <- st_centroid(st_as_sfc(st_bbox(chm_30)))


# Combine centroids with labels

centroids <- rbind(
  st_sf(label = "A", geometry = chm_3_centroid),
  st_sf(label = "B", geometry = chm_15_centroid),
  st_sf(label = "C", geometry = chm_30_centroid)
)


# Convert sf to data frame with coordinates

obs_points_jitter = cbind(obs_points, st_coordinates(obs_points))


# Add jitter to x and y (adjust amount as needed)

set.seed(42)  # For reproducibility
obs_points_jitter$x_jit = obs_points_jitter$X + rnorm(nrow(obs_points_jitter), mean = 0, sd = 500)  # 500 meters jitter
obs_points_jitter$y_jit = obs_points_jitter$Y + rnorm(nrow(obs_points_jitter), mean = 0, sd = 500)



# --- plot main map  ---

# sample_plot = ggplot() +
#   
#   # Bavaria highlight
#   geom_sf(data = bavaria, fill= NA, color = "black", alpha = 0.5) +
#   
#   # Observation points colored by management type
#   geom_point(data = obs_points_jitter, aes(x = x_jit, y = y_jit, color = mngt_type), size = 0.6) +
#   #geom_sf(data = obs_points, aes(color = mngt_type), size = 0.3) +
#   scale_color_manual(values = mngt_colors, name = "Management type") +
#   
#   # Add CHM inset points
#   geom_sf(data = centroids, color = "white", fill = "lightblue", shape = 21, size = 12, alpha = 0.7) +
#   geom_sf(data = centroids, color = "grey60", fill = "darkblue", shape = 21, size = 8, alpha = 1) +
#   geom_sf_text(data = centroids, aes(label = label), fontface = "bold", size = 4 , color = "white") + #,nudge_x = 5000, nudge_y = 5000
#   # scale
#   annotation_scale(location = "bl", width_hint = 0.5) +
#   # Zoom to Bavaria area
#   coord_sf(
#     xlim = c(bbox["xmin"] - padding_x, bbox["xmax"] + padding_x),
#     ylim = c(bbox["ymin"] - padding_y, bbox["ymax"] + padding_y)
#   ) +
#   theme_void() +
#   theme(
#     legend.position = c(0, 0.35),   # Middle-left inside plot
#     legend.justification = c(0, 0.5), 
#     legend.background = element_rect(fill = alpha("white", 0.7), color = NA)  # Optional: semi-transparent legend background
#   ) +
#   guides(color = guide_legend(override.aes = list(size = 4))) 
# sample_plot

#------------------ set aside on top

sample_plot = ggplot() +
  geom_sf(data = bavaria, fill= NA, color = "black", alpha = 0.5) +
  
  geom_point(
    data = subset(obs_points_jitter, mngt_type != "Set aside"),
    aes(x = x_jit, y = y_jit, color = mngt_type),
    size = 0.6
  ) +
  geom_point(
    data = subset(obs_points_jitter, mngt_type == "Set aside"),
    aes(x = x_jit, y = y_jit, color = mngt_type),
    size = 0.9
  ) +
  
  scale_color_manual(values = mngt_colors, name = "Management type") +
  
  geom_sf(data = centroids, color = "white", fill = "lightblue", shape = 21, size = 12, alpha = 0.7) +
  geom_sf(data = centroids, color = "grey60", fill = "darkblue", shape = 21, size = 8, alpha = 1) +
  geom_sf_text(data = centroids, aes(label = label), fontface = "bold", size = 4 , color = "white") +
  annotation_scale(location = "bl", width_hint = 0.5) +
  coord_sf(
    xlim = c(bbox["xmin"] - padding_x, bbox["xmax"] + padding_x),
    ylim = c(bbox["ymin"] - padding_y, bbox["ymax"] + padding_y)
  ) +
  theme_void() +
  theme(
    legend.position = c(0, 0.35),
    legend.justification = c(0, 0.5),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))
sample_plot

#-------------------

chm_3_df <- as.data.frame(chm_3, xy = TRUE)
chm_15_df <- as.data.frame(chm_15, xy = TRUE)
chm_30_df <- as.data.frame(chm_30, xy = TRUE)

# Rename value columns
names(chm_3_df)[3] <- "height"
names(chm_15_df)[3] <- "height"
names(chm_30_df)[3] <- "height"


# Extract min and max from each dataset
min_height <- min(c(min(chm_3_df$height, na.rm = TRUE),
                    min(chm_15_df$height, na.rm = TRUE),
                    min(chm_30_df$height, na.rm = TRUE)))

max_height <- max(c(max(chm_3_df$height, na.rm = TRUE),
                    max(chm_15_df$height, na.rm = TRUE),
                    max(chm_30_df$height, na.rm = TRUE)))



# --- CHM plots ---
# ------------------------------------------------------------

# --- CHM 3 Plot with Label ---


# communal forest - pine forest 3 years post disturbance

center_x <- 501850
center_y <- 5545400

# Define zoom extent (50 meters in each direction)
x_min <- center_x - 50
x_max <- center_x + 50
y_min <- center_y - 50
y_max <- center_y + 50

# Define padding from edges (adjust as needed)
padding_x <- 5
padding_y <- 5

# Position label slightly inside top-left corner
label_x <- x_min + padding_x
label_y <- y_max - padding_y

# Update ggplot with zoom window

p_chm3 <- ggplot() +
  geom_raster(data = chm_3_df, aes(x = x, y = y, fill = height)) +
  scale_fill_viridis_c(option = "D", limits = c(min_height, max_height)) +
  annotate("text", x = label_x, y = label_y, label = "A", size = 8, fontface = "bold", hjust = 0, vjust = 1, color = "white") +
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
  theme_void()

p_chm3


# --------- CHM 15 Plot with Label
# small private forest - beech 15 years post disturbance

center_x <- 768880
center_y <- 5294840

# Define zoom extent (50 meters in each direction)
x_min <- center_x - 50
x_max <- center_x + 50
y_min <- center_y - 50
y_max <- center_y + 50

# Define padding from edges (adjust as needed)
padding_x <- 5
padding_y <- 5

# Position label slightly inside top-left corner
label_x <- x_min + padding_x
label_y <- y_max - padding_y

# Update ggplot with zoom window

p_chm15 <- ggplot() +
  geom_raster(data = chm_15_df, aes(x = x, y = y, fill = height)) +
  scale_fill_viridis_c(option = "D", limits = c(min_height, max_height)) +
  labs(fill = "Height [m]") +   # <-- legend title
  annotate("text", x = label_x, y = label_y, label = "B", size = 8, fontface = "bold", hjust = 0, vjust = 1, color = "white") +
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
  theme_void()

p_chm15


# --------- CHM 30 Plot with Label
# large private forest - spruce 30 years post disturbance


center_x <- 652450
center_y <- 5342300

# Define zoom extent (50 meters in each direction)
x_min <- center_x - 50
x_max <- center_x + 50
y_min <- center_y - 50
y_max <- center_y + 50

# Define padding from edges (adjust as needed)
padding_x <- 5
padding_y <- 5

# Position label slightly inside top-left corner
label_x <- x_min + padding_x
label_y <- y_max - padding_y

# Update ggplot with zoom window

p_chm30 <- ggplot() +
  geom_raster(data = chm_30_df, aes(x = x, y = y, fill = height)) +
  scale_fill_viridis_c(option = "D", limits = c(min_height, max_height)) +
  annotate("text", x = label_x, y = label_y, label = "C", size = 8, fontface = "bold", hjust = 0, vjust = 1, color = "white") +
  annotation_scale(location = "bl", width_hint = 0.5) +
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
  theme_void()

p_chm30

# ------------------------------------------------------------



# -- sort CHM plots

p_chm3_legend <- p_chm3 + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "grey60", fill = NA, size = 1))

p_chm15_nolegend <- p_chm15 + 
  theme(legend.position = "right",
        panel.border = element_rect(color = "grey60", fill = NA, size = 1))

p_chm30_nolegend <- p_chm30 + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "grey60", fill = NA, size = 1))



# --- Stack plots with spacing

# Stack plots and add legend beside the middle plot

chm_column <- plot_grid(
  p_chm3_legend, 
  NULL, 
  plot_grid(p_chm15_nolegend), 
  NULL, 
  p_chm30_nolegend,
  ncol = 1,
  rel_heights = c(1, 0.02, 1, 0.02, 1),
  align = "v"
)

chm_column


# Place overview_p on top of sample_plot

sample_plot_inset = ggdraw() +
  draw_plot(sample_plot) +  # Base plot
  draw_plot(overview_p, x = 0.69, y = 0.62, width = 0.4, height = 0.4)  # Inset map in top-right

sample_plot_inset


# Combine with overview map, both columns with matching height

final_plot <- (sample_plot_inset | chm_column) + 
  plot_layout(widths = c(2, 1))


final_plot


ggsave("03_work/results/overview_map.tif", plot = final_plot, dpi = 300, width = 10, height = 6)





# ----------------------------------------------
#  plot with separate sample plots per mngt
# ----------------------------------------------

# Vector of management types you want to plot

mngt_types <- c("Municipal forest", "Small private", "Large private", 
                "Set aside", "Federal forest", "State forest", "Other")


# Function to create sample_plot per category

create_mngt_plot <- function(mngt_name, label_letter = NULL) {
  
  obs_subset <- obs_points_jitter %>% dplyr::filter(mngt_type == mngt_name)
  
  p <- ggplot() +
    geom_sf(data = bavaria, fill = NA, color = "black", alpha = 0.5) +
    geom_point(data = obs_subset, aes(x = x_jit, y = y_jit, color = mngt_type), size = 0.5) +
    scale_color_manual(values = mngt_colors, guide = "none") +
    coord_sf(xlim = c(bbox["xmin"] - padding_x, bbox["xmax"] + padding_x),
             ylim = c(bbox["ymin"] - padding_y, bbox["ymax"] + padding_y)) +
    theme_void() +
    theme(panel.border = element_rect(color = "grey60", fill = NA, size = 1),
          plot.title = element_text(hjust = 0.5, size = 10)) +
    ggtitle(mngt_name)
  
  return(p)
}


p_communal <- create_mngt_plot("Municipal forest")
p_small <- create_mngt_plot("Small private")
p_large <- create_mngt_plot("Large private")

p_communal = p_communal +
  geom_sf(data = centroids[centroids$label == "A", ], color = "grey60", fill = "darkblue", shape = 21, size = 6, alpha = 1) +
  geom_sf_text(data = centroids[centroids$label == "A", ], aes(label = label), fontface = "bold", size = 3, color = "white")+ theme(plot.margin = margin(0, 0, 0, 0))
p_small = p_small +
  geom_sf(data = centroids[centroids$label == "B", ], color = "grey60", fill = "darkblue", shape = 21, size = 6, alpha = 1) +
  geom_sf_text(data = centroids[centroids$label == "B", ], aes(label = label), fontface = "bold", size = 3, color = "white")+ theme(plot.margin = margin(0, 0, 0, 0))
p_large = p_large +
  geom_sf(data = centroids[centroids$label == "C", ], color = "grey60", fill = "darkblue", shape = 21, size = 6, alpha = 1) +
  geom_sf_text(data = centroids[centroids$label == "C", ], aes(label = label), fontface = "bold", size = 3, color = "white")+ theme(plot.margin = margin(0, 0, 0, 0))


p_set_aside <- create_mngt_plot("Set aside")+ theme(plot.margin = margin(0, 0, 0, 0))
p_federal <- create_mngt_plot("Federal forest")+ theme(plot.margin = margin(0, 0, 0, 0))
p_state <- create_mngt_plot("State forest") +
  ggspatial::annotation_scale(
    location  = "br",                 # bottom-right
    width_hint = 0.5,
    height     = unit(0.2, "cm"),
    pad_x      = unit(0.2, "cm"),
    pad_y      = unit(0.05, "cm")   
  ) +
  theme(plot.margin = margin(0, 0, 0, 0))

p_other <- create_mngt_plot("Other")+ theme(plot.margin = margin(0, 0, 0, 0))


p_chm3_legend <- p_chm3 + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "grey60", fill = NA, size = 1))+ theme(plot.margin = margin(0, 0, 0, 0))

p_chm15_nolegend <- p_chm15 + 
  guides(
    fill = guide_colorbar(
      title = "Height [m]",
      title.position = "left",   
      title.hjust = 1,           
      barwidth = unit(4, "cm")   
    )
  ) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(margin = margin(r = 6)),  # small gap between title and bar
    panel.border = element_rect(color = "grey60", fill = NA, size = 1),
    plot.margin = margin(0, 0, 0, 0)
  )

p_chm30_nolegend <- p_chm30 + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "grey60", fill = NA, size = 1))+ theme(plot.margin = margin(0, 0, 0, 0))


# Row 1: CHM plots + horizontal legend
row1 <- plot_grid(p_chm3_legend, p_chm15_nolegend, p_chm30_nolegend, ncol = 3, align = "hv")


# Row 2: First three management types
row2 <- plot_grid(p_communal, p_small, p_large, ncol = 3, align = "hv")

# Row 3: Next three management types
row3 <- plot_grid(p_set_aside, p_federal, p_state, ncol = 3, align = "hv")

# Row 4: Other + overview
row4 <- plot_grid(p_other, overview_p,  ncol = 3,  align = "hv")

# Final grid
final_layout <- plot_grid(
  row1,
  row2,
  row3,
  row4,
  ncol = 1,
  rel_heights = c(1.12, 1, 1, 1)  # First row slightly taller
)

final_layout


ggsave("03_work/results/supporting_information/overview_map_per_mngt.tif", plot = final_layout, dpi = 300, width = 6, height = 8)

