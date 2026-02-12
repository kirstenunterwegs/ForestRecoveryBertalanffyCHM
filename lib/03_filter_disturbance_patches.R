##############################################################
#
# Filtering and processing of forest disturbance data
#
#############################################################


# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
#
# This script prepares and filters forest disturbance patches across Bavaria 
# for subsequent analysis of natural disturbance recovery and canopy height dynamics.
#
# Specifically, it:
#   1. Merges all annual disturbance patch rasters and assigns unique patch IDs 
#      across years.
#   2. Masks patches to the Bavarian administrative boundary.
#   3. Excludes anthropogenic disturbances (e.g., harvest) based on disturbance
#      agent classification (Senf & Seidl, 2021).
#   4. Removes strictly protected areas (designated after 1986) not part of the 
#      study’s sampling design.
#   5. Masks orthophoto survey areas to include only those captured within the 
#      vegetation period (June–September) and with at least two temporal observations.
#   6. Converts filtered disturbance rasters to polygon features, calculates patch 
#      areas, and subsets to minimum patch sizes (0.1 ha and 0.25 ha).
#   7. Generates buffer layers for:
#        - surrounding intact forest canopy (outer 500 m),
#        - inner patch edges (10 m exclusion zone),
#      and creates a raster-based inner-edge mask.
#
# Output:
#   - Merged and filtered disturbance rasters 
#   - Disturbance patch polygons by area threshold 
#   - Buffer layers and edge masks 
#
# Notes:
#   - All spatial data are processed in EPSG:3035 (ETRS89 / LAEA Europe) 
#     unless otherwise stated.
#   - Some outputs are reprojected to EPSG:25832 (ETRS89 / UTM zone 32N) 
#     for integration with local datasets.
# -------------------------------------------------------------------------

# --- load libraries

library(terra)
library(sf)

# set directory

setwd("~/data") # or how you name your directory

################################################################
# --- load and merge disturbance patches into one layer ---
################################################################

# assign an individual ID to each disturbance patch
# disturbance patches derived from Senf & Seidl 2021 (https://doi.org/10.1111/gcb.15679)

folder_path <- "02_dataRaw/disturbances/patches/"

# List all .tif files
tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)

# Initialize an empty list to store relabeled rasters
relabeled_rasters <- list()

# Initialize the starting value for unique IDs
current_max_id <- 0

# Loop through each patch file and relabel patch id, to get unique ID across all years

for (tif_file in tif_files) {
  
  print(tif_file)
  
  raster <- rast(tif_file)
  
  # Get patch IDs
  unique_ids <- unique(values(raster))
  unique_ids <- unique_ids[!is.na(unique_ids)]
  
  # Create a relabeling map
  new_ids <- seq(current_max_id + 1, current_max_id + length(unique_ids))
  id_map <- setNames(new_ids, unique_ids)
  
  # Relabel the raster
  relabeled_raster <- classify(raster, data.frame(unique_ids, new_ids), others = NA)
  
  # Update the current max ID
  current_max_id <- max(new_ids)
  
  # Append to the list of relabeled rasters
  relabeled_rasters[[length(relabeled_rasters) + 1]] <- relabeled_raster
}

# Merge all relabeled rasters
merged_raster <- do.call(merge, relabeled_rasters)

# Write the combined raster to a new file
writeRaster(merged_raster, "03_work/data_processed/disturbances/patches/distrubance_patches_merges_unique_ids.tif", overwrite = TRUE)


################################
# --- mask to Bavaria --- (aoi)
###############################

bavaria <- vect("02_dataRaw/admin/bayern_Grenze.gpkg")
bavaria_3035 <- project(bavaria,"EPSG:3035")

patches <- rast("03_work/data_processed/disturbances/patches/distrubance_patches_merges_unique_ids.tif")
crs(patches) <- "EPSG:3035"

patches_bavaria <- mask(crop(patches, bavaria_3035), bavaria_3035)

writeRaster(patches_bavaria, "03_work/data_processed/disturbances/patches/disturbance_patches_merged_unique_ids_bavaria.tif", overwrite = TRUE)
#patches_bavaria <- rast("03_work/data_processed/disturbances/patches/disturbance_patches_merged_unique_ids_bavaria.tif")

# overall disturbed area in bavaria
n_non_na <- global(!is.na(patches_bavaria), "sum", na.rm = TRUE)[1, 1]
n_non_na * 0.09 # 482,403 ha



####################################################################
# --- subset to only disturbances initialized by natural agents ---
####################################################################

# load agent attribution from Senf & Seidl (2021). Global Change Biology

agent <- rast("02_dataRaw/disturbances/agent_classes_germany.tif")
crs(agent) <- "EPSG:3035"
agent_bavaria <- mask(crop(agent, bavaria_3035), bavaria_3035) # mask to Bavaria

# set all disturbances classified as "harvest" to NA to exclude those

patches_bavaria[agent_bavaria == 3] <- NA
writeRaster(patches_bavaria, "03_work/data_processed/disturbances/patches/disturbance_patches_merged_unique_ids_bavaria_natural.agent.tif", overwrite = TRUE)

# natural disturbed area in bavaria
n_non_na <- global(!is.na(patches_bavaria), "sum", na.rm = TRUE)[1, 1]
n_non_na * 0.09 # 170,362.8 ha



#####################################################################################
# --- mask out strictly protected areas not in or confounding our study design  ---
#####################################################################################

# ------------------------------------------------------------------------------------------
# We remove all protected areas designated after the onset of the disturbance mapping
# to exclude zones with potentially confounding or biasing management or conservation measures 
# affecting forest recovery. 
# This includes:
#   - Strict forest reserves designated after 1986
#   - Biosphere reserve core zones (excluded)
#   - National park buffer zones
# ------------------------------------------------------------------------------------------

# --- strict forest reserves

reserves <- vect("02_dataRaw/set_aside_forests/protected_area_all/Old_Reserves_230110.gpkg")
reserves.out <- reserves[reserves$ausweisung > 1986] # ausweisung (german) = designation

reserves.in <- reserves[reserves$ausweisung < 1986]
writeVector(reserves.in, "02_dataRaw/set_aside_forests/selected_set_aside/strict_forest_reserves_pre_1986.gpkg")

# --- biosphere reserves core zones (excluding BDG area) 

bio.reserve <- vect("02_dataRaw/set_aside_forests/protected_area_all/biospherereserve_roehn_core_zones.gpkg")

# --- national park buffer zones

bdg <- vect("02_dataRaw/set_aside_forests/protected_area_all/nps/bdg/npb_zonierung_22_epsg25832.shp")
bdg.buffer <- bdg[bdg$zone_txt == "Pflegezone"]

baywald <- vect("02_dataRaw/set_aside_forests/protected_area_all/nps/baywald/Nationalpark_(Alt-_und_Erweiterungsgebiet).shp")
baywald.buffer <- baywald[baywald$TEXT == "Erweiterungsgebiet"]
baywald.buffer <- project(baywald.buffer, "EPSG:25832")


# merge protected area masking layers 

mask_protected_areas <- rbind(reserves.out, bio.reserve, bdg.buffer,baywald.buffer)

#writeVector(mask_protected_areas, "02_dataRaw/set_aside_forests/protected_area_all/mask_protected_areas.gpkg")
#mask_protected_areas <- vect("02_dataRaw/set_aside_forests/protected_area_all/mask_protected_areas.gpkg")

# create mask layer

mask_layer_protected <- erase(bavaria, mask_protected_areas)
mask_layer_protected_3035 <- project(mask_layer_protected, "EPSG:3035")# reproject mask layer

# mask disturbance patches

patches_bavaria_protected.masked <- mask(patches_bavaria, mask_layer_protected_3035, touches=FALSE)

writeRaster(patches_bavaria_protected.masked, "03_work/data_processed/disturbances/patches/disturbance_patches_merged_unique_ids_bavaria_natural.agent_protected.masked.tif", overwrite = TRUE)


# natural disturbed area in Bavaria filtered for protected areas not in our design

n_non_na <- global(!is.na(patches_bavaria_protected.masked), "sum", na.rm = TRUE)[1, 1]
n_non_na * 0.09 # 170,362.8 ha



#######################################################################################################
# --- mask out Orthophoto survey dates pre June and post September to avoid distortions in the CHMs ---
#######################################################################################################

survey2017 <- vect("02_dataRaw/DOP_survey_dates/2017/MD_UTM32_BildFluguebersicht_2017.shp")
survey2018 <- vect("02_dataRaw/DOP_survey_dates/2018/MD_UTM32_BildFluguebersicht_2018.shp")
survey2019 <- vect("02_dataRaw/DOP_survey_dates/2019/MD_UTM32_BildFluguebersicht_2019.shp")
survey2020 <- vect("02_dataRaw/DOP_survey_dates/2020/MD_UTM32_BildFluguebersicht_2020.shp")
survey2021 <- vect("02_dataRaw/DOP_survey_dates/2021/MD_UTM32_BildFluguebersicht_2021.shp")
survey2022.2023 <- vect("02_dataRaw/DOP_survey_dates/2022_2023/Befliegungsübersicht_2022_2023.shp")

# subset flight areas to those months, where survey was in June, July, August or September

subset_survey2017 <- survey2017[grep("^\\d{4}(06|07|08|09)", as.character(survey2017$UA)), ]
subset_survey2018 <- survey2018[grep("^\\d{4}/(06|07|08|09)", as.character(survey2018$FLUGDATUM)), ]
subset_survey2019 <- survey2019[grep("^\\d{4}/(06|07|08|09)", as.character(survey2019$FLUGDATUM)), ]
subset_survey2020 <- survey2020[grep("^\\d{4}/(06|07|08|09)", as.character(survey2020$FLUGDATUM)), ]
subset_survey2021 <- survey2021[grep("^\\d{4}/(06|07|08|09)", as.character(survey2021$FLUGDATUM)), ]
subset_survey2022.2023 <- survey2022.2023[grep("^\\d{4}/(06|07|08|09)", as.character(survey2022.2023$flugdatum)), ]


# ---create masking layer with at least two observations to subset distrubances

# combine all subsetted survey layers into one list

all_layers_list <- list(subset_survey2017, subset_survey2018, subset_survey2019,
                        subset_survey2020, subset_survey2021, subset_survey2022.2023)

# Find all pairwise intersections to ensure we have at least two observations

overlaps <- list()


for (i in seq_along(all_layers_list)) {
  for (j in seq_along(all_layers_list)) {
    if (i < j) { # Avoid redundant comparisons and self-comparisons
      overlap <- intersect(all_layers_list[[i]], all_layers_list[[j]])
      if (!is.null(overlap)) {
        overlaps <- append(overlaps, list(overlap))
      }
    }
  }
}

# Combine all overlaps into a single SpatVector
if (length(overlaps) > 0) {
  combined_overlaps <- do.call(rbind, overlaps)
  print(combined_overlaps)
} else {
  print("No overlaps found.")
}

names(combined_overlaps) <- make.unique(names(combined_overlaps)) 

#writeVector(combined_overlaps, "03_work/data_processed/survey_date/areas_min.two_surveys.gpkg", overwrite=TRUE)
#combined_overlaps <- vect("03_work/data_processed/survey_date/areas_min.two_surveys.shp")

combined_overlaps.3035 <- project(combined_overlaps, "EPSG:3035")

# mask disturbance patches

patches_bavaria_protected_survey.masked <- mask(patches_bavaria_protected.masked, combined_overlaps.3035, touches=FALSE)

#writeRaster(patches_bavaria_protected_survey.masked, "03_work/data_processed/disturbances/patches/disturbance_patches_merged_unique_ids_bavaria_natural.agent_protected_survey.masked.tif", overwrite = TRUE)
#patches_bavaria_protected_survey.masked <- rast("03_work/data_processed/disturbances/patches/disturbance_patches_merged_unique_ids_bavaria_natural.agent_protected_survey.masked.tif")

# natural disturbed area in bavaria filtered for protected areas not in our design and masked by survey date

n_non_na <- global(!is.na(patches_bavaria_protected_survey.masked), "sum", na.rm = TRUE)[1, 1]
n_non_na * 0.09 # 170,362.8 ha


#######################################################################################################
# --- create disturbance patch polygons and filter by size and create distrubance patch edge masks
#######################################################################################################

# Note: (if outcome is 0 geometry, restart R and try again from here)

disturbance.patch.polygons <- as.polygons(patches_bavaria_protected_survey.masked, aggregate = T, values=T, na.rm=TRUE)


# calculate area of each patch

disturbance.patch.polygons$area.ha <- expanse(disturbance.patch.polygons, unit="ha", transform=TRUE)
sum(disturbance.patch.polygons$area.ha) # 133689.1 ha - 48473 patches

writeVector(disturbance.patch.polygons, "03_work/data_processed/disturbances/patches/disturbance_patches_filtered_all.gpkg", overwrite=T)


# subset to patches larger than 0.25 ha  

disturbance.patch.polygons_sub_0.25 <- disturbance.patch.polygons[disturbance.patch.polygons$area.ha > 0.25]

sum(disturbance.patch.polygons_sub_0.25$area.ha)# 133346.7 ha - 45675 patches (~ 6% filtered out)
writeVector(disturbance.patch.polygons_sub_0.25, "03_work/data_processed/disturbances/patches/disturbance_patches_filtered_0.25ha.shp", overwrite=T)  


# -- load disturbance polygons for further processing

disturbance.patch.polygons_sub <- vect("03_work/data_processed/disturbances/patches/disturbance_patches_filtered_0.25ha.gpkg")
disturbance.patch.polygons.25832 <- project(disturbance.patch.polygons_sub, "EPSG:25832")

writeVector(disturbance.patch.polygons.25832,"03_work/data_processed/disturbances/patches/disturbance_patches_filtered_0.25ha_25832.gpkg" )
disturbance.patch.polygons_sub <- vect("03_work/data_processed/disturbances/patches/disturbance_patches_filtered_0.25ha_25832.gpkg")



#-- buffer polygons to extract surrounding (intact) forest canopy and mask edge area of distrubance patches on the inside

disturbance_sf <- st_as_sf(disturbance.patch.polygons_sub)

# outer buffer 500m 

outer_buffer_500 <- st_buffer(disturbance_sf, dist = 500)
outer_buffer_500 <- vect(outer_buffer_500)
writeVector(outer_buffer_500, "03_work/data_processed/disturbances/patches/buffer/outer_buffer_500.gpkg") 

outer_buffer_area <- erase(outer_buffer_500, disturbance.patch.polygons_sub)
writeVector(outer_buffer_area, "03_work/data_processed/disturbances/patches/buffer/outer_buffer_area_500.gpkg") 

# inner buffer 10 m

inner_buffer_10 <- st_buffer(disturbance_sf, dist = -10)
inner_buffer_10 <- vect(inner_buffer_10)
writeVector(inner_buffer_10, "03_work/data_processed/disturbances/patches/buffer/inner_buffer_10.gpkg")  


# --- remove inner pixel row of disturbance patches to create a patch edge mask ---

disturbance.patches.raster <- rast("03_work/data_processed/disturbances/patches/disturbance_patches_merged_unique_ids_bavaria_natural.agent_protected_survey.masked.tif")

edge_mask <- ifel(!is.na(disturbance.patches.raster), 1, NA)
edge_mask_focal <- focal(edge_mask, w=3, na.policy= "all", fun="mean")

writeRaster(edge_mask_focal, "03_work/data_processed/disturbances/patches/buffer/inner_patch_edge_mask.tif", overwrite = T)

edge_mask_polygon <- as.polygons(edge_mask_focal)
edge_mask_polygon <- project(edge_mask_polygon, "EPSG:25832")
writeVector(edge_mask_polygon, "03_work/data_processed/disturbances/patches/buffer/inner_patch_edge_mask.gpkg") 

# area of inner patch edge mask
edge_mask_polygon$area.ha <- expanse(edge_mask_polygon, unit="ha", transform=TRUE)
sum(edge_mask_polygon$area.ha)# 32068.08 ha


rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

