###########################################################################
#
# Canopy Height Model Tile Processing & Extraction for Disturbance Patches 
#
##########################################################################


# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
#
# This script prepares, filters, and extracts CHM data for all disturbance 
# patches in Bavaria for the CHMs across the years 2017–2023.
#
# Specifically, it:
#   1. Builds polygon grids for all nDOM/CHM tiles and stores their file paths.
#   2. Filters annual orthophoto survey areas to the vegetation period 
#      (June–September) and masks tile grids accordingly.
#   3. Identifies all CHM tiles intersecting each disturbance patch and its 
#      buffers, then loads and merges the required tiles.
#   4. Crops and masks CHM rasters to patch boundaries, outer forested buffers, 
#      and inner buffer zones (10 m and 30 m), excluding disturbed or non-forest areas.
#   5. Runs all patch–year extractions in parallel
#
# Output:
#   - Annual tile-grid layers (masked to vegetation-period surveys)
#   - CHM rasters for patch interiors, buffer areas, and inner-edge masks
#   - Error log and summary of extracted patch–year combinations
#
# Notes:
#   - CHM tiles are processed in EPSG:25832; forest/disturbance layers are 
#     temporarily projected from EPSG:3035 where needed.
# -------------------------------------------------------------------------

# ---libraries

library(terra)
library(tictoc)
library(parallel)
library(dplyr)
library(stringr)

# set directory

setwd("~/NAS/Projects/ForestRecovery/")



###########################################################################################
# --- create a vector file grid of all CHM tiles with respective file path as attribute ---
###########################################################################################

years <- c("2017", "2018", "2019", "2020", "2021", "2022", "2023") # years for which CHM is available

tile_paths <- c("/data/public/Projects/DANK/ForestRecovery/02_dataRaw/nDOM_lwf/nDOM_17_neu", # paths for server access to all CHMs
                "/data/public/Projects/DANK/ForestRecovery/02_dataRaw/nDOM_lwf/nDOM_18",
                "/data/public/Projects/DANK/ForestRecovery/02_dataRaw/nDOM_lwf/nDOM_19",
                "/data/public/Projects/DANK/ForestRecovery/02_dataRaw/nDOM_lwf/nDOM_20",
                "/data/public/Projects/DANK/ForestRecovery/02_dataRaw/nDOM_lwf/nDOM_21",
                "/data/public/Projects/DANK/ForestRecovery/02_dataRaw/nDOM_lwf/nDOM_22",
                "/data/public/Projects/DANK/ForestRecovery/02_dataRaw/nDOM_lwf/nDOM_23")

# extent extraction function 

extract.ext = function(file_tile){
  tile = rast(file_tile)
  crs(tile) <- "EPSG:25832"
  ext_tile = as.polygons(ext(tile), crs = terra::crs(tile))
  ext_tile$file_tile = file_tile
  return(ext_tile)  
}

# loop thorough each year

for (i in seq_along(years)) {
  year <- years[i]
  dir_tiles <- tile_paths[i]
  
  cat("Processing year:", year, "\n")
  
  # list files for the current year
  tic("List files")
  files_tiles_YYYY <- list.files(dir_tiles, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  toc()
  
  # process all tiles for the current year
  ext_tiles_YYYY <- vector(mode = "list", length = length(files_tiles_YYYY))
  
  tic(paste("Process all tiles for year", year))
  for (j in 1:length(files_tiles_YYYY)) {
    file_tile <- files_tiles_YYYY[j]
    cat(j, "Processing tile:", file_tile, "\n")
    ext_tiles_YYYY[[j]] <- extract.ext(file_tile)
  }
  toc()
  
  # create single vector grid of scans
  grid_YYYY <- vect(ext_tiles_YYYY)
  
  # Write the output 
  output_path <- paste0("03_work/data_processed/CHM/tile_grids/grid_", year, ".gpkg") 
  writeVector(grid_YYYY, output_path)
  
  cat("Finished processing for year:", year, "\n")
  
}



##############################################################
# ------ mask all CHM survey grids to vegetation period -----
##############################################################

# ------------------------------------------------------------------------------
# Only survey dates between June - September for each year should be considered
# to avoid CHMs based on Orthophotos collected outside of vegetation period
# as this can create bias in the generated CHMs 
# ------------------------------------------------------------------------------


# load all survey grids and subset them to those areas where survey was in June, July, August or September

survey2017 <- vect("02_dataRaw/DOP_survey_dates/2017/MD_UTM32_BildFluguebersicht_2017.shp")
survey2018 <- vect("02_dataRaw/DOP_survey_dates/2018/MD_UTM32_BildFluguebersicht_2018.shp")
survey2019 <- vect("02_dataRaw/DOP_survey_dates/2019/MD_UTM32_BildFluguebersicht_2019.shp")
survey2020 <- vect("02_dataRaw/DOP_survey_dates/2020/MD_UTM32_BildFluguebersicht_2020.shp")
survey2021 <- vect("02_dataRaw/DOP_survey_dates/2021/MD_UTM32_BildFluguebersicht_2021.shp")
survey2022.2023 <- vect("02_dataRaw/DOP_survey_dates/2022_2023/Befliegungsübersicht_2022_2023.shp")

subset_survey2017 <- survey2017[grep("^\\d{4}(06|07|08|09)", as.character(survey2017$UA)), ]
subset_survey2018 <- survey2018[grep("^\\d{4}/(06|07|08|09)", as.character(survey2018$FLUGDATUM)), ]
subset_survey2019 <- survey2019[grep("^\\d{4}/(06|07|08|09)", as.character(survey2019$FLUGDATUM)), ]
subset_survey2020 <- survey2020[grep("^\\d{4}/(06|07|08|09)", as.character(survey2020$FLUGDATUM)), ]
subset_survey2021 <- survey2021[grep("^\\d{4}/(06|07|08|09)", as.character(survey2021$FLUGDATUM)), ]
subset_survey2022.2023 <- survey2022.2023[grep("^\\d{4}/(06|07|08|09)", as.character(survey2022.2023$flugdatum)), ]


# define the years and file paths

years <- 2017:2023
subsets <- list(
  "2017" = subset_survey2017,
  "2018" = subset_survey2018,
  "2019" = subset_survey2019,
  "2020" = subset_survey2020,
  "2021" = subset_survey2021,
  "2022" = subset_survey2022.2023,
  "2023" = subset_survey2022.2023
)

output_folder <- "03_work/data_processed/CHM/tile_grids/masked_grids" 


# loop through each year

for (year in years) {
  
  grid_file <- sprintf(grid_path_template, year)
  output_file <- file.path(output_folder, paste0("grid_", year, ".gpkg"))
  
  # Load grid file and corresponding subset survey layer for masking
  grid <- vect(grid_file)
  mask_subset <- subsets[[as.character(year)]]
  
  # Mask the grid
  masked_grid <- mask(grid, mask_subset)
  
  
  # Save the masked grid
  writeVector(masked_grid, output_file, filetype = "GPKG", overwrite = TRUE)
  
}


#################################################################################
# --- Extract CHMs for each year per disturbance patch - parallel version ---
#################################################################################


# load all files and paths

patches <- vect("03_work/data_processed/disturbances/patches/disturbance_patches_filtered_0.25ha_25832.gpkg")
years <- 2017:2023

base_path <- "03_work/data_processed/CHM/tile_grids/masked_grids" # where tile grids are stored
patches_file <- "03_work/data_processed/disturbances/patches/disturbance_patches_filtered_0.25ha_25832.gpkg"
buffers_file <- "03_work/data_processed/disturbances/patches/buffer/outer_buffer_area_500.gpkg"
innerbuffer30_file <- "03_work/data_processed/disturbances/patches/buffer/inner_patch_edge_mask.gpkg"
innerbuffer10_file <- "03_work/data_processed/disturbances/patches/buffer/inner_buffer_10.gpkg"
forest_file <- "02_dataRaw/disturbances/germany/forestcover_germany.tif"
dist_all_file <- "02_dataRaw/disturbances/disturbance_year_1986-2020_germany.tif"




# --- define the function to process patches ----

process_patches <- function(patch_indices, patches_file, buffers_file, innerbuffer10_file, innerbuffer30_file, forest_file, dist_all_file, years, base_path) {
  
  require(terra)  # Load terra package in kernel
  
  # load layers locally in each kernel
  patches <- vect(patches_file)  
  buffers <- vect(buffers_file)  
  innerbuffer10 <- vect(innerbuffer10_file)
  innerbuffer30 <- vect(innerbuffer30_file)  
  forest <- rast(forest_file)
  dist_all <- rast(dist_all_file)
  
  error_log <- data.frame(patch_index = integer(), error_message = character(), stringsAsFactors = FALSE)  # Initialize a log for errors
  
  for (i in patch_indices) {
    
    tryCatch({
      
      # get the patch
      buffer <- buffers[i]
      patch <- patches[i]
      print(paste0("Processing patch: ", i))
      
      # loop through years
      for (year in years) {
        print(paste0("Processing year: ", year))
        
        # load grid for the specific year
        grid_file <- file.path(base_path, paste0("grid_", year, ".gpkg"))
        grid_YYYY <- vect(grid_file)
        
        # find tiles intersecting the patch buffer
        files_tiles_patch <- grid_YYYY[is.related(grid_YYYY, buffer, "intersects")]$file_tile
        
        if (length(files_tiles_patch) == 0) {
          next
        }
        
        # merge tiles
        
        if (length(files_tiles_patch) == 1) {
          tiles_patch <- rast(files_tiles_patch)  # load single tile
          crs(tiles_patch) <- "EPSG:25832"
        } else if (length(files_tiles_patch) > 1) {
          tiles_patch <- lapply(files_tiles_patch, rast)
          tiles_patch <- do.call(merge, tiles_patch)  # merge multiple tiles
          crs(tiles_patch) <- "EPSG:25832"
        }
        
        # define patch id for writing the rasters
        patch_id <- as.numeric(patch[1, ][[1]])
        
        # outer buffer - mask to buffer area around patch and exclude non-forested and disturbed pixels  - not used in analysis!
        
        buffer.3035 <- project(buffer, "EPSG:3035")
        
        chm_buffer <- mask(crop(tiles_patch, buffer), buffer)
        forest_cropped <- crop(forest, buffer.3035)
        forest_resampled <- project(forest_cropped, chm_buffer, method = "near") #reproject and resample for masking
        chm_buffer <- mask(chm_buffer, forest_resampled, maskvalues = 0)
        
        dist_all_cropped <- crop(dist_all, buffer.3035)
        dist_all_resampled <- project(dist_all_cropped, chm_buffer, method = "near") #reproject and resample for masking
        chm_buffer <- mask(chm_buffer, dist_all_resampled, inverse = TRUE)
        
        writeRaster(chm_buffer, paste0("03_work/data_processed/CHM/disturbance_patches/outer_buffer/","patch_", patch_id, "_" , year, "_outer_buffer", ".tif"), overwrite = TRUE)
        
        # patch processing
        
        chm_patch <- mask(crop(tiles_patch, patch), patch, touches = F)
        writeRaster(chm_patch,  paste0("03_work/data_processed/CHM/disturbance_patches/patch/", "patch_", patch_id, "_", year, ".tif"), overwrite = TRUE)
        
        # mask to inner buffer 10 m # - not used in analysis but can be trialed to retain more disturbed area with the risk of more comission errors in disturbed area !
        
        chm_patch_inner10 <- mask(chm_patch, innerbuffer10)
        if (any(!is.na(values(chm_patch_inner10)))) {
          writeRaster(chm_patch_inner10, paste0("03_work/data_processed/CHM/disturbance_patches/inner_buffer_10/","patch_", patch_id, "_", year, "_10m_innerbuffer", ".tif"), overwrite = TRUE)
        }
        
        # mask to inner buffer 30 m
        
        chm_patch_inner30 <- mask(chm_patch, innerbuffer30, touches = F)
        #if (any(!is.na(values(chm_patch_inner30)))) { # excludes all raster with NAs pnly
        if (sum(!is.na(values(chm_patch_inner30))) > 700) { # excludes all raster with pixel values covering less than a Landsat scene (due to reprojection, sometimes border pixels are left over)
          writeRaster(chm_patch_inner30, paste0("03_work/data_processed/CHM/disturbance_patches/inner_buffer_30/","patch_", patch_id, "_", year, "_30m_innerbuffer", ".tif"), overwrite = TRUE)
        }
        
        print(paste0("Finished patch ", i, " for year ", year))
      }
    }, error = function(e) {
      # append error to data frame
      error_log <<- rbind(error_log, data.frame(patch_index = i, error_message = e$message, stringsAsFactors = FALSE))
      print(paste0("Error processing patch: ", i, " - ", e$message))
    })
  }
  
  return(error_log)  # return error log
}



# split patches into chunks for parallel processing

n_cores <- 30 # Number of cores to use -> with 8 cores ~ 28 h
patch_indices <- split(1:length(patches), cut(1:length(patches), n_cores, labels = FALSE))


# --- Parallel processing ---

tic("Parallel processing")

cl <- makeCluster(n_cores)

# export necessary variables and functions
clusterExport(cl, c("patch_indices", "years", "base_path", "patches_file", "buffers_file", 
                    "innerbuffer10_file", "innerbuffer30_file", "forest_file", "dist_all_file", "process_patches"))

# load terra on all workers
clusterEvalQ(cl, library(terra))

# process in parallel
result <- parLapply(cl, patch_indices, function(indices) {
  process_patches(
    patch_indices = indices,
    patches_file = patches_file,
    buffers_file = buffers_file,
    innerbuffer10_file = innerbuffer10_file,
    innerbuffer30_file = innerbuffer30_file,
    forest_file = forest_file,
    dist_all_file = dist_all_file,
    years = years,
    base_path = base_path
  )
})

all_errors <- do.call(rbind, result) # 20 patches had errors

stopCluster(cl)

toc()


################# finish CHM extraction ########################################


# calculate amount of patches which are still included: 

folder_path <- "03_work/data_processed/CHM/disturbance_patches/inner_buffer_30/"

# List all files
file_list <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)

# Extract patch_id and year from filenames
file_metadata <- tibble(
  filename = basename(file_list),
  patch_id = str_extract(filename, "(?<=patch_)\\d+"),
  year = str_extract(filename, "_\\d{4}_") %>% str_remove_all("_")
)

patch_summary <- file_metadata %>%
  group_by(patch_id) %>%
  summarise(
    file_count = n(),
    unique_years = n_distinct(year),
    .groups = "drop"
  ) %>%
  mutate(years_match_files = file_count == unique_years)



