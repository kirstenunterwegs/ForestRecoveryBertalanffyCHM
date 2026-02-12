################################################################################
#
# Extract recovery data from CHM & build modelling datasets
#
################################################################################


# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
#
# This script extract canopy height recovery metrics for disturbance patches
# from CHM time series, links them to environmental and management covariates,
# and aggregates everything to a Landsat-resolution modelling dataset.
#
# Specifically, it:
#   1. Links CHM patch raster to disturbance information,
#       management, topography, climate, soil, and species layers.
#   2. Extracts multi-year CHM height trajectories per disturbance patch,
#      adds patch size and surrounding forest structure metrics.
#   3. Aggregates CHM and covariates to Landsat resolution by disturbance patch
#      and year, and filters out edge / low-coverage pixels.
#   4. Cleans missing or inconsistent covariates (disturbance, climate, height,
#      buffer metrics, patch size) and restricts to high-severity disturbances.
#   5. Splits data into first/last CHM observations per pixel. 
#   6. Creates spatially thinned subsets (250–1000 m grids) to reduce spatial
#      autocorrelation for model fitting.
#   7. Check correlation between different predictor variables
#
# Output:
#   - Pixel- and Landsat-scale recovery tables (RData / RDS)
#   - Spatially thinned subsets for modelling
#
# Notes:
#   - Parallel extraction is used; runtime depends strongly on amount of
#     cores used
# -------------------------------------------------------------------------


# ---libraries

library(terra)
library(data.table)
library(dplyr)
library(stringr)
library(parallel)
library(tictoc)
library(ggplot2)
library(corrplot)

# set directory

setwd("~/data") # or how you name your directory


# ------------------------------------------------------------
# prepare directories and meta data frames for the extraction
# ------------------------------------------------------------

# define directories to load data in extraction

dir_dist_patches <- "03_work/data_processed/CHM/disturbance_patches/inner_buffer_30/" # removed inner edge pixel of disturbance patch
dir_dist_size <- "03_work/data_processed/CHM/disturbance_patches/patch/" # original patches for patch size calculation
dir_dist_buffer <- "03_work/data_processed/CHM/disturbance_patches/outer_buffer/" # 500 m buffer around disturbance
dir_dist <- "02_dataRaw/disturbances/"
dir_mngt <- "03_work/data_processed/management/ownership_update.shp"
dir_topo <- "03_work/data_processed/topography/topography_bavaria.tif"
dir_spec <- "02_dataRaw/tree_species/species_class_sum_INT1U.tif"
dir_clim <- "03_work/data_processed/climate/climate_bavaria.tiff"
dir_eco <- "02_dataRaw/ecoregion/siteconditions_raster.tif"
dir_soil <- "02_dataRaw/soil/Shape/uebk25by20230116.shp"


# --- List all files ---

# -- Patches -- (inner buffered patches - 30m)

#file_list <- list.files(dir_dist_patches, pattern = "\\.tif$", full.names = TRUE)
#saveRDS(file_list, "03_work/data_processed/CHM/disturbance_patches/list_inner_buffer_30.rds" )
file_list <- readRDS("03_work/data_processed/CHM/disturbance_patches/list_inner_buffer_30.rds")


# get metadta to load patch and corresponding buffer area
patch_metadata <- data.table(
  path = file_list,
  patch_id = str_extract(basename(file_list), "(?<=patch_)\\d+"),
  year = str_extract(basename(file_list), "_\\d{4}_") %>% str_remove_all("_")
)


# -- Patches -- original size
 
# file_list_patch <- list.files(dir_dist_size, pattern = "\\.tif$", full.names = TRUE)
# saveRDS(file_list_patch, "03_work/data_processed/CHM/disturbance_patches/file_list_patch.rds" )
file_list_patch <- readRDS("03_work/data_processed/CHM/disturbance_patches/file_list_patch.rds")


patch_size_metadata <- data.table(
  path = file_list_patch,
  patch_id = str_extract(basename(file_list_patch), "(?<=patch_)\\d+"),
  year = str_extract(basename(file_list_patch), "(?<=_)\\d{4}(?=\\.tif)")  
  )


# -- Buffer --

#file_list_buffer <- list.files(dir_dist_buffer, pattern = "\\.tif$", full.names = TRUE)
#saveRDS(file_list_buffer, "03_work/data_processed/CHM/disturbance_patches/list_outer_buffer.rds" )
file_list_buffer <- readRDS("03_work/data_processed/CHM/disturbance_patches/list_outer_buffer.rds")


buffer_metadata <- data.table(
  path = file_list_buffer,
  patch_id = str_extract(basename(file_list_buffer), "(?<=patch_)\\d+"),
  year = str_extract(basename(file_list_buffer), "_\\d{4}_") %>% str_remove_all("_")
)




# ----------------------------------
# define extraction function
# ----------------------------------



extract.recovery = function(patch_metadata, buffer_metadata, patch_size_metadata, id, dir_dist, dir_mngt, dir_topo, dir_spec, dir_clim, dir_eco, dir_soil){
  
  # load chms per patch for all years
  print( paste0("ID:" , id))
  patch_files <- patch_metadata[patch_id == id, path]
  chms = rast(patch_files)
  names_chms = as.numeric(str_extract(patch_files, "(?<=patch_)\\d+_\\d{4}") %>% str_extract("(?<=_)[^_]+$")) # change name to patch_id, years
  names(chms) = names_chms
  
  # Add disturbance information
  
  files_dist = list.files(dir_dist, pattern = ".tif", full.names = T)
  dist = rast(files_dist)
  names_dist =  gsub("^.*/([^_]+)_([^_]+).*\\.tif$", "\\1_\\2", files_dist)
  names(dist) = names_dist
  
  ext_chms = as.polygons(ext(chms), crs = terra::crs(chms)) # get extent of disturbance patch
  
  ext_3035 = project(ext_chms, terra::crs(dist))
  dist_sub = crop(dist, ext_3035)
  dist_sub$IDlandsat = dist_sub$agent_classes
  values(dist_sub$IDlandsat) = 1:ncell(dist_sub)
  dist_sub = project(dist_sub, subset(chms,1), method = "near")
   
  # remove areas that are NA in the chms
  dist_sub = ifel(all(is.na(chms)), NA, dist_sub)
  
  
   # add management - ownership information
  management <-  vect(dir_mngt)
  mngt_sub <- crop(management, ext_chms)
  
  # if mngt_sub has no geometry, assign 6 as mngt: 
  if (nrow(mngt_sub) == 0) {
    mngt_rast <- chms[[1]]
    mngt_rast[!is.na(chms[[1]])] <- 6
    names(mngt_rast) <- "mngt_type"
  } else { 
  mngt_rast <- rasterize(mngt_sub, dist_sub, field = "mngt_type", fun = "max")
  mngt_rast = ifel(all(is.na(chms)), NA, mngt_rast)
  
    # fill up missing management information (ownership layer doesn't cover all forest/disturbed areas)
    unique_mngt <- unique(mngt_rast$mngt_type)
    if (length(unique_mngt$mngt_type) == 1) {       # when there is only one mngt type present, assign this one
      mngt_rast[!is.na(chms[[1]]) & is.na(mngt_rast)] <- unique_mngt
    } else {                              # when there is more than one mngt type, assign 6 as managed without ownership specifics
      mngt_rast[!is.na(chms[[1]]) & is.na(mngt_rast)] <- 6
    }
  } 
    
  # add topography 
  topo <- rast(dir_topo)
  topo_sub <- crop(topo, ext_chms)
  topo_sub = resample(topo_sub, subset(chms,1), method = "near")
  topo_sub = ifel(all(is.na(chms)), NA, topo_sub)
  
  # add tree species/ forest type
  spec <- rast(dir_spec)
  spec_sub <- crop(spec, ext_3035)
  spec_sub = project(spec_sub, subset(chms,1), method = "near")
  spec_sub = ifel(all(is.na(chms)), NA, spec_sub)
  names(spec_sub) <- "species"
  
  # add climate
  climate <- rast(dir_clim)
  climate_sub <- crop(climate, ext_chms)
  climate_sub = resample(climate_sub, subset(chms,1), method = "near")
  climate_sub = ifel(all(is.na(chms)), NA, climate_sub)
  
  ecoregion <- rast(dir_eco)
  eco_sub <- crop(ecoregion, ext_chms)
  eco_sub = resample(eco_sub, subset(chms,1), method = "near")
  eco_sub = ifel(all(is.na(chms)), NA, eco_sub)
  names(eco_sub) <- "ecoregion"
  
  # add soil
  soil <- vect(dir_soil)
  soil_sub <- crop(soil, ext_chms)
  soil_sub <- rasterize(soil_sub, dist_sub, field = "LEG_EINH", fun = "max")
  soil_sub = ifel(all(is.na(chms)), NA, soil_sub)
  names(soil_sub) <- "soil"
  
  # convert to data.table
  chms_dist = c(chms, dist_sub, mngt_rast, topo_sub, spec_sub, climate_sub, eco_sub, soil_sub)
  recovery = as.data.table(chms_dist, xy = T, cells = T)
  recovery = melt(recovery, id.vars = c("cell","x","y", names(dist_sub), "mngt_type", names(topo_sub), "species", names(climate_sub), "ecoregion",  "soil"), 
                  variable.name = "chm_year", value.name = "h")
  recovery$chm_year = as.integer(as.character(recovery$chm_year))
  recovery$patch = id 
  
  # add dist patch size
  patch_files_full <- patch_size_metadata[patch_id == id, path]
  full_patch = rast(patch_files_full[1]) # load one whole patch
  patch_ha <- (sum(!is.na(values(full_patch))))/10000 # calculate the area of the patch
  recovery$patch_ha =patch_ha
  
  # load disturbance patch buffer area (already masked to fores area and by disturbed area)
  
  buffer_file <-  buffer_metadata[patch_id == id][which.max(as.numeric(year)), path] # path to most recent CHM of surroundings
  buffer = rast(buffer_file)
  buffer <- clamp(buffer, lower = 0, upper = 50, values = TRUE) # set neg val to 0 and very high val to 50 m
  
  
  #--- create IDlandsat to aggregate values per pixel before getting global stats
  
        ext_buffer = as.polygons(ext(buffer), crs = terra::crs(buffer))
        ext_dist = project(ext_buffer, terra::crs(dist))
        dist_sub = crop(dist, ext_dist)
        dist_sub$IDlandsat = dist_sub$agent_classes
        values(dist_sub$IDlandsat) = 1:ncell(dist_sub)
        dist_sub = project(dist_sub, buffer, method = "near")

  # aggregate CHM metrics per landsatID in buffer
        
     buffer_landsat_agg <- zonal(buffer, dist_sub$IDlandsat, fun = "mean", as.raster=T)

     recovery[, `:=`(
       ha_forest_buffer = global(buffer_landsat_agg, fun = "notNA")[1, 1]/ 10000 , #sum(!is.na(values(buffer))), # 1m2 res of CHM pixel
       b_mean_height = global(buffer_landsat_agg, fun = "mean", na.rm = TRUE)[1, 1],
       b_perc_75 = global(buffer_landsat_agg, fun = function(x) quantile(x, probs = 0.75, na.rm = TRUE))[1, 1],
       b_perc_95 = global(buffer_landsat_agg, fun = function(x) quantile(x, probs = 0.95, na.rm = TRUE))[1, 1],
       b_perc_99 = global(buffer_landsat_agg, fun = function(x) quantile(x, probs = 0.99, na.rm = TRUE))[1, 1],
       b_max_height = global(buffer_landsat_agg, fun = "max", na.rm = TRUE)[1, 1],
       b_sd_height = global(buffer_landsat_agg, fun = "sd", na.rm = TRUE)[1, 1]
     )]
  
  return(recovery)
}



# ---------------------------------------
# parallel extraction of recovery data
# ---------------------------------------

# ---- Paralleled function 

parallel_extract_recovery <- function(id, patch_metadata, buffer_metadata, patch_size_metadata, dir_dist, dir_mngt, dir_topo, dir_spec, dir_clim, dir_eco, dir_soil) {
  tryCatch({
    # Run the extraction function
    extract.recovery(patch_metadata, buffer_metadata, patch_size_metadata, id, dir_dist, dir_mngt, dir_topo, dir_spec, dir_clim, dir_eco, dir_soil)
  }, error = function(e) {
    # Log the error message with the corresponding patch ID
    return(data.table(patch_id = id, error_message = e$message))
  })
}


# Get IDs to process

ids <- unique(patch_metadata$patch_id) # get all patch IDs
ids_b <- unique(buffer_metadata$patch_id) # get all patch IDs where we extracted buffer area
ids_size <- unique(patch_size_metadata$patch_id)

# there are more IDs in buffer than in patch_metadata, as those ones are the masked ones, so small/thin patches are excluded

ids <- intersect(ids, ids_b) 
ids <- intersect(ids, ids_size)

# define number of cores and split ids into different thread categories

n_cores <- 40 # Number of cores to use -> with 8 cores ~ 28 h  / with 40 cores ~ 8 hours
id_indices <- split(ids, cut(1:length(ids), n_cores, labels = FALSE))


tic("Start data extraction for recovery dt")

# -----------------------
# --- Run in parallel ---
# -----------------------

cl <- makeCluster(n_cores)

clusterExport(cl, c("extract.recovery", "patch_metadata", "buffer_metadata", "patch_size_metadata", "dir_dist", "parallel_extract_recovery", "dir_mngt", "dir_topo", "dir_spec", "dir_clim", "dir_eco", "dir_soil"))

clusterEvalQ(cl,{library(terra) # Ensure packages is loaded on all workers
  library(data.table) 
  library(dplyr) 
  library(stringr)
}) 

recovery_list <- parLapply(cl, ids, function(id) {
  parallel_extract_recovery(id, patch_metadata, buffer_metadata, patch_size_metadata, dir_dist, dir_mngt, dir_topo, dir_spec, dir_clim, dir_eco, dir_soil)
})
stopCluster(cl)

toc()

# Separate successful results from errors
results <- Filter(function(x) !("error_message" %in% names(x)), recovery_list)
error_logs <- Filter(function(x) "error_message" %in% names(x), recovery_list)


# Bind & store results
recovery.dt <- rbindlist(results, fill = TRUE)
save(recovery.dt, file = "03_work/data_processed/analysis.dt/recovery_all.dt.RData")

# Save errors for debugging
if (length(error_logs) > 0) {
  error_log_dt <- rbindlist(error_logs, fill = TRUE)
  fwrite(error_log_dt, "03_work/data_processed/analysis.dt/recovery_all_errors.csv")
}

#remove(recovery.dt) 


# %%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Within the data repository we provide a version of the full recovery.dt, 
# but without x,y coordinates as we are not allowed to share spatially 
# explicit forest ownership data.
# But you can reproduce the filtering and subsetting steps which follow 
# without the coordinates.
#
# We provide the spatial subsets produced later on within the repository.
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- data table provided in the repository: 
recovery.dt <- recovery.dt[, !c("x", "y")]
save(recovery.dt, file = "03_work/data_processed/analysis.dt/recovery_all.dt_noXY.RData")
# ---

# ---- this is a very large data table, it takes several minutes to load! -----

load("03_work/data_processed/analysis.dt/recovery_all.dt_noXY.RData")

# ---------------------------------------------------------------------
# Aggregate recovery dt to Landsat pixel resolution of disturbance map
# ---------------------------------------------------------------------

#error_log_dt <- read.csv("03_work/data_processed/analysis.dt/recovery_all_errors.csv")
#load("03_work/data_processed/analysis.dt/recovery_all.dt.RData")


# --- create data set with first and last survey per pixel ---

# subset to first and last chm survey years
recovery.dt[, `:=`(chm_minyear = min(chm_year), chm_maxyear = max(chm_year)),by = .(cell, patch)]
recovery_fullperiod = recovery.dt[chm_year == chm_minyear | chm_year == chm_maxyear]
recovery.dt$chm_maxyear = NULL
recovery.dt$chm_minyear = NULL
recovery_fullperiod$chm_maxyear = NULL
recovery_fullperiod$chm_minyear = NULL

save(recovery_fullperiod, file = "03_work/data_processed/analysis.dt/recovery_all_fullperiod.RData")


# --- aggregate to Landsat pixel resolution ---

#load("03_work/data_processed/analysis.dt/recovery_all_fullperiod.RData")

tic("start landsat aggregation")

recovery_landsat = recovery_fullperiod[,list(
  
    x = mean(x)
  , y = mean(y)
  
  # height metrics
  , h_mean = mean(ifelse(is.finite(h), h, NA), na.rm = TRUE)
  , h_sd = sd(ifelse(is.finite(h), h, NA), na.rm = TRUE)
  , h_g75 = quantile(ifelse(is.finite(h), h, NA), probs = 0.75, na.rm = TRUE)
  , h_q95 = quantile(ifelse(is.finite(h), h, NA), probs = 0.95, na.rm = TRUE)
  , h_q99 = quantile(ifelse(is.finite(h), h, NA), probs = 0.99, na.rm = TRUE)
  
  # climate
  , prec = mean(ifelse(is.finite(prec), prec, NA), na.rm = TRUE)
  , temp = mean(ifelse(is.finite(temp), temp, NA), na.rm = TRUE)
  # topography
  , elevation = mean(ifelse(is.finite(elevation), elevation, NA), na.rm = TRUE)
  , slope = mean(ifelse(is.finite(slope), slope, NA), na.rm = TRUE)
  , aspect = mean(ifelse(is.finite(aspect), aspect, NA), na.rm = TRUE)
  
  # rank occurrence of different classes per Landsat pixel and chose most prevalent:
  # (for categorical variables)
  , mngt_type = names(which.max(table(mngt_type))) 
  , ecoregion = names(which.max(table(ecoregion)))
  , soil = names(which.max(table(soil)))
  , species = names(which.max(table(species)))
  , ncells = .N  # amount of entries which go into aggregation, so I can filter out border & overlap pixel areas
  
  # Landsat pixel level information
  , agent_classes = unique(agent_classes)
  , disturbance_severity = unique(disturbance_severity)
  , disturbance_year = unique(disturbance_year)
  , patch_ha = unique(patch_ha)
  
  # Surrounding height in buffer information
  , ha_forest_buffer = unique(ha_forest_buffer)
  , b_mean_height = unique(b_mean_height)
  , b_perc_75 = unique(b_perc_75)
  , b_perc_95 = unique(b_perc_95)
  , b_perc_99 = unique(b_perc_99)
  , b_max_height = unique(b_max_height)
  , b_sd_height = unique(b_sd_height)
  
), by = .(IDlandsat, patch, chm_year)] 

toc()

save(recovery_landsat, file = "03_work/data_processed/analysis.dt/recovery_all_landsat_hmean_hmax.RData")

load("03_work/data_processed/analysis.dt/recovery_all_landsat_hmean_hmax.RData")

# -- filter: remove all Landsat pixel with less than 500 entries per pixel to remove boundary and reprojection inaccuracies

recovery_landsat_filtered <- recovery_landsat[ncells > 800]

# fill up species with no species information with dummy code 99 for "other species"

recovery_landsat_filtered <- recovery_landsat_filtered[, species := as.numeric(species)]
recovery_landsat_filtered <- recovery_landsat_filtered[is.na(species), species := 99]

recovery_landsat_filtered <- recovery_landsat_filtered[, mngt_type := as.numeric(mngt_type)]
recovery_landsat_filtered <- recovery_landsat_filtered[, ecoregion := as.numeric(ecoregion)]


# ---------------------------------------
# ---- check NAs in various columns -----
# ---------------------------------------

# --- missing disturbance information 

na_disturbance_year <- recovery_landsat_filtered %>%  # 663 patches with missing disturbance information, :2165 pixels
  filter(is.na(disturbance_year))

  # delete these entries
  # pattern of missing disturbance info not specific, but random in the landscape
  recovery_landsat_filtered <- na.omit(recovery_landsat_filtered, cols = "disturbance_year")

    
# --- missing climate data
  
  na_temp <- recovery_landsat_filtered %>%  # 51 patches without or missing climate information
  filter(is.na(temp))
  
  # patches crossing the border of climate data layer, so we just fill up those patches with information from the rest of the patches
  ids <- unique(na_temp$patch)
  
  # Calculate mean temp and prec per patch, excluding NAs where there is climate information in the patch (21 of 51 patches)
  mean_climate <- recovery_landsat_filtered %>%
    filter(patch %in% ids) %>%
    filter(is.finite(temp) & is.finite(prec)) %>%  # Remove rows with NaN or NA in temp or prec
    group_by(patch) %>%
    summarise(
      mean_temp = mean(temp[is.finite(temp)], na.rm = TRUE),
      mean_prec = mean(prec[is.finite(prec)], na.rm = TRUE)
    )

  recovery_landsat_filtered <- recovery_landsat_filtered %>% #  join with the original dataframe
    left_join(mean_climate, by = "patch")
  
  # Assign the calculated means to rows where temp or prec is NA
  recovery_landsat_filtered <- recovery_landsat_filtered %>%
    mutate(
      temp = ifelse(is.na(temp), mean_temp, temp),
      prec = ifelse(is.na(prec), mean_prec, prec)
    ) %>%
    select(-mean_temp, -mean_prec)  # Remove the temporary mean columns
  
  # left with 460 pixels without climate data
  # Delete those data entries
  recovery_landsat_filtered <- na.omit(recovery_landsat_filtered, cols = "temp")
  

# --- missing height entries  
  
  na_h <- recovery_landsat_filtered %>%  # 148 pixels with no height data, as first or last CHM doesn't cover whole patch area (10 patches)
  filter(is.na(h_q95))
  
  recovery_landsat_filtered <- na.omit(recovery_landsat_filtered, cols = "h_q95")

  # check missing surrounding buffer entries
  na_b_mean_height  <- recovery_landsat_filtered %>%  # 459 pixels without buffer information (37 patches)
    filter(is.na(b_mean_height))

  recovery_landsat_filtered <- na.omit(recovery_landsat_filtered, cols = "b_mean_height")

  # identify pixels with patch size = 0, as those had problems with patch size extraction: 
  # Get unique patch IDs where patch ha is below minimum patch size
  
  patch_ids_to_update <- unique(recovery_landsat_filtered[patch_ha < 0.25, patch])
  
  # Loop through patch IDs and update patch_ha
  
  for (id in patch_ids_to_update) {
    
    print(id)
    
    # Get raster file path
    patch_file_path <- patch_size_metadata[patch_id == id, path]
    
    if (length(patch_file_path) >= 1 && file.exists(patch_file_path[1])) {
      
      # Load the raster
      full_patch <- tryCatch({
        rast(patch_file_path[2])  # 2 in this case, as 1 did not exist for original extraction!
      }, error = function(e) {
        message(paste("Could not load raster for patch", id))
        return(NULL)
      })
      
      if (!is.null(full_patch)) {
        # Calculate patch size in ha
        patch_ha_new <- sum(!is.na(values(full_patch))) / 10000
        
        # Update patch_ha in your main data.table
        recovery_landsat_filtered[patch == id & patch_ha < 0.25, patch_ha := patch_ha_new]
      }
    }
  }
  

# --- filter out observations with disturbance severity < 0.75 (remove less severe disturbances)
  
recovery_landsat_filtered_sev <- recovery_landsat_filtered[disturbance_severity > 0.75]


save(recovery_landsat_filtered_sev, file = "03_work/data_processed/analysis.dt/recovery_all_landsat_filtered_sev0.75__hmean_hmax.RData")


# -----------------------------------------------
# separate into first and last observation year
# -----------------------------------------------

load("03_work/data_processed/analysis.dt/recovery_all_landsat_filtered_sev0.75__hmean_hmax.RData")


# Separate into two data tables based on the chm_year (earlier vs. later)

recovery_landsat_filtered_sev_firstyr <- recovery_landsat_filtered_sev %>%
  group_by(patch, IDlandsat) %>%
  filter(chm_year == min(chm_year)) %>%  
  ungroup()

recovery_landsat_filtered_sev_lastyr <- recovery_landsat_filtered_sev %>%
  group_by(patch, IDlandsat) %>%
  filter(chm_year == max(chm_year)) %>% 
  ungroup()

save(recovery_landsat_filtered_sev_firstyr, file = "03_work/data_processed/analysis.dt/recovery_all_landsat_filtered_sev0.75__hmean_hmax_firstyr.RData")
save(recovery_landsat_filtered_sev_lastyr, file = "03_work/data_processed/analysis.dt/recovery_all_landsat_filtered_sev0.75__hmean_hmax_lastyr.RData")


# ----------------------------------------------------------
# further data cleaning and processing pre modelling 
# ----------------------------------------------------------

setDT(recovery_landsat_filtered_sev_firstyr)

# remove species group 99 (other, where we do not have any species classification by Blickensdörfer)

recovery_landsat_filtered_sev_firstyr = recovery_landsat_filtered_sev_firstyr[species != 99] # remove 24440 obs

# remove all observation from the disturbance year and the year after, to avoid error in disturbance year classification

recovery_landsat_filtered_sev_firstyr[, t := chm_year - disturbance_year]
recovery_landsat_filtered_sev_firstyr[t < 0, t := 250]
recovery_landsat_filtered_sev_firstyr = recovery_landsat_filtered_sev_firstyr[t > 1] #remove 12002 observations

dt_sev0.75_clean_hmean_hmax = recovery_landsat_filtered_sev_firstyr


saveRDS(dt_sev0.75_clean_hmean_hmax, "03_work/data_processed/analysis.dt/dt_sev0.75_clean_hmean_hmax.rds")
#dt_sev0.75_clean_hmean_hmax = readRDS("03_work/data_processed/analysis.dt/dt_sev0.75_clean_hmean_hmax.rds")


# --------------------------------------------------------------
# create spatial subsets to account for spatial autocorrelation
# --------------------------------------------------------------

# convert centroid to point layer for random draw of observation per grid cell

obs_points <- vect(dt_sev0.75_clean_hmean_hmax, geom = c("x", "y"), crs = "EPSG:25832")


# ---- spatial sub sampling for 1000 m grid

# raster grid covering the extent
r <- rast(ext(obs_points), resolution = 1000, crs = "EPSG:25832")

# give each 4x4 km cell an ID
r[] <- 1:ncell(r)

# extract grid ID for each point
grid_ids = terra::extract(r, obs_points)
dt_sev0.75_clean_hmean_hmax[, grid_id := grid_ids[[2]]]

hist(dt_sev0.75_clean_hmean_hmax$grid_id)

# Sample one random obs per grid cell
set.seed(789)
dt_sev0.75_sampled <- dt_sev0.75_clean_hmean_hmax[, .SD[sample(.N, 1)], by = grid_id]

hist(dt_sev0.75_sampled$mngt_type)
dt_sev0.75_sampled %>% group_by(mngt_type) %>% count()

saveRDS(dt_sev0.75_sampled, "03_work/data_processed/analysis.dt/spatial_sub1000.rds")



# ---- spatial sub sampling for 500 m grid

r <- rast(ext(obs_points), resolution = 500, crs = "EPSG:25832")

# give each 4x4 km cell an ID
r[] <- 1:ncell(r)

# extract grid ID for each point
grid_ids = terra::extract(r, obs_points)
dt_sev0.75_clean_hmean_hmax[, grid_id := grid_ids[[2]]]

hist(dt_sev0.75_clean_hmean_hmax$grid_id)

# Sample one random obs per grid cell
set.seed(789)
dt_sev0.75_sampled <- dt_sev0.75_clean_hmean_hmax[, .SD[sample(.N, 1)], by = grid_id]

hist(dt_sev0.75_clean_hmean_hmax$mngt_type)
hist(dt_sev0.75_sampled$mngt_type)

dt_sev0.75_sampled %>% group_by(mngt_type) %>% count()


saveRDS(dt_sev0.75_sampled, "03_work/data_processed/analysis.dt/spatial_sub500.rds")



# ---- spatial sub sampling for 250 m grid

r <- rast(ext(obs_points), resolution = 250, crs = "EPSG:25832")

# give each 4x4 km cell an ID
r[] <- 1:ncell(r)

# extract grid ID for each point
grid_ids = terra::extract(r, obs_points)
dt_sev0.75_clean_hmean_hmax[, grid_id := grid_ids[[2]]]

hist(dt_sev0.75_clean_hmean_hmax$grid_id)

# Sample one random obs per grid cell
set.seed(789)
#set.seed(437) # second subset for trial
dt_sev0.75_sampled <- dt_sev0.75_clean_hmean_hmax[, .SD[sample(.N, 1)], by = grid_id]

hist(dt_sev0.75_clean_hmean_hmax$mngt_type)
hist(dt_sev0.75_sampled$mngt_type)

dt_sev0.75_sampled %>% group_by(mngt_type) %>% count()


saveRDS(dt_sev0.75_sampled, "03_work/data_processed/analysis.dt/spatial_sub250.rds")


# -------------------------------------------------
# ---  explore correlations between variables:
# -------------------------------------------------


recovery_dt = as.data.table(readRDS("03_work/data_processed/analysis.dt/spatial_sub250.rds"))
recovery_dt = recovery_landsat_filtered_sev_firstyr

recovery_landsat_filtered_sev_firstyr_sub <- recovery_dt[, c("prec", "temp", "elevation", "slope", "aspect", "patch_ha")]

recovery_landsat_filtered_sev_firstyr_sub <- na.omit(recovery_landsat_filtered_sev_firstyr_sub, cols = "ecoregion")
recovery_landsat_filtered_sev_firstyr_sub[] <- lapply(recovery_landsat_filtered_sev_firstyr_sub, function(x) {
  if (is.factor(x) || is.character(x)) {
    as.numeric(factor(x))  
  } else {
    as.numeric(x)  
  }
})

str(recovery_landsat_filtered_sev_firstyr_sub)

cor(recovery_landsat_filtered_sev_firstyr_sub)
cor_matrix <- cor(recovery_landsat_filtered_sev_firstyr_sub)
corrplot(cor_matrix, method = "circle") 
corrplot(cor_matrix,
         method = "circle",
         type = "upper",     
         addCoef.col = "black",  
         number.cex = 0.7,       
         tl.col = "black",       
         tl.cex = 0.8)   

cor_plot = corrplot(cor_matrix,
                    method = "circle",
                    type = "upper",     
                    addCoef.col = "black",  
                    number.cex = 0.7,       
                    tl.col = "black",       
                    tl.cex = 0.8)  


tiff("03_work/data_processed/results/supporting_information/corr_plot_sub250.tif", width = 12, height = 6, units = "in", res = 300)
corrplot(cor_matrix,
         method = "circle",
         type = "upper",
         addCoef.col = "black",
         number.cex = 1.2,    
         tl.col = "black",
         tl.cex = 1.1)        

dev.off()


