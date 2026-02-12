##############################################################################
#
# Prepare auxiliary data layer for the data extraction of recovering patches
#
#############################################################################


# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
#
# This script prepares all necessary data layers we use in the combined
# data extraction og height, environmental and management information for 
# the recovering distrubance patches
#
# Specifically, it:
#   1. Modifies the forest ownership map and adds set aside forest areas.
#   2. Crops and masks the topography layers to Bavaria.
#   3. Crops and masks the climate layers to Bavaria.
#   4. Modifies the soil encoding to generate coarser soil categories,
#       which we use in the model fitting
#
# -------------------------------------------------------------------------


# ---libraries

library(terra)
library(dplyr)
library(stringr)

# set directory

setwd("~/data") # or how you name your directory



#%%%%%%%%%%%%%%%%% (1) management %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# load set asides

# strict forest reserves in Bavaria designated before 1986
reserves <- vect("02_dataRaw/set_aside_forests/selected_set_aside/strict_forest_reserves_pre_1986.gpkg") 
# Berchtesgaden and Bavarian forest National Park
bdg <- vect("02_dataRaw/set_aside_forests/selected_set_aside/bdg_forest.gpkg")
baywald <- vect("02_dataRaw/set_aside_forests/selected_set_aside/forest.bay_wald_forest.gpkg")

setasides <- rbind(reserves, bdg, baywald)


# load ownership map 

ownership <- vect("02_dataRaw/ownership/Forstliche_bersichtskarte_lyrx.shp")

unique(ownership$besitzart)


# simplify ownership df and add ID 

ownership$besitzart[ownership$besitzart == "Sonstiger Staatswald SO"] <- "Staatswald ST" # collapse state forest categories
ownership$mngt_type <- as.integer(factor(ownership$besitzart, levels = unique(ownership$besitzart)))

ownership.df <- as.data.frame(ownership)


#  replace ownership with set asides where it covers the map 

ownership_masked <- mask(ownership, setasides, inverse = TRUE) # mask out set aside areas


# change attribute table setup for setasides to merge with the masked ownership vector layer

setasides_transformed <- setasides[, c("objectid", "besitzart", "name", "aelf", "gemeinde", "waldges")]
names(setasides_transformed) <- c("code", "besitzart", "globalid", "st_area_sh", "st_length_", "mngt_type")
setasides_transformed$mngt_type <- 0 # all set aside forests get mngt ID 0!

ownership_update <- rbind(ownership_masked, setasides_transformed)

writeVector(ownership_update, "03_work/data_processed/management/ownership_update.shp")



#%%%%%%%%%%%%%%%%% (2) topography %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elev <- rast("02_dataRaw/topography/dem.germany.tif")
slope <- rast("02_dataRaw/topography/slope.tif")
aspect <- rast("02_dataRaw/topography/aspect.tif")

bavaria <- vect("02_dataRaw/admin/bayern_Grenze.gpkg")

elev.bavaria <- crop(elev, bavaria)
slope.bavaria <- crop(slope, bavaria)
aspect.bavaria <- crop(aspect, bavaria)

topo.stack <- c(elev.bavaria, slope.bavaria, aspect.bavaria)
names(topo.stack) <- c("elevation", "slope", "aspect")

writeRaster(topo.stack, "03_work/data_processed/topography/topography_bavaria.tif", overwrite = T)



#%%%%%%%%%%%%%%%%% (3) climate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp <- rast("02_dataRaw/dwd/multi_annual_air_temp_mean_1991_2020_17_25832.tiff")
precp <- rast("02_dataRaw/dwd/multi_annual_precipitation_1991-2020_17_25832.tiff")

temp.bavaria <- crop(temp, bavaria)
precp.bavaria <- crop(precp, bavaria)

climate.stack <- c(temp.bavaria, precp.bavaria)
names(climate.stack) <- c("temp", "prec")

writeRaster(climate.stack, "03_work/data_processed/climate/climate_bavaria.tiff", overwrite = T)



#%%%%%%%%%%%%%%%%% (4) soil code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

soil <- vect("02_dataRaw/soil/Shape/uebk25by20230116.shp")

soil_df <- as.data.frame(soil)[, c("LEG_EINH", "LEG_TEXT")]

# table with unique combinations of LEG_EINH and LEG_TEXT
soil_unique <- soil_df %>%
  distinct(LEG_EINH, LEG_TEXT) %>%
  rename(soil = LEG_EINH, soil_type = LEG_TEXT)

write.csv(soil_unique, "03_work/data_processed/codes/soil_code.csv", row.names = FALSE)
#soil_code <- read.csv("03_work/data_processed/codes/soil_code.csv")

# relabel/cluster soil types

soil_code <- soil_code %>%
  mutate(soil_id = soil) %>%
  mutate(soil = as.integer(str_extract(soil_id, "\\d+")),
         soil = as.numeric(soil)) 

write.csv(soil_code, "03_work/data_processed/codes/soil_code_clustered.csv", row.names = FALSE)
