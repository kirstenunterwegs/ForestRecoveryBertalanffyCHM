################################################################################
# Model fit of the modified Bertalanffy model with 250 m grid subset 
# Student error distribution
################################################################################


###### load libraries

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(stringr))

###### set dir

setwd("~/data") # or how you name your directory

# path = dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(path)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### load data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

sub250 = as.data.table(readRDS("03_work/data_processed/analysis.dt/spatial_sub250.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### prepare data for model fit ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


sub250[, slope_std := scale(slope)]
sub250[, aspect_std := scale(aspect)]
sub250[, elevation_std := scale(elevation)]

sub250[, t_alt := t]
sub250[, t_alt := fifelse(t == 250, 35, t)]
sub250[, log10t := log10(t_alt)] 
sub250[, patch_ha_log := log(patch_ha)]
sub250[, patch_ha_log_std := scale(patch_ha_log)]


# --- rename management and species

# note: in the process of writing the manuscript, we relabeled some management types (communal forest = municipal forest)

# Store original values
sub250[, mngt_id := mngt_type]
sub250[, species_id := species]

# relabel
sub250[, mngt_type := factor(mngt_type,
                             levels = 0:6,
                             labels = c("Set aside", "Large private", "Small private",
                                        "State forest", "Communal forest", "Federal forest", "Other"))]
species_map <- c(`2` = "Birch", `3` = "Beech", `4` = "Douglas fir", `5` = "Oak", 
                 `6` = "Alder", `8` = "Spruce", `9` = "Pine", `10` = "Larch", 
                 `14` = "Fir", `16` = "ODH", `17` = "ODL", `99` = "Other")
sub250[, species := species_map[as.character(species)]]
sub250[, species := factor(species)]

# --- re-classify soil to a higher level grouping by removing fine classification

sub250[, soil_id := soil]
sub250[, soil := as.numeric(str_extract(soil_id, "\\d+"))]


#%%%%%%%%%%%%%%%%%%%%%%%%%#
##### model fit####
#%%%%%%%%%%%%%%%%%%%%%%%%%#

cat("==== Starting model fit ====\n")

# Record start time
start_time <- Sys.time()

#fitting the model inside tryCatch to handle failures 
fit_success <- FALSE

try({
  fit_fe_sub250_cmdLOG_student14 = brm(  
  bf(
    h_mean ~ exp(hmaxstarLOG) * (1 - exp(-exp(gLOG) * t))^3 + exp(hOffsetLOG),
    gLOG ~ 1 + slope_std + elevation_std + aspect_std + patch_ha_log_std + (1 | mngt_type) + (1 | species) + (1 | mngt_type : species) + (1 | soil),       
    hOffsetLOG ~ 1 + slope_std + elevation_std + aspect_std + patch_ha_log_std +  (1 | mngt_type) + (1 | species) + (1 | mngt_type : species) + (1 | soil), 
    hmaxstarLOG ~ 1 + slope_std + elevation_std + aspect_std + patch_ha_log_std +  (1 | mngt_type) + (1 | species) + (1 | mngt_type : species) + (1 | soil),
    # Note: as the student distribution cannot shift its mode, we do not model the change in distribution here across time
    nl = TRUE
  ),
  data = sub250, 
  family = "student",
  backend = "cmdstanr", 
  prior = c(
    # fixed
    #g
    prior(normal(-3, 1), nlpar = gLOG, class = "b")  
    ,prior(normal(0, 0.5), nlpar = gLOG, class = "b", coef = "slope_std") # slope 
    ,prior(normal(0, 0.5), nlpar = gLOG, class = "b", coef = "elevation_std") # temperature 
    ,prior(normal(0, 0.5), nlpar = gLOG, class = "b", coef = "aspect_std") # aspect
    ,prior(normal(0, 0.5), nlpar = gLOG, class = "b", coef = "patch_ha_log_std") # patch size
    # hOffsetLOG
    ,prior(normal(1,1), nlpar = hOffsetLOG, class = "b") 
    ,prior(normal(0, 0.5), nlpar = hOffsetLOG, class = "b", coef = "slope_std") # slope 
    ,prior(normal(0, 0.5), nlpar = hOffsetLOG, class = "b", coef = "elevation_std") # temperature
    ,prior(normal(0, 0.5), nlpar = hOffsetLOG, class = "b", coef = "aspect_std") # aspect
    ,prior(normal(0, 0.5), nlpar = hOffsetLOG, class = "b", coef = "patch_ha_log_std") # patch size
    # hmax
    ,prior(normal(3,1), nlpar = hmaxstarLOG, class = "b") 
    ,prior(normal(0, 0.5), nlpar = hmaxstarLOG, class = "b", coef = "slope_std") # slope 
    ,prior(normal(0, 0.5), nlpar = hmaxstarLOG, class = "b", coef = "elevation_std") # temperature
    ,prior(normal(0, 0.5), nlpar = hmaxstarLOG, class = "b", coef = "aspect_std") # aspect
    ,prior(normal(0, 0.5), nlpar = hmaxstarLOG, class = "b", coef = "patch_ha_log_std") # patch size
    # sd
    ,prior(normal(0, 0.5), lb = 0, nlpar = gLOG, class = "sd")   # --- g
    ,prior(normal(0, 0.5), nlpar = gLOG, group = mngt_type,lb = 0, class = "sd") 
    ,prior(normal(0, 0.5), nlpar = gLOG, group = species, lb = 0, class = "sd")
    ,prior(normal(0, 0.5), nlpar = gLOG, group = mngt_type:species, lb = 0, class = "sd")
    ,prior(normal(0, 0.5), nlpar = gLOG, group = soil, lb = 0, class = "sd")
    ,prior(normal(0, 0.5), lb = 0, nlpar = hOffsetLOG, class = "sd")  # --- hOffset
    ,prior(normal(0, 0.5), nlpar = hOffsetLOG, group = mngt_type, lb = 0, class = "sd") 
    ,prior(normal(0, 0.5), nlpar = hOffsetLOG, group = species, lb = 0,class = "sd") 
    ,prior(normal(0, 0.5), nlpar = hOffsetLOG, group = mngt_type:species, lb = 0,class = "sd") 
    ,prior(normal(0, 0.5), nlpar = hOffsetLOG, group = soil, lb = 0,class = "sd") 
    ,prior(normal(0, 0.5), lb = 0, nlpar = hmaxstarLOG, class = "sd")  # --- hmax
    ,prior(normal(0, 0.5), nlpar = hmaxstarLOG, group = species, lb = 0,class = "sd") 
    ,prior(normal(0, 0.5), nlpar = hmaxstarLOG, group = mngt_type, lb = 0,class = "sd") 
    ,prior(normal(0, 0.5), nlpar = hmaxstarLOG, group = mngt_type:species, lb = 0,class = "sd") 
    ,prior(normal(0, 0.5), nlpar = hmaxstarLOG, group = soil, lb = 0,class = "sd") 

  ),
  iter = 3000, warmup = 1500, chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = 444,
  silent = FALSE
)

# save model object
save(fit_fe_sub250_cmdLOG_student14, file = "03_work/models/model_fits/fit_fe_sub250_cmdLOG_student14.RData")

fit_success <- TRUE
}, silent = TRUE)

# Record end time
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")


# logging
if (fit_success ) {
  cat("==== Model fit completed successfully ====\n")
  cat(sprintf("Runtime: %.2f minutes\n", as.numeric(duration)))
} else {
  cat("==== Model fit failed or model file not found ====\n")
  cat("Check model specification, data input, or CmdStan backend.\n")
}

