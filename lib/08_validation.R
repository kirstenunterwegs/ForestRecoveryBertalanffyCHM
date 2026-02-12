################################################################################
#
# Validation of model performance with later height observations 
#
################################################################################


# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
#
# This script evaluates the brms model fits by comparing predictions
# from early post-disturbance data to later CHM observations for the same pixels.
#
# Specifically, it:
#   1. Loads teh models and performs basic checks.
#   2. Reconstructs the training dataset and standardizes all predictors; 
#      applies identical transformations to the validation (later-year) data.
#   3. Matches first-year and last-year observations by pixel ID and patch.
#   4. Generates posterior height predictions at t1 and t2,
#      computes posterior Δh from full draws and from aggregated means.
#   5. Builds a validation table containing observed vs. predicted Δh and Δt.
#   6. Computes R², RMSE, and bias 
#      produces diagnostic plots for all data and for cases with Δt ≥ 4.
#
# Output:
#   - Predicted vs. observed growth plots for reporting
# -------------------------------------------------------------------------


# ---libraries

library(data.table)
library(ggplot2)
library(brms)
library(patchwork)
library(ggdist)
library(stringr)

# set directory

setwd("~/NAS/Projects/ForestRecovery/")


# -----------------
#  Load fit 
# -----------------


load("03_work/models/model_fits/fit_fe_sub250_cmdLOG_asymlap27.RData")
fit =  fit_fe_sub250_cmdLOG_asymlap27 

# # student distribution model
# load("03_work/models/model_fits/fit_fe_sub250_cmdLOG_student14.RData")
# fit = fit_fe_sub250_cmdLOG_student14


# ---------------------
# pp checks
# ---------------------

pp_check_full = pp_check(fit_fe_sub250_cmdLOG_asymlap27)
ggsave("03_work/data_processed/results/supporting_information/pp_check_full_asyml27.tif", plot = pp_check_full, dpi = 300, width = 12, height = 6)
pp_check_lim = pp_check(fit_fe_sub250_cmdLOG_asymlap27) + xlim(-10,40)
ggsave("03_work/data_processed/results/supporting_information/pp_check_lim_asyml27.tif", plot = pp_check_lim, dpi = 300, width = 12, height = 6)


#pp_check_full = pp_check(fit_fe_sub250_cmdLOG_student14)
#ggsave("03_work/data_processed/results/supporting_information/pp_check_full_stduent14.tif", plot = pp_check_full, dpi = 300, width = 12, height = 6)
#pp_check_lim = pp_check(fit_fe_sub250_cmdLOG_student14) + xlim(-10,40)
#ggsave("03_work/data_processed/results/supporting_information/pp_check_lim_stduent14.tif", plot = pp_check_lim, dpi = 300, width = 12, height = 6)


#--------------------------------------------------------------
# Load, standardize and compare original data with model data
#--------------------------------------------------------------

recovery_landsat_filtered_sev_firstyr = as.data.table(readRDS("03_work/data_processed/analysis.dt/spatial_sub250.rds"))


# now we compare whether this is the same data set as the one used to fit the model
recovery_landsat_filtered_sev_firstyr_fitted = as.data.table(fit$data)

# check nb entries 
nrow(recovery_landsat_filtered_sev_firstyr_fitted) == nrow(recovery_landsat_filtered_sev_firstyr) 

# check response distribution 
hist(recovery_landsat_filtered_sev_firstyr_fitted$h_mean)
hist(recovery_landsat_filtered_sev_firstyr$h_mean)

# check t
hist(recovery_landsat_filtered_sev_firstyr_fitted$t)
hist(recovery_landsat_filtered_sev_firstyr$t)

#--- standardize predictors

# we need to save the means and sds for normalizing the lastyr validation data

slope_mean = mean(recovery_landsat_filtered_sev_firstyr$slope)
slope_sd = sd(recovery_landsat_filtered_sev_firstyr$slope)

aspect_mean = mean(recovery_landsat_filtered_sev_firstyr$aspect)
aspect_sd = sd(recovery_landsat_filtered_sev_firstyr$aspect)

temp_mean = mean(recovery_landsat_filtered_sev_firstyr$temp)
temp_sd = sd(recovery_landsat_filtered_sev_firstyr$temp)

prec_mean = mean(recovery_landsat_filtered_sev_firstyr$prec)
prec_sd = sd(recovery_landsat_filtered_sev_firstyr$prec)

h_elevation_mean = mean(recovery_landsat_filtered_sev_firstyr$elevation)
h_elevation_sd = sd(recovery_landsat_filtered_sev_firstyr$elevation)


recovery_landsat_filtered_sev_firstyr[
  , `:=`(
    slope_std = (slope - slope_mean) / slope_sd
    , aspect_std = (aspect - aspect_mean) / aspect_sd
    , prec_std = (prec - prec_mean) / prec_sd
    , temp_std = (temp - temp_mean) / temp_sd
    , elevation_std = (elevation - h_elevation_mean) / h_elevation_sd
  )
]

recovery_landsat_filtered_sev_firstyr[, logt := log10(t)] 
recovery_landsat_filtered_sev_firstyr[, tfact := as.character(ifelse(t < 30, t, 30))]
recovery_landsat_filtered_sev_firstyr[, patch_ha_log := log(patch_ha)]

recovery_landsat_filtered_sev_firstyr[, t_alt := t]
recovery_landsat_filtered_sev_firstyr[, t_alt := fifelse(t == 250, 35, t)]
recovery_landsat_filtered_sev_firstyr[, log10t := log10(t_alt)] 

patch_ha_log_mean = mean(recovery_landsat_filtered_sev_firstyr$patch_ha_log)
patch_ha_log_sd  = sd(recovery_landsat_filtered_sev_firstyr$patch_ha_log)

recovery_landsat_filtered_sev_firstyr[, patch_ha_log_std := (patch_ha_log - patch_ha_log_mean) / patch_ha_log_sd]


# --- rename mngt and species

# Store original values
recovery_landsat_filtered_sev_firstyr[, mngt_id := mngt_type]
recovery_landsat_filtered_sev_firstyr[, species_id := species]

recovery_landsat_filtered_sev_firstyr[, mngt_type := factor(mngt_type,
                                                            levels = 0:6,
                                                            labels = c("Set aside", "Large private", "Small private",
                                                                       "State forest", "Communal forest", "Federal forest", "Other"))]
species_map <- c(`2` = "Birch", `3` = "Beech", `4` = "Douglas fir", `5` = "Oak", 
                 `6` = "Alder", `8` = "Spruce", `9` = "Pine", `10` = "Larch", 
                 `14` = "Fir", `16` = "ODH", `17` = "ODL", `99` = "Other")
recovery_landsat_filtered_sev_firstyr[, species := species_map[as.character(species)]]
recovery_landsat_filtered_sev_firstyr[, species := factor(species)]

recovery_landsat_filtered_sev_firstyr[, soil_id := soil]
recovery_landsat_filtered_sev_firstyr[, soil := as.numeric(str_extract(soil_id, "\\d+"))]



# --- check correspondence of predictor variables ---

hist(recovery_landsat_filtered_sev_firstyr_fitted$slope_std)
hist(recovery_landsat_filtered_sev_firstyr$slope_std)

hist(recovery_landsat_filtered_sev_firstyr_fitted$aspect_std)
hist(recovery_landsat_filtered_sev_firstyr$aspect_std)

hist(recovery_landsat_filtered_sev_firstyr_fitted$elevation_std)
hist(recovery_landsat_filtered_sev_firstyr$elevation_std)

hist(recovery_landsat_filtered_sev_firstyr_fitted$log10t)
hist(recovery_landsat_filtered_sev_firstyr$log10t)

# check categorical layers 

setorder(recovery_landsat_filtered_sev_firstyr_fitted, mngt_type)
recovery_landsat_filtered_sev_firstyr_fitted[,.N,mngt_type]

setorder(recovery_landsat_filtered_sev_firstyr, mngt_type)
recovery_landsat_filtered_sev_firstyr[,.N,mngt_type]

setorder(recovery_landsat_filtered_sev_firstyr_fitted, species)
recovery_landsat_filtered_sev_firstyr_fitted[,.N, species]

setorder(recovery_landsat_filtered_sev_firstyr, species)
recovery_landsat_filtered_sev_firstyr[,.N, species]

#-------------------------------------------------------------------
# Load, standardize and compare original data with validation data
#------------------------------------------------------------------

# we repeat the exact same thing, except for scaling / standardization
# and, obviously, we don't expect distributions to be exactly the same

load("03_work/data_processed/analysis.dt/recovery_all_landsat_filtered_sev0.75__hmean_hmax_lastyr.Rdata") 
recovery_landsat_filtered_sev_lastyr = as.data.table(recovery_landsat_filtered_sev_lastyr)

# remove species group 99 (other, where we do not have any species classification by Blickensdörfer)
recovery_landsat_filtered_sev_lastyr = recovery_landsat_filtered_sev_lastyr[species != 99] 

# remove all observation from the disturbance year and the year after, to avoid error in disturbance year classification
recovery_landsat_filtered_sev_lastyr[, t := chm_year - disturbance_year]
recovery_landsat_filtered_sev_lastyr[t < 0, t := 250]
recovery_landsat_filtered_sev_lastyr[, patch_ha_log := log(patch_ha)]


# --- rename mngt and species

# Store original values
recovery_landsat_filtered_sev_lastyr[, mngt_id := mngt_type]
recovery_landsat_filtered_sev_lastyr[, species_id := species]

recovery_landsat_filtered_sev_lastyr[, mngt_type := factor(mngt_type,
                                                           levels = 0:6,
                                                           labels = c("Set aside", "Large private", "Small private",
                                                                      "State forest", "Communal forest", "Federal forest", "Other"))]
species_map <- c(`2` = "Birch", `3` = "Beech", `4` = "Douglas fir", `5` = "Oak", 
                 `6` = "Alder", `8` = "Spruce", `9` = "Pine", `10` = "Larch", 
                 `14` = "Fir", `16` = "ODH", `17` = "ODL", `99` = "Other")
recovery_landsat_filtered_sev_lastyr[, species := species_map[as.character(species)]]
recovery_landsat_filtered_sev_lastyr[, species := factor(species)]


### - for models build with data subsets,subset last yr observations for the same pixels: 

recovery_landsat_filtered_sev_lastyr = recovery_landsat_filtered_sev_lastyr[
  recovery_landsat_filtered_sev_firstyr[, .(IDlandsat, patch)],
  on = .(IDlandsat, patch),
  nomatch = 0
]

# now we compare whether this is the same data set as the one used to fit the model
recovery_landsat_filtered_sev_firstyr_fitted = as.data.table(fit$data)

# check nb entries 
nrow(recovery_landsat_filtered_sev_firstyr_fitted) == nrow(recovery_landsat_filtered_sev_lastyr)

# check response distribution
min(recovery_landsat_filtered_sev_firstyr$h_mean)
min(recovery_landsat_filtered_sev_lastyr$h_mean) # height Δ of -31 -> some mining happening at few pixels

hist(recovery_landsat_filtered_sev_firstyr$h_mean, breaks = seq(-50,50))
hist(recovery_landsat_filtered_sev_lastyr$h_mean, breaks = seq(-50,50)) 
# height distribution shift is due to disturbances between both observation points + growth

# expected differences in t
hist(recovery_landsat_filtered_sev_firstyr$t)
hist(recovery_landsat_filtered_sev_lastyr$t)

# --- standardize predictors for newdata

recovery_landsat_filtered_sev_lastyr[
  , `:=`(
    slope_std = (slope - slope_mean) / slope_sd
    , aspect_std = (aspect - aspect_mean) / aspect_sd
    , temp_std = (temp - temp_mean) / temp_sd
    , prec_std = (prec - prec_mean) / prec_sd
    , patch_ha_log_std = (patch_ha_log - patch_ha_log_mean) / patch_ha_log_sd
    , elevation_std = (elevation - h_elevation_mean) / h_elevation_sd
  )
]

recovery_landsat_filtered_sev_lastyr[, logt := log10(t)] 
recovery_landsat_filtered_sev_lastyr[, t_alt := t]
recovery_landsat_filtered_sev_lastyr[, t_alt := fifelse(t == 250, 35, t)]
recovery_landsat_filtered_sev_lastyr[, log10t := log10(t_alt)] 
recovery_landsat_filtered_sev_lastyr[, soil_id := soil]
recovery_landsat_filtered_sev_lastyr[, soil := as.numeric(str_extract(soil_id, "\\d+"))]


# check correspondence of predictor variables

hist(recovery_landsat_filtered_sev_firstyr_fitted$slope_std)
hist(recovery_landsat_filtered_sev_lastyr$slope_std)

hist(recovery_landsat_filtered_sev_firstyr_fitted$aspect_std)
hist(recovery_landsat_filtered_sev_lastyr$aspect_std)

# check categorical layers 

setorder(recovery_landsat_filtered_sev_firstyr, mngt_type)
recovery_landsat_filtered_sev_firstyr_fitted[,.N,mngt_type]

setorder(recovery_landsat_filtered_sev_lastyr, mngt_type)
recovery_landsat_filtered_sev_lastyr[,.N,mngt_type]

setorder(recovery_landsat_filtered_sev_firstyr, species)
recovery_landsat_filtered_sev_firstyr_fitted[,.N, species]

setorder(recovery_landsat_filtered_sev_lastyr, species)
recovery_landsat_filtered_sev_lastyr[,.N, species]

setorder(recovery_landsat_filtered_sev_firstyr, soil)
recovery_landsat_filtered_sev_firstyr[,.N, soil]

setorder(recovery_landsat_filtered_sev_lastyr, soil)
recovery_landsat_filtered_sev_lastyr[,.N, soil]



#%%%%%%%%%%%%%%%%%##
#### Validation ####
#%%%%%%%%%%%%%%%%%##


recovery_landsat_filtered_sev_lastyr = recovery_landsat_filtered_sev_lastyr[as.numeric(t) >= 2]
recovery_landsat_filtered_sev_firstyr = recovery_landsat_filtered_sev_firstyr[
  recovery_landsat_filtered_sev_lastyr[, .(IDlandsat, patch)],
  on = .(IDlandsat, patch),
  nomatch = 0
]


olddata = recovery_landsat_filtered_sev_firstyr[, c("IDlandsat","patch","h_mean","t","slope_std","aspect_std", "elevation_std", "patch_ha_log_std", "soil", "species", "mngt_type","log10t", "logt")]     
newdata = recovery_landsat_filtered_sev_lastyr[, c("IDlandsat","patch","h_mean","t","slope_std","aspect_std","elevation_std", "patch_ha_log_std", "soil", "species", "mngt_type","log10t", "logt")]  
colnames_data = colnames(olddata)
colnames(olddata) = colnames_data
colnames(newdata) = colnames_data


head(olddata)
head(newdata)


#------------ validation calculating the dh_pred before aggregating: 


# Ensure same chunking

olddata[, chunk := floor(1:nrow(olddata)/1000)]
newdata[, chunk := floor(1:nrow(newdata)/1000)]
chunks = unique(olddata$chunk)  # assuming newdata and olddata are aligned

# Allocate space for posterior differences

hpred1 = vector(mode = "list",length = length(chunks))
hpred2 = vector(mode = "list",length = length(chunks))



for(i in seq_along(chunks)){
  chunk_current = chunks[i]
  cat("Chunk", chunk_current, "\n")
  
  olddata_current = olddata[chunk == chunk_current]
  newdata_current = newdata[chunk == chunk_current]
  
  # Ensure rows match!
  stopifnot(nrow(olddata_current) == nrow(newdata_current))
  
  # Get posterior linear predictions for asymmetric Laplace model (linpred can account for asymmetry in model prediction)
  hpred_old = posterior_linpred(fit, newdata = olddata_current, re_formula = NULL) #NA to remove random effects
  hpred_new = posterior_linpred(fit, newdata = newdata_current, re_formula = NULL)
  
  # Get posterior predictions for student model with epred
  # hpred_old = posterior_epred(fit, newdata = olddata_current, re_formula = NULL) 
  # hpred_new = posterior_epred(fit, newdata = newdata_current, re_formula = NULL)
  
  
  hpred_old = colMeans(hpred_old)
  hpred_new = colMeans(hpred_new)
  
  hpred1[[i]] = hpred_old
  hpred2[[i]] = hpred_new
}

hpred1 = unlist(hpred1)
hpred2 = unlist(hpred2)

recovery_landsat_filtered_sev_firstyr$hpred1 = hpred1
recovery_landsat_filtered_sev_lastyr$hpred2 = hpred2


# visualize the distribution of predicted height values (for distributions we split to better see them)

g10 = ggplot(data = recovery_landsat_filtered_sev_firstyr[t < 10], aes(x = t, y = hpred1)) + stat_halfeye(aes(x = as.factor(t))) + theme_classic()
g20 = ggplot(data = recovery_landsat_filtered_sev_firstyr[t >= 10 & t < 20], aes(x = t, y = hpred1)) + stat_halfeye(aes(x = as.factor(t))) + theme_classic()
g30 = ggplot(data = recovery_landsat_filtered_sev_firstyr[t >= 20 & t < 30], aes(x = t, y = hpred1)) + stat_halfeye(aes(x = as.factor(t))) + theme_classic()

g10 + g20 + g30 + plot_layout(ncol = 1)

# create combined data table for validation

validation_berta = merge(
  recovery_landsat_filtered_sev_firstyr[,c("IDlandsat","patch","slope_std","aspect_std","elevation_std", "patch_ha_log_std", "soil","mngt_type","species","disturbance_year","chm_year","t","h_mean","hpred1")] 
  , recovery_landsat_filtered_sev_lastyr[,c("IDlandsat","patch","chm_year","t","h_mean","hpred2")]
  , by = c("IDlandsat","patch")
  , suffixes = c("1", "2"),
  all = FALSE  # only matches are retained, so observations from last obs year, where no disturbances happened in between
)

validation_berta[
  , `:=`(
    dt = t2 - t1
    , dt2 = chm_year2 - chm_year1
    , dh_emp = h_mean2 - h_mean1
    , dh_pred_agg = hpred2 - hpred1 # height growth difference for col mean values per draw
  )
]


# --- save validation dt
save(validation_berta, file = "03_work/analysis/validation/validation_berta_dt_asyml_27_reNULL.RData")
#save(validation_berta, file = "03_work/analysis/validation/validation_berta_dt_student14_reNULL.RData")
# ----

#load("03_work/analysis/validation/validation_berta_dt_asyml_27_reNULL.RData")


# --- raw validation ---

(
  gval = ggplot(validation_berta, aes(x = dh_pred_agg, y = dh_emp)) +
    geom_point(alpha = 0.1, size = 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm") +
    labs(x = "Predicted delta Height growth", y = "Observed delta Height growth", 
         title = "Predicted vs Observed Height growth (no filtering)") +
    theme_classic() + xlim(c(-45,25)) + ylim(c(-45,25))
)


# --- stats ----

# overall

validation_berta[ 
  , list(
    growth = mean(dh_pred_agg)
    , r2 = round(cor(dh_pred_agg, dh_emp)^2,2)
    , rmse = round(sqrt(mean((dh_pred_agg - dh_emp)^2)),3)
    , me = round(mean(dh_pred_agg - dh_emp),3) # bias
  )
]

# by management

validation_berta[
  , list(
    growth = mean(dh_pred_agg)
    , r2 = round(cor(dh_pred_agg, dh_emp)^2,2)
    , rmse = round(sqrt(mean((dh_pred_agg - dh_emp)^2)),3)
    , me = round(mean(dh_pred_agg - dh_emp),3) # bias
  )
  , mngt_type
]


# validation only for expected large positive growth and without disturbances
validation_berta_growthstrict = validation_berta[dt >= 4]

(
  gval_growthstrict = ggplot(validation_berta_growthstrict, aes(x = dh_pred_agg, y = dh_emp)) +
    geom_point(alpha = 0.1, size = 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm") +
    labs(x = "Predicted delta Height growth", y = "Observed delta Height growth", 
         title = "Predicted vs Observed Height growth (only growth >= 4 years)") +
    theme_classic() #+ xlim(c(-45,25)) + ylim(c(-45,25))
)


# visualize predicted growth versus observed growth

ggplot(data = validation_berta_growthstrict, aes(x = hpred1, y = h_mean1)) + geom_point(alpha = 0.1, size = 0.25) + theme_classic() + geom_abline(intercept = 0, slope = 1) + geom_smooth(method = "lm") + ylim(c(0,30)) + xlim(c(0,30))
ggplot(data = validation_berta_growthstrict, aes(x = hpred2, y = h_mean2)) + geom_point(alpha = 0.1, size = 0.25) + theme_classic() + geom_abline(intercept = 0, slope = 1) + geom_smooth(method = "lm") + ylim(c(0,30)) + xlim(c(0,30))
ggplot(data = validation_berta_growthstrict, aes(x = dh_pred_agg, y = dh_emp)) + geom_point(alpha = 0.1, size = 0.25) + theme_classic() + geom_abline(intercept = 0, slope = 1) + geom_smooth(method = "lm") #+ ylim(c(-5,10)) + xlim(c(-5,10))

# ggplot(validation_berta_growthstrict, aes(x = hpred2-hpred1, y = h_mean2-h_mean1)) +
#   geom_hex(bins = 80) +
#   scale_fill_viridis_c(option = "D") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   theme_classic()


# --- stats ---- for strict growth pixels

# overall

validation_berta_growthstrict[
  , list(
    growth = mean(dh_pred_agg)
    , r2 = round(cor(dh_pred_agg, dh_pred_agg)^2,2)
    , rmse = round(sqrt(mean((dh_pred_agg - dh_emp)^2)),3)
    , me = round(mean(dh_pred_agg - dh_emp),3) # bias
  )
]

# by management

validation_berta_growthstrict[
  , list(
    growth = mean(dh_pred_agg)
    , r2 = round(cor(dh_pred_agg, dh_emp)^2,2)
    , rmse = round(sqrt(mean((dh_pred_agg - dh_emp)^2)),3)
    , me = round(mean(dh_pred_agg - dh_emp),3) # bias
  )
  , mngt_type
]


# ---------------------------------------------
# plot validation 
# ---------------------------------------------


# --- all observations

stats_all <- validation_berta[
  , list(
    growth = mean(dh_pred_agg),
    r2 = round(cor(dh_pred_agg, dh_emp)^2, 2),
    rmse = round(sqrt(mean((dh_pred_agg - dh_emp)^2)), 3),
    me = round(mean(dh_pred_agg - dh_emp), 3)
  )
]

# Create the annotation text

annotation_text_all <- paste0(
  "Mean Predicted Growth: ", round(stats_all$growth, 2), "\n",
  "R²: ", stats_all$r2, "\n",
  "RMSE: ", stats_all$rmse, "\n",
  "Bias (ME): ", stats_all$me
)

# Build the plot with annotation

gval <- ggplot(validation_berta, aes(x = dh_pred_agg, y = dh_emp)) +
  geom_point(alpha = 0.1, size = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted delta Height growth",
    y = "Observed delta Height growth",
    title = "Predicted vs Observed Height growth (no filtering)"
  ) +
  annotate(
    "text",
    x = -45,  
    y = 25,   
    hjust = 0, vjust = 1,
    label = annotation_text_all,
    size = 3.5
  ) +
  theme_classic() +
  xlim(c(-45, 25)) +
  ylim(c(-45, 25))


# --- observations where we expect growth (strict growth)


# Compute statistics

stats_growthstrict <- validation_berta_growthstrict[
  , list(
    growth = mean(dh_pred_agg),
    r2 = round(cor(dh_pred_agg, dh_emp)^2, 2),
    rmse = round(sqrt(mean((dh_pred_agg - dh_emp)^2)), 3),
    me = round(mean(dh_pred_agg - dh_emp), 3)
  )
]

# Create the annotation text

annotation_text <- paste0(
  "Mean Predicted Growth: ", round(stats_growthstrict$growth, 2), "\n",
  "R²: ", stats_growthstrict$r2, "\n",
  "RMSE: ", stats_growthstrict$rmse, "\n",
  "Bias (ME): ", stats_growthstrict$me
)

# Build the plot with annotation

gval_growthstrict <- ggplot(validation_berta_growthstrict, aes(x = dh_pred_agg, y = dh_emp)) +
  geom_point(alpha = 0.1, size = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted delta Height growth",
    y = "Observed delta Height growth",
    title = "Predicted vs Observed Height growth (only growth >= 4 years)"
  ) +
  annotate(
    "text",
    x = min(validation_berta_growthstrict$dh_pred_agg, na.rm = TRUE), 
    y = max(validation_berta_growthstrict$dh_emp, na.rm = TRUE),
    hjust = 0, vjust = 1,
    label = annotation_text,
    size = 3.5
  ) +
  theme_classic()



gcombo = gval + gval_growthstrict

ggsave(plot = gcombo, filename = "03_work/analysis/validation/gcombo_asyml_27_agg.png", width = 14, height = 7)
#ggsave(plot = gcombo, filename = "03_work/analysis/validation/gcombo_student14_agg.png", width = 14, height = 7) 


