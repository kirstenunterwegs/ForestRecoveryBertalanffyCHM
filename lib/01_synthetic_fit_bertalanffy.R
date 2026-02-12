############################################################################
#
# Check model's ability to disentangle g (heigt growth rate) and 
# hOffset (height offset = disturbance legacies) using synthetic data
#
############################################################################

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
# 
# We use the Bertalanffy growth model to estimate three key parameters:
#   - g:       growth rate
#   - hOffset: initial height offset (representing legacy effects)
#   - hmax:    maximum attainable height (asymptote)
#
# The basic model formulation is:
#   h = (hmax - hOffset) * (1 - exp(-g * t))^3
#
# Since disturbances leave legacies in the forest, growth does not start 
# from zero. To capture this, we include an additional height offset term 
# (hOffset), representing the initial height at the start of observation:
#   h = (hmax - hOffset) * (1 - exp(-g * t))^3 + hOffset
#
# There ia an interplay between hOffset (starting height) 
# and g (growth rate) in determining height (h) to reach a certain height (h).
#
# We generate synthetic height data (h) using different combinations of g, 
# hOffset, and hmax to test whether the model can reliably distinguish between 
# these parameters, even when random noise is added to simulate measurement 
# uncertainty (as typically found in CHM-derived height data).
#
# The simulated data cover growth from 0–35 years and include mature, 
# undisturbed stands approximated at 150 years. (and tested other time windows)
#
# -------------------------------------------------------------------------
# Notes
# -------------------------------------------------------------------------
# In many simulations, we use the term (hmax - hOffset) directly.
# In some cases, we simplify this by defining a new parameter:
#   hmax_star = hmax - hOffset
#
# The reason is that the surrounding canopy height (used as a proxy for hmax) 
# can be biased. Therefore, we later estimate hmax directly within the model. 
# Whenever this modification is applied, it is explicitly labeled as such.
#
# -------------------------------------------------------------------------
# Implementation note
# -------------------------------------------------------------------------
# Before each model fitting routine (or "suite"), new synthetic data are 
# generated. This allows immediate experimentation with different parameter 
# variations of interest.
# -------------------------------------------------------------------------


# --- libraries 

# data analysis & wrangling
library(data.table)
library(brms)
library(tictoc)
library(tidyr)
library(dplyr)
library(stringr)

# plotting
library(ggplot2)
library(patchwork)
library(viridis)


# --- set working directory 

setwd("~/data") # or how you name your directory



################# Single simulation ########################


n_synthetic <- 800  # number of observations
n_fixed <- 20       # Number of rows with t = 150
n_random <- n_synthetic - n_fixed # number of observation in other time window


# prepare different hmax for different disturbance patch IDs to generate synthetic data; patches represent disturbance patches

selected_patches = c("100603", "100174", "100169", "1007344", "1002348", "1010147", "1010292", "1009814", "1007813", "100030", "100894", "1010434", "1009612",
                     "1010869", "1008131", "1008392", "1009341", "1010831", "1007956", "1009902", "1007126", "1005823", "100167", "1007497", "1009106", "1002333",
                     "1007860", "1009058", "1002598", "1007231", "101011", "1008853", "1004322", "1010448", "1010735", "100399", "1010605", "1007710", "1005271",
                     "1009749")

patch_hmax_mapping <- tibble(
  patch = selected_patches,
  hmax = c(31.52689, 29.60121, 31.10225, 16.64266, 30.91800, 32.62781, 30.14731, 29.48504, 30.38545, 31.12933, 23.68818, 32.59286, 26.38217, 34.31919, 26.10335,
           30.07681, 29.34345, 30.10626, 40.38669, 29.38597, 30.38545, 30.05920, 36.28186, 27.50942, 31.43734, 31.53233, 28.01574, 24.54890, 30.24358, 30.34470,
           28.09915, 31.21759, 26.95848, 29.72182, 27.35790, 29.56066, 30.91459, 33.10059, 27.50942, 30.07681)
)


# set time ranges where we have observation data:

synthetic_data <- tibble(
  t = runif(n_synthetic, min = 0, max = 150),  # Generate a sequence of time points
  # t = c(runif(n_random, min = 0, max = 35), rep(150, n_fixed)),  # Mix random values with fixed 150s
  # t = c(runif(n_random, min = 0, max = 15), rep(150, n_fixed)),  # shortly after disturbance
  # t = c(runif(n_random, min = 20, max = 35), rep(150, n_fixed)),  # later after disturbance
  # t = c(runif(n_synthetic, min = 20, max = 35)),  # later after disturbance
  patch = sample(selected_patches, n_synthetic, replace = TRUE)  # Randomly assign patches to rows
) %>%
  left_join(patch_hmax_mapping, by = "patch") 


# -- without different management classes

g = 0.08 # g between 0.03 and 0.09 is realistic
# we add noise to the height estimates ad varying degrees; we subtract the 2 from hmax, as we set this as the hOffset
synthetic_data$h = (synthetic_data$hmax - 2) * (1 - exp(-g * synthetic_data$t))^3 + rnorm(nrow(synthetic_data), 2, 6) # mean error 2m with SD 6 m; alternative: rnorm(nrow(synthetic_data), 0, 3)


# -- with different management classes, mimicking different effects on growth rates and hOffset per management type

patch_mngt_mapping <- tibble(
  patch = selected_patches,
  mngt = sample(0:1, length(selected_patches), replace = TRUE)  # Assign unique mngt to each patch
)
synthetic_data <- synthetic_data %>%
  left_join(patch_mngt_mapping, by = "patch")


synthetic_data$g <- ifelse(synthetic_data$mngt == 1, 0.03, 0.05)# Set g based on mngt value
synthetic_data$h <- (synthetic_data$hmax - 2 ) * (1 - exp(-synthetic_data$g * synthetic_data$t))^3 + rnorm(nrow(synthetic_data), 2, 6) 
synthetic_data$g <- NULL


# -- fit the model 

tic("start fit")

synthetic_fit =  brms::brm(
  bf(
    h ~ (hmax-hOffset) * (1 - exp(-g * t))^3 + hOffset
    , g ~ 1  + (1 | patch) + mngt
    , hOffset ~ 1 + (1 | patch) + mngt
    #, sigma ~ 1 + (1 | patch) # potentially let sigma vary by patch
    # , sigma ~ 1 
    , nl = TRUE
  )
  , data = synthetic_data
  , family = "student" 
  , prior = c(
    prior(normal(0.05, 0.02), lb = 0, nlpar = g, class = "b"),   
    prior(normal(5,5), lb = 0, nlpar = hOffset, class = "b"),
    prior(normal(0.025, 0.025), lb = 0, nlpar = g, class = "sd"),   
    prior(normal(10,10), lb = 0, nlpar = hOffset, class = "sd")
  )
  , iter = 2500, warmup = 1000, chains = 4, cores = cores
  , control = list(adapt_delta = 0.99, max_treedepth = 10)
  , seed = 444
  , silent = F
)

toc()

save(synthetic_fit, file = "03_work/models/synthetic_fit_gaussian_sigma_v1.RData")
load("03_work/models/synthetic_fit_gaussian_sigma_v1.RData")


# check the model:

print(synthetic_fit, digits = 3)
plot(synthetic_fit)

# extract model estimates

g = fixef(synthetic_fit)["g_Intercept", "Estimate"]
ranef_g = ranef(synthetic_fit)$patch[, , "g_Intercept"]
g_patch = g + ranef_g[, "Estimate"]
hist(g_patch)
mean(g_patch)


# --- check g estimate with management differentiation

g_patch_df <- tibble(g_patch = g_patch, patch = rownames(ranef_g))
g_patch_df <- left_join(g_patch_df, patch_mngt_mapping, by = "patch")

g_patch_df %>% group_by(mngt) %>% summarise(g_mean = mean(g_patch))

ggplot(g_patch_df, aes(x = g_patch, fill = factor(mngt))) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  labs(title = "Distribution of g values per Management Type",
       x = "g value", y = "Frequency", fill = "Management Type (mngt)") +
  scale_fill_manual(values = c("blue", "red"), labels = c("mngt = 0", "mngt = 1")) +
  theme_minimal()


# --- check hOffset estimate with management differentiation

hOffset = fixef(synthetic_fit)["hOffset_Intercept", "Estimate"]
ranef_hOffset = ranef(synthetic_fit)$patch[, , "hOffset_Intercept"]
hOffset_patch = hOffset + ranef_hOffset[, "Estimate"]
hist(hOffset_patch)
mean(hOffset_patch)


# %%%%%%% Multiple repetitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

################################################################################
#
# Synthetic model fits (without estimating hmax)
# We add different noise levels to the h values (measured later by the CHM)
# and went through different combinations of time windows of measurements, 
# growth rates (g) and noise levels 
#
################################################################################


# Define variations

t_ranges <- list(
  #"0-20" = function(n) runif(n, min = 0, max = 20),
  #"15-30" = function(n) runif(n, min = 15, max = 30),
  "0-35" = function(n) runif(n, min = 15, max = 30),
  "0-150" = function(n) runif(n, min = 0, max = 150),
  "0-35,150" = function(n, fixed) c(runif(n - fixed, min = 0, max = 35), rep(150, fixed))
)

g_values <- c(0.02, 0.05, 0.09)  # Growth rates
h_noise_levels <- c(1, 5, 10)  # Noise standard deviations

# Expand grid of all parameter combinations
param_grid <- expand.grid(
  t_range = names(t_ranges),
  g = g_values,
  h_noise = h_noise_levels
)

# Set constants
n_synthetic <- 600
n_fixed <- 20
n_random <- n_synthetic - n_fixed
num_replicates <- 3  # Number of repetitions of model fits

# prepare different hmax for different patch IDs to generate synthetic data; patches represent disturbance patches

selected_patches = c("100603", "100174", "100169", "1007344", "1002348", "1010147", "1010292", "1009814", "1007813", "100030", "100894", "1010434", "1009612",
                     "1010869", "1008131", "1008392", "1009341", "1010831", "1007956", "1009902", "1007126", "1005823", "100167", "1007497", "1009106", "1002333",
                     "1007860", "1009058", "1002598", "1007231", "101011", "1008853", "1004322", "1010448", "1010735", "100399", "1010605", "1007710", "1005271",
                     "1009749")

patch_hmax_mapping <- tibble(
  patch = selected_patches,
  hmax = c(31.52689, 29.60121, 31.10225, 16.64266, 30.91800, 32.62781, 30.14731, 29.48504, 30.38545, 31.12933, 23.68818, 32.59286, 26.38217, 34.31919, 26.10335,
           30.07681, 29.34345, 30.10626, 40.38669, 29.38597, 30.38545, 30.05920, 36.28186, 27.50942, 31.43734, 31.53233, 28.01574, 24.54890, 30.24358, 30.34470,
           28.09915, 31.21759, 26.95848, 29.72182, 27.35790, 29.56066, 30.91459, 33.10059, 27.50942, 30.07681)
)


# Create an empty data frame to store results
results_df <- tibble()


# Start the repeated simulations
for (replicate_id in 1:num_replicates) {
  for (i in seq_len(nrow(param_grid))) {
    
    # Extract parameters
    t_range_name <- param_grid$t_range[i]
    g_val <- param_grid$g[i]
    h_noise <- param_grid$h_noise[i]
    
    # Generate synthetic data based on t range
    if (t_range_name == "0-35,150") {
      synthetic_t <- t_ranges[[t_range_name]](n_synthetic, n_fixed)
    } else {
      synthetic_t <- t_ranges[[t_range_name]](n_synthetic)
    }
    
    synthetic_data <- tibble(
      t = synthetic_t,  
      patch = sample(selected_patches, n_synthetic, replace = TRUE)
    ) %>%
      left_join(patch_hmax_mapping, by = "patch") %>%
      mutate(h = (hmax - 10) * (1 - exp(-g_val * t))^3 + rnorm(n_synthetic, 10, h_noise))  # Apply noise first 1m, second time with 10 m Offset + noise
    
    # Fit Bayesian model
    tic(paste(
      "Replicate:", replicate_id, 
      "- t_range:", t_range_name,
      "- g =", g_val, 
      "- h_noise =", h_noise
    ))
    
    synthetic_fit <- brm(
      bf(
        h ~ (hmax - hOffset) * (1 - exp(-g * t))^3 + hOffset,
        g ~ 1 + (1 | patch),
        hOffset ~ 1 + (1 | patch),
        nl = TRUE
      ),
      data = synthetic_data,
      family = "student",
      prior = c(
        prior(normal(0.05, 0.02), lb = 0, nlpar = g, class = "b"),   
        prior(normal(5,5), lb = 0, nlpar = hOffset, class = "b"),
        prior(normal(0.025, 0.025), lb = 0, nlpar = g, class = "sd"),   
        prior(normal(10,10), lb = 0, nlpar = hOffset, class = "sd")
      ),
      iter = 2500, warmup = 1000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.95, max_treedepth = 10),
      seed = 444,
      silent = FALSE
    )
    toc()
    
    # Extract estimates
    g_estimate <- fixef(synthetic_fit)["g_Intercept", "Estimate"]
    ranef_g <- ranef(synthetic_fit)$patch[, , "g_Intercept"]
    g_patch <- g_estimate + ranef_g[, "Estimate"]
    
    hOffset_estimate <- fixef(synthetic_fit)["hOffset_Intercept", "Estimate"]
    ranef_hOffset <- ranef(synthetic_fit)$patch[, , "hOffset_Intercept"]
    hOffset_patch <- hOffset_estimate + ranef_hOffset[, "Estimate"]
    
    # Convert to tibble for storing
    result_entry <- tibble(
      replicate = replicate_id,
      t_range = t_range_name,
      g_input = g_val,
      h_noise_input = h_noise,
      patch = names(g_patch),
      g_patch_estimate = g_patch,
      hOffset_patch_estimate = hOffset_patch
    )
    
    # Append to results_df
    results_df <- bind_rows(results_df, result_entry)
  }
}

df.hOffset_10 <- results_df



# Convert results_df to long format for plotting

plot_data <- results_df %>%
  pivot_longer(cols = c(g_patch_estimate, hOffset_patch_estimate),
               names_to = "parameter",
               values_to = "estimate") %>%
  mutate(
    parameter = factor(parameter, levels = c("g_patch_estimate", "hOffset_patch_estimate"),
                       labels = c("g_patch", "hOffset_patch")),
    g = factor(g_input),
    h_noise = factor(h_noise_input),
    #t_range = factor(t_range, levels = c("0-20", "15-30" ,"0-35", "0-150", "0-35,150"))
    t_range = factor(t_range, levels = c("0-35", "0-150", "0-35,150"))
  ) %>%
  
  mutate(g_label = case_when(
    g == 0.02 ~ "g 0.02",
    g == 0.05 ~ "g 0.05",
    g == 0.09 ~ "g 0.09",
    TRUE ~ as.character(g)
  )) %>%
  
  mutate(h_noise_label = factor(case_when(
    h_noise == 1 ~ "h noise 1",
    h_noise == 5 ~ "h noise 5",
    h_noise == 10 ~ "h noise 10",
    TRUE ~ as.character(h_noise)), 
   levels = c("h noise 1", "h noise 5", "h noise 10")
  ))


# -- plot g

g_patch_plot <- ggplot(filter(plot_data, parameter == "g_patch"), aes(x = estimate, fill = t_range)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  geom_vline(aes(xintercept = as.numeric(g_input)), linetype = "dashed", size = 1) +
  facet_wrap(~g_label + h_noise_label) + #, scales = "free"
  theme_minimal() +
  labs(title = "Distribution of g_patch across different model configurations; hOffset = 10 m",
       x = "g_patch Estimate", y = "Count", fill = "t_range") +
  theme(legend.position = "bottom")


ggsave("03_work/analysis/synthetic_simulation/g_patch_hist_hOffset_10m.tiff", plot = g_patch_plot, width = 10, height = 6, dpi = 300)

# -- plot hOffset

hOffset_patch_plot <- ggplot(filter(plot_data, parameter == "hOffset_patch"), aes(x = estimate, fill = t_range)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  #geom_vline(xintercept = 1, linetype = "dashed", size = 1) +  # Reference line at 1
  geom_vline(xintercept = 10, linetype = "dashed", size = 1) +  # Reference line at 10
  facet_wrap(~g_label + h_noise_label) + # , scales = "free"
  theme_minimal() +
  labs(title = "Distribution of hOffset_patch across different model configurations; hOffset = 10m", # change hOffset label depending on teh level you chose
       x = "hOffset_patch Estimate", y = "Count", fill = "t_range") +
  theme(legend.position = "bottom")

ggsave("03_work/analysis/synthetic_simulation/hOffset_patch_hist_hOffset_10m.tiff", plot = hOffset_patch_plot, width = 10, height = 6, dpi = 300)


hOffset_patch_plot_sub <- ggplot(filter(plot_data, parameter == "hOffset_patch" & t_range == "0-35,150"), aes(x = estimate, fill = t_range)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  #geom_vline(xintercept = 1, linetype = "dashed", size = 1) +  # Reference line at 1
  geom_vline(xintercept = 10, linetype = "dashed", size = 1) +  # Reference line at 10
  facet_wrap(~g_label + h_noise_label) + # , scales = "free"
  theme_minimal() +
  labs(title = "Distribution of hOffset_patch across different model configurations; hOffset = 10m",
       x = "hOffset_patch Estimate", y = "Count", fill = "t_range") +
  theme(legend.position = "bottom")

ggsave("03_work/analysis//synthetic_simulation/hOffset_patch_hist_hOffset_10m_t0-15.150.tiff", plot = hOffset_patch_plot_sub, width = 10, height = 6, dpi = 300)



################################################################################
#
#  Add noise to hmax
# Note: we ended up not using hmax extracted from the surrounding environment
#       around disturbance patch, as it is a very bad predictor for the real 
#       hmax; we still keep this here as we tested the implications
#
################################################################################


# Define the parameter grid

g_values <- c(0.03, 0.06, 0.09)  # Growth rates
h_noise_levels <- c(1, 5,10)  # Noise standard deviations
hmax_noise_levels <- c(1, 5, 10) 
param_grid <- expand.grid(g = g_values, h_noise = h_noise_levels, hmax_noise = hmax_noise_levels)

n_synthetic <- 800
n_fixed <- 20
n_random <- n_synthetic - n_fixed
num_replicates <- 3  # Number of repetitions of model fits

# prepare different hmax for different patch IDs to generate synthetic data; patches represent distrubance patches

selected_patches = c("100603", "100174", "100169", "1007344", "1002348", "1010147", "1010292", "1009814", "1007813", "100030", "100894", "1010434", "1009612",
                     "1010869", "1008131", "1008392", "1009341", "1010831", "1007956", "1009902", "1007126", "1005823", "100167", "1007497", "1009106", "1002333",
                     "1007860", "1009058", "1002598", "1007231", "101011", "1008853", "1004322", "1010448", "1010735", "100399", "1010605", "1007710", "1005271",
                     "1009749")

patch_hmax_mapping <- tibble(
  patch = selected_patches,
  hmax = c(31.52689, 29.60121, 31.10225, 16.64266, 30.91800, 32.62781, 30.14731, 29.48504, 30.38545, 31.12933, 23.68818, 32.59286, 26.38217, 34.31919, 26.10335,
           30.07681, 29.34345, 30.10626, 40.38669, 29.38597, 30.38545, 30.05920, 36.28186, 27.50942, 31.43734, 31.53233, 28.01574, 24.54890, 30.24358, 30.34470,
           28.09915, 31.21759, 26.95848, 29.72182, 27.35790, 29.56066, 30.91459, 33.10059, 27.50942, 30.07681)
)

# Create an empty dataframe to store results
results_df <- tibble()

# Start the repeated simulations

for (replicate_id in 1:num_replicates) {
  for (i in seq_len(nrow(param_grid))) {
    
    g_val <- param_grid$g[i]
    h_noise <- param_grid$h_noise[i]
    hmax_noise <- param_grid$hmax_noise[i]
    
    # Generate synthetic data
    
    synthetic_data <- tibble(
      t = c(runif(n_random, min = 0, max = 35), rep(150, n_fixed)),  # Time points
      patch = sample(selected_patches, n_synthetic, replace = TRUE)
    ) %>%
      left_join(patch_hmax_mapping, by = "patch") %>%
      mutate(h = (hmax - 1) * (1 - exp(-g_val * t))^3 + rnorm(n_synthetic, 1, h_noise)) %>% # apply noise/ create hOffset
      mutate(hmax = hmax + rnorm(n_synthetic, 0, hmax_noise)) # apply noise on hmax
    
    # Fit Bayesian model
    
    tic(paste("Replicate:", replicate_id, "- Fitting model for g =", g_val, "and h_noise =", h_noise))
    synthetic_fit <- brm(
      bf(
        h ~ (hmax - hOffset) * (1 - exp(-g * t))^3 + hOffset,
        g ~ 1 + (1 | patch),
        hOffset ~ 1 + (1 | patch),
        nl = TRUE
      ),
      data = synthetic_data,
      family = "student",  # Changed from "student" to "gaussian"
      prior = c(
        prior(normal(0.05, 0.02), lb = 0, nlpar = g, class = "b"),
        prior(normal(5,5), lb = 0, nlpar = hOffset, class = "b"),
        prior(normal(0.025, 0.025), lb = 0, nlpar = g, class = "sd"),
        prior(normal(10,10), lb = 0, nlpar = hOffset, class = "sd")
      ),
      iter = 2500, warmup = 1000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.95, max_treedepth = 10),
      seed = 444,
      silent = FALSE
    )
    toc()
    
    # Extract estimates
    
    g_estimate <- fixef(synthetic_fit)["g_Intercept", "Estimate"]
    ranef_g <- ranef(synthetic_fit)$patch[, , "g_Intercept"]
    g_patch <- g_estimate + ranef_g[, "Estimate"]
    
    hOffset_estimate <- fixef(synthetic_fit)["hOffset_Intercept", "Estimate"]
    ranef_hOffset <- ranef(synthetic_fit)$patch[, , "hOffset_Intercept"]
    hOffset_patch <- hOffset_estimate + ranef_hOffset[, "Estimate"]
    
    # Convert to tibble for storing
    
    result_entry <- tibble(
      replicate = replicate_id,
      g_input = g_val,
      h_noise_input = h_noise,
      hmax_noise_input = hmax_noise,
      patch = names(g_patch),
      g_patch_estimate = g_patch,
      hOffset_patch_estimate = hOffset_patch
    )
    
    # Append to results_df
    
    results_df <- bind_rows(results_df, result_entry)
  }
}

df_false_hmax <- results_df

# Save results as RDS

saveRDS(df_false_hmax, "03_work/analysis/synthetic_simulation/models/sensitivity_synthetic_sim_reps_falsehmax.rds")
#results_df <- readRDS("03_work/analysis/synthetic_simulation/models/sensitivity_synthetic_sim_reps_falsehmax.rds")


# Convert results_df to long format for plotting

plot_data <- results_df %>%
  pivot_longer(cols = c(g_patch_estimate, hOffset_patch_estimate),
               names_to = "parameter",
               values_to = "estimate") %>%
  mutate(
    parameter = factor(parameter, levels = c("g_patch_estimate", "hOffset_patch_estimate"),
                       labels = c("g_patch", "hOffset_patch")),
    g = factor(g_input),
    h_noise = factor(h_noise_input),
    hmax_noise_input = factor(hmax_noise_input)
  ) %>%
  
  mutate(g_label = case_when(
    g == 0.03 ~ "g 0.03",
    g == 0.06 ~ "g 0.06",
    g == 0.09 ~ "g 0.09",
    TRUE ~ as.character(g)
  )) %>%
  
  mutate(h_noise_label = factor(case_when(
    h_noise == 1 ~ "h noise 1",
    h_noise == 5 ~ "h noise 5",
    h_noise == 10 ~ "h noise 10",
    TRUE ~ as.character(h_noise)), 
    levels = c("h noise 1", "h noise 5", "h noise 10")
  ))


# g

g_patch_hmax <- ggplot(filter(plot_data, parameter == "g_patch"), aes(x = estimate, fill = hmax_noise_input)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  geom_vline(aes(xintercept = as.numeric(g_input)), linetype = "dashed", size = 1) +
  facet_wrap(~g_label + h_noise_label, scales = "free") + #, scales = "free"
  theme_minimal() +
  labs(title = "Distribution of g_patch across different hmax noise level, hOffset = 1, t= 0-35; 150, scales=free",
       x = "g_patch Estimate", y = "Count") +
  theme(legend.position = "bottom")

ggsave("03_work/analysis/synthetic_simulation/synthetic_simulation/g_estimation_hmax_noise_sim_freescales.tiff", plot = g_patch_hmax, width = 10, height = 6, dpi = 300)


# hOffset

hOffset_patch_hmax <-ggplot(filter(plot_data, parameter == "hOffset_patch"), aes(x = estimate, fill = hmax_noise_input)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +  # Reference line at 1
  #geom_vline(xintercept = 10, linetype = "dashed", size = 1) +  # Reference line at 10
  facet_wrap(~g_label + h_noise_label, scales = "free") + # , scales = "free"
  theme_minimal() +
  labs(title = "Distribution of hOffset_patch across different different hmax noise level, hOffset = 1, t= 0-35; 150, scales = free",
       x = "hOffset_patch Estimate", y = "Count") +
  theme(legend.position = "bottom")


ggsave("03_work/analysis/synthetic_simulation/synthetic_simulation/hOffset_estimation_hmax_noise_sim_freescales.tiff", plot = hOffset_patch_hmax, width = 10, height = 6, dpi = 300)




################################################################################
#
#  Estimate hmaxstar with hmax
# Modification we pointed to in the Setup description:
# We estimate hmaxstar ( = hmax - hOffset) with the (measured) hmax, as we do 
# not trust the surrounding hmax measured to be a reliable asymptote for the model.
# Hence we generate a seperate parameter, which is estimated by the model.
#
################################################################################


# Define variations

t_ranges <- list(
  #"0-20" = function(n) runif(n, min = 0, max = 20),
  #"15-30" = function(n) runif(n, min = 15, max = 30),
  #"0-35" = function(n) runif(n, min = 15, max = 30),
  #"0-150" = function(n) runif(n, min = 0, max = 150),
  "0-35,150" = function(n, fixed) c(runif(n - fixed, min = 0, max = 35), rep(150, fixed))
)

g_values <- c(0.05)  # Growth rates
h_noise_levels <- c(1, 3, 5)  # Noise standard deviations

# --- alternative simulation setup
# g_values <- c(0.02, 0.05, 0.09)  # Growth rates
# h_noise_levels <- c(1, 5, 10)  # Noise standard deviations

# Expand grid of all parameter combinations

param_grid <- expand.grid(
  t_range = names(t_ranges),
  g = g_values,
  h_noise = h_noise_levels
)

# Set constants

n_synthetic <- 600
n_fixed <- 20
n_random <- n_synthetic - n_fixed
num_replicates <- 3  # Number of repetitions


# prepare different hmax for different patch IDs

selected_patches = c("100603", "100174", "100169", "1007344", "1002348", "1010147", "1010292", "1009814", "1007813", "100030", "100894", "1010434", "1009612",
                     "1010869", "1008131", "1008392", "1009341", "1010831", "1007956", "1009902", "1007126", "1005823", "100167", "1007497", "1009106", "1002333",
                     "1007860", "1009058", "1002598", "1007231", "101011", "1008853", "1004322", "1010448", "1010735", "100399", "1010605", "1007710", "1005271",
                     "1009749")

patch_hmax_mapping <- tibble(
  patch = selected_patches,
  hmax = c(31.52689, 29.60121, 31.10225, 16.64266, 30.91800, 32.62781, 30.14731, 29.48504, 30.38545, 31.12933, 23.68818, 32.59286, 26.38217, 34.31919, 26.10335,
           30.07681, 29.34345, 30.10626, 40.38669, 29.38597, 30.38545, 30.05920, 36.28186, 27.50942, 31.43734, 31.53233, 28.01574, 24.54890, 30.24358, 30.34470,
           28.09915, 31.21759, 26.95848, 29.72182, 27.35790, 29.56066, 30.91459, 33.10059, 27.50942, 30.07681)
)


# Create an empty data frame to store results
results_df <- tibble()


# Start the repeated simulations

for (replicate_id in 1:num_replicates) {  # number of repeated simulations
  for (i in seq_len(nrow(param_grid))) {  # go through each parameter & noise combination
    
    # Extract parameters
    t_range_name <- param_grid$t_range[i]
    g_val <- param_grid$g[i]
    h_noise <- param_grid$h_noise[i]
    
    # Generate synthetic data based on t range
    synthetic_t <- as.numeric(NA)
    
    if (t_range_name == "0-35,150") {
      synthetic_t <- t_ranges[[t_range_name]](n_synthetic, n_fixed) # create the observation range we have
    } else {
      synthetic_t <- t_ranges[[t_range_name]](n_synthetic) # continuous observation period
    }
    
    synthetic_data <- tibble(
      t = synthetic_t,  
      patch = sample(selected_patches, n_synthetic, replace = TRUE)
    ) %>%
      left_join(patch_hmax_mapping, by = "patch") %>%
      mutate(h = (hmax-10) * (1 - exp(-g_val * t))^3 + rnorm(n_synthetic, 10, h_noise))  
    
    
    # Fit Bayesian model
    
    tic(paste(
      "Replicate:", replicate_id, 
      "- t_range:", t_range_name,
      "- g =", g_val, 
      "- h_noise =", h_noise
    ))
    
    synthetic_fit3 <- brm(  
      bf(
        h ~ hmaxstar * (1 - exp(-g * t))^3 + hOffset,
        g ~ 1 + (1 | patch),
        hOffset ~ 1 + (1 | patch),
        hmaxstar ~ 1 + hmax, # using patch surrounding hmax as predictor for hmaxstar
        nl = TRUE
      ),
      data = synthetic_data,
      family = "student",
      prior = c(
        prior(normal(0.05, 0.02), lb = 0, nlpar = g, class = "b"),   
        prior(normal(5,5), lb = 0, nlpar = hOffset, class = "b"),
        prior(normal(-5,5), nlpar = hmaxstar, class = "b", coef = "Intercept"), # intercept
        prior(normal(1,1), nlpar = hmaxstar, class = "b", coef = "hmax"), # slope
        prior(normal(0.025, 0.025), lb = 0, nlpar = g, class = "sd"),   
        prior(normal(10,10), lb = 0, nlpar = hOffset, class = "sd")
        # prior(normal(5,5), lb = 0, nlpar = hmaxstar, class = "sd")
      ),
      iter = 2500, warmup = 1000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 10),
      seed = 444,
      silent = FALSE
    )
    
    toc()
    
    get_prior(synthetic_fit3)
    
    # Extract estimates
    
    g_estimate <- fixef(synthetic_fit3)["g_Intercept", "Estimate"]
    ranef_g <- ranef(synthetic_fit3)$patch[, , "g_Intercept"]
    g_patch <- g_estimate + ranef_g[, "Estimate"]
    
    hOffset_estimate <- fixef(synthetic_fit3)["hOffset_Intercept", "Estimate"]
    ranef_hOffset <- ranef(synthetic_fit3)$patch[, , "hOffset_Intercept"]
    hOffset_patch <- hOffset_estimate + ranef_hOffset[, "Estimate"]
    
    hmaxstar_intercept <- fixef(synthetic_fit3)["hmaxstar_Intercept", "Estimate"]
    hmaxstar_hmax <- fixef(synthetic_fit3)["hmaxstar_hmax", "Estimate"]
    hmax <- mean(patch_hmax_mapping$hmax) * hmaxstar_hmax + hmaxstar_intercept + hOffset_estimate 
   
    
    # Convert to tibble for storing
    
    result_entry <- tibble(
      replicate = replicate_id
      ,t_range = t_range_name
      ,g_input = g_val
      ,h_noise_input = h_noise
      ,patch = names(g_patch)
      ,g_patch_estimate = g_patch
      ,hOffset_patch_estimate = hOffset_patch
      ,hmaxstar_estimate = hmax
    )
    
    # Append to results_df
    
    results_df <- bind_rows(results_df, result_entry)
  }
}

hmax_fit2 <- results_df

saveRDS(hmax_fit2, "03_work/analysis/synthetic_simulation/models/sensitivity_synthetic_sim_hmax_fit2.rds")
#results_df <- readRDS("03_work/analysis/synthetic_simulation/models/sensitivity_synthetic_sim_hmax_fit2.rds")


# calculate hmax (hmaxstar + hOffset)

results_df$hmax_patch <- results_df$hmaxstar_estimate + results_df$hOffset_patch_estimate

# Convert results_df to long format for plotting

plot_data <- results_df %>%
  pivot_longer(cols = c(g_patch_estimate, hOffset_patch_estimate, hmax_patch),
               names_to = "parameter",
               values_to = "estimate") %>%
  mutate(
    parameter = factor(parameter, levels = c("g_patch_estimate", "hOffset_patch_estimate", "hmax_patch"),
                       labels = c("g_patch", "hOffset_patch", "hmax_patch")),
    g = factor(g_input),
    h_noise = factor(h_noise_input),
    t_range = factor(t_range, levels = c("0-150", "0-35,150"))
  ) %>%
  
  mutate(g_label = case_when(
    g == 0.05 ~ "g 0.05",
    TRUE ~ as.character(g)
  )) %>%
  
  mutate(h_noise_label = factor(case_when(
    h_noise == 1 ~ "h noise 1",
    h_noise == 3 ~ "h noise 3",
    h_noise == 5 ~ "h noise 5",
    TRUE ~ as.character(h_noise)),
    levels = c("h noise 1", "h noise 3", "h noise 5")
  ))


# g

g_patch_plot <- ggplot(filter(plot_data, parameter == "g_patch"), aes(x = estimate, fill = t_range)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  geom_vline(aes(xintercept = as.numeric(g_input)), linetype = "dashed", size = 1) +
  facet_wrap(~g_label + h_noise_label) + #, scales = "free"
  theme_minimal() +
  labs(title = "Distribution of g_patch; hOffset = 10 m",
       x = "g_patch Estimate", y = "Count", fill = "t_range") +
  theme(legend.position = "bottom")
g_patch_plot

# hOffset

hOffset_patch_plot <- ggplot(filter(plot_data, parameter == "hOffset_patch"), aes(x = estimate, fill = t_range)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  #geom_vline(xintercept = 1, linetype = "dashed", size = 1) +  # Reference line at 1
  geom_vline(xintercept = 10, linetype = "dashed", size = 1) +  # Reference line at 10
  facet_wrap(~g_label + h_noise_label) + # , scales = "free"
  theme_minimal() +
  labs(title = "Distribution of hOffset_patch; hOffset = 10m",
       x = "hOffset_patch Estimate", y = "Count", fill = "t_range") +
  theme(legend.position = "bottom")
hOffset_patch_plot

# hmax

hmax_patch_plot <- ggplot(filter(plot_data, parameter == "hmax_patch"), aes(x = estimate, fill = t_range)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  #geom_vline(xintercept = 1, linetype = "dashed", size = 1) +  # Reference line at 1
  # geom_vline(xintercept = 10, linetype = "dashed", size = 1) +  # Reference line at 10
  facet_wrap(~g_label + h_noise_label) + # , scales = "free"
  theme_minimal() +
  labs(title = "Distribution of hmax_patch estimation; hOffset = 10m",
       x = "hOffset_patch Estimate", y = "Count", fill = "t_range") +
  theme(legend.position = "bottom")
hmax_patch_plot

hist(patch_hmax_mapping$hmax)



################################################################################
#
#  Vary g and hOffset by two management types with fixed h noise level 2 m
#
################################################################################

# Define the grid for g and h_noise

g_diff_values <- seq(0, 0.04, by = 0.02)  # g difference between mngt 1 and mngt 0 (from no difference to 0.04)
h_offset_diff_values <- seq(0, 10, by = 4)  # Varying Offset from 0 to 10

# Define the synthetic data settings

n_synthetic <- 800
n_fixed <- 20       # Number of rows with t = 150
n_random <- n_synthetic - n_fixed
num_replicates <- 3  # Number of repetitions

# prepare different hmax for different patch IDs to generate synthetic data; patches represent distrubance patches

selected_patches = c("100603", "100174", "100169", "1007344", "1002348", "1010147", "1010292", "1009814", "1007813", "100030", "100894", "1010434", "1009612",
                     "1010869", "1008131", "1008392", "1009341", "1010831", "1007956", "1009902", "1007126", "1005823", "100167", "1007497", "1009106", "1002333",
                     "1007860", "1009058", "1002598", "1007231", "101011", "1008853", "1004322", "1010448", "1010735", "100399", "1010605", "1007710", "1005271",
                     "1009749")

patch_hmax_mapping <- tibble(
  patch = selected_patches,
  hmax = c(31.52689, 29.60121, 31.10225, 16.64266, 30.91800, 32.62781, 30.14731, 29.48504, 30.38545, 31.12933, 23.68818, 32.59286, 26.38217, 34.31919, 26.10335,
           30.07681, 29.34345, 30.10626, 40.38669, 29.38597, 30.38545, 30.05920, 36.28186, 27.50942, 31.43734, 31.53233, 28.01574, 24.54890, 30.24358, 30.34470,
           28.09915, 31.21759, 26.95848, 29.72182, 27.35790, 29.56066, 30.91459, 33.10059, 27.50942, 30.07681)
)


# Create an empty dataframe to store results

results_df <- tibble()

# Start the repeated simulations

for (replicate_id in 1:num_replicates) {
  for (g_diff in g_diff_values) {  # Loop over g differences
    for (h_offset_diff in h_offset_diff_values) {  # Loop over h_noise_diff levels
      
      # Create the synthetic data
      
      synthetic_data <- tibble(
        t = c(runif(n_random, min = 0, max = 35), rep(150, n_fixed)),  # Generate a sequence of t
        patch = sample(selected_patches, n_synthetic, replace = TRUE)  # Randomly assign patches to rows
      ) %>%
        left_join(patch_hmax_mapping, by = "patch")
      
      patch_mngt_mapping <- tibble(
        patch = selected_patches,
        mngt = sample(0:1, length(selected_patches), replace = TRUE)  # Assign unique mngt to each patch
      )
      synthetic_data <- synthetic_data %>%
        left_join(patch_mngt_mapping, by = "patch")
      
      # Assign different g and hOffset  based on mngt 
      
      synthetic_data$g <- ifelse(synthetic_data$mngt == 1, 0.03 + g_diff, 0.05)  # g difference for mngt = 1 and mngt = 0
      
      synthetic_data$h <- synthetic_data$hmax * (1 - exp(-synthetic_data$g * synthetic_data$t))^3 + 
        rnorm(nrow(synthetic_data), mean = ifelse(synthetic_data$mngt == 1, 2 + h_offset_diff, 2), sd = 2)  # Varying the mean; here h noise 2 m
      
      # Fit the model
      
      tic(paste("Replicate:", replicate_id, "- Fitting model for g_diff =", g_diff, "and h_offset_diff =", h_offset_diff))
      
      synthetic_fit = brms::brm(
        bf(
          h ~ (hmax - hOffset) * (1 - exp(-g * t))^3 + hOffset,
          g ~ 1 + (1 | patch) + mngt,
          hOffset ~ 1 + (1 | patch) + mngt,
          nl = TRUE
        ),
        data = synthetic_data,
        family = "student",  
        prior = c(
          prior(normal(0.05, 0.02), lb = 0, nlpar = g, class = "b"),
          prior(normal(5,5), lb = 0, nlpar = hOffset, class = "b"),
          prior(normal(0.025, 0.025), lb = 0, nlpar = g, class = "sd"),
          prior(normal(10,10), lb = 0, nlpar = hOffset, class = "sd")
        ),
        iter = 2500, warmup = 1000, chains = 4, cores = 4,
        control = list(adapt_delta = 0.99, max_treedepth = 10),
        seed = 444,
        silent = FALSE
      )
      
      toc()
      
      # Extract estimates
      
      g_estimate <- fixef(synthetic_fit)["g_Intercept", "Estimate"]
      g_mngt1 <- fixef(synthetic_fit)["g_mngt", "Estimate"]
      
      ranef_g <- ranef(synthetic_fit)$patch[, , "g_Intercept"]
      ranef_g <- as.data.frame(ranef_g)
      ranef_g$patch <- rownames(ranef_g)
      ranef_g <- merge(ranef_g,  unique(synthetic_data[,c("patch","mngt")]), by = "patch", all.x = F)
      ranef_g$g_patch <- ifelse(ranef_g$mngt==0, (ranef_g$Estimate + g_estimate), (ranef_g$Estimate +g_estimate + g_mngt1))
      
      hOffset_estimate <- fixef(synthetic_fit)["hOffset_Intercept", "Estimate"]
      hOffset_mngt1 <- fixef(synthetic_fit)["hOffset_mngt", "Estimate"]
      
      ranef_hOffset <- ranef(synthetic_fit)$patch[, , "hOffset_Intercept"]
      ranef_hOffset <- as.data.frame(ranef_hOffset)
      ranef_hOffset$patch <- rownames(ranef_hOffset)
      ranef_hOffset <- merge(ranef_hOffset,  unique(synthetic_data[,c("patch","mngt")]), by = "patch", all.x = F)
      ranef_hOffset$hOffset_patch <- ifelse(ranef_hOffset$mngt==0, (ranef_hOffset$Estimate + hOffset_estimate), (ranef_hOffset$Estimate +hOffset_estimate + hOffset_mngt1))
      
      # Convert to tibble for storing
      
      result_entry <- tibble(
        replicate = replicate_id,
        g_diff = g_diff,
        h_offset_diff = h_offset_diff,
        patch = names(g_patch),
        g_patch_estimate = ranef_g$g_patch,
        hOffset_patch_estimate = ranef_hOffset$hOffset_patch,
        mngt = ranef_g$mngt,
      )
      
      # Append to results_df
      
      results_df <- bind_rows(results_df, result_entry)
    }
  }
}

df_g_offset_diff <- results_df

# -- Save results as RDS

saveRDS(df_g_offset_diff, "03_work/analysis/synthetic_simulation/models/sensitivity_synthetic_sim_reps_mngt_diff_hnoise2.rds")
#results_df <- readRDS("03_work/analysis/synthetic_simulation/models/sensitivity_synthetic_sim_reps_mngt_diff_hnoise2.rds")


# -- Convert results_df to long format for plotting

plot_data <- results_df %>%
  
  pivot_longer(cols = c(g_patch_estimate, hOffset_patch_estimate),
               names_to = "parameter",
               values_to = "estimate") %>%
  mutate(
    parameter = factor(parameter, levels = c("g_patch_estimate", "hOffset_patch_estimate"),
                       labels = c("g_patch_estimate", "hOffset_patch_estimate" )),
    g_diff = as.numeric(g_diff),
    h_offset_diff = as.numeric(h_offset_diff),
    mngt = factor(mngt)
  )  %>%
  
  mutate(g_diff_label = case_when(
    g_diff == 0 ~ "g0: 0.05 ; g1: 0.03", # ifelse(synthetic_data$mngt == 1, 0.03 + g_diff, 0.05) (diff -0.02)
    g_diff == 0.02 ~ "g1 & g2: 0.05",  #  (diff 0, as I add 0.2 onto the 0.3, which by default create the difference in g)
    g_diff == 0.04 ~ "g0: 0.05 ; g1: 0.07", # (diff +0.02)
    TRUE ~ as.character(g_diff)
  )) %>%
  
  mutate(hOffset_diff_label = case_when(
    h_offset_diff == 0 ~ "hOffset 0 & 1: 2m",
    h_offset_diff == 4 ~ "hOffset0: 2m; hOffset1: 6m",
    h_offset_diff == 8 ~ "hOffset0: 2m; hOffset1: 10m",
    TRUE ~ as.character(h_offset_diff)
  ))

# Define your label order explicitly (to ensure facets follow this order)

g_diff_levels <- c("g1 & g2: 0.05",
                   "g0: 0.05 ; g1: 0.03",
                   "g0: 0.05 ; g1: 0.07")

hOffset_levels <- c("hOffset 0 & 1: 2m",
                    "hOffset0: 2m; hOffset1: 6m",
                    "hOffset0: 2m; hOffset1: 10m")

plot_data <- plot_data %>%
  mutate(
    g_diff_label = factor(g_diff_label, levels = g_diff_levels),
    hOffset_diff_label = factor(hOffset_diff_label, levels = hOffset_levels)
  )

# -- plot g

g_mngt <- ggplot(filter(plot_data, str_detect(parameter, "g_patch_estimate")), aes(x = estimate, fill = mngt)) +
  geom_histogram(alpha = 0.6, bins = 50, position = "identity", colour = "black") +
  geom_vline(aes(xintercept = ifelse(mngt == 0, 0.05, 0.03 + g_diff), color = mngt), linetype = "dashed", size = 1) +
  facet_wrap(~ g_diff_label + hOffset_diff_label) + #, scales = "free"
  scale_fill_manual(values = c( "0" = "#3B82F6", "1" = "#F59E0B"),name = "Management") +
  scale_color_manual(values = c( "0" = "#3B82F6", "1" = "#F59E0B"),name = "Management") +
  theme_minimal() +
  labs(x = "g Estimate", y = "Count", fill = "mngt") +
  theme(legend.position = "bottom",
    strip.text = element_text(size = 13),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 12))                 

g_mngt

ggsave("03_work/analysis/synthetic_simulation/g_estimation_mngt_hnoise2_paperversion.tiff", plot = g_mngt, width = 11, height = 6.5, dpi = 300)



# -- plot hOffset

hOffset_mngt <- ggplot(filter(plot_data, str_detect(parameter, "hOffset_patch_estimate")), aes(x = estimate, fill = mngt)) +
  geom_histogram(alpha = 0.6, bins = 50, position = "identity", colour = "black") +
  geom_vline(aes(xintercept = ifelse(mngt == 0, 2, 2 + h_offset_diff), color = mngt), linetype = "dashed", size = 1) +
  facet_wrap(~ g_diff_label + hOffset_diff_label) + # ,scales = "free"
  scale_fill_manual(values = c( "0" = "#3B82F6", "1" = "#F59E0B"),name = "Management") +
  scale_color_manual(values = c( "0" = "#3B82F6", "1" = "#F59E0B"),name = "Management") +
  theme_minimal() +
  labs(x = "hOffset_patch Estimate", y = "Count", fill = "mngt") + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 13),  
        axis.title = element_text(size = 15),  
        axis.text = element_text(size = 12))  

hOffset_mngt

ggsave("03_work/analysis/synthetic_simulation/hOffset_estimation_mngt_hnoise2_paperversion.tiff", plot = hOffset_mngt, width = 10, height = 6, dpi = 300)


