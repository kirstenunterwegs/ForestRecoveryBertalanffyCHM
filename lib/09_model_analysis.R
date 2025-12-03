#####################################################################
#
# Model Evaluation & Recovery Simulation
#
#####################################################################


# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
#
# This script evaluates hierarchical Bayesian canopy-height growth models
# (asymmetric Laplace / student-t), summarises fixed and random effects,
# and generates management-focused model diagnostics and predictions.
#
# Specifically, it:
#   1. Loads the brms models and extracts posterior samples for
#      g, hOffset, and hmaxstar 
#   2. Summarises fixed vs. random effects as % change and plots component-wise
#      effect sizes for all predictors (growth rate, disturbance legacies,
#      asymptotic height).
#   3. Derives posterior distributions of g, hOffset and hmaxstar per
#      management type and forest type, and visualises them with half-eye plots.
#   4. Simulates canopy-height trajectories and plots mean growth curves
#      including observed data and highlighted time windows.
#   5. Quantifies time to reach a given canopy-height recovery threshold (e.g.
#      5 m) per management type and decomposes recovery into contributions from
#      disturbance legacies (hOffset) vs. subsequent growth (g).
#   6. Exports summary tables of fixed effects and random-effect SDs for all
#      model variants for use in the supplementary material.
#
# Output:
#   - Effect size panel plots (fixed vs. random, g / hOffset / hmaxstar)
#   - Posterior summaries per management type and species
#   - Simulated growth-curve plots with distributions at key time windows
#   - Recovery time and contribution plots for different management types
#   - Tables of fixed and random effects for supplementary material
# -------------------------------------------------------------------------


# ---libraries

library(brms)
library(data.table)

library(tidybayes)
library(posterior)
library(tidyr)
library(bayesplot)
library(dplyr)
library(stringr)

#for plotting
library(MetBrewer)
library(patchwork)
library(ggplot2)
library(ggridges)
library(sjPlot)


# set directory

setwd("~/NAS/Projects/ForestRecovery/")


# --------------------------------------
# load the model
# -------------------------------------


load("03_work/models/model_fits/fit_fe_sub250_cmdLOG_asymlap27.RData")
fit_fe_sub250_cmdLOG_asymlap27

# # student model
# load("03_work/models/model_fits/fit_fe_sub250_cmdLOG_student14.RData")
# fit_fe_sub250_cmdLOG_student14
#
# # model with spatial subset 500 m grid
# load("03_work/models//model_fits/fit_fe_sub250_cmdLOG_asymlap27_sub500.RData")
# fit_fe_sub250_cmdLOG_asymlap27_sub500
# 
# # model with spatial subset 1000 m grid
# load("03_work/models//model_fits/fit_fe_sub250_cmdLOG_asymlap27_sub1000.RData")
# fit_fe_sub250_cmdLOG_asymlap27_sub1000


# Set model to analyze

#---------------------------------------
model = fit_fe_sub250_cmdLOG_asymlap27
#model = fit_fe_sub250_cmdLOG_student14

# -- model checks
#model = fit_fe_sub250_cmdLOG_asymlap27_sub500
#model = fit_fe_sub250_cmdLOG_asymlap27_sub1000
#---------------------------------------


posterior_samples <- as_draws_df(model)

# Extract management type random effects 

mngt_labels <- c(
  "State.forest" = "State forest",
  "Other" = "Other",
  "Federal.forest" = "Federal forest",
  "Set.aside" = "Set aside",
  "Communal.forest" = "Municipal forest",
  "Small.private" = "Small private",
  "Large.private" = "Large private"
)

# --- growth rate per management type
posterior_g_mngt <- posterior_samples %>%
  select(starts_with("r_mngt_type__gLOG[")) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(mngt_type = gsub("r_mngt_type__gLOG\\[(.*),Intercept\\]", "\\1", parameter)) %>%
  mutate(g = posterior_samples$b_gLOG_Intercept + value) %>%
  select(g, mngt_type) %>%
  mutate(g = exp(g),
         mngt_type = recode(mngt_type, !!!mngt_labels)) 


# --- hOffset per management type
posterior_hOffset_mngt <- posterior_samples %>%
  select(starts_with("r_mngt_type__hOffsetLOG[")) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(mngt_type = gsub("r_mngt_type__hOffsetLOG\\[(.*),Intercept\\]", "\\1", parameter)) %>%
  mutate(hOffset = value + posterior_samples$b_hOffsetLOG_Intercept) %>%
  select(hOffset, mngt_type)%>%
  mutate(hOffset = exp(hOffset),
         mngt_type = recode(mngt_type, !!!mngt_labels)) 

# --- hmax per management type
posterior_hmax_mngt <- posterior_samples %>%
  select(starts_with("r_mngt_type__hmaxstarLOG[")) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(mngt_type = gsub("r_mngt_type__hmaxstarLOG\\[(.*),Intercept\\]", "\\1", parameter)) %>%
  mutate(hmaxstar = value + posterior_samples$b_hmaxstarLOG_Intercept) %>%
  select(hmaxstar, mngt_type)%>%
  mutate(hmaxstar = exp(hmaxstar),
         mngt_type = recode(mngt_type, !!!mngt_labels)) 



# ------------------------------------------
# define colour scheme for management types
# ------------------------------------------

mngt_colors <- c(
  "State forest" = "#A00E00",
  "Other"     = "#beaed4"  ,      
  "Federal forest" = "#DC7400" , 
  "Set aside"    = "#F6C200",
  "Municipal forest"     = "#529A6F",
  "Small private"    = "#132B69",
  "Large private" = "#e7298a"
)



# -----------------------------------------
# Model estimates for all predictors
# -----------------------------------------

### ------ version with posterior draws ------- ###

# --- Extract posterior draws ---

draws <- as_draws_df(model)


# model components to extract

components <- c("gLOG", "hOffsetLOG", "hmaxstarLOG")


# --- fixed Effects (exclude intercepts) ---

fixed_effects <- draws %>%
  select(matches("^b_")) %>%
  select(-matches("Intercept")) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "Draw") %>%
  mutate(
    Component = str_extract(Parameter, paste(components, collapse = "|")),
    Variable = str_remove(Parameter, paste0("^b_", Component, "_")),
    Component = str_remove(Component, "LOG"),
    Type = "Fixed",
    Draw_exp = (exp(Draw) - 1) * 100
  )


# summarise fixed effects

fixed_summary <- fixed_effects %>%
  group_by(Component, Variable, Type) %>%
  summarise(
    Estimate = median(Draw_exp),
    Q2.5 = quantile(Draw_exp, 0.025),
    Q97.5 = quantile(Draw_exp, 0.975),
    .groups = "drop"
  )


# --- random Effects (exclude intercepts) ---

r_cols <- grep("^r_", names(draws), value = TRUE)

random_effects <- draws %>%
  select(all_of(r_cols)) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "Draw") %>%
  filter(str_detect(Parameter, paste(components, collapse = "|"))) %>%
  mutate(
    Component = str_extract(Parameter, paste(components, collapse = "|")),
    Group = str_extract(Parameter, "(?<=r_)[^\\[]+"),  # group name: species, soil, etc.
    Level = str_extract(Parameter, "(?<=\\[)[^,]+"),   
    Component = str_remove(Component, "LOG"),
    Variable = paste0("sd_", Group),
    Type = "Random",
    Draw_exp = (exp(Draw) - 1) * 100
  )


# summarise SD across group levels per draw, then summarise SDs across draws

random_summary <- random_effects %>%
  group_by(Component, Variable, DrawIndex = rep(1:(n() / length(unique(Level))), each = length(unique(Level)))) %>%
  summarise(Draw_sd = sd(Draw_exp), .groups = "drop") %>%
  group_by(Component, Variable) %>%
  summarise(
    Estimate = median(Draw_sd),
    Q2.5 = quantile(Draw_sd, 0.025),
    Q97.5 = quantile(Draw_sd, 0.975),
    Type = "Random",
    .groups = "drop"
  )%>%
  mutate(
    Variable = Variable %>%
      str_remove("__gLOG|__hOffsetLOG|__hmaxstarLOG")
  )


# --- rename variables for plotting ---

var_labels <- c(
  "elevation_std" = "Elevation",
  "slope_std" = "Slope",
  "patch_ha_log_std" = "Disturbance size",
  "aspect_std" = "Aspect",
  "sd_mngt_type" = "Management type",
  "sd_mngt_type:species" = "Management type\nx Forest type",
  "sd_soil" = "Soil",
  "sd_species" = "Forest type"
)


# define function to prepare data for plotting

prep_plot_data <- function(df, component_name, effect_type, fixed_levels = NULL) {
  out <- df %>%
    dplyr::filter(Component == component_name, Type == effect_type) %>%
    dplyr::arrange(Estimate) %>%
    dplyr::mutate(Variable = dplyr::recode(Variable, !!!var_labels))
  if (is.null(fixed_levels)) {
    out <- out %>% dplyr::mutate(Variable = factor(Variable, levels = unique(Variable)))
  } else {
    out <- out %>% dplyr::mutate(Variable = factor(Variable, levels = fixed_levels))
  }
  out
}


# prepare data for fixed effects

plot_g_fixed     <- prep_plot_data(fixed_summary, "g", "Fixed")
fixed_order     <- levels(plot_g_fixed$Variable) # lock in order of fixed effects for the growth rate
plot_h_fixed    <- prep_plot_data(fixed_summary, "hOffset",  "Fixed", fixed_levels = fixed_order)
plot_hmax_fixed <- prep_plot_data(fixed_summary, "hmaxstar", "Fixed", fixed_levels = fixed_order)


# prepare data for random effects

plot_g_random     <- prep_plot_data(random_summary, "g", "Random")
random_order <- levels(plot_g_random$Variable)
plot_h_random     <- prep_plot_data(random_summary, "hOffset", "Random", fixed_levels = random_order)
plot_hmax_random  <- prep_plot_data(random_summary, "hmaxstar", "Random", fixed_levels = random_order)


# function to get min/max x-axis limits per component

get_x_limits <- function(fixed_df, random_df, component_name) {
  fixed_vals <- fixed_df %>% filter(Component == component_name) %>% select(Q2.5, Q97.5)
  random_vals <- random_df %>% filter(Component == component_name) %>% select(Q2.5, Q97.5)
  all_vals <- bind_rows(fixed_vals, random_vals)
  range(all_vals)
}

xlims_g <- get_x_limits(fixed_summary, random_summary, "g")
xlims_h <- get_x_limits(fixed_summary, random_summary, "hOffset")
xlims_hmax <- get_x_limits(fixed_summary, random_summary, "hmaxstar")


# function to create plots

make_effect_plot <- function(data, component_label, effect_label,
                             show_x = TRUE, show_y = FALSE,
                             xlim = NULL,
                             title_size = 18, y_text_size = 17,
                             y_title_size = 17,
                             panel_tag = NULL) {
  p <- ggplot(data, aes(x = Estimate, y = Variable, color = Type)) +
    geom_vline(xintercept = 0, color = "darkgrey", linewidth = 0.3) +
    geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5), linewidth = 1.2, size = 1) +
    coord_cartesian(xlim = xlim) +
    theme_minimal(base_size = 14) +
    labs(
      x = if (show_x) "Effect Size (% Change)" else NULL,
      y = if (show_y) effect_label else NULL,
      title = component_label
    ) +
    scale_color_manual(values = c("Fixed" = "black", "Random" = "grey40")) + #cadetblue
    theme(
      axis.text.y      = element_text(size = y_text_size),
      axis.text.x      = element_text(size = 16),
      axis.title.x     = element_text(size = if (show_x) 18 else 0),
      axis.title.y     = element_text(size = if (show_y) y_title_size else 0),
      plot.title       = element_text(size = title_size, hjust = 0),
      legend.position  = "none",
      panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.8),
      panel.background = element_blank()
    ) +
    { if (!show_x) theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title.x = element_blank()) else NULL } +
    { if (!show_y) theme(axis.text.y = element_blank(),
                         axis.ticks.y = element_blank(),
                         axis.title.y = element_blank()) else NULL }
  
  # add panel label
  if (!is.null(panel_tag)) {
    p <- p + annotate(
      "text",
      x = xlim[1] + diff(xlim) * 0.01,      
      y = Inf,                              
      label = panel_tag,
      hjust = 0.4, vjust = 2,
      size = 6
    )
  }
  
  return(p)
}


# --- plot growth rate variable estimates

pg_random    <- make_effect_plot(plot_g_random,    "Growth rate", "Random effects",
                                 show_x = FALSE, show_y = TRUE,  xlim = xlims_g, panel_tag = "(a)")
pg_fixed     <- make_effect_plot(plot_g_fixed,     "", "Fixed effects",
                                 show_x = TRUE, show_y = TRUE,  xlim = xlims_g, panel_tag = "(d)")


# --- plot disturbance legacies estimates

ph_random    <- make_effect_plot(plot_h_random,    "Disturbance legacies", "",
                                 show_x = FALSE, show_y = FALSE, xlim = xlims_h, panel_tag = "(b)")
ph_fixed     <- make_effect_plot(plot_h_fixed,     "", "",
                                 show_x = TRUE, show_y = FALSE, xlim = xlims_h, panel_tag = "(e)")


# --- plot asymptotic height estimates

phmax_random <- make_effect_plot(plot_hmax_random, "Asymptotic canopy height", "",
                                 show_x = FALSE, show_y = FALSE, xlim = xlims_hmax, panel_tag = "(c)")
phmax_fixed  <- make_effect_plot(plot_hmax_fixed,  "", "",
                                 show_x = TRUE, show_y = FALSE, xlim = xlims_hmax, panel_tag = "(f)")



overall_effects <-
  (pg_random | ph_random | phmax_random) /
  (pg_fixed  | ph_fixed  | phmax_fixed) 

overall_effects

ggsave("03_work/results/all_model_effects_transformed_asyml27_panels.tif", plot = overall_effects, dpi = 300, width = 14, height = 8)
ggsave("03_work/results/all_model_effects_transformed_asyml27_panels.jpg", plot = overall_effects, dpi = 300, width = 14, height = 8)

# ggsave("03_work/results/all_model_effects_transformed_student14_panels.tif", plot = overall_effects, dpi = 300, width = 15, height = 8)
# ggsave("03_work/results/all_model_effects_transformed_student14_panels.jpg", plot = overall_effects, dpi = 300, width = 15, height = 8)

#ggsave("03_work/results/results_control/all_model_effects_asyml27_sub500.tif", plot = overall_effects, dpi = 300, width = 14, height = 8)
#ggsave("03_work/results/results_control/all_model_effects_asyml27_sub1000.tif", plot = overall_effects, dpi = 300, width = 14, height = 8)



# --------------------------------------------------
# Plot posterior distributions per management type
# --------------------------------------------------


# set theme for posterior plots

theme_posterior <- theme_classic() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.1, color = "grey75"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 17),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 15, color = "grey30"),
    axis.ticks.y  = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_blank(),
    panel.grid.minor.x = element_line(color = "grey80", linewidth = 0.4, linetype = "dotted")
    )
  

# get n observation per group: 

n_by_mngt <- model$data %>%
  count(mngt_type, name = "n")%>%
  mutate(mngt_type = recode(mngt_type,
                       "Other" = "Other",
                       "Communal forest" = "Municipal forest"))


# ---- Summary for g ----

g_summary <- posterior_g_mngt %>%
  group_by(mngt_type) %>%
  summarise(g_median = median(g, na.rm = TRUE)) %>%
  left_join(n_by_mngt, by = "mngt_type") %>%
  arrange(g_median) %>%
  mutate(
    mngt_label = paste0(mngt_type, " (n = ", n, ")"),
    mngt_type = factor(mngt_type, levels = mngt_type),
    mngt_label = factor(mngt_label, levels = mngt_label)
  )


posterior_g_mngt <- posterior_g_mngt %>%
  left_join(g_summary %>% select(mngt_type, mngt_label), by = "mngt_type") %>%
  mutate(
    mngt_label = factor(mngt_label, levels = levels(g_summary$mngt_label))
  )


# apply same to hOffset, preserving order from g_summary

posterior_hOffset_mngt <- posterior_hOffset_mngt %>%
  left_join(g_summary %>% select(mngt_type, mngt_label), by = "mngt_type") %>%
  mutate(
    mngt_label = factor(mngt_label, levels = levels(g_summary$mngt_label))
  )


# apply same to hmaxstar, preserving order from g_summary
posterior_hmax_mngt <- posterior_hmax_mngt %>%
  left_join(g_summary %>% select(mngt_type, mngt_label), by = "mngt_type") %>%
  mutate(
    mngt_type = factor(mngt_type, levels = levels(g_summary$mngt_type)),
    mngt_label = factor(mngt_label, levels = levels(g_summary$mngt_label))
  )


# global medians for orientation vertical line

median_g <- median(posterior_g_mngt$g, na.rm = TRUE)
median_hOffset <- median(posterior_hOffset_mngt$hOffset, na.rm = TRUE)
median_hmax <- median(posterior_hmax_mngt$hmaxstar, na.rm = TRUE)



# ---- plot g ----

p_g_mngt <- ggplot(posterior_g_mngt, aes(x = mngt_label, y = g, fill = mngt_type, color = mngt_type)) +
  stat_halfeye(alpha = 0.3, fill_type = "segments", show.legend = FALSE) +
  stat_interval(.width = 0.95, size = 1.2, alpha = 0.4, show.legend = FALSE) +
  stat_interval(.width = 0.80, size = 2.0, alpha = 0.7, show.legend = FALSE) +
  stat_interval(.width = 0.50, size = 3.0, alpha = 1, show.legend = FALSE) +
  stat_summary(geom = "point", fun = median, size = 2.5, shape = 21, fill = "black", color = "black", stroke = 0.6) +
  geom_hline(yintercept = median_g, color = "grey30", linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = mngt_colors) +
  scale_color_manual(values = mngt_colors) +
  coord_flip(ylim = c(0.045, 0.09), clip = "off") +
  labs(x = "Management Type", y = "Height growth rate") +
  theme_posterior+
  theme(    axis.title.y = element_text(margin = margin(r = 15)))

# ---- plot hOffset ----

p_hOffset_mngt <- ggplot(posterior_hOffset_mngt, aes(x = mngt_label, y = hOffset, fill = mngt_type, color = mngt_type)) +
  stat_halfeye(alpha = 0.3, fill_type = "segments", show.legend = FALSE) +
  stat_interval(.width = 0.95, size = 1.2, alpha = 0.4, show.legend = FALSE) +
  stat_interval(.width = 0.80, size = 2.0, alpha = 0.7, show.legend = FALSE) +
  stat_interval(.width = 0.50, size = 3.0, alpha = 1, show.legend = FALSE) +
  stat_summary(geom = "point", fun = median, size = 2.5, shape = 21, fill = "black", color = "black", stroke = 0.6) +
  geom_hline(yintercept = median_hOffset, color = "grey30", linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = mngt_colors) +
  scale_color_manual(values = mngt_colors) +
  coord_flip(ylim = c(0, 2), clip = "off") + # asym laplace
  #coord_flip(ylim = c(0, 4), clip = "off") + # asym laplace sub 1000
  #coord_flip(ylim = c(1.5, 7), clip = "off") + #student
  labs(x = "", y = "Disturbance legacies [m]") +
  theme_posterior+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# ---- plot hmaxstar ----

p_hmax_mngt <- ggplot(posterior_hmax_mngt, aes(x = mngt_label, y = hmaxstar, fill = mngt_type, color = mngt_type)) +
  stat_halfeye(alpha = 0.3, fill_type = "segments", show.legend = FALSE) +
  stat_interval(.width = 0.95, size = 1.2, alpha = 0.4, show.legend = FALSE) +
  stat_interval(.width = 0.80, size = 2.0, alpha = 0.7, show.legend = FALSE) +
  stat_interval(.width = 0.50, size = 3.0, alpha = 1, show.legend = FALSE) +
  stat_summary(geom = "point", fun = median, size = 2.5, shape = 21,
               fill = "black", color = "black", stroke = 0.6) +
  geom_hline(yintercept = median_hmax, color = "grey30", linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = mngt_colors) +
  scale_color_manual(values = mngt_colors) +
  coord_flip(ylim = c(13, 25), clip = "off") + # asym laplace
  #coord_flip(ylim = c(9, 23), clip = "off") + # student
  labs(x = "", y = "Asymptotic canopy \nheight [m]") +
  theme_posterior +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# ---- Combine plots ----

p_mngt_estimates <- p_g_mngt + p_hOffset_mngt + p_hmax_mngt
p_mngt_estimates


ggsave("03_work/results/mngt_g_hOffset_hmaxstar_asyml_27.tif", plot = p_mngt_estimates, dpi = 300, width = 12, height = 6)
#ggsave("03_work/results/mngt_g_hOffset_hmaxstar_student14.tif", plot = p_mngt_estimates, dpi = 300, width = 12, height = 6)

#ggsave("03_work/results/results_control/mngt_g_hOffset_asyml_27_sub500.tif", plot = p_mngt_estimates, dpi = 300, width = 12, height = 6)
#ggsave("03_work/results/results_control/mngt_g_hOffset_asyml_27_sub1000.tif", plot = p_mngt_estimates, dpi = 300, width = 12, height = 6)

# plots with only g and hOffset (no asymptote)
# ggsave("03_work/results/mngt_g_hOffset_asyml_27.tif", plot = p_mngt_estimates, dpi = 300, width = 12, height = 6)
# ggsave("03_work/results/mngt_g_hOffset_student14.tif", plot = p_mngt_estimates, dpi = 300, width = 12, height = 6)


# --- summaries median per mngt

hOffset_summary <- posterior_hOffset_mngt %>%
  group_by(mngt_type) %>%
  summarise(
    hOffset_median = median(hOffset, na.rm = TRUE)
  )


hmax_summary <- posterior_hmax_mngt %>%
  group_by(mngt_type) %>%
  summarise(
    hmax_median = median(hmaxstar, na.rm = TRUE)
  )



#------------------------------------------
# Simulate growth curves from the model
#------------------------------------------

# Set fixed values
t_seq <- round(seq(0, 125, length.out = 100),0)  # time sequence

hmax_val <- exp(fixef(model)["hmaxstarLOG_Intercept", "Estimate"])

# Get levels from the model data
species_levels = unique(model$data$species)
mngt_levels = unique(model$data$mngt_type)
mngt_species_levels = unique(model$data$`mngt_type:species`)
soil_levels = unique(model$data$soil)


# Create new data grid

new_data <- expand_grid(
  t = t_seq,
  species = species_levels,
  mngt_type = mngt_levels,
  soil = soil_levels
) %>%
  mutate(
    slope_std = 0,
    aspect_std = 0,
    elevation_std = 0,
    patch_ha_log_std = 0,
  )

setDT(new_data)

new_data$logt = log10(new_data$t)
new_data = new_data[t>1]
new_data[, t_alt := t]
new_data[, t_alt := fifelse(t == 250, 35, t)]
new_data[, log10t := log10(t_alt)]


# identify valid mngt:species combinations in the data for model fit and filter for those in the new data dt

valid_combos <- model$data %>% distinct(mngt_type, soil, species) #
new_data_filtered <- new_data %>% semi_join(valid_combos, by = c("mngt_type", "soil", "species")) #
setDT(new_data_filtered)


# chunking for efficiency

new_data_filtered[, chunk := floor(1:nrow(new_data_filtered) / 1000)]
chunks = unique(new_data_filtered$chunk)
hpred = vector(mode = "list", length = length(chunks))


# List to store prediction draws
draws_list <- vector("list", length(chunks))


# Collect posterior draws for each chunk

for (i in seq_along(chunks)) {
  chunk_id <- chunks[i]
  cat("Processing chunk", chunk_id, "\n")
  chunk_data <- new_data_filtered[chunk == chunk_id]
  
  draws <- posterior_linpred(model, newdata = chunk_data, re_formula = ~ (1 | mngt_type)) 
  #draws <- posterior_epred(model, newdata = chunk_data, re_formula = ~ (1 | mngt_type)) # for student model
  
  hpred_draws = colMeans(draws)

  draws_list[[i]] <- hpred_draws
}

hpred= unlist(draws_list)

new_data_filtered$hpred = hpred


#save(new_data_filtered, file = "03_work/analysis/growth_curve_simulations/intermediate_results/predictions_asyml27_growth_curves_mngt_nointeract.RData")
#save(new_data_filtered, file = "03_work/analysis/growth_curve_simulations/predictions_student14_growth_curves_mngtonly.RData")

#load("03_work/analysis/growth_curve_simulations/predictions_asyml27_reNULL_growth_curves_mngtonly.RData")
#load("03_work/analysis/growth_curve_simulations/predictions_student14_growth_curves_mngtonly.RData")


plot_summary = new_data_filtered %>% 
  group_by(t, mngt_type) %>%
  summarise(
    h_median = median(hpred, na.rm = TRUE),
    h_mean = mean(hpred, na.rm=TRUE),
    .groups = "drop"
  )%>%
  mutate(
    mngt_type = recode(mngt_type,
                       "Other" = "Other",
                       "Communal forest" = "Municipal forest"
    )
  )



# -----------------------------------------
# plot growth curves 
# -----------------------------------------

# prepare data for plotting:

# only plot until 90 years, as here asymptote already reached
plot_summary_sub = plot_summary %>% filter(t <91)

# prepare raw data points to add to the plot
set.seed(12) 

# subsetting observations to x samples per mngt type - t combination for display:

raw_points_sub <- model$data %>%
  mutate(
    mngt_type = recode(mngt_type,
                       "Other"            = "Other",
                       "Communal forest"  = "Municipal forest"),
    t = ifelse(t == 250, 89, t)
  ) %>%
  select(t, h_mean, mngt_type) %>%
  mutate(mngt_type = factor(mngt_type, levels = names(mngt_colors))) %>%
  group_by(mngt_type, t) %>%
  # sample up to 40 obs per (mngt_type, t)
  group_modify(~ dplyr::slice_sample(.x, n = min(40, nrow(.x)))) %>% 
  ungroup() %>%
  mutate(mngt_type = factor(mngt_type, levels = levels(plot_summary_sub$mngt_type)))


# --- plot growth curves ---

growth_curve_mngt = ggplot(plot_summary_sub, aes(x = t, y = h_median, color = mngt_type)) +
  
  # Raw data points
  geom_jitter(
    data = raw_points_sub,
    aes(x = t, y = h_mean, color = mngt_type),
    inherit.aes = FALSE,
    alpha = 0.3, size = 0.3,
    width = 0.4) +
  
  # Rectangles for highlighted time windows
  annotate("rect", xmin = 3, xmax = 7, ymin = 0, ymax = 4, fill = "lightgrey", alpha = 0.7, color = "darkgrey", linewidth=0.8) +
  annotate("rect", xmin = 18, xmax = 22, ymin = 7, ymax = 11, fill = "lightgrey", alpha = 0.7, color = "darkgrey", linewidth=0.8) +
  annotate("rect", xmin = 48, xmax = 52, ymin = 16, ymax = 20, fill = "lightgrey", alpha = 0.7, color = "darkgrey", linewidth=0.8) +
  
  
  #---- Retangles for student version
  #Rectangles for highlighted time windows
  # annotate("rect", xmin = 2, xmax = 8, ymin = 2, ymax = 6, fill = "lightgrey", alpha = 0.5, color = "darkgrey") +
  # annotate("rect", xmin = 22, xmax = 28, ymin = 9, ymax = 13, fill = "grey", alpha = 0.5, color = "darkgrey") +
  # annotate("rect", xmin = 47, xmax = 53, ymin = 14, ymax = 18, fill = "darkgrey", alpha = 0.5, color = "darkgrey") +
  #----------------


  # Main data layers
  geom_line(size = 1) +
  
  
  # # --- Rectangle labels in top left corner of each
  annotate("text", x = 3.2, y = 3.8, label = "A", hjust = 0, vjust = 1, size = 5) +
  annotate("text", x = 18.2, y = 10.8, label = "B", hjust = 0, vjust = 1, size = 5) +
  annotate("text", x = 48.2, y = 19.8, label = "C", hjust = 0, vjust = 1, size = 5) +
  
 
   #---- Labels for student version
  # annotate("text", x = 2.2, y = 5.8, label = "A", hjust = 0, vjust = 1, size = 5) +
  # annotate("text", x = 22.2, y = 12.8, label = "B", hjust = 0, vjust = 1, size = 5) +
  # annotate("text", x = 47.2, y = 17.8, label = "C", hjust = 0, vjust = 1, size = 5) +
  #----------------
  
  # geom_ribbon(aes(ymin = h_lower, ymax = h_upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = mngt_colors) +
  
  labs(
    x = "Time since disturbance (years)",
    y = "Canopy height (m)",
    color = "Management type",
    fill = "Management type"
  ) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_continuous(
    expand = c(0.01, 0.01),
    breaks = c(20, 40, 60, 80, 89),                  # where to put the ticks
    labels = c("20", "40", "60", "80", "> 100")         # what to show there
  ) +
  theme_classic() +
  theme(
    legend.position = "top", #c(0.90, 0.03)
    legend.justification = c("right", "bottom"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 11),
    legend.text  = element_text(size = 12),
    legend.title  = element_text(size = 14)) +

  guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8)))

growth_curve_mngt 


# ---  plot distributions of simulated heights per defined time window ---


#---- adjust management labels: 

new_data_filtered = new_data_filtered %>% 
  mutate( mngt_type = recode(mngt_type,
                       "Other" = "Other",
                       "Communal forest" = "Municipal forest" ))


#--- subplot A

draws_filtered <- new_data_filtered %>%
  filter(t ==5)


p_A <- ggplot(draws_filtered, aes(x = factor(mngt_type), y = hpred, fill = mngt_type)) +
  geom_boxplot(alpha = 0.5, width = 1) +  # Width = 1 fills entire slot
  scale_fill_manual(values = mngt_colors) +
  scale_y_continuous(position = "right") +
  labs(x = NULL, y = "Predicted Height (m)", title = "A") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.7),
    panel.grid.minor.y = element_line(color = "grey80", linewidth = 0.4, linetype = "dotted"))

p_A


#--- subplot B

draws_filtered <- new_data_filtered %>%
  filter(t ==20)


p_B <- ggplot(draws_filtered, aes(x = factor(mngt_type), y = hpred, fill = mngt_type)) +
  geom_boxplot(alpha = 0.5, width = 1) +  # Width = 1 fills entire slot
  scale_fill_manual(values = mngt_colors) +
  scale_y_continuous(
    position = "right",
    breaks = seq(1, 30, by = 1),               
    minor_breaks = seq(0.5, 30, by = 0.5)  
  ) + 
  labs(x = NULL, y = "Predicted Height (m)", title = "B") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.7),
    panel.grid.minor.y = element_line(color = "grey80", linewidth = 0.4, linetype = "dotted"))

p_B


#--- subplot C

draws_filtered <- new_data_filtered %>%
  filter(t == 49)


p_C <- ggplot(draws_filtered, aes(x = factor(mngt_type), y = hpred, fill = mngt_type)) +
  geom_boxplot(alpha = 0.5, width = 1) +  # Width = 1 fills entire slot
  scale_fill_manual(values = mngt_colors) +
  scale_y_continuous(
    position = "right",
    breaks = seq(1, 30, by = 1),              
    minor_breaks = seq(0.5, 30, by = 0.5)     
  ) + 
  labs(x = NULL, y = "Predicted Height (m)", title = "C") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.7),
    panel.grid.minor.y = element_line(color = "grey80", linewidth = 0.4, linetype = "dotted"))

p_C


# --- formatting for Legend

growth_curve_mngt <- growth_curve_mngt +
  guides(
    color = guide_legend(
      nrow = 1,                    # force into one line
      title.position = "top",   
      title.hjust = 0.5,        
      byrow = TRUE,
      override.aes = list(size = 4, alpha = 0.8)
    ),
    fill = guide_legend(nrow = 1, byrow = TRUE)
  )  +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 13, margin = margin(l = 2, r = 10)), 
    legend.box.spacing = unit(0.3, "cm"),   
    legend.spacing.x  = unit(0.5, "cm"),    
    legend.key.width  = unit(0.8, "cm")     
  )

# Extract legend from the plot
legend_mngt <- cowplot::get_legend(growth_curve_mngt)

# remove legend from the growth_curve_mngt panel itself
growth_curve_mngt_nolegend <- growth_curve_mngt +
  theme(legend.position = "none")

# combine plots (growth curve + subplots A/B/C)
combined_plot_core <- growth_curve_mngt_nolegend + p_C / p_B / p_A +
  plot_layout(ncol = 2, widths = c(2, 1))

# assemble final figure with legend at top
final_combined_plot <- cowplot::plot_grid(
  legend_mngt,
  combined_plot_core,
  ncol = 1,
  rel_heights = c(0.1, 1)  # adjust vertical space for legend
)

final_combined_plot

ggsave("03_work/results/growth_curves_mngt_distributions_sub_asymlaplace.tif", plot = final_combined_plot, dpi = 300, width = 13, height = 7)
#ggsave("03_work/results/growth_curves_mngt_distributions_sub_student14.jpg", plot = final_combined_plot, dpi = 300, width = 13, height = 7)



###--- calculate IQR per management type for the different time steps (and overall) ---

variation_summary <- new_data_filtered %>%
  filter(t<91) %>% # only include until 90, as growth curve flattens of and I display data until this timestep
  group_by(mngt_type, t) %>%
  summarise(
    IQR = IQR(hpred, na.rm = TRUE),
    SD = sd(hpred, na.rm = TRUE),
    Mean = mean(hpred, na.rm = TRUE),
    CV = SD / Mean
  ) %>%
  ungroup()

variation_summary_overall <- variation_summary %>%
  group_by(mngt_type) %>%
  summarise(
    mean_IQR = mean(IQR),
    mean_SD = mean(SD),
    mean_CV = mean(CV)
  )

  

# --------------------------------------------------------
# Contribution of g and hOffset to recovery 5 m
# --------------------------------------------------------

h_recovery <- 5 # recovery threshold 5 m


posterior_g_mngt$.draw <- 1:nrow(posterior_g_mngt)
posterior_hOffset_mngt$.draw <- 1:nrow(posterior_hOffset_mngt)



posterior_hmaxstar_mngt <- posterior_samples %>%
  select(starts_with("r_mngt_type__hmaxstarLOG[")) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(
    mngt_type = gsub("r_mngt_type__hmaxstarLOG\\[(.*),Intercept\\]", "\\1", parameter),
    hmaxstar_log = posterior_samples$b_hmaxstarLOG_Intercept + value,
    hmaxstar = exp(hmaxstar_log)
  ) %>%
  select(hmaxstar, mngt_type) %>%
  mutate(
    mngt_type = recode(mngt_type, !!!mngt_labels)
  )


posterior_hmaxstar_mngt$.draw <- 1:nrow(posterior_hmaxstar_mngt)


posterior_contrib_mngt <- posterior_g_mngt %>%
  left_join(posterior_hOffset_mngt, by = c(".draw", "mngt_type")) %>%
  left_join(posterior_hmaxstar_mngt, by = c(".draw", "mngt_type")) %>%
  mutate(
    h_recovery = h_recovery, 
    t_recovery = -(1/g) * log(1- ((h_recovery - hOffset)/ hmaxstar)^ (1/3)),
    h_growth = h_recovery - hOffset,
    hOffset_contrib = hOffset / h_recovery,
    growth_contrib = h_growth / h_recovery
  ) %>%
  select(.draw, mngt_type, hOffset, g, hOffset_contrib, growth_contrib, h_recovery, t_recovery)


# when hOffset > h_recovery, that set it's contribution = 1, growth contribution = 0 & t_recovery = 0

posterior_contrib_mngt$hOffset_contrib = ifelse(posterior_contrib_mngt$hOffset >= h_recovery, 1 , posterior_contrib_mngt$hOffset_contrib)
posterior_contrib_mngt$growth_contrib = ifelse(posterior_contrib_mngt$hOffset >= h_recovery, 0 , posterior_contrib_mngt$growth_contrib)
posterior_contrib_mngt$t_recovery = ifelse(posterior_contrib_mngt$hOffset >= h_recovery, 0 , posterior_contrib_mngt$t_recovery)


posterior_contrib_mngt_long <- posterior_contrib_mngt %>%
  pivot_longer(cols = c("hOffset_contrib", "growth_contrib"),
               names_to = "component", values_to = "proportion")


# mean & median contributions

mean_hOffset_overall = mean(posterior_contrib_mngt$hOffset_contrib)
median_hOffset_overall = median(posterior_contrib_mngt$hOffset_contrib)
mean_g_overall = mean(posterior_contrib_mngt$growth_contrib)
median_g_overall = median(posterior_contrib_mngt$growth_contrib)


# Calculate median contributions and median time to recovery for each management type

posterior_contrib_mngt_medians <- posterior_contrib_mngt %>%
  group_by(mngt_type) %>%
  summarise(
    median_growth_contrib = median(growth_contrib, na.rm = TRUE),
    median_hOffset_contrib = median(hOffset_contrib, na.rm = TRUE),
    median_t_recovery = median(t_recovery, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c("median_growth_contrib", "median_hOffset_contrib"),
               names_to = "component", values_to = "median_proportion")


# --- combined plot ---

# Reorder based on median recovery time

mngt_order <- posterior_contrib_mngt %>%
  group_by(mngt_type) %>%
  summarise(median_t_recovery = median(t_recovery, na.rm = TRUE)) %>%
  arrange(median_t_recovery) %>%
  pull(mngt_type)


# Apply factor levels based on order

posterior_contrib_mngt$mngt_type <- factor(posterior_contrib_mngt$mngt_type, levels = mngt_order)
posterior_contrib_mngt_medians$mngt_type <- factor(posterior_contrib_mngt_medians$mngt_type, levels = mngt_order)
posterior_contrib_mngt_long$mngt_type <- factor(posterior_contrib_mngt_long$mngt_type, levels = mngt_order)

mean_recovery_time_mngt = mean(posterior_contrib_mngt$t_recovery, na.rm=T)
median_recovery_time_mngt = median(posterior_contrib_mngt$t_recovery, na.rm=T)


# recovery time plot

recov_time_mngt = ggplot(posterior_contrib_mngt, aes(x = mngt_type, y = t_recovery, fill = mngt_type)) +
  scale_fill_manual(values = mngt_colors) +
  scale_color_manual(values = mngt_colors) +
  geom_violin(alpha = 0.5, color = "grey30", linewidth = 0.1) +
  geom_boxplot(
    aes(x = mngt_type, y = t_recovery),
    width = 0.1,
    alpha = 0.2,
    color = "grey30",
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  )+
  labs( x = "Management type",
    y =  "Time to recovery (years)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_flip()
recov_time_mngt


posterior_contrib_mngt_medians = posterior_contrib_mngt_medians %>%
  mutate(median_proportion = median_proportion * 100)


# contribution plot

median_mngt_g_hOffset_contrib_flip = ggplot(posterior_contrib_mngt_medians, aes(x = mngt_type, y = median_proportion, fill = component)) +
  geom_bar(stat = "identity", width = 0.6, color = "black",) +  
  scale_fill_manual(
    values = c("white", "grey30"),  
    labels = c("Height Growth","Distrubance legacies" ),
    guide = guide_legend(ncol = 1)) +
  labs(
    x = NULL,
    y = "Median contribution to recovery (%)",
    fill = "Component"
  ) +
  theme_tidybayes() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 14),  
    legend.text  = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12)
  ) +
  coord_flip()#+


recov_contrib_combined = (recov_time_mngt + median_mngt_g_hOffset_contrib_flip) +
  plot_layout(ncol = 2, widths = c(1, 0.5))

recov_contrib_combined



ggsave("03_work/results/recov_time_contrib_5m_asyml_27.tif", plot = recov_contrib_combined, dpi = 300, width = 12, height = 6)
#ggsave("03_work/results/recov_time_contrib_5m_student14.tif", plot = recov_contrib_combined, dpi = 300, width = 12, height = 6)
#ggsave("03_work/results_control/recov_time_contrib_5m_asyml_27_sub500.tif", plot = recov_contrib_combined, dpi = 300, width = 12, height = 6)
#ggsave("03_work/results_control/recov_time_contrib_5m_asyml_27_sub1000.tif", plot = recov_contrib_combined, dpi = 300, width = 12, height = 6)



# -------------------------------------------------
# print model tables for supplementary material 
# -------------------------------------------------


# load all models if not done yet:

load("03_work/models/model_fits/fit_fe_sub250_cmdLOG_asymlap27.RData")
fit_fe_sub250_cmdLOG_asymlap27

# student model
load("03_work/models/model_fits/fit_fe_sub250_cmdLOG_student14.RData")
fit_fe_sub250_cmdLOG_student14

load("03_work/models//model_fits/fit_fe_sub250_cmdLOG_asymlap27_sub500.RData")
fit_fe_sub250_cmdLOG_asymlap27_sub500

load("03_work/models//model_fits/fit_fe_sub250_cmdLOG_asymlap27_sub1000.RData")
fit_fe_sub250_cmdLOG_asymlap27_sub1000


# chose model to create table for:

fit <- fit_fe_sub250_cmdLOG_asymlap27
#fit <- fit_fe_sub250_cmdLOG_student14
#fit = fit_fe_sub250_cmdLOG_asymlap27_sub500
#fit = fit_fe_sub250_cmdLOG_asymlap27_sub1000


vc  <- VarCorr(fit, summary = TRUE)

# -------- SDs (keep Estimate, Est.Error, Q2.5, Q97.5) ----------

get_vc_sd_df <- function(vc_list) {
  out <- list()
  for (g in names(vc_list)) {
    sds <- vc_list[[g]]$sd
    if (is.null(sds)) next
    
    # If a plain numeric vector (rare), wrap it to the expected columns
    if (is.null(dim(sds))) {
      df <- data.frame(
        Term      = if (!is.null(names(sds))) names(sds) else rep("Intercept", length(sds)),
        Estimate  = as.numeric(sds),
        Est.Error = NA_real_,
        Q2.5      = NA_real_,
        Q97.5     = NA_real_,
        check.names = FALSE
      )
    } else {
      # matrix with the 4 columns od estimates
      keep <- intersect(colnames(sds), c("Estimate","Est.Error","Q2.5","Q97.5"))
      df <- as.data.frame(sds[, keep, drop = FALSE], check.names = FALSE)
      df$Term <- rownames(sds)
      # reorder columns
      df <- df[, c("Term","Estimate","Est.Error","Q2.5","Q97.5")]
    }
    
    df$Group <- g
    out[[length(out) + 1]] <- df[, c("Group","Term","Estimate","Est.Error","Q2.5","Q97.5")]
  }
  do.call(rbind, out)
}


vc_sd_df <- get_vc_sd_df(vc)


# ---- Fixed effects

fx <- fixef(fit, probs = c(.025, .975))
fixed_tab <- data.frame(
  Term     = rownames(fx),
  Estimate = fx[, "Estimate"],
  SE       = fx[, "Est.Error"],
  `2.5%`   = fx[, "Q2.5"],
  `97.5%`  = fx[, "Q97.5"],
  row.names = NULL, check.names = FALSE
)

# ---- Render nice tables

tab_df(fixed_tab, title = "Fixed effects (posterior mean & 95% CI)")
tab_df(vc_sd_df,  title = "Random-effects SDs (posterior mean, SE, 95% CI)")

if (!is.null(vc_cor)) tab_df(vc_cor, title = "Random-effects correlations (VarCorr)")


################################################ END