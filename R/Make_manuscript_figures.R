rm(list = ls())
library(INLA)
library(INLAutils)
library(brinla)
library(tidyverse)
library(GGally)
library(raster)
library(readxl)
library(ggpubr)
library(patchwork)
library(ggthemes)
library(ztable)
library(magrittr)
library(viridis)

# -----------------------------------------------------------------------------------------------------------------
# Data import -----------
# -----------------------------------------------------------------------------------------------------------------
source("source/misc_functions.R")
load('model_data/regression_data.RData')
load("results/traits_model_list.RData")

figure_data <- model_data
col_list <- c('BG_imp', 'rain_imp', 'temp_imp', 'dry_imp', 'HFP_imp', 
              'rain_breadth', 'temp_breadth', 'dry_breadth', 'HFP_breadth', 'BG_breadth',
              'avg.r', 'HWI', 'Mass')

for (i in 1:length(col_list)){
      figure_data[, col_list[i]] <- unscale_fun(figure_data[, col_list[i]])
}

relative_scale <- scale_y_continuous(breaks = log(c(seq(0.5, 3, by = 0.5), 4, 5, 10, 20, 30, 40)), labels = c(seq(0.5, 3, by = 0.5), 4, 5, 10, 20, 30, 40))

# -----------------------------------------------------------------------------------------------------------------
## Overview of range changes -----------
# -----------------------------------------------------------------------------------------------------------------
cells_change_log <- ggplot(model_data, aes(x = reorder(species, log_prop_change))) +
      geom_bar(aes(y = log_prop_change, col = "Total range change"), stat = 'identity', alpha = 0.3) +
      geom_point(aes(y = log_relative_colonisation, col = "Colonisation"), pch = '-', size = 7) +
      geom_point(aes(y = log_relative_extinction, col = "Extinction"), pch = '-', size = 7) +
      geom_segment(aes(xend = reorder(species,  log_relative_colonisation), y = log_relative_extinction, yend = log_relative_colonisation), alpha = .3, lty = 2, lwd = .3) +
      geom_hline(aes(yintercept = 0), alpha = .5) + 
      theme_classic() +
      xlab("Species") +
      ylab("Relative change factor") +
      theme(legend.title = element_blank(),
            legend.position = c(.35, .9),
            legend.key = element_blank(),
            text = element_text(size = 16),
            legend.text = element_text(size = 14),
            legend.background = element_blank()) +
      scale_x_discrete(labels = element_blank()) +
      scale_colour_colorblind() +
      scale_y_continuous(breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40))

cor_2_rsq <- round(rsq(model_data$log_prop_change, model_data$log_relative_colonisation), digits = 2)
cor_3_rsq <- round(rsq(model_data$log_prop_change, model_data$log_relative_extinction), digits = 2)
coeff <- 2
cor_plot <- ggplot(model_data, aes(log_prop_change)) +
      geom_point(aes(y = log_relative_colonisation), col = "#000000") +
      geom_point(aes(y = log_relative_extinction*coeff), col = "#E69F00") +
      theme_classic() +
      xlab("Total range change") +
      ylab("Colonisation") +
      scale_colour_colorblind(guide = "none") +
      geom_text(aes(x = Inf, y = -Inf, hjust = 1.8, vjust = -2, label = paste0('R²= ', cor_2_rsq)), col = "#000000", size = 5)  +
      geom_text(aes(x = Inf, y = -Inf, hjust = 3, vjust = -2, label = paste0('R²= ', cor_3_rsq)), col = "#E69F00", size = 5) +
      scale_y_continuous(name = "Colonisation", breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40),
                         sec.axis = sec_axis(~./coeff, name="Extinction",
                                             breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)))  +
      scale_x_continuous(breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)) + 
      theme(axis.line.y.right = element_line(color = "#E69F00"), 
            axis.ticks.y.right = element_line(color = "#E69F00"),
            axis.text.y.right = element_text(color = "#E69F00"),
            axis.title.y.right = element_text(color = "#E69F00"),
            text = element_text(size = 16))

load('model_output/Ardeotis_kori_figure_files.RData')
preds <- rbind(preds_1, preds_2) %>% 
      mutate(ind = ifelse(ind == 1, "1980-1999", "2000-2020"))
median_plot_kori <- ggplot() +
      geom_tile(data = preds, aes(x = x, y = y, fill = all_pred)) +
      geom_polygon(data = TZ_no_lakes, aes(x=long, y=lat, group=group), 
                   fill = NA, color = "black") +
      geom_point(data = atlas_sp_simple,
                 aes(x = x, y = y, col = as.factor(presence), alpha = as.factor(presence), shape = as.factor(presence)), pch = 15, cex = 3) +
      geom_point(data = ebird_sp_simple %>% filter(presence == 0),
                 aes(x = x, y = y, shape = as.factor(presence)), col = "orangered", alpha = 0.3) +
      geom_point(data = ebird_sp_simple %>% filter(presence == 1),
                 aes(x = x, y = y, shape = as.factor(presence)), col = "black") +
      coord_sf() +
      facet_wrap(.~ ind, ncol = 1) + 
      ggtitle('Kori bustard') +
      scale_fill_viridis(name = "P(occ)", na.value="transparent") +
      scale_color_manual("Occurrence", values = c("orangered", "black")) +
      scale_shape_manual("Data source", values = c(15, 16), labels = c("Atlas", "eBird")) +
      scale_alpha_manual(values = c(0.3, 1), guide = "none") +
      xlab(element_blank()) +
      ylab(element_blank()) +
      theme_bw() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 18),
            strip.text.x = element_text(size = 18),
            text = element_text(size = 16))

load('model_output/Tockus_deckeni_figure_files.RData')
preds <- rbind(preds_1, preds_2) %>% 
      mutate(ind = ifelse(ind == 1, "1980-1999", "2000-2020"))
median_plot_hornbill <- ggplot() +
      geom_tile(data = preds, aes(x = x, y = y, fill = all_pred)) +
      geom_polygon(data = TZ_no_lakes, aes(x=long, y=lat, group=group), 
                   fill = NA, color = "black") +
      geom_point(data = atlas_sp_simple,
                 aes(x = x, y = y, col = as.factor(presence), alpha = as.factor(presence), shape = as.factor(presence)), pch = 15, cex = 3) +
      geom_point(data = ebird_sp_simple %>% filter(presence == 0),
                 aes(x = x, y = y, shape = as.factor(presence)), col = "orangered", alpha = 0.3) +
      geom_point(data = ebird_sp_simple %>% filter(presence == 1),
                 aes(x = x, y = y, shape = as.factor(presence)), col = "black") +
      coord_sf() +
      facet_wrap(.~ ind, ncol = 1) + 
      ggtitle("Von der Decken's hornbill") +
      scale_fill_viridis(name = "P(occorrence)", na.value="transparent") +
      scale_color_manual("Occurrence", values = c("orangered", "black")) +
      scale_shape_manual("Data source", values = c(15, 16), labels = c("Atlas", "eBird")) +
      scale_alpha_manual(values = c(0.3, 1), guide = "none") +
      xlab(element_blank()) +
      ylab(element_blank()) +
      theme_bw() +
      theme(legend.position = "right",
            plot.title = element_text(hjust = 0.5, size = 18),
            strip.text.x = element_text(size = 18),
            text = element_text(size = 16),
            legend.text = element_text(size = 14)) + 
      guides(shape = guide_legend(override.aes = list(size = 3), nrow = 2),
             color = guide_legend(override.aes = list(size = 3), nrow = 2))


combined_range_plots <- (cells_change_log / cor_plot) | median_plot_kori | median_plot_hornbill + plot_annotation(tag_levels = 'A') 

ggsave(filename = "figures/final_figures/combined_range_plots_R.png", plot = combined_range_plots, width = 38, height = 20,
       units = "cm", dpi = 300)


# -----------------------------------------------------------------------------------------------------------------
## Model output - forest plot ---------
# -----------------------------------------------------------------------------------------------------------------
prop_change_forest <- make_forest_plot(model_prop_change, "Sensitivity", "#0072B2") + 
      ggtitle("Total range change")

loss_forest <- make_forest_plot(model_loss_log, "Sensitivity", "#0072B2") + 
      ggtitle("Extinction") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank())

make_forest_plot(model_col_log, "Sensitivity", "#0072B2")  + 
      ggtitle("Colonisation") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank()) 

prop_change_forest_niche <- make_forest_plot(model_prop_change_breadth, "Niche breadth", "#E69F00") + 
      ggtitle("Total range change") +
      xlab(element_blank()) 

loss_forest_niche <- make_forest_plot(model_loss_log_breadth, "Niche breadth", "#E69F00") + 
      ggtitle("Extinction") +
      xlab(element_blank()) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank())

col_forest_niche <- make_forest_plot(model_col_log_breadth, "Niche breadth", "#E69F00") + 
      ggtitle("Colonisation") +
      xlab(element_blank()) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank())

forest_patch <- prop_change_forest_niche + loss_forest_niche + col_forest_niche + 
      prop_change_forest + loss_forest + col_forest 
forest_annotated <- forest_patch + plot_annotation(tag_levels = 'A')

ggsave(filename = "figures/final_figures/forest_plots_comb.png", plot = forest_annotated, width = 22, height = 20,
       units = "cm", dpi = 400)

# -----------------------------------------------------------------------------------------------------------------
## Selected plots showing raw data with model results ---------
# -----------------------------------------------------------------------------------------------------------------
rain_breadth_ext <- make_raw_data_plots(model_loss_log_breadth, figure_data, "rain_breadth", newscale.y = relative_scale, custom_unit_step = 100)
mass_sens_ext <- make_raw_data_plots(model_loss_log, figure_data, "Mass", newscale.y = relative_scale, custom_unit_step = 10)
dry_sens_col <-  make_raw_data_plots(model_col_log, figure_data, "dry_imp", newscale.y = relative_scale)
hot_sens_ext <-  make_raw_data_plots(model_loss_log, figure_data, "temp_imp", newscale.y = relative_scale)

raw_data_plots <- (rain_breadth_ext[[1]] + ylab("Meaningful extinction") | mass_sens_ext[[1]] + ylab("Meaningful extinction")) / (dry_sens_col[[1]] + ylab("Meaningful colonisation") | hot_sens_ext[[1]] + ylab("Meaningful extinction")) + plot_annotation(tag_levels = 'A') 

ggsave(filename = "figures/final_figures/raw_data_plots.png", plot = raw_data_plots, width = 25, height = 15,
       units = "cm", dpi = 400)

# -----------------------------------------------------------------------------------------------------------------
## Supplementary table of species data ------------
# -----------------------------------------------------------------------------------------------------------------
species_table <- figure_data %>% 
      dplyr::select(species, log_prop_change, relative_extinction, relative_colonisation, everything()) %>% 
      mutate(Mass = exp(Mass)) %>% 
      dplyr::select(-area_chance_transitions, -sum_dist_80s, -sum_dist_20s, -log_prop_change, -log_relative_extinction, -log_relative_colonisation,
                    -cells_lost, -cells_colonised)

options(ztable.type="viewer")
z = ztable(species_table) %>% makeHeatmap(margin = 2)
z

# -----------------------------------------------------------------------------------------------------------------
## Supplementary figures: All raw data plots --------------
# -----------------------------------------------------------------------------------------------------------------
list_change_sens <- make_raw_data_plots(model_prop_change, figure_data, newscale.y = relative_scale)
plot_change_sens <- wrap_plots(list_change_sens, ncol = 3) + plot_annotation(tag_levels = 'A', title = "Total range change - Sensitivity") 
png(filename = "figures/final_figures/plot_change_sens.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_change_sens
dev.off()

list_loss_sens <- make_raw_data_plots(model_loss_log, figure_data, newscale = relative_scale)
plot_loss_sens <- wrap_plots(list_loss_sens) + plot_annotation(tag_levels = 'A', title = "Meaningful extinction - Sensitivity") 
png(filename = "figures/final_figures/plot_loss_sens.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_loss_sens
dev.off()

list_col_sens <- make_raw_data_plots(model_col_log, figure_data, newscale = relative_scale)
plot_col_sens <- wrap_plots(list_col_sens) + plot_annotation(tag_levels = 'A', title = "Meaningful colonisation - Sensitivity") 
png(filename = "figures/final_figures/plot_col_sens.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_col_sens
dev.off()

list_change_spec <- make_raw_data_plots(model_prop_change_breadth, figure_data, newscale = relative_scale)
plot_change_spec <- wrap_plots(list_change_spec) + plot_annotation(tag_levels = 'A', title = "Total range change - Niche breadth") 
png(filename = "figures/final_figures/plot_change_spec.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_change_spec
dev.off()

list_loss_spec <- make_raw_data_plots(model_loss_log_breadth, figure_data, newscale = relative_scale)
plot_loss_spec <- wrap_plots(list_loss_spec) + plot_annotation(tag_levels = 'A', title = "Meaningful extinction - Niche breadth") 
png(filename = "figures/final_figures/plot_loss_spec.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_loss_spec
dev.off()

list_col_spec <- make_raw_data_plots(model_col_log_breadth, figure_data, newscale = relative_scale)
plot_col_spec <- wrap_plots(list_col_spec) + plot_annotation(tag_levels = 'A', title = "Meaningful colonisation - Niche breadth") 
png(filename = "figures/final_figures/plot_col_spec.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_col_spec
dev.off()
