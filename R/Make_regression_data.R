rm(list = ls())
library(tidyverse)
library(readxl)
library(INLA)

# -----------------------------------------------------------------------------------------------------------------
# Data import & preparation
# -----------------------------------------------------------------------------------------------------------------
source("source/misc_functions.R")
species_results <- read.csv('model_output/species_results.csv') %>% dplyr::select(-X)
traits <- readxl::read_xlsx('model_data/ELEData/TraitData/AVONET2_eBird.xlsx', sheet = "AVONET2_eBird")
colours <-  read.csv('model_data/plumage_lightness.csv') %>%
      mutate(species = gsub(pattern = "_", replacement = " ", species)) %>% 
      filter(sex == "F")
# Higher reflectance = lighter colors = like hotter areas, lower heat load 
# Darker birds are limited from occupying areas with high temperatures

species_results$species[species_results$species == "Estrilda erythronotos"] <- "Brunhilda erythronotos"
species_results$species[species_results$species == "Riparia cincta"] <- "Neophedina cincta"
species_results$species[species_results$species == "Sylvia boehmi"] <- "Curruca boehmi"
species_results$species[species_results$species == "Sylvia communis"] <- "Curruca communis"
species_results$species[species_results$species == "Turdoides aylmeri"] <- "Argya aylmeri"

species_results$temp_breadth <- species_results$temp_max - species_results$temp_min
species_results$rain_breadth <- species_results$rain_max - species_results$rain_min
species_results$dry_breadth <- species_results$dry_max - species_results$dry_min
species_results$HFP_breadth <- species_results$HFP_max - species_results$HFP_min
species_results$BG_breadth <- species_results$BG_max - species_results$BG_min

# Add trait data from Tobias et al. 2022
df_merged <- merge(species_results, traits, all.x = TRUE, all.y = FALSE, by.x = "species", by.y = "Species2") %>% 
      mutate(Trophic.Level = factor(Trophic.Level, levels = c("Omnivore", "Carnivore", "Herbivore")), 
             Primary.Lifestyle = as.factor(Primary.Lifestyle))
df_merged <- merge(df_merged, colours, all.x = TRUE,by = "species")

df_merged$Migratory_ability[df_merged$Migration == 1] <- "low"
df_merged$Migratory_ability[df_merged$Migration == 2] <- "moderate"
df_merged$Migratory_ability[df_merged$Migration == 3] <- "high"
df_merged <- df_merged %>% 
      mutate(Migratory_ability = as.factor(Migratory_ability))

all_species <- df_merged %>% 
      mutate(BG_imp = scale(BG_imp),
             rain_imp = scale(rain_imp),
             temp_imp = scale(temp_imp),
             dry_imp = scale(dry_imp),
             HFP_imp = scale(HFP_imp),
             Mass = scale(log(Mass)),   # We don't expect linear relationship with mass
             HWI = scale(`Hand-Wing.Index`),
             avg.r = scale(avg.r),
             rain_breadth = scale(rain_breadth),
             temp_breadth = scale(temp_breadth),
             dry_breadth = scale(dry_breadth),
             HFP_breadth = scale(HFP_breadth),
             BG_breadth = scale(BG_breadth),
             log_relative_colonisation = log((cells_colonised)/area_chance_transitions),
             log_relative_extinction = log((cells_lost)/area_chance_transitions),
             relative_colonisation = cells_colonised/area_chance_transitions,
             relative_extinction = cells_lost/area_chance_transitions,
             log_prop_change = log(sum_dist_20s/sum_dist_80s),
             prop_change = sum_dist_20s/sum_dist_80s,
             auc_mean = (auc_80s + auc_20s)/2) %>% 
      mutate(Migratory_ability = factor(Migratory_ability, levels = c("low", "moderate", "high")),
             Primary.Lifestyle = factor(Primary.Lifestyle, levels = c("Insessorial", "Terrestrial", "Aerial", "Generalist")),
             Trophic.Level = factor(Trophic.Level, levels = c("Carnivore", "Herbivore", "Omnivore"))) %>% 
      dplyr::select(species, relative_colonisation, relative_extinction, area_chance_transitions, sum_dist_80s, log_prop_change, 
                    log_relative_colonisation, log_relative_extinction, prop_change, cells_lost, cells_colonised, sum_dist_20s,
                    Trophic.Level, Primary.Lifestyle, Migratory_ability,
                    BG_imp, rain_imp, temp_imp,
                    dry_imp, HFP_imp, HWI, avg.r,
                    rain_breadth, temp_breadth, dry_breadth,
                    HFP_breadth, BG_breadth, auc_80s, auc_20s, auc_mean, Mass)

model_data <-  all_species %>% 
      filter(species != "Gyps africanus",        # only scavenger
             species != "Struthio camelus")  %>%     # Mass outlier
      drop_na()

# -----------------------------------------------------------------------------------------------------------------
# Linear combinations for INLA model
# -----------------------------------------------------------------------------------------------------------------
BG_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                            BG_imp = seq(min(model_data$BG_imp), max(model_data$BG_imp), len = 100))
names(BG_lc) <- paste0("BG_imp", 1:100)
rain_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              rain_imp = seq(min(model_data$rain_imp), max(model_data$rain_imp), len = 100))
names(rain_lc) <- paste0("rain_imp", 1:100)
temp_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              temp_imp = seq(min(model_data$temp_imp), max(model_data$temp_imp), len = 100))
names(temp_lc) <- paste0("temp_imp", 1:100)
dry_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             dry_imp = seq(min(model_data$dry_imp), max(model_data$dry_imp), len = 100))
names(dry_lc) <- paste0("dry_imp", 1:100)
HFP_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             HFP_imp = seq(min(model_data$HFP_imp), max(model_data$HFP_imp), len = 100))
names(HFP_lc) <- paste0("HFP_imp", 1:100)
Mass_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              Mass = seq(min(model_data$Mass), max(model_data$Mass), len = 100))
names(Mass_lc) <- paste0("Mass", 1:100)
HWI_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              HWI = seq(min(model_data$HWI), max(model_data$HWI), len = 100))
names(HWI_lc) <- paste0("HWI", 1:100)
reflect_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                 avg.r = seq(min(model_data$avg.r, na.rm = TRUE), max(model_data$avg.r, na.rm = TRUE), len = 100))
names(reflect_lc) <- paste0("avg.r", 1:100)
Trophic_level_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
                                       Trophic.LevelHerbivore = c(0, 1, 0),
                                       Trophic.LevelOmnivore = c(0, 0, 1))
names(Trophic_level_lc) <- paste0("Trophic.Level", 1:3)
Migratory_ability_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
                                           Migratory_abilitymoderate = c(0, 1, 0),
                                           Migratory_abilityhigh = c(0, 0, 1))
names(Migratory_ability_lc) <- paste0("Migratory_ability", 1:3)
Primary.Lifestyle_lc <- inla.make.lincombs("(Intercept)" = rep(1, 4),
                                           Primary.LifestyleTerrestrial = c(0, 1, 0, 0),
                                           Primary.LifestyleAerial = c(0, 0, 1, 0),
                                           Primary.LifestyleGeneralist = c(0, 0, 0, 1))
names(Primary.Lifestyle_lc) <- paste0("Primary.Lifestyle", 1:4)

all_lc_sens <- c(BG_lc, rain_lc, temp_lc, dry_lc, HFP_lc, HWI_lc, Mass_lc, reflect_lc, 
                 Trophic_level_lc, Migratory_ability_lc, Primary.Lifestyle_lc)


BG_lc_breadth <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                  BG_breadth = seq(min(model_data$BG_breadth), max(model_data$BG_breadth), len = 100))
names(BG_lc_breadth) <- paste0("BG_breadth", 1:100)
rain_lc_breadth <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                    rain_breadth = seq(min(model_data$rain_breadth), max(model_data$rain_breadth), len = 100))
names(rain_lc_breadth) <- paste0("rain_breadth", 1:100)
temp_lc_breadth <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                    temp_breadth = seq(min(model_data$temp_breadth), max(model_data$temp_breadth), len = 100))
names(temp_lc_breadth) <- paste0("temp_breadth", 1:100)
dry_lc_breadth <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                   dry_breadth = seq(min(model_data$dry_breadth), max(model_data$dry_breadth), len = 100))
names(dry_lc_breadth) <- paste0("dry_breadth", 1:100)
HFP_lc_breadth <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                   HFP_breadth = seq(min(model_data$HFP_breadth), max(model_data$HFP_breadth), len = 100))
names(HFP_lc_breadth) <- paste0("HFP_breadth", 1:100)


all_lc_breadth <- c(BG_lc_breadth, rain_lc_breadth, temp_lc_breadth, dry_lc_breadth, HFP_lc_breadth, HWI_lc, Mass_lc, reflect_lc, 
                  Trophic_level_lc, Migratory_ability_lc, Primary.Lifestyle_lc)

save(model_data, all_species, all_lc_sens, all_lc_breadth, file = 'model_data/regression_data.RData')

