suppressPackageStartupMessages(library(INLA, quietly=TRUE))
library(tidyverse)
library(raster)
library(pROC)

modelData_list <- list.files(path = "/users/jhw538/scratch/TZ_spatio_temporal_INLA/model_output_final/good/", pattern = ".RData",
                             full.names = TRUE, recursive = TRUE)

species_results <- data.frame(matrix(NA, nrow = length(modelData_list), ncol = 24))
colnames(species_results) <- c('species', 'rsq_random',
                               'rain_imp', 'temp_imp', 'dry_imp', 'BG_imp', 'HFP_imp',
                               'sum_dist_80s', 'sum_dist_20s', 'cells_colonised', 'cells_lost', 'area_chance_transitions', 
                               'rain_min', 'rain_max', 'temp_min', 'temp_max', 'BG_min', 'BG_max',
                               'dry_min', 'dry_max', 'HFP_min', 'HFP_max', 'auc_80s', 'auc_20s')

add_quantiles <- function(df, column, new_var){
      med_value <- median(df[, column], na.rm = TRUE)
      df[, new_var] <- NA
      df[, new_var][df[, column] <= 0] <- "lower"
      df[, new_var][df[, column] > 0 & df[, column] <= med_value] <- "middle"
      df[, new_var][df[, column] > med_value] <- "upper"
      return(df)
}

for (j in 1:length(modelData_list)) {
      load(modelData_list[j])
      print(paste0("Starting file: ", modelData_list[j]))
      species_results[j, 'species'] <- species
      
      pred.index <- inla.stack.index(stack = integrated_stack, tag = "pred")$data
      lincomb.index.temp <- grep('TZ_lc_noTempMax', rownames(model$summary.lincomb.derived))
      lincomb.index.rain <- grep('TZ_lc_noRain', rownames(model$summary.lincomb.derived))
      lincomb.index.dry <- grep('TZ_lc_noDryspell', rownames(model$summary.lincomb.derived))
      lincomb.index.BG <- grep('TZ_lc_noBG', rownames(model$summary.lincomb.derived))
      lincomb.index.HFP <- grep('TZ_lc_noHFP', rownames(model$summary.lincomb.derived))
      lincomb.index.allFixed <- grep('TZ_lc_allFixed', rownames(model$summary.lincomb.derived))
      
      IP_df <- data.frame(IP_sp) %>% dplyr::select(date_index, x, y)
      
      IP_df$pred_median <- model$summary.fitted.values[pred.index, "0.5quant"]
      IP_df$pred_sd <- model$summary.fitted.values[pred.index, "sd"]
      IP_df$no_temp_median <- model$summary.lincomb.derived[lincomb.index.temp, "0.5quant"]
      IP_df$no_rain_median <- model$summary.lincomb.derived[lincomb.index.rain, "0.5quant"]
      IP_df$no_dry_median <- model$summary.lincomb.derived[lincomb.index.dry, "0.5quant"]
      IP_df$no_BG_median <- model$summary.lincomb.derived[lincomb.index.BG, "0.5quant"]
      IP_df$no_HFP_median <- model$summary.lincomb.derived[lincomb.index.HFP, "0.5quant"]
      IP_df$allFixed_median <- model$summary.lincomb.derived[lincomb.index.allFixed, "0.5quant"]
      
      pred_data <- data.frame(all_pred = IP_df$pred_median,
                              no_temp_median = IP_df$no_temp_median,  
                              no_rain_median = IP_df$no_rain_median,
                              no_dry_median = IP_df$no_dry_median,
                              no_BG_median = IP_df$no_BG_median,
                              no_HFP_median = IP_df$no_HFP_median,
                              allFixed_median = IP_df$allFixed_median,
                              sd = IP_df$pred_sd,
                              ind = rep(c(1,2), each = length(IP_df$pred_median)/2))
      
      predcoordsGroup <- do.call(rbind, list(predcoords, predcoords))
      pred_data_coords <- cbind(pred_data, predcoordsGroup) %>% 
            dplyr::select(x, y, everything())
      
      preds_1 <- pred_data_coords %>% 
            filter(ind == 1) %>% 
            rasterFromXYZ(.) %>% 
            mask(., TZ_no_lakes) %>% 
            as.data.frame(., xy = T) %>% 
            mutate(ind = "1980-1999") %>% 
            filter(!if_all(-ind, is.na))
      
      preds_2 <- pred_data_coords %>% 
            filter(ind == 2) %>% 
            rasterFromXYZ(.) %>% 
            mask(., TZ_no_lakes) %>% 
            as.data.frame(., xy = T) %>% 
            mutate(ind = "2000-2020") %>% 
            filter(!if_all(-ind, is.na))
      
      all_pred_df <- rbind(preds_1, preds_2)
      print("median done")
      
      xmedian <- list()
      for (k in 1:2) {
            xmedian[[k]] <- inla.mesh.project(
                  projgrid,  model$summary.random$shared.field$`0.5quant`[index_set$shared.field.group == k])
      }
      
      xy.inGroup <- c(xy.in,xy.in)
      xmedian <- unlist(xmedian)
      xmedian <- xmedian[xy.inGroup]
      
      dataObj <- data.frame(random_shared = xmedian,
                            ind = rep(c(1,2), each = length(xmedian)/2))
      
      spatObj_coords <- cbind(dataObj, predcoordsGroup) %>% 
            dplyr::select(x, y, everything())
      
      spatObj_1 <- spatObj_coords %>% 
            filter(ind == 1) %>% 
            rasterFromXYZ(.) %>% 
            mask(., TZ_no_lakes) %>% 
            as.data.frame(., xy = T) %>% 
            mutate(ind = "1980-1999") %>% 
            filter(!if_all(-ind, is.na))
      
      spatObj_2 <- spatObj_coords %>% 
            filter(ind == 2) %>% 
            rasterFromXYZ(.) %>% 
            mask(., TZ_no_lakes) %>% 
            as.data.frame(., xy = T) %>% 
            mutate(ind = "2000-2020") %>% 
            filter(!if_all(-ind, is.na))
      
      random_shared_df <- rbind(spatObj_1, spatObj_2)
      
      xmedian2 <- list()
      for (k in 1:2) {
            xmedian2[[k]] <- inla.mesh.project(
                  projgrid,  model$summary.random$eBird.field$`0.5quant`[index_set_eBird$eBird.field.group == k])
      }
      
      xmedian2 <- unlist(xmedian2)
      xmedian2 <- xmedian2[xy.inGroup]
      
      dataObj2 <- data.frame(random_eBird = xmedian2,
                             ind = rep(c(1,2), each = length(xmedian2)/2))
      spatObj2_coords <- cbind(dataObj2, predcoordsGroup) %>% 
            dplyr::select(x, y, everything())
      
      spatObj2_1 <- spatObj2_coords %>% 
            filter(ind == 1) %>% 
            rasterFromXYZ(.) %>% 
            mask(., TZ_no_lakes) %>% 
            as.data.frame(., xy = T) %>% 
            mutate(ind = "1980-1999") %>% 
            filter(!if_all(-ind, is.na))
      
      spatObj2_2 <- spatObj2_coords %>% 
            filter(ind == 2) %>% 
            rasterFromXYZ(.) %>% 
            mask(., TZ_no_lakes) %>% 
            as.data.frame(., xy = T) %>% 
            mutate(ind = "2000-2020") %>% 
            filter(!if_all(-ind, is.na))
      
      random_eBird_df <- rbind(spatObj2_1, spatObj2_2)
      
      print("random done")
      all_merged_1 <- merge(all_pred_df, random_eBird_df, by = c("ind", "x", "y"))
      all_merged <- merge(all_merged_1, random_shared_df, by = c("ind", "x", "y"))
      all_merged$all_preds_est <- all_merged$random_eBird + all_merged$random_shared + all_merged$allFixed_median
      all_merged$all_preds_no_rain <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_rain_median
      all_merged$all_preds_no_temp <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_temp_median
      all_merged$all_preds_no_dry <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_dry_median
      all_merged$all_preds_no_BG <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_BG_median
      all_merged$all_preds_no_HFP <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_HFP_median
      all_merged$all_random <- all_merged$random_eBird + all_merged$random_shared
      
      all_merged[, 4:20] <- lapply(all_merged[, 4:20], inla.link.cloglog, inverse = TRUE) 
      
      rsq <- function(x, y) summary(lm(y~x))$r.squared
      
      species_results[j, 'rsq_random'] <- rsq(all_merged$all_preds_est, all_merged$all_random)
      
      species_results[j, 'rain_imp'] <- 1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_rain)
      species_results[j, 'temp_imp'] <- 1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_temp)
      species_results[j, 'dry_imp'] <- 1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_dry)
      species_results[j, 'BG_imp'] <- 1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_BG)
      species_results[j, 'HFP_imp'] <- 1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_HFP)
      
      original_values <- data.frame(orig_values = unlist(all.seq)) %>% 
            mutate(covariate = sub("*.seq\\d+", "", rownames(.)),
                   sequence = as.numeric(gsub("\\D", "", rownames(.))))
      
      scale_params <- median(model$summary.lincomb.derived$`0.5quant`[grep("TZ_ann_rain|TZ_BG|TZ_dryspell|TZ_max_temp|TZ_HFP", 
                                                                           rownames(model$summary.lincomb.derived))])
      effect_combs <- data.frame(covariate = gsub('[[:digit:]]+', '', sub("*_lc\\d+", "", rownames(model$summary.lincomb.derived))),
                                 sequence = as.numeric(gsub("\\D", "", rownames(model$summary.lincomb.derived))),
                                 quant_05 = inla.link.cloglog((model$summary.lincomb.derived$`0.5quant` - scale_params) - 0.36, inverse = TRUE),
                                 quant_0025 = inla.link.cloglog((model$summary.lincomb.derived$`0.025quant` - scale_params) - 0.36, inverse = TRUE),
                                 quant_0975 = inla.link.cloglog((model$summary.lincomb.derived$`0.975quant` - scale_params) - 0.36, inverse = TRUE))
      effect_combs_main <- effect_combs %>% filter(covariate %in% c("TZ_ann_rain", "TZ_BG", "TZ_dryspell", "TZ_max_temp", "TZ_HFP"))
      effect_combs_m <- merge(original_values, effect_combs_main)
      
      df <- data.frame(covariate = unique(effect_combs_m$covariate))
      
      min_max <- effect_combs_m %>% 
            group_by(covariate) %>% 
            filter(quant_05 > 0.5) %>% 
            summarize(min = min(orig_values),
                      max = max(orig_values))
      df_merged <- merge(df, min_max, all.x = TRUE)
      
      species_results[j, 'rain_min'] <- df_merged$min[df_merged$covariate == "TZ_ann_rain"]
      species_results[j, 'rain_max'] <- df_merged$max[df_merged$covariate == "TZ_ann_rain"]
      species_results[j, 'temp_min'] <- df_merged$min[df_merged$covariate == "TZ_max_temp"]
      species_results[j, 'temp_max'] <- df_merged$max[df_merged$covariate == "TZ_max_temp"]
      species_results[j, 'BG_min'] <- df_merged$min[df_merged$covariate == "TZ_BG"]
      species_results[j, 'BG_max'] <- df_merged$max[df_merged$covariate == "TZ_BG"]
      species_results[j, 'dry_min'] <- df_merged$min[df_merged$covariate == "TZ_dryspell"]
      species_results[j, 'dry_max'] <- df_merged$max[df_merged$covariate == "TZ_dryspell"]
      species_results[j, 'HFP_min'] <- df_merged$min[df_merged$covariate == "TZ_HFP"]
      species_results[j, 'HFP_max'] <- df_merged$max[df_merged$covariate == "TZ_HFP"]
      
      dist_80s <- pred_data$all_pred[pred_data$ind == "1980-1999"]  
      dist_20s <- pred_data$all_pred[pred_data$ind == "2000-2020"]  
      
      dist_80s_P <- inla.link.cloglog(dist_80s, inverse = TRUE)
      dist_80s_P <- dist_80s_P[!is.na(dist_80s_P)]
      dist_20s_P <- inla.link.cloglog(dist_20s, inverse = TRUE)
      dist_20s_P <- dist_20s_P[!is.na(dist_20s_P)]
      
      # Number of cells colonised and lost
      cells_colonised <- (1-dist_80s_P) * dist_20s_P
      cells_lost <- dist_80s_P * (1-dist_20s_P)
      
      species_results[j, 'cells_colonised'] <- sum(cells_colonised)
      species_results[j, 'cells_lost'] <- sum(cells_lost)
      
      species_results[j, 'sum_dist_80s'] <- sum(dist_80s_P)
      species_results[j, 'sum_dist_20s'] <- sum(dist_20s_P)
      
      species_results[j, 'area_chance_transitions'] <- sum(dist_80s_P * (1-dist_80s_P))
      
      atlas_sp_simple <- as.data.frame(atlas_sp)  %>% 
            dplyr::select("presence", "date_index", "x", 'y') %>% 
            mutate(ind = ifelse(date_index == 1, '1980-1999', '2000-2020'))
      
      ebird_sp_simple <- as.data.frame(ebird_sp)  %>% 
            dplyr::select("presence", "date_index", "x", 'y') %>% 
            mutate(ind = ifelse(date_index == 1, '1980-1999', '2000-2020'))
      
      print(paste0(species, " eBird presence 1980-1999: ", sum(ebird_sp_simple$presence[ebird_sp_simple$ind == '1980-1999'], na.rm = TRUE)))
      print(paste0(species, " eBird presence 2000-2020: ", sum(ebird_sp_simple$presence[ebird_sp_simple$ind == '2000-2020'], na.rm = TRUE)))
      print(paste0(species, " Atlas presence 1980-1999: ", sum(atlas_sp_simple$presence[atlas_sp_simple$ind == '1980-1999'], na.rm = TRUE)))
      print(paste0(species, " Atlas presence 2000-2020: ", sum(atlas_sp_simple$presence[atlas_sp_simple$ind == '2000-2020'], na.rm = TRUE)))
      
      all_obs <- rbind(atlas_sp_simple, ebird_sp_simple)
      all_obs_80s <- all_obs %>% filter(ind == '1980-1999') 
      all_obs_20s <- all_obs %>% filter(ind == '2000-2020')
      
      all_obs_data <- all_pred_df %>% 
            dplyr::select(x, y, ind, all_pred, sd) %>% 
            mutate(pred_median = inla.link.cloglog(all_pred, inv = TRUE),
                   pred_sd = sd)
      
      pred_data_80s <- all_obs_data %>% filter(ind == '1980-1999') 
      pred_data_20s <- all_obs_data %>% filter(ind == '2000-2020')
      
      pred_raster_80s <- rasterFromXYZ(pred_data_80s[, c("x", "y", "pred_median")])
      pred_raster_20s <- rasterFromXYZ(pred_data_20s[, c("x", "y", "pred_median")])
      
      all_obs_80s_spdf <- SpatialPointsDataFrame(coords = data.frame(all_obs_80s$x, all_obs_80s$y), data = data.frame(presence = all_obs_80s$presence))
      all_obs_20s_spdf <- SpatialPointsDataFrame(coords = data.frame(all_obs_20s$x, all_obs_20s$y), data = data.frame(presence = all_obs_20s$presence))
      
      ext_80s <- extract(pred_raster_80s, all_obs_80s_spdf)
      ext_20s <- extract(pred_raster_20s, all_obs_20s_spdf)
      
      all_obs_80s$pred <- ext_80s
      all_obs_20s$pred <- ext_20s
      
      all_obs_80s <- all_obs_80s %>% 
            filter(!is.na(pred))
      all_obs_20s <- all_obs_20s %>% 
            filter(!is.na(pred))
      
      auc_80s <- round(as.numeric(auc(all_obs_80s$presence, all_obs_80s$pred)), digits = 3)
      auc_20s <- round(as.numeric(auc(all_obs_20s$presence, all_obs_20s$pred)), digits = 3)
      
      species_results[j, 'auc_80s'] <- auc_80s
      species_results[j, 'auc_20s'] <- auc_20s
      
      print(paste0("Finished ", j, " of ", length(modelData_list), ": ", species))
      
}

write.csv(species_results, '/users/jhw538/scratch/TZ_spatio_temporal_INLA/model_output_final/species_results.csv')
