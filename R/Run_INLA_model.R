rm(list = ls())
args = commandArgs(trailingOnly = TRUE)
options(scipen=999)

library(tidyverse)
library(rgdal)
library(raster)
library(rgeos)
library(ggpubr)
suppressPackageStartupMessages(library(INLA, quietly=TRUE))
library(ggthemes)
library(terra)
library(tidyterra)
library(viridis)
library(patchwork)

sapply(list.files(path = "source/", pattern="*.R", full.names = T), source, .GlobalEnv)

load("model_data/TZ_INLA_model_file_temporal_V2.RData")

# ------------------------------------------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------------------------------------------
species = c('Eremopterix leucopareia')  

max.edge = 0.35
multiplier <- 0.5

# Gaussian priors for fixed effects. Change from default due to cloglog link
fixed_mean = 0
fixed_precision = 1

# ------------------------------------------------------------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------------------------------------------------------------
ebird_full <- ebird_full %>%
      mutate(date_index = ifelse(date > '2000-01-01',2,1))

ebird_filtered <- ebird_full %>% 
      filter(APPROVED == 1,  
             `ALL SPECIES REPORTED` == 1,         
             `EFFORT DISTANCE KM` < 15,
             `DURATION MINUTES` >= 5,
             `DURATION MINUTES` <= 5*60,
             `NUMBER OBSERVERS` <= 10)   

ebird_filtered <- ebird_filtered %>% 
      group_by(LATITUDE, LONGITUDE, `SAMPLING EVENT IDENTIFIER`, `DURATION MINUTES`, 
               `EFFORT DISTANCE KM`, `NUMBER OBSERVERS`, `OBSERVATION DATE`, `LOCALITY`, date_index) %>% 
      summarise(occurrence = ifelse(species %in% `SCIENTIFIC NAME`, TRUE, FALSE)) %>% 
      ungroup()  %>%
      group_by(LATITUDE, LONGITUDE, `DURATION MINUTES`, 
               `EFFORT DISTANCE KM`, `NUMBER OBSERVERS`, `OBSERVATION DATE`, `LOCALITY`) %>% 
      slice_head() %>%    
      ungroup() %>% 
      rename(duration_minutes = `DURATION MINUTES`,
             effort_distance_km = `EFFORT DISTANCE KM`,
             number_observers = `NUMBER OBSERVERS`,
             x = LONGITUDE,
             y = LATITUDE) %>% 
      mutate(ebird_effort = number_observers*duration_minutes)  

ebird_sp <- SpatialPointsDataFrame(
      coords = ebird_filtered[, c("x", "y")],
      data = data.frame(presence = ebird_filtered$occurrence, 
                        ebird_effort = ebird_filtered$ebird_effort,
                        duration_minutes = ebird_filtered$duration_minutes, 
                        number_observers = ebird_filtered$number_observers, 
                        date_index = ebird_filtered$date_index),
      proj4string = crs(proj))

# Only include observation data points for Tanzania
in_sp <- rgeos::gIntersection(ebird_sp, TZ_outline)
ebird_sp <- ebird_sp[in_sp, ]

atlas_full <- atlas_full %>%
      filter(time_period != 'x') %>%
      mutate(date_index = ifelse(time_period == '20s',2,1))

atlas_filtered <- atlas_full %>% 
      mutate(Scientific = trimws(Scientific, which = 'both')) %>% 
      filter(Scientific == species) %>% 
      rename(x = Long,
             y = Lat) %>% 
      mutate(presence = ifelse(occurrence == 1, TRUE, FALSE)) %>% 
      dplyr::select(-V1); if(is_empty(atlas_filtered$presence)){print("ERROR: No Atlas data available")}

atlas_sp <- SpatialPointsDataFrame(
      coords = data.frame(atlas_filtered[, c("x", "y")]),
      data = data.frame(presence = atlas_filtered$presence, effort = atlas_filtered$effort,
                        date_index = atlas_filtered$date_index),
      proj4string = crs(proj))

# Log transform effort, then range 0-1
range01 <- function(x){(x - min(x))/(max(x) - min(x))}
ebird_sp$duration_minutes <- range01(log(ebird_sp$duration_minutes+1))
atlas_sp$effort <- range01(log(atlas_sp$effort+1))

filtered_covs <- temporal_variables[, c('TZ_ann_rain_1980s_1.s', 'TZ_ann_rain_2000s_1.s',
                                        'TZ_ann_rain_1980s_2.s', 'TZ_ann_rain_2000s_2.s',
                                        'TZ_max_temp_1980s_1.s', 'TZ_max_temp_2000s_1.s',
                                        'TZ_max_temp_1980s_2.s', 'TZ_max_temp_2000s_2.s',
                                        'TZ_dryspell_1980s_1.s', 'TZ_dryspell_2000s_1.s',
                                        'TZ_dryspell_1980s_2.s', 'TZ_dryspell_2000s_2.s',
                                        'TZ_BG_1980s_1.s', 'TZ_BG_2000s_1.s',
                                        'TZ_BG_1980s_2.s', 'TZ_BG_2000s_2.s',
                                        'indicator_90s', 'indicator_2010s',
                                        'TZ_HFP_1980s_1.s', 'TZ_HFP_2000s_1.s',
                                        'TZ_HFP_1980s_2.s', 'TZ_HFP_2000s_2.s')]

# Sample the covariates spatially closest to the occurrence data points
Nearest_covs_ebird <- GetNearestCovariate(ebird_sp, filtered_covs)
Nearest_covs_atlas <- GetNearestCovariate(atlas_sp, filtered_covs)

# Add sampled covariates to occurrence data
ebird_sp@data[, names(Nearest_covs_ebird@data)] <- Nearest_covs_ebird@data
ebird_sp <- as(ebird_sp, 'data.frame')

# Combine covariates from different time periods into single variable, different times identified by separate date_index
ebird_sp <- ebird_sp %>% mutate(annual_rain_1 = ifelse(date_index == 1, TZ_ann_rain_1980s_1.s, TZ_ann_rain_2000s_1.s),
                                annual_rain_2 = ifelse(date_index == 1, TZ_ann_rain_1980s_2.s, TZ_ann_rain_2000s_2.s),
                                hottest_temp_1 = ifelse(date_index == 1, TZ_max_temp_1980s_1.s, TZ_max_temp_2000s_1.s),
                                hottest_temp_2 = ifelse(date_index == 1, TZ_max_temp_1980s_2.s, TZ_max_temp_2000s_2.s),
                                max_dryspell_1 = ifelse(date_index == 1, TZ_dryspell_1980s_1.s, TZ_dryspell_2000s_1.s),
                                max_dryspell_2 = ifelse(date_index == 1, TZ_dryspell_1980s_2.s, TZ_dryspell_2000s_2.s),
                                BG_1 = ifelse(date_index == 1, TZ_BG_1980s_1.s, TZ_BG_2000s_1.s),
                                BG_2 = ifelse(date_index == 1, TZ_BG_1980s_2.s, TZ_BG_2000s_2.s),
                                indicator = ifelse(date_index == 1, indicator_90s, indicator_2010s),
                                HFP_1 = ifelse(date_index == 1, TZ_HFP_1980s_1.s, TZ_HFP_2000s_1.s),
                                HFP_2 = ifelse(date_index == 1, TZ_HFP_1980s_2.s, TZ_HFP_2000s_2.s))

# Make spdf, add intercept, so model can distinguish eBird presence from Atlas presence
ebird_sp <- SpatialPointsDataFrame(coords = ebird_sp[, c("x", "y")],
                                   data = ebird_sp,
                                   proj4string = crs(proj))
ebird_sp@data[, 'ebird_intercept'] <- 1
ebird_sp$presence <- as.numeric(ebird_sp$presence)

atlas_sp@data[, names(Nearest_covs_atlas@data)] <- Nearest_covs_atlas@data
atlas_sp <- as(atlas_sp, 'data.frame')

atlas_sp <- atlas_sp %>% mutate(annual_rain_1 = ifelse(date_index == 1, TZ_ann_rain_1980s_1.s, TZ_ann_rain_2000s_1.s),
                                annual_rain_2 = ifelse(date_index == 1, TZ_ann_rain_1980s_2.s, TZ_ann_rain_2000s_2.s),
                                hottest_temp_1 = ifelse(date_index == 1, TZ_max_temp_1980s_1.s, TZ_max_temp_2000s_1.s),
                                hottest_temp_2 = ifelse(date_index == 1, TZ_max_temp_1980s_2.s, TZ_max_temp_2000s_2.s),
                                max_dryspell_1 = ifelse(date_index == 1, TZ_dryspell_1980s_1.s, TZ_dryspell_2000s_1.s),
                                max_dryspell_2 = ifelse(date_index == 1, TZ_dryspell_1980s_2.s, TZ_dryspell_2000s_2.s),
                                BG_1 = ifelse(date_index == 1, TZ_BG_1980s_1.s, TZ_BG_2000s_1.s),
                                BG_2 = ifelse(date_index == 1, TZ_BG_1980s_2.s, TZ_BG_2000s_2.s),
                                indicator = ifelse(date_index == 1, indicator_90s, indicator_2010s),
                                HFP_1 = ifelse(date_index == 1, TZ_HFP_1980s_1.s, TZ_HFP_2000s_1.s),
                                HFP_2 = ifelse(date_index == 1, TZ_HFP_1980s_2.s, TZ_HFP_2000s_2.s))

atlas_sp <- SpatialPointsDataFrame(coords = atlas_sp[, c('x', 'y')],
                                   data = atlas_sp,
                                   proj4string = crs(proj))
atlas_sp@data[,'atlas_intercept'] <- 1
atlas_sp$presence <- as.numeric(atlas_sp$presence)

print(species)
print(paste0('eBird presence records 1980s: ', length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 1, ]))))
print(paste0('eBird presence records 2000s: ', length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 2, ]))))
print(paste0('Atlas presence records 1980s: ', length(rownames(atlas_sp@data[atlas_sp@data$presence == TRUE & as.data.frame(atlas_sp)$date_index == 1, ]))))
print(paste0('Atlas presence records 2000s: ', length(rownames(atlas_sp@data[atlas_sp@data$presence == TRUE & as.data.frame(atlas_sp)$date_index == 2, ]))))

# Plot the raw data
plot_80s <- ggplot() +
      geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group),
                   colour = "white", fill = NA) +
      geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 0 & as.data.frame(atlas_sp)$date_index == 1, ],
                 aes(x = x, y = y, col = "Absence"), alpha = 0.8, pch = 15, cex = 2) +
      geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 1 & as.data.frame(atlas_sp)$date_index == 1, ],
                 aes(x = x, y = y, col = "Presence"), alpha = 0.8, pch = 15, cex = 3) +
      geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 0 & as.data.frame(ebird_sp)$date_index == 1, ],
                 aes(x = x, y = y, fill = "Absence"), pch = 21, alpha = 0.5) +
      geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 1 & as.data.frame(ebird_sp)$date_index == 1, ],
                 aes(x = x, y = y, fill = "Presence"), pch = 23, cex = 2) +
      theme_void() +
      coord_equal() +
      scale_fill_discrete(name = 'eBird') +
      scale_color_discrete(name = 'Atlas') +
      ggtitle(paste0("1980-1999, eBird presence: ", length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 1, ]))))

plot_2000s <- ggplot() +
      geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group),
                   colour = "white", fill = NA) +
      geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 0 & as.data.frame(atlas_sp)$date_index == 2, ],
                 aes(x = x, y = y, col = "Absence"), alpha = 0.8, pch = 15, cex = 2) +
      geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 1 & as.data.frame(atlas_sp)$date_index == 2, ],
                 aes(x = x, y = y, col = "Presence"), alpha = 0.8, pch = 15, cex = 3) +
      geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 0 & as.data.frame(ebird_sp)$date_index == 2, ],
                 aes(x = x, y = y, fill = "Absence"), pch = 21, alpha = 0.5) +
      geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 1 & as.data.frame(ebird_sp)$date_index == 2, ],
                 aes(x = x, y = y, fill = "Presence"), pch = 23, cex = 2) +
      theme_void() +
      coord_equal() +
      scale_fill_discrete(name = 'eBird') +
      scale_color_discrete(name = 'Atlas') +
      ggtitle(paste0("2000-2020, eBird presence: ", length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 2, ]))))

raw_plots_comb <- ggpubr::ggarrange(plot_80s, plot_2000s, nrow = 2, common.legend = TRUE, legend = "right")
raw_plots_comb

# ------------------------------------------------------------------------------------------------------------------------
# Mesh
# ------------------------------------------------------------------------------------------------------------------------
Meshpars <- list(max.edge = c(max.edge, max.edge*4), 
                 offset = c(max.edge, max.edge*5), 
                 cutoff = max.edge/2)

Mesh <- MakeSpatialRegion(
      data = NULL,
      bdry = ROI,    
      meshpars = Meshpars,
      proj = proj
)
Mesh$mesh$crs <- proj

# ------------------------------------------------------------------------------------------------------------------------
# SPDE model on the mesh
# ------------------------------------------------------------------------------------------------------------------------
# Penalized complexity priors for SPDE. Pick dynamically based on spatial range of observations 
eBird.prior.range = c(round(multiplier*diff(range(ebird_filtered$y[ebird_filtered$occurrence == 1])), digits = 2), 0.5)
eBird.prior.sigma = c(round(multiplier*diff(range(ebird_filtered$y[ebird_filtered$occurrence == 1]))/1.2, digits = 2), 0.5)

prior_range = c(round(multiplier*diff(range(atlas_filtered$y[atlas_filtered$presence == 1])), digits = 2), 0.5)
prior_sigma = c(round(multiplier*diff(range(atlas_filtered$y[atlas_filtered$presence == 1]))/1.5, digits = 2), 0.5)

pcspde <- inla.spde2.pcmatern(
      mesh = Mesh$mesh,
      alpha = 2,
      prior.range = prior_range,   
      prior.sigma = prior_sigma)

eBird_spde <- inla.spde2.pcmatern(
      mesh = Mesh$mesh,
      alpha = 2,
      prior.range = prior_range,
      prior.sigma = prior_sigma)

# ------------------------------------------------------------------------------------------------------------------------
# Space-time index set
# ------------------------------------------------------------------------------------------------------------------------
# Set the number of time periods
n_time_layers = 2

index_set <- inla.spde.make.index(name ='shared.field',
                                  n.spde = pcspde$n.spde,
                                  n.group = n_time_layers)

# Make separate random spatial field for eBird data
index_set_eBird <- inla.spde.make.index(name ='eBird.field',
                                        n.spde = eBird_spde$n.spde,
                                        n.group = n_time_layers)

# ------------------------------------------------------------------------------------------------------------------------
# Projection matrices
# ------------------------------------------------------------------------------------------------------------------------
projmat_eBird <- inla.spde.make.A(mesh = Mesh$mesh,
                                  loc = as.matrix(ebird_sp@coords),
                                  group = ebird_sp$date_index)

projmat_atlas <- inla.spde.make.A(mesh = Mesh$mesh,
                                  loc = as.matrix(atlas_sp@coords),
                                  group = atlas_sp$date_index)

# ------------------------------------------------------------------------------------------------------------------------
# Estimation stacks
# ------------------------------------------------------------------------------------------------------------------------
stk.eBird <- inla.stack(data = list(resp = ebird_sp@data[, 'presence']),
                        A = list(1, projmat_eBird, projmat_eBird),
                        tag = 'eBird',
                        effects = list(ebird_sp@data, eBird.field = index_set_eBird, shared.field = index_set))

stk.atlas <- inla.stack(data = list(resp = atlas_sp@data[, 'presence']),
                        A = list(1, projmat_atlas),
                        tag = 'atlas',
                        effects = list(atlas_sp@data, shared.field = index_set))

# ------------------------------------------------------------------------------------------------------------------------
# Prediction data
# ------------------------------------------------------------------------------------------------------------------------
# Set grid locations for predictions
Nxy.scale <- 0.1  # about 10km resolution

Boundary <- Mesh$mesh$loc[Mesh$mesh$segm$int$idx[, 2], ]
Nxy.size <- c(diff(range(Boundary[, 1])), diff(range(Boundary[, 2])))
Nxy <- round(Nxy.size / Nxy.scale)

projgrid <- inla.mesh.projector(Mesh$mesh,
                                xlim = range(Boundary[, 1]),
                                ylim = range(Boundary[, 2]),
                                dims = Nxy)

# Get the index of points on the grid within the boundary
xy.in <- splancs::inout(projgrid$lattice$loc, Boundary)

# Select only points on the grid that fall within the boundary
predcoords <- projgrid$lattice$loc[which(xy.in), ]
colnames(predcoords) <- c('x','y')
Apred <- projgrid$proj$A[which(xy.in), ]

# Construct the prediction data, which is spatial points by temporal layers
spatial_points <- cbind(c(Mesh$mesh$loc[,1]), c(Mesh$mesh$loc[,2]))
spatial_points <- projgrid$lattice$loc[which(xy.in), ]
colnames(spatial_points) <- c("x", "y")

NearestPredCovs <- GetNearestCovariate(points = spatial_points, covs = temporal_variables)

spatiotemporal_points <- rbind(spatial_points, spatial_points)
SP_Points_data <- data.frame(date_index = rep(c(1,2), each = nrow(spatiotemporal_points)/2))
SP_Points_data[, 'annual_rain_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_ann_rain_1980s_1.s, NearestPredCovs@data$TZ_ann_rain_2000s_1.s)
SP_Points_data[, 'annual_rain_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_ann_rain_1980s_2.s, NearestPredCovs@data$TZ_ann_rain_2000s_2.s)

SP_Points_data[, 'hottest_temp_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_max_temp_1980s_1.s, NearestPredCovs@data$TZ_max_temp_2000s_1.s)
SP_Points_data[, 'hottest_temp_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_max_temp_1980s_2.s, NearestPredCovs@data$TZ_max_temp_2000s_2.s)

SP_Points_data[, 'max_dryspell_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_dryspell_1980s_1.s, NearestPredCovs@data$TZ_dryspell_2000s_1.s)
SP_Points_data[, 'max_dryspell_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_dryspell_1980s_2.s, NearestPredCovs@data$TZ_dryspell_2000s_2.s)

SP_Points_data[, 'BG_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_BG_1980s_1.s, NearestPredCovs@data$TZ_BG_2000s_1.s)
SP_Points_data[, 'BG_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_BG_1980s_2.s, NearestPredCovs@data$TZ_BG_2000s_2.s)

SP_Points_data[, 'indicator'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$indicator_90s, NearestPredCovs@data$indicator_2010s)

SP_Points_data[, 'HFP_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_HFP_1980s_1.s, NearestPredCovs@data$TZ_HFP_2000s_1.s)
SP_Points_data[, 'HFP_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_HFP_1980s_2.s, NearestPredCovs@data$TZ_HFP_2000s_2.s)

# Add effort variables as constantly high effort across the projection grid
SP_Points_data[, 'duration_minutes'] <- 0.9
SP_Points_data[, 'effort'] <- 0.9

SP_Points_data <- cbind(SP_Points_data, spatiotemporal_points)

IP_sp <- sp::SpatialPointsDataFrame(coords = spatiotemporal_points, data = SP_Points_data, proj4string = proj)
IP_sp@data$Intercept <- 1
IP_sp@data[, c("x", "y")] <- IP_sp@coords

projmat.pred <- inla.spde.make.A(mesh = Mesh$mesh,
                                 loc = IP_sp@coords,
                                 group = IP_sp@data$date_index)

# ------------------------------------------------------------------------------------------------------------------------
# Prediction stack
# ------------------------------------------------------------------------------------------------------------------------
stk.pred <- inla.stack(tag='pred',
                       data = list(resp = NA),
                       A = list(1, projmat.pred),
                       effects = list(IP_sp@data,
                                      shared.field = index_set))  # Do not include separate eBird random field here

# ------------------------------------------------------------------------------------------------------------------------
# Integration stack
# ------------------------------------------------------------------------------------------------------------------------
coordnames <- colnames(temporal_variables@coords)
Names <- colnames(temporal_variables@data)

Names  <- c(coordnames, Names)

#get mesh triangle centroids
Points <- cbind(c(Mesh$mesh$loc[,1]), c(Mesh$mesh$loc[,2]))

#set column names of centroids to cordinate names from covariate data
colnames(Points) <- coordnames

#get value of nearest covariate to each mesh centroids
NearestCovs <- GetNearestCovariate(points=Points, covs=temporal_variables)

#set intercept to 1 for each mesh element
NearestCovs$Intercept <- rep(1,nrow(NearestCovs))

#add coordinate to data part of spatialpoints object if they are given
NearestCovs@data[, colnames(NearestCovs@coords)] <- NearestCovs@coords
covs_duplicated <- rbind(NearestCovs@data, NearestCovs@data)
covs_duplicated$index <- rep(c(1,2), each = nrow(covs_duplicated)/2)

# Projector matrix for integration points.
projmat.ip <- Matrix::Diagonal(n = n_time_layers*pcspde$n.spde)  # from mesh to integration points

stk.ip <- inla.stack(tag = "ip",
                     data = list(resp = NA, e = Mesh$w),  # 0 for count, NA for incidental, e is area of mesh polygon
                     A = list(1, projmat.ip, projmat.ip),
                     effects = list(covs_duplicated,
                                    shared.field = index_set,
                                    eBird.field = index_set_eBird))

integrated_stack <- inla.stack(stk.eBird, stk.atlas, stk.ip, stk.pred)

# ------------------------------------------------------------------------------------------------------------------------
# Model formula
# ------------------------------------------------------------------------------------------------------------------------
sdres <- sd(integrated_stack$effects$data$presence, na.rm = T)
h.spec <- list(rho = list(prior="pc.prec", param = c(0.5, 0.01)))

# At each time point, spatial locations are linked through the spde.
# Across time, the process evolves according to an AR(1) process.
form <- resp ~ 0 + ebird_intercept + atlas_intercept + 
      annual_rain_1 + annual_rain_2 + hottest_temp_1 + hottest_temp_2 + 
      max_dryspell_1 + max_dryspell_2 + 
      BG_1 + BG_2 + indicator + HFP_1 + HFP_2 +
      x + y +
      duration_minutes + effort + 
      f(shared.field, model = pcspde, group = shared.field.group,  # In each year, spatial locations are linked by spde
        control.group = list(model = 'ar1', hyper = h.spec)) +  # Across time, process evolves according to AR(1) process
      f(eBird.field, model = eBird_spde, group = eBird.field.group,
        control.group = list(model = 'ar1', hyper = h.spec))

# ------------------------------------------------------------------------------------------------------------------------
## Linear combinations -----------------
# ------------------------------------------------------------------------------------------------------------------------
# Make linear combinations, used to derive the variation explained by individual fixed effects in the model
TZ_lc_allFixed <-  inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                      atlas_intercept = rep(1, NROW(SP_Points_data)),
                                      annual_rain_1 = SP_Points_data$annual_rain_1,
                                      annual_rain_2 = SP_Points_data$annual_rain_2,
                                      max_dryspell_1 = SP_Points_data$max_dryspell_1,
                                      max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                      hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                      hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                      HFP_1 = SP_Points_data$HFP_1,
                                      HFP_2 = SP_Points_data$HFP_2,
                                      BG_1 = SP_Points_data$BG_1,
                                      BG_2 = SP_Points_data$BG_2,
                                      duration_minutes = SP_Points_data$duration_minutes,
                                      effort = SP_Points_data$effort,
                                      x = SP_Points_data$x,
                                      y = SP_Points_data$y,
                                      indicator = SP_Points_data$indicator)
names(TZ_lc_allFixed) <- paste0("TZ_lc_allFixed", 1:NROW(SP_Points_data))
TZ_lc_noTempMax <-  inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                       atlas_intercept = rep(1, NROW(SP_Points_data)),
                                       annual_rain_1 = SP_Points_data$annual_rain_1,
                                       annual_rain_2 = SP_Points_data$annual_rain_2,
                                       max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                       max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                       HFP_1 = SP_Points_data$HFP_1,
                                       HFP_2 = SP_Points_data$HFP_2,
                                       BG_1 = SP_Points_data$BG_1,
                                       BG_2 = SP_Points_data$BG_2,
                                       duration_minutes = SP_Points_data$duration_minutes,
                                       effort = SP_Points_data$effort,
                                       x = SP_Points_data$x,
                                       y = SP_Points_data$y,
                                       indicator = SP_Points_data$indicator)
names(TZ_lc_noTempMax) <- paste0("TZ_lc_noTempMax", 1:NROW(SP_Points_data))
TZ_lc_noRain <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                   atlas_intercept = rep(1, NROW(SP_Points_data)),
                                   max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                   max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                   hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                   hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                   HFP_1 = SP_Points_data$HFP_1,
                                   HFP_2 = SP_Points_data$HFP_2,
                                   BG_1 = SP_Points_data$BG_1,
                                   BG_2 = SP_Points_data$BG_2,
                                   duration_minutes = SP_Points_data$duration_minutes,
                                   effort = SP_Points_data$effort,
                                   x = SP_Points_data$x,
                                   y = SP_Points_data$y,
                                   indicator = SP_Points_data$indicator)
names(TZ_lc_noRain) <- paste0("TZ_lc_noRain", 1:NROW(SP_Points_data))
TZ_lc_noDryspell <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                       atlas_intercept = rep(1, NROW(SP_Points_data)),
                                       annual_rain_1 = SP_Points_data$annual_rain_1,
                                       annual_rain_2 = SP_Points_data$annual_rain_2,
                                       hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                       hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                       HFP_1 = SP_Points_data$HFP_1,
                                       HFP_2 = SP_Points_data$HFP_2,
                                       BG_1 = SP_Points_data$BG_1,
                                       BG_2 = SP_Points_data$BG_2,
                                       duration_minutes = SP_Points_data$duration_minutes,
                                       effort = SP_Points_data$effort,
                                       x = SP_Points_data$x,
                                       y = SP_Points_data$y,
                                       indicator = SP_Points_data$indicator)
names(TZ_lc_noDryspell) <- paste0("TZ_lc_noDryspell", 1:NROW(SP_Points_data))
TZ_lc_noBG <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                 atlas_intercept = rep(1, NROW(SP_Points_data)),
                                 annual_rain_1 = SP_Points_data$annual_rain_1,
                                 annual_rain_2 = SP_Points_data$annual_rain_2,
                                 max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                 max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                 hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                 hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                 HFP_1 = SP_Points_data$HFP_1,
                                 HFP_2 = SP_Points_data$HFP_2,
                                 duration_minutes = SP_Points_data$duration_minutes,
                                 effort = SP_Points_data$effort,
                                 x = SP_Points_data$x,
                                 y = SP_Points_data$y)
names(TZ_lc_noBG) <- paste0("TZ_lc_noBG", 1:NROW(SP_Points_data))
TZ_lc_noHFP <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                  atlas_intercept = rep(1, NROW(SP_Points_data)),
                                  annual_rain_1 = SP_Points_data$annual_rain_1,
                                  annual_rain_2 = SP_Points_data$annual_rain_2,
                                  max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                  max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                  hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                  hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                  BG_1 = SP_Points_data$BG_1,
                                  BG_2 = SP_Points_data$BG_2,
                                  duration_minutes = SP_Points_data$duration_minutes,
                                  effort = SP_Points_data$effort,
                                  x = SP_Points_data$x,
                                  y = SP_Points_data$y,
                                  indicator = SP_Points_data$indicator)
names(TZ_lc_noHFP) <- paste0("TZ_lc_noHFP", 1:NROW(SP_Points_data))

lc_combined <- c(all_lc, TZ_lc_noRain, TZ_lc_noTempMax, TZ_lc_noDryspell, TZ_lc_noBG, TZ_lc_noHFP, TZ_lc_allFixed)

# ------------------------------------------------------------------------------------------------------------------------
# INLA model
# ------------------------------------------------------------------------------------------------------------------------
# Model output string
model_name <- paste0(sub(" ", "_", species), "_r", 
                     prior_range[1], "_", prior_range[2], "_s", prior_sigma[1], "_", prior_sigma[2], ".RData")

# Gaussian priors
C.F. <- list(
      mean = fixed_mean,
      prec = fixed_precision   # Precision for all fixed effects except intercept
)

model <- inla(form, family = "binomial", control.family = list(link = "cloglog"), 
              lincomb = all_lc,
              data = inla.stack.data(integrated_stack), 
              verbose = FALSE,
              control.predictor = list(A = inla.stack.A(integrated_stack), 
                                       link = NULL, compute = TRUE), 
              control.fixed = C.F.,
              E = inla.stack.data(integrated_stack)$e, 
              control.compute = list(waic = FALSE, dic = FALSE, cpo = TRUE))

logCPO_vect = log(model$cpo$cpo[model$cpo$cpo != 0])
logCPO_vect = logCPO_vect[is.finite(logCPO_vect)]

logCPO = round(-sum(logCPO_vect, na.rm = T), digits = 2)

save(model, all.seq, integrated_stack, SP_Points_data, predcoords,
     TZ_no_lakes, TZ_outline, projgrid, index_set, index_set_eBird,
     xy.in, species, prior_range, prior_sigma, logCPO, max.edge, 
     IP_sp, file = paste0('model_output/', model_name))

# ------------------------------------------------------------------------------------------------------------------------
# Effect plots 
# ------------------------------------------------------------------------------------------------------------------------
original_values <- data.frame(orig_values = unlist(all.seq)) %>% 
      mutate(covariate = sub("*.seq\\d+", "", rownames(.)),
             sequence = as.numeric(gsub("\\D", "", rownames(.))))

# Scale effect plots so that flat trends lie around the median
scale_params <- median(model$summary.lincomb.derived$`0.5quant`[grep("TZ_HFP|TZ_ann_rain|TZ_BG|TZ_dryspell|TZ_max_temp|TZ_HFP", 
                                                                     rownames(model$summary.lincomb.derived))])
effect_combs <- data.frame(covariate = gsub('[[:digit:]]+', '', sub("*_lc\\d+", "", rownames(model$summary.lincomb.derived))),
                           sequence = as.numeric(gsub("\\D", "", rownames(model$summary.lincomb.derived))),
                           quant_05 = inla.link.cloglog((model$summary.lincomb.derived$`0.5quant` - scale_params) - 0.36, inverse = TRUE),
                           quant_0025 = inla.link.cloglog((model$summary.lincomb.derived$`0.025quant` - scale_params) - 0.36, inverse = TRUE),
                           quant_0975 = inla.link.cloglog((model$summary.lincomb.derived$`0.975quant` - scale_params) - 0.36, inverse = TRUE))
effect_combs_main <- effect_combs %>% 
      filter(covariate %in% c("TZ_ann_rain", "TZ_BG", "TZ_dryspell", "TZ_max_temp", "TZ_HFP"))
effect_combs_m <- merge(original_values, effect_combs_main)

facet_labels <- c(
      TZ_ann_rain = "Annual rainfall",
      TZ_BG = "Bareground cover",
      TZ_dryspell = "Longest dryspell duration",
      TZ_max_temp = "Hottest temperature",
      TZ_HFP = "Human footprint"
)

# Make effect plots, using linear combinations
effects_plot <- ggplot(effect_combs_m) +
      geom_line(aes(x = orig_values, y = quant_05)) +
      geom_line(aes(x = orig_values, y = quant_0025), lty = 2, alpha = .5) +
      geom_line(aes(x = orig_values, y = quant_0975), lty = 2, alpha = .5) +
      facet_wrap(~ covariate, scale = 'free_x', labeller = as_labeller(facet_labels)) +
      theme_few() +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlab("") +
      ylab("Probability of occurrence"); effects_plot

ggsave(plot = effects_plot, filename = paste0("figures/effects_", sub(" ", "_", species), "_E", max.edge, "_r", 
                                              prior_range[1], "_", prior_range[2], "_s", 
                                              prior_sigma[1], "_", prior_sigma[2], 
                                              ".png"))

# ------------------------------------------------------------------------------------------------------------------------
# Prediction plots
# ------------------------------------------------------------------------------------------------------------------------
TZ_no_lakes_vect <- vect(TZ_no_lakes)

pred.index <- inla.stack.index(stack = integrated_stack, tag = "pred")$data

IP_df <- data.frame(SP_Points_data) %>% dplyr::select(date_index, x, y)

IP_df$pred_median <- model$summary.fitted.values[pred.index, "0.5quant"]
IP_df$pred_sd <- model$summary.fitted.values[pred.index, "sd"]

pred_data <- data.frame(all_pred = inla.link.cloglog(IP_df$pred_median, inv = T),
                        sd = IP_df$pred_sd,
                        ind = rep(c(1,2), each = length(IP_df$pred_median)/2))
predcoordsGroup <- do.call(rbind, list(predcoords, predcoords))
pred_data_coords <- cbind(pred_data, predcoordsGroup) %>% 
      dplyr::select(x, y, everything())

preds_1 <- pred_data_coords %>% 
      filter(ind == 1) %>% 
      rast(., type = "xyz") %>% 
      mask(., TZ_no_lakes_vect)
      
preds_2 <- pred_data_coords %>% 
      filter(ind == 2) %>% 
      rast(., type = "xyz") %>% 
      mask(., TZ_no_lakes_vect)

median_plot_1 <- ggplot() +
      geom_spatraster(data = preds_1$all_pred) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "black") +
      coord_sf() +
      ggtitle("1980-1999") +
      scale_fill_viridis(name = "P(occurrence)", na.value="transparent") +
      theme_bw()

median_plot_2 <- ggplot() +
      geom_spatraster(data = preds_2$all_pred) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "black") +
      coord_sf() +
      ggtitle("2000-2020") +
      scale_fill_viridis(name = "P(occurrence)", na.value="transparent") +
      theme_bw() 

sd_plot_1 <- ggplot() +
      geom_spatraster(data = preds_1$sd) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "black") +
      coord_sf() +
      ggtitle("1980-1999") +
      scale_fill_viridis(name = "SD", na.value="transparent") +
      theme_bw()

sd_plot_2 <- ggplot() +
      geom_spatraster(data = preds_2$sd) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "black") +
      coord_sf() +
      ggtitle("2000-2020") +
      scale_fill_viridis(name = "SD", na.value="transparent") +
      theme_bw()

# ------------------------------------------------------------------------------------------------------------------------
# Posterior median of the space-time random field
# ------------------------------------------------------------------------------------------------------------------------
xmedian <- list()
for (j in 1:2) {
      xmedian[[j]] <- inla.mesh.project(
            projgrid,  model$summary.random$shared.field$`0.5quant`[index_set$shared.field.group == j])
}

xy.inGroup <- c(xy.in,xy.in)
xmedian <- unlist(xmedian)
xmedian <- xmedian[xy.inGroup]

dataObj <- data.frame(random_shared = xmedian,
                      ind = rep(c(1,2), each = length(xmedian)/2))
shared_random <- cbind(dataObj, predcoordsGroup) %>% 
      dplyr::select(x, y, everything())

shared_1 <- shared_random %>% 
      filter(ind == 1) %>% 
      rast(., type = "xyz") %>% 
      mask(., TZ_no_lakes_vect)

shared_2 <- shared_random %>% 
      filter(ind == 2) %>% 
      rast(., type = "xyz") %>% 
      mask(., TZ_no_lakes_vect)

random_plot_1 <- ggplot() +
      geom_spatraster(data = shared_1$random_shared) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "black") +
      coord_sf() +
      ggtitle("1980-1999 - shared random field") +
      scale_fill_viridis(name = "Median", na.value="transparent") +
      theme_bw()

random_plot_2 <- ggplot() +
      geom_spatraster(data = shared_2$random_shared) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "black") +
      coord_sf() +
      ggtitle("2000-2020 - shared random field") +
      scale_fill_viridis(name = "Median", na.value="transparent") +
      theme_bw()

xmedian2 <- list()
for (j in 1:2) {
      xmedian2[[j]] <- inla.mesh.project(
            projgrid,  model$summary.random$eBird.field$`0.5quant`[index_set_eBird$eBird.field.group == j])
}

xmedian2 <- unlist(xmedian2)
xmedian2 <- xmedian2[xy.inGroup]

dataObj2 <- data.frame(random_eBird = xmedian2,
                       ind = rep(c(1,2), each = length(xmedian2)/2))

eBird_random <- cbind(dataObj2, predcoordsGroup) %>% 
      dplyr::select(x, y, everything())

eBird_random_1 <- eBird_random %>% 
      filter(ind == 1) %>% 
      rast(., type = "xyz") %>% 
      mask(., TZ_no_lakes_vect)

eBird_random_2 <- eBird_random %>% 
      filter(ind == 2) %>% 
      rast(., type = "xyz") %>% 
      mask(., TZ_no_lakes_vect)

eBird_random_plot_1 <- ggplot() +
      geom_spatraster(data = eBird_random_1$random_eBird) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "black") +
      coord_sf() +
      ggtitle("1980-1999 - eBird random field") +
      scale_fill_viridis(name = "Median", na.value="transparent") +
      theme_bw()

eBird_random_plot_2 <- ggplot() +
      geom_spatraster(data = eBird_random_1$random_eBird) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "black") +
      coord_sf() +
      ggtitle("2000-2020 - eBird random field") +
      scale_fill_viridis(name = "Median", na.value="transparent") +
      theme_bw()

maps_combined <- (median_plot_1 | median_plot_2 | sd_plot_1 | sd_plot_2) / 
      (random_plot_1 | random_plot_2 | eBird_random_plot_1 | eBird_random_plot_2) / effects_plot + 
      plot_annotation(title = paste0(species, ", ", logCPO),
                      theme = theme(plot.title = element_text(size = 22)))

pdf(paste0("figures/output_", sub(" ", "_", species), "_E", max.edge, "_r",
           prior_range[1], "_", prior_range[2], "_s",
           prior_sigma[1], "_", prior_sigma[2], ".pdf"),
    width = 20, height = 20)
maps_combined
dev.off()


