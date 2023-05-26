sortBase <- function(vec, nsplines = 2, method, user_cp_quantiles = NULL, vecname, display_plot) {
      #' Thin-plate spline basis function, modified from code
      #' provided in Crainiceanu, C., Ruppert, D. & Wand, M.P. Bayesian analysis for
      #' penalized spline regression using WinBUGS. J. Stat. Soft. 14, 1?24(2005).
      #'
      #' This function computes a thin-plate spline basis for a given numeric vector
      #' using a specified number of basis functions (nsplines).
      #' The basis functions are centered at control points determined by either the quantiles
      #' or kmeans clusters of the unique values in the input vector.
      #'
      #' @param vec A numeric vector for which the thin-plate spline basis is computed
      #' @param nsplines Integer specifying the number of basis functions (default: 2)
      #' @param methods Character string specifying the control point placement method
      #' @param user_cp_quantiles Optional numeric vector specifying user-defined quantiles for placing control points
      #' @return A matrix with the same number of rows as the input vector and nsplines columns,
      #'         representing the thin-plate spline basis for the input vector.

      if (!is.null(user_cp_quantiles)){
            control_points <- quantile(unique(vec), probs = user_cp_quantiles, na.rm = TRUE)
            plot_title <- paste0(vecname, ", Cp quantiles provided by user")
      } else {
            if (method == "quantile"){
                  control_points <- quantile(unique(vec), seq(0, 1, length = (nsplines+2))[-c(1, (nsplines+2))], na.rm = TRUE)
                  label = "quantiles"
            } else if (method == "kmeans"){
                  clusters <- kmeans(unique(na.omit(vec)), centers = nsplines)
                  control_points <- sort(clusters$centers)
                  label = "kmeans"
            } else if (method == "gmm"){
                  require(mixtools)
                  fit <- normalmixEM(unique(na.omit(vec)), k = nsplines)
                  control_points <- fit$mu
                  label = "Gaussian mixture model"
            }
            plot_title <- paste0(vecname, ", Cp based on ", label)
      }
      
      if(display_plot == TRUE){
            plot(density(na.omit(vec)), main = plot_title)
            abline(v = control_points, col = "red")
      }

      # Define the thin-plate spline penalty matrix
      zFE       <- cbind(rep(1, length(vec) ), vec)

      z_K <- (abs(outer(vec, control_points, "-"))) ^ 3
      OMEGA.all <- (abs(outer(control_points, control_points, "-"))) ^ 3
      svd.OMEGA.all  <- svd(OMEGA.all)
      sqrt.OMEGA.all <- t(svd.OMEGA.all$v %*% (t(svd.OMEGA.all$u) *
                                                     sqrt(svd.OMEGA.all$d)))
      z.out     <- t(solve(sqrt.OMEGA.all, t(z_K)))
      return(z.out)
}

prepareSplineData <- function(df, vector, nsplines = 2, method = "quantile", user_cp_quantiles = NULL, display_plot = TRUE){
      #' Prepare spline data for a given vector in a data frame
      #'
      #' This function computes the spline basis for a given vector in a data frame using
      #' either natural splines or the custom sortBase method. The resulting spline basis
      #' is added to the data frame and prediction vectors are assigned to the global environment.
      #'
      #' @param df A data frame containing the input vector
      #' @param vector The vector within the data frame for which the spline basis is computed
      #' @param nsplines Integer specifying the number of basis functions (default: 2)
      #' @param methods Character string specifying the control point placement method (default: "quantile")
      #' @param user_cp_quantiles Optional numeric vector specifying user-defined quantiles for placing control points
      #' @return A modified data frame with added columns for each of the computed spline basis functions
      
      # Create prediction vector, as sequence of 100 values, spanning full range of values. 
      sequence <- seq(min(vector, na.rm = T), max(vector, na.rm = T), length = 100)
      # Save this vector, we use it later to label the x axis in the model effect plots
      assign(paste0(deparse(substitute(vector)), ".seq"), sequence, envir = .GlobalEnv)
      
      # Add prediction vector to beginning of original values, scale all to mean 0 and sd 1
      var.s <- c(scale(c(sequence, vector)))  
      
      n_data <- NROW(vector)/2
      vec_name <- sub("^.*\\$", "", deparse(substitute(vector)))
      
      z.var.s <- scale(sortBase(vec = var.s, nsplines = nsplines, method, user_cp_quantiles, vecname = vec_name, display_plot = display_plot))
      for (i in 1:nsplines) {
            assign(paste0(deparse(substitute(vector)), "_", i, ".s"), c(z.var.s[1:100, i]), envir = .GlobalEnv)
            df[, paste0(deparse(substitute(vector)), "_1980s_", i, ".s")] <- c(z.var.s[101:(100+n_data), i])
            df[, paste0(deparse(substitute(vector)), "_2000s_", i, ".s")] <- c(z.var.s[(n_data+101):NROW(z.var.s), i])
      }
      return(df)
}

make_forest_plot <- function(inla_model, label, col){
      #' This function creates a forest plot based on the coefficients of an INLA model. The forest plot visualizes the posterior estimates of the model's fixed effects. Depending on the label argument, the function can be used to create a forest plot for different purposes, such as visualizing the importance of different species traits for niche breadth or sensitivity.
      #'
      #' @param inla_model An INLA model object from which the coefficients are extracted.
      #' @param label A character string indicating the type of the forest plot, either "Niche breadth" or "Sensitivity".
      #' @param col A color to be used for the background of the plot.
      #'
      #' @return A ggplot object representing the forest plot. Each point in the plot represents a fixed effect in the model, with horizontal error bars showing the 95% credible interval for the effect, and the color indicating whether the effect is statistically significant.
      #'
      #' @details The function first extracts the coefficients of the fixed effects from the INLA model, then relabels the coefficients for better interpretability in the plot. Depending on the label argument, it adjusts the labels appropriately. The function then creates a ggplot object representing the forest plot, where each point corresponds to a fixed effect, the horizontal position of the point represents the median of the posterior distribution of the effect, and the horizontal error bar shows the 95% credible interval. The color of the point and the error bar indicates whether the effect is statistically significant, based on whether the 95% credible interval excludes zero.
      #'
      #' The function uses the ggplot2 package for plotting, and thus the returned object can be further customized using ggplot2 functions.
      
      model_fixed <- data.frame(inla_model$summary.fixed) %>% 
            mutate(ID = rownames(.))
      
      if(label == "Niche breadth"){
            model_fixed <- model_fixed %>% 
                  mutate(ID = replace(ID, ID == "(Intercept)", "Intercept"),
                         ID = replace(ID, ID == "Trophic.LevelOmnivore", "Trophic level: Omnivore"),
                         ID = replace(ID, ID == "Trophic.LevelHerbivore", "Trophic level: Herbivore"),
                         ID = replace(ID, ID == "Primary.LifestyleGeneralist", "Primary lifestyle: Generalist"),
                         ID = replace(ID, ID == "Primary.LifestyleAerial", "Primary lifestyle: Aerial"),
                         ID = replace(ID, ID == "Primary.LifestyleTerrestrial", "Primary lifestyle: Terrestrial"),
                         ID = replace(ID, ID == "Migratory_abilityhigh", "Migratory ability: High"),
                         ID = replace(ID, ID == "Migratory_abilitymoderate", "Migratory ability: Moderate"),
                         ID = replace(ID, ID == "BG_breadth", "Bare ground cover"),
                         ID = replace(ID, ID == "rain_breadth", "Annual rainfall"),
                         ID = replace(ID, ID == "temp_breadth", "Hottest temperature"),
                         ID = replace(ID, ID == "dry_breadth", "Dry spell duration"),
                         ID = replace(ID, ID == "HFP_breadth", "Human footprint"),
                         ID = replace(ID, ID == "HWI", "Hand-wing Index"),
                         ID = replace(ID, ID == "Mass", "Body mass"),
                         ID = replace(ID, ID == "avg.r", "Dorsal reflectance"))
      }
      
      if(label == "Sensitivity"){
            model_fixed <- model_fixed %>% 
                  mutate(ID = replace(ID, ID == "(Intercept)", "Intercept"),
                         ID = replace(ID, ID == "Trophic.LevelOmnivore", "Trophic level: Omnivore"),
                         ID = replace(ID, ID == "Trophic.LevelHerbivore", "Trophic level: Herbivore"),
                         ID = replace(ID, ID == "Primary.LifestyleGeneralist", "Primary lifestyle: Generalist"),
                         ID = replace(ID, ID == "Primary.LifestyleAerial", "Primary lifestyle: Aerial"),
                         ID = replace(ID, ID == "Primary.LifestyleTerrestrial", "Primary lifestyle: Terrestrial"),
                         ID = replace(ID, ID == "Migratory_abilityhigh", "Migratory ability: High"),
                         ID = replace(ID, ID == "Migratory_abilitymoderate", "Migratory ability: Moderate"),
                         ID = replace(ID, ID == "BG_imp", "Bare ground cover"),
                         ID = replace(ID, ID == "rain_imp", "Annual rainfall"),
                         ID = replace(ID, ID == "temp_imp", "Hottest temperature"),
                         ID = replace(ID, ID == "dry_imp", "Dry spell duration"),
                         ID = replace(ID, ID == "HFP_imp", "Human footprint"),
                         ID = replace(ID, ID == "HWI", "Hand-wing Index"),
                         ID = replace(ID, ID == "Mass", "Body mass"),
                         ID = replace(ID, ID == "avg.r", "Dorsal reflectance"))
      }
      
      model_fixed$ID <- factor(model_fixed$ID, levels = c("Trophic level: Omnivore", "Trophic level: Herbivore", 
                                                          "Primary lifestyle: Generalist", "Primary lifestyle: Aerial", "Primary lifestyle: Terrestrial",
                                                          "Migratory ability: High", "Migratory ability: Moderate", "Hand-wing Index", "Body mass", "Dorsal reflectance",
                                                          "Bare ground cover", "Annual rainfall", "Hottest temperature", "Dry spell duration", "Human footprint",
                                                          "Intercept"))
      model_fixed$significant <- ifelse((model_fixed$X0.025quant > 0 & model_fixed$X0.975quant > 0)|(model_fixed$X0.025quant < 0 & model_fixed$X0.975quant < 0), "yes", "no")
      
      forest_plot <- ggplot() + 
            annotate("text", y = 1.6, x = max(inla_model$summary.fixed$`0.975quant`) - 0.1*max(inla_model$summary.fixed$`0.975quant`), 
                     col = col, label = label,
                     hjust = 0,
                     angle = 90,
                     size = 3.5) +
            geom_rect(aes(ymin = 1.5, ymax = 6.5, xmin = min(inla_model$summary.fixed$`0.025quant`), 
                          xmax = max(inla_model$summary.fixed$`0.975quant`)),
                      fill = col, col = col, alpha = .3) +
            geom_point(data = model_fixed, aes(y = ID, x = X0.5quant, col = significant)) +
            geom_errorbar(data = model_fixed, aes(y = ID, xmin = X0.025quant, xmax = X0.975quant, col = significant), width = 0.1) +
            geom_vline(aes(xintercept = 0), lty = 2, alpha = .3) +
            xlab("Posterior estimates") +
            ylab(element_blank()) +
            theme_minimal() +
            scale_y_discrete(limits=rev) +
            scale_colour_manual(values = c("black", "#D55E00"), guide = "none") 
      return(forest_plot)
}

rsq <- function(x, y) summary(lm(y~x))$r.squared

unscale_fun <- function(attr_column, value = NULL){
      #' This function unscales a scaled variable back to its original scale. Scaling is a common pre-processing step in machine learning and statistics, which involves subtracting the mean and dividing by the standard deviation. This function does the reverse operation, multiplying by the standard deviation and adding the mean.
      #'
      #' @param attr_column A numeric vector that has been scaled. This vector should have the attributes 'scaled:scale' (the standard deviation of the original variable) and 'scaled:center' (the mean of the original variable), which are typically added by the scale() function in R.
      #' @param value An optional numeric value to be unscaled. If provided, this function will unscale the provided value instead of the entire vector.
      #'
      #' @return A numeric vector or single value that has been unscaled back to the original scale of the variable.
      #'
      #' @details This function checks if a specific value is provided to be unscaled. If not, it unscales the entire vector. This can be used to interpret the coefficients of a model that was fitted with scaled predictors, as it allows to convert these coefficients back to the original scale of the predictors.
      
      if(is.null(value)){
            unscaled <- attr_column * attr(attr_column, 'scaled:scale') + attr(attr_column, 'scaled:center')
      } else {
            unscaled <- value * attr(attr_column, 'scaled:scale') + attr(attr_column, 'scaled:center')
      }
      return(unscaled)
}

get_slope_unscaled <- function(model, data, covar, custom_unit_step, log_transformed = NULL){
      #' This function is used to extract and adjust the slope of a given covariate from a log-linear model, making it interpretable on the original scale of the data. It can handle models where the response variable and certain covariates were log-transformed prior to modeling.
      #'
      #' @param model An object representing the fitted model, from which the coefficients are to be extracted.
      #' @param data The original dataset used in the model, needed to obtain the scale of the covariates.
      #' @param covar A string representing the name of the covariate for which the slope is to be extracted and adjusted.
      #' @param custom_unit_step An optional numerical value representing a custom unit step size to be used in adjusting the slope. If NULL, the function will calculate an appropriate unit step size based on the scale of the covariate in the original data.
      #' @param log_transformed A vector of character strings representing the names of covariates that were log-transformed in the model. This is used to adjust the slopes of these covariates.
      #'
      #' @return The function does not return a value but prints the adjusted slope (expressed as a percent change in the response per unit change in the covariate on the original scale) and the unit step size used in adjusting the slope.
      #'
      #' @details The function operates by extracting the model's coefficients and adjusting the slope of the specified covariate. If the covariate was log-transformed in the model (as indicated by covar being present in log_transformed), the function adjusts the slope to express it in terms of a unit change on the original scale of the covariate. The adjustment involves dividing the log-transformed slope by the value representing a one-unit step on the original scale of the covariate. This value is obtained using the `unscale_fun` function.
      #' 
      #' The function also calculates an appropriate unit step size based on the scale of the covariate in the original data, unless a custom unit step size is specified by the user. The slope is then adjusted for this unit step size
      
      # Extract coefficients from model
      model_fixed <- data.frame(model$summary.fixed) %>% 
            mutate(ID = rownames(.))
      
      # Calculate slope 
      slope <- model_fixed$mean[model_fixed$ID == covar] * 100
      
      # Adjust for log transformation
      if(covar %in% log_transformed) {
            slope <- (exp(model_fixed$mean[model_fixed$ID == covar] / unscale_fun(value = 1, data[, covar])) - 1) * 100
      }
      
      # Calculate unit step
      unit_step <- round(unscale_fun(value = 1, data[, covar]), digits = 2)
      
      if(is.null(custom_unit_step)){
            new_unit_step = 10^floor(log10(abs(unit_step))) 
      } else {
            new_unit_step = custom_unit_step
      }
      
      # Adjust slope for new unit step
      new_slope = (slope * new_unit_step) / unit_step
      
      # Format output
      new_slope <- sprintf("%.1f%%", new_slope)  
      
      # Print result
      print(sprintf("ΔX: %s, ΔY: %s", new_unit_step, new_slope))
}

make_raw_data_plots <- function(inla_model, data, covariate = NULL, newscale.y = NULL, custom_unit_step = NULL){
      #' This function generates a list of ggplot2 plots of raw data for each covariate in the INLA model against the response variable. It also adds estimated lines and uncertainty ribbons for numeric covariates, or point estimates for factor covariates, derived from the INLA model.
      #'
      #' @param inla_model An INLA model object.
      #' @param data The data frame that was used to fit the INLA model.
      #' @param covariate A character string specifying a single covariate to plot. If NULL, the function will generate plots for all covariates in the INLA model.
      #' @param newscale.y A custom y scale, created with scale_y_continuous() or similar. If NULL, the default y scale is used.
      #' @param custom_unit_step An optional custom unit step for the get_slope_unscaled function, which adjusts the slope of the estimate for the unit step.
      #'
      #' @return A list of ggplot2 objects, each being a plot for one covariate.
      #'
      #' @details The function generates plots using ggplot2. For numeric covariates, it creates a scatter plot of the raw data with an estimated line and uncertainty ribbons derived from the INLA model. For factor covariates, it creates a boxplot of the raw data with point estimates derived from the INLA model. If the covariate is "Mass", the function applies a custom x scale. In addition, the function generates labels for x and y axes, and the plot title is the output of the get_slope_unscaled function for the given covariate.
      #'
      #' The helper function create_label() retrieves a human-readable label for a given variable name, and the helper function create_plot() generates a ggplot2 object for a given covariate.
      
      # Helper functions
      create_label <- function(var_name, label_df) {
            label_df$label_name[label_df$var_name == var_name]
      }
      
      create_plot <- function(covariate) {
            if (is.numeric(data[, covariate])) {
                  p <- ggplot(data, aes(x = get(covariate))) +
                        geom_point(aes(y = get(response))) +
                        geom_line(data = lincomb_data[grepl(covariate, rownames(lincomb_data)) == TRUE, ],
                                  aes(x = seq(min(data[, covariate], na.rm = TRUE), 
                                              max(data[, covariate], na.rm = TRUE), 
                                              len = 100),
                                      y = X0.5quant), 
                                  lty = 2, alpha = 0.5, col = "red") +
                        geom_ribbon(data = lincomb_data[grepl(covariate, rownames(lincomb_data)) == TRUE, ],
                                    aes(x = seq(min(data[, covariate], na.rm = TRUE), 
                                                max(data[, covariate], na.rm = TRUE), 
                                                len = 100), 
                                        ymin = X0.025quant, 
                                        ymax = X0.975quant), 
                                    alpha = 0.1) +
                        theme_linedraw() +
                        xlab(create_label(covariate, xlab_df)) +
                        ylab(create_label(response, ylab_df)) + 
                        geom_rug(aes(x = get(covariate), y = get(response)), 
                                 col = rgb(.5, 0, 0, alpha = .4)) +
                        ggtitle(get_slope_unscaled(inla_model, data, covariate, custom_unit_step, "Mass"))
            } else {
                  p <- ggplot(data, aes(x = get(covariate))) +
                        geom_boxplot(aes(y = get(response))) +
                        geom_point(data = lincomb_data[grepl(covariate, rownames(lincomb_data)) == TRUE, ], 
                                   aes(x = levels(data[, covariate]), y = X0.5quant), 
                                   col = "red", cex = 2) +
                        theme_linedraw() +
                        xlab(create_label(covariate, xlab_df)) +
                        ylab(create_label(response, ylab_df))
            }
            if(covariate == "Mass") p <- p + scale_x_continuous(breaks = log(c(seq(4, 10, by = 2),
                                                                               15, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000, 10000)), 
                                                                labels = c(seq(4, 10, by = 2),
                                                                           15, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000, 10000))
            if(!is.null(newscale.y)) p <- p + newscale.y
            print(p)
      }
      
      lincomb_data <- data.frame(inla_model$summary.lincomb.derived) %>% 
            mutate(ID = rownames(.)) %>% 
            dplyr::select(-mode, -kld, -mean, -sd)
      xlab_df <- data.frame(var_name = c('BG_imp', 'rain_imp', 'temp_imp', 'dry_imp', 'HFP_imp', 
                                         'rain_breadth', 'temp_breadth', 'dry_breadth', 'HFP_breadth', 'BG_breadth',
                                         'HWI', 'avg.r', 'Mass',
                                         'Trophic.Level', 'Primary.Lifestyle', 'Migratory_ability'), 
                            label_name = c('Bare ground sensitivity', 'Annual rainfall sensitivity', 'Hottest temperature sensitivity', 'Dry spell duration sensitivity', 'Human footprint sensitivity', 
                                           'Annual rainfall niche breadth', 'Hottest temperature niche breadth', 'Dry spell duration niche breadth', 'Human footprint niche breadth', 'Bare ground niche breadth',
                                           'Hand-wing index', 'Average dorsal reflectance', 'Body mass',
                                           'Trophic level', 'Primary lifestyle', 'Migratory ability'))
      
      ylab_df <- data.frame(var_name = c('relative_colonisation', 'relative_extinction', 'log_prop_change'), 
                            label_name = c('Meaningful colonisation', 'Meaningful extinction', 'Proportional change'))
      
      response <- as.character(inla_model$.args$formula[2])
      
      covariates <- if(is.null(covariate)) {
            strsplit(as.character(inla_model$.args$formula[3]), split = " \\+ ")[[1]]
      } else {
            covariate
      }
      
      plot_list <- lapply(covariates, create_plot)
      
      return(plot_list)
}



