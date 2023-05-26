rm(list = ls())
library(tidyverse)
library(INLA)
library(MASS)

# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
source("source/misc_functions.R")
load('model_data/regression_data.RData')

#  Get sample sizes by trait levels
model_data %>% 
      group_by(Primary.Lifestyle) %>% 
      summarise(N = n())
model_data %>% 
      group_by(Trophic.Level) %>% 
      summarise(N = n())
model_data %>% 
      group_by(Migratory_ability) %>% 
      summarise(N = n())
mean(model_data$auc_mean)
min(model_data$auc_mean)
max(model_data$auc_mean)

# -----------------------------------------------------------------------------------------------------------------
# Statistical test
# -----------------------------------------------------------------------------------------------------------------
form_prop_change <- log_prop_change ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + HWI + Mass + avg.r

model_prop_change <- inla(form_prop_change, data = model_data, 
                          family = "T",  
                          lincomb = all_lc_sens,
                          control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                          control.predictor=list(compute=TRUE))

hist(model_prop_change$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_prop_change$cpo$pit))),
       model_prop_change$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_prop_change$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))

form_col_log <- log_relative_colonisation ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + HWI + Mass + avg.r

model_col_log <- inla(form_col_log, data = model_data, 
                         family = "T",  
                         lincomb = all_lc_sens,
                         control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                         control.predictor=list(compute=TRUE))

hist(model_col_log$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_col_log$cpo$pit))), model_col_log$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_col_log$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))

form_loss_log <- log_relative_extinction ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + HWI + Mass + avg.r

model_loss_log <- inla(form_loss_log, data = model_data, 
                              family = "T", 
                              lincomb = all_lc_sens,
                              control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                              control.predictor=list(compute=TRUE))
summary(model_loss_log)
hist(model_loss_log$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_loss_log$cpo$pit))), model_loss_log$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_loss_log$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))

form_prop_change_breadth <- log_prop_change ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      rain_breadth + temp_breadth + dry_breadth + HFP_breadth + BG_breadth + HWI + Mass + avg.r

model_prop_change_breadth <- inla(form_prop_change_breadth, data = model_data, 
                                 family = "T", 
                                 lincomb = all_lc_breadth,
                                 control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                 control.predictor=list(compute=TRUE))

hist(model_prop_change_breadth$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_prop_change_breadth$cpo$pit))), model_prop_change_breadth$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_prop_change_breadth$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))

form_col_log_breadth <- log_relative_colonisation ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      rain_breadth + temp_breadth + dry_breadth + HFP_breadth + BG_breadth + HWI + Mass + avg.r

model_col_log_breadth <- inla(form_col_log_breadth, data = model_data, 
                                 family = "T",  
                                 lincomb = all_lc_breadth,
                                 control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                 control.predictor=list(compute=TRUE))

hist(model_col_log_breadth$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_col_log_breadth$cpo$pit))), model_col_log_breadth$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_col_log_breadth$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))

form_loss_log_breadth <- log_relative_extinction ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      rain_breadth + temp_breadth + dry_breadth + HFP_breadth + BG_breadth + HWI + Mass + avg.r

model_loss_log_breadth <- inla(form_loss_log_breadth, data = model_data, 
                                      family = "T",  
                                      lincomb = all_lc_breadth,
                                      control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                      control.predictor=list(compute=TRUE))

hist(model_loss_log_breadth$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_loss_log_breadth$cpo$pit))), model_loss_log_breadth$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_loss_log_breadth$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))

save(model_prop_change, model_col_log, model_loss_log, 
     model_prop_change_breadth, model_col_log_breadth,model_loss_log_breadth, file = "results/traits_model_list.RData")





