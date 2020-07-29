# Example code to predict temperature-specific R0 curves for schistosomiasis with and without snail control
# Supplement to "Interventions can shift the thermal optimum for disease: Evidence from human schistosomiasis" by K.H. Nguyen, P.H. Boersch-Supan, V.J. Harwood, J.R. Rohr (2020)
# https://github.com/pboesu/schistosomiasis

library(dplyr)
library(readr)

#implementation of the R_0 function from Gao, S., Liu, Y., Luo, Y. & Xie, D. Control problems of a mathematical model for schistosomiasis transmission dynamics. Nonlinear Dynamics 63, 503-512, doi:10.1007/s11071-010-9818-z (2010).
gao_R_0 <- function(params){
  with(params, {  
    #calculate compound parameters
    gamma = k * Lambda_1 * gamma_1 / M_0
    d_1 = mu_1 + delta_1 + eta
    d_2 = mu_2 + theta
    d_3 = mu_2 + delta_2 + theta
    d_4 = mu_4 + tau
    delta = gamma_2 * Lambda_2
    #calculate R_0
    R_0 = sqrt((beta_1 * beta_2 * gamma * delta)/ (d_1 * d_2 * d_3 * d_4 * mu_1 * mu_3))
    R_0
  })
}

#read in predicted thermal performance curves and constants.
clean_parameters_wide <- read_rds('data/clean_parameters_wide.rds')

#set control parameters (additional mortality on snails and cercaria; human treatment rate) to zero
clean_parameters_wide %>% mutate(eta = 0,
                                 theta = 0,
                                 tau = 0) -> clean_parameters_wide_no_control 

#set snail control to 20%
clean_parameters_wide %>% mutate(eta = 0,
                                 theta = 0.2,
                                 tau = 0) -> clean_parameters_wide_snail_control 
# we can now loop over the rows of this data frame to predict an R_0 value for each temperature

#prepare a vector for the results that is the same size as there are rows in the data
R_0_temp_no_control <- numeric(nrow(clean_parameters_wide))
R_0_temp_snail_control <- numeric(nrow(clean_parameters_wide))

for (i in 1:nrow(clean_parameters_wide)){ #for each row
  R_0_temp_no_control[i] <- gao_R_0(clean_parameters_wide_no_control[i, ]) #apply the R_0 function to that row
  R_0_temp_snail_control[i] <- gao_R_0(clean_parameters_wide_snail_control[i, ]) #apply the R_0 function to that row
}

#plot normalised R0 curves
plot(R_0_temp_no_control/max(R_0_temp_no_control) ~ clean_parameters_wide_no_control$temperature, type = 'l', main='', xlab='Temperature (°C)', ylab = latex2exp::TeX('R_{0}/max(R_{0})'), lwd = 2)
lines(R_0_temp_snail_control/max(R_0_temp_snail_control) ~ clean_parameters_wide$temperature, type = 'l', col = 'darkgreen',lty = 1, lwd = 2)
legend('topright', lwd = c(2,2), col = c('black','darkgreen'), legend = c('no control','snail control'), bty = 'n', cex = 0.8)
