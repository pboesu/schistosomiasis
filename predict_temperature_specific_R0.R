# Example code to predict temperature-specific R0 curves for schistosomiasis with and without snail control
# Supplement to "Interventions can shift the thermal optimum for parasitic disease transmission" by K.H. Nguyen, P.H. Boersch-Supan, R.B. Hartman, S.Y. Mendiola, V.J. Harwood, D.J. Civitello, J.R. Rohr (2021)
# https://github.com/pboesu/schistosomiasis

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

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

#read in predicted thermal performance curves
clean_parameters_long <- read_csv('data/trait_predictions.csv', comment = "#")

#plot predicted thermal performance curves (simplified version of Figure 2)
ggplot(clean_parameters_long, aes(x = temperature, y = lci)) +
  geom_line(lty=3) +
  geom_line(aes(x = temperature, y = uci), lty=3) +
  facet_wrap(~ variable, scales = 'free_y') +
  geom_line(aes(x = temperature, y = value), col = 'red') +
  ylab("1/day") +
  theme_classic()

#reshape data to wide format and add temperature-invariant parameters
clean_parameters_wide <- clean_parameters_long %>%
  select(-uci, -lci, -parameter_type, -units) %>% 
  tidyr::pivot_wider(values_from = value, names_from = variable) %>%
  mutate(beta_1 = 4.06e-09, #Temperature-invariant parameters from Gao et al. 2010 Table 1
         delta_1 = 0.0039,
         k = 300,
         Lambda_1 = 8000,
         M_0 = 1e+06,
         mu_1 = 3.84e-05)

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

#plot normalised R0 curves (simplified version of Figure 3A)
plot(R_0_temp_no_control/max(R_0_temp_no_control) ~ clean_parameters_wide_no_control$temperature, type = 'l', main='', xlab='Temperature (°C)', ylab = latex2exp::TeX('R_{0}/max(R_{0})'), lwd = 2)
lines(R_0_temp_snail_control/max(R_0_temp_snail_control) ~ clean_parameters_wide$temperature, type = 'l', col = 'darkgreen',lty = 1, lwd = 2)
legend('topright', lwd = c(2,2), col = c('black','darkgreen'), legend = c('no control','snail control'), bty = 'n', cex = 0.8)
