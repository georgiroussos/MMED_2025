library(deSolve)
library(dplyr)
library(ggplot2)
library(ellipse)


# Import data -------------------------------------------------------------

#df_C <- read.csv("MMED_C.csv", skip = 1) %>%
# mutate(Species = "C") %>%
#  rename(time = category,
#        N = dissected,
#       I = Infected) %>%
#  mutate(prop_inf = I / N)

df_obs <- read.csv("C:\\Users\\davry\\OneDrive\\Desktop\\mmed\\MMED Project\\MMED_V_CI.csv") 


# Parameters --------------------------------------------------------------
disease_params <- function(lambda = 0.03,
                           gamma = 0.26,
                           mu = 0.0242,
                           phi = exp(-0.0242 * 9)) {
  return(as.list(environment()))
}

init <- c(S = 1000, E = 0, I = 0)
tseqMonth <- 0:15

# SEI model ---------------------------------------------------------------

SImod <- function(tt, yy, parms) with(c(parms, as.list(yy)), {
  N <- S + E + I
  deriv <- rep(NA, 3)
  deriv[1] <- -phi * lambda * S
  deriv[2] <- phi * (lambda * S - gamma * E)
  deriv[3] <- phi * gamma * E
  return(list(deriv))
})

#SEI model as defined by John 


# Simulation --------------------------------------------------------------

simEpidemic <- function(init, tseq = tseqMonth, modFunction = SImod, parms = disease_params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms = parms))
  return(simDat)
}



# Negative log-likelihood ------------------------------------------------
nllikelihood <- function(parms, obsDat = df_obs) {
  simDat <- simEpidemic(init, tseq = obsDat$time, parms = parms)
  pt <- simDat$I / (simDat$S + simDat$E + simDat$I + 1e-10) #avoid 0 denominator
  pt <- pmin(pmax(pt, 1e-10), 1 - 1e-10)  # restrict probabilities
  nlls <- -dbinom(obsDat$I, size = obsDat$N, prob = pt, log = TRUE)
  return(sum(nlls))
}
#before: nllikelihood() just used obs$I/obs$N (not model output!)



# Parameter substitution function -----------------------------------------
subsParms <- function(fit.params, fixed.params = disease_params()) {
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_', '', nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })
}


# Object function ---------------------------------------------------------

objFXN <- function(fit.params, fixed.params = disease_params(), obsDat = df_obs) {
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat)
}



# Initial parameter guesses ----------------------------------------------
guess.params <- c(log_lambda = log(0.05), log_gamma = log(0.3))


# Optimisations -----------------------------------------------------------

optim.sann_V <- optim(par = guess.params,
                    fn = objFXN,
                    method = "SANN",
                    control = list(trace = 3, maxit = 150),
                    fixed.params = disease_params(),
                    obsDat = df_obs,
                    # add hessian = T
                    hessian = T)

optim.final_V <- optim(par = optim.sann_V$par,
                     fn = objFXN,
                     method = "Nelder-Mead",
                     control = list(trace = 3, maxit = 800, reltol = 1e-7),
                     hessian = TRUE,
                     fixed.params = disease_params(),
                     obsDat = df_obs)

MLEfits_V <- optim.final_V$par
exp(unname(MLEfits_V))  

# Plot --------------------------------------------------------------------

fitted_df <- simEpidemic(init, tseq = tseqMonth, parms = subsParms(MLEfits_V, disease_params())) %>%
  mutate(predicted = I / (S + E + I)) %>%
  select(time, predicted)

df_plot <- df_obs %>%
  mutate(prop_inf = I / N) %>%
  left_join(fitted_df, by = "time")  # merge on time

T.vivax_fitted_SEI_model_MLE <- ggplot(df_plot, aes(x = time)) +
  geom_point(aes(y = prop_inf, shape = "Observed", color = "Observed"), size = 20) +
  geom_line(aes(y = prop_inf, color = "Observed"), linewidth = 5) +
  geom_point(aes(y = predicted, shape = "Predicted", color = "Predicted"), 
             size = 20, stroke = 1.2, fill = "white") +
  geom_line(aes(y = predicted, color = "Predicted"), linetype = "dashed", linewidth = 5) +
  geom_errorbar(aes(y = prop_inf, ymin = Lower_CI, ymax = Upper_CI, color = "Observed"), width = 0.2, linewidth = 4, alpha = 0.9) +
  scale_shape_manual("Data Type", values = c("Observed" = 16, "Predicted" = 21)) +
  scale_color_manual("Data Type", values = c("Observed" = "black", "Predicted" = "blue")) +
  labs(title = "T.vivax fitted SEI model (MLE)",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  coord_cartesian(ylim = c(0, 0.25)) + 
  theme_bw(base_size = 150) + theme(legend.position = c(0.02, 0.98),legend.justification = c("left", "top"),legend.background = element_rect(fill = "white", color = "black"), legend.text = element_text(size = 75),
                                    legend.title = element_text(size = 80),) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 21), 
                                                  fill = c(NA, "white"))))

ggsave("C:\\Users\\davry\\OneDrive\\Desktop\\mmed\\MMED Project\\graphs final\\T.vivax_fitted_SEI_model_MLE.png", 
       T.vivax_fitted_SEI_model_MLE, 
       device = "png", 
       width = 130, height = 112.5, units = "cm", dpi = 300,
       limitsize = FALSE)

# Contour plots -----------------------------------------------------------

fisherInfMatrix_V <- solve(optim.final_V$hessian)  

ell_V <- ellipse(fisherInfMatrix_V, centre = MLEfits_V, level = 0.95)
ell_exp_V <- exp(ell_V)
print(ell_exp_V)

lambda_hat_V <- exp(MLEfits_V['log_lambda'])
gamma_hat_V <- exp(MLEfits_V['log_gamma'])
gamma_hat_V
lambda_hat_V

png("C:\\Users\\davry\\OneDrive\\Desktop\\mmed\\MMED Project\\graphs final\\MLE_ConfidenceRegion_SEI_test.png", 
    width = 150, height = 120.5, units = "cm", res = 300)

par(cex.lab = 10,    
    cex.axis = 7,  
    cex.main = 8,
    mar = c(25, 25, 25, 20),  
    mgp = c(15, 2.5, 0))  

MLE_ConfidenceRegion_SEI_test <- plot(1, 1, type = 'n', log = 'xy',
     xlim = c(0.005,0.035), ylim = c(0.05,5),
     las = 1,
     xlab = expression(lambda~"(/ nine days)"), ylab = expression(gamma~"(/ nine days)"),
     main = "MLE and 95% Confidence Region for Model Parameters\nfor SEI ", bty = "n", cex = 20)

#Vivax Plot
points(lambda_hat_V, gamma_hat_V, pch = 16, cex = 4, col = 'blue')
lines(ell_exp_V, col = 'blue', lwd = 5)

#Congolense plot
points(lambda_hat_C, gamma_hat_C, pch = 16, cex = 4, col = 'red')
lines(ell_exp_C, col = 'red', lwd = 5)

#Add legend 
legend("topright", 
       legend = c("MLE (T. congolense)", 
                  "95% CI (T. congolense)", 
                  "MLE (T. vivax)", 
                  "95% CI (T. vivax)"),
       col = c("red", "red", "blue", "blue"),
       pch = c(16, NA, 16, NA),   
       lty = c(NA, 1, NA, 1),     
       lwd = 5,
       pt.cex = 5,
       bty = "n",
       cex = 7 )

dev.off()

# Parameter CIs -----------------------------------------------------------

fisherInfMatrix_V <- solve(optim.final_V$hessian)
log_se <- sqrt(diag(fisherInfMatrix))


log_lambda_hat <- MLEfits["log_lambda"]
log_gamma_hat <- MLEfits["log_gamma"]

ci_log_lambda <- c(log_lambda_hat - 1.96 * log_se["log_lambda"],
                   log_lambda_hat + 1.96 * log_se["log_lambda"])

ci_log_gamma <- c(log_gamma_hat - 1.96 * log_se["log_gamma"],
                  log_gamma_hat + 1.96 * log_se["log_gamma"])

ci_lambda <- exp(ci_log_lambda)
ci_gamma <- exp(ci_log_gamma)

lambda_hat <- exp(log_lambda_hat)
gamma_hat <- exp(log_gamma_hat)
lambda_hat 
gamma_hat 

cat("lambda_hat:", lambda_hat, "\n")
cat("95% CI for lambda:", round(ci_lambda, 5), "\n\n")

cat("gamma_hat:", gamma_hat, "\n")
cat("95% CI for gamma:", round(ci_gamma, 5), "\n")
