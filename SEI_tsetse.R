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

df_obs <- read.csv("MMED_C_CI.csv") 


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

optim.sann <- optim(par = guess.params,
                    fn = objFXN,
                    method = "SANN",
                    control = list(trace = 3, maxit = 150),
                    fixed.params = disease_params(),
                    obsDat = df_obs,
                    # add hessian = T
                    hessian = T)

optim.final <- optim(par = optim.sann$par,
                     fn = objFXN,
                     method = "Nelder-Mead",
                     control = list(trace = 3, maxit = 800, reltol = 1e-7),
                     hessian = TRUE,
                     fixed.params = disease_params(),
                     obsDat = df_obs)

MLEfits <- optim.final$par
exp(unname(MLEfits))  

# Plot --------------------------------------------------------------------

fitted_df <- simEpidemic(init, tseq = tseqMonth, parms = subsParms(MLEfits, disease_params())) %>%
  mutate(predicted = I / (S + E + I)) %>%
  select(time, predicted)

df_plot <- df_obs %>%
  mutate(prop_inf = I / N) %>%
  left_join(fitted_df, by = "time")  # merge on time

ggplot(df_plot, aes(x = time)) +
  geom_point(aes(y = prop_inf), color = "black", size = 2) +
  geom_line(aes(y = prop_inf), color = "black") +
  geom_point(shape = 21, aes(y = predicted), color = "darkorange", fill = "white", size = 2, stroke = 1.2) +
  geom_line(aes(y = predicted), color = "darkorange", linetype = "dashed", linewidth = 1.2) +
  geom_errorbar(aes(y = prop_inf, ymin = Lower_CI, ymax = Upper_CI ))+
  labs(title = "Fitted SEI Model",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  theme_bw()

# Contour plots -----------------------------------------------------------

fisherInfMatrix <- solve(optim.final$hessian)  

ell <- ellipse(fisherInfMatrix, centre = MLEfits, level = 0.95)
ell_exp <- exp(ell)
ell_exp

lambda_hat <- exp(MLEfits['log_lambda'])
gamma_hat <- exp(MLEfits['log_gamma'])
gamma_hat
lambda_hat

plot(1, 1, type = 'n', log = 'xy',
     xlim = c(0.005,0.035), ylim = c(0.05,2.5),
     las = 1,
     xlab = expression(lambda), ylab = expression(gamma),
     main = "-log(likelihood) contours", bty = "n")

points(lambda_hat, gamma_hat, pch = 16, cex = 2, col = 'black')
lines(ell_exp)
legend("topright", c('MLE', '95% Confidence Region'), lty = c(NA, 1), pch = c(16, NA),
       col = c('black', 'black'), bg = 'white', bty = 'n')

# Parameter CIs -----------------------------------------------------------

fisherInfMatrix <- solve(optim.final$hessian)
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


# AIC



log_likelihood <- -nllikelihood(parms = disease_params(lambda = lambda_hat  , gamma = gamma_hat), obsDat = df_obs)

num_pars <- 2
traditional_AIC <- -2 *log_likelihood + 2*(num_pars)
traditional_AIC
