library(deSolve)
library(dplyr)
library(ggplot2)
library(ellipse)

# Import data -------------------------------------------------------------

df_obs <- read.csv("MMED_C.csv", skip = 1) %>%
  mutate(Species = "C") %>%
  rename(time = category,
         N = dissected,
         I = Infected) %>%
  mutate(prop_inf = I / N)


# Parameters --------------------------------------------------------------
disease_params <- function(lambda = 0.03,
                           gamma = 0.26,
                           mu = 0.0242,
                           phi = exp(-0.0242 * 3)) {
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
  geom_point(aes(y = prop_inf), color = "black") +
  geom_line(aes(y = predicted), color = "blue") +
  labs(title = "Fitted SEI Model",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  theme_minimal()



# Contour plots -----------------------------------------------------------

fisherInfMatrix <- solve(optim.sann$hessian)

## Initialize plot of parameters
plot(1,1, type = 'n', log = 'xy',
     ## xlim = range(alpha.seq), ylim = range(Beta.seq),
     xlim = c(0.03,4), ylim = c(.0015,0.05),
     las = 1,
     xlab = expression(lambda), ylab = expression(gamma),
     main = "-log(likelihood) contours", bty = "n")

## Add true parameter values to the plot #### Here
#with(true_pars, points(lambda, gamma, pch = 16, cex = 2, col = 'red'))
## Add MLE to the plot
points(exp(MLEfits['log_lambda']), exp(MLEfits['log_gamma']), pch = 16, cex = 2, col = 'black')
## Add 95% contour ellipse from Hessian
lines(exp(ellipse(fisherInfMatrix, centre = MLEfits, level = .95)))
#### Here
#legend("topright", c('truth', 'MLE', '95% Confidence Region'), lty = c(NA, NA, 1), pch = c(16,16, NA),
#col = c('red', 'black', 'black'), bg='white', bty = 'n')
legend("topright", c('MLE', '95% Confidence Region'), lty = c(NA, 1), pch = c(16, NA),
       col = c('black', 'black'), bg='white', bty = 'n')

