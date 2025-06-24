library(tidyverse)

# Data prep ---------------------------------------------------------------

df_C <- read.csv("/Users/georgiaroussos/Desktop/MMED_C.csv", skip =1)
df_C$Species <- "C"

df_V <- read.csv("/Users/georgiaroussos/Desktop/MMED_V.csv", skip =1)
df_V$Species <- "V"

df_comb <- rbind(df_C, df_V)


df_C <- df_C %>%
  rename(time = category,
         N = dissected,
         I = Infected) %>%
  mutate(prop_inf = I / N)


# Disease parameters ------------------------------------------------------


disease_params <- function(lambda = 0.9 
                           , tau =0.1 
){
  return(as.list(environment()))
}


# Catalytic model ---------------------------------------------------------

catalytic_func <- function(lambda, a, tau) { #as a function of lambda (rate of infection); 
  #a (age) and tau (delay of infection)
  infections <- 1 - exp(-lambda * (a - tau)) #forgot multiplication *
  tibble(time = a, prevalence = infections)
}

# MLE ---------------------------------------------------------------------

#probability for binomial

nllikelihood <- function(parms = disease_params(), obsDat = df_C) {
  lambda <- parms$lambda
  tau <- parms$tau
  
  #pt <- obsDat$I/obsDat$N
  pt <- 1 - exp(-lambda * pmax((obsDat$time - tau), 0))  # avoids negative time # without this specification causes issues
 #catalytic models for infection assumes probability of infection is P(a)=1-exp(-lambda*(a-t))
  
  nlls <- -dbinom(obsDat$I, size = obsDat$N, prob = pt, log = TRUE)
  
  return(sum(nlls))
}


#substitution for log-transformed parameters
subsParms <- function(fit.params, fixed.params = disease_params()) {
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
    for (nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for (nm in loggedParms) assign(gsub('log_', '', nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })
}

guess.params <- c(log_lambda = log(0.9), log_tau = log(0.2)) #FIXME 
subsParms(guess.params, disease_params())


## Make likelihood a function of fixed and fitted parameters.
objFXN <- function(fit.params ## paramters to fit
                   , fixed.params =disease_params() ## fixed paramters
                   , obsDat=df_C) {
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat) ## then call likelihood
}

objFXN(guess.params, disease_params())



# Optimisation ------------------------------------------------------------

trace <- 3 #check

init.pars <- c(log_lambda = log(0.9), log_tau = log(0.2))

optim.vals <- optim(
  par = init.pars,
  fn = objFXN,
  fixed.params = disease_params(),
  obsDat = df_C,
  control = list(trace = trace, maxit = 150),
  method = "SANN"
)

best_pars <- exp(unname(optim.vals$par))
print(paste("Estimated lambda:", round(best_pars[1], 4)))
print(paste("Estimated tau:", round(best_pars[2], 4)))

best_pars 

# Simulation --------------------------------------------------------------


fitted_values <- catalytic_func(lambda = best_pars[1], a = df_C$time, tau = best_pars[2]) %>%
  mutate(N = df_C$N,
         observed_infections = df_C$I,
         fitted_infections = prevalence * N)


ggplot(fitted_values, aes(x = time)) +
  geom_point(aes(y = observed_infections), color = "black", size = 2) +
  geom_line(aes(y = fitted_infections), color = "blue", linewidth = 1) +
  labs(title = "Catalytic Model",
       x = "Age",
       y = "Number Infected") +
  theme_minimal()
