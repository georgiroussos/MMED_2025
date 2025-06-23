
# Data prep ---------------------------------------------------------------

a<- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

df_C <- read.csv("/Users/georgiaroussos/Desktop/MMED_C.csv", skip =1)
df_C$Species <- "C"

df_V <- read.csv("/Users/georgiaroussos/Desktop/MMED_V.csv", skip =1)
df_V$Species <- "V"

df_comb <- rbind(df_C, df_V)


disease_params <- function(lambda = 0.9 ## transmission coefficient when prevalence is 0
                           , tau =0.1 ## 60 year natural life expectancy
){
  return(as.list(environment()))
}

y0 <- df_C$Infected[1]
init <- c(y0) 


# Catalytic model ---------------------------------------------------------

catalytic_func <- function(lambda, a, tau) { #as a function of lambda (rate of infection); 
  #a (age) and tau (delay of infection)
  
  infections = 100*(1-exp(-lambda(a-tau)))
  
  tibble(time = a, I = infections)
}


# MLE ---------------------------------------------------------------------

#probability for binomial


nllikelihood <- function(parms = disease_params(), obsDat=df_C) {
  
  p_t <- df_C$Infected / df_C$dissected
  
  nlls <- -dbinom(df_C$Infected, size = df_C$dissected, prob = p_t, log = TRUE)
  
  return(sum(nlls))
}


subsParms <- function(fit.params, fixed.params=disease_params()) {
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
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


init.pars <- c(log_lambda = log(0.9), log_tau = log(0.2))
trace <- 3 #check

optim.vals <- optim(par = init.pars
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = df_C
                    , control = list(trace = trace, maxit = 150)
                    , method = "SANN")
exp(unname(optim.vals$par))




# Simulation --------------------------------------------------------------

simEpidemic <- function(init, tseq = seq(0, 15, 1), modFunction=catalytic_func, parms = disease_params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  simDat$I <- rowSums(simDat[, Infected])
  simDat$P <- with(simDat, Infected/dissected)
  return(simDat)
}

trueParms <- disease_params() 


fitDat <- simEpidemic(init, parms = subsParms(optim.vals$par, trueParms))


