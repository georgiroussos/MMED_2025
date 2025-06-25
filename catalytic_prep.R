library(tidyverse)
library(ellipse)

# Data prep ---------------------------------------------------------------

df_C <- read.csv("MMED_C.csv", skip = 1) %>%
  mutate(Species = "C") %>%
  rename(time = category,
         N = dissected,
         I = Infected) %>%
  mutate(prop_inf = I / N)


# To produce the confidence intervals


c_binom <- mapply(function(x,n) binom.test(x,n)$conf.int, df_C$I, df_C$N)


CI_df <- data.frame(time = 0:15, Lower_CI = c(c_binom[1,]) , Upper_CI = c(c_binom[2,]) )

# merge CI with dataframe
df_C_CI <- left_join(df_C, CI_df, by = "time")

# only for first time
# write.csv(df_C_CI, "MMED_C_CI.csv")

# Catalytic model ---------------------------------------------------------

catalytic_func <- function(lambda, a, tau) {
  infections <- 1 - exp(-lambda * pmax((a - tau), 0)) #avoid negative time
  return(infections)
}

# MLE ---------------------------------------------------------------------

neg_log_likelihood <- function(params, data) {
  lambda <- params[1]
  tau <- params[2]
  
  pt <- catalytic_func(lambda, data$time, tau)
  
  # avoid log(0)
  pt <- pmin(pmax(pt, 1e-10), 1 - 1e-10)
  
  -sum(dbinom(data$I, size = data$N, prob = pt, log = TRUE))
}


init_params <- c(lambda = 0.05, tau = 5)

optim.vals <- optim(par = init_params,
                    fn = neg_log_likelihood,
                    data = df_C,
                    method = "SANN",
                    lower = c(0.0001, 0),
                    upper = c(1, max(df_C$time)),
                    # add that we want the hessian matrix
                    hessian = T)


# Results -----------------------------------------------------------------


optim.vals$par  # best estimates 

lambda_hat <- optim.vals$par[1]
tau_hat <- optim.vals$par[2]

optim.vals$par[1]
optim.vals$par[2]

#df_fit <- df_C %>%
#  mutate(predicted = catalytic_func(lambda_hat, time, tau_hat))

df_fit = data.frame(time = 0:15, predicted = catalytic_func(lambda_hat, 0:15, tau_hat))

df_fit <- left_join(df_fit, df_C_CI, by = "time")


ggplot(df_fit, aes(x = time)) +
  geom_point(aes(y = prop_inf), color = "black", size = 2) +
  geom_line(aes(y = prop_inf), color = "black") +
  geom_point(shape = 21, aes(y = predicted), color = "darkorange", fill = "white", size = 2, stroke = 1.2) +
  geom_line(aes(y = predicted), color = "darkorange", linetype = "dashed", linewidth = 1.2) +
  geom_errorbar(aes(y = prop_inf, ymin = Lower_CI, ymax = Upper_CI ))+
  labs(title = "Catalytic Model",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  theme_bw()

# Adding and plotting parameter CIs

MLEfits <- optim.vals$par
fisherInfMatrix <- solve(optim.vals$hessian)


plot(1,1, type = 'n', log = 'xy',
     ## xlim = range(alpha.seq), ylim = range(Beta.seq),
     xlim = c(0.005,0.035), ylim = c(.5,2),
     las = 1,
     xlab = expression(lambda), ylab = expression(tau),
     main = "-log(likelihood) contours", bty = "n")
## Add true parameter values to the plot #### Here
#with(true_pars, points(lambda_h, tau_h, pch = 16, cex = 2, col = 'red'))
## Add MLE to the plot
points(MLEfits['lambda'], MLEfits['tau'], pch = 16, cex = 2, col = 'black')
## Add 95% contour ellipse from Hessian
lines(ellipse(fisherInfMatrix, centre = MLEfits, level = .95))


#### Here
#legend("topright", c('truth', 'MLE', '95% Confidence Region'), lty = c(NA, NA, 1), pch = c(16,16, NA),
#col = c('red', 'black', 'black'), bg='white', bty = 'n')
legend("topright", c('MLE', '95% Confidence Region'), lty = c(NA, 1), pch = c(16, NA),
       col = c('black', 'black'), bg='white', bty = 'n')


prop_sigma<-sqrt(diag(fisherInfMatrix))
prop_sigma<-diag(fisherInfMatrix)
upper<-optim.vals$par+1.96*prop_sigma
lower<-optim.vals$par-1.96*prop_sigma
interval<-data.frame(value=optim.vals $par, upper=upper, lower=lower)