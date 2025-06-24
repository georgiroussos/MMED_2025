library(tidyverse)

# Data prep ---------------------------------------------------------------

df_C <- read.csv("/Users/georgiaroussos/Desktop/MMED_C.csv", skip = 1) %>%
  mutate(Species = "C") %>%
  rename(time = category,
         N = dissected,
         I = Infected) %>%
  mutate(prop_inf = I / N)

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
                    upper = c(1, max(df_C$time)))


# Results -----------------------------------------------------------------


optim.vals$par  # best estimates 

lambda_hat <- optim.vals$par[1]
tau_hat <- optim.vals$par[2]

optim.vals$par[1]
optim.vals$par[2]

df_fit <- df_C %>%
  mutate(predicted = catalytic_func(lambda_hat, time, tau_hat))


ggplot(df_fit, aes(x = time)) +
  geom_point(aes(y = prop_inf), color = "black") +
  geom_line(aes(y = predicted), color = "blue") +
  labs(title = "Catalytic Model",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  theme_minimal()

