library(ggplot2)
library(scales)
library(tidyverse)
library(ellipse)

# Data prep ---------------------------------------------------------------

df_C <- read.csv("C:\\Users\\davry\\OneDrive\\Desktop\\mmed\\MMED Project\\MMED_C.csv", skip = 1) %>%
mutate(Species = "C") %>%
 rename(time = category,
       N = dissected,
      I = Infected) %>%
 mutate(prop_inf = I / N)


# To produce the confidence intervals


# c_binom <- mapply(function(x,n) binom.test(x,n)$conf.int, df_C$I, df_C$N)


# CI_df <- data.frame(time = 0:15, Lower_CI = c(c_binom[1,]) , Upper_CI = c(c_binom[2,]) )

# merge CI with dataframe
# df_C_CI <- left_join(df_C, CI_df, by = "time")

# only for first time
# write.csv(df_C_CI, "MMED_C_CI.csv")

# import data again
df_C <- read.csv("C:\\Users\\davry\\OneDrive\\Desktop\\mmed\\MMED Project\\MMED_C_CI.csv")


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

optim.vals_C <- optim(par = init_params,
                    fn = neg_log_likelihood,
                    data = df_C,
                    method = "SANN",
                    lower = c(0.0001, 0),
                    upper = c(1, max(df_C$time)),
                    # add that we want the hessian matrix
                    hessian = T)


# Results -----------------------------------------------------------------


optim.vals_C$par  # best estimates 

lambda_hat <- optim.vals_C$par[1]
tau_hat <- optim.vals_C$par[2]

optim.vals_C$par[1]
optim.vals_C$par[2]

#df_fit <- df_C %>%
#  mutate(predicted = catalytic_func(lambda_hat, time, tau_hat))

df_fit = data.frame(time = 0:15, predicted = catalytic_func(lambda_hat, 0:15, tau_hat))

df_fit <- left_join(df_fit, df_C, by = "time")


T.congolense_fitted_catalytic_model_MLE <- ggplot(df_fit, aes(x = time)) +
  geom_point(aes(y = prop_inf, shape = "Observed", color = "Observed"), size = 20) +
  geom_line(aes(y = prop_inf, color = "Observed"), linewidth = 5) +
  geom_point(aes(y = predicted, shape = "Predicted", color = "Predicted"), 
             size = 20, stroke = 1.2, fill = "white") +
  geom_line(aes(y = predicted, color = "Predicted"), linetype = "dashed", linewidth = 5) +
  geom_errorbar(aes(y = prop_inf, ymin = Lower_CI, ymax = Upper_CI, color = "Observed"), width = 0.2, linewidth = 4, alpha = 0.9) +
  scale_shape_manual("Data Type", values = c("Observed" = 16, "Predicted" = 21)) +
  scale_color_manual("Data Type", values = c("Observed" = "black", "Predicted" = "red")) +
  labs(title = "T.congolense fitted Catalytic model \n (MLE)",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  coord_cartesian(ylim = c(0, 0.25)) + 
  theme_bw(base_size = 150) + theme(legend.position = c(0.02, 0.98),legend.justification = c("left", "top"),legend.background = element_rect(fill = "white", color = "black"), legend.text = element_text(size = 75),
                                    legend.title = element_text(size = 80),) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 21), 
                                                  fill = c(NA, "white"))))

ggsave("C:\\Users\\davry\\OneDrive\\Desktop\\mmed\\MMED Project\\graphs final\\T.congolense_fitted_catalytic_model_MLE.png", 
       T.congolense_fitted_catalytic_model_MLE, 
       device = "png", 
       width = 130, height = 112.5, units = "cm" , dpi = 300, 
       limitsize = FALSE)


# Adding and plotting parameter CIs

MLEfits_C <- optim.vals_C$par
fisherInfMatrix_C <- solve(optim.vals_C$hessian)

png("C:\\Users\\davry\\OneDrive\\Desktop\\mmed\\MMED Project\\graphs final\\MLE_ConfidenceRegion.png", 
    width = 130, height = 112.5, units = "cm", res = 300)

# Adjust graphical parameters
par(
  cex.lab = 5,          
  cex.axis = 4.5,       
  cex.main = 4,         
  mar = c(8, 8, 5, 3),  
  mgp = c(5, 1.5, 0),   
  tcl = -0.5,           
  las = 1               
)

# Plot setup
plot(1, 1, type = 'n', log = 'xy',
     xlim = c(0.005, 0.035), 
     ylim = c(0.1, 2),
     xlab = expression(lambda~"(/ nine days)"), 
     ylab = expression(tau~"(/ nine days)"),
     main = "MLE and 95% Confidence Region for Model Parameters\nfor Catalytic model",
     bty = "n")

# Add MLE points and confidence ellipses
points(MLEfits_C['lambda'], MLEfits_C['tau'], pch = 16, cex = 4, col = 'blue') 
lines(ellipse(fisherInfMatrix_C, centre = MLEfits_C, level = 0.95), col = 'blue', lwd = 5) 

points(MLEfits_V['lambda'], MLEfits_V['tau'], pch = 16, cex = 4, col = 'red') 
lines(ellipse(fisherInfMatrix_V, centre = MLEfits_V, level = 0.95), col = 'red', lwd = 5) 

# Add legend
legend("topright", 
       legend = c("MLE (T. congolense)", "95% CI (T. congolense)", 
                  "MLE (T. vivax)", "95% CI (T. vivax)"),
       col = c("blue", "blue", "red", "red"),
       pch = c(16, NA, 16, NA),
       lty = c(NA, 1, NA, 1),
       lwd = 5,
       pt.cex = 4,
       bty = "n",
       cex = 4)

dev.off()


prop_sigma<-sqrt(diag(fisherInfMatrix))
prop_sigma<-diag(fisherInfMatrix)
upper<-optim.vals$par+1.96*prop_sigma
lower<-optim.vals$par-1.96*prop_sigma
interval<-data.frame(value=optim.vals $par, upper=upper, lower=lower)
