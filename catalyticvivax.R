# Load libraries
library(ggplot2)

# Data
age <- 0:15
y_obs <- c(
  0.0003711952487, 0.004826678368, 0.01164079823, 0.02031144211,
  0.02876480541, 0.03907074974, 0.04829545455, 0.05719237435,
  0.0670995671, 0.07714285714, 0.0830449827, 0.09251101322,
  0.09947643979, 0.1095890411, 0.1176470588, 0.1290322581
)

# Catalytic model function
catalytic_model <- function(pars, age) {
  lambda <- pars[1]
  tau <- pars[2]
  
  y_pred <- ifelse(age >= tau,
                   1 - exp(-lambda * (age - tau)),
                   0.0)
  return(y_pred)
}

# Residual sum of squares function
residuals_fn <- function(pars) {
  y_pred <- catalytic_model(pars, age)
  sum((y_pred - y_obs)^2)
}

# Initial guess and bounds
initial_guess <- c(lambda = 0.1, tau = 3)
lower_bounds <- c(0.0001, 1)
upper_bounds <- c(0.1, 50)

# Least squares optimization (bounded)
fit <- optim(
  par = initial_guess,
  fn = residuals_fn,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  hessian = TRUE
)

# Extract fitted values
lambda_fit <- fit$par["lambda"]
tau_fit <- fit$par["tau"]

cat(sprintf("Estimated λ = %.5f\n", lambda_fit))
cat(sprintf("Estimated τ = %.2f days\n", tau_fit))

# Generate predicted values
y_fit <- catalytic_model(fit$par, age)



# Create dataframe for plotting
df_vc <- data.frame(
  Age = age,
  Observed = y_obs,
  Fitted = y_fit)

library(tidyverse)

dataset4 <- read.csv("MMED_V_CI.csv")

dataset4 <- dataset4 %>%
  rename(Age = time)

df_Vc <- left_join(df_vc,dataset2, by = "Age" )


library(ggplot2)
plot2<-ggplot(df_Vc, aes(x = Age)) +
  geom_point(aes(y = Observed, shape = "Observed", color = "Observed"), size = 20) +
  geom_line(aes(y = Observed, color = "Observed"),linewidth=5) +
  geom_point(aes(y = Fitted, shape = "Predicted", color = "Predicted"), 
             size = 20, stroke = 1.2, fill = "white") +
  geom_line(aes(y = Fitted, color = "Predicted"), linetype = "dashed", linewidth = 5) +
  geom_errorbar(aes(y=Observed, ymin = Lower_CI, ymax = Upper_CI, color = "Observed"),width = 0.2, linewidth = 4, alpha = 0.9)+
  #geom_errorbar(aes(y = prop_inf, ymin = Lower_CI, ymax = Upper_CI, color = "Observed"), width = 0.2) +
  scale_shape_manual("Data Type", values = c("Observed" = 16, "Predicted" = 21)) +
  scale_color_manual("Data Type", values = c("Observed" = "black", "Predicted" = "blue")) +
  labs(title = "T.vivax fitted catalytic model (LSE)",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  theme(plot.title = element_text(size = 75, face = "bold"))+
  coord_cartesian(ylim=c(0,0.25))+
  theme_bw(base_size = 150) +
  theme(legend.position = c(0.02, 0.98),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = "white", color = "black"), 
        legend.text = element_text(size = 75),
  legend.title = element_text(size = 80)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 21), 
                                                  fill = c(NA, "white"))))

plot2 <- plot2 +
  theme(text = element_text(size = 80)) 
plot2
# Save to working directory with proper plot name
ggsave("Ttvivaxcat_plot.png", plot = plot2, device = "png",
       width = 130, height = 112.5, units = "cm", dpi = 300,limitsize = FALSE)











ggplot(df_vc, aes(x = Age)) +
  geom_point(aes(y = Observed, shape = "Observed", color = "Observed"), size = 2.5) +
  geom_line(aes(y = Observed, color = "Observed")) +
  geom_point(aes(y = Fitted, shape = "Predicted", color = "Predicted"), 
             size = 2.5, stroke = 1.2, fill = "white") +
  geom_line(aes(y = Fitted, color = "Predicted"), linetype = "dashed", linewidth = 1.5) +
  #geom_errorbar(aes(y = prop_inf, ymin = Lower_CI, ymax = Upper_CI, color = "Observed"), width = 0.2) +
  scale_shape_manual("Data Type", values = c("Observed" = 16, "Predicted" = 21)) +
  scale_color_manual("Data Type", values = c("Observed" = "black", "Predicted" = "red")) +
  labs(title = "T.vivax fitted Catalytic model (LSE)",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = c(16, 21), 
                                                  fill = c(NA, "white"))))



