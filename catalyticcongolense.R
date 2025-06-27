# Load libraries
library(ggplot2)

# Data
age <- 0:15
y_obs <- c(
  0, 0.002218279, 0.006024096, 0.012527634,
  0.019281332, 0.027057497, 0.036262204, 0.04475043,
  0.049568966, 0.06010929, 0.074829932, 0.07826087,
  0.094972067, 0.104895105, 0.10483871, 0.117021277
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
df_cc <- data.frame(
  Age = age,
  Observed = y_obs,
  Fitted = y_fit)
  
library(tidyverse)

dataset3 <- read.csv("MMED_C_CI.csv")

dataset3 <- dataset3 %>%
  rename(Age = time)

df_Cc <- left_join(df_cc,dataset3, by = "Age" )




library(ggplot2)
plot3<-ggplot(df_Cc, aes(x = Age)) +
  geom_point(aes(y = Observed, shape = "Observed", color = "Observed"), size = 20) +
  geom_line(aes(y = Observed, color = "Observed"),linewidth=5) +
  geom_point(aes(y = Fitted, shape = "Predicted", color = "Predicted"), 
             size = 20, stroke = 1.2, fill = "white") +
  geom_line(aes(y = Fitted, color = "Predicted"), linetype = "dashed", linewidth = 5) +
  geom_errorbar(aes(y=Observed, ymin = Lower_CI, ymax = Upper_CI, color = "Observed"),width = 0.2, linewidth = 4, alpha = 0.9)+
  #geom_errorbar(aes(y = prop_inf, ymin = Lower_CI, ymax = Upper_CI, color = "Observed"), width = 0.2) +
  scale_shape_manual("Data Type", values = c("Observed" = 16, "Predicted" = 21)) +
  scale_color_manual("Data Type", values = c("Observed" = "black", "Predicted" = "red")) +
  labs(title = "T.congolense fitted catalytic model (LSE)",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  #theme(plot.title = element_text(size = 50, face = "bold", title.position="topleft"))+
  coord_cartesian(ylim=c(0,0.25))+
  theme_bw(base_size = 130) +
  theme(legend.position = c(0.02, 0.98),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = "white", color = "black"), 
        legend.text = element_text(size = 75),
        legend.title = element_text(size = 80)
      ) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 21), 
                                                  fill = c(NA, "white"))))

plot3
# Save to working directory with proper plot name
ggsave("3tcongcat_plot.png", plot = plot3, device = "png",
       width = 130, height = 112.5, units = "cm", dpi = 300,limitsize = FALSE)




























 
  
 