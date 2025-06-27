# Load required packages
library(ggplot2)

# Data: Age and observed proportions
age <- 0:15
y_obs <- c(
  0, 0.002218279, 0.006024096, 0.012527634,
  0.019281332, 0.027057497, 0.036262204, 0.04475043,
  0.049568966, 0.06010929, 0.074829932, 0.07826087,
  0.094972067, 0.104895105, 0.10483871, 0.117021277
)

# SEI model function (returns all timesteps)
sei_model <- function(params, time_steps = length(y_obs)) {
  lambda <- params[1]
  gamma <- params[2]
  psi <- 0.804
  
  S <- numeric(time_steps)
  E <- numeric(time_steps)
  I <- numeric(time_steps)
  
  S[1] <- 1000
  E[1] <- 0
  I[1] <- 0
  ratio_I <- numeric(time_steps)
  ratio_I[1] <- I[1] / (S[1] + E[1] + I[1])
  
  for (k in 2:time_steps) {
    S[k] <- psi * ((1 - lambda) * S[k - 1])
    E[k] <- psi * (lambda * S[k - 1] + (1 - gamma) * E[k - 1])
    I[k] <- psi * (gamma * E[k - 1] + I[k - 1])
    
    total_k <- S[k] + E[k] + I[k]
    ratio_I[k] <- I[k] / total_k
  }
  
  return(ratio_I)
}

# LSE objective function (sum of squared errors)
sei_lse <- function(params) {
  pred <- sei_model(params)
  sum((pred - y_obs)^2)
}

# Initial guess and bounds
initial_guess <- c(0.05, 0.4)
lower_bounds <- c(0.0001, 0.0001)
upper_bounds <- c(10, 10)

# Optimization using bounded L-BFGS-B method
fit <- optim(
  par = initial_guess,
  fn = sei_lse,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  hessian = TRUE
)

# Extract fitted parameters
lambda_hat <- fit$par[1]
gamma_hat <- fit$par[2]
cat(sprintf("Estimated λ = %.5f\n", lambda_hat))
cat(sprintf("Estimated γ = %.5f\n", gamma_hat))

# Simulate fitted model
I_fit <- sei_model(fit$par)

# Create data frame for plotting
df_sc <- data.frame(
  Age = age,
  Observed = y_obs,
  Fitted = I_fit
)


dataset1 <- read.csv("MMED_C_CI.csv")

dataset1 <- dataset1 %>%
  rename(Age = time)

df_Cs <- left_join(df_sc,dataset1, by = "Age" )




library(ggplot2)
plot1<-ggplot(df_Cs, aes(x = Age)) +
  geom_point(aes(y = Observed, shape = "Observed", color = "Observed"), size = 20) +
  geom_line(aes(y = Observed, color = "Observed"),linewidth=5) +
  geom_point(aes(y = Fitted, shape = "Predicted", color = "Predicted"), 
             size = 20, stroke = 1.2, fill = "white") +
  geom_line(aes(y = Fitted, color = "Predicted"), linetype = "dashed", linewidth = 5) +
  geom_errorbar(aes(y=Observed, ymin = Lower_CI, ymax = Upper_CI, color = "Observed"),width = 0.2, linewidth = 4, alpha = 0.9)+
  #geom_errorbar(aes(y = prop_inf, ymin = Lower_CI, ymax = Upper_CI, color = "Observed"), width = 0.2) +
  scale_shape_manual("Data Type", values = c("Observed" = 16, "Predicted" = 21)) +
  scale_color_manual("Data Type", values = c("Observed" = "black", "Predicted" = "red")) +
  labs(title = "T.congolense fitted SEI model (LSE)",
       x = "Ovarian Age Category",
       y = "Proportion Infected") +
  #theme(plot.title = element_text(size = 50, face = "bold", title.position="topleft"))+
  coord_cartesian(ylim=c(0,0.25))+
  theme_bw(base_size = 150) +
  theme(legend.position = c(0.02, 0.98),legend.justification = c("left", "top"),
        legend.background = element_rect(fill = "white", color = "black"), 
        legend.text = element_text(size = 75),
        legend.title = element_text(size = 80)
  ) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 21), 
                                                  fill = c(NA, "white"))))

plot1
# Save to working directory with proper plot name
ggsave("3tcongSEI_plot.png", plot = plot1, device = "png",
       width = 130, height = 112.5, units = "cm", dpi = 300,limitsize = FALSE)



