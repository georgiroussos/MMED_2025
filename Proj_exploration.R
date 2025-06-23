
df_C <- read.csv("/Users/georgiaroussos/Desktop/MMED_C.csv", skip = 1) %>% 
  mutate(prop_infected = Infected/dissected, 
         species = "Congolense")

df_V <- read.csv("/Users/georgiaroussos/Desktop/MMED_V.csv", skip = 1) %>% 
  mutate(prop_infected = Infected/dissected, 
         species = "Vivax")

df_comb <- rbind(df_C, df_V)

ggplot() +
  geom_line(data = df_comb, aes(x = category, y = prop_infected, group = species),
            color = "black", alpha = 0.2) +
  geom_point(data = df_comb, aes(x = category, y = prop_infected, group = species, color = species),
             size = 2,  fshape = 21) +
  #scale_color_brewer(palette = "Paired") +   # Color scale for lines, points, and error bars
  #scale_fill_brewer(palette = "Paired") + 
  labs(title = "Proportion of Tstetse Flies Infected",
       x = "Age Category (t)",
       y = "Proportion Infected", 
       color = "Species") +
  theme_bw()

