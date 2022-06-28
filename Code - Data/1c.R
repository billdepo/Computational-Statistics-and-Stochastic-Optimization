# 1.c - Simulated samples

set.seed(42) # generate reproducible results across runs

# --- IMPORTS ---
library(ggplot2)

# --- FUNCTIONS ---

# T = (X1+X2+...Xn)^2 / n
T_statistical_function <- function(sample) {
  sample_size = length(sample)
  summ = sum(sample)
  T = summ^2 / sample_size
  return(T)
}

# --- MAIN ---
n = 80 # sample points X1, X2, ..., Xn
num_simulations = 10000 # number of times to simulate n values and subsequently T values
T_values = c() # save the simulated T values

for (i in 1:num_simulations){
  sample = runif(n, min=0, max=1) # sample from uniform distribution n values
  T = T_statistical_function(sample) # calculate T = (X1+...+Xn)^2 / n
  T_values = append(T_values, T) # add T values to the list of T values
}

simulated_mean_value = mean(T_values)
simulated_std_value = sd(T_values)
real_mean_value = 1/3 + (n-1)/4 # analytical calculation of E(T)


# --- PLOTTING ---

# convert T values list to df to use ggplot
df <- data.frame(T = T_values)

#optimal bin width hopt=3.491*s*n^(-1/3)
opt_h <- 3.491 * simulated_std_value * num_simulations^(-1/3) 

figure <- ggplot(data=df, aes(x=T)) +
  geom_histogram(
    binwidth = opt_h, # bin width = optimal
    aes(y = ..density..), 
    color = "white", # histogram bins borders
    fill = "indianred" # histogram bins fill color
  ) +
  geom_density(size = 1.1, 
               color = "red") + # simulated T values dist. curve
  labs(title = "Simulated T values distribution histogram",
       subtitle = 'Gaussian kernel') # title

figure # show figure


# --- PRINTING RESULTS ---

# print(paste("Simulated mean of T:", simulated_mean_value))
# print(paste("Real mean of T:", real_mean_value))
# print(paste("Deviation of simulated mean from theoretical:",
#             round(100*abs((simulated_mean_value-real_mean_value)/real_mean_value), 2), "%"))
# print(paste("Simulated standard deviation of T:", simulated_std_value))


cat(sprintf("Simulated mean of T: %.4f
Real mean of T: %.4f
Deviation of simulated mean from theoretical: %.2f %%
Simulated standard deviation of T: %.4f",
simulated_mean_value,
real_mean_value,
round(100*abs((simulated_mean_value-real_mean_value)/real_mean_value), 2),
simulated_std_value))
