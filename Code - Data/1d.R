# 1.d - Bootstrap-Jackknife

set.seed(42) # generate reproducible results across runs

# --- IMPORTS ---
library(ggplot2)

# --- FUNCTIONS ---

bootstrap_sampling <- function(original_sample){
  n = length(original_sample)
  bootstrap_sample = sample(original_sample, size=n, replace = TRUE)
  return(bootstrap_sample)
}

jackknife_sampling <- function(original_sample, exclude_index){
  jackknife_sample = original_sample[-exclude_index] # remove item with selected index
  return(jackknife_sample)
}

# T = (X1+X2+...Xn)^2 / n
T_statistical_function <- function(sample) {
  sum = 0
  sample_size = length(sample)
  for(i in 1:sample_size){
    X_i = sample[i]
    sum = sum + X_i
  }
  T = sum^2 / sample_size
  return(T)
}

# --- MAIN ---
B = 10000 # number of bootstrap samples to generate
original_sample = readRDS("data1.rds") # read data file
n = length(original_sample) # sample length

# Bootstrap method
bootstrap_T_values = c() # save the T value for each bootstrap sample

for (i in 1:B){
  bootstrap_sample = bootstrap_sampling(original_sample) # generate bootstrap sample
  bootstrap_T_value =  T_statistical_function(bootstrap_sample) # calculate T value for each bootstrap sample
  bootstrap_T_values = append(bootstrap_T_values, bootstrap_T_value) # keep those values
}

bootstrap_se = sd(bootstrap_T_values) # bootstrap T standard error

# Jackknife method
jackknife_T_values = c() # save the T value for each jackknife sample
original_sample_T_value = T_statistical_function(original_sample)

for (i in 1:n){ # remove 1 item each time for the original sample of size n
  jackknife_sample = jackknife_sampling(original_sample, i) # remove sample i 
  jackknife_T_value = T_statistical_function(jackknife_sample) # calculate T 
  jackknife_T_values = append(jackknife_T_values, jackknife_T_value) # save them
}

jackknife_mean_value = mean(jackknife_T_values)

# calculate jackknife estimator
jackknife_estimator = n * original_sample_T_value - (n-1) * jackknife_mean_value

# standard error calculation for jackknife estimator
summ = 0
for (i in 1:n){
  summ = summ + (jackknife_T_values[i] - jackknife_mean_value)^2
}
  
jackknife_estimator_variance = ((n-1)/(n)) * summ
jackknife_estimator_se = sqrt(jackknife_estimator_variance) # jackknife T se


# --- PLOTTING ---

# 1c reimplementation

n = 80 # sample points X1, X2, ..., Xn
num_simulations = B # number of times to simulate n values and subsequently T values
T_values = c() # save the simulated T values

for (i in 1:num_simulations){
  sample = runif(n, min=0, max=1) # sample from uniform distribution n values
  T = T_statistical_function(sample) # calculate T = (X1+...+Xn)^2 / n
  T_values = append(T_values, T) # add T values to the list of T values
}
simulated_std_value = sd(T_values)

uniform_df <- data.frame(T_value = T_values)
uniform_opt_h <- 3.491 * simulated_std_value * B^(-1/3) #optimal bin width hopt=3.491*s*n^(-1/3)


# Bootstrap

bootstrap_opt_h <- 3.491 * bootstrap_se * B^(-1/3) # hopt=3.491*s*n^(-1/3)

# convert accepted samples list to df to use ggplot
bootstrap_df <- data.frame(T_value = bootstrap_T_values)

# combine two dfs in one
uniform_df$Simulation = 'Uniform' # create new column denoting to which method the T value is about
bootstrap_df$Simulation = 'Bootstrap'
df_combined = rbind(uniform_df, bootstrap_df)

# plot the bootstrap and uniform T values histoframs in a single figure
figure = ggplot(data=df_combined, aes(T_value, fill = Simulation)) + 
  geom_histogram( # histogram bootstrap
    data = bootstrap_df,
    binwidth = bootstrap_opt_h,
    kernel = 'gaussian',
    aes(y = ..density..), 
    alpha = 0.2
  ) +
  geom_histogram( # histogram uniform
    data = uniform_df,
    binwidth = bootstrap_opt_h,
    kernel = 'gaussian',
    aes(y = ..density..), 
    alpha = 0.2
  ) +
  theme(legend.position = c(0.9, 0.7)) +
  labs(title = "T values histograms",
       subtitle = "Bootstrap vs Uniform simulation, Gaussian kernel")

figure # show figure



# --- PRINTING RESULTS ---
cat(sprintf("Bootstrap T standard error: %.4f", bootstrap_se))
cat(sprintf("Jackknife T standard error: %.4f", jackknife_estimator_se))

# print(paste("Bootstrap T standard error:", bootstrap_se))
# print(paste("Jackknife T standard error:", jackknife_estimator_se))



