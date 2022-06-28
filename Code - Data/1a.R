# 1.a - Rejection Sampling
# Sample from N(0,1) using utilizing a Cauchy proposal distribution

set.seed(42) # generate reproducible results across runs


# --- IMPORTS ---
library(ggplot2)


# --- FUNCTIONS ---

# returns Cauchy sample derived from transformed inverse sample
cauchy_sampling <- function(uniform_sample) {
  transformed_sample <- tan(pi * (uniform_sample-1/2))
  return(transformed_sample)
}

# f(x)
f <- function(x) {
  return((1/sqrt(2*pi))*exp(-x^2/2))
}

# g(x)
g <- function(x) {
  return(1/(pi*(1+x^2)))
}

# acceptance probability 
acceptance_prob <- function(M, sample) {
  probability <- f(sample)/(M*g(sample)) # acceptance probability
  return(probability)
}

# --- MAIN ---
M = sqrt(2*pi/exp(1)) # optimal M constant multiplier
num_total_samples = 0 # counter of total samples generated
num_accepted_samples = 0 # counter of samples accepted (<= total)
required_samples = 1000 # number of accepted samples to be generated
accepted_samples = c() # save accepted samples' values

while (num_accepted_samples < required_samples){
  num_total_samples = num_total_samples + 1 # calculate all samples generated
  uniform_sample = runif(1, min=0, max=1) # sample from uniform distribution
  cauchy_sample = cauchy_sampling(uniform_sample) # transform to get cauchy sample
  accept_probability = acceptance_prob(M, cauchy_sample) # acceptance prob
  uniform_probability = runif(1, min=0, max=1) # gerete number in [0,1]
  if (accept_probability >= uniform_probability){
    num_accepted_samples  = num_accepted_samples + 1
    accepted_samples = append(accepted_samples, cauchy_sample)
  }
}


# --- PLOTTING ---

opt_h <- 3.491 * 1 * required_samples^(-1/3) # Normal dist, optimal bin width

# convert accepted samples list to df to use ggplot
df <- data.frame(x = accepted_samples)

figure <- ggplot(data=df, aes(x)) +
    geom_histogram(
      binwidth = opt_h, # bin width = optimal
      aes(y = ..density..),
      color = "white", # histogram bins borders
      fill = "indianred" # histogram bins fill color
    ) +
    geom_density(size=1.1,  # plot simulated distribution curve
                 kernel="gaussian",
                 mapping = aes(color="Simulated"),
                 key_glyph = draw_key_path) +
    stat_function(fun = f, # plot N(0,1) pdf curve  
                  size=1.1, 
                  mapping = aes(color="Standard normal"),
                  key_glyph = draw_key_path) + 
    scale_color_manual(name = "Pdf:",
                     values = c("red", "darkblue"), # Color specification
                     labels = c("Simulated - Gaussian kernel", "Standard normal")) +
    theme(legend.position = c(0.8, 0.7)) +
    labs(title = "Simulated standard normal values histogram",
         subtitle = "Method: Rejection sampling") # title

figure # show figure


# show h(x)=f(x)/g(x) that is used to select the best M
x = seq(-5, 5, 0.05) # create a sequence from -5 to 5 with step 0.05
h <- function(x) f(x)/g(x) # define h function
h_values = h(x) # caluclate the values in order to find max
max_h = max(h_values) # get M as the maximum value of h
M_func <- function(x) max_h # define M constant function

figure2 <- ggplot(data = data.frame(x=x), mapping = aes(x=x)) +
  stat_function(fun=h, mapping = aes(color="h(x)"), size=1.1) +
  stat_function(fun=M_func, mapping = aes(color="optimal M"), size=1.1) +
  scale_color_manual(name = "Function:",
                     values = c("darkblue", "indianred"), # Color specification
                     labels = c("h(x)", "optimal M")) +
  theme(legend.position = c(0.9, 0.7)) +
  labs(title = "Graph of h(x) and optimal M value") # title


figure2 # show figure


# show f(x) and M_optimal * g(x) in a single plot 
G <- function(x) max_h * g(x) # G = M*g(x)

figure3 <- ggplot(data = data.frame(x=x), mapping = aes(x=x)) +
  stat_function(fun=f, mapping = aes(color="f(x)"), size=1.1) +
  stat_function(fun=G, mapping = aes(color="M*g(x)"), size=1.1) +
  scale_color_manual(name = "Function:",
                     values = c("darkblue", "indianred"), # Color specification
                     labels = c("f(x)", "M*g(x)")) +
  theme(legend.position = c(0.9, 0.7)) +
  labs(title = "Graph of f(x) and M*g(x) with optimal M") # title


figure3 # show figure


# --- PRINTING RESULTS ---
# Results of simulation
deviation = round(100*abs(round(100*num_accepted_samples/num_total_samples,2)-
                            round(100*1/M, 2))/round(100*1/M, 2), 2)
accpt_prob = round(100*num_accepted_samples/num_total_samples,2)

# print(paste("SIMULATION RESULTS"))
# print(paste("Total number of samples generated:", num_total_samples))
# print(paste("Total number of samples accepted:", num_accepted_samples))
# print(paste("Simulation acceptance probability:", accpt_prob, "%"))
# print(paste("Theoretical acceptance probability:", round(100*1/M, 2), "%"))
# 
# print(paste("Deviation of simulation acceptance probability from theoretical:",
#             deviation, "%"))
# print(paste("Estimation for mean from accepted samples:", mean(accepted_samples)))
# print(paste("Estimation for standard deviation from accepted samples:", 
#             sd(accepted_samples)))
#   


cat(sprintf("SIMULATION RESULTS
Total number of samples generated: %d
Total number of samples accepted: %d
Simulation acceptance probability (SAP): %.2f %%
Theoretical acceptance probability (TAP): %.2f %%
Deviation of SAP from TAP: %.2f %%
Estimation for mean from accepted sample points: %.4f
Estimation for standard deviation from accepted sample points: %.4f",
num_total_samples,
num_accepted_samples,
accpt_prob,
round(100*1/M, 2),
deviation,
mean(accepted_samples),
sd(accepted_samples)))
