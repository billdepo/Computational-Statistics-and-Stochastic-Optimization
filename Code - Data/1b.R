# 1.b - Rao-Blackwellization

set.seed(42) # generate reproducible results across runs

# --- FUNCTIONS ---

rao_blackwellization_sample <- function(n, r, p) {
  #generate n random values that follow gamma distribution Gamma(r,(1-p)/p)
  gamma_sample = rgamma(n, shape=r, rate=(1-p)/p) # generate gamma sample
  conditional_sample = c() # save the final conditional samples
  for (i in 1:n){
    gamma_sample_point = gamma_sample[i]
    conditional_sample_point = rpois(1, lambda=gamma_sample_point) #generate one conditional sample point from poisson
    conditional_sample = append(conditional_sample, conditional_sample_point)
  }
  return(conditional_sample)
}

# --- MAIN ---
n = 5000 # number of sample points
r = 1 # number of successes to observe occurring  
p = 0.5 # probability of success in a trial
# then the random number of failures we have seen, X, will have the negative binomial distribution:
# X ~ NB(r,p)

# draw a sample from negative binomial 
neg_bionomial_sample = rnbinom(n, size = r, prob = p) # draw n nbinomially distributed values
mean(neg_bionomial_sample)
var(neg_bionomial_sample)

# draw a sample using the rao-blackwellization method
conditional_sample = rao_blackwellization_sample(n, r, p)
mean(conditional_sample)
var(conditional_sample)


# EDW PREPEI NA PARW VARIANCE TOY EKTIMITI??

# 
# 
# hist(y_rnbinom,                                          # Plot of randomly drawn nbinom density
#      breaks = 100,
#      main = "")
# 

