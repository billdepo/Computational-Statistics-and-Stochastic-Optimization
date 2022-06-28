# 1.b - Rao-Blackwellization
library(ggplot2)

set.seed(42) # generate reproducible results across runs


log_seq <- function(start_point, end_point, num_points, log_base = 10) {
  sequence <- seq(log(start_point, log_base), 
                  log(end_point, log_base), length.out=num_points)
  return(log_base ** sequence)
}

n = 5000 # number of sample points for sample X = (X1, X2, ..., Xn)
r = 1 # number of successes to observe occurring
p = 0.5 # probability of success in a trial
std_est_sample = c()
rao_est_sample = c()

std_mean = c()
rao_mean = c()

std_var = c()
rao_var = c()

std_var_theor = c()
rao_var_theor = c()


#m = seq(2, 3000, 40)  # number of experiments (number of samples taken)
m = unique(floor(log_seq(2, 5000, num=40))) # selected a logarithmic range 

for (i in m){
  
  # take a sample point from NB(r,p) and add it to the standard estim. sample
  x_nb = replicate(i, rnbinom(n=n, size=r, prob=p)) # m samples of size n each

  # take a sample point lambda_i from Gamma(r,(1-p)/p) first and for that
  # lambda_i, get a sample point from Poisson(lambda_i). This is X_i
  lambda_i = replicate(i, rgamma(n=n, shape=r, scale = (1-p)/p))

  # calculate means
  std_est = colMeans(x_nb) # mean of columns 1 x m, this is the std estimator
  std_cur_mean = mean(std_est) # mean of all the estimator's values
  std_mean = append(std_mean, std_cur_mean)
  
  rao_est = colMeans(lambda_i) # mean of columns 1 x m, this is rao estimator
  rao_cur_mean = mean(rao_est)
  rao_mean = append(rao_mean, rao_cur_mean)
  
  # calculate variances
  std_cur_var = var(std_est)
  std_var = append(std_var, std_cur_var)
  
  rao_cur_var = var(rao_est)
  rao_var = append(rao_var, rao_cur_var) 
  
}


# create a data frame with all the values
df = data.frame(sims=m, # num of simulations ie sample size
                std_means=std_mean,
                rao_means=rao_mean,
                std_vars=std_var,
                rao_vars=rao_var,
                real_mean = 1,
                std_vars_theor = 2/n,
                rao_vars_theor = 1/n)


colors <- c("Rao estimator mean" = "darkblue", 
            "Standard estimator mean" = "indianred", 
            "Real mean value" = "royalblue1")

figure1 <- ggplot(data = df, aes(x=sims)) +
  geom_line(aes(y=std_means, color='Standard estimator mean'), size=0.9) +
  geom_line(aes(y=rao_means, color='Rao estimator mean'), size=0.9) +
  geom_line(aes(y=real_mean, color='Real mean value'), linetype="dashed", size=0.9) +
  labs(
    title = "Mean of Standard and Rao-Blackwellized estimators",
    subtitle = "Simulation for variable number of samples of size n=5000 each",
    x = 'Number of samples',
    y = 'Mean',
    color = 'Legend') +
  scale_color_manual(values = colors) +
  theme(legend.position = c(0.79, 0.20),
        legend.key.size = unit(0.5, 'cm'))

figure1 # show figure

colors <- c("Rao estimator variance" = "darkblue", 
            "Standard estimator variance" = "indianred",
            "Theoretical rao estimator variance" = "deepskyblue4",
            "Theoretical standard estimator variance" = "magenta")

figure2 <- ggplot(data = df, aes(x=sims)) +
  geom_line(aes(y=std_vars, color='Standard estimator variance'), size=0.9) +
  geom_line(aes(y=rao_vars, color='Rao estimator variance'), size=0.9) +
  geom_line(aes(y=rao_vars_theor, color='Theoretical rao estimator variance'), linetype="dashed", size=0.9) +
  geom_line(aes(y=std_vars_theor, color='Theoretical standard estimator variance'), linetype="dashed", size=0.9) +
  labs(
    title = "Variance of Standard and Rao-Blackwellized estimators",
    subtitle = "Simulation for variable number of samples of size n=5000 each",
    x = 'Number of samples',
    y = 'Variance',
    color = 'Legend') +
  #scale_y_continuous(trans='log10') +
  scale_color_manual(values = colors) +
  theme(legend.position = c(0.73, 0.81),
        legend.key.size = unit(0.35, 'cm'))

figure2 # show figure

#------------------------------------------
# NOW SIMULATE FOR VARIABLE SAMPLE SIZE n AND FIXED NUMBER OF SIMULATIONS 


k = 300 # number of simulated samples
std_est_sample2 = c()
rao_est_sample2 = c()

std_mean2 = c()
rao_mean2 = c()

std_var2 = c()
rao_var2 = c()

std_var_theor2 = c()
rao_var_theor2 = c()
sample_size = c()


sequence =  unique(floor(log_seq(2, 10000, num=150))) # selected a logarithmic range 

for (i in sequence){
  # take k samples of size i from NB(r,p)
  x_nb2 = replicate(k, rnbinom(n=i, size=r, prob=p)) #sample of size i
  
  # take k samples of size i of L_i points from Gamma(r,(1-p)/p) 
  lambda_i2 = replicate(k, rgamma(n=i, shape=r, scale = (1-p)/p))
  
  # calculate means
  std_est2 = colMeans(x_nb2) # mean of columns 1 x m, this is the std estimator
  std_cur_mean2 = mean(std_est2) # mean of all the estimator's values
  std_mean2 = append(std_mean2, std_cur_mean2)
  
  rao_est2 = colMeans(lambda_i2) # mean of columns 1 x m, this is rao estimator
  rao_cur_mean2 = mean(rao_est2)
  rao_mean2 = append(rao_mean2, rao_cur_mean2)
  
  # calculate variances
  std_cur_var2 = var(std_est2)
  std_var2 = append(std_var2, std_cur_var2)
  
  rao_cur_var2 = var(rao_est2)
  rao_var2 = append(rao_var2, rao_cur_var2) 
  
  print(paste('Done for sample size i=', i))
}

for (i in sequence){
  std_var_theor2 = append(std_var_theor2, 2/i)
  rao_var_theor2 = append(rao_var_theor2, 1/i)
  sample_size = append(sample_size, i)
}

# create a data frame with all the values
df2 = data.frame(sims=sample_size, # num of simulations ie sample size
                std_means=std_mean2,
                rao_means=rao_mean2,
                std_vars=std_var2,
                rao_vars=rao_var2,
                real_mean = 1,
                std_vars_theor = std_var_theor2,
                rao_vars_theor = rao_var_theor2)


colors <- c("Rao estimator mean" = "darkblue", 
            "Standard estimator mean" = "indianred", 
            "Real mean value" = "royalblue1")

figure3 <- ggplot(data = df2, aes(x=sims)) +
  geom_line(aes(y=std_means, color='Standard estimator mean'), size=1.1) +
  geom_line(aes(y=rao_means, color='Rao estimator mean'), size=1.1) +
  geom_line(aes(y=real_mean, color='Real mean value'), linetype="dashed", size=0.9) +
  labs(
    title = "Mean of Standard and Rao-Blackwellized estimators",
    subtitle = "Simulation for k=300 samples of variable size",
    x = 'Sample size n',
    y = 'Mean value',
    color = 'Legend') +
  scale_color_manual(values = colors,
                     guide=guide_legend(override.aes=list(linetype=c(1,1,2), lwd=c(1,1,0.5)))) +
  theme(legend.position = c(0.9, 0.85),
        legend.key.size = unit(0.5, 'cm'))

figure3 # show figure


colors <- c("Rao estimator variance" = "darkblue", 
            "Standard estimator variance" = "indianred",
            "Theoretical rao estimator variance" = "deepskyblue4",
            "Theoretical standard estimator variance" = "magenta")

figure4 <- ggplot(data = df2, aes(x=sims)) +
  geom_line(aes(y=std_vars, color='Standard estimator variance'), size=1.1) +
  geom_line(aes(y=rao_vars, color='Rao estimator variance'), size=1.1) +
  geom_line(aes(y=rao_vars_theor, color='Theoretical rao estimator variance'), linetype="dashed", size=0.9) +
  geom_line(aes(y=std_vars_theor, color='Theoretical standard estimator variance'), linetype="dashed", size=0.9) +
  labs(
    title = "Variance of Standard and Rao-Blackwellized estimators ",
    subtitle = "Simulation for k=300 samples of variable size",
    x = 'Sample size n',
    y = 'Variance value',
    color = 'Legend') +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = colors,
                     guide=guide_legend(override.aes=list(linetype=c(1,1,2,2), lwd=c(1,1,0.5,0.5))))+
  theme(legend.position = c(0.85, 0.85),
        legend.key.size = unit(0.5, 'cm'))

figure4 # show figure
