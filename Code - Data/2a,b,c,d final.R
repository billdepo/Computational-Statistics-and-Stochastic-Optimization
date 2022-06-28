# 2 - pdf estimation with Epanechnikov kernel

set.seed(42) # generate reproducible results across runs

# --- IMPORTS ---

library(ggplot2)

# --- FUNCTIONS ---

epanechnikov_kernel <- function(u){
  if ((u <= 1) && (u >= -1)){ #if |u|<=1
    kernel = (3/4) * (1-u^2) # Epanechnikov K(u) calculation
  } else {
    kernel = 0 # if |u|>1
  }
  return(kernel)
}

estimated_pdf <- function(x, sample, h){
  n = length(sample)
  summ = 0
  for (i in 1:n){ # calculate Σ[Κ((x-x_i)/h)]
    summ = summ + epanechnikov_kernel((x-sample[i])/h)
  }
  pdf = (1/(n * h)) * summ
  return(pdf)
}


# below function is used in 2b and 2c so that a single pdf function 
#includes the kernel implementation inside it too
full_pdf <- function(x, sample=faithful$eruptions, h=0.208){
  n = length(sample)
  summ = 0
  for (i in 1:n){ # calculate Σ[Κ((x-x_i)/h)]
    u = (x-sample[i])/h
    if ((u <= 1) & (u >= -1)){ #if |u|<=1
      kernel = (3/4) * (1-u^2) # Epanechnikov K(u) calculation
    } else {
      kernel = 0 # if |u|>1
    }
    summ = summ + kernel
  }
  pdf = (1/(n * h)) * summ
  return(pdf)
}


# --- MAIN ---
#2a
data = faithful$eruptions # load Old Faithful geyser data
n = length(data) # sample size

x_seq = seq(-2,2,0.1)
ep_kernel = c()
for (i in x_seq){
  ep_kernel =  append(ep_kernel, epanechnikov_kernel(i))
}

#plot the epanechnikov kernel
figure1 <- ggplot(data = data.frame(x=x_seq), aes(x=x_seq, y=ep_kernel)) +
  geom_line(color='indianred', aes( y=ep_kernel), size=1.1) +
  labs(title = "Graph of the epanechnikov kernel function", x="u", y = "K(u)") # title

figure1 # show figure

cat(sprintf("Max value in data: %.2f
Min value in data: %.2f
Max-min range in data: %.2f
Mean value of data: %.3f
Standard deviation of data: %.3f",
max(data),
min(data),
abs(max(data)-min(data)),
mean(data),
sd(data)))


# print(paste("Max value in data:", max(data)))
# print(paste("Min value in data:", min(data)))
# print(paste("Mean value in data:", mean(data)))
# print(paste("Standard deviation in data", sd(data)))
# print(paste("Max-Min range in data:", abs(max(data)-min(data))))

# we use the above printed results to select an appropriate range for h 

#h_seq = seq(0.05, sd(data), by=0.01) 

log_seq <- function(start_point, end_point, num_points, log_base = 10) {
  sequence <- seq(log(start_point, log_base), 
                  log(end_point, log_base), length.out=num_points)
  return(log_base ** sequence)
}

h_seq = log_seq(0.01, sd(data), num=30) # selected a logarithmic range for h

cv_log_likelihood_seq = c()

for (h_index in 1:length(h_seq)){ # iterate over each h value in the range
  h = h_seq[h_index] # get the selected h value
  # cross-validation method, remove 1 observation each time
  cv_log_likelihood = 0 # take cross validated log likelihood
  for (i in 1:n){ # iterate over all observations
    x_i = data[i] # take the i observation
    cv_data = data[-i] # exclude i observation from sample
    cv_pdf_on_i = estimated_pdf(x_i, cv_data, h)
    cv_log_likelihood = cv_log_likelihood + log(cv_pdf_on_i) # log() = ln()
  }
  # print(paste("h=", h, " Log likelihood=", cv_log_likelihood))
  cv_log_likelihood_seq = append(cv_log_likelihood_seq, cv_log_likelihood)
}

# plot(x=h_seq[100:300], y=cv_log_likelihood_seq[100:300], type="l") # plot log likelihood function


#plot the log likelihood vs h
df_ll = data.frame(x=h_seq, y=cv_log_likelihood_seq)
figure2 <- ggplot(data = df_ll, aes(x=x,y=y)) +
  geom_line(color='indianred', size=1.1) +
  labs(title = "Cross-validated log-likelihood vs h",
       x="h", y = "LL(h)") # title

figure2 # show figure


# repeat the process for the narrowed h range
# to narrow it down we find the max and second max value
ll_max = max(cv_log_likelihood_seq) # max LL(h)
idx_max = which(cv_log_likelihood_seq == ll_max) # index of max
ll_2max = max(cv_log_likelihood_seq[cv_log_likelihood_seq != ll_max]) # 2nd max
idx_2max = which(cv_log_likelihood_seq == ll_2max) # index of max
# define new h range to perform a seach for the optimal h
h_start = min(h_seq[idx_max], h_seq[idx_2max])
h_end = max(h_seq[idx_max], h_seq[idx_2max])

#h_seq = seq(h_start, h_end, by=0.001) 
h_seq = log_seq(h_start, h_end, num=30) # selected a logarithmic range for h

cv_log_likelihood_seq = c()

for (h_index in 1:length(h_seq)){
  h = h_seq[h_index]
  cv_log_likelihood = 0 
  for (i in 1:n){ 
    x_i = data[i] 
    cv_data = data[-i] 
    cv_pdf_on_i = estimated_pdf(x_i, cv_data, h)
    cv_log_likelihood = cv_log_likelihood + log(cv_pdf_on_i)
  }
  cv_log_likelihood_seq = append(cv_log_likelihood_seq, cv_log_likelihood)
}

# plot again to get a more zoomed and accurate picture around the max
df_ll = data.frame(x=h_seq, y=cv_log_likelihood_seq)
figure3 <- ggplot(data = df_ll, aes(x=x,y=y)) +
  geom_line(color='indianred', size=1.1) +
  xlim(h_start,h_end) +
  ylim(min(cv_log_likelihood_seq), max(cv_log_likelihood_seq)) +
  labs(title = "Cross-validated log-likelihood vs h",
       x="h", y = "LL(h)") # title

figure3 # show figure


max_log_likelihood = max(cv_log_likelihood_seq)
index_max_log_likelihood = which.max(cv_log_likelihood_seq) # index of max
h_optimal = h_seq[index_max_log_likelihood] # corresponding h for max log likelihood

print(paste('The optimal h =', h_optimal))

#plot(density(data, kernel = "epanechnikov")) # plot pdf using h that R calculates itself


# background color only for the plot region
plot.new() 
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "#ebebeb") # background color of ggplot


# add a new plot
par(new = TRUE)


plot(density(data, bw=h_optimal, kernel = "epanechnikov"),
     col='indianred',
     lwd=3,
     main="Estimated f(x) pdf using density function",
     sub='h=optimal=0.208, kernel=Epanechnikov',
     xlab='x (mins)') # plot pdf using optimal calculated h
grid(col = "white")
polygon(density(data, bw=h_optimal, kernel = "epanechnikov"),
        col = rgb(205/256, 92/256, 92/256, alpha = 0.6)) # color under curve


#2b
# get density function's x points (default = 512 points)
dens = density(data, bw=h_optimal, kernel = "epanechnikov")

x_seq = sort(dens$x)
pdf_values = c()
for (i in 1:length(x_seq)){
  x_i = x_seq[i]
  pdf = full_pdf(x_i) # estimated pdf value for each of the x points
  pdf_values = append(pdf_values, pdf) # add to a list
}


# background color only for the plot region
plot.new() 
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "#ebebeb") # background color of ggplot


# add a new plot
par(new = TRUE, cex=0.8)

plot(density(data, bw=h_optimal, kernel = "epanechnikov"),
     col='indianred',
     ylim=c(0,max(pdf_values)), 
     lwd=3,
     main='Estimated pdf for optimal h \n 
          Custom vs R density function implementation',
     sub='h=optimal=0.208, kernel=Epanechnikov',
     xlab='x (mins)') # plot pdf using optimal calculated h
grid(col="white")
polygon(density(data, bw=h_optimal, kernel = "epanechnikov"),
        col = rgb(205/256, 92/256, 92/256, alpha = 0.6)) # color under curve
par(new=TRUE)
lines(x_seq, pdf_values, type="l", col='darkblue', lwd=3)
polygon(c(min(x_seq), x_seq, min(x_seq)), c(0, pdf_values, 0),
        col = rgb(55/256, 198/256, 255/256, alpha = 0.2)) # color under curve
legend(
  'topright',
  title = "Implementation",
  col = c('indianred', 'darkblue'),
  lty=1:1,
  legend = c('density \n function', 'custom')
)



# plot(density(data, bw=h_optimal, kernel = "epanechnikov"), 
#      col='indianred',
#      ylim=c(0,max(pdf_values)), 
#      xlab='t (mins)',
#      main='Estimated pdf f for optimal h (0.21) - custom vs R density() implementation')
# polygon(density(data, bw=h_optimal, kernel = "epanechnikov"),
#         col = rgb(205/256, 92/256, 92/256, alpha = 0.6)) # color under curve
# lines(x_seq, pdf_values, type="l", col='blue')
# legend(
#   'topright',
#   title = "Implementation",
#   col = c('indianred', 'black'),
#   lty=1:1,
#   legend = c('R density()', 'custom')
# )

#2c
probability = integrate(Vectorize(full_pdf), lower =3.5, upper = Inf) 
# print(paste("Probability that the time interval between successive volcano \n
# explosions is more than 3.5 minutes:", probability$value, "with ", 
# probability$abs.error, "absolute error"))


cat(sprintf("Probability that the time interval between successive volcano
explosions is more than 3.5 minutes: %.5f with %.5f absolute error",
probability$value,
probability$abs.error))

#2d
set.seed(42)
sample_size = 250 # number of sample points to simulate from estimated f
f_simul_values = c() # list of simulated f values


# for (i in 1:sample_size){
#   # Step 1: randomly select an observation from the initial sample with 1/n probability (n=size of sample=272)
#   random_obs = sample(data, 1) # we call this a random variable X
#   
#   # Step 2: simulate a random variable Y from the kernel's distribution
#   # here we need to generate U1,U2,U3 ~ Uniform(-1,1)
#   # and then if |U3|>=|U2| & |U3|>=|U1|: keep U2, else U3
#   u_values = runif(3, min=-1, max=1) # generate three uniform values Uniform(-1,1)
#   u1 = u_values[1]
#   u2 = u_values[2]
#   u3 = u_values[3]
#   
#   if ( (abs(u3) >= abs(u2)) && (abs(u3) >= abs(u1)) ){
#     u = u2
#   } else{
#     u = u3
#   }
#   # now u is the value of the random variable Y
#   
#   
#   # Step 3: the random sample point from estimated f is Z = X + h * Y
#   f_sample = random_obs + h * u
#   f_simul_values = append(f_simul_values, f_sample)
# }

hY <- apply(matrix(runif(3*sample_size, -h_optimal, h_optimal), 3), 2, median)
x = sample(data, sample_size, replace = TRUE)
f_simul_values = x + hY # Z = X + h*Y sampling


simul_probability = length(which(f_simul_values > 3.5)) / sample_size
# print(paste("Simulated probability that the time interval between successive volcano explosions is more than 3.5 minutes:", simul_probability))

cat(sprintf("Simulated probability that the time interval between successive
volcano explosions is more than 3.5 minutes: %.5f",
simul_probability))
