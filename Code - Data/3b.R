# 3b. EM for unknown parameter(s) estimation

set.seed(42) # generate reproducible results across runs

# --- IMPORTS ---

library(ggplot2)

# --- FUNCTIONS ---

EM_algorithm <- function(x, p, l, tolerance=1e-10){
  n = length(x) # length of sample / number of observations
  m = length(p) # number of populations/mixture distributions 
  W = matrix(,nrow=n, ncol=m) # create empty matrix for w_ij
  # save in a list the p1,p2,l1,l2 values we find over iterations
  p1_values = c(p[1])
  p2_values = c(p[2])
  l1_values = c(l[1])
  l2_values = c(l[2])
  diffs = c()
  
  repeat{ # repeat for as long as the condition concerning tolerance is satisfied
    
    # E-step - recalculate W matrix values  (wij) 
    for (i in 1:n){ # iterate over rows (observations)
      
      # calculate denominator
      denominator = 0 # this is the sum we see in the wij equation's denominator
      for (j in 1:m){ # iterate over columns (populations)
        denominator = denominator + p[j] * dpois(x[i], lambda=l[j]) 
      }
      
      # calculate wij
      for (j in 1:m){ # iterate over columns (populations)
        W[i,j] = p[j] * dpois(x[i], lambda=l[j]) / denominator # update matrix
      }
    }
    
    # M-step - update p and l values
    for (j in 1:m){
      
      p[j] = 0 # new estimations of p1, p2
      l[j] = 0 # new estimations of l1, l2
      for (i in 1:n){
        p[j] = p[j] + W[i, j]
        l[j] = l[j] + W[i, j] * x[i]
      }
      l[j] = l[j] / p[j]
      p[j] = p[j] / n
    }
    p1_values = append(p1_values, p[1])
    p2_values = append(p2_values, p[2])
    
    l1_values = append(l1_values, l[1])
    l2_values = append(l2_values, l[2])
    
    # repeat condition to exit
    l1_last_two_values = tail(l1_values, 2)
    l2_last_two_values = tail(l2_values, 2)
    difference = (l1_last_two_values[2] - l1_last_two_values[1])^2
    difference = difference + (l2_last_two_values[2] - l2_last_two_values[1])^2
    if (i==1){ # for the first iteration keep two diffs
      diffs = append(diffs, difference)
    }
    diffs = append(diffs, difference)
    if (difference <= tolerance){
      break # exit the repeat loop
    }
  } # end of repeat loop
  
  diffs = append(diffs, difference) # duplicate last diff
  return(list(p1=p1_values, p2=p2_values, 
              l1=l1_values, l2=l2_values, diffs=diffs))
}


pi_plots <- function(result){
  
  df = data.frame(p1=result$p1, p2=result$p2, 
                  l1=result$l1, l2=result$l2, diffs=result$diffs)
  
  # p1, p2 plots
  colors <- c("π1" = "darkblue", "π2" = "indianred")
  
  ggplot(data=df, aes(x=1:nrow(df))) +
    geom_line(aes(y=p1, color='π1'), size = 1) + 
    #geom_point(aes (y=p1, color='π1'), size = 3) +
    geom_line(aes(y=p2, color='π2'), size = 1) +
    #geom_point(aes (y=p2, color='π2'), size = 3) +
    labs(
      title = 'π1 and π2 values over EM iterations',
      x = 'iteration',
      y = 'p-value',
      color = 'Legend') +
    scale_color_manual(values = colors) +
    theme(legend.position = c(0.9, 0.85))
}

lambda_plots <- function(result){
  
  df = data.frame(p1=result$p1, p2=result$p2, 
                  l1=result$l1, l2=result$l2, diffs=result$diffs)
  
  # l1, l2 plots
  colors <- c("λ1" = "darkblue", "λ2" = "indianred")
  
  ggplot(data=df, aes(x=1:nrow(df))) +
    geom_line(aes(y=l1, color='λ1'), size = 1) + 
    #geom_point(aes(y=l1, color='λ1'), size = 3) + 
    geom_line(aes(y=l2, color='λ2'), size = 1) + 
    #geom_point(aes(y=l2, color='λ2'), size = 3) +  
    labs(
      title = 'λ1 and λ2 values over EM iterations',
      x = 'iteration',
      y = 'λ value',
      color = 'Legend') +
    scale_color_manual(values = colors) +
    theme(legend.position = c(0.9, 0.85))
}

convergence_plot <- function(result){
  df = data.frame(p1=result$p1, p2=result$p2, 
                  l1=result$l1, l2=result$l2, diffs=result$diffs)
  
  # difference (convergence criterion) plot
  conv_limit <- function(x) 1e-10 # define M constant function
  
  ggplot(data=df, aes(x=1:nrow(df))) +
    geom_line(aes(y=diffs, color='difference'), size = 1) + 
    scale_y_continuous(trans='log10') +
    stat_function(fun=conv_limit, mapping = aes(color="1e-10"), size=1.1) +
    scale_color_manual(name = "Function:",
                       values = c("darkblue", "indianred"), # Color specification
                       labels = c("difference", "1e-10")) +
    theme(legend.position = c(0.9, 0.7)) +
    labs(
      title = 'Convergence criterion value over EM iterations',
      x = 'iteration',
      y = 'criterion value')
  
}

print_results <- function(result){
  
  num_iterations = length(result$p1)
  
  p1=tail(result$p1, 1)
  p2=tail(result$p2, 1)
  l1=tail(result$l1, 1)
  l2=tail(result$l2, 1)
  
  # print(paste('Number of iterations before convergence:', num_iterations))
  # print(paste('λ1 =', l1, 'λ2 =', l2))
  # print(paste('π1 =', p1, 'π2 =', p2))
  
  cat(sprintf("Number of iterations before convergence: %d
  λ1: %.3f
  λ2: %.3f
  π1: %.3f
  π2: %.3f",
              num_iterations,
              l1,
              l2,
              p1,
              p2))
  
}

# --- MAIN ---
#3b

sample = c(2,7,3,9) # data observations

# initialize p and lambda values
p1 = runif(1, min=0, max=1) # random value for p1 in [0,1] and then p2=1-p1
p = c(p1, 1-p1) # p = [p1, p2]

# lambda = [lambda1, lambda2] random initialization with a value between [1,10]
l = runif(2, min=1, max=20)  

result = EM_algorithm(sample, p, l) # EM algorithm

# ---  PRINT RESULTS --- 

print_results(result)

# ---  PLOTS --- 

pi_plots(result) # display plot of converging values p1, p2
lambda_plots(result) # display plot of converging values l1, l2
convergence_plot(result) # convergence criterion values in each iteration


# run with the values that were found theoreticall with MLE in 2a
p = c(0.5, 0.5)
l = c(2.5, 8)

result = EM_algorithm(sample, p, l)

# ---  PRINT RESULTS --- 

print_results(result)
