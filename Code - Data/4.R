#4 - diabetes prediction models

set.seed(4) # generate reproducible results across runs


# --- IMPORTS ---
library(ggplot2)
#install.packages("lars") # install lars package just once and then comment out
library(lars)
#install glmnet package, run the below once and then comment out
# install.packages("glmnet", repos = "https://cran.us.r-project.org")

library(glmnet)


# --- FUNCTIONS ---

# convert integer to binary list
int2bin <- function(integer, num_of_bits=10) {
  binary_vector = rev(as.numeric(intToBits(integer)))
  return(binary_vector[-(1:(length(binary_vector) - num_of_bits))])
}

# create string that says which variables are used
concat_string <- function(a_list) {
  full_string = ""
  for (i in 1:length(a_list)) {
    full_string = paste(full_string, a_list[i]) # add variable name
    if (i!=length(a_list)){
      full_string = paste(full_string, " + ") # add sign symbol between
    }
  }
  return(full_string)
}

# --- MAIN ---
# read data
data("diabetes")
attach(diabetes)

# create columns for all variables, first column=y=target variable
dataset = data.frame(cbind(y = diabetes$y, diabetes$x))

# display useful info for the dataset and variables
str(dataset) # variables and their types
summary(dataset) # min, mean, max, median etc. per variable
anyNA(dataset) # check for NA values -> 0 NA values


# 4a
dataset_columns = colnames(dataset) # 11 column names incl. target y
indep_vars = dataset_columns[-1] # keep 10 independent variables' names only
print(indep_vars)
num_models = 2^length(indep_vars) # number of models 2^10

list_of_ = c(1)
models_df = data.frame(
  id = c(0:1023)
)
models_df$num_of_vars = NA
models_df$vars_used = NA
models_df$BIC = NA

# add the null model
models_df$vars_used[1] = 1 
models_df$num_of_vars[1] = 0 
null_model = lm( y~1, data = dataset["y"]) # null model
models_df$BIC[1] = BIC(null_model) # its BIC


for (i in 2:num_models){ # 2^10 = 1,024 models in total, null already indexed
  variables_mask = int2bin(i-1) # binary mask for variables selection
  
  # select the appropriate variables each time
  variables_selected = indep_vars[which(variables_mask==1)] # select vars
  models_df$num_of_vars[i] = sum(variables_mask) # update df on number of vars
  models_df$vars_used[i]=concat_string(variables_selected) # update df on vars used
  variables_selected = append(variables_selected, "y") # add y to selected vars
  sub_dataset = dataset[variables_selected] # dataset with selected vars and y
  
  # for each set of variables selected, use them all to fit linear model
  # this is done through the sub_dataset usage (limited dataset each time)
  model = lm( y~., data = sub_dataset)
  models_df$BIC[i] = BIC(model) # update df with BIC value for each model
}

best_model_index = which.min(models_df$BIC)
best_model_variables = models_df[best_model_index, ]$vars_used
best_model_BIC = models_df[best_model_index, ]$BIC

print("Best model with full enumeration - BIC criterion ")
print(paste("Variables:", best_model_variables))
print(paste("BIC:", best_model_BIC))

# now that we have found the best model, build it and name it M1
M1 = lm(y~sex+bmi+map+hdl+ltg, data=dataset)
summary(M1) # print info on the fitted model


# 4b

fit <- glmnet(x=dataset[, c('age', 'sex', 'bmi', 'map', 'tc', 
                            'ldl', 'hdl', 'tch', 'ltg', 'glu')],
              y=dataset[, 'y'],
              alpha=1) # fit is an object of class glmnet
# plot(fit, label=TRUE) # visualize coefficients vs L1 norm
plot(fit, xvar='lambda', label=TRUE) # visualize coefficients vs log(lambda) values
title("Coefficients vs log(lambda)", line = 2.5, adj = 0)
# coef(fit, s = 20) # show coefficients of model with selected lambda

# cross-validation stochasticness -> setting seed for reproducible results
set.seed(4)
cvfit <- cv.glmnet(x, as.matrix(y), alpha=1, type.measure = "mse")
plot(cvfit, ylab='CV-MSE') # CV-MSE vs log(lambda)
title("CV-MSE vs log(lambda)", line = 2.5, adj = 0)


# lambda.min is the value of λ that gives minimum mean cross-validated error, 
# while lambda.1se is the value of λ that gives the most regularized model such 
# that the cross-validated error is within one standard error of the minimum

# M2
l_min = cvfit$lambda.min # lambda that gives minimum mean cv error
M2_coeffs = coef(cvfit, s = "lambda.min") # check the coefficients
print("Coefficients of fitted Lasso model:")
print(M2_coeffs)

M2 = lm(y~sex+bmi+map+tc+hdl+ltg+glu, data=dataset) # define model at mean cv error

# M3
l_1se = cvfit$lambda.1se
M3_coeffs = coef(cvfit, s = "lambda.1se") # check the coefficients
print("Coefficients of fitted Lasso model:")
print(M3_coeffs)
# we can see age, sex, tc, ldl, tch, glu coefficients reached 0
# and we will exclude them from the model

M3 = lm(y~bmi+map+hdl+ltg, data=dataset) # define model at 1 se with some vars removed
#predict(M3, newx = x[1:2,], type = "response") # make a pred


# 4c

RMSE_func <- function(y_preds, y_real){
  diff = y_preds - y_real
  diff_sq = diff**2
  RMSE = sqrt(sum(diff_sq)/length(diff_sq))
  return(RMSE)  
}

set.seed(4)

folds_num = 5
dataset_shuffled = dataset[sample(nrow(dataset)), ] # shuffle data randomly
folds = cut(seq(1, nrow(dataset_shuffled)), breaks=folds_num, labels=FALSE) # 5 folds creation

RMSE1_seq = c() # keep the RMSE value on the test folds for the 5 folds tested
RMSE2_seq = c()
RMSE3_seq = c()

# performing 5-fold cv
for(i in 1:folds_num){
  # segment your data by fold using the which() function 
  test_idx = which(folds==i, arr.ind=TRUE) # get indexes for test set
  test_set = dataset_shuffled[test_idx, ] # get the test fold for test set
  train_set = dataset_shuffled[-test_idx, ] # get the rest folds for train set
  
  M1_preds = predict(M1, test_set) # predictions of M1 on test set
  RMSE1 = RMSE_func(M1_preds, test_set[['y']])
  RMSE1_seq = append(RMSE1_seq, RMSE1)
  
  M2_preds = predict(M2, test_set) # predictions of M2 on test set
  RMSE2 = RMSE_func(M2_preds, test_set[['y']])
  RMSE2_seq = append(RMSE2_seq, RMSE2)
  
  M3_preds = predict(M3, test_set) # predictions of M3 on test set
  RMSE3 = RMSE_func(M3_preds, test_set[['y']])
  RMSE3_seq = append(RMSE3_seq, RMSE3)
  
}

RMSE1_mean = mean(RMSE1_seq)
RMSE2_mean = mean(RMSE2_seq)
RMSE3_mean = mean(RMSE3_seq)

cat