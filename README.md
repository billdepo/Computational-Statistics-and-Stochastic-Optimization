#Computational-Statistics-and-Stochastic-Optimization
Semester project in the context of the <a href="https://dsml.ece.ntua.gr/">DSML NTUA</a> MSc's course <a href="https://dsml.ece.ntua.gr/studies/courses/ypologistike-statistike-kai-stochastike-beltistopoiese">"Computational Statistics and Stochastic Optimization"</a>.

# Introduction
The project employs a range of statistical and stochastic methods in order to solve a range of problems by yielding the power of the programming language R and whenever possible verify simulation results theoretically. We briefly describe below the four assignments this report deals with.

# Exercise 1
Here, we explore the rejection sampling method to sample from a known distribution. Next, we experiment with an expected value estimator and its Rao-Blackwellized version that is used instead of the original estimator to reduce the estimator's variance while maintining the property of being unbiased (same expected value as the parameter for estimation). Moreover, values from a statistical function of random variables following the uniform distribution are simulated and a histogram of the resulting distribution is generated. Finally, we calculate the function’s standard error using Bootstrap and Jackknife methods and generate a histogram to compare with the previous method.

# Exercise 2
Here we examine the ”Old Faithful geyser” dataset, which contains the times in minutes of consecutive explosions of a volcano in the USA. As a first step, we estimate the probability density function (pdf) of the data using an Epanechnikov kernel and determine the optimal kernel width by maximizing the cross-validated likelihood. The estimated pdf is plotted in two different ways and the results are compared. Lastly, the probability for two consecutive explosions to be more than 3.5 minutes apart is calculated by integrating an appropriate region under the estimated pdf curve as well as by simulating values from the estimated pdf and the two results are then compared.

# Exercise 3
This exercise demonstrates how Maximum Likelihood is used for parameter estima- tion and subsequently the Expectation Maximization algorithm is employed for the same task when some information is missing.

# Exercise 4
Finally, exercise four explores the full space of possible models in the feature/variable selection problem given a diabetes prediction dataset containing ten variables by utilizing the BIC criterion. In addition to that, the Lasso method is applied for the same task and cross-validation is used to calculate the lambda penalty parameter of the method. Subsequently, based on two different parameter values, two models are specified by the inclusion of appropriate predictors to the models, and finally we compare all three models on their predictive power using 5-fold cross-validation and each model’s Root Mean Square Error (RMSE) on the validation fold each time.