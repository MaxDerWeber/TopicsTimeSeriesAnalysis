## Assignment 6 - 

## Load packages ----
library(urca)
install.packages("sandwich")
library(sandwich)
library(tseries)


## Question 11 ----

# define counters
i1 <- 0
i2 <- 0
i3 <- 0

for (s in 1:1000){
  # Data Generation
  
  a <- 1
  beta_0 <- 1 
  
  epsilon <- rnorm(1000,0,1)
  e <- epsilon[3:1000] - 5/6*epsilon[2:999] + 1/6*epsilon[1:998]
  u <- cumsum(a*e) # works only for a = 1
  x <- beta_0 + u
  #plot(1:998,x, type = 'l') 
  
  ## Testing ----
  
  # Long-Run Variances with Andrew's
  
  omega2 <- lrvar(x, type = "Andrews", prewhite = F, adjust = F, kernel = "Quadratic Spectral", bw = length(x)^(4/5))
  omega3 <- lrvar(x, type = "Andrews", prewhite = F, adjust = F, kernel = "Quadratic Spectral", bw = length(x)^(1/3)) 
  # macht es Sinn, dass die LongRun Variance kleiner ist bei kleinerer BW? normalerweise steigt doch die Varianz mit sinkender BW?
  
  
  #  testing strategy 1): Dickey-Fuller without corrections
  #  t1 <- as.numeric(adf.test(x,alternative = "stationary", k = 0)$statistic)  # k defines lag order
  #  ?adf.test -> tseries package
  #  beta0 ist nicht =0 !! schauen ob das in einer der functions spezifiziert werden muss
  
  # ur.df urca-package
  #  t1new <- ur.df(x, type = "none")  #type= "none" means no drift; hier kein lag_order nur lags ?
  #  summary(t1new) # mit @ auf die einzelnen slots zugreifen (statt mit [] oder $ wie sonst üblich)
  #  t1new@teststat
  #  t1new@testreg$call #komische Regression
  
  
  
  x_lag1 <- x[1:length(x)-1]
  x_t <- x[2:length(x)]            
  ols <- lm(formula = x_t ~ x_lag1 )
  
  summary(ols)
  res_output <- ols$coefficients
  a_hat <- res_output[2]
  beta0_hat <- res_output[1]
  sigma_a1_hat <- (summary(ols)$coefficients[2,2])  #?? stimmt das ?? here Standard Error for the formular
  #  gamma0 <- acf(e, type = "covariance")$acf[1] # autocovar of e should be the variance of the errors(=residuals) in the ols?
  gamma0 <- (summary(ols)$sigma) ^ 2
  
  t1 <- as.numeric((a_hat -1) / sigma_a1_hat)
  
  # testing strategy 2) & 3): Phillips-Perron correction (1)
  t2 <- as.numeric((sqrt(gamma0 / (omega2^2)) * t1) - (((omega2^2 - gamma0)/(2*omega2)) * ((length(x) * sigma_a1_hat) / gamma0))) # length(x) equals T
  t3 <- as.numeric((sqrt(gamma0 / (omega3^2)) * t1) - (((omega3^2 - gamma0)/(2*omega3)) * ((length(x) * sigma_a1_hat) / gamma0))) # length(x) equals T
  
  # Fill  t1, t2, t3
  
  
  if(t1 < -2.68){i1 <- i1+1}    
  if(t2 < -2.68){i2 <- i2+1}    
  if(t3 < -2.68){i3 <- i3+1}    
  
}

P1 <- 1/1000 * i1   # P1 = 0.970
P2 <- 1/1000 * i2   # P2 = 1
P3 <- 1/1000 * i3   # P3 = 0.925


# All Test don't seem to perform very well, as they reject H0 most of the time, telling us, 
# that a1 is lower than 1 and x is stationary, which it isn't.
# The best performing test is strategy 3. With my little knowledge of non-parametrics,
# I suppose it is because of the smaller bandwidth, which results in more accurate omega_hat.
# q = T 4/5 (here 251.1886) seems unresonably large bandwidth to me.
# I want to add, that I excpected different results for P3 (something closer 0),
# as a correct PP correction should take care of the omitted short term variance in the ols estimation
# which is not corrected in the standard Dickey-Fuller test (P1).

