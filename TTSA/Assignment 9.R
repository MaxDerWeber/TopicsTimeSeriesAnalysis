## Assignment 9


# Question 15 ----

b <- 0.25
epsilon1 <- rnorm(1001)       # standard normal with n=1001 (to have 1000 periods with epsilon_t-1)
epsilon2 <- rt(1001, df = 5)  # t-distribution with 5 deg. of freedom, n = 1001
epsilon3 <- runif(1001)       # standard uniform distribution, n = 1001


x1 <- epsilon1[2:1001] + b * epsilon1[1:1000]
x2 <- epsilon2[2:1001] + b * epsilon2[1:1000]
x3 <- epsilon3[2:1001] + b * epsilon3[1:1000]


#am Anfang der loop x aus d1,x2,x3 auswählen für die functions

p_gram <- function(b,x,j){
  lambda = ((2*pi*j)/T)
  1/(2*pi*length(x)) *  ( (sum(x[i,t] * cos(lambda*t)))^2 + (sum(x[i,t] * sin(lambda*t)))^2 )
}

f_x = function(b,j) {
  lambda <- ((2*pi*j)/T)
  (1+b^2+2*b*cos(lambda))*(sigma/(2*pi))
  }


M <- ((length(x1-1)/2))
Whittle_estim <- function(b,x){
  for (j in 1:M){
    result <- result + log(f_x(b,j))-(p_gram(x,j)/f_x(b,j))
  }
}



output = matrix(nrow = 1000, ncol = 3, data = NA) #Contains estimates for b
output2 = matrix(nrow = 1000, ncol = 3, data = NA)#Contains variances 

x <- x1
for(s in 1:1000){
  #simulate data
  epsilon1 <- rnorm(1001)       # standard normal with n=1001 (to have 1000 periods with epsilon_t-1)
  epsilon2 <- rt(1001, df = 5)  # t-distribution with 5 deg. of freedom, n = 1001
  epsilon3 <- runif(1001)       # standard uniform distribution, n = 1001
  
  x1 <- epsilon1[2:1001] + b * epsilon1[1:1000]
  x2 <- epsilon2[2:1001] + b * epsilon2[1:1000]
  x3 <- epsilon3[2:1001] + b * epsilon3[1:1000]
  
  # estimation
  s <- 1
  sigma <-  var(x1)
  output[s,1] <- optimize(interval = c(0,1),maximum = TRUE, Whittle_estim(x = x1))$maximum #bug
  output[s,2] <- optimize(interval = c(0,1),maximum = TRUE, Whittle_estim(x = x2))$maximum #bug
  output[s,3] <- optimize(interval = c(0,1),maximum = TRUE, Whittle_estim(x = x3))$maximum #bug

  output2[s,1] <- var(x1[s])
  output2[s,1] <- var(x2[s])
  output2[s,1] <- var(x3[s])
  
  
  # Means and Biases
  
  # MSE
  
  # MLE
}







# Question 16 ----


m <- 2000
t <- 1000
epsilon <- rnorm(m+t)

Y <- matrix(data = NA, nrow = 5, ncol = 1000)
y_raw <- matrix(data = NA, nrow = 5, ncol = 3000)
corr <- matrix(data = NA, nrow = 5, ncol = 26)
D <- c(-0.45, -0.25, 0.0001, 0.25, 0.45)
psi <- matrix(data = NA, 5, 160)

for (i in 1:5){
#  y_raw <- fracdiff.sim(m+t, d = D[i], innov = epsilon )$series 

  # the gamma function can't handly very small numbers thats why I set the max j to 160 (formula from slides 157/159)
  d <- D[i]
  
  for (j in 1:160){
    psi[i,j] = gamma(j+d)/(gamma(j+1)*gamma(d)) # rows = d-values, columns = j_1:j_160
  }

  for (t2 in (1:3000)){
    y_raw[i,t2] <- sum(psi[i,1:160]*epsilon[(t2+160):(t2+1)])
  }

  Y[i,] <- y_raw[i,(m+1):(m+t)]
   
  acf_output <- acf(Y[i,],lag.max = 25, type = 'correlation',plot = F, na.action = na.pass)
  corr[i,] <- acf_output$acf
}

## plots

par(mfrow = c(3,2))
plot(Y[1,], type = 'l', col = 1, main = 'd = -0.45', ylim = c(-3,3))
plot(Y[2,], type = 'l', col = 2, main = 'd = -0.25', ylim = c(-3,3))
plot(Y[3,], type = 'l', col = 3, main = 'd = 0.0001', ylim = c(-3,3))
plot(Y[4,], type = 'l', col = 4, main = 'd = 0.25', ylim = c(-3,3))
plot(Y[5,], type = 'l', col = 8, main = 'd = 0.45', ylim = c(-3,3))
plot(1,1,legend("bottomright", c('d = -0.45', 'd = -0.25', 'd = 0', 'd = 0.25', 'd = 0.45'), col = c(1,2,3,4,8), lwd = c(2,2,2,2,2)))

# We can see that the persistance increases with d, as there are more persistant trends recognizable.
# Especially the trend around t=600 in this sample.


par(mfrow = c(3,2))
plot(corr[1,], type = 'h', col = 1, lwd = 3, main = 'd = -0.45')
abline(h = 0)
plot(corr[2,], type = 'h', col = 2, lwd = 3, main = 'd = -0.25')
abline(h = 0)
plot(corr[3,], type = 'h', col = 3, lwd = 3, main = 'd = 0')
abline(h = 0)
plot(corr[4,], type = 'h', col = 4, lwd = 3, main = 'd = 0.25')
abline(h = 0)
plot(corr[5,], type = 'h', col = 8, lwd = 3, main = 'd = 0.45', ylim = c(-0.1,1))
abline(h = 0)
plot(1,1,legend("bottomright", c('d = -0.45', 'd = -0.25', 'd = 0', 'd = 0.25', 'd = 0.45'), col = c(1,2,3,4,8), lwd = c(2,2,2,2,2)))

# It seems as if a smaller negative d-value produces more countercyclical autocorrelations and lower persistance.
# On the other hand, the larger the positive d-value,the higher the autocorrelations and therefore the persistance of the WN. 
# All autocorrelations go off with increasing h.
# This finding reflects our theoretical insight from slide 161 that autocorrelation grows with phi when d > 0.






