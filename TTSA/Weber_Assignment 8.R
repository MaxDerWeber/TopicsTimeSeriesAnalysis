## Assignment 8 Question 14

# plotting spectral density ----

sigma_squared <- 1 #from what is given
counter <- 1
L <- seq(from = 0,to = pi, by = 0.01)


spectral_density <- function(lambda){
  sigma_squared / (2*((2*q+1)^2)*pi) * (2*q +1 + (2 * sum((2*q+1-h)*cos(lambda*h))))
}


f_x_lambda <- matrix(data = NA, nrow = 3, ncol = 315) # empty matrix for function points


for (q in c(1,3,10)){
  h <- 1:(2*q)
  L <- seq(from = 0,to = pi, by = 0.01)
  
  f_x_lambda[counter,]  <- unlist(lapply(L,spectral_density))
  
  counter <- counter + 1
}

# plot the functions
plot(L, f_x_lambda[1,], type = 'l', lwd = 3, xlab = 'lambda',ylab = 'f (lambda)', main = 'spectral density')
points(L, f_x_lambda[2,], type = 'l', col = 'blue', lwd = 3)
points(L, f_x_lambda[3,], type = 'l', col = 'red', lwd = 3)
legend('topright', legend = c('q = 1', 'q = 3', 'q = 10'), col = c('black', 'blue', 'red'), lwd = c(3,3,3))


# calculations ----
integrant <- function(l){
  f <- sigma_squared / (2*((2*q+1)^2)*pi) * (2*q +1 + (2 * sum((2*q+1-h)*cos(l*h))))
  return(f)
}


s_lambda0 <- matrix(data = 0, nrow = 6, ncol = 315)
counter <- 1
counter_l <- 1

for (q in c(1,3,10)){
  for(l in L){
    h <- 1:(2*q)
    s_lambda0[counter,(l+0.01)*100] <- 2 * integrate(Vectorize(integrant), lower = 0, upper = l)$value
    counter_l <- counter_l + 1
  }
  counter <- counter + 1
}

s_lambda0[4,] <- s_lambda0[1,1:315]/s_lambda0[1,315]
s_lambda0[5,] <- s_lambda0[2,1:315]/s_lambda0[2,315]
s_lambda0[6,] <- s_lambda0[3,1:315]/s_lambda0[3,315]

plot(L, s_lambda0[4,], type = 'p', lwd = 2, xlab = 'lambda',ylab = 'spectral density', main = 'variance ratio')
points(L, s_lambda0[5,], type = 'p', col = 'blue', lwd = 2)
points(L, s_lambda0[6,], type = 'p', col = 'red', lwd = 2)
legend('bottomright', legend = c('q = 1', 'q = 3', 'q = 10'), col = c('black', 'blue', 'red'), lwd = c(3,3,3))

