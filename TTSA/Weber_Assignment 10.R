## Assignment 10 Q17

t <- 1000
m <- 2000

# Local Whittle and log PR estimation for different m(alpha)
alpha <- seq(from = 0.2, to = 0.8, by = 0.1)
m_alpha <- length(y) ^ alpha 

# function to minimize LW estimator

LW_R_d <- function(d, y){
  M <- m_alpha[i] #loopen counter
  first_part <- matrix(data = NA, ncol = 1, nrow = M)
  last_part <- matrix(data = NA, ncol = 1, nrow = M)
  
  for (j in 1:M){
    lambda <- 2 * pi * j / length(y)
    part1 <- (1/(2*pi*length(y))) 
    part2 <- sum(y * cos(lambda*(1:length(y))))
    part3 <- sum(y * sin(lambda*(1:length(y))))
    denominator <- lambda^(-2*d)
    first_part_point <- part1 * (part2^2 + part3^2) / denominator
    
    first_part[j] <- first_part_point
    last_part[j] <- log(lambda)
  }
  R_LW <- log(1/M * sum(first_part)) - (2*d/M) * sum(last_part)
  return(R_LW)
}




### Monte Carlo
test_statLW <- matrix(NA, nrow = 1000, ncol = 7) #matrix containing the 7 test statistics for all simulations
test_statPR <- matrix(NA, nrow = 1000, ncol = 7)
d_LW_est <- matrix(NA, nrow = 1000, ncol = 7)
d_PR_est <- matrix(NA, nrow = 1000, ncol = 7)

for (s in 1:1000){
  # Generation of {x} as ARMA process with p being 0.5, q being 0.3 and epsilon being generated as before
  x <- matrix(NA, nrow = t+m, ncol = 1)
  epsilon <- rnorm(t+m)
  
  for (i in 2:3000) {
    x[1] <- epsilon[1]
    x[i]= 0.5 * x[i-1] + epsilon[i] + 0.3 * epsilon[i-1]
  }
  
  # Generation of the ARFIMA coefficients
  psi <- matrix(NA, nrow = 3000, ncol=1)
  d <- 0.25
  
  for (j in 2:3000) {
    psi[1] <- 1 
    psi[j] <- ((j - 2 + d) / (j - 1)) * psi[j - 1]
  }
  
  # Generation of {y} 
  y <- matrix(NA, nrow = 2999, ncol = 1)
  
  for (tt in 1:3000) {
    ti <- 1:tt
    y[tt] <- sum(x[tt - ti + 1] * psi[ti])
  }
  
  
  
  
  # plugging different BW for both estimators
  for (i in 1:7){
    # LW estimation
    outputLW <- optimize(LW_R_d, y = y, interval = c(-1,1), maximum = F)$minimum
    d_LW_est[s,i] <- outputLW
    test_statLW[s,i] <- 2 * sqrt(M) * (outputLW - 0.25) #slide 177
   
    
    # PR estimation
    M <- m_alpha[i] #loopen counter
    ts_data <- y[1:M]
    I_j <- matrix(NA, ncol = 1, nrow = M)
    R_j <- matrix(NA, ncol = 1, nrow = M)
    
    for (j in 1:M){
      lambda <- 2 * pi * j / length(y)
      part1 <- (1/(2*pi*T)) 
      part2 <- sum(y * cos(lambda*(1:length(y))))
      part3 <- sum(y * sin(lambda*(1:length(y))))
      I_j[j] <- log(part1 * (part2^2 + part3^2))    #log(i_y)
      R_j[j] <- -log(4*sin(lambda/2)^2)  #is this sqared sine correct?
      R_j_mean <- mean(R_j)
      
      estimator_PR <- sum((R_j - R_j_mean) * I_j) / sum((R_j - R_j_mean)^2)
      d_PR_est[s,i] <- estimator_PR
      test_statPR[s,i] <- (sqrt(M) * (0.25 - estimator_PR)) / sqrt(pi^2/24)
    }
  }
  print(s)
}


# Coverage Probabilities estimated with our MC Simulation

cov_prob_LW_m1 <- sum(test_statLW[,1]<=1.96 & test_statLW[,1] >= -1.96)/1000
cov_prob_LW_m2 <- sum(test_statLW[,2]<=1.96 & test_statLW[,2] >= -1.96)/1000
cov_prob_LW_m3 <- sum(test_statLW[,3]<=1.96 & test_statLW[,3] >= -1.96)/1000
cov_prob_LW_m4 <- sum(test_statLW[,4]<=1.96 & test_statLW[,4] >= -1.96)/1000
cov_prob_LW_m5 <- sum(test_statLW[,5]<=1.96 & test_statLW[,5] >= -1.96)/1000
cov_prob_LW_m6 <- sum(test_statLW[,6]<=1.96 & test_statLW[,6] >= -1.96)/1000
cov_prob_LW_m7 <- sum(test_statLW[,7]<=1.96 & test_statLW[,7] >= -1.96)/1000

cov_prob_PR_m1 <- sum(test_statPR[,1]<=1.96 & test_statPR[,1] >= -1.96)/1000
cov_prob_PR_m2 <- sum(test_statPR[,2]<=1.96 & test_statPR[,2] >= -1.96)/1000
cov_prob_PR_m3 <- sum(test_statPR[,3]<=1.96 & test_statPR[,3] >= -1.96)/1000
cov_prob_PR_m4 <- sum(test_statPR[,4]<=1.96 & test_statPR[,4] >= -1.96)/1000
cov_prob_PR_m5 <- sum(test_statPR[,5]<=1.96 & test_statPR[,5] >= -1.96)/1000
cov_prob_PR_m6 <- sum(test_statPR[,6]<=1.96 & test_statPR[,6] >= -1.96)/1000
cov_prob_PR_m7 <- sum(test_statPR[,7]<=1.96 & test_statPR[,7] >= -1.96)/1000

### results

rbind(cov_prob_LW_m1, cov_prob_LW_m2, cov_prob_LW_m3, cov_prob_LW_m4, cov_prob_LW_m5, cov_prob_LW_m6, cov_prob_LW_m7)
# cov_prob_LW_m1 = 0.000
# cov_prob_LW_m2 = 0.931
# cov_prob_LW_m3 = 0.944
# cov_prob_LW_m4 = 0.972
# cov_prob_LW_m5 = 0.988
# cov_prob_LW_m6 = 0.700
# cov_prob_LW_m7 = 0.000

diag(var(d_LW_est))
# 0.000000000 0.057256227 0.018967798 0.006833253 0.002485172 0.001029985 0.000499923

rbind(cov_prob_PR_m1, cov_prob_PR_m2, cov_prob_PR_m3, cov_prob_PR_m4, cov_prob_PR_m5, cov_prob_PR_m6, cov_prob_PR_m7)
# cov_prob_PR_m1 = 0.680
# cov_prob_PR_m2 = 0.827
# cov_prob_PR_m3 = 0.885
# cov_prob_PR_m4 = 0.916
# cov_prob_PR_m5 = 0.926
# cov_prob_PR_m6 = 0.569
# cov_prob_PR_m7 = 0.000


diag(var(d_PR_est))
#0.3661428749 0.0781791093 0.0262250668 0.0102247053 0.0040218497 0.0017393581 0.0007378615



# the variance drops with the bandwidth 
# medium bw seems to be the best in terms of estimation and overall MSE

