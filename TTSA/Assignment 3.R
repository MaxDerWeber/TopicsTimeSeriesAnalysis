### TTSA Assignment 3 -- Question 6

# import and defining time series ----

BIP_season_adjusted <- read_delim("~/Desktop/TTSA DATA.txt", 
                      + "\t", escape_double = FALSE, col_types = cols(observation_date = col_date(format = "%Y-%m-%d")), 
                      + trim_ws = TRUE)

BIP_data <- BIP_season_adjusted$CPMNACSCAB1GQDE
BIP_ts <- ts(BIP_data, start = c(1991,1), frequency = 4)
plot(BIP_ts)

# generate the logarithmic quarterly gdp growth series r_t ----

log_BIP_ts <- log(BIP_ts)
length(BIP_season_adjusted$CPMNACSCAB1GQDE)
log_x_t <- log_BIP_ts[2:108]
log_x_tminus1 <- log_BIP_ts[1:107]

r_t1 <- log_x_t - log_x_tminus1
r_t <- ts(r_t1, start = c(1991,2), frequency = 4)
plot(r_t, type='l')
#View(r_t)

# fit the ARMA model using AIC ----
aic <- matrix(rep(0,16),4,4)
rownames(aic) <- c("AR0", "AR1", "AR2", "AR3")
colnames(aic) <- c("MA0", "MA1", "MA2", "MA3")

for (i in 1:4) {
  for (j in 1:4) {
    aic[i,j] <- AIC(arima(r_t, order = c(i-1,0,j-1)))
  }
}
# aic:

#MA0       MA1       MA2       MA3
#AR0 -720.8904 -725.8781 -724.4745 -724.5310
#AR1 -726.9301 -725.6936 -722.0653 -722.9290
#AR2 -725.4215 -722.9675 -721.7433 -724.0187
#AR3 -724.5110 -722.5836 -721.6223 -722.3127

min(aic) # -> ARMA(0,1) is the best according to AIC as it produces the lowest value (-726.9301)

# evaluate model ----
ar1ma0
ar1ma0 <- arima(r_t, order = c(1,0,0))
ar1ma0_residuals <- ts(ar1ma0$residuals,start = c(1991,2), end = c(2017,3), frequency = 4)
fitted_ar1ma0 <- r_t - ar1ma0_residuals

plot(r_t, col = 'red')
points(ar1ma0_residuals, type = 'l', col = 'blue')
points(fitted_ar1ma0, type = 'l', col = 'black')

acf(ar1ma0_residuals, type = "correlation") #the first value is negligable
## the residuals are free of serial correlation

