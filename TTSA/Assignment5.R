## Assignment 5 ----

##Question 9 ----

install.packages("forecast")
library(forecast)

# r_t from Assignment 3

prediction <- matrix(rep(0,28),4,7)
rownames(prediction) <- c("p = 0:","p = 1:","p = 2:","p = 3:") 
colnames(prediction) <- c("1","2","3","4","5","6","  true forecast")
forecast_errors <- matrix(rep(0,28),4,7)
rownames(forecast_errors) <- c("p = 0:","p = 1:","p = 2:","p = 3:") 
colnames(forecast_errors) <- c("1","2","3","4","5","6","  true forecast")

for (p in 1:4){
  for (i in 1:7){
    # one-step ahead predicition for the last 6 periods as well as a true forecst for one period 
    fit <- arima(r_t[1:(100+i)], order = c((p-1),0,0))
    output <- predict(fit, n.ahead = 1)
    prediction[p, i] <- round(output$pred[1],5)
#    output <- forecast(fit, 1)
#    prediction[p,i] <- output$mean #does the same thing
    
    # fill error matrix (last column is irrelevant as there is only forecasted value)
    forecast_errors[p,i] <- round(r_t[101+i] - prediction[p, i],5)
  }
}
  

prediction
forecast_errors

rmse0 <- sqrt(1/6*sum(forecast_errors[1,1:6]^2))
rmse1 <- sqrt(1/6*sum(forecast_errors[2,1:6]^2))
rmse2 <- sqrt(1/6*sum(forecast_errors[3,1:6]^2))
rmse3 <- sqrt(1/6*sum(forecast_errors[4,1:6]^2))
RMSE <- cbind(c(rmse0,rmse1, rmse2, rmse3))

RMSE

## AR(1) produces the lowest RMSE of 0.003753149