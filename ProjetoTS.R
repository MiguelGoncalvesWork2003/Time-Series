library(astsa)
library(xts)
library(forecast)
library(fpp)
library(ggplot2)
library(lubridate)
library(urca)
library(tseries)
library(Quandl)
library(fable)
library(scales)
library(zoo)
library(dplyr)
library(tsbox)
ts <- read.csv(file.choose())
ts_clean <- ts[c('date','wind_united_states')]

ts_clean$date <- as.character(ts_clean$date)
ts_clean$date <- mdy(paste(ts_clean$date))
ts_clean$date <- as.Date(paste(ts_clean$date, "01", sep = "-"), format = "%Y-%m-%d")


ts_clean1 <- ts(ts_clean$wind_united_states, start = c(year(min(ts_clean$date)), month(min(ts_clean$date))), frequency = 12)


ts_clean1 <- window(ts_clean1, start = c(2001, 1), end = c(2022, 12))


plot(ts_clean1, main = "Subsetted Time Series (2001-2022)", xlab = "Year", ylab = "Wind United States Value")


#train set e test set


train_set <- window(ts_clean1, start = c(2001, 1), end = c(2020, 12))


test_set <- window(ts_clean1, start = c(2021, 1), end = c(2022, 12))


plot(train_set, main = "Training Set (2001-2020)", xlab = "Year", ylab = "Wind United States Value")


acf(train_set, main = "ACF of Training Set (2001-2020)")


pacf(train_set, main = "PACF of Training Set (2001-2020)")

#claramente tem trend e seasonality e é heterocedastic
#primeiro analise ao wind_united_states


par(mfrow = c(3, 1))
lambda = BoxCox.lambda(train_set)
lambda
plot(train_set, main = expression(train_set), ylab="")
plot(log(train_set), main = expression(log(train_set)), ylab = "")
plot(BoxCox(train_set,lambda), main = expression(BoxCox(train_set)), ylab = "")

par(mfrow = c(1, 1))
# Define the start year dynamically
#start_year <- start(train_set)[1] + (length(train_set) - length(wind)) / 12

# Create the ts object with adjusted start year
#wind <- ts(wind, start = c(start_year, 1), frequency = 12)

wind <- ts(BoxCox(train_set, lambda), start = start(train_set), frequency = frequency(train_set))



tsplot(wind, main = expression("BoxCox(wind)"))

nsdiffs(wind) #deu 1 entao D=1 no modelo multiplicativo

ndiffs(diff(wind,12)) #deu 1 entao d=1 no modelo multiplicativo

ndiffs(diff(diff(wind,12)))

wind_diff <- diff(diff(wind,12))


par(mfrow = c(2, 1))
plot(wind, main = expression("BoxCox(wind)"))
plot(wind_diff, main = expression(paste(Delta, Delta[12], "BoxCox(wind)")))




#Training and test set

#plots todos para começar a analisar o modelo SARIMA

maxlag <- 48
par(mfrow=c(3,2), mar=c(3,3,4,2))

plot(wind, main = expression("BoxCox(wind)"))
plot(wind_diff, main = expression(paste(Delta, Delta[12], "BoxCox(wind)")))

Acf(wind, type='correlation', lag=maxlag, ylab="", main=expression(paste("ACF for BoxCox(wind)")))
Acf(wind_diff, type='correlation', lag=maxlag, na.action=na.omit, ylab="", main=expression(paste("ACF for ", Delta, Delta[12], "BoxCox(wind)")))

Acf(wind, type='partial', lag=maxlag, na.action=na.omit, ylab="", main=expression(paste("PACF for BoxCox(wind)")))
Acf(wind_diff, type='partial', lag=maxlag, na.action=na.omit, ylab="", main=expression(paste("PACF for ", Delta,Delta[12], "BoxCox(wind)")))

#ACF e PACF que é supost analisar e daqui tirar o modelo
#acf2(wind_diff, max.lag=maxlag, main=expression(paste("ACF for ", Delta, Delta[12], "wind")))


ggtsdisplay(wind_diff, lag.max = 48, theme = theme_light(), main = "Stationary training data")

#--------------------------------------------------------------------------------------------------------

M1 = sarima(wind, 2, 1, 1, 3, 1, 0, 12)

M1 = sarima(wind, 0, 1, 1, 0, 1, 1, 12)
#mau pelo Ljung-box

M2 = sarima(wind, 2, 1, 0, 2, 1, 0, 12)
#nao é o melhor

M3 = sarima(wind, 2, 1, 1, 2, 1, 1, 12)
#ajsutando para sarima(wind, 2, 1, 1, 0, 1, 1, 12) igual ao M6
M3= sarima(wind, 2, 1, 1, 0, 1, 1, 12)
#bom

M4 = sarima(wind, 0, 1, 1, 2, 1, 0, 12)
#muito bom

M5 = sarima(wind, 2, 1, 0, 0, 1, 1, 12)
#nao é bom

M6 = sarima(wind, 2, 1, 1, 0, 1, 1, 12)
#bom

M7 = sarima(wind, 2, 1, 1, 3, 1, 0, 12)
#nao muito bom

M9= sarima(wind, 9, 1, 0, 0, 1, 1, 12)
#--------------------------------------------------------------------------------------------------------
#Forecast 12 meses para M4
# Subset the time series (exclude the last observation)

# Fit the ARIMA model on the training subset
M4.fit_plots <- Arima(train_set, order = c(0, 1, 1), seasonal = c(2, 1, 0))
M4.fct.v1 <- forecast(M4.fit_plots, h = 24, level = 95)

# Plot the forecast along with the observed data
autoplot(ts_c(Observed = ts_clean1, `Fixed training (24-step ahead forecast)` = M4.fct.v1$mean),
        ylab = NULL, size = 1, main = "Forecasts from ARIMA(0,1,1)(2,1,0)[12]") + scale_color_manual(values = c("black",
                                                                                                                "tomato")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),
                                                                                                                                                   legend.position = c(0.76, 0.2)) + autolayer(M4.fct.v1, showgap = F)

# Auxiliar: Obtaining the SE from forecasts
aux_M4 = sarima.for(train_set, n.ahead = 24, p = 0, d = 1, q = 1, P = 2, D = 1, Q = 0,S = 12, plot = F)
test_set
# Creating a forecast table for M4
fct.tab_M4 = data.frame(
  Month = as.yearmon(time(test_set)), 
  Observed = as.numeric(test_set),  # Converter test_set para numérico
  Forecast = round(M4.fct.v1$mean, 2), 
  Error = round((as.numeric(test_set) - M4.fct.v1$mean), 2), 
  Error_PCT = percent(
    (as.vector((as.numeric(test_set) - M4.fct.v1$mean) / as.numeric(test_set))), 
    accuracy = 0.01
  ), 
  Covered_by_CI = between(
    as.numeric(test_set), 
    (M4.fct.v1$mean - 1.96 * aux_M4$se), 
    (M4.fct.v1$mean + 1.96 * aux_M4$se)
  )
)

# Display the forecast table for M4
knitr::kable(
  fct.tab_M4, 
  digits = 2, 
  align = "lccccc", 
  caption = "Forecasts (24-months ahead) for M4"
)

#Accuracy Measures
accuracy_M4 = rbind(accuracy(M4.fct.v1$mean, test_set))
rownames(accuracy_M4) <- c("Fixed training (24-step ahead forecast)")

knitr::kable(accuracy_M4, digits = 4, align = "ccccccc", caption = "Accuracy measures (24-step ahead forecast) M4")




#Forecast 12 meses para M6
#Model on the training subset
M6.fit_plots <- Arima(train_set, order = c(2, 1, 1), seasonal = c(0, 1, 1))
M6.fct.v1 <- forecast(M6.fit_plots, h = 24, level = 95)


# Plot the forecast along with the observed data
autoplot(ts_c(Observed = ts_clean1, `Fixed training (24-step ahead forecast)` = M6.fct.v1$mean),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(2,1,1)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "tomato")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),
                                                                                                                                                    legend.position = c(0.76, 0.2)) + autolayer(M6.fct.v1, showgap = F)


# Auxiliar: Obtaining the SE from forecasts
aux_M6 = sarima.for(train_set, n.ahead = 25, p = 2, d = 1, q = 1, P = 0, D = 1, Q = 1,S = 12, plot = F)


# Criação da tabela de previsões para o modelo M6
fct.tab_M6 = data.frame(
  Month = as.yearmon(time(test_set)), 
  Observed = as.numeric(test_set), 
  Forecast = round(M6.fct.v1$mean, 2), 
  Error = round((as.numeric(test_set) - M6.fct.v1$mean), 2), 
  Error_PCT = scales::percent(
    (as.vector((as.numeric(test_set) - M6.fct.v1$mean) / as.numeric(test_set))), 
    accuracy = 0.01
  ), 
  Covered_by_CI = dplyr::between(
    as.numeric(test_set), 
    (M6.fct.v1$mean - 1.96 * aux_M6$se), 
    (M6.fct.v1$mean + 1.96 * aux_M6$se)
  )
)

# Exibição da tabela
knitr::kable(fct.tab_M6, digits = 2, align = "lccccc", caption = "Forecasts (24-months ahead for Model M6)")


#Accuracy Measures M6
accuracy_M6 = rbind(accuracy(M6.fct.v1$mean, test_set))
rownames(accuracy_M6) <- c("Fixed training (24-step ahead forecast)")

knitr::kable(accuracy_M6, digits = 4, align = "ccccccc", caption = "Accuracy measures (24-step ahead forecast) M6")


#Forecast 12 meses para M9
#Model on the training subset
M9.fit_plots <- Arima(train_set, order = c(2, 1, 1), seasonal = c(0, 1, 1))
M9.fct.v1 <- forecast(M9.fit_plots, h = 24, level = 95)


# Plot the forecast along with the observed data
autoplot(ts_c(Observed = ts_clean1, `Fixed training (24-step ahead forecast)` = M9.fct.v1$mean),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(9,1,0)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "tomato")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),
                                                                                                                                                    legend.position = c(0.76, 0.2)) + autolayer(M9.fct.v1, showgap = F)


# Auxiliar: Obtaining the SE from forecasts
aux_M9 = sarima.for(train_set, n.ahead = 25, p = 2, d = 1, q = 1, P = 0, D = 1, Q = 1,S = 12, plot = F)


# Criação da tabela de previsões para o modelo M6
fct.tab_M9 = data.frame(
  Month = as.yearmon(time(test_set)), 
  Observed = as.numeric(test_set), 
  Forecast = round(M9.fct.v1$mean, 2), 
  Error = round((as.numeric(test_set) - M6.fct.v1$mean), 2), 
  Error_PCT = scales::percent(
    (as.vector((as.numeric(test_set) - M6.fct.v1$mean) / as.numeric(test_set))), 
    accuracy = 0.01
  ), 
  Covered_by_CI = dplyr::between(
    as.numeric(test_set), 
    (M9.fct.v1$mean - 1.96 * aux_M9$se), 
    (M9.fct.v1$mean + 1.96 * aux_M9$se)
  )
)

# Exibição da tabela
knitr::kable(fct.tab_M9, digits = 2, align = "lccccc", caption = "Forecasts (24-months ahead for Model M6)")


#Accuracy Measures M9
accuracy_M9 = rbind(accuracy(M9.fct.v1$mean, test_set))
rownames(accuracy_M9) <- c("Fixed training (24-step ahead forecast)")

knitr::kable(accuracy_M9, digits = 4, align = "ccccccc", caption = "Accuracy measures (24-step ahead forecast) M9")






#--------------------------------------------------------------------------------------------------------
#4.2. Forecasting 12 months with 1-step ahead - Interaction

# converting the ts object into a DF
ts_df = ts_df(ts_clean1)

# forecasting 12 times 1-step ahead for M4
output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time <= as.Date("2020-12-01") %m+% months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, model = M4.fit_plots)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), value = tail(fitted(fit),
                                                                            1))
  
  output = bind_rows(output, aux)
}

# converting output to TS
M4.fct.v2 = ts_ts(output)


# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Fixed training (1-step ahead forecast)` = M4.fct.v2),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(0,1,1)(2,1,0)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "gold")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),
                                                                                                                                                  legend.position = c(0.76, 0.2))


# forecasting 12 times 1-step ahead for M6
output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time <= as.Date("2020-12-01") %m+% months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, model = M6.fit_plots)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), value = tail(fitted(fit),
                                                                            1))
  
  output = bind_rows(output, aux)
}

# converting output to TS
M6.fct.v2 = ts_ts(output)

# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Fixed training (1-step ahead forecast)` = M6.fct.v2),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(2,1,1)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "gold")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),
                                                                                                                                                  legend.position = c(0.76, 0.2))
#M9
ts_df = ts_df(ts_clean1)

# forecasting 12 times 1-step ahead for M4
output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time <= as.Date("2020-12-01") %m+% months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, model = M9.fit_plots)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), value = tail(fitted(fit),
                                                                            1))
  
  output = bind_rows(output, aux)
}

# converting output to TS
M9.fct.v2 = ts_ts(output)


# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Fixed training (1-step ahead forecast)` = M9.fct.v2),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(9,1,0)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "gold")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),
                                                                                                                                                  legend.position = c(0.76, 0.2))


#--------------------------------------------------------------------------------------------------------
#4.3. Forecasting 12 months with 1-step ahead - Expanding windows

# converting the ts object into a DF
ts_df = ts_df(ts_clean1)

# forecasting M1 specification with expanding windows
output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time <= as.Date("2020-12-01") %m+% months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, order = c(0, 1, 1), seasonal = c(2, 1, 0))
  
  M4.fct.expw = forecast(fit, h = 1, level = 95)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), model_specification = as.character(fit),
                  value = M4.fct.expw$mean)
  
  output = bind_rows(output, aux)
}

# converting output to TS
M4.fct.v3 = ts_ts(output)


# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Expanding training (1-step ahead forecast)` = M4.fct.v3),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(0,1,1)(2,1,0)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "purple")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),legend.position = c(0.76, 0.2))

#para M6
# converting the ts object into a DF
ts_df = ts_df(ts_clean1)

# forecasting M1 specification with expanding windows
output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time <= as.Date("2020-12-01") %m+% months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, order = c(2, 1, 1), seasonal = c(0, 1, 1))
  
  M6.fct.expw = forecast(fit, h = 1, level = 95)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), model_specification = as.character(fit),
                  value = M6.fct.expw$mean)
  
  output = bind_rows(output, aux)
}

# converting output to TS
M6.fct.v3 = ts_ts(output)


# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Expanding training (1-step ahead forecast)` = M6.fct.v3),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(2,1,1)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "purple")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank())

#para M6
# converting the ts object into a DF
ts_df = ts_df(ts_clean1)

# forecasting M1 specification with expanding windows
output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time <= as.Date("2020-12-01") %m+% months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, order = c(2, 1, 1), seasonal = c(0, 1, 1))
  
  M6.fct.expw = forecast(fit, h = 1, level = 95)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), model_specification = as.character(fit),
                  value = M6.fct.expw$mean)
  
  output = bind_rows(output, aux)
}

# converting output to TS
M6.fct.v3 = ts_ts(output)


# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Expanding training (1-step ahead forecast)` = M6.fct.v3),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(2,1,1)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "purple")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank())



#para M9
# converting the ts object into a DF
ts_df = ts_df(ts_clean1)

# forecasting M1 specification with expanding windows
output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time <= as.Date("2020-12-01") %m+% months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, order = c(9, 1, 0), seasonal = c(0, 1, 1))
  
  M9.fct.expw = forecast(fit, h = 1, level = 95)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), model_specification = as.character(fit),
                  value = M9.fct.expw$mean)
  
  output = bind_rows(output, aux)
}

# converting output to TS
M9.fct.v3 = ts_ts(output)


# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Expanding training (1-step ahead forecast)` = M9.fct.v3),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(9,1,0)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "purple")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank())


#--------------------------------------------------------------------------------------------------------
#4.4
#M4
ts_df = ts_df(ts_clean1)


output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time >= as.Date("2000-01-01") %m+% months(i) & time <= as.Date("2020-12-01") %m+%
             months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, order = c(0, 1, 1), seasonal = c(2, 1, 0))
  
  M4.fct.recw = forecast(fit, h = 1, level = 95)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), model_specification = as.character(fit),
                  value = M4.fct.recw$mean)
  
  output = bind_rows(output, aux)
}

# converting output to TS
M4.fct.v4 = ts_ts(output)

# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Recursive training (1-step ahead forecast)` = M4.fct.v4),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(0,1,1)(2,1,0)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),
                                                                                                                                                             legend.position = c(0.76, 0.2))

#M6
ts_df = ts_df(ts_clean1)


output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time >= as.Date("2000-01-01") %m+% months(i) & time <= as.Date("2020-12-01") %m+%
             months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, order = c(2, 1, 1), seasonal = c(0, 1, 1))
  
  M6.fct.recw = forecast(fit, h = 1, level = 95)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), model_specification = as.character(fit),
                  value = M6.fct.recw$mean)
  
  output = bind_rows(output, aux)
}

# converting output to TS
M6.fct.v4 = ts_ts(output)

# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Recursive training (1-step ahead forecast)` = M6.fct.v4),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(2,1,1)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),
                                                                                                                                                             legend.position = c(0.76, 0.2))
#M9
ts_df = ts_df(ts_clean1)


output = data.frame()
for (i in 0:23) {
  df_select = ts_df %>%
    filter(time >= as.Date("2000-01-01") %m+% months(i) & time <= as.Date("2020-12-01") %m+%
             months(i))
  
  ts_select = ts(df_select$value, freq = 12, end = c(year(tail(df_select$time,
                                                               1)), month(tail(df_select$time, 1))))
  
  fit = Arima(ts_select, order = c(9, 1, 0), seasonal = c(0, 1, 1))
  
  M9.fct.recw = forecast(fit, h = 1, level = 95)
  
  aux = bind_cols(time = as.Date("2021-01-01") %m+% months(i), model_specification = as.character(fit),
                  value = M9.fct.recw$mean)
  
  output = bind_rows(output, aux)
}

# converting output to TS
M9.fct.v4 = ts_ts(output)

# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Recursive training (1-step ahead forecast)` = M9.fct.v4),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(9,1,0)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                                                 "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(), legend.background = element_blank(),
                                                                                                                                                             legend.position = c(0.76, 0.2))
#plotting all forecasts
#M4
# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Fixed training (24-step ahead forecast)` = M4.fct.v1$mean,
              `Fixed training (1-step ahead forecast)` = M4.fct.v2, `Expanding training (1-step ahead forecast)` = M4.fct.v3,
              `Recursive training (1-step ahead forecast)` = M4.fct.v4), ylab = NULL, size = 1,
         main = "Forecasts from ARIMA(0,1,1)(2,1,0)[12]") + scale_color_manual(values = c("black",
                                                                                          "tomato", "gold", "purple", "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(),
                                                                                                                                                                  legend.background = element_blank(), legend.position = c(0.76, 0.2))
# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Fixed training (24-step ahead forecast)` = M6.fct.v1$mean,
              `Fixed training (1-step ahead forecast)` = M6.fct.v2, `Expanding training (1-step ahead forecast)` = M6.fct.v3,
              `Recursive training (1-step ahead forecast)` = M6.fct.v4), ylab = NULL, size = 1,
         main = "Forecasts from ARIMA(2,1,1)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                          "tomato", "gold", "purple", "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(),
                                                                                                                                                                  legend.background = element_blank(), legend.position = c(0.76, 0.2))
# forecasts plot
autoplot(ts_c(Observed = ts_clean1, `Fixed training (24-step ahead forecast)` = M9.fct.v1$mean,
              `Fixed training (1-step ahead forecast)` = M9.fct.v2, `Expanding training (1-step ahead forecast)` = M9.fct.v3,
              `Recursive training (1-step ahead forecast)` = M9.fct.v4), ylab = NULL, size = 1,
         main = "Forecasts from ARIMA(9,1,0)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                          "tomato", "gold", "purple", "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(),
                                                                                                                                                                  legend.background = element_blank(), legend.position = c(0.76, 0.2))
#5 Assessing quality
#M4

autoplot(ts_c(Observed = test_set, `Fixed training (24-step ahead forecast)` = M4.fct.v1$mean,
              `Fixed training (1-step ahead forecast)` = M4.fct.v2, `Expanding training (1-step ahead forecast)` = M4.fct.v3,
              `Recursive training (1-step ahead forecast)` = M4.fct.v4), ylab = NULL, size = 1,
         main = "Forecasts from ARIMA(0,1,1)(2,1,0)[12]") + scale_color_manual(values = c("black",
                                                                                          "tomato", "gold", "purple", "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(),
                                                                                                                                                                  legend.background = element_blank(), legend.position = c(0.76, 0.8))
autoplot(ts_c(Observed = test_set, `Fixed training (24-step ahead forecast)` = M6.fct.v1$mean,
              `Fixed training (1-step ahead forecast)` = M6.fct.v2, `Expanding training (1-step ahead forecast)` = M6.fct.v3,
              `Recursive training (1-step ahead forecast)` = M6.fct.v4), ylab = NULL, size = 1,
         main = "Forecasts from ARIMA(2,1,1)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                          "tomato", "gold", "purple", "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(),
                                                                                                                                                                  legend.background = element_blank(), legend.position = c(0.76, 0.8))
autoplot(ts_c(Observed = test_set, `Fixed training (24-step ahead forecast)` = M9.fct.v1$mean,
              `Fixed training (1-step ahead forecast)` = M9.fct.v2, `Expanding training (1-step ahead forecast)` = M9.fct.v3,
              `Recursive training (1-step ahead forecast)` = M9.fct.v4), ylab = NULL, size = 1,
         main = "Forecasts from ARIMA(9,1,0)(0,1,1)[12]") + scale_color_manual(values = c("black",
                                                                                          "tomato", "gold", "purple", "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(),
                                                                                                                                                                  legend.background = element_blank(), legend.position = c(0.76, 0.8))
print(M6.fct.v1$mean)
print(M6.fct.v2)
print(M6.fct.v3)
print(M6.fct.v4)
#nao conseguimos ver o 3o porque o 4o é quase identico

#Plotting the forecast errors
autoplot(ts_c(`Fixed training (24-step ahead forecast)` = test_set - M4.fct.v1$mean,
              `Fixed training (1-step ahead forecast)` = test_set - M4.fct.v2, `Expanding training (1-step ahead forecast)` = test_set -
                M4.fct.v3, `Recursive training (1-step ahead forecast)` = test_set - M4.fct.v4),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(0,1,1)(2,1,0)[12]") + scale_color_manual(values = c("tomato",
                                                                                                                 "gold", "purple", "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(),
                                                                                                                                                                               legend.background = element_blank(), legend.position = c(0.3, 0.8)) + geom_hline(yintercept = 0,
                                                                                                                                                                                                                                                                color = "black", size = 1)  
autoplot(ts_c(`Fixed training (24-step ahead forecast)` = test_set - M6.fct.v1$mean,
              `Fixed training (1-step ahead forecast)` = test_set - M6.fct.v2, `Expanding training (1-step ahead forecast)` = test_set -
                M6.fct.v3, `Recursive training (1-step ahead forecast)` = test_set - M6.fct.v4),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(2,1,1)(0,1,1)[12]") + scale_color_manual(values = c("tomato",
                                                                                                                 "gold", "purple", "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(),
                                                                                                                                                                               legend.background = element_blank(), legend.position = c(0.3, 0.8)) + geom_hline(yintercept = 0,
                                                                                                                                                                                                                                                                color = "black", size = 1)  
autoplot(ts_c(`Fixed training (24-step ahead forecast)` = test_set - M9.fct.v1$mean,
              `Fixed training (1-step ahead forecast)` = test_set - M9.fct.v2, `Expanding training (1-step ahead forecast)` = test_set -
                M9.fct.v3, `Recursive training (1-step ahead forecast)` = test_set - M9.fct.v4),
         ylab = NULL, size = 1, main = "Forecasts from ARIMA(9,1,0)(0,1,1)[12]") + scale_color_manual(values = c("tomato",
                                                                                                                 "gold", "purple", "darkolivegreen4")) + theme_light() + theme(legend.title = element_blank(),
                                                                                                                                                                               legend.background = element_blank(), legend.position = c(0.3, 0.8)) + geom_hline(yintercept = 0,
                                                                                                                                                                                                                                                                color = "black", size = 1)  
accuracy = rbind(accuracy(M4.fct.v1$mean, test_set), accuracy(M4.fct.v2, test_set),
                 accuracy(M4.fct.v3, test_set), accuracy(M4.fct.v4, test_set))
rownames(accuracy) <- c("Fixed training (12-step ahead forecast)", "Fixed training (1-step ahead forecast)",
                        "Expanding training (1-step ahead forecast)", "Recursive training (1-step ahead forecast)")

knitr::kable(accuracy, digits = 4, align = "ccccccc", caption = "Accuracy measures (24-months forecasts)")

accuracy = rbind(accuracy(M6.fct.v1$mean, test_set), accuracy(M6.fct.v2, test_set),
                 accuracy(M6.fct.v3, test_set), accuracy(M6.fct.v4, test_set))
rownames(accuracy) <- c("Fixed training (12-step ahead forecast)", "Fixed training (1-step ahead forecast)",
                        "Expanding training (1-step ahead forecast)", "Recursive training (1-step ahead forecast)")

knitr::kable(accuracy, digits = 4, align = "ccccccc", caption = "Accuracy measures (24-months forecasts)")

accuracy = rbind(accuracy(M9.fct.v1$mean, test_set), accuracy(M9.fct.v2, test_set),
                 accuracy(M9.fct.v3, test_set), accuracy(M9.fct.v4, test_set))
rownames(accuracy) <- c("Fixed training (12-step ahead forecast)", "Fixed training (1-step ahead forecast)",
                        "Expanding training (1-step ahead forecast)", "Recursive training (1-step ahead forecast)")

knitr::kable(accuracy, digits = 4, align = "ccccccc", caption = "Accuracy measures (24-months forecasts)")

