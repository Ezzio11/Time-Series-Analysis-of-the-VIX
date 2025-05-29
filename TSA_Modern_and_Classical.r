################################################################################
### PACKAGE MANAGEMENT

# Install required packages 
install.packages(c("tseries", "forecast", "lmtest", "seastests", 
                   "readr", "zoo", "xts", "ggplot2", "bizdays"))

# Load required libraries
library(tseries)      # Time series analysis
library(forecast)     # Forecasting functions
library(lmtest)       # Statistical tests
library(seastests)    # Seasonality tests
library(readr)        # Data reading
library(zoo)          # Time series handling
library(xts)          # Extended time series
library(ggplot2)      # Data visualization
library(bizdays)      # Business days handling 

# Model comparison function
calc_metrics <- function(actual, predicted) {
  mse <- mean((actual - predicted)^2)
  mae <- mean(abs(actual - predicted))
  rmse <- sqrt(mse)
  mape <- mean(abs(actual - predicted)/actual) * 100
  return(c(MSE = mse, MAE = mae, RMSE = rmse, MAPE = mape))
}

################################################################################
### DATA LOADING AND PREPARATION

# Load dataset
VIXCLS <- read_csv("VIXCLS.csv", 
                   col_types = cols(
                     observation_date = col_date(format = "%d-%m-%y"), 
                     VIXCLS = col_number()
                   ))

# View data structure
View(VIXCLS)
attach(VIXCLS)

# Create time series vector
VIXCLS_val <- VIXCLS$VIXCLS

################################################################################
### TIME SERIES INITIALIZATION

# Create daily time series (252 trading days/year)
VIXCLS_ts <- ts(VIXCLS_val, frequency = 252, start = c(1990, 1))

# Initial plot
ts.plot(VIXCLS_ts, 
        main = "VIX Time Series", 
        ylab = "Daily VIX Value", 
        xlab = "Days")

################################################################################
### SEASONALITY AND TRANSFORMATION

# Check for seasonality
seastests::kw(VIXCLS_ts, freq = 252)  # p-value < 0.01 indicates seasonality

# Log-transform to stabilize variance
VIXCLS_ts <- log(VIXCLS_ts)
ts.plot(VIXCLS_ts, 
        main = "Log-Transformed VIX", 
        ylab = "log(VIX)", 
        xlab = "Days")

################################################################################
### ANNUAL AGGREGATION

# Create annual aggregated series
Annual_VIX_ts <- aggregate.ts(VIXCLS_ts, nfrequency = 1, FUN = mean)
ts.plot(Annual_VIX_ts, 
        main = "Annual VIX Averages", 
        ylab = "Annual Mean VIX (log)", 
        xlab = "Years")

################################################################################
### STATIONARITY ANALYSIS

# ACF and PACF plots
acf(Annual_VIX_ts, lag.max = 10, main = "ACF of Annual VIX")
pacf(Annual_VIX_ts, lag.max = 10, main = "PACF of Annual VIX")

################################################################################
### ARIMA MODELING

# Model 1: ARIMA(1,0,1)
Model1 <- arima(Annual_VIX_ts, order = c(1,0,1), method = "ML")
summary(Model1)
coeftest(Model1)
Residuals <- residuals(Model1)
acf(Residuals)
pacf(Residuals)
Box.test(Residuals, type = "Ljung", lag = 10, fitdf = 2)
Box.test(Residuals, type = "Box-Pierce", lag = 10, fitdf = 2)
accuracy(Model1)
f1 <- forecast(Model1, h = 10)
print(f1)
plot(f1)

# Model 2: ARIMA(1,0,0)
Model2 <- arima(Annual_VIX_ts, order = c(1,0,0), method = "ML")
summary(Model2)
coeftest(Model2)
Residuals <- residuals(Model2)
acf(Residuals)
pacf(Residuals)
Box.test(Residuals, type = "Ljung", lag = 10, fitdf = 2)
Box.test(Residuals, type = "Box-Pierce", lag = 10, fitdf = 2)
accuracy(Model2)
f2 <- forecast(Model2, h = 10)
print(f2)
plot(f2)

# Automated ARIMA selection
auto.arima(Annual_VIX_ts, trace = TRUE)

################################################################################
### TIME SERIES DECOMPOSITION

# Re-create daily time series for decomposition
vix_daily_ts <- ts(VIXCLS$VIXCLS, frequency = 252, start = c(1990, 1))

# Additive decomposition
decomp <- decompose(vix_daily_ts, type = "additive")
plot(decomp)

# Component extraction
trend_comp <- decomp$trend
seasonal_comp <- decomp$seasonal

# Detrended series
detrended_vix <- vix_daily_ts - trend_comp
plot(detrended_vix, main = "Detrended VIX Series")

# Deseasonalized series
deseason_vix <- vix_daily_ts - seasonal_comp
plot(deseason_vix, main = "Deseasonalized VIX Series")

# Moving average smoothing
ma_smoothed <- ma(vix_daily_ts, 5, centre = TRUE)

################################################################################
### TREND ANALYSIS

# Prepare clean data for regression
time_index <- seq_along(deseason_vix)
valid_points <- !is.na(deseason_vix)
clean_vix <- deseason_vix[valid_points]
clean_time <- time_index[valid_points]

# Linear trend model
linear_model <- lm(clean_vix ~ clean_time)
summary(linear_model)

# Quadratic trend model
quad_model <- lm(clean_vix ~ clean_time + I(clean_time^2))
summary(quad_model)

# Cubic trend model
cubic_model <- lm(clean_vix ~ clean_time + I(clean_time^2) + I(clean_time^3))
summary(cubic_model)

# Calculate metrics
metrics_linear <- calc_metrics(clean_vix, predict(linear_model))
metrics_quad <- calc_metrics(clean_vix, predict(quad_model))
metrics_cubic <- calc_metrics(clean_vix, predict(cubic_model))

metrics <- rbind(Linear = metrics_linear,
                 Quadratic = metrics_quad,
                 Cubic = metrics_cubic)

print(metrics)

################################################################################
### FORECASTING WITH CUBIC TREND MODEL (BEST MODEL)

# Create business calendar that excludes weekends
create.calendar("MyCal", weekdays = c("saturday", "sunday"))

# Last known date
last_date <- max(observation_date)

n_forecast <- 252  # number of future trading periods

# Generate future trading dates starting from the next business day
future_dates <- bizseq(from = last_date + 1, to = last_date + 365, cal = "MyCal")[1:n_forecast]

# Initialize forecast list
forecast_list <- vector("list", n_forecast)

# Loop to generate forecasts
for (i in 1:n_forecast) {
  next_period <- data.frame(clean_time = length(clean_time) + i)
  forecast_i <- predict(cubic_model, newdata = next_period, interval = "predict")
  final_forecast_i <- forecast_i + seasonal_comp[(i - 1) %% length(seasonal_comp) + 1]
  
  forecast_list[[i]] <- data.frame(
    date = future_dates[i],
    fit = final_forecast_i[,"fit"],
    lwr = final_forecast_i[,"lwr"],
    upr = final_forecast_i[,"upr"]
  )
}

# Combine list into a single data frame
forecast_df <- do.call(rbind, forecast_list)

# Print the forecasts
print(forecast_df)

# Plot original vs prediction
original_df <- data.frame(
  date = observation_date,
  vix = clean_vix
)

ggplot() +
  geom_line(data = original_df, aes(x = date, y = vix), color = "black", alpha = 0.6) +
  geom_ribbon(data = forecast_df, aes(x = date, ymin = lwr, ymax = upr), fill = "lightblue", alpha = 0.4) +
  geom_line(data = forecast_df, aes(x = date, y = fit), color = "blue", size = 1) +
  labs(
    title = "Original VIX + Forecast (Cubic Model)",
    x = "Date",
    y = "VIX"
  ) +
  theme_minimal()
################################################################################
### FINAL: YEARLY AGGREGATION & PLOTTING

# Modify forecast horizon to 10 years (2520 trading days)
n_forecast <- 2520

# Generate future trading dates for 10 years
end_date <- offset(last_date + 1, n_forecast - 1, "MyCal")
future_dates <- bizseq(from = last_date + 1, 
                       to = end_date, 
                       cal = "MyCal")

# Initialize forecast list
forecast_list <- vector("list", n_forecast)

# Loop to generate 2520 daily forecasts
for (i in 1:n_forecast) {
  next_period <- data.frame(clean_time = length(clean_time) + i)
  forecast_i <- predict(cubic_model, newdata = next_period, interval = "predict")
  
  # Seasonal adjustment using modulo for yearly cycle (252 trading days/year)
  seasonal_adj <- seasonal_comp[(i - 1) %% 252 + 1]
  final_forecast_i <- forecast_i + seasonal_adj
  
  forecast_list[[i]] <- data.frame(
    date = future_dates[i],
    fit = final_forecast_i[,"fit"],
    lwr = final_forecast_i[,"lwr"],
    upr = final_forecast_i[,"upr"]
  )
}

# Combine list into a single data frame
forecast_df <- do.call(rbind, forecast_list)

# --------------------------------------------------------------------------------
# Prepare historical yearly data (ensure this is done BEFORE any plot)
historical_yearly <- data.frame(
  year = time(Annual_VIX_ts),
  vix = exp(Annual_VIX_ts)  # Back-transform from log scale
)

# --------------------------------------------------------------------------------
# Group forecasts by year
forecast_df$year <- format(forecast_df$date, "%Y")

# Aggregate yearly forecast means
yearly_forecast <- aggregate(cbind(fit, lwr, upr) ~ year, 
                             data = forecast_df, 
                             FUN = mean)

# Convert year to numeric for plotting
yearly_forecast$year <- as.numeric(yearly_forecast$year)

# Round for cleaner table output
yearly_forecast[, 2:4] <- round(yearly_forecast[, 2:4], 2)

# --------------------------------------------------------------------------------
# Print tables
cat("Historical Yearly Averages:\n")
print(historical_yearly)

cat("\n\nForecasted Yearly Averages:\n")
print(yearly_forecast)

# --------------------------------------------------------------------------------
# Plot historical vs forecasted yearly averages
ggplot() +
  geom_line(data = historical_yearly, aes(x = year, y = vix), 
            color = "darkblue", size = 1) +
  geom_line(data = yearly_forecast, aes(x = year, y = fit), 
            color = "red", size = 1) +
  geom_ribbon(data = yearly_forecast, 
              aes(x = year, ymin = lwr, ymax = upr), 
              fill = "orange", alpha = 0.3) +
  labs(title = "Historical vs Forecasted Yearly VIX Averages",
       x = "Year", y = "VIX") +
  theme_minimal()
