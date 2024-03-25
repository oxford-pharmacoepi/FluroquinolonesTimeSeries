library(dplyr)
library(readr)
library(lubridate)
library(forecast)
library(parallel)
library(ggplot2)
library(lmtest)


# Function to read and preprocess data
preprocess_data <- function(data, age_group) {
  data %>%
    filter(
      outcome_cohort_name == "fluroquinolones",
      analysis_interval == "months",
      denominator_sex == "Both",
      denominator_age_group == age_group
    ) %>%
    mutate(month = month(incidence_start_date)) %>%
    mutate(
      intrv = ifelse(incidence_start_date >= "2019-04-01", 1, 0),
      quarter = factor(quarter(incidence_start_date)),
      intervention_2015 = if_else(incidence_start_date >= as.Date("2015-01-01"), 1, 0),
      month = month(incidence_start_date)
    ) %>%
    mutate(quarter1 = as.factor(if_else(quarter == "1", 1, 0))) %>%
    mutate(quarter2 = as.factor(if_else(quarter == "2", 1, 0))) %>%
    mutate(quarter3 = as.factor(if_else(quarter == "3", 1, 0))) %>%
    mutate(quarter4 = as.factor(if_else(quarter == "4", 1, 0))) %>%
    filter(incidence_start_date!= "2019-04-01")

  #
}


model_parameters <- list()
fit_and_save_arima <- function(x, cdm_name) {

  age_group <- unique(x$age_group)

  x <- x[min(which(x$n_events != 0)):dim(x)[1], ]

  x <- x[!is.na(x$person_years), ]

  x[is.na(x$n_events), ]$n_events <- 5

  x[x$n_events == 0, ]$n_events <- 5

  step <- as.numeric(x$incidence_start_date >= as.Date("2019-04-01"))

  ramp <- append(rep(0, sum(step == 0)), seq(1, sum(step != 0), 1))

  x$time <- 1:dim(x)[1]


  ts_x <- ts(x$n_events, frequency = 12, start = c(
    year(min(x$incidence_start_date)),
    month(min(x$incidence_start_date))
  ))


  # model_arima <- auto.arima(ts_x, seasonal=TRUE, xreg=cbind(step, ramp),
  #                           max.d=1, max.D=1, stepwise=FALSE, trace=TRUE)

  model_arima <- auto.arima(ts_x,
    seasonal = TRUE, xreg = cbind(step, ramp),
    max.p = 15, max.q = 15,
    max.P = 2, max.Q = 2,
    max.order = 15,
    allowdrift = TRUE,
    allowmean = FALSE, approximation = FALSE,
    ic = "aic",
    max.d = 1, max.D = 1, stepwise = FALSE, trace = TRUE
  )

  checkresiduals(model_arima)
  p <- as.numeric(arimaorder(model_arima)[1])
  d <- as.numeric(arimaorder(model_arima)[2])
  q <- as.numeric(arimaorder(model_arima)[3])
  P <- as.numeric(arimaorder(model_arima)[4])
  D <- as.numeric(arimaorder(model_arima)[5])
  Q <- as.numeric(arimaorder(model_arima)[6])
  period <- as.numeric(arimaorder(model_arima)[7])

  aic_value <- model_arima$aic

  model_results_df <- data.frame(
    p = p, d = d, q = q,
    P = P, D = D, Q = Q,
    period = period,
    AIC = aic_value )

  write.csv(model_results_df, paste0("models_param/arima_parameter_", cdm_name, "_", age_group, ".csv"), row.names = FALSE)




  if (is.na(period)) {
    model_arima_nointrv <- Arima(window(ts_x, end = c(2019, 3)), order = c(p, d, q))
  } else {
    model_arima_nointrv <- Arima(window(ts_x, end = c(2019, 3)),
      order = c(p, d, q),
      seasonal = list(order = c(P, D, Q), period = period)
    )
  }

  coef_test_results_with_intrv <- matrix(coeftest(model_arima), ncol = 4)
  coef_df_with_intrv <- data.frame(
    Parameter = row.names(coeftest(model_arima)),
    coef_test_results_with_intrv,
    Model = "With Intervention"
  )
  names(coef_df_with_intrv)[2:5] <- c("Estimate", "Std. Error", "Statistics", "P-value")

  # Assuming 'model_arima_nointrv' is your fitted model for 'Without intervention'
  coef_test_results_without_intrv <- matrix(coeftest(model_arima_nointrv), ncol = 4)
  coef_df_without_intrv <- data.frame(
    Parameter = row.names(coeftest(model_arima_nointrv)),
    coef_test_results_without_intrv,
    Model = "Remove Intervention"
  )
  names(coef_df_without_intrv)[2:5] <- c("Estimate", "Std. Error", "Statistics", "P-value")

  # Combine the coefficient data frames
  coef_df_combined <- rbind(coef_df_with_intrv, coef_df_without_intrv)

  # Save the combined coefficients to CSV
  csv_filename <- paste0("models_coef/arima_coef_", cdm_name, "_", age_group, ".csv")
  write.csv(coef_df_combined, csv_filename, row.names = FALSE)


  fitted_nointrv <- fitted(model_arima_nointrv)
  res_se <- sd(resid(model_arima_nointrv))
  res_se_x <- sd(resid(model_arima_nointrv))
  # Assuming a 95% confidence interval
  alpha <- 0.05
  z_value <- qnorm(1 - alpha / 2)
  # Calculate the confidence intervals
  low_b <- fitted(model_arima_nointrv) - z_value * res_se
  up_b <- fitted(model_arima_nointrv) + z_value * res_se


  forecast_nointrv <- forecast::forecast(model_arima_nointrv,
    h = length(ts_x) - length(fitted_nointrv)
  )

  start_date <- min(x$incidence_start_date)
  end_date <- start_date %m+% months(length(fitted_nointrv) + length(forecast_nointrv$mean) - 1) # Adjust for the total length
  date_sequence <- seq(start_date, end_date, by = "month")

  # Step 2: Combine the series (fitted and forecasted values)
  combined_series <- c(fitted_nointrv, forecast_nointrv$mean)
  combined_95_lower <- c(low_b, forecast_nointrv$lower[, "95%"])
  combined_95_upper <- c(up_b, forecast_nointrv$upper[, "95%"])


  # Step 3: Create a Data Frame
  pred_nointrv <- data.frame(
    Date = date_sequence,
    Forecast = combined_series* 100000/x$person_years,
    Lower_95_CI = combined_95_lower* 100000/x$person_years,
    Upper_95_CI = combined_95_upper* 100000/x$person_years,
    Model = "Remove intervention"
  )

  fitted_intrv <- fitted(model_arima)
  res_se_intrv <- sd(resid(model_arima))
  # Assuming a 95% confidence interval
  alpha <- 0.05
  z_value <- qnorm(1 - alpha / 2)
  # Calculate the confidence intervals
  low_b_intrv <- fitted(model_arima) - z_value * res_se_intrv
  up_b_intrv <- fitted(model_arima) + z_value * res_se_intrv

  pred_intrv <- data.frame(
    Date = date_sequence,
    Forecast = as.numeric(fitted_intrv)* 100000/x$person_years,
    Lower_95_CI = as.numeric(low_b_intrv)* 100000/x$person_years,
    Upper_95_CI = as.numeric(up_b_intrv)* 100000/x$person_years,
    Model = "With intervention"
  )

  pred_combined <- rbind(pred_intrv, pred_nointrv)

  write.csv(pred_combined, paste0("arima_pred/arima_prediction_", cdm_name, "_", age_group, ".csv"), row.names = FALSE)



  model_res_arima <- residuals(model_arima)
  # Create and save ARIMA fitted, validated and predicted plot
  png_filename <- paste0("Diag/ARIMA_diag_", cdm_name, "_", age_group, ".png")
  png(png_filename)
  par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
  acf(as.vector(model_res_arima), main = "")
  pacf(as.vector(model_res_arima), main = "")

  title(paste0("ARIMA diagnostics ", cdm_name, " ", age_group), outer = TRUE)
  dev.off() # Close the device

  png_qq <- paste0("qq_", unique(x$cdm_name), "_", age_group, ".png")


  if (any(model_res_arima != 0)) {
    png(png_qq)
    par(mfrow = c(1, 1))
    qqnorm(model_res_arima, main = paste0("Q-Q Plot of ARIMA Residuals ", unique(x$cdm_name)))
    qqline(model_res_arima)
    dev.off() # Close the device
  }

  standardized_res <- scale(model_res_arima)

  sp <- paste0("scatterplot_", unique(x$cdm_name), "_", age_group, ".png")
  png(sp)
  print(
    ggplot(data.frame(Time = x$incidence_start_date, Residuals = as.vector(standardized_res)), aes(x = Time, y = Residuals)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "Time vs. Standardized Residuals", x = "Time", y = "Standardized Residuals")
  )
  dev.off()


}

data_young <- preprocess_data(incidenceTable, "19 to 59")
data_old <- preprocess_data(incidenceTable, "60 to 150")

# Combine and tag data with age_group before splitting by cdm_name
tag_and_split_data <- function(data, age_group) {
  data <- mutate(data, age_group = age_group) # Tagging data with age_group
  split(data, data$cdm_name) # Splitting by cdm_name
}

data_young_tagged <- tag_and_split_data(data_young, "19 - 59")
data_old_tagged <- tag_and_split_data(data_old, "60 - 150")

# Combine tagged data into a single list for processing


process_function <- function(name, list_element) {
  # Extract the first (and only) dataframe from the list element
  data_subset <- list_element

  # Call the analysis function
  fit_and_save_arima(data_subset, name)
}
# Setup parallel processing
cl1 <- makeCluster(detectCores() - 1)
clusterExport(cl1, varlist = c("fit_and_save_arima", "data_young_tagged", "process_function"), envir = environment())
clusterEvalQ(cl1, {
  library(forecast)
    library(dplyr)
    library(readr)
    library(vcd)
    library(tsModel)
    library(foreign)
    library("lmtest")
    library("Epi")
    library("splines")
    library(ggplot2)
    library(forecast)
    library(lubridate)
    library(pbs)
    library(parallel)
  })

# Iterate over combined_data names for parallel execution
parLapply(cl1, names(data_young_tagged), function(name) {
  process_function(name, data_young_tagged[[name]])
})
stopCluster(cl1)



cl2 <- makeCluster(detectCores() - 1)

clusterExport(cl2, varlist = c("fit_and_save_arima", "data_old_tagged", "process_function"), envir = environment())
clusterEvalQ(cl2, {
  library(forecast)
  library(dplyr)
  library(readr)
  library(vcd)
  library(tsModel)
  library(foreign)
  library("lmtest")
  library("Epi")
  library("splines")
  library(ggplot2)
  library(forecast)
  library(lubridate)
  library(pbs)
  library(parallel)
})
parLapply(cl2, names(data_old_tagged), function(name) {
  process_function(name, data_old_tagged[[name]])
})
# Clean up
stopCluster(cl2)
