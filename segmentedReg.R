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
library(here)
library(parallel)

incidenceTable <- read_csv(here("incidenceTable (4).csv"))

model_parameters <- list()
fit_and_save_segreg <- function(x, cdm_name){
  age_group <- unique(x$age_group)


  x <- x[min(which(x$n_events_sum != 0)):dim(x)[1],]

  x <- x[!is.na(x$person_years_sum),]

  x[is.na(x$n_events_sum),]$n_events_sum <- 5

  x[x$n_events_sum==0,]$n_events_sum<-5

  x$time <- 1:dim(x)[1]


  # previous models tested
  # model <- glm(n_events ~ offset(log(person_years)) + intervention_2015 +
  #                intrv*time + quarter4, family=quasipoisson, data=x)


  # Model with lagged term, spline, intervention and intervention*time
  model_with_intrv <- glm(n_events_sum ~ offset(log(person_years_sum)) + lag(n_events_sum) +
                            intrv*time, family=quasipoisson, data=x)
  model_with_intrv <- update(model_with_intrv, . ~ . + pbs(as.numeric(x$quarter), df=4))


  x_without_intrv <- x %>% dplyr::filter(incidence_start_date <= as.Date("2019-03-01"))

  model_without_intrv <- glm(n_events_sum ~ offset(log(person_years_sum)) + lag(n_events_sum) +
                               time + pbs(as.numeric(quarter), df=4), family=quasipoisson, data=x_without_intrv)

  # Test for simpler model: Model with only intervention and intervention*time
  model_simple <- glm(n_events_sum ~ offset(log(person_years_sum)) + intrv*time, family=quasipoisson, data=x)
  # Calculate the critical value for the 95% confidence interval
  critical_value <- qnorm(0.975)



  # Calculate predictions and SEs for the model with intervention
  pred_simple <- predict(model_simple, type="response", se.fit = TRUE)
  # Ensure element-wise operation for the calculation of CI
  ci_simple_lower <- pred_simple$fit - critical_value * pred_simple$se.fit
  ci_simple_upper <- pred_simple$fit + critical_value * pred_simple$se.fit



  # Calculate predictions and SEs for the model with intervention
  pred_with_intrv <- predict(model_with_intrv, type="response", se.fit = TRUE)
  # Ensure element-wise operation for the calculation of CI
  ci_with_intrv_lower <- pred_with_intrv$fit - critical_value * pred_with_intrv$se.fit
  ci_with_intrv_upper <- pred_with_intrv$fit + critical_value * pred_with_intrv$se.fit

  # Calculate predictions and SEs for the model without intervention

  pred_without_intrv <- predict(model_without_intrv, type="response", se.fit = TRUE, newdata = x %>% dplyr::select(c("person_years_sum",
                                                                                                              "n_events_sum",
                                                                                                              "quarter",
                                                                                                              "time")))
  ci_without_intrv_lower <- pred_without_intrv$fit - critical_value * pred_without_intrv$se.fit
  ci_without_intrv_upper <- pred_without_intrv$fit + critical_value * pred_without_intrv$se.fit


  # Combine forecasts with CIs
  forecasts_simple <- data.frame(
    Date = x$incidence_start_date,
    Forecast = c(pred_simple$fit) * 100000/x$person_years_sum,  #lagged term included so first term no pred
    Lower_95_CI = c(ci_simple_lower) * 100000/x$person_years_sum,
    Upper_95_CI = c(ci_simple_upper) * 100000/x$person_years_sum,
    Model = "Step and slope only",
    age_group = age_group,
    cdm_name = cdm_name
  )

  # Combine forecasts with CIs
  forecasts_with_intrv <- data.frame(
    Date = x$incidence_start_date,
    Forecast = c(NA, pred_with_intrv$fit) * 100000/x$person_years_sum,  #lagged term included so first term no pred
    Lower_95_CI = c(NA, ci_with_intrv_lower) * 100000/x$person_years_sum,
    Upper_95_CI = c(NA, ci_with_intrv_upper) * 100000/x$person_years_sum,
    Model = "With Intervention",
    age_group = age_group,
    cdm_name = cdm_name
  )

  forecasts_without_intrv <- data.frame(
    Date = x$incidence_start_date,
    Forecast = pred_without_intrv$fit * 100000/x$person_years_sum,#lagged term included so remove first obs, as no pred
    Lower_95_CI = ci_without_intrv_lower * 100000/x$person_years_sum,
    Upper_95_CI = ci_without_intrv_upper * 100000/x$person_years_sum,
    Model = "Remove Intervention",
    age_group = age_group,
    cdm_name = cdm_name
  )

  forecasts <- rbind(forecasts_with_intrv, forecasts_without_intrv, forecasts_simple)
  write.csv(forecasts, paste0("segreg_pred/segreg_prediction_", cdm_name, "_",
                              age_group, ".csv"), row.names = FALSE)

  # Extract and combine diagnostics for both models, including coefficients and p-values
  diagnostics <- data.frame(
    Model = c(rep("With Intervention", nrow(summary(model_with_intrv)$coefficients)),
              rep("Remove Intervention", nrow(summary(model_without_intrv)$coefficients))),
    Parameter = c(rownames(summary(model_with_intrv)$coefficients), rownames(summary(model_without_intrv)$coefficients)),
    Estimate = c(summary(model_with_intrv)$coefficients[, "Estimate"], summary(model_without_intrv)$coefficients[, "Estimate"]),
    "Std. Error" = c(summary(model_with_intrv)$coefficients[, "Std. Error"], summary(model_without_intrv)$coefficients[, "Std. Error"]),
    Statistics = c(summary(model_with_intrv)$coefficients[, "t value"], summary(model_without_intrv)$coefficients[, "t value"]),
    "P-value" = c(summary(model_with_intrv)$coefficients[, "Pr(>|t|)"], summary(model_without_intrv)$coefficients[, "Pr(>|t|)"])
  )

  diagnostics$deviance <- c(rep(deviance(model_with_intrv),
                           nrow(summary(model_with_intrv)$coefficients)),
                       rep(deviance(model_without_intrv),
                           nrow(summary(model_without_intrv)$coefficients)))

  csv_filename <- paste0("models_coef/segreg_coef_", cdm_name, "_", age_group, ".csv")
  write.csv(diagnostics, csv_filename, row.names = FALSE)

  residuals <- resid(model_with_intrv, type = "deviance")
  png_filename <- paste0("Diag/SegReg_diag_", cdm_name, "_", age_group, ".png")
  png(png_filename)
  par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
  acf(residuals, main = "")
  pacf(residuals, main = "")
  title(paste0("Segemented regression diagnostics ", cdm_name, " ", age_group), outer = TRUE)
  dev.off() # Close the device


  # save(residuals, file = paste0("residuals_", cdm_name, "_", age_group, ".RData"))
  residuals_data <- data.frame(incidence_start_date = x$incidence_start_date, residuals = c(NA, residuals))

  # Save this dataframe
  save(residuals_data, file = paste0("residuals_with_dates_", cdm_name, "_", age_group, ".RData"))




}


preprocess_data_seg <- function(data, age_group) {
  data <- data %>%
    filter(outcome_cohort_name=="fluroquinolones") %>%
    filter(analysis_interval=="months")%>%
    mutate(quarter=quarter(incidence_start_date))%>%
    mutate(year=year(incidence_start_date))%>%
    filter(denominator_sex=="Both")%>%
    filter(denominator_age_group == age_group) %>%
    group_by(year, quarter, cdm_name)%>%
    summarise(
      n_events_sum = sum(n_events),
      person_years_sum = sum(person_years),
      .groups = 'drop'
    ) %>%
    mutate(incidence_start_date = ymd(paste(year, quarter * 3 - 2, "01", sep = "-")))%>%
    arrange(incidence_start_date)%>%
    mutate(time=row_number())%>%
    mutate(intrv = ifelse(incidence_start_date >= "2019-04-01", 1, 0)) %>%
    filter(incidence_start_date!= as.Date("2019-04-01"))

  return(data)
}
data_young <- preprocess_data_seg(incidenceTable, "19 to 59")
data_old <- preprocess_data_seg(incidenceTable, "60 to 150")

# Combine and tag data with age_group before splitting by cdm_name
tag_and_split_data <- function(data, age_group) {
  data <- mutate(data, age_group = age_group) # Tagging data with age_group
  split(data, data$cdm_name) # Splitting by cdm_name
}

data_young_tagged_segreg <- tag_and_split_data(data_young, "19 - 59")
data_old_tagged_segreg <- tag_and_split_data(data_old, "60 - 150")

# Combine tagged data into a single list for processing


process_function_segreg <- function(name, list_element) {
  # Extract the first (and only) dataframe from the list element
  data_subset <- list_element

  # Call the analysis function
  fit_and_save_segreg(data_subset, name)
}

cl1 <- makeCluster(detectCores() - 1)
clusterExport(cl1, varlist = c("fit_and_save_segreg", "data_young_tagged_segreg",
                               "process_function_segreg"), envir = environment())
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
parLapply(cl1, names(data_young_tagged_segreg), function(name) {
  process_function_segreg(name, data_young_tagged_segreg[[name]])
})
stopCluster(cl1)



cl2 <- makeCluster(detectCores() - 1)

clusterExport(cl2, varlist = c("fit_and_save_segreg", "data_old_tagged_segreg",
                               "process_function_segreg"), envir = environment())
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
parLapply(cl2, names(data_old_tagged_segreg), function(name) {
  process_function_segreg(name, data_old_tagged_segreg[[name]])
})
# Clean up
stopCluster(cl2)
