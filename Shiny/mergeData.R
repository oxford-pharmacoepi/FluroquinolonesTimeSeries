library(here)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(png)
library(tidyverse)

extract_info_and_combine_files <- function(directory) {
  # Define the mapping from filename keywords to age groups
  age_group_mapping <- list(
    '60 - 150' = '60 and over',
    '19 - 59' = '19 - 59'
    # Add other mappings as necessary
  )
  
  # Get all CSV files in the directory
  files <- list.files(path = directory, pattern = "\\.csv$", full.names = TRUE)
  
  # Iterate over files, read data, extract info, and bind rows
  combined_df <- files %>%
    map_df(function(file_path) {
      # Extract model name and age keyword from the filename
      parts <- str_split(basename(file_path), "_")[[1]]  # Split the filename into parts
      model_name <- parts[3]  # Assuming 'modelX' is always in the third position
      age_keyword <- parts[length(parts)]  # Assuming age-related keyword is always last
      age_keyword <- str_replace(age_keyword, "\\.csv", "")  # Remove the .csv extension
      
      # Map the age keyword to an age group
      age_group <- age_group_mapping[[age_keyword]] %>% 
        if_else(is.null(.), 'Unknown', .)  # Default to 'Unknown' if not found
      
      # Read the CSV file
      df <- read_csv(file_path)
      df <- df %>%
        rename_with(~ gsub("Std..Error", "Std. Error", .)) %>%
        rename_with(~ gsub("P.value", "P-value", .)) %>%
        rename_with(~ gsub("deviance", "Deviance", .))
      
      # Add the extracted info as new columns
      df <- df %>%
        mutate(cdm_name = model_name, age_group = age_group)
      
      return(df)
    })
  
  return(combined_df)
}

arima_para <- extract_info_and_combine_files("data_gold/arima_para")
arima_coef <- extract_info_and_combine_files("data_gold/arima_coef")

segreg_coef <- extract_info_and_combine_files("data_gold/segreg_coef")



# --------- fitted values -----------

extract_info_and_combine_arima_pred <- function(directory) {
  # Define the mapping from filename keywords to age groups
  age_group_mapping <- list(
    '60 - 150' = '60 and over',
    '19 - 59' = '19 - 59'
    # Add other mappings as necessary
  )
  
  # Get all CSV files in the directory
  files <- list.files(path = directory, pattern = "\\.csv$", full.names = TRUE)
  
  # Iterate over files, read data, extract info, and bind rows
  combined_df <- files %>%
    map_df(function(file_path) {
      # Extract model name and age keyword from the filename
      parts <- str_split(basename(file_path), "_")[[1]]  # Split the filename into parts
      model_name <- parts[3]  # Assuming 'modelX' is always in the third position
      age_keyword <- parts[length(parts)]  # Assuming age-related keyword is always last
      age_keyword <- str_replace(age_keyword, "\\.csv", "")  # Remove the .csv extension
      
      # Map the age keyword to an age group
      age_group <- age_group_mapping[[age_keyword]] %>% if_else(is.null(.), 'Unknown', .)  # Default to 'Unknown' if not found
      
      # Read the CSV file
      df <- read_csv(file_path)
      # Add the extracted info as new columns
      df <- df %>%
        mutate(cdm_name = model_name, age_group = age_group)
      
      return(df)
    })
  
  return(combined_df)
}
arima_pred <- extract_info_and_combine_arima_pred("data_gold/arima_pred")

files <- list.files(path = "data_gold/segreg_pred", full.names = TRUE)
segreg_pred <- lapply(files, read.csv)%>%
  bind_rows() %>% mutate(age_group = if_else(age_group == "60 - 150", "60 and over",
                                             "19 - 59"))

incidenceTable_4 <- read_csv("data_gold/incidenceTable (4).csv")
incidenceTable_4 <- incidenceTable_4 %>% mutate(age_group = case_when(denominator_age_group == "19 to 59" ~ "19 - 59",
                                                                      denominator_age_group == "60 to 150" ~ "60 and over",
                                                                      TRUE ~ .data$denominator_age_group)) %>%
  filter(age_group %in% c("60 and over", "19 - 59"),
         denominator_sex == "Both") %>%
  select(age_group,
         cdm_name,
         incidence_100000_pys,
         incidence_100000_pys_95CI_lower,
         incidence_100000_pys_95CI_upper,
         incidence_start_date) %>%
  mutate(cdm_name = if_else(cdm_name=="project_3619", "East Scotland data", cdm_name))%>%
  mutate(cdm_name = if_else(cdm_name=="CPRDAurumFull", "CPRD Aurum", cdm_name))

arima_pred <- arima_pred %>% mutate(Date = as.Date(Date))  %>%
  mutate(cdm_name = if_else(cdm_name=="project", "East Scotland data", cdm_name)) %>% 
  mutate(cdm_name = if_else(cdm_name=="CPRDAurumFull", "CPRD Aurum", cdm_name)) %>%
  left_join(incidenceTable_4, by = c("cdm_name",
                                     "age_group",
                                     "Date" = "incidence_start_date")) %>% 
  filter(Date <= as.Date("2022-01-01"))

segreg_pred <- segreg_pred %>% mutate(Date = as.Date(Date)) %>% 
  mutate(cdm_name = if_else(cdm_name=="project_3619", "East Scotland data", cdm_name)) %>%
  mutate(cdm_name = if_else(cdm_name=="CPRDAurumFull", "CPRD Aurum", cdm_name)) %>%
  left_join(incidenceTable_4, by =
              c("cdm_name",
                "age_group",
                "Date" = "incidence_start_date"))%>%
  filter(Date <= as.Date("2022-01-01"))

save(
  segreg_pred, arima_pred,
  segreg_coef, arima_coef,
  arima_para,
  file = here("mergeData.Rdata")
)
