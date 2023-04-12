
# Purpose -----------------------------------------------------------------
# Identify potentially duplicate enrolled subjects based on:
# age
# dob
# sex
# ethnicity
# race
# height
# weight
# temperature
# systolic_bp
# diastolic_bp

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(haven)

# Data --------------------------------------------------------------------
data <- purrr::map(c('vs', 'dm'),
                   ~haven::read_xpt(paste0('data/', .x, '.xpt'))
) %>%
  `names<-`(c('vs', 'dm'))

dm <- data$dm %>%
  select(USUBJID, SEX, RACE, ETHNIC, ARM) %>%
  filter(ARM != 'Screen Failure') %>%
  rowwise() %>%
  mutate(YOB = sample(1950:2004, 1)) %>%
  ungroup()

vs <- data$vs %>%
  filter(EPOCH == 'SCREENING' & !is.na(VSORRES)) %>%
  group_by(USUBJID, EPOCH, VSTESTCD) %>%
  mutate(VSSEQ = row_number()) %>%
  ungroup() %>%
  filter(VSSEQ == 1) %>%
  select(USUBJID, VSTESTCD, VSSTRESN)

patient_data <- left_join(dm %>% select(-ARM),
                          vs,
                          by = c('USUBJID')) %>%
  pivot_wider(names_from = VSTESTCD,
              values_from = VSSTRESN) %>%
  bind_rows(., slice(., rep(1, 1)))

# Analysis ---------------------------------------------------------------------
# Define the variables and their expected variances
variances <- list(
  SYSBP = 2,
  DIABP = 3,
  PULSE = 2,
  TEMP = 1,
  HEIGHT = 1,
  WEIGHT = 3
)

# Define the quantiles for each variable
quantiles <- list(
  SYSBP = 0.90,
  DIABP = 0.95,
  PULSE = 0.95,
  TEMP = 0.95,
  HEIGHT = 0.99,
  WEIGHT = 0.90
)

# Define a function to calculate the Mahalanobis distance for a single variable
mahalanobis_var <- function(x, var, q) {
  dist <- mahalanobis(x = x, center = mean(x), cov = var(x), inverted = TRUE)
  threshold <- qchisq(p = q, df = 1)  # Assume 1 degree of freedom
  duplicates <- which(dist < threshold)
  return(duplicates)
}

# Create a list of the input vectors
input_list <- list(
  data = patient_data %>% select(names(variances)),
  variances = variances,
  quantiles = quantiles
)

# Use purrr to apply the function to each variable and combine the results
duplicates <- purrr::pmap(input_list, mahalanobis_var) %>%
  reduce(union)


# Calculate Mahalanobis distance for each variable separately
dists_systolic_bp <- mahalanobis(x = data$systolic_bp,
                                 center = mean(data$systolic_bp),
                                 cov = var(data$systolic_bp),
                                 inverted = TRUE)

dists_diastolic_bp <-
  mahalanobis(
    x = data$diastolic_bp,
    center = mean(data$diastolic_bp),
    cov = var(data$diastolic_bp),
    inverted = TRUE
  )

dists_total_cholesterol <-
  mahalanobis(
    x = data$total_cholesterol,
    center = mean(data$total_cholesterol),
    cov = var(data$total_cholesterol),
    inverted = TRUE
  )

# Define different thresholds for each variable
threshold_age <- qchisq(p = 0.95, df = 1)  # Age has 1 degree of freedom
threshold_systolic_bp <- qchisq(p = 0.90, df = 1)  # Systolic blood pressure has higher variance, so we use a lower threshold
threshold_diastolic_bp <- qchisq(p = 0.95, df = 1)  # Diastolic blood pressure has lower variance, so we use a higher threshold
threshold_total_cholesterol <- qchisq(p = 0.95, df = 1)  # Total cholesterol has similar variance to age, so we use a similar threshold

# Identify potential duplicates for each variable based on its own threshold
duplicates_age <- which(dists_age < threshold_age)
duplicates_systolic_bp <- which(dists_systolic_bp < threshold_systolic_bp)
duplicates_diastolic_bp <- which(dists_diastolic_bp < threshold_diastolic_bp)
duplicates_total_cholesterol <- which(dists_total_cholesterol < threshold_total_cholesterol)

# Combine duplicate patient IDs from all variables
duplicates <- unique(c(duplicates_age, duplicates_systolic_bp, duplicates_diastolic_bp, duplicates_total_cholesterol))
