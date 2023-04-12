
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
library(tidymodels)
library(mclust)


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
            by = c('USUBJID')
            ) %>%
  pivot_wider(names_from = VSTESTCD,
              values_from = VSSTRESN) %>%
  bind_rows(., slice(., rep(1, 1)))

# Model ------------------------------------------------------------------------
# Split the patient data into training and testing sets
patient_split <- initial_split(patient_data, prop = 0.8, strata = "SEX")
patient_train <- training(patient_split)
patient_test <- testing(patient_split)

# Preprocess the patient data using a recipe
patient_recipe <- recipe(USUBJID ~ ., data = patient_train) %>%
  step_rm(USUBJID) %>%
  step_dummy(SEX, RACE, ETHNIC) %>%
  step_center(all_numeric(), -all_outcomes()) %>%
  step_scale(all_numeric(), -all_outcomes())

# Fit a k-means clustering model to the preprocessed patient data
patient_kmeans <- patient_recipe %>%
  prep() %>%
  bake(patient_train) %>%
  mclust::Mclust(method = "EM") %>%
  predict()

# Add the cluster labels to the patient data
patient_train <- patient_train %>%
  mutate(cluster = patient_kmeans$classification)

# Identify potential duplicate patients within each cluster
potential_duplicates <- patient_train %>%
  group_by(cluster) %>%
  filter(n() > 1) %>%
  group_by(YOB, SEX, RACE, ETHNIC) %>% #, HEIGHT, WEIGHT, TEMP, SYSBP, DIABP) %>%
  summarize(n_duplicates = n(),
            .groups = 'keep') %>%
  filter(n_duplicates > 1)

# Print the potential duplicate patients
potential_duplicates

# Plot the data colored by cluster assignment
P <- ggplot(potential_duplicates,
            aes(x = X1, y = X2, color = cluster)) +
  geom_point() +
  scale_color_discrete(name = "Cluster") +
  theme_BW()

P