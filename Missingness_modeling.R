library(tidyverse)
library(rsample)


framingham_processed <- readRDS("Data/framingham_processed.rds")


glimpse(framingham_processed)
dim(framingham_processed)


# Train/test split stratified by ten_year_chd

set.seed(123)

split <- initial_split(framingham_processed, 
                       prop = 0.8, 
                       strata = ten_year_chd)

train <- training(split)
test  <- testing(split)

dim(train)
dim(test)

train |> count(ten_year_chd) |> mutate(prop = n / sum(n))
test  |> count(ten_year_chd) |> mutate(prop = n / sum(n))



# continuous variables to scale
continuous_vars <- c("age", "cigs_per_day", "tot_chol",
                     "pulse_pressure", "bmi", "heart_rate", "glucose")

train_means <- train |> 
  summarise(across(all_of(continuous_vars), ~mean(.x, na.rm = TRUE)))

train_sds <- train |> 
  summarise(across(all_of(continuous_vars), ~sd(.x, na.rm = TRUE)))

# Scale training set
train_scaled <- train |>
  mutate(across(all_of(continuous_vars),
                ~(.x - train_means[[cur_column()]]) / train_sds[[cur_column()]]))

# Scale test set using train parameters
test_scaled <- test |>
  mutate(across(all_of(continuous_vars),
                ~(.x - train_means[[cur_column()]]) / train_sds[[cur_column()]]))


# verification
train_scaled |>
  summarise(across(all_of(continuous_vars),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE))))


# checking missingness after dropping co-missingness
train_scaled |>
  summarise(across(everything(), ~sum(is.na(.))))









