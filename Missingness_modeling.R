library(tidyverse)
library(rsample)
library(rjags)
library(coda)



# Data loading and splitting -------------------------------------

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


saveRDS(train, "Data/train.rds")
saveRDS(test,  "Data/test.rds")



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
framingham_processed |>
  summarise(across(everything(), ~sum(is.na(.))))




# modeling missiness in total_chol ----------------------------------


train |>
  ggplot(aes(x = tot_chol)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 40, 
                 fill = "steelblue", 
                 color = "white",
                 alpha = 0.7) +
  geom_density(color = "darkred", linewidth = 0.8) +
  geom_vline(aes(xintercept = mean(tot_chol, na.rm = TRUE)),
             color = "darkred", linetype = "dashed") +
  labs(title = "Distribution of Total Cholesterol",
       subtitle = "Training set — original scale",
       x = "Total Cholesterol (mg/dL)",
       y = "Density") +
  theme_minimal(base_size = 11)

ggsave("Figures/total_chol_distribution_for_imputation.png")


# Complete cases for model fitting (train only)

predictors_train <- train_scaled |>
  mutate(
    sex              = as.integer(sex),
    education        = as.integer(education),
    current_smoker   = as.integer(current_smoker),
    bp_meds          = as.integer(bp_meds),
    prevalent_stroke = as.integer(prevalent_stroke),
    prevalent_hyp    = as.integer(prevalent_hyp),
    diabetes         = as.integer(diabetes),
    ten_year_chd     = as.integer(ten_year_chd)
  ) |>
  select(-tot_chol) |> # remove that column
  as.matrix()


complete_rows  <- complete.cases(predictors_train) & !is.na(train$tot_chol)


y_fit          <- train$tot_chol[complete_rows]
predictors_fit <- predictors_train[complete_rows, ]

length(y_fit) # 2912
dim(predictors_fit) # 2912 x 14


# Train rows with missing tot_chol

train_miss_idx        <- which(is.na(train$tot_chol))
predictors_train_miss <- predictors_train[train_miss_idx, ]


train_miss_idx

# Test rows with missing tot_chol

predictors_test <- test_scaled |>
  mutate(
    sex              = as.integer(sex),
    education        = as.integer(education),
    current_smoker   = as.integer(current_smoker),
    bp_meds          = as.integer(bp_meds),
    prevalent_stroke = as.integer(prevalent_stroke),
    prevalent_hyp    = as.integer(prevalent_hyp),
    diabetes         = as.integer(diabetes),
    ten_year_chd     = as.integer(ten_year_chd)
  ) |>
  select(-tot_chol) |>
  as.matrix()

test_miss_idx        <- which(is.na(test$tot_chol))
predictors_test_miss <- predictors_test[test_miss_idx, ]


test_miss_idx



predictors_test_miss

# model definition

chol_model_string <- "
model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + inprod(X[i,], beta[])
  }

  # Priors
  beta0 ~ dnorm(0, 0.01)
  
  for (j in 1:K) {
    beta[j] ~ dnorm(0, 0.01)
  }
  
  tau   ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}
"

jags_data <- list(
  y = y_fit,
  X = predictors_fit,
  N = length(y_fit),
  K = ncol(predictors_fit)
)


chol_model <- jags.model(
  textConnection(chol_model_string),
  data     = jags_data,
  n.chains = 3,
  n.adapt  = 3000
)


update(chol_model, n.iter = 10000)

posterior_totchol <- coda.samples(
  chol_model,
  variable.names = c("beta0", "beta", "sigma"),
  n.iter = 50000,
  thin   = 20
)


cat("\n--- Posterior Summary ---\n")
print(summary(posterior_totchol))

cat("\nGelman-Rubin Diagnostic:\n")
print(gelman.diag(posterior_totchol))

cat("\nEffective Sample Sizes:\n")
print(effectiveSize(posterior_totchol))


posterior_matrix <- as.matrix(posterior_totchol)

params_to_plot <- c("beta0", "beta[2]", "beta[6]", "beta[7]", 
                    "beta[9]", "beta[13]",
                    "beta[10]", "sigma")

param_labels <- c("beta0 — Intercept", "beta[2] — Age", "beta[6] — BP Meds", 
                  "beta[7] — Prev. Stroke", "beta[9] — Diabetes", 
                  "beta[13] — CHD",
                  "beta[10] — BMI", "sigma")

posterior_df <- as.data.frame(posterior_matrix[, params_to_plot])
colnames(posterior_df) <- param_labels

posterior_df |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = param_labels)) |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.7, color = "white") +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Posterior Distributions — Total Cholesterol Imputation Model",
       x = "Value", y = "Density") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

ggsave("Figures/Posterior_distributions_tot_chol.png")



# imputing missing data in train and test split

beta0_samples <- posterior_matrix[, "beta0"]
beta_samples  <- posterior_matrix[, paste0("beta[", 1:14, "]")]
sigma_samples <- posterior_matrix[, "sigma"]

n_samples <- length(beta0_samples)

n_samples


# Posterior predictive for train missing
n_train_miss <- nrow(predictors_train_miss)
ppc_train    <- matrix(NA, nrow = n_samples, ncol = n_train_miss)

for (s in 1:n_samples) {
  mu_train       <- beta0_samples[s] + predictors_train_miss %*% beta_samples[s, ]
  ppc_train[s, ] <- rnorm(n_train_miss, mean = mu_train, sd = sigma_samples[s])
}


# Posterior predictive for test missing
n_test_miss <- nrow(predictors_test_miss)
ppc_test    <- matrix(NA, nrow = n_samples, ncol = n_test_miss)

for (s in 1:n_samples) {
  mu_test       <- beta0_samples[s] + predictors_test_miss %*% beta_samples[s, ]
  ppc_test[s, ] <- rnorm(n_test_miss, mean = mu_test, sd = sigma_samples[s])
}


# posterior summary
imputed_train_totchol <- round(colMeans(ppc_train))
imputed_test_totchol  <- round(colMeans(ppc_test))

ci_train <- apply(ppc_train, 2, quantile, probs = c(0.025, 0.975))
ci_test  <- apply(ppc_test,  2, quantile, probs = c(0.025, 0.975))

cat("Imputed train tot_chol values:\n")
print(imputed_train_totchol)
cat("\n95% Credible Intervals (train):\n")
print(round(ci_train, 1))

cat("\nImputed test tot_chol values:\n")
print(imputed_test_totchol)
cat("\n95% Credible Intervals (test):\n")
print(round(ci_test, 1))


# posterior predictive distribution for some missing samples

plot_data <- bind_rows(
  as.data.frame(ppc_train[, 1:3]) |>
    setNames(paste0("Train Obs ", train_miss_idx[1:3])) |>
    pivot_longer(everything(), names_to = "observation", values_to = "value") |>
    mutate(set = "Train"),
  
  as.data.frame(ppc_test[, 1:3]) |>
    setNames(paste0("Test Obs ", test_miss_idx[1:3])) |>
    pivot_longer(everything(), names_to = "observation", values_to = "value") |>
    mutate(set = "Test")
)

mean_labels <- plot_data |>
  group_by(observation) |>
  summarise(mean_val = mean(value)) |>
  mutate(label = paste0("Mean = ", round(mean_val)))

plot_data |>
  ggplot(aes(x = value, fill = set)) +
  geom_density(alpha = 0.7, color = "white") +
  facet_wrap(~observation, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  geom_vline(data = mean_labels,
             aes(xintercept = mean_val),
             linetype = "dashed", color = "black") +
  geom_text(data = mean_labels,
            aes(x = mean_val, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5,
            size = 2.8, color = "black", inherit.aes = FALSE) +
  labs(title = "Posterior Predictive Distributions (Sample of missing Total Cholesterol obs)",
       subtitle = "Dashed line = posterior mean",
       x = "Total Cholesterol (mg/dL)",
       y = "Density",
       fill = "Dataset") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("Figures/posterior_predictive_distribution_missing_tot_chol.png")


train$tot_chol[train_miss_idx] <- imputed_train_totchol
test$tot_chol[test_miss_idx]   <- imputed_test_totchol


cat("Missing tot_chol in train:", sum(is.na(train$tot_chol)), "\n")
cat("Missing tot_chol in test:",  sum(is.na(test$tot_chol)),  "\n")

saveRDS(train, "Data/train.rds")
saveRDS(test,  "Data/test.rds")



