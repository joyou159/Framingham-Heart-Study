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

# Scaling

continuous_vars <- c("age", "cigs_per_day", "tot_chol",
                     "pulse_pressure", "bmi", "heart_rate", "glucose")

log_vars        <- c("glucose", "bmi", "pulse_pressure", "tot_chol", "heart_rate")

cat_vars        <- c("sex", "education", "current_smoker", "bp_meds",
                     "prevalent_stroke", "prevalent_hyp", "diabetes", "ten_year_chd")

train <- train |>
  mutate(across(all_of(cat_vars), ~ as.integer(as.factor(.x)) - 1L))

test <- test |>
  mutate(across(all_of(cat_vars), ~ as.integer(as.factor(.x)) - 1L))

# Verify: all should show 0 and 1
train |> select(all_of(cat_vars)) |> 
  summarise(across(everything(), list(min=min, max=max)))

train_means <- train |>
  mutate(across(all_of(log_vars), log)) |>
  summarise(across(all_of(continuous_vars), ~mean(.x, na.rm = TRUE)))

train_sds <- train |>
  mutate(across(all_of(log_vars), log)) |>
  summarise(across(all_of(continuous_vars), ~sd(.x, na.rm = TRUE)))

# Scale training set — log transform then z-score
train_scaled <- train |>
  mutate(across(all_of(log_vars), log)) |>
  mutate(across(all_of(continuous_vars),
                ~(.x - train_means[[cur_column()]]) / train_sds[[cur_column()]]))

# Scale test set using train parameters
test_scaled <- test |>
  mutate(across(all_of(log_vars), log)) |>
  mutate(across(all_of(continuous_vars),
                ~(.x - train_means[[cur_column()]]) / train_sds[[cur_column()]]))

# Verification — training means should be ~0, sds ~1
train_scaled |>
  summarise(across(all_of(continuous_vars),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE))))

# Missingness summary
framingham_processed |>
  summarise(across(everything(), ~sum(is.na(.))))

train_scaled |>
  select(all_of(continuous_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "steelblue", colour = "white", linewidth = 0.2) +
  geom_density(colour = "firebrick", linewidth = 0.7) +
  facet_wrap(~variable, scales = "free", ncol = 4) +
  labs(title = "Continuous predictors (log-transformed & scaled)",
       x = NULL, y = "Density") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave("Figures/cont_distributions_log_and_scaled.png")


# tot_chol imputation model ----------------------------------------

predictors_train <- train_scaled |>
  select(-tot_chol) |>
  as.matrix()


# Complete cases only for model fitting
complete_rows  <- complete.cases(predictors_train) & !is.na(train_scaled$tot_chol)

y_fit          <- train_scaled$tot_chol[complete_rows]   # log-transformed + scaled
predictors_fit <- predictors_train[complete_rows,]

predictors_fit

length(y_fit)
dim(predictors_fit)
colnames(predictors_fit)  

train_miss_idx        <- which(is.na(train_scaled$tot_chol))
predictors_train_miss <- predictors_train[train_miss_idx, ]

train_miss_idx

# Predictor matrix — test set
predictors_test <- test_scaled |>
  select(-tot_chol) |>
  as.matrix()

# Test rows with missing tot_chol
test_miss_idx        <- which(is.na(test_scaled$tot_chol))
predictors_test_miss <- predictors_test[test_miss_idx, ]

test_miss_idx

# JAGS model

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
  n.iter = 20000,
  thin   = 5
  )

# Diagnostics

cat("\n--- Posterior Summary ---\n")
print(summary(posterior_totchol))

cat("\nGelman-Rubin Diagnostic:\n")
print(gelman.diag(posterior_totchol))

cat("\nEffective Sample Sizes:\n")
print(effectiveSize(posterior_totchol))


posterior_matrix <- as.matrix(posterior_totchol)

params_to_plot <- c("beta[2]", "beta[6]", "beta[1]",
                    "beta[10]", "beta[11]", "beta[3]")

labels <- c("Age", "BP Meds", "Sex",
            "BMI", "Heart Rate", "Education")

posterior_df <- as.data.frame(posterior_matrix[, params_to_plot])
colnames(posterior_df) <- labels

posterior_df |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = labels)) |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.75, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Posterior Distributions — Significant Predictors of Total Cholesterol",
    subtitle = "95% credible interval excludes zero | Dashed line at zero for reference",
    x        = "Posterior value",
    y        = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/posterior_distributions_log_tot_chol.png")


# imputing missing data in train and test split

beta0_samples <- posterior_matrix[, "beta0"]
beta_samples  <- posterior_matrix[, paste0("beta[", 1:14, "]")]
sigma_samples <- posterior_matrix[, "sigma"]
n_samples     <- length(beta0_samples)
cat("Number of posterior samples:", n_samples, "\n")

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

ppc_train_orig <- exp(ppc_train * train_sds[["tot_chol"]] + train_means[["tot_chol"]])
ppc_test_orig  <- exp(ppc_test  * train_sds[["tot_chol"]] + train_means[["tot_chol"]])



# posterior summary
imputed_train_totchol <- round(colMeans(ppc_train_orig))
imputed_test_totchol  <- round(colMeans(ppc_test_orig))

ci_train <- apply(ppc_train_orig, 2, quantile, probs = c(0.025, 0.975))
ci_test  <- apply(ppc_test_orig,  2, quantile, probs = c(0.025, 0.975))

cat("Imputed train tot_chol values (original scale):\n")
print(imputed_train_totchol)
cat("\n95% Credible Intervals — train (original scale):\n")
print(round(ci_train, 1))

cat("\nImputed test tot_chol values (original scale):\n")
print(imputed_test_totchol)
cat("\n95% Credible Intervals — test (original scale):\n")
print(round(ci_test, 1))

# posterior predictive distribution for some missing samples

plot_data <- bind_rows(
  as.data.frame(ppc_train_orig[, 1:3]) |>
    setNames(paste0("Train Obs ", train_miss_idx[1:3])) |>
    pivot_longer(everything(), names_to = "observation", values_to = "value") |>
    mutate(set = "Train"),
  
  as.data.frame(ppc_test_orig[, 1:3]) |>
    setNames(paste0("Test Obs ", test_miss_idx[1:3])) |>
    pivot_longer(everything(), names_to = "observation", values_to = "value") |>
    mutate(set = "Test")
)

mean_labels <- plot_data |>
  group_by(observation, set) |>
  summarise(mean_val = mean(value), .groups = "drop") |>
  mutate(label = paste0("Mean = ", round(mean_val)))

plot_data |>
  ggplot(aes(x = value, fill = set)) +
  geom_density(alpha = 0.7, color = "white", linewidth = 0.3) +
  geom_vline(data = mean_labels,
             aes(xintercept = mean_val),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_text(data = mean_labels,
            aes(x = mean_val, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5,
            size = 2.8, color = "grey20", inherit.aes = FALSE) +
  facet_wrap(~observation, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  labs(
    title    = "Posterior Predictive Distributions — Missing Total Cholesterol",
    subtitle = "Sample of 3 train and 3 test observations | Dashed line = posterior mean",
    x        = "Total Cholesterol (mg/dL)",
    y        = "Density",
    fill     = "Dataset"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )
ggsave("Figures/posterior_predictive_distribution_missing_tot_chol.png")


train$tot_chol[train_miss_idx] <- imputed_train_totchol
test$tot_chol[test_miss_idx]   <- imputed_test_totchol


cat("Missing tot_chol in train:", sum(is.na(train$tot_chol)), "\n")
cat("Missing tot_chol in test:",  sum(is.na(test$tot_chol)),  "\n")

saveRDS(train, "Data/train.rds")
saveRDS(test,  "Data/test.rds")

# modeling missingness in BMI -------------------------

train <- readRDS("Data/train.rds")
test  <- readRDS("Data/test.rds")

cat("Missing tot_chol in train:", sum(is.na(train$tot_chol)), "\n")
cat("Missing tot_chol in test:",  sum(is.na(test$tot_chol)),  "\n")

# Remaining missingness
train |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0)


# data preprocessing 

train_means_bmi <- train |>
  mutate(across(all_of(log_vars), log)) |>
  summarise(across(all_of(continuous_vars), ~mean(.x, na.rm = TRUE)))

train_sds_bmi <- train |>
  mutate(across(all_of(log_vars), log)) |>
  summarise(across(all_of(continuous_vars), ~sd(.x, na.rm = TRUE)))

train_scaled_bmi <- train |>
  mutate(across(all_of(log_vars), log)) |>
  mutate(across(all_of(continuous_vars),
                ~(.x - train_means_bmi[[cur_column()]]) / train_sds_bmi[[cur_column()]]))

test_scaled_bmi <- test |>
  mutate(across(all_of(log_vars), log)) |>
  mutate(across(all_of(continuous_vars),
                ~(.x - train_means_bmi[[cur_column()]]) / train_sds_bmi[[cur_column()]]))

predictors_train_bmi <- train_scaled_bmi |>
  select(-bmi) |>
  as.matrix()

predictors_test_bmi <- test_scaled_bmi |>
  select(-bmi) |>
  as.matrix()

# Complete cases for fitting
complete_rows_bmi  <- complete.cases(predictors_train_bmi) & !is.na(train_scaled_bmi$bmi)
y_fit_bmi          <- train_scaled_bmi$bmi[complete_rows_bmi]
predictors_fit_bmi <- predictors_train_bmi[complete_rows_bmi, ]

colnames(predictors_fit_bmi)  


cat("y_fit_bmi length:", length(y_fit_bmi), "\n")
cat("predictors_fit_bmi dim:", dim(predictors_fit_bmi), "\n")

# Missing indices
train_miss_idx_bmi        <- which(is.na(train_scaled_bmi$bmi))
predictors_train_miss_bmi <- predictors_train_bmi[train_miss_idx_bmi, ]

test_miss_idx_bmi        <- which(is.na(test_scaled_bmi$bmi))



cat("Train missing BMI:", length(train_miss_idx_bmi), "\n")
cat("Test missing BMI:",  length(test_miss_idx_bmi),  "\n") # = 0 (no imputation needed)


# JAGS model — BMI imputation
bmi_model_string <- "
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

jags_data_bmi <- list(
  y = y_fit_bmi,
  X = predictors_fit_bmi,
  N = length(y_fit_bmi),
  K = ncol(predictors_fit_bmi)
)

bmi_model <- jags.model(
  textConnection(bmi_model_string),
  data     = jags_data_bmi,
  n.chains = 3,
  n.adapt  = 3000
)

update(bmi_model, n.iter = 10000)

posterior_bmi <- coda.samples(
  bmi_model,
  variable.names = c("beta0", "beta", "sigma"),
  n.iter = 20000,
  thin   = 5
)

# Diagnostics
cat("\n--- Posterior Summary ---\n")
print(summary(posterior_bmi))

cat("\nGelman-Rubin Diagnostic:\n")
print(gelman.diag(posterior_bmi))

cat("\nEffective Sample Sizes:\n")
print(effectiveSize(posterior_bmi))



# plotting

posterior_matrix_bmi <- as.matrix(posterior_bmi)

params_to_plot_bmi <- c("beta[8]", "beta[4]", "beta[1]",
                        "beta[10]", "beta[3]", "beta[9]")

labels_bmi <- c("Prev. Hypertension", "Current Smoker", "Sex",
                "Total Cholesterol",  "Education",      "Diabetes")

posterior_df_bmi <- as.data.frame(posterior_matrix_bmi[, params_to_plot_bmi])
colnames(posterior_df_bmi) <- labels_bmi

posterior_df_bmi |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = labels_bmi)) |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.75, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Posterior Distributions — Significant Predictors of BMI",
    subtitle = "95% credible interval excludes zero | Dashed line at zero for reference",
    x        = "Posterior value",
    y        = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/posterior_distributions_log_bmi.png")


# data imputation
beta0_samples_bmi <- posterior_matrix_bmi[, "beta0"]
beta_samples_bmi  <- posterior_matrix_bmi[, paste0("beta[", 1:14, "]")]
sigma_samples_bmi <- posterior_matrix_bmi[, "sigma"]
n_samples_bmi     <- length(beta0_samples_bmi)
cat("Number of posterior samples:", n_samples_bmi, "\n")

n_train_miss_bmi <- nrow(predictors_train_miss_bmi)
ppc_train_bmi    <- matrix(NA, nrow = n_samples_bmi, ncol = n_train_miss_bmi)

for (s in 1:n_samples_bmi) {
  mu_train_bmi       <- beta0_samples_bmi[s] + predictors_train_miss_bmi %*% beta_samples_bmi[s, ]
  ppc_train_bmi[s, ] <- rnorm(n_train_miss_bmi, mean = mu_train_bmi, sd = sigma_samples_bmi[s])
}

ppc_train_bmi_orig <- exp(ppc_train_bmi * train_sds_bmi[["bmi"]] + train_means_bmi[["bmi"]])

imputed_train_bmi <- round(colMeans(ppc_train_bmi_orig), 1)
ci_train_bmi      <- apply(ppc_train_bmi_orig, 2, quantile, probs = c(0.025, 0.975))

cat("Imputed train BMI values (original scale):\n")
print(imputed_train_bmi)
cat("\n95% Credible Intervals — train (original scale):\n")
print(round(ci_train_bmi, 1))

cat("\nImputed BMI range:", range(imputed_train_bmi), "\n")
cat("Observed BMI range:", range(train$bmi, na.rm = TRUE), "\n")



set.seed(123)
sample_idx <- sample(1:n_train_miss_bmi, 6)

plot_data_bmi <- as.data.frame(ppc_train_bmi_orig[, sample_idx]) |>
  setNames(paste0("Train Obs ", train_miss_idx_bmi[sample_idx])) |>
  pivot_longer(everything(), names_to = "observation", values_to = "value")

mean_labels_bmi <- plot_data_bmi |>
  group_by(observation) |>
  summarise(mean_val = mean(value), .groups = "drop") |>
  mutate(label = paste0("Mean = ", round(mean_val, 1)))

plot_data_bmi |>
  ggplot(aes(x = value, fill = observation)) +
  geom_density(alpha = 0.75, color = "white", linewidth = 0.3) +
  geom_vline(data = mean_labels_bmi,
             aes(xintercept = mean_val),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_text(data = mean_labels_bmi,
            aes(x = mean_val, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5,
            size = 2.8, color = "grey20", inherit.aes = FALSE) +
  facet_wrap(~observation, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Posterior Predictive Distributions — Missing BMI (Train)",
    subtitle = "Random sample of 6 observations | Dashed line = posterior mean",
    x        = "BMI (kg/m²)",
    y        = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/posterior_predictive_distribution_missing_bmi.png")

# saving imputed values
train$bmi[train_miss_idx_bmi] <- imputed_train_bmi

cat("Missing BMI in train:", sum(is.na(train$bmi)), "\n")
cat("Missing BMI in test:",  sum(is.na(test$bmi)),  "\n")

saveRDS(train, "Data/train.rds")
saveRDS(test,  "Data/test.rds")

# modeling missingness in cigs/day -------------------------

train <- readRDS("Data/train.rds")
test  <- readRDS("Data/test.rds")


cat("\nRemaining missingness in train:\n")
train |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |> print()




