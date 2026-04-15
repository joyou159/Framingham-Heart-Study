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

# How many smokers report 0?
train |>
  filter(current_smoker == 1) |>
  count(cigs_per_day == 0)

train |>
  filter(current_smoker == 1) |>
  ggplot(aes(x = cigs_per_day)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Distribution of cigs/day among smokers",
       x = "Cigarettes per day", y = "Count") +
  theme_minimal()

# Non-smokers: assign 0 directly
train$cigs_per_day[train$current_smoker == 0 & is.na(train$cigs_per_day)] <- 0
test$cigs_per_day[test$current_smoker  == 0 & is.na(test$cigs_per_day)]  <- 0

cat("Missing cigs_per_day after zero-fill:\n")
cat("Train:", sum(is.na(train$cigs_per_day)), "\n")
cat("Test: ", sum(is.na(test$cigs_per_day)),  "\n")

train |>
  filter(is.na(cigs_per_day)) |>
  count(current_smoker)

test |>
  filter(is.na(cigs_per_day)) |>
  count(current_smoker)



# Subset to smokers only
train_smokers <- train |> filter(current_smoker == 1)
test_smokers  <- test  |> filter(current_smoker == 1)


# Scale continuous predictors
continuous_vars_cigs <- c("age", "tot_chol", "pulse_pressure", 
                          "bmi", "heart_rate", "glucose")

log_vars_cigs <- c("glucose", "bmi", "pulse_pressure", "tot_chol", "heart_rate")

train_means_cigs <- train_smokers |>
  mutate(across(all_of(log_vars_cigs), log)) |>
  summarise(across(all_of(continuous_vars_cigs), ~mean(.x, na.rm = TRUE)))

train_sds_cigs <- train_smokers |>
  mutate(across(all_of(log_vars_cigs), log)) |>
  summarise(across(all_of(continuous_vars_cigs), ~sd(.x, na.rm = TRUE)))

train_smokers_scaled <- train_smokers |>
  mutate(across(all_of(log_vars_cigs), log)) |>
  mutate(across(all_of(continuous_vars_cigs),
                ~(.x - train_means_cigs[[cur_column()]]) / train_sds_cigs[[cur_column()]]))

test_smokers_scaled <- test_smokers |>
  mutate(across(all_of(log_vars_cigs), log)) |>
  mutate(across(all_of(continuous_vars_cigs),
                ~(.x - train_means_cigs[[cur_column()]]) / train_sds_cigs[[cur_column()]]))


# JAGS data preparation
predictors_train_cigs <- train_smokers_scaled |>
  select(-cigs_per_day, -current_smoker) |>  
  as.matrix()

predictors_test_cigs <- test_smokers_scaled |>
  select(-cigs_per_day, -current_smoker) |>
  as.matrix()


complete_rows_cigs  <- complete.cases(predictors_train_cigs) & 
  !is.na(train_smokers_scaled$cigs_per_day)

y_fit_cigs          <- train_smokers_scaled$cigs_per_day[complete_rows_cigs]
predictors_fit_cigs <- predictors_train_cigs[complete_rows_cigs, ]

cat("y_fit_cigs length:", length(y_fit_cigs), "\n")
cat("predictors_fit_cigs dim:", dim(predictors_fit_cigs), "\n")
cat("Columns:", colnames(predictors_fit_cigs), "\n")

train_miss_idx_cigs        <- which(is.na(train_smokers_scaled$cigs_per_day))
predictors_train_miss_cigs <- predictors_train_cigs[train_miss_idx_cigs, ]

test_miss_idx_cigs        <- which(is.na(test_smokers_scaled$cigs_per_day))
predictors_test_miss_cigs <- predictors_test_cigs[test_miss_idx_cigs, ]

cat("Train missing cigs_per_day:", length(train_miss_idx_cigs), "\n")
cat("Test missing cigs_per_day:",  length(test_miss_idx_cigs),  "\n")

colnames(predictors_fit_cigs)

# model definition in JAGS
cigs_model_string <- "
model {
  # Likelihood — Negative Binomial
  for (i in 1:N) {
    y[i] ~ dnegbin(p[i], r)
    p[i] <- r / (r + lambda[i])
    log(lambda[i]) <- beta0 + inprod(X[i,], beta[])
  }

  # Priors
  beta0 ~ dnorm(0, 0.01)

  for (j in 1:K) {
    beta[j] ~ dnorm(0, 0.01)
  }

  # Dispersion parameter
  r ~ dgamma(0.01, 0.01)
}
"

jags_data_cigs <- list(
  y = as.integer(y_fit_cigs),
  X = predictors_fit_cigs,
  N = length(y_fit_cigs),
  K = ncol(predictors_fit_cigs)
)

cigs_model <- jags.model(
  textConnection(cigs_model_string),
  data     = jags_data_cigs,
  n.chains = 3,
  n.adapt  = 3000
)

update(cigs_model, n.iter = 10000)

posterior_cigs <- coda.samples(
  cigs_model,
  variable.names = c("beta0", "beta", "r"),
  n.iter = 15000,
  thin   = 5
)

# Diagnostics
cat("\n--- Posterior Summary ---\n")
print(summary(posterior_cigs))

cat("\nGelman-Rubin Diagnostic:\n")
print(gelman.diag(posterior_cigs))

cat("\nEffective Sample Sizes:\n")
print(effectiveSize(posterior_cigs))

# Six significant predictors (95% CI excludes zero)
params_to_plot_cigs <- c("beta[1]", "beta[2]", "beta[3]",
                         "beta[10]", "beta[11]", "beta[12]")

labels_cigs <- c("Sex", "Age", "Education",
                 "Heart Rate", "Glucose", "Ten Year CHD")

posterior_df_cigs <- as.data.frame(posterior_matrix_cigs[, params_to_plot_cigs])
colnames(posterior_df_cigs) <- labels_cigs

posterior_df_cigs |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = labels_cigs)) |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.75, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Posterior Distributions — Significant Predictors of Cigs/Day",
    subtitle = "Negative Binomial model | 95% CI excludes zero | Dashed line at zero",
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

ggsave("Figures/posterior_distributions_cigs_per_day.png")

#  Dispersion parameter r alone
data.frame(r = posterior_matrix_cigs[, "r"]) |>
  ggplot(aes(x = r)) +
  geom_density(fill = "coral", alpha = 0.75, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = mean(posterior_matrix_cigs[, "r"]),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  annotate("text",
           x    = mean(posterior_matrix_cigs[, "r"]),
           y    = Inf,
           label = paste0("Mean = ", round(mean(posterior_matrix_cigs[, "r"]), 2)),
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey20") +
  labs(
    title    = "Posterior Distribution — Dispersion Parameter (r)",
    subtitle = "Negative Binomial model | Dashed line = posterior mean",
    x        = "r",
    y        = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/posterior_distribution_dispersion_r_cigs.png")



# Posterior samples
beta0_samples_cigs <- posterior_matrix_cigs[, "beta0"]
beta_samples_cigs  <- posterior_matrix_cigs[, paste0("beta[", 1:13, "]")]
r_samples_cigs     <- posterior_matrix_cigs[, "r"]
n_samples_cigs     <- length(beta0_samples_cigs)
cat("Number of posterior samples:", n_samples_cigs, "\n")

# Posterior predictive — train missing
n_train_miss_cigs  <- nrow(predictors_train_miss_cigs)
ppc_train_cigs     <- matrix(NA, nrow = n_samples_cigs, ncol = n_train_miss_cigs)

for (s in 1:n_samples_cigs) {
  lambda_train        <- exp(beta0_samples_cigs[s] + predictors_train_miss_cigs %*% beta_samples_cigs[s, ])
  ppc_train_cigs[s, ] <- rnbinom(n_train_miss_cigs, size = r_samples_cigs[s], mu = lambda_train)
}

# Posterior predictive — test missing
n_test_miss_cigs  <- nrow(predictors_test_miss_cigs)
ppc_test_cigs     <- matrix(NA, nrow = n_samples_cigs, ncol = n_test_miss_cigs)

for (s in 1:n_samples_cigs) {
  lambda_test        <- exp(beta0_samples_cigs[s] + predictors_test_miss_cigs %*% beta_samples_cigs[s, ])
  ppc_test_cigs[s, ] <- rnbinom(n_test_miss_cigs, size = r_samples_cigs[s], mu = lambda_test)
}

# Posterior summary
imputed_train_cigs <- round(colMeans(ppc_train_cigs))
imputed_test_cigs  <- round(colMeans(ppc_test_cigs))

ci_train_cigs <- apply(ppc_train_cigs, 2, quantile, probs = c(0.025, 0.975))
ci_test_cigs  <- apply(ppc_test_cigs,  2, quantile, probs = c(0.025, 0.975))

cat("Imputed train cigs/day:\n")
print(imputed_train_cigs)
cat("\n95% Credible Intervals — train:\n")
print(round(ci_train_cigs, 1))

cat("\nImputed test cigs/day:\n")
print(imputed_test_cigs)
cat("\n95% Credible Intervals — test:\n")
print(round(ci_test_cigs, 1))


set.seed(123)
sample_idx_train <- sample(1:n_train_miss_cigs, 3)
sample_idx_test  <- sample(1:n_test_miss_cigs,  3)

plot_data_cigs <- bind_rows(
  as.data.frame(ppc_train_cigs[, sample_idx_train]) |>
    setNames(paste0("Train Obs ", train_miss_idx_cigs[sample_idx_train])) |>
    pivot_longer(everything(), names_to = "observation", values_to = "value") |>
    mutate(set = "Train"),
  
  as.data.frame(ppc_test_cigs[, sample_idx_test]) |>
    setNames(paste0("Test Obs ", test_miss_idx_cigs[sample_idx_test])) |>
    pivot_longer(everything(), names_to = "observation", values_to = "value") |>
    mutate(set = "Test")
)

mean_labels_cigs <- plot_data_cigs |>
  group_by(observation, set) |>
  summarise(mean_val = mean(value), .groups = "drop") |>
  mutate(label = paste0("Mean = ", round(mean_val, 1)))

plot_data_cigs |>
  ggplot(aes(x = value, fill = set)) +
  geom_histogram(bins = 30, alpha = 0.75, color = "white") +
  geom_vline(data = mean_labels_cigs,
             aes(xintercept = mean_val),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_text(data = mean_labels_cigs,
            aes(x = mean_val, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5,
            size = 2.8, color = "grey20", inherit.aes = FALSE) +
  facet_wrap(~observation, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  labs(
    title    = "Posterior Predictive Distributions — Missing Cigs/Day",
    subtitle = "Negative Binomial model | 3 train and 3 test observations | Dashed line = posterior mean",
    x        = "Cigarettes per Day",
    y        = "Count",
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

ggsave("Figures/posterior_predictive_distribution_missing_cigs_per_day.png")

train$cigs_per_day[train$current_smoker == 1 & is.na(train$cigs_per_day)] <- imputed_train_cigs
test$cigs_per_day[test$current_smoker == 1 & is.na(test$cigs_per_day)] <- imputed_test_cigs

cat("Missing cigs_per_day in train:", sum(is.na(train$cigs_per_day)), "\n")
cat("Missing cigs_per_day in test:",  sum(is.na(test$cigs_per_day)),  "\n")


train <- train |> select(-current_smoker) # no longer needed
test  <- test  |> select(-current_smoker) 

saveRDS(train, "Data/train.rds")
saveRDS(test,  "Data/test.rds")

# modeling missingness in bp_meds ----------------------------
train <- readRDS("Data/train.rds")
test  <- readRDS("Data/test.rds")

colnames(train)

train |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |> print()

test |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |> print()


cat("Cross-tabulation: prevalent_hyp vs bp_meds\n")
table(train$prevalent_hyp, train$bp_meds, useNA = "always")

prop.table(table(train$prevalent_hyp, train$bp_meds), margin = 1)


train |>
  filter(!is.na(bp_meds)) |>
  group_by(prevalent_hyp) |>
  summarise(
    n         = n(),
    n_bp_meds = sum(bp_meds),
    prop_bp   = mean(bp_meds)
  )


train |>
  filter(!is.na(bp_meds)) |>
  mutate(
    prevalent_hyp = factor(prevalent_hyp, labels = c("No Hypertension", "Hypertension")),
    bp_meds       = factor(bp_meds,       labels = c("No BP Meds", "On BP Meds"))
  ) |>
  ggplot(aes(x = prevalent_hyp, fill = bp_meds)) +
  geom_bar(position = "fill", alpha = 0.85, color = "white") +
  scale_fill_manual(values = c("No BP Meds" = "steelblue", "On BP Meds" = "tomato")) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title    = "BP Medication Use by Hypertension Status",
    subtitle = "All BP medication users have prevalent hypertension",
    x        = NULL,
    y        = "Proportion",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/bp_meds_vs_prevalent_hyp.png")


# Non-hypertensives → 0 (deterministic from perfect separation)
train$bp_meds[train$prevalent_hyp == 0 & is.na(train$bp_meds)] <- 0
test$bp_meds[test$prevalent_hyp   == 0 & is.na(test$bp_meds)]  <- 0

# Hypertensives → sample from observed proportion (9.46%)
set.seed(123)
n_train_hyp_miss <- sum(train$prevalent_hyp == 1 & is.na(train$bp_meds))
n_test_hyp_miss  <- sum(test$prevalent_hyp  == 1 & is.na(test$bp_meds))

cat("Hypertensives with missing bp_meds — train:", n_train_hyp_miss, "\n")
cat("Hypertensives with missing bp_meds — test:",  n_test_hyp_miss,  "\n")

# Store imputed values
imputed_train_bp <- rbinom(n_train_hyp_miss, size = 1, prob = 0.0946)
imputed_test_bp  <- rbinom(n_test_hyp_miss,  size = 1, prob = 0.0946)

cat("\nImputed bp_meds values — train:\n")
print(imputed_train_bp)
cat("n(0):", sum(imputed_train_bp == 0), "| n(1):", sum(imputed_train_bp == 1), "\n")

cat("\nImputed bp_meds values — test:\n")
print(imputed_test_bp)
cat("n(0):", sum(imputed_test_bp == 0), "| n(1):", sum(imputed_test_bp == 1), "\n")

# Assign
train$bp_meds[train$prevalent_hyp == 1 & is.na(train$bp_meds)] <- imputed_train_bp
test$bp_meds[test$prevalent_hyp   == 1 & is.na(test$bp_meds)]  <- imputed_test_bp

# Verify
cat("\nMissing bp_meds in train:", sum(is.na(train$bp_meds)), "\n")
cat("Missing bp_meds in test:",   sum(is.na(test$bp_meds)),  "\n")

# Save
saveRDS(train, "Data/train.rds")
saveRDS(test,  "Data/test.rds")



# modeling missingness in education ----------------------------
train <- readRDS("Data/train.rds")
test  <- readRDS("Data/test.rds")

colnames(train)

train |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |> print()

test |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |> print()


table(train$education, useNA = "always")
table(test$education,  useNA = "always")


# Proportions
train |>
  filter(!is.na(education)) |>
  count(education) |>
  mutate(prop = round(n / sum(n), 3))

# Missingness rate
cat("Missing education in train:", sum(is.na(train$education)), "/", nrow(train), "\n")
cat("Missing education in test:",  sum(is.na(test$education)),  "/", nrow(test),  "\n")


# Scaling
continuous_vars_edu <- c("age", "cigs_per_day", "tot_chol",
                         "pulse_pressure", "bmi", "heart_rate", "glucose")

log_vars_edu <- c("glucose", "bmi", "pulse_pressure", "tot_chol", "heart_rate")

train_means_edu <- train |>
  mutate(across(all_of(log_vars_edu), log)) |>
  summarise(across(all_of(continuous_vars_edu), ~mean(.x, na.rm = TRUE)))

train_sds_edu <- train |>
  mutate(across(all_of(log_vars_edu), log)) |>
  summarise(across(all_of(continuous_vars_edu), ~sd(.x, na.rm = TRUE)))

train_scaled_edu <- train |>
  mutate(across(all_of(log_vars_edu), log)) |>
  mutate(across(all_of(continuous_vars_edu),
                ~(.x - train_means_edu[[cur_column()]]) / train_sds_edu[[cur_column()]]))

test_scaled_edu <- test |>
  mutate(across(all_of(log_vars_edu), log)) |>
  mutate(across(all_of(continuous_vars_edu),
                ~(.x - train_means_edu[[cur_column()]]) / train_sds_edu[[cur_column()]]))

# Design matrices
predictors_train_edu <- train_scaled_edu |>
  select(-education) |>
  as.matrix()

predictors_test_edu <- test_scaled_edu |>
  select(-education) |>
  as.matrix()

#  Complete cases
complete_rows_edu  <- complete.cases(predictors_train_edu) & !is.na(train_scaled_edu$education)
y_fit_edu          <- train_scaled_edu$education[complete_rows_edu] + 1L  # JAGS needs 1-indexed
predictors_fit_edu <- predictors_train_edu[complete_rows_edu, ]

cat("y_fit_edu length:", length(y_fit_edu), "\n")
cat("predictors_fit_edu dim:", dim(predictors_fit_edu), "\n")
cat("Education levels:", table(y_fit_edu), "\n")

# Missing indices
train_miss_idx_edu        <- which(is.na(train_scaled_edu$education))
predictors_train_miss_edu <- predictors_train_edu[train_miss_idx_edu, ]

test_miss_idx_edu        <- which(is.na(test_scaled_edu$education))
predictors_test_miss_edu <- predictors_test_edu[test_miss_idx_edu, ]

cat("Train missing education:", length(train_miss_idx_edu), "\n")
cat("Test missing education:",  length(test_miss_idx_edu),  "\n")


# model definition 

edu_model_string <- "
model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dcat(p[i, ])

    # Cumulative probabilities
    for (k in 1:K_ord) {
      logit(Q[i, k]) <- alpha[k] - inprod(X[i,], beta[])
    }

    # Cell probabilities
    p[i, 1] <- Q[i, 1]
    for (k in 2:K_ord) {
      p[i, k] <- Q[i, k] - Q[i, k-1]
    }
    p[i, K_ord + 1] <- 1 - Q[i, K_ord]
  }

  # Priors — ordered cutpoints
  alpha[1] ~ dnorm(0, 0.01)
  for (k in 2:K_ord) {
    alpha[k] <- alpha[k-1] + exp(delta[k])
    delta[k] ~ dnorm(0, 0.01)
  }

  # Priors — coefficients
  for (j in 1:K_pred) {
    beta[j] ~ dnorm(0, 0.01)
  }
}
"
jags_data_edu <- list(
  y     = as.integer(y_fit_edu),
  X     = predictors_fit_edu,
  N     = length(y_fit_edu),
  K_ord  = 3,                      
  K_pred = ncol(predictors_fit_edu)
)


edu_model <- jags.model(
  textConnection(edu_model_string),
  data     = jags_data_edu,
  n.chains = 3,
  n.adapt  = 1000 
)

update(edu_model, n.iter = 3000)

posterior_edu <- coda.samples(
  edu_model,
  variable.names = c("alpha", "beta"),
  n.iter = 10000, 
  thin   = 3         
)


cat("\n--- Posterior Summary ---\n")
print(summary(posterior_edu))

cat("\nGelman-Rubin Diagnostic:\n")
print(gelman.diag(posterior_edu))

cat("\nEffective Sample Sizes:\n")
print(effectiveSize(posterior_edu))


colnames(predictors_fit_edu)


posterior_matrix_edu <- as.matrix(posterior_edu)

# Significant predictors only
params_to_plot_edu <- c("alpha[1]", "alpha[2]", "alpha[3]",
                        "beta[2]", "beta[6]", "beta[8]",
                        "beta[9]", "beta[10]", "beta[13]")

labels_edu <- c("Cutpoint 1", "Cutpoint 2", "Cutpoint 3",
                "Age", "Prev. Hypertension", "Total Cholesterol",
                "BMI", "Heart Rate", "Pulse Pressure")

# Build data frame
posterior_df_edu <- as.data.frame(posterior_matrix_edu[, params_to_plot_edu])
colnames(posterior_df_edu) <- labels_edu

# Separate cutpoints and betas for clean plotting
# Plot 1 — significant betas
posterior_df_edu |>
  select(-starts_with("Cutpoint")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = labels_edu[4:9])) |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.75, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Posterior Distributions — Significant Predictors of Education",
    subtitle = "Ordinal logistic model | 95% CI excludes zero | Dashed line at zero",
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

ggsave("Figures/posterior_distributions_education.png")

# Plot 2 — cutpoints
posterior_df_edu |>
  select(starts_with("Cutpoint")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = labels_edu[1:3])) |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.75, color = "white", linewidth = 0.3) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title    = "Posterior Distributions — Ordinal Cutpoints",
    subtitle = "Thresholds separating education categories 1→2→3→4",
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

ggsave("Figures/posterior_distributions_education_cutpoints.png")


table(train$education)

# Visualize cutpoints on latent scale
cutpoints <- c(alpha1 = -0.24, alpha2 = 1.06, alpha3 = 2.16)

# Cumulative probabilities at each cutpoint
cum_probs <- plogis(cutpoints)

# Cell probabilities
cell_probs <- c(
  cum_probs[1],
  cum_probs[2] - cum_probs[1],
  cum_probs[3] - cum_probs[2],
  1 - cum_probs[3]
)

cat("Cumulative probabilities:\n")
print(round(cum_probs, 3))
cat("\nCell probabilities:\n")
print(round(cell_probs, 3))

# plot latent scale
x <- seq(-5, 5, length.out = 1000)
y <- dlogis(x)

plot_df <- data.frame(x = x, y = y)

ggplot(plot_df, aes(x = x, y = y)) +
  geom_line(linewidth = 0.8, color = "grey30") +
  
  # Shaded regions for each category
  geom_ribbon(data = subset(plot_df, x <= cutpoints[1]),
              aes(ymin = 0, ymax = y), fill = "#E41A1C", alpha = 0.4) +
  geom_ribbon(data = subset(plot_df, x > cutpoints[1] & x <= cutpoints[2]),
              aes(ymin = 0, ymax = y), fill = "#377EB8", alpha = 0.4) +
  geom_ribbon(data = subset(plot_df, x > cutpoints[2] & x <= cutpoints[3]),
              aes(ymin = 0, ymax = y), fill = "#4DAF4A", alpha = 0.4) +
  geom_ribbon(data = subset(plot_df, x > cutpoints[3]),
              aes(ymin = 0, ymax = y), fill = "#984EA3", alpha = 0.4) +
  
  # Cutpoint lines
  geom_vline(xintercept = cutpoints, linetype = "dashed",
             color = "grey20", linewidth = 0.5) +
  
  # Cutpoint labels
  annotate("text", x = cutpoints[1], y = max(y) * 1.05,
           label = paste0("α_1 = ", cutpoints[1]),
           hjust = 1.1, size = 3.5, color = "grey20") +
  annotate("text", x = cutpoints[2], y = max(y) * 1.05,
           label = paste0("α_2 = ", cutpoints[2]),
           hjust = 1.1, size = 3.5, color = "grey20") +
  annotate("text", x = cutpoints[3], y = max(y) * 1.05,
           label = paste0("α_3 = ", cutpoints[3]),
           hjust = 1.1, size = 3.5, color = "grey20") +
  
  # Category labels with probabilities
  annotate("text", x = -2.5,  y = 0.05,
           label = paste0("Some HS\n", round(cell_probs[1]*100, 1), "%"),
           size = 3, color = "#E41A1C", fontface = "bold") +
  annotate("text", x = 0.4,   y = 0.05,
           label = paste0("HS Grad\n", round(cell_probs[2]*100, 1), "%"),
           size = 3, color = "#377EB8", fontface = "bold") +
  annotate("text", x = 1.6,   y = 0.15,
           label = paste0("Some College\n", round(cell_probs[3]*100, 1), "%"),
           size = 3, color = "#4DAF4A", fontface = "bold") +
  annotate("text", x = 3.5,   y = 0.05,
           label = paste0("College+\n", round(cell_probs[4]*100, 1), "%"),
           size = 3, color = "#984EA3", fontface = "bold") +
  
  labs(
    title    = "Latent Scale — Ordinal Cutpoints for Education",
    subtitle = "Shaded areas represent probability mass in each education category",
    x        = "Latent scale",
    y        = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/education_latent_scale_cutpoints.png")


# imputation of education data

# Posterior samples
alpha_samples <- posterior_matrix_edu[, c("alpha[1]", "alpha[2]", "alpha[3]")]
beta_samples_edu  <- posterior_matrix_edu[, paste0("beta[", 1:13, "]")]
n_samples_edu     <- nrow(posterior_matrix_edu)
cat("Number of posterior samples:", n_samples_edu, "\n")

# Function to compute category probabilities
get_cat_probs <- function(alpha, lin_pred) {
  K <- length(alpha) + 1
  cum_p <- plogis(alpha - lin_pred)
  p <- c(cum_p[1],
         cum_p[2] - cum_p[1],
         cum_p[3] - cum_p[2],
         1 - cum_p[3])
  return(p)
}

# Posterior predictive — train missing
n_train_miss_edu  <- nrow(predictors_train_miss_edu)
ppc_train_edu     <- matrix(NA, nrow = n_samples_edu, ncol = n_train_miss_edu)

for (s in 1:n_samples_edu) {
  for (i in 1:n_train_miss_edu) {
    lin_pred <- sum(predictors_train_miss_edu[i, ] * beta_samples_edu[s, ])
    probs    <- get_cat_probs(alpha_samples[s, ], lin_pred)
    ppc_train_edu[s, i] <- sample(1:4, size = 1, prob = probs)
  }
}

# Posterior predictive — test missing
n_test_miss_edu  <- nrow(predictors_test_miss_edu)
ppc_test_edu     <- matrix(NA, nrow = n_samples_edu, ncol = n_test_miss_edu)

for (s in 1:n_samples_edu) {
  for (i in 1:n_test_miss_edu) {
    lin_pred <- sum(predictors_test_miss_edu[i, ] * beta_samples_edu[s, ])
    probs    <- get_cat_probs(alpha_samples[s, ], lin_pred)
    ppc_test_edu[s, i] <- sample(1:4, size = 1, prob = probs)
  }
}

# Posterior summary — most frequent category (mode)
imputed_train_edu <- apply(ppc_train_edu, 2, function(x) as.integer(names(which.max(table(x)))))
imputed_test_edu  <- apply(ppc_test_edu,  2, function(x) as.integer(names(which.max(table(x)))))

# Back to 0-indexed
imputed_train_edu <- imputed_train_edu - 1L
imputed_test_edu  <- imputed_test_edu  - 1L

cat("Imputed train education:\n")
print(table(imputed_train_edu))

cat("\nImputed test education:\n")
print(table(imputed_test_edu))

# Sanity check — should match observed proportions roughly
cat("\nObserved train proportions:\n")
print(prop.table(table(train$education)))

cat("\nImputed train proportions:\n")
print(prop.table(table(imputed_train_edu)))


# posterior predictive distribution

# Random sample of 3 train and 3 test
set.seed(123)
sample_idx_train_edu <- sample(1:n_train_miss_edu, 3)
sample_idx_test_edu  <- sample(1:n_test_miss_edu,  3)

# Compute posterior category probabilities
prob_train_df <- lapply(sample_idx_train_edu, function(i) {
  counts <- table(factor(ppc_train_edu[, i], levels = 1:4))
  props  <- as.numeric(counts) / n_samples_edu
  data.frame(
    observation = paste0("Train Obs ", train_miss_idx_edu[i]),
    category    = factor(0:3, labels = c("Some HS", "HS Grad",
                                         "Some College", "College+")),
    probability = props,
    set         = "Train"
  )
}) |> bind_rows()

prob_test_df <- lapply(sample_idx_test_edu, function(i) {
  counts <- table(factor(ppc_test_edu[, i], levels = 1:4))
  props  <- as.numeric(counts) / n_samples_edu
  data.frame(
    observation = paste0("Test Obs ", test_miss_idx_edu[i]),
    category    = factor(0:3, labels = c("Some HS", "HS Grad",
                                         "Some College", "College+")),
    probability = props,
    set         = "Test"
  )
}) |> bind_rows()

prob_df <- bind_rows(prob_train_df, prob_test_df)

# Plot
prob_df |>
  ggplot(aes(x = category, y = probability, fill = set)) +
  geom_col(alpha = 0.85, color = "white") +
  geom_text(aes(label = paste0(round(probability * 100, 1), "%")),
            vjust = -0.4, size = 2.8, color = "grey20") +
  facet_wrap(~observation, ncol = 3) +
  scale_fill_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title    = "Posterior Predictive Histogram — Missing Education",
    subtitle = "3 train and 3 test observations | Bars show posterior probability per category",
    x        = NULL,
    y        = "Posterior Probability",
    fill     = "Dataset"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  )

ggsave("Figures/posterior_predictive_uncertainty_education.png")

# completing the imputation
train$education[train_miss_idx_edu] <- imputed_train_edu
test$education[test_miss_idx_edu]   <- imputed_test_edu


cat("Missing education in train:", sum(is.na(train$education)), "\n")
cat("Missing education in test:",  sum(is.na(test$education)),  "\n")


cat("\nObserved train proportions:\n")
print(prop.table(table(train$education)))

saveRDS(train, "Data/train.rds")
saveRDS(test,  "Data/test.rds")

# modeling missingness in glucose ----------------------------

train <- readRDS("Data/train.rds")
test  <- readRDS("Data/test.rds")

colnames(train)

train |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |> print()

test |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |> print()


# Scaling
continuous_vars_gluc <- c("age", "cigs_per_day", "tot_chol",
                          "pulse_pressure", "bmi", "heart_rate", "glucose")

log_vars_gluc <- c("bmi", "pulse_pressure", "tot_chol", "heart_rate", "glucose")

train_means_gluc <- train |>
  mutate(across(all_of(log_vars_gluc), log)) |>
  summarise(across(all_of(continuous_vars_gluc), ~mean(.x, na.rm = TRUE)))

train_sds_gluc <- train |>
  mutate(across(all_of(log_vars_gluc), log)) |>
  summarise(across(all_of(continuous_vars_gluc), ~sd(.x, na.rm = TRUE)))

train_scaled_gluc <- train |>
  mutate(across(all_of(log_vars_gluc), log)) |>
  mutate(across(all_of(continuous_vars_gluc),
                ~(.x - train_means_gluc[[cur_column()]]) / train_sds_gluc[[cur_column()]]))

test_scaled_gluc <- test |>
  mutate(across(all_of(log_vars_gluc), log)) |>
  mutate(across(all_of(continuous_vars_gluc),
                ~(.x - train_means_gluc[[cur_column()]]) / train_sds_gluc[[cur_column()]]))

# Design matrices
predictors_train_gluc <- train_scaled_gluc |>
  select(-glucose) |>
  as.matrix()

predictors_test_gluc <- test_scaled_gluc |>
  select(-glucose) |>
  as.matrix()

# Complete cases
complete_rows_gluc  <- complete.cases(predictors_train_gluc) & !is.na(train_scaled_gluc$glucose)
y_fit_gluc          <- train_scaled_gluc$glucose[complete_rows_gluc]
predictors_fit_gluc <- predictors_train_gluc[complete_rows_gluc, ]

cat("y_fit_gluc length:", length(y_fit_gluc), "\n")
cat("predictors_fit_gluc dim:", dim(predictors_fit_gluc), "\n")
cat("Columns:", colnames(predictors_fit_gluc), "\n")

# Missing indices
train_miss_idx_gluc        <- which(is.na(train_scaled_gluc$glucose))
predictors_train_miss_gluc <- predictors_train_gluc[train_miss_idx_gluc, ]

test_miss_idx_gluc        <- which(is.na(test_scaled_gluc$glucose))
predictors_test_miss_gluc <- predictors_test_gluc[test_miss_idx_gluc, ]

cat("\nTrain missing glucose:", length(train_miss_idx_gluc), "\n")
cat("Test missing glucose:",  length(test_miss_idx_gluc),  "\n")


# JAGS model for glucose missingnes
gluc_model_string <- "
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

jags_data_gluc <- list(
  y = y_fit_gluc,
  X = predictors_fit_gluc,
  N = length(y_fit_gluc),
  K = ncol(predictors_fit_gluc)
)

gluc_model <- jags.model(
  textConnection(gluc_model_string),
  data     = jags_data_gluc,
  n.chains = 3,
  n.adapt  = 3000
)

update(gluc_model, n.iter = 10000)

posterior_gluc <- coda.samples(
  gluc_model,
  variable.names = c("beta0", "beta", "sigma"),
  n.iter = 20000,
  thin   = 5
)

# Diagnostics
cat("\n--- Posterior Summary ---\n")
print(summary(posterior_gluc))

cat("\nGelman-Rubin Diagnostic:\n")
print(gelman.diag(posterior_gluc))

cat("\nEffective Sample Sizes:\n")
print(effectiveSize(posterior_gluc))

colnames(predictors_fit_gluc)

# posterior distribution 
posterior_matrix_gluc <- as.matrix(posterior_gluc)

params_to_plot_gluc <- c("beta[4]", "beta[8]", "beta[11]",
                         "beta[12]", "beta[13]", "sigma")

labels_gluc <- c("Cigs Per Day", "Diabetes", "Heart Rate",
                 "Ten Year CHD", "Pulse Pressure", "Sigma")

posterior_df_gluc <- as.data.frame(posterior_matrix_gluc[, params_to_plot_gluc])
colnames(posterior_df_gluc) <- labels_gluc

posterior_df_gluc |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = labels_gluc)) |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.75, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Posterior Distributions — Significant Predictors of Glucose",
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

ggsave("Figures/posterior_distributions_glucose.png")

# posterior predictive distribution for missing values

# Posterior samples
beta0_samples_gluc <- posterior_matrix_gluc[, "beta0"]
beta_samples_gluc  <- posterior_matrix_gluc[, paste0("beta[", 1:13, "]")]
sigma_samples_gluc <- posterior_matrix_gluc[, "sigma"]
n_samples_gluc     <- length(beta0_samples_gluc)
cat("Number of posterior samples:", n_samples_gluc, "\n")

# Posterior predictive — train missing
n_train_miss_gluc  <- nrow(predictors_train_miss_gluc)
ppc_train_gluc     <- matrix(NA, nrow = n_samples_gluc, ncol = n_train_miss_gluc)

for (s in 1:n_samples_gluc) {
  mu_train           <- beta0_samples_gluc[s] + predictors_train_miss_gluc %*% beta_samples_gluc[s, ]
  ppc_train_gluc[s, ] <- rnorm(n_train_miss_gluc, mean = mu_train, sd = sigma_samples_gluc[s])
}

# Posterior predictive — test missing
n_test_miss_gluc  <- nrow(predictors_test_miss_gluc)
ppc_test_gluc     <- matrix(NA, nrow = n_samples_gluc, ncol = n_test_miss_gluc)

for (s in 1:n_samples_gluc) {
  mu_test            <- beta0_samples_gluc[s] + predictors_test_miss_gluc %*% beta_samples_gluc[s, ]
  ppc_test_gluc[s, ] <- rnorm(n_test_miss_gluc, mean = mu_test, sd = sigma_samples_gluc[s])
}

# Back-transform: scaled log → original scale
ppc_train_gluc_orig <- exp(ppc_train_gluc * train_sds_gluc[["glucose"]] + train_means_gluc[["glucose"]])
ppc_test_gluc_orig  <- exp(ppc_test_gluc  * train_sds_gluc[["glucose"]] + train_means_gluc[["glucose"]])

# Posterior summary 
imputed_train_gluc <- round(colMeans(ppc_train_gluc_orig))
imputed_test_gluc  <- round(colMeans(ppc_test_gluc_orig))

ci_train_gluc <- apply(ppc_train_gluc_orig, 2, quantile, probs = c(0.025, 0.975))
ci_test_gluc  <- apply(ppc_test_gluc_orig,  2, quantile, probs = c(0.025, 0.975))

cat("Imputed train glucose (original scale):\n")
print(imputed_train_gluc)
cat("\n95% Credible Intervals — train:\n")
print(round(ci_train_gluc, 1))

cat("\nImputed test glucose (original scale):\n")
print(imputed_test_gluc)
cat("\n95% Credible Intervals — test:\n")
print(round(ci_test_gluc, 1))


# Select 3 train and 3 random test
set.seed(123)
high_val_idx     <- 139  # imputed value = 153
other_train_idx  <- sample(setdiff(1:n_train_miss_gluc, high_val_idx), 2)
sample_idx_train_gluc <- c(high_val_idx, other_train_idx)
sample_idx_test_gluc  <- sample(1:n_test_miss_gluc, 3)

# Build plot data
plot_data_gluc <- bind_rows(
  as.data.frame(ppc_train_gluc_orig[, sample_idx_train_gluc]) |>
    setNames(paste0("Train Obs ", train_miss_idx_gluc[sample_idx_train_gluc],
                    ifelse(sample_idx_train_gluc == high_val_idx, " ★", ""))) |>
    pivot_longer(everything(), names_to = "observation", values_to = "value") |>
    mutate(set = "Train"),
  
  as.data.frame(ppc_test_gluc_orig[, sample_idx_test_gluc]) |>
    setNames(paste0("Test Obs ", test_miss_idx_gluc[sample_idx_test_gluc])) |>
    pivot_longer(everything(), names_to = "observation", values_to = "value") |>
    mutate(set = "Test")
)

mean_labels_gluc <- plot_data_gluc |>
  group_by(observation, set) |>
  summarise(mean_val = mean(value), .groups = "drop") |>
  mutate(label = paste0("Mean = ", round(mean_val, 1)))

# Plot
plot_data_gluc |>
  ggplot(aes(x = value, fill = set)) +
  geom_density(alpha = 0.75, color = "white", linewidth = 0.3) +
  geom_vline(data = mean_labels_gluc,
             aes(xintercept = mean_val),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_text(data = mean_labels_gluc,
            aes(x = mean_val, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5,
            size = 2.8, color = "grey20", inherit.aes = FALSE) +
  facet_wrap(~observation, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  labs(
    title    = "Posterior Predictive Distributions — Missing Glucose",
    subtitle = "★ High value observation (diabetic) | Dashed line = posterior mean",
    x        = "Glucose (mg/dL)",
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

ggsave("Figures/posterior_predictive_distribution_missing_glucose.png")


# adding the imputed values

train$glucose[train_miss_idx_gluc] <- imputed_train_gluc
test$glucose[test_miss_idx_gluc]   <- imputed_test_gluc

cat("Missing glucose in train:", sum(is.na(train$glucose)), "\n")
cat("Missing glucose in test:",  sum(is.na(test$glucose)),  "\n")

cat("\nRemaining missingness in train:\n")
train |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |> print()

cat("\nRemaining missingness in test:\n")
test |> summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0) |> print()


saveRDS(train, "Data/train.rds")
saveRDS(test,  "Data/test.rds")


cat("Train dimensions:", dim(train), "\n")
cat("Test dimensions:",  dim(test),  "\n")