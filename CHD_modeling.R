library(tidyverse)
library(rjags)
library(coda)


# Data loading and some sanity checks ----------------------------------------

train <- readRDS("Data/train.rds")
test  <- readRDS("Data/test.rds")

# Verify no missingness
cat("Missing in train:", sum(is.na(train)), "\n")
cat("Missing in test:",  sum(is.na(test)),  "\n")

colnames(train)

# Define variable groups
continuous_vars <- c("age", "cigs_per_day", "tot_chol",
                     "pulse_pressure", "bmi", "heart_rate", "glucose")

cat_vars <- c("sex", "education", "bp_meds", "prevalent_stroke",
              "prevalent_hyp", "diabetes", "ten_year_chd")

#  Continuous — train vs test overlay
bind_rows(
  train |> select(all_of(continuous_vars)) |> mutate(set = "Train"),
  test  |> select(all_of(continuous_vars)) |> mutate(set = "Test")
) |>
  pivot_longer(-set, names_to = "variable", values_to = "value") |>
  ggplot(aes(x = value, fill = set, color = set)) +
  geom_density(alpha = 0.4, linewidth = 0.6) +
  facet_wrap(~variable, scales = "free", ncol = 4) +
  scale_fill_manual(values  = c("Train" = "steelblue", "Test" = "tomato")) +
  scale_color_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  labs(
    title    = "Distribution of Continuous Variables — Original Scale",
    subtitle = "Train vs Test overlay",
    x        = NULL,
    y        = "Density",
    fill     = "Dataset",
    color    = "Dataset"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/continuous_distributions_original_scale_after_imputation.png")

# Categorical — train vs test side by side
bind_rows(
  train |> select(all_of(cat_vars)) |> mutate(set = "Train"),
  test  |> select(all_of(cat_vars)) |> mutate(set = "Test")
) |>
  pivot_longer(-set, names_to = "variable", values_to = "value") |>
  mutate(value = factor(value)) |>
  ggplot(aes(x = value, fill = set)) +
  geom_bar(position = "dodge", alpha = 0.85, color = "white") +
  facet_wrap(~variable, scales = "free", ncol = 4) +
  scale_fill_manual(values = c("Train" = "steelblue", "Test" = "tomato")) +
  labs(
    title    = "Distribution of Categorical Variables",
    subtitle = "Train vs Test comparison",
    x        = NULL,
    y        = "Count",
    fill     = "Dataset"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/categorical_distributions_after_imputation.png")

cat("Train — prevalent_stroke:\n")
train |> count(prevalent_stroke) |> mutate(prop = round(n / sum(n), 3)) |> print()

cat("\nTest — prevalent_stroke:\n")
test |> count(prevalent_stroke) |> mutate(prop = round(n / sum(n), 3)) |> print()



# setting up data for JAGS modeling

continuous_vars <- c("age", "cigs_per_day", "tot_chol",
                     "pulse_pressure", "bmi", "heart_rate", "glucose")
log_vars        <- c("glucose", "bmi", "pulse_pressure", "tot_chol", "heart_rate")
cat_vars        <- c("sex", "education", "bp_meds",
                     "prevalent_stroke", "prevalent_hyp", "diabetes", "ten_year_chd")


train_means <- train |>
  mutate(across(all_of(log_vars), log)) |>
  summarise(across(all_of(continuous_vars), ~ mean(.x, na.rm = TRUE)))

train_sds <- train |>
  mutate(across(all_of(log_vars), log)) |>
  summarise(across(all_of(continuous_vars), ~ sd(.x, na.rm = TRUE)))


# log-transform and scaling
train_scaled <- train |>
  mutate(across(all_of(log_vars), log)) |>
  mutate(across(all_of(continuous_vars),
                ~ (.x - train_means[[cur_column()]]) / train_sds[[cur_column()]]))

test_scaled <- test |>
  mutate(across(all_of(log_vars), log)) |>
  mutate(across(all_of(continuous_vars),
                ~ (.x - train_means[[cur_column()]]) / train_sds[[cur_column()]]))


Y_train <- train_scaled$ten_year_chd
Y_test  <- test_scaled$ten_year_chd

predictor_vars <- c("age", "cigs_per_day", "tot_chol",
                    "pulse_pressure", "bmi", "heart_rate", "glucose",
                    "sex", "education", "bp_meds",
                    "prevalent_stroke", "prevalent_hyp", "diabetes")


X_train <- train_scaled |>
  select(all_of(predictor_vars)) |>
  as.matrix()

X_test <- test_scaled |>
  select(all_of(predictor_vars)) |>
  as.matrix()


N_train <- nrow(X_train)
N_test  <- nrow(X_test)
P       <- ncol(X_train)

cat("X_train:", N_train, "x", P, "\n")
cat("X_test: ", N_test,  "x", P, "\n")
cat("Columns:", colnames(X_train), "\n")


# check the determinant of X'X for train and test
det(t(X_train) %*% X_train)
det(t(X_test) %*% X_test)



# Feature selection using SSVS ----------------------------------
ssvs_model_string <- "
model {

  for (i in 1:N) {
    Y[i] ~ dbern(pi[i])
    logit(pi[i]) <- alpha + inprod(X[i,], beta[])
  }

  for (j in 1:P) {
    beta[j] ~ dnorm(0, tau / lambda[j])
    
    gamma[j] ~ dbern(0.5)
    
    lambda[j] <- gamma[j] * c1 + (1 - gamma[j]) * c0
  }

  alpha ~ dnorm(0, 0.01)

  tau ~ dgamma(0.1, 0.1)

  c0 <- 0.01
  c1 <- 10
}
"

jags_data <- list(
  Y = Y_train,
  X = X_train,
  N = N_train,
  P = P
)

inits <- function() list(
  alpha = 0,
  beta  = rnorm(P, 0, 0.1),
  gamma = rbinom(P, 1, 0.5),
  tau   = 1
)
params <- c("alpha", "beta", "gamma", "tau")

jags_model <- jags.model(
  textConnection(ssvs_model_string),
  data     = jags_data,
  inits    = inits,
  n.chains = 3,
  n.adapt  = 3000
)

update(jags_model, n.iter = 10000)

ssvs_samples <- coda.samples(
  jags_model,
  variable.names = params,
  n.iter         = 20000,
  thin           = 5
)

# Diagnostics
cat("\n--- Posterior Summary ---\n")
print(summary(ssvs_samples))

cat("\nGelman-Rubin Diagnostic:\n")
gelman.diag(ssvs_samples)

cat("\nEffective Sample Sizes:\n")
print(effectiveSize(ssvs_samples))


gamma_samples <- as.matrix(ssvs_samples)[, grep("^gamma", colnames(as.matrix(ssvs_samples)))]
pip <- colMeans(gamma_samples)
pip
data.frame(variable = colnames(X_train), PIP = pip)

# plotting pip 
data.frame(variable = colnames(X_train), PIP = pip) |>
  arrange(PIP) |>
  mutate(variable = factor(variable, levels = variable)) |>
  ggplot(aes(x = PIP, y = variable)) +
  geom_col(fill = "steelblue", alpha = 0.85, width = 0.6) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "tomato", linewidth = 0.5) +
  geom_text(aes(label = round(PIP, 3)), hjust = -0.15, size = 3.5, color = "grey30") +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "Posterior Inclusion Probabilities (PIP)",
    subtitle = "Dashed line at 0.5 threshold",
    x        = "PIP",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/pip_barplot.png", width = 7, height = 5)


# looking at the most common two sampled models
samples_mat <- as.matrix(ssvs_samples)
gamma_samples <- samples_mat[, grep("^gamma", colnames(samples_mat))]
gamma_str <- apply(gamma_samples, 1, paste0, collapse = "")
model_freq <- sort(table(gamma_str), decreasing = TRUE)

top2 <- names(model_freq)[1:2]
top2


decode_model <- function(model_str, var_names) {
  as.logical(as.numeric(strsplit(model_str, "")[[1]])) |> 
    (\(x) var_names[x])()
}

lapply(top2, decode_model, var_names = colnames(X_train))


model_prop <- sort(model_freq / sum(model_freq), decreasing = TRUE)

data.frame(
  model = sapply(top2, decode_model, var_names = colnames(X_train),
                 simplify = FALSE) |>
    sapply(\(x) paste(x, collapse = " + ")),
  proportion = as.numeric(model_prop[top2])
)



# plotting the beta posteriror distribution for all of the top 5 included variables
top5_vars <- data.frame(variable = colnames(X_train), PIP = pip) |>
  arrange(desc(PIP)) |>
  slice(1:5)

top5_vars


beta_samples <- as.matrix(ssvs_samples)[, grep("^beta", colnames(as.matrix(ssvs_samples)))]
colnames(beta_samples) <- colnames(X_train)

beta_top5_df <- beta_samples[, top5_vars$variable] |>
  as.data.frame() |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  left_join(top5_vars, by = "variable") |>
  mutate(variable = factor(variable, levels = top5_vars$variable),  # sorted by PIP
         label = paste0(variable, "\n(PIP = ", round(PIP, 3), ")"))

ggplot(beta_top5_df, aes(x = value, fill = variable, color = variable)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2")  +
  scale_color_brewer(palette = "Set2") +
  labs(
    title    = "Posterior Distribution of β - Top 5 Variables by PIP",
    subtitle = "Vertical dashed line at zero",
    x        = "β value",
    y        = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(color = "grey40", size = 10)
  )


ggsave("Figures/beta_posterior_top5_pip.png", width = 10, height = 6)




