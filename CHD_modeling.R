library(tidyverse)
library(rjags)
library(coda)
library(caret)
library(pROC)

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



# Model 1: age + sex + cigs_per_day + prevalent_hyp + pulse_pressure ---------------------------------

X_m1_train <- X_train[, c("age", "sex", "cigs_per_day", "prevalent_hyp", "pulse_pressure")]

mod1_string <- "
model {
  for (i in 1:n) {
    Y[i]          ~ dbern(pi[i])
    logit(pi[i]) <- beta[1]
                   + X[i,1]*beta[2] + X[i,2]*beta[3]
                   + X[i,3]*beta[4] + X[i,4]*beta[5] + X[i,5]*beta[6]
    like[i]      <- dbin(Y[i], pi[i], 1)
  }
  for (j in 1:6) { beta[j] ~ dnorm(0, 0.01) }
}
"

data1 <- list(Y = Y_train, X = X_m1_train, n = N_train)

inits1 <- function() list(beta = rnorm(6, 0, 0.1))

model1 <- jags.model(
  textConnection(mod1_string),
  data     = data1,
  inits    = inits1,
  n.chains = 3,
  n.adapt  = 3000
)

update(model1, n.iter = 10000)

samps1 <- coda.samples(
  model1,
  variable.names = c("beta", "like"),
  n.iter         = 20000,
  thin           = 5
)

# --- DIC ---
DIC1 <- dic.samples(model1, n.iter = 20000, thin = 5)

# --- WAIC ---
like_cols <- grep("^like", colnames(samps1[[1]]))
like1     <- rbind(samps1[[1]][, like_cols],
                   samps1[[2]][, like_cols],
                   samps1[[3]][, like_cols])
fbar1 <- colMeans(like1)
Pw1   <- sum(apply(log(like1), 2, var))
WAIC1 <- -2 * sum(log(fbar1)) + 2 * Pw1

cat("=== Model 1: age + sex + cigs_per_day + prevalent_hyp + pulse_pressure ===\n")
cat("DIC:\n");  print(DIC1)
cat("WAIC:", round(WAIC1, 2), "\n")
cat("pW (WAIC):", round(Pw1, 2), "\n")

# --- Diagnostics ---
beta_samps1 <- samps1[, grep("^beta", colnames(samps1[[1]]))]

cat("\nGelman-Rubin:\n")
print(gelman.diag(beta_samps1))

cat("\nEffective Sample Sizes:\n")
print(effectiveSize(beta_samps1))

cat("\nPosterior Summary:\n")
print(summary(beta_samps1))


png("Figures/trace_plots_model_1_beta_1to3.png", width = 1400, height = 1200, res = 120)
par(mfrow = c(3, 2))
for (j in 1:3) {
  traceplot(beta_samps1[, j], main = paste0("beta[", j, "] trace"))
  densplot(beta_samps1[, j],  main = paste0("beta[", j, "] density"))
}
dev.off()

png("Figures/trace_plots_model_1_beta_4to6.png", width = 1400, height = 1200, res = 120)
par(mfrow = c(3, 2))
for (j in 4:6) {
  traceplot(beta_samps1[, j], main = paste0("beta[", j, "] trace"))
  densplot(beta_samps1[, j],  main = paste0("beta[", j, "] density"))
}
dev.off()


# Model 2: age + sex + cigs_per_day + prevalent_hyp + pulse_pressure + glucose ----------

X_m2_train <- X_train[, c("age", "sex", "cigs_per_day", "prevalent_hyp", "pulse_pressure", "glucose")]

mod2_string <- "
model {
  for (i in 1:n) {
    Y[i]          ~ dbern(pi[i])
    logit(pi[i]) <- beta[1]
                   + X[i,1]*beta[2] + X[i,2]*beta[3]
                   + X[i,3]*beta[4] + X[i,4]*beta[5] + X[i,5]*beta[6]
                   + X[i,6]*beta[7]
    like[i]      <- dbin(Y[i], pi[i], 1)
  }
  for (j in 1:7) { beta[j] ~ dnorm(0, 0.01) }
}
"

data2 <- list(Y = Y_train, X = X_m2_train, n = N_train)

inits2 <- function() list(beta = rnorm(7, 0, 0.1))

model2 <- jags.model(
  textConnection(mod2_string),
  data     = data2,
  inits    = inits2,
  n.chains = 3,
  n.adapt  = 3000
)

update(model2, n.iter = 10000)

samps2 <- coda.samples(
  model2,
  variable.names = c("beta", "like"),
  n.iter         = 20000,
  thin           = 5
)

# --- DIC ---
DIC2 <- dic.samples(model2, n.iter = 20000, thin = 5)

# --- WAIC ---
like_cols2 <- grep("^like", colnames(samps2[[1]]))
like2      <- rbind(samps2[[1]][, like_cols2],
                    samps2[[2]][, like_cols2],
                    samps2[[3]][, like_cols2])
fbar2 <- colMeans(like2)
Pw2   <- sum(apply(log(like2), 2, var))
WAIC2 <- -2 * sum(log(fbar2)) + 2 * Pw2

cat("=== Model 2: age + sex + cigs_per_day + prevalent_hyp + pulse_pressure + glucose ===\n")
cat("DIC:\n");  print(DIC2)
cat("WAIC:", round(WAIC2, 2), "\n")
cat("pW (WAIC):", round(Pw2, 2), "\n")

# --- Diagnostics ---
beta_samps2 <- samps2[, grep("^beta", colnames(samps2[[1]]))]

cat("\nGelman-Rubin:\n")
print(gelman.diag(beta_samps2))

cat("\nEffective Sample Sizes:\n")
print(effectiveSize(beta_samps2))

cat("\nPosterior Summary:\n")
print(summary(beta_samps2))


png("Figures/trace_plots_model_2_beta_1to3.png", width = 1400, height = 1200, res = 120)
par(mfrow = c(3, 2))
for (j in 1:3) {
  traceplot(beta_samps2[, j], main = paste0("beta[", j, "] trace"))
  densplot(beta_samps2[, j],  main = paste0("beta[", j, "] density"))
}
dev.off()

png("Figures/trace_plots_model_2_beta_4to6.png", width = 1400, height = 1200, res = 120)
par(mfrow = c(3, 2))
for (j in 4:6) {
  traceplot(beta_samps2[, j], main = paste0("beta[", j, "] trace"))
  densplot(beta_samps2[, j],  main = paste0("beta[", j, "] density"))
}
dev.off()

png("Figures/trace_plots_model_2_beta_7.png", width = 1400, height = 500, res = 120)
par(mfrow = c(1, 2))
traceplot(beta_samps2[, 7], main = "beta[7] trace")
densplot(beta_samps2[, 7],  main = "beta[7] density")
dev.off()





# Comparison between model_1 and model_2 --------------------------
cat("\n=== Model Comparison ===\n")
cat("        WAIC     pW\n")
cat("Model1:", round(WAIC1, 2), "  ", round(Pw1, 2), "\n")
cat("Model2:", round(WAIC2, 2), "  ", round(Pw2, 2), "\n")
cat("\nDIC Model 1:\n"); print(DIC1)
cat("\nDIC Model 2:\n"); print(DIC2)



# Extract beta samples for both models
beta_names_m1 <- c("intercept", "age", "sex", "cigs_per_day", 
                   "prevalent_hyp", "pulse_pressure")

beta_names_m2 <- c("intercept", "age", "sex", "cigs_per_day",
                   "prevalent_hyp", "pulse_pressure", "glucose")

# Combine all chains
beta_mat1 <- rbind(samps1[[1]], samps1[[2]], samps1[[3]])[, grep("^beta", colnames(samps1[[1]]))]
beta_mat2 <- rbind(samps2[[1]], samps2[[2]], samps2[[3]])[, grep("^beta", colnames(samps2[[1]]))]

colnames(beta_mat1) <- beta_names_m1
colnames(beta_mat2) <- beta_names_m2

# Model 1
as.data.frame(beta_mat1) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = beta_names_m1)) |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.5, linewidth = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Posterior Distributions — Model 1",
    subtitle = "age + sex + cigs_per_day + prevalent_hyp + pulse_pressure",
    x        = "β value", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/beta_posterior_model1.png", width = 10, height = 6)

# Model 2
as.data.frame(beta_mat2) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = beta_names_m2)) |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.5, linewidth = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Posterior Distributions — Model 2",
    subtitle = "age + sex + cigs_per_day + prevalent_hyp + pulse_pressure + glucose",
    x        = "β value", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(color = "grey40", size = 10)
  )

ggsave("Figures/beta_posterior_model2.png", width = 10, height = 6)



# Posterior Predictive on Test Set -------------------------------------

beta_mat2 <- do.call(rbind, lapply(samps2, function(ch) ch[, grep("^beta", colnames(ch))]))

X_m2_test <- X_test[, c("age", "sex", "cigs_per_day", "prevalent_hyp", "pulse_pressure", "glucose")]

# For each posterior draw, compute pi for every test observation
# beta_mat2: (n_samples x 7), X_m2_test: (N_test x 6)
# Design matrix with intercept
X_test_design <- cbind(1, X_m2_test)  # N_test x 7

# Matrix of predicted probabilities: n_samples x N_test
logit_mat  <- beta_mat2 %*% t(X_test_design)
pi_mat     <- 1 / (1 + exp(-logit_mat))   # n_samples x N_test

dim(pi_mat) # samples x N_test

pi_mean    <- colMeans(pi_mat)


# Results Reporting — Model 2 ------------------------------------

roc_obj    <- roc(Y_test, pi_mean)
best_thresh <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")

Y_hat_opt <- ifelse(pi_mean >= best_thresh$threshold, 1, 0)

cm_opt <- confusionMatrix(
  factor(Y_hat_opt, levels = c(0, 1)),
  factor(Y_test,    levels = c(0, 1)),
  positive = "1"
)
print(cm_opt)

# Summary Metrics Table
metrics_df <- data.frame(
  Metric = c("AUC", "Optimal Threshold", "Accuracy", "Balanced Accuracy",
             "Sensitivity", "Specificity",
             "Pos Pred Value", "Neg Pred Value", "Kappa"),
  Value = round(c(
    auc(roc_obj),
    best_thresh$threshold,
    cm_opt$overall["Accuracy"],
    cm_opt$byClass["Balanced Accuracy"],
    cm_opt$byClass["Sensitivity"],
    cm_opt$byClass["Specificity"],
    cm_opt$byClass["Pos Pred Value"],
    cm_opt$byClass["Neg Pred Value"],
    cm_opt$overall["Kappa"]
  ), 4)
)
print(metrics_df)

# ROC Curve
ggroc(roc_obj, color = "steelblue", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.3, y = 0.2,
           label = paste0("AUC = ", round(auc(roc_obj), 4)),
           size = 4.5, color = "steelblue") +
  labs(
    title    = "ROC Curve — Model 2",
    subtitle = paste0("Optimal threshold (Youden's J) = ", round(best_thresh$threshold, 4)),
    x        = "Specificity",
    y        = "Sensitivity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    panel.grid.minor = element_blank()
  )
ggsave("Figures/roc_model2.png", width = 7, height = 5)

# Confusion Matrix Heatmap
cm_df <- as.data.frame(cm_opt$table)
colnames(cm_df) <- c("Predicted", "Actual", "Freq")

ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Freq), size = 6, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title    = "Confusion Matrix — Model 2",
    subtitle = paste0("Threshold = ", round(best_thresh$threshold, 4),
                      " (Youden's J)"),
    x        = "Actual",
    y        = "Predicted",
    fill     = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    panel.grid       = element_blank(),
    legend.position  = "right"
  )
ggsave("Figures/confusion_matrix_model2.png", width = 6, height = 5)


# Predicted Probability Distribution by True Class
data.frame(pi_mean = pi_mean, Y_test = factor(Y_test, labels = c("No CHD", "CHD"))) |>
  ggplot(aes(x = pi_mean, fill = Y_test)) +
  geom_density(alpha = 0.5, linewidth = 0.7) +
  geom_vline(xintercept = best_thresh$threshold,
             linetype = "dashed", color = "grey30", linewidth = 0.6) +
  scale_fill_manual(values = c("No CHD" = "steelblue", "CHD" = "tomato")) +
  annotate("text", x = best_thresh$threshold + 0.02, y = 8,
           label = paste0("t = ", round(best_thresh$threshold, 4)),
           size = 3.5, color = "grey30", hjust = 0) +
  labs(
    title    = "Posterior Predictive Probabilities by True Class — Model 2",
    subtitle = "Dashed line at optimal threshold",
    x        = "Posterior Mean P(Y = 1)",
    y        = "Density",
    fill     = "True Class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10)
  )
ggsave("Figures/pred_prob_by_class_model2.png", width = 8, height = 5)
