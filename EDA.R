library(tidyverse)      # data wrangling + ggplot2
library(skimr)          # rich summary statistics
library(DataExplorer)   # automated EDA
library(naniar)         # missing data analysis
library(VIM)            # missing data visualization
library(corrplot)       # correlation matrices
library(GGally)         # pair-plots
library(ggmosaic)       # mosaic plots for categorical vars
library(patchwork)      # combining plots
library(mice)           # missing data patterns (before imputation)
library(janitor)        # clean column names
library(e1071)
library(kableExtra)
library(gt)
library(modelsummary)
library(patchwork)


install.packages("")

# Dataset loading and investigation ---------------------------------------


df <- read.csv("framingham.csv") %>% clean_names()  # janitor: lowercase snake_case

dim(df)
head(df)
str(df)
glimpse(df)

map(df, function(x) length(unique(x)))

binary_vars <- c("sex", "current_smoker", "bp_meds",
                     "prevalent_stroke", "prevalent_hyp", "diabetes")

ordinal_vars <- c("education")

continuous_vars <- c("age", "cigs_per_day", "tot_chol", "sys_bp",
                     "dia_bp", "bmi", "heart_rate", "glucose")

target_var <- "ten_year_chd"


# Convert convert discrete vars into factors


framingham <- df |>
  mutate(
    # binary -> factor
    male = factor(male, labels = c("Female", "Male")),
    current_smoker   = factor(current_smoker, labels = c("No", "Yes")),
    bp_meds = factor(bp_meds, labels = c("No", "Yes")),
    prevalent_stroke = factor(prevalent_stroke, labels = c("No", "Yes")),
    prevalent_hyp = factor(prevalent_hyp, labels = c("No", "Yes")),
    diabetes = factor(diabetes, labels = c("No", "Yes")),
    
    # ordinal -> ordered factor
    education = factor(education, 
                       levels = 1:4,
                       labels  = c("Some HS", "HS Grad", "Some College", "College+"),
                       ordered = TRUE),
    
    ten_year_chd = factor(ten_year_chd, labels = c("No", "Yes"))
  ) |> 
  rename(sex = male)


glimpse(framingham)


# Summary Statistics ------------------------------------------------------

skim(framingham) 

datasummary_skim(
  framingham,
  type      = "numeric",
  fun_numeric = list(
    N        = \(x) sum(!is.na(x)),
    Mean     = \(x) round(mean(x,na.rm = TRUE), 2),
    SD       = \(x) round(sd(x,na.rm = TRUE), 2),
    Median   = \(x) round(median(x,na.rm = TRUE), 2),
    IQR      = \(x) round(IQR(x,na.rm = TRUE), 2),
    P5       = \(x) round(quantile(x, 0.05, na.rm = TRUE), 2),
    P95      = \(x) round(quantile(x, 0.95, na.rm = TRUE), 2),
    Skewness = \(x) round(moments::skewness(x, na.rm = TRUE), 2),
    N_Miss   = \(x) as.integer(sum(is.na(x)))
  ),
  output = "Figures/numeric_summary.png"
)

cat_summary <- framingham |>
  select(where(is.factor)) |>
  mutate(across(everything(), as.character)) |> 
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") |>
  group_by(Variable) |>
  summarise(
    N        = sum(!is.na(Value)),
    N_Miss   = sum(is.na(Value)),
    N_Levels = n_distinct(Value, na.rm = TRUE),
    Mode     = names(sort(table(Value), decreasing = TRUE))[1],
    Mode_N   = max(table(Value)),
    `Mode %` = round(max(table(Value)) / sum(!is.na(Value)) * 100, 1)
  )

cat_summary |>
  gt::gt() |>
  gt::tab_header(
    title    = "Categorical Variables — Descriptive Summary",
    subtitle = "Framingham Heart Study"
  ) |>
  gt::opt_stylize(style = 6, color = "blue") |>
  gt::gtsave("Figures/table_categorical.png")

# variables Analysis ---------------------------------------

# target variable analysis

chd_summary <- framingham |>
  count(ten_year_chd) |>
  mutate(pct = n / sum(n) * 100)


p_chd_bar <- ggplot(chd_summary,
                    aes(x = ten_year_chd, y = n, fill = ten_year_chd)) +
  geom_col(width = 0.5, alpha = 0.9) + # to represent the value in each cell
  geom_text(aes(label = paste0("n = ", n, "\n(", round(pct, 1), "%)")),
            vjust = -0.3, fontface="bold", size=4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "10-Year CHD Risk: Class Distribution",
       subtitle = paste0("Imbalance ratio ~ 1 : ", round(chd_summary$n[1] / chd_summary$n[2], 1)),
       x = "10-Year CHD", y = "Count") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

p_chd_bar

ggsave("Figures/target_variable_overview.png", plot = p_chd_bar)


# Continuous Variables 

framingham |>
  select(all_of(continuous_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  drop_na(value) |>                                   
  mutate(variable = factor(variable,
                           levels = continuous_vars,
                           labels = c("Age (years)", "Cigarettes/Day",
                                      "Total Cholesterol", "Systolic BP",
                                      "Diastolic BP", "BMI",
                                      "Heart Rate", "Glucose"))) |>
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 35, fill = "steelblue",
                 color = "white", alpha = 0.75) +
  geom_density(color = "darkred", linewidth = 0.9) +
  facet_wrap(~variable, scales = "free", ncol = 4) +
  labs(title    = "Distributions of Continuous Predictors",
       subtitle = "Histogram with kernel density overlay",
       x = NULL, y = "Density") +
  theme_minimal(base_size = 11)


ggsave("Figures/Distributions_of_Continuous_var.png")


framingham |>
  select(all_of(continuous_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  drop_na(value) |
  mutate(
    variable = factor(
      variable,
      levels = continuous_vars,
      labels = c("Age", "Cigarettes/Day", "Total Chol.",
                 "Systolic BP", "Diastolic BP",
                 "BMI", "Heart Rate", "Glucose")
    )
  ) |>
  ggplot(aes(x = "", y = value)) +
  geom_jitter(width = 0.2, alpha = 0.07, color = "steelblue", size = 0.5, na.rm = TRUE) +
  geom_boxplot(alpha = 0.5, fill = "lightyellow",
               outlier.shape = NA, linewidth = 0.7, na.rm = TRUE) +
  stat_summary(fun = \(x) mean(x, na.rm = TRUE),
               geom = "point", shape = 18, size = 3, color = "darkred", na.rm = TRUE) +
  facet_wrap(~variable, scales = "free_y", ncol = 4) +
  labs(title    = "Boxplots of Continuous Predictors",
       subtitle = "Red diamond = mean; jittered points show raw data",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11)

ggsave("Figures/Boxplots_of_continuous_vars.png")


skew_df <- framingham |>
  select(all_of(continuous_vars)) |>
  summarise(across(everything(), ~skewness(., na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "skewness") |>
  mutate(
    direction = ifelse(skewness > 0, "Right-skewed", "Left-skewed"),
    variable  = fct_reorder(variable, skewness)
  )

ggplot(skew_df, aes(x = variable, y = skewness, fill = direction)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed",
             color = "gray40", linewidth = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(title    = "Skewness of Continuous Predictors",
       subtitle = "Dashed lines at ±1 — beyond this, transformation may help",
       x = NULL, y = "Skewness", fill = NULL) +
  theme_minimal(base_size = 13)

ggsave("Figures/skewness_plot.png")

# All binary/categorical vars in one panel

cat_plot_data <- framingham |>
  select(all_of(binary_vars), education) |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  drop_na(value) |>
  count(variable, value) |>
  group_by(variable) |>
  mutate(
    pct = n / sum(n) * 100,
    variable = recode(variable,
                      "male"             = "Sex",
                      "current_smoker"   = "Current Smoker",
                      "bp_meds"          = "BP Medication",
                      "prevalent_stroke" = "Prevalent Stroke",
                      "prevalent_hyp"    = "Hypertension",
                      "diabetes"         = "Diabetes",
                      "education"        = "Education Level"
    )
  )

ggplot(cat_plot_data,
       aes(x = value, y = pct, fill = value)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
            vjust = -0.4, size = 3.2) +
  facet_wrap(~variable, scales = "free_x", ncol = 4) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand  = expansion(mult = c(0, 0.15))) +
  labs(title    = "Categorical Predictor Distributions",
       x = NULL, y = "Percentage") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x     = element_text(angle = 30, hjust = 1))

ggsave("Figures/Categorical_distribution.png")


# Bivariate investigation -------------------------------------------------

# continuous variables
framingham |>
  select(all_of(continuous_vars), ten_year_chd) |>
  pivot_longer(-ten_year_chd, names_to = "variable", values_to = "value") |>
  mutate(variable = factor(variable,
                           levels = continuous_vars,
                           labels = c("Age", "Cigarettes/Day", "Total Chol.",
                                      "Systolic BP", "Diastolic BP",
                                      "BMI", "Heart Rate", "Glucose"))) |>
  ggplot(aes(x = value, fill = ten_year_chd, color = ten_year_chd)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values  = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  facet_wrap(~variable, scales = "free", ncol = 4) +
  labs(title    = "Continuous Predictors by CHD Status",
       x = NULL, y = "Density",
       fill = "10-yr CHD", color = "10-yr CHD") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("Figures/continous_predictors_by_CHD.png")


# Age and systolic BP show strongest association
framingham |>
  mutate(chd_num = as.integer(ten_year_chd) - 1) |>
  select(all_of(continuous_vars), chd_num) |>
  summarise(across(all_of(continuous_vars),
                   ~cor(., chd_num, use = "pairwise.complete.obs"))) |>
  pivot_longer(everything(),
               names_to  = "variable",
               values_to = "r") |>
  mutate(
    variable = fct_reorder(variable, abs(r)),
    dir      = ifelse(r > 0, "Positive", "Negative")
  ) |>
  ggplot(aes(x = variable, y = r, fill = dir)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(title    = "Point-Biserial Correlation with 10-Year CHD",
       x = NULL, y = "Pearson r with CHD", fill = "Direction") +
  theme_minimal(base_size = 13)

ggsave("Figures/continous_predictors_by_CHD_pearson_r.png")



# categorical variables
framingham |>
  select(all_of(binary_vars), education, ten_year_chd) |>
  mutate(across(c(all_of(binary_vars), education), as.character)) |>
  pivot_longer(-ten_year_chd) |>
  drop_na(value) |>                                                   
  mutate(name = recode(name,
                       "male"             = "Sex",
                       "current_smoker"   = "Current Smoker",
                       "bp_meds"          = "BP Medication",
                       "prevalent_stroke" = "Prevalent Stroke",
                       "prevalent_hyp"    = "Hypertension",
                       "diabetes"         = "Diabetes",
                       "education"        = "Education Level"
  )) |>
  ggplot(aes(x = value, fill = ten_year_chd)) +
  geom_bar(position = "fill", alpha = 0.9) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  facet_wrap(~name, scales = "free_x", ncol = 4) +
  labs(title = "CHD Proportion by Categorical Predictors",
       x = NULL, y = "Proportion", fill = "10-yr CHD") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

ggsave("Figures/categorical_predictors_by_CHD_proportion.png")






# Missing Data ------------------------------------------------------------

# Missingness bar chart

glucose_miss <- miss_var_summary(framingham) |>
  filter(variable == "glucose") |>
  mutate(pct_miss = as.numeric(pct_miss)) |>
  pull(pct_miss) |>
  round(1)

p_miss_bar <- miss_var_summary(framingham) |>
  filter(n_miss > 0) |>
  mutate(
    pct_miss = as.numeric(pct_miss),
    variable = fct_reorder(variable, pct_miss)
  ) |>
  ggplot(aes(x = variable, y = pct_miss)) +
  geom_col(fill = "tomato", alpha = 0.85) +
  geom_text(aes(label = paste0(round(pct_miss, 1), "%\n(n=", n_miss, ")")),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 12),
                     labels = scales::percent_format(scale = 1)) +
  labs(title    = "Missing Data by Variable",
       subtitle = paste0("Glucose accounts for the majority of missingness (",
                         glucose_miss, "%)"),
       x = NULL, y = "% Missing") +
  theme_minimal(base_size = 13)


p_miss_bar

ggsave("Figures/missing_proportions.png", plot = p_miss_bar)



# missingness heatmap
p_miss_heat <- vis_miss(framingham, sort_miss = TRUE) +
  labs(title = "Missingness Pattern Across All Rows") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_miss_heat

ggsave("Figures/missing_pattern.png", plot = p_miss_heat)


# co-occurrence of missingness (Are the same rows missing multiple variables?)
co_missingness <- gg_miss_upset(framingham,
              nsets = 6,
              nintersects = 15) 
co_missingness

# Is missingness in glucose related to CHD?
chd_relation_to_glocuse_missingness <-  framingham |>
  mutate(glucose_missing = is.na(glucose)) |>
  group_by(glucose_missing, ten_year_chd) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(glucose_missing) |>
  mutate(pct = n / sum(n) * 100) |>
  ggplot(aes(x = glucose_missing, y = pct, fill = ten_year_chd)) +
  geom_col(position = "stack", alpha = 0.9) +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_x_discrete(labels = c("FALSE" = "Glucose Observed",
                              "TRUE"  = "Glucose Missing")) +
  labs(title    = "CHD Rate: Glucose Missing vs Observed",
       subtitle = "Missing Completely at Random? (W.r.t CHD)",
       x = NULL, y = "Proportion (%)", fill = "CHD") +
  theme_minimal(base_size = 13)

chd_relation_to_glocuse_missingness
ggsave("Figures/Glucose_missing_vs_observed_chd.png", plot = chd_relation_to_glocuse_missingness)

chisq.test(table(framingham$ten_year_chd, is.na(framingham$glucose)))

summary(glm(is.na(glucose) ~ ten_year_chd, 
            data = framingham, 
            family = binomial()))


# ── continuous panel (your existing plot) ─────────────────────────────────
p_continuous <- framingham |>
  mutate(glucose_missing = factor(is.na(glucose),
                                  levels = c(FALSE, TRUE),
                                  labels = c("Observed", "Missing"))) |>
  select(glucose_missing, age, sys_bp, bmi, tot_chol,
         heart_rate, cigs_per_day, dia_bp) |>
  pivot_longer(-glucose_missing) |>
  group_by(glucose_missing, name) |>
  summarise(mean = mean(value, na.rm = TRUE),
            se   = sd(value,   na.rm = TRUE) / sqrt(n()),
            .groups = "drop") |>
  ggplot(aes(x = glucose_missing, y = mean,
             color = glucose_missing, group = 1)) +
  geom_pointrange(aes(ymin = mean - 1.96 * se,
                      ymax = mean + 1.96 * se), size = 0.7) +
  geom_line(color = "gray60", linetype = "dashed") +
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  labs(title    = "Continuous Covariate Profiles by Glucose Missingness",
       x = NULL, y = "Mean ± 95% CI") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")



p_continuous

ggsave("Figures/continous_covar_profiles_by_glucose_missingness.png", plot = p_continuous)



p_categorical <- framingham |>
  mutate(glucose_missing = factor(is.na(glucose),
                                  levels = c(FALSE, TRUE),
                                  labels = c("Observed", "Missing"))) |>
  select(glucose_missing, all_of(binary_vars), education) |>
  mutate(across(all_of(binary_vars), as.character),
         education = as.character(education)) |>
  pivot_longer(-glucose_missing, values_drop_na = TRUE) |>  # ← drop NA values
  group_by(glucose_missing, name, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(glucose_missing, name) |>
  mutate(pct = n / sum(n) * 100) |>
  ggplot(aes(x = value, y = pct, fill = glucose_missing)) +
  geom_col(position = "dodge", alpha = 0.85) +
  facet_wrap(~name, scales = "free_x", ncol = 4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Categorical Covariate Profiles by Glucose Missingness",
       x = NULL, y = "Proportion (%)",
       fill = "Glucose") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

p_categorical


ggsave("Figures/categorical_covar_profiles_by_glucose_missingness.png", plot = p_categorical)

