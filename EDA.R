library(tidyverse)      
library(skimr)          
library(DataExplorer)   
library(naniar)         
library(VIM)            
library(corrplot)       
library(GGally)         
library(ggmosaic)       
library(patchwork)      
library(janitor)        
library(modelsummary)
library(moments)
library(ggcorrplot)


# Data loading and investigation ---------------------------------------


df <- read.csv("framingham.csv") %>% clean_names()  # janitor: lowercase snake_case

dim(df)
head(df)
str(df)
glimpse(df)

map(df, function(x) length(unique(x)))

# variables names
binary_vars <- c("sex", "current_smoker", "bp_meds",
                 "prevalent_stroke", "prevalent_hyp", "diabetes")

ordinal_vars <- c("education")

continuous_vars <- c("age", "cigs_per_day", "tot_chol", "sys_bp",
                     "dia_bp", "bmi", "heart_rate", "glucose")

target_var <- "ten_year_chd"



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


# for numeric variables
datasummary_skim(
  framingham,
  type = "numeric",
  fun_numeric = list(
    Mean   = \(x) round(mean(x, na.rm = TRUE), 2),
    SD     = \(x) round(sd(x, na.rm = TRUE), 2),
    Median = \(x) round(median(x, na.rm = TRUE), 2),
    IQR    = \(x) round(IQR(x, na.rm = TRUE), 2),
    Min    = \(x) round(min(x, na.rm = TRUE), 2),
    Max    = \(x) round(max(x, na.rm = TRUE), 2),
    P5     = \(x) round(quantile(x, 0.05, na.rm = TRUE), 2),
    P95    = \(x) round(quantile(x, 0.95, na.rm = TRUE), 2),
    Skewness = \(x) round(skewness(x, na.rm = TRUE), 2),
    N_Miss = \(x) as.integer(sum(is.na(x)))
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
  ) |>
  gt::opt_stylize(style = 6, color = "blue") |>
  gt::gtsave("Figures/table_categorical.png")


# Univariate investigation ---------------------------------------

# target variable 
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


# histograms with kernel density estimation
framingham |>
  select(all_of(continuous_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  drop_na(value) |>        # drop missing values (manually)                           
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


# boxplots
framingham |>
  select(all_of(continuous_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  drop_na(value) |>                                   # drop missing values
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
  geom_jitter(width = 0.2, alpha = 0.07, color = "steelblue", size = 0.5) +
  geom_boxplot(alpha = 0.5, fill = "lightyellow",
               outlier.shape = NA, linewidth = 0.7) +
  stat_summary(fun = mean, geom = "point",
               shape = 18, size = 3, color = "darkred") +
  facet_wrap(~variable, scales = "free_y", ncol = 4) +
  labs(title    = "Boxplots of Continuous Predictors",
       subtitle = "Red diamond = mean; jittered points show raw data",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11)

ggsave("Figures/Boxplots_of_continuous_vars.png")



# comparing skewness
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
       subtitle = "Dashed lines at +1 or -1",
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
  drop_na(all_of(continuous_vars)) |> 
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

ggsave("Figures/continous_predictors_by_CHD_point_biserial_corr.png")


# categorical variables
framingham |>
  select(all_of(binary_vars), education, ten_year_chd) |>
  drop_na(all_of(binary_vars)) |> 
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




# correlation matrix
framingham_cor <- framingham |>
  mutate(
    sex            = as.integer(sex),
    current_smoker = as.integer(current_smoker),
    bp_meds        = as.integer(bp_meds),
    prevalent_stroke = as.integer(prevalent_stroke),
    prevalent_hyp  = as.integer(prevalent_hyp),
    diabetes       = as.integer(diabetes),
    education      = as.integer(education),
    ten_year_chd   = as.integer(ten_year_chd)
  )

cor_matrix <- cor(framingham_cor, use = "pairwise.complete.obs")

ggcorrplot(
  cor_matrix,
  type       = "upper",
  lab        = TRUE,
  lab_size   = 3,
  method     = "square",
  colors     = c("tomato", "white", "steelblue"),
  outline.color = "white",
  tl.cex     = 10,
  tl.srt     = 45
) +
  labs(title    = "Correlation Matrix",
       subtitle = "Blue = positive | Red = negative correlation") +
  theme(plot.title    = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        legend.title  = element_text(size = 10))

ggsave("Figures/correlation_matrix.png", width = 12, height = 10, dpi = 150)



# Missing Data ------------------------------------------------------------


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


# missingness pattern 
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
png("Figures/co_missingness.png", width = 2400, height = 1600, res = 300)
print(co_missingness)  # explicit print() to be safe
dev.off()


# dropping co-missingness (observations has missingness > 1 missing value)

glimpse(framingham)


# First check how many rows have more than 1 missing
framingham |>
  mutate(n_missing = rowSums(is.na(across(everything())))) |>
  count(n_missing)


framingham |>
  mutate(n_missing = rowSums(is.na(across(everything())))) |>
  filter(n_missing > 1 | is.na(heart_rate)) |>
  count(ten_year_chd) |>
  mutate(prop = n / sum(n))



framingham_clean <- framingham |>
  filter(!is.na(heart_rate)) |>
  mutate(n_missing = rowSums(is.na(across(everything())))) |>
  filter(n_missing <= 1) |>
  select(-n_missing)

glimpse(framingham_clean)


# preprocessing (feature engineering + outlier treatment) ----------------------------------------

# Feature engineering

framingham_clean <- framingham_clean |>
  mutate(pulse_pressure = sys_bp - dia_bp) |>
  select(-sys_bp, -dia_bp)

glimpse(framingham_clean)


framingham_clean |>
  summarise(
    min = min(pulse_pressure),
    max = max(pulse_pressure),
    mean = mean(pulse_pressure),
    sd = sd(pulse_pressure)
  )



# outlier detection
continuous_vars <- c("age", "cigs_per_day", "tot_chol", "pulse_pressure",
                     "bmi", "heart_rate", "glucose")

framingham_clean |>
  select(all_of(continuous_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  group_by(variable) |>
  summarise(
    P1  = quantile(value, 0.01, na.rm = TRUE),
    P99 = quantile(value, 0.99, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE)
  )


# glucose vs diabetes (leave all of them, they are clinically plausible)
framingham_clean |>
  filter(glucose > 200) |>
  select(glucose, diabetes, ten_year_chd)


dim(framingham_clean)


framingham_clean <- framingham_clean |>
  filter(is.na(tot_chol) | tot_chol <= 500)

# Confirm number of observations
nrow(framingham_clean)



saveRDS(framingham_clean, "Data/framingham_processed.rds")

