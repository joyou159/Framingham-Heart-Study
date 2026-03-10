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

install.packages("")

# Dataset loading and investigation ---------------------------------------


df <- read.csv("framingham.csv") %>% clean_names()  # janitor: lowercase snake_case

dim(df)
head(df)
str(df)
glimpse(df)

map(df, function(x) length(unique(x)))

binary_vars <- c("male", "current_smoker", "bp_meds",
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
  )


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

# Target Variable Analysis ---------------------------------------

# Class imbalance

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

# CHD rate by sex

p_chd_sex <- framingham |>
  group_by(male, ten_year_chd) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(male) |>
  mutate(pct = n / sum(n) * 100) |>
  filter(ten_year_chd == "Yes") |>
  ggplot(aes(x = male, y = pct, fill = male)) +
  geom_col(width = 0.5, alpha = 0.9) +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
            vjust = -0.4, fontface = "bold", size = 4) +
  scale_fill_manual(values = c("orchid3", "steelblue")) +
  scale_y_continuous(limits = c(0, 25),
                     labels = scales::percent_format(scale = 1)) +
  labs(title = "CHD Rate by Sex",
       x = NULL, y = "% with CHD") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

p_chd_sex

p_chd_bar + p_chd_sex +
  plot_annotation(title = "Target Variable Overview",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))

ggsave("Figures/target_variable_overview.png")

