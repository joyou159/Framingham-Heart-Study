library(tidyverse)

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



