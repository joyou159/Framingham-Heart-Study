library(tidyverse)      
library(skimr)          
library(DataExplorer)   
library(naniar)         
library(VIM)            
library(corrplot)       
library(GGally)         
library(ggmosaic)       
library(patchwork)      
library(mice)           
library(janitor)        


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

