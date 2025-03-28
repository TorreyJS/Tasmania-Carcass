################################################################################
# Bounded variable analysis -- use for variables like percent consumption, etc
# glmm models
# Both within AND between site analyses (summer 2023)
################################################################################
### Load libraries ---
library(readxl);library(forcats);library(glmmTMB); library(emmeans)
library(DHARMa); library(plotrix);library(dplyr); library(tidyr)
library(ggplot2); library(car);library(ggsci); library(nlme)

### Data prep ---- 
# this code reads in data, renames "rep" to "block", creates "plot" and "site_rep"
# columns, makes day a factor, changes 100% consumed values to 99.9%, creates a 
# consumed proportion column (decimal instead of percent), and adds a time column (integer)

mydata <- read.xlsx("Data/2023_carcass persistence pct.xlsx") |> janitor::clean_names() |> 
  mutate(block = as.character(rep),
         plot = paste(site, trt, rep, sep = "_"),
         site_rep = paste(site, block, sep = "_"),
         day_fct = as.factor(days)) |> 
  mutate(consumed_pct = case_when(
    overall == 100 ~ 99.9,
    .default = overall)) |> 
  mutate(consumed_prop = consumed_pct/100) |> 
  mutate(time = as.integer(as.factor(days)))

# remove Takone E1, because cage was breached
mydata <- mydata %>%
  filter(plot != "T_E_1")

### filter each site into separate objects

## First, analyzing each site separately
data_H <- mydata |>  filter(site == "A")
data_N <- mydata |>  filter(site == "B")
data_L <- mydata |>  filter(site == "G")
data_M <- mydata |>  filter(site == "T")

### Main analysis ----

# single site - repeat for each site

# full model 
m_H1 <- glmmTMB(consumed_prop ~ trt*day_fct + (1|block) +
                  cs(day_fct + 0 | plot),
                family = beta_family(), 
                data = data_H) # throws a warning message

plot(simulateResiduals(m_a1))

# reduced model
m_H2 <- glmmTMB(consumed_prop ~ trt*day_fct + (1|block),
                family = beta_family(), 
                data = data_H)

plot(simulateResiduals(m_H2))

car::Anova(m_H2, type = 3)

(H_emm <- emmeans(m_H2, ~ trt | day_fct, type = "response"))
pairs(H_emm)


## combined sites:

m_combo <- glmmTMB(consumed_prop ~ site*trt*day_fct + (1|block),
                   family = beta_family(),
                   disp = ~ day_fct,
                   data = mydata)

plot(simulateResiduals(m_combo))

car::Anova(m_combo, type = 3) 

(combo_emm <- emmeans(m_combo, ~ site, type = "response"))
pairs(combo_emm)


################################################################################
### Consumption: Days of persistence
pers <- read.xlsx("data/day removed both yrs.xlsx")

pers_summ <- pers %>% filter(season == "summer" & trtctrl == "T") |>
  mutate(filterID = paste0(site, treatment, rep))

# remove TE1
pers_summ <- pers_summ %>%
  filter(filterID != "TE1")

pers_summ <- pers_summ %>%
  mutate(Site = case_when(site == "B" ~ "No Devils",
                          site == "A" ~ "High Devils",
                          site == "G" ~ "Low Devils",
                          site == "T" ~ "Moderate Devils"),
         Access = case_when(treatment == "E" ~ "Excluded",
                            treatment == "S" ~ "Open"))
pers_summ$Site <- as.factor(pers_summ$Site)

### can do stats on the pers_summ dataframe now

## STATS: lme with replicate as random
persistence <- lme(day_remove ~ Site*Access,
                   random = ~ 1 | rep,
                   data = pers_summ, na.action = na.exclude)

plot(persistence)
qqnorm(persistence) 

# anova
car:: Anova(persistence, type = "marginal") 

# Calculate emmeans
eres <- emmeans(persistence, ~ Site)
pairs(eres)  
