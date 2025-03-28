################################################################################
# Carcass data analysis: COMPARING SUMMER AND WINTER (end of experiment)
# input file: timefinal_bothyears.xlsx
################################################################################
library(readxl);library(nlme);library(emmeans);library(dplyr);library(tidyr); 
library(openxlsx);library(car); library(lme4);library(glmmTMB)

seasonal <- read_xlsx("Data/timefinal_bothyears.xlsx")|> 
  rename(plot = "trt1", splitplot = "trt2") |> 
  mutate(block = as.character(block),
         time_fct = as.character(time),
         year_fct = as.character(year),
         sample_site = paste(block, plot, splitplot, sep = "_"))|> 
  arrange(site, block, plot, splitplot, time) 

str(seasonal)
# this model structure is appropriate for all unbounded variables:
# columns 14, 15, 17-20, 29, 30, 40
colnames(seasonal)[30] <- "bac_shandiv"

hist(seasonal$ITS_shandiv)
# no transformation: pH, 16S shandiv, ITS shandiv
# log transformation: N, NH4, NO3, C, EC, rk

mod3 <- lme(log(ugN_gdry) ~ year*plot*splitplot, 
            random = ~ 1|site/block,
            data = seasonal, na.action = na.exclude)

plot(mod3)

Anova(mod3, type = "3", test.statistic = "F")


emm1 <- emmeans(mod3, ~ year|splitplot, type = "response")
emm1
pairs(emm1)

# 2. DAYS OF PERSISTENCE:
day <- read_xlsx("data/day removed both yrs.xlsx", sheet =1)
# site names must match
day <- day %>% mutate(site = case_when(site == "G" ~ "goulds",
                                       site == "T" ~ "takone",
                                       site == "B" ~ "bruny",
                                       site == "A" ~ "arthur",
                                       TRUE ~ site),
                      rep = factor(rep))

hist(day$day_remove) # no transformation needed
str(day)
unique(day$rep)
mod3 <- lme(day_remove ~ season*treatment, 
            random = ~ 1|site/rep,
            data = day, na.action = na.exclude)

plot(mod3)

Anova(mod3, type = "3", test.statistic = "F")


emm1 <- emmeans(mod3, ~ season, type = "response")
emm1
pairs(emm1)


# 3. PERCENT CONSUMPTION (use glmm)
hist(seasontrt$pct.consumed)

seasontrt <- seasontrt %>% 
  mutate(consumed2 = case_when(pct.consumed == 100 ~ 99,
                               pct.consumed == 0 ~ 1,
                               TRUE ~ pct.consumed),
         prop_consumed = consumed2/100)

mod2 <- glmmTMB(prop_consumed ~ plot*year + (1|site/block),
                family = beta_family(), 
                data = seasontrt) # throws a warning message

car::Anova(mod2, type = 3)

emm2 <- emmeans(mod2, ~ plot, type = "response")
emm2
pairs(emm2)

# 4. Percent BONE CONSUMPTION
sk <- read.csv("data/skeleton.csv") # this is only from day 30, summer
# value 'bone' is amount of bone CONSUMED
# I want to know about SITE and TREATMENT diffs
sk <- sk %>%
  mutate(bone2 = case_when(bone == 100 ~ 99,
                        bone == 0 ~ 1,
                        TRUE ~ bone),
         prop_consumed = bone2/100,
         rep = factor(rep))

mod1 <- glmmTMB(prop_consumed ~ site*trt + (1|rep),
                family = beta_family(),
                data = sk)

car::Anova(mod1, type = 3)

emm3 <- emmeans(mod1, ~ site, type = "response")
emm3
pairs(emm3)
# H,M>C; M>L
