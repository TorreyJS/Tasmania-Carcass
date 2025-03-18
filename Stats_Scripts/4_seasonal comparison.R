################################################################################
# Carcass data analysis: COMPARING SUMMER AND WINTER (end of experiment)
# input file: timefinal_bothyears.xlsx
################################################################################
library(readxl);library(nlme);library(emmeans);library(dplyr);library(tidyr); 
library(openxlsx);library(car); library(lme4)

seasonal <- read_xlsx("data/timefinal_bothyears.xlsx")|> 
  rename(plot = "trt1", splitplot = "trt2") |> 
  mutate(block = as.character(block),
         time_fct = as.character(time),
         year_fct = as.character(year),
         sample_site = paste(block, plot, splitplot, sep = "_"))|> 
  arrange(site, block, plot, splitplot, time) 


########################## carcasses ################################################
# 1. CARCASS CONSUMPTION 
# do not consider where splitplot = C (all NA values)
seasontrt <- check %>% filter(trt2 == "T")
seasontrt$year <- as.factor(seasontrt$year)

# Fit the models

### 1. linear model:
mod3 <- lm(pct.consumed ~ year*trt1*site, 
           data = seasontrt, na.action = na.exclude)

Anova(mod3)
emmeans(mod3, ~ site) 
emmeans(mod3, ~ year)
emmeans(mod3, ~ trt1) 
#emmeans(mod3, pairwise ~ site|year_fct)

# get values:
yumsum <- seasontrt %>% 
  group_by(year) %>%
  summarise(mean_cons = mean(pct.consumed, na.rm =T),
            SE_cons = std.error(pct.consumed))

# 2. DAYS OF PERSISTENCE:
# days of persistence didn't come through in this document (probably an issue while compiling)
day <- read_xlsx("data/day removed both yrs.xlsx", sheet =1)

day <- day %>% filter(trtctrl == "T")
str(day) # change year to factor or character so it will play nicely in model
day$year <- as.factor(day$year)
# site is a mess
day <- day %>% mutate(site = case_when(site == "G" ~ "low",
                                       site == "T" ~ "mod",
                                       site == "B" ~ "no",
                                       site == "A" ~ "high",
                                       TRUE ~ site))

###### Comparing ONLY time final, and treatment sites!

# Linear model:
mod4 <- lm(day_remove ~ year*treatment*site, 
           data = day, na.action = na.exclude)
plot(mod4)
Anova(mod4)
emmeans(mod4, ~ site)
emmeans(mod4, ~ year) 
emmeans(mod4, ~ treatment) 

# get values:
daysum <- day %>% 
  group_by(year) %>%
  summarise(mean_day = mean(day_remove, na.rm =T),
            SE_day = std.error(day_remove))

##############################nutrients#########################################

####### 2: SOIL NUTRIENTS
###### Comparing ONLY time final, and treatment sites!

#### follow this process for EACH variable:
# ammonium, nitrate, TDN, DOC, pH, EC
# Linear model:
modamm <- lm(ugNH4_gdry ~ year*plot*site, 
             data = seasontrt, na.action = na.exclude)
Anova(modamm) # site***, plot**, year***
emmeans(modamm, ~ year) # more ammonium (~3x) in summer than winter!


############################## microbes ########################################
# Just alpha diversity:

# SHANNON DIVERSITY:

colnames(seasontrt)[29] <- "bac_shandiv"
# there is a SINGLE row with NA (2022 hi E4) - need to omit this (row 48)
strt <- seasontrt[seasontrt$Number != "48",]
str(strt)
xtabs(~ year + trt1 + site, data = strt)
### the issue is there are no "E" values for High Devils 2022, since all cages 
# failed so were or changed to "O"; for analysis, these need to be changed back 
# to "E", mention in the text that the cages failed
strt <- strt %>% 
  mutate(trt1 = case_when(Number %in% c(26, 38, 43, 45) ~ "E", TRUE~trt1))
xtabs(~ year + trt1 + site, data = strt) # great, that solved that

# Linear model --- do for BOTH 16S and ITS alpha diversity
shanbac <- lm(bac_shandiv ~ year*trt1*site, # or ITS_shandiv
              data = strt, na.action = na.exclude)
Anova(shanbac) 
emmeans(shanbac, pairwise ~ year|site) 
emmeans(shanbac, ~ year) 


