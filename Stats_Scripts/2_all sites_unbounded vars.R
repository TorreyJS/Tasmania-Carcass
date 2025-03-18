################################################################################
# Carcass data analysis stage 2: BETWEEN site analysis
# input file: emmeans data from stage 1 (emmeans_within sites_unbounded.xlsx)
# This script compares across ALL SITES, so only interactions with SITE are of interest
# note that to avoid NaNs, adjust = "fdr" is necessary. Do this for ALL COMPARISONS for consistency
################################################################################
library(dplyr);library(tidyr);library(nlme);library(emmeans);library(performance)
library(ggplot2)

stage1 <- read_xlsx("outputs/emmeans_within sites_unbounded_no TE1.xlsx") |> 
  mutate(time_fct = as.factor(time_fct)) |> 
  rename(carcass_loc = plot, 
         sampling_loc = splitplot) |> 
  mutate(time_ = as.numeric(time_fct))

str(stage1)

# Repeat the process below for each variable of interest:

### example: TOTAL DISSOLVED NITROGEN: 
s1_n <- stage1 |> dplyr::filter(response_var == "ugN_gdry") # 80 observations

s1_t <- s1_n %>% filter(sampling_loc == "T")
s1_c <- s1_n %>% filter(sampling_loc == "C")

ct <- inner_join(s1_t, s1_c, by = c("time_fct","carcass_loc","response_var","site","time_"))

ct <- ct %>%
  mutate(emmean = emmean.x/emmean.y,
         SE = rowMeans(cbind(SE.x, SE.y), na.rm = TRUE))

hist(ct$emmean) # if log distribution--log transform data in model (if not, delete 'log()')

# only want interactions including site, since this is for ACROSS SITE comparisons
m1_n <- gls(log(emmean) ~ site*carcass_loc +
              site*time_fct,
            correlation = corAR1(form = ~1|time_),
            weights = ~ I(1/SE), 
            data = ct)

shapiro.test(m1_n$residuals) # 
plot(m1_n) 
qqPlot(m1_n$residuals)
hist(m1_n$residuals)

car::Anova(m1_n, type = 3) 

emmeans_results <- emmeans(m1_n, ~ site | carcass_loc , type = "response")
pairs(emmeans_results, adjust = "fdr")
