################################################################################
# Summer carcass data analysis stage 1: UNBOUNDED variables at EACH site 
# input file: data with outliers checked (none removed), negative PO4 values set to 0
# input file name: 2023Carcass_NoOuts.csv
# output file: all_emmeans and all_anova
# all_emmeans used for 4_stage 2 stats comparing across sites!
# this code corresponds to approximately lines 422-429 of MS
################################################################################
library(readxl);library(nlme);library(emmeans);library(dplyr);library(tidyr); 
library(openxlsx);library(DataExplorer) 


mydata <- read_xlsx("data/2023Carcass_NoOuts_062624.xlsx") |> 
  rename(plot = "trt1", splitplot = "trt2", shandiv_16S = "16s_shandiv") |> 
  mutate(block = as.character(block),
         time_fct = as.character(time),
         sample_site = paste(block, plot, splitplot, sep = "_"),
         po4_plus1 = as.numeric(ugPO4_gdry) + 0.001,
         filterID = paste0(site, plot, block))|> 
  arrange(site, block, plot, splitplot, time) #|> janitor::clean_names()

## remove TE1 - cage failed
mydata <- mydata %>%
  filter(filterID != "TE1")

table(mydata$sample_site, mydata$time_fct, mydata$site)

# split plot design
# main plot = "plot" cage enclosure/stake ("E"/"S")
# split plot = "splitplot" control/treatment ("C"/"T"), sample location

plot_histogram(mydata)
plot_boxplot(mydata, by = 'site')

table(mydata$site, mydata$year)  # 120 for all sites but Bruny (50)


### filter each site into separate objects

## analyzing each site separately
data3H <- mydata %>%  filter(site == "A") # Salmon River (High Devils)
data3N <- mydata %>%  filter(site == "B") # Bruny Island (No Devils)
data3L <- mydata %>%  filter(site == "G") # Blue Tier (Low Devils)
data3M <- mydata %>%  filter(site == "T") # West Takone (Moderate Devils)


### Analysis ---- 

# repeated measures structure over time correlation handled with corAR1:
corr_str1 = corAR1(form = ~ time|block/plot/sample_site, value = 0.2, fixed = FALSE)

#### First, define function to fit a log transformed OR untransformed model to each variable
lm_func <- function(yvar, log_mod = 1, df) {
  # Define formula based on whether or not log is used (1 or 0)
  if (log_mod) {
    f <- as.formula(paste0("log(", yvar, ") ~ time_fct*plot*splitplot"))
    print(paste("Data log transformed for", yvar))
  } else {
    f <- as.formula(paste0(yvar, " ~ time_fct*plot*splitplot"))
  }
  
  # Fit the model using lme
  m0 <- lme(f,
            random = ~ 1 | block/plot/sample_site,
            corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
            data = df, na.action = na.exclude)
  
  # Print the fitted model summary
  print(summary(m0))
  
  # Print diagnostic plots to console to assess fit
  print(plot(m0, main = paste(yvar)))
  print(qqnorm(m0, main = paste(yvar)))
  
  # Perform ANOVA and calculate emmeans
  a <- anova(m0, type = "marginal")
  
  # Calculate emmeans
  e <- emmeans(m0, ~ splitplot*plot | time_fct, type = 'link')
  
  # Manually back-transform from log scale to original scale
  if (log_mod) {
    emmeans_df <- as.data.frame(summary(e))
    emmeans_df$emmean <- exp(emmeans_df$emmean)
    # Apply exponential transformation to SE using correct formula
    emmeans_df$SE <- emmeans_df$emmean * sqrt(exp(emmeans_df$SE^2) - 1)
    emmeans_df$lower.CL <- exp(emmeans_df$lower.CL)  # Adjust confidence limits
    emmeans_df$upper.CL <- exp(emmeans_df$upper.CL)
  } else {
    emmeans_df <- as.data.frame(summary(e))
  }
  
  # Add the response variable name to the dataframes
  emmeans_df$response_var <- yvar
  anova_df <- as.data.frame(a)
  anova_df$response_var <- yvar
  
  # Round p-values to 4 decimal places in ANOVA dataframe
  anova_df$`p-value` <- round(anova_df$`p-value`, 4)
  
  return(list(model = m0, emmeans = emmeans_df, anova = anova_df))
}

#### Next, pull out ONLY unbounded variables (bounded variables analyzed separately)
# only unbounded variables taken over all time points fit the first model structure:
unbound <- colnames(mydata)[c(13,14,16:19,28,29,39,49)]
# determine whether each variable is log or not (by looking at histogram)
# then print the decision (1 = log, 0 = no transform) to a list:
vars_log_df <- data.frame(
  variable = unbound,
  log_mod = c(1,1,1,1,1,1,1,0,0,1), # change these for each SITE
  stringsAsFactors = FALSE)

# Initialize empty dataframes to store results
all_emmeans <- data.frame()
all_anova <- data.frame()

# Loop through each response variable and apply the function (do this 4x)
for (i in 1:nrow(vars_log_df)) {
  response_var <- vars_log_df$variable[i]
  log_mod <- vars_log_df$log_mod[i]
  
  result <- lm_func(response_var, log_mod = log_mod, df = data3a)
  
  # Combine the results into a single dataframe
  all_emmeans <- rbind(all_emmeans, result$emmeans)
  all_anova <- rbind(all_anova, result$anova)
}

# After each loop, save EACH output to the correct file:
high_emmeans <- all_emmeans
high_emmeans$site <- "High Devils"
high_anova <- all_anova
high_anova$site <- "High Devils"

no_emmeans <- all_emmeans
no_emmeans$site <- "No Devils"
no_anova <- all_anova
no_anova$site <- "No Devils"

low_emmeans <- all_emmeans
low_emmeans$site <- "Low Devils"
low_anova <- all_anova
low_anova$site <- "Low Devils"

mod_emmeans <- all_emmeans
mod_emmeans$site <- "Moderate Devils"
mod_anova <- all_anova
mod_anova$site <- "Moderate Devils"


#### Once all emmeans results are saved, add row names  
no_anova <- data.frame(row_name = rownames(no_anova), no_anova, row.names = NULL)
low_anova <- data.frame(row_name = rownames(low_anova), low_anova, row.names = NULL)
mod_anova <- data.frame(row_name = rownames(mod_anova), mod_anova, row.names = NULL)
high_anova <- data.frame(row_name = rownames(high_anova), high_anova, row.names = NULL)

## then combine:
allsite_anova <- rbind(no_anova, low_anova, mod_anova, high_anova)
sig_anova <- allsite_anova[allsite_anova$p.value<0.05,] # filter just significant effects

allsite_emmeans <- rbind(no_emmeans, low_emmeans, mod_emmeans, high_emmeans)

#### Save the output files:
write.xlsx(allsite_anova, "outputs/anova_within sites_unbounded_no TE1.xlsx")
write.xlsx(sig_anova, "outputs/anova_within sites_sig only_no TE1.xlsx")
write.xlsx(allsite_emmeans, "outputs/emmeans_within sites_unbounded_no TE1.xlsx")




######### SINGLE VARIABLE ANALYSIS -- for emmeans comparisons OR determining if log/not #######
# repeated measures structure over time correlation handled with corAR1:
corr_str1 = corAR1(form = ~ time|block/plot/sample_site, value = 0.2, fixed = FALSE)

hist(data3H$ugC_gdry) # 

# Fit the model using lme
m0 <- lme((ugC_gdry) ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = data3H, na.action = na.exclude)

m1 <- lme(log(ugC_gdry) ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = data3H, na.action = na.exclude)

AIC(m0, m1)

# Print diagnostic plots to console to assess fit of chosen model
qqnorm(m1) # 
shapiro.test(m1$residuals) 

# Calculate emmeans
eres <- emmeans(m1, ~ splitplot | time_fct, type = 'link') # can change main effects 
pairs(eres)  

