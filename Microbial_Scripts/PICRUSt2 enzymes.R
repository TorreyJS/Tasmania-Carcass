################################################################################
# PICRUSt2.0 pathway results
# 2/13/25
# Ornithine and Lysine Decarboxylase and C, N, P enzymes
################################################################################
library(dplyr);library(readr);library(openxlsx);library(tidyr)
library(stringr);library(emmeans);library(plotrix);library(nlme)
library(ggplot2)


# bring in KEGG data: 
kegg <- read_tsv("Data/pred_metagenome_unstrat.tsv.gz")
colnames(kegg)[1] <-"keggfxn"

pwys <- kegg[,1]

################## ORNITHINE DECARBOXYLASE

### to get RELATIVE ABUNDANCE:
# Step 1: Calculate total reads per sample
total_reads_per_sample <- colSums(kegg[,-1])  # Exclude KEGG function column (1)

# Step 2: Extract only K01581 row
ornithine_reads <- kegg %>%
  filter(keggfxn == "K01581") %>%  # Assuming the KEGG column is named "KO"
  pivot_longer(-keggfxn, names_to = "sample", values_to = "reads")

# Step 3: Compute relative abundance
ornithine_reads <- ornithine_reads %>%
  mutate(total_reads = total_reads_per_sample[sample],
         rel_abundance = (reads / total_reads)*100)

meta$sample <- as.character(meta$number)

# match with metadata
orn1 <- left_join(ornithine_reads, meta, by="sample")

#colnames(orn1)[1] <- "sum"
# quick look: are they significantly different?
ornmod <- aov(rel_abundance ~ trt, data = orn1)
anova(ornmod)

emmeans(ornmod, pairwise ~ trt)
# ornithine DOWN under carcasses
orn1$sum <- as.numeric(orn1$sum)


############ Do the stats right: stage 1
str(orn1)
# change rep to "block"
# need plot: E/S--first part of trt
# need splitplot: T/C
# need sample_site: block_plot_splitplot
# need time_fct

orn_stat <- orn1 %>%
  mutate(block = as.factor(rep),
         plot = str_sub(trt, 1, 1),
         splitplot = case_when(trt == "E" ~ "T",
                               trt == 'EC' ~ "C",
                               trt == "S" ~ "T",
                               trt == "SC" ~ "C"),
         sample_site = paste(block, plot, splitplot, sep = "_"),
         time_fct = as.factor(time),
         rel_abundance = as.numeric(rel_abundance),
         TREAT = case_when(trt == "EC" ~ "Control",
                           trt == "E" ~ "Excluded",
                           trt == "SC" ~ "Control",
                           trt == "S" ~ "Open"))
enzO <- orn_stat[c(5:16)]

# for quantifying:
orn.agg <- orn_stat %>%
  group_by(trt) %>%
  summarise(mean = mean(rel_abundance, na.rm = T),
            SE = std.error(rel_abundance))
# percent change = -24% for E, -14.3% for S

# split into sites for stage 1:
ornA <- orn_stat %>% filter(site == "A") # medium devils
ornB <- orn_stat %>% filter(site == "B") # no devils
ornG <- orn_stat %>% filter(site == "G") # low devils
ornT <- orn_stat %>% filter(site == "T") # high devils

corr_str1 = corAR1(form = ~ time|block/plot/sample_site, value = 0.2, fixed = FALSE)

hist(ornT$rel_abundance) 

# remove any outliers:
#ornA2 <- ornA %>% filter(sample != "1014" & sample != "1200")
#ornG2 <- ornG %>% filter(sample != "1012")

# Fit the model using lme
m0 <- lme(rel_abundance ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = ornT, na.action = na.exclude)

m1 <- lme(log(rel_abundance) ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = ornT, na.action = na.exclude)

AIC(m0, m1)
# A: m0 lower for both data sets
# B: m0
# G: m0
# T: m1

# Print diagnostic plots to console to assess fit
qqnorm(m0) 

res_df <- ornG %>%
  mutate(residuals = residuals(m0)) %>%
  arrange(desc(abs(residuals)))  # Sort by highest absolute residuals

# View top residuals
head(res_df, 10) 
#A: sample 1014 is outlier, 1200 too 
#G: 1012

shapiro.test(m0$residuals) #

#anova
anova(m0, type = "marginal") 
# A/low: time:split** (p<0.01, F4,73=5.3)
# B: time:split** (p<0.01, F4,22=5.4)
# G: time:split** (p<0.01, F4,78=8.9)
# T: time:split*** (P,0.001, F4,78=6.9)

# Calculate emmeans
eres <- emmeans(m0, ~ splitplot| time_fct, type = 'link')
pairs(eres)  
# A (low): C>T at time 1, 3
# B (no): C>T at time 2
# G (high): C>T at time 2, 3, 4
# T (mod): C>T at time 1

## get emmeans to write to table
e = emmeans(m1, ~ splitplot*plot | time_fct, type = 'response') 
edf <- as.data.frame(summary(e))

hiemm <- edf
hiemm$site <- "High"
# 
modemm <- edf
modemm$site <- "Mod"
# 
lowemm <- edf
lowemm$site <- "Low"
# 
ctrlemm <- edf
ctrlemm$site <- "Ctrl"

emmeans_ornithine <- rbind(hiemm, modemm, lowemm, ctrlemm)
write.xlsx(emmeans_ornithine, "outputs/emmeans_ornithine.xlsx")

####### Stage 2: #########
# only want interactings including site, since this is for ACROSS SITE comparisons
hist(emmeans_ornithine$response) # fairly normal?
# sampling_loc = splitplot, carcass_loc = plot
emmeans_ornithine$time <- as.numeric(emmeans_ornithine$time_fct)

ornmod1 <- gls(response ~ site*splitplot + 
                 site*plot +
                 site*time_fct,
               correlation = corAR1(form = ~1|time),
               weights = ~ I(1/SE), 
               data = emmeans_ornithine)

ornmod2 <- gls(log(response) ~ site*splitplot + 
                 site*plot +
                 site*time_fct,
               correlation = corAR1(form = ~1|time),
               weights = ~ I(1/SE), 
               data = emmeans_ornithine)

AIC(ornmod1, ornmod2) # ornmod1 is better

plot(ornmod1) # 2 extrememly distinct groups...
anova(ornmod1, type = "marginal")
# site:time** (p<0.01, F12,52=3.0), site:split*** (p<0.001, F3,52 = 7.0); time** (p<0.05, F4,52=3.5), split**, site***

emmeans_results <- emmeans(ornmod1, ~ splitplot | time_fct | site)
pairs(emmeans_results, adjust = "fdr")
# no diffs between consecutive times at any site. 
# T>C at every single site
# E not diff from S
# FOR PLOT: ctrl A, exc B, open B



################## LYSINE DECARBOXYLASE ######
### to get RELATIVE ABUNDANCE:
# Step 1: Calculate total reads per sample
total_reads_per_sample <- colSums(kegg[,-1])  # Exclude KEGG function column (1)

# Step 2: Extract only K01581 row
lysinecarboxylase <- kegg %>%
  filter(keggfxn == "K01582") %>%  # Assuming the KEGG column is named "KO"
  pivot_longer(-keggfxn, names_to = "sample", values_to = "reads")

# Step 3: Compute relative abundance
lysine_reads <- lysinecarboxylase %>%
  mutate(total_reads = total_reads_per_sample[sample],
         rel_abundance = (reads / total_reads)*100)


enzL <- lysine_reads[c(2,5)]
enzO$number <- as.character(enzO$number)
carcenz <- full_join(enzO, enzL, by = c("number" = "sample"))
ce2 <- carcenz[c(2,1,13)]
colnames(ce2)[2:3] <- c("ODC_relabund","Lysine_relabund")

write.xlsx(ce2, "outputs/carcass enzymes.xlsx")

# match with metadata
lysine1 <- full_join(lys, meta, by = "number")

colnames(lysine1)[1] <- "sum"
# are they significantly different?
lysmod <- aov(sum ~ trt, data = lysine1)
anova(lysmod)

emmeans(lysmod, pairwise ~ trt)
# lysine up under carcasses--no significant differences


################################################################################
##### C, N, and P enzymes:
# bring in the pathways of interest:
interest <- read.xlsx("Data/pathways of interest.xlsx", sheet = 2)

################## acid phosphatases: 
aphos <- interest[8:11,4]

### to get RELATIVE ABUNDANCE:
# Step 1: Calculate total reads per sample
total_reads_per_sample <- colSums(kegg[,-1])  # Exclude KEGG function column (1)

# Step 2: Extract only genes of interest:
aphos_reads <- kegg %>%
  filter(keggfxn %in% aphos) %>%  
  pivot_longer(-keggfxn, names_to = "sample", values_to = "reads")  %>%
  mutate(total_reads = total_reads_per_sample[sample],
         rel_abundance = (reads / total_reads)*100)

meta <- read.xlsx("data/2023_metadata.xlsx")
meta$sample <- as.character(meta$number)

# match with metadata
ap1 <- left_join(aphos_reads, meta, by="sample")

acid <- ap1 %>% filter(keggfxn == "K01078" | keggfxn == "K09474")
alk <- ap1 %>% filter(keggfxn == "K01113" | keggfxn == "K01077")
# phoC <- ap1 %>% filter(keggfxn == "K09474")

# are they significantly different?
ornmod <- aov(rel_abundance ~ trt, data = alk)
anova(ornmod)
# all together = no diff
# just acid: p<0.05, F3,796 = 3.5
# just alkaline: p=0.88

# continue with just acid phosphatases
acid2 <- acid %>% group_by(sample, site, time, trt, rep) %>%
  summarise(sum_abund = sum(rel_abundance))

apmod <- aov(sum_abund ~ trt, data = acid2)
anova(apmod) # p<0.001, F3,396 = 7.1

emmeans(apmod, pairwise ~ trt)
# acid phosphatase lower under E than EC, lower at E than S, not diff btwn S than SC or at ctrls


############ Do the stats right: stage 1
colnames(acid2)
# change rep to "block"
# need plot: E/S--first part of trt
# need splitplot: T/C
# need sample_site: block_plot_splitplot
# need time_fct

phos_stat <- acid2 %>%
  mutate(block = as.factor(rep),
         plot = str_sub(trt, 1, 1),
         splitplot = case_when(trt == "E" ~ "T",
                               trt == 'EC' ~ "C",
                               trt == "S" ~ "T",
                               trt == "SC" ~ "C"),
         sample_site = paste(block, plot, splitplot, sep = "_"),
         time_fct = as.factor(time),
         rel_abundance = as.numeric(sum_abund),
         TREAT = case_when(trt == "EC" ~ "Control",
                           trt == "E" ~ "Excluded",
                           trt == "SC" ~ "Control",
                           trt == "S" ~ "Open"))

# for quantifying:
phos.agg <- phos_stat %>%
  group_by(trt) %>%
  summarise(mean = mean(sum_abund, na.rm = T),
            SE = std.error(sum_abund))
# percent change = -17% for E, -4.5% for S

# split into sites for stage 1:
phosA <- phos_stat %>% filter(site == "A") # medium devils
phosB <- phos_stat %>% filter(site == "B") # no devils
phosG <- phos_stat %>% filter(site == "G") # low devils
phosT <- phos_stat %>% filter(site == "T") # high devils

corr_str1 = corAR1(form = ~ time|block/plot/sample_site, value = 0.2, fixed = FALSE)

hist(phosT$rel_abundance) 

# remove any outliers:
phosA2 <- phosA %>% filter(sample != "1060")
phosG2 <- phosG %>% filter(sample != "1360")
phosT2 <- phosT %>% filter(sample != "1247")

# Fit the model using lme
m0 <- lme(rel_abundance ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = phosB, na.action = na.exclude)

m1 <- lme(log(rel_abundance) ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = phosB, na.action = na.exclude)

AIC(m0, m1)

# Print diagnostic plots to console to assess fit
qqnorm(m0) 


phosT$residuals <- residuals(m0)  # Assign residuals directly
res_df <- phosT[order(-abs(phosT$residuals)), ]  

# # View top residuals
head(res_df, 10)
#A: sample 1060 is outlier
#G: samples 1360, 1012
# T: 1247


shapiro.test(m0$residuals) #

#anova
anova(m0, type = "marginal") 

# Calculate emmeans
eres <- emmeans(m0, ~ splitplot | time_fct)
pairs(eres)  

## get emmeans to write to table
e = emmeans(m0, ~ splitplot*plot | time_fct, type = 'response') 
edf <- as.data.frame(summary(e))

#hiemm <- edf
#hiemm$site <- "High"
# 
#modemm <- edf
#modemm$site <- "Mod"
# 
#lowemm <- edf
#lowemm$site <- "Low"
# 
#ctrlemm <- edf
#ctrlemm$site <- "Ctrl"

emmeans_phos <- rbind(hiemm, modemm, lowemm, ctrlemm)
#write.xlsx(emmeans_phos, "outputs/emmeans_acid phosphatase.xlsx")

####### Stage 2: #########
# only want interactings including site, since this is for ACROSS SITE comparisons
hist(emmeans_phos$emmean) # some left skew
# sampling_loc = splitplot, carcsas_loc = plot
emmeans_phos$time <- as.numeric(emmeans_phos$time_fct)

phosmod1 <- gls(emmean ~ site*splitplot + 
                  site*plot +
                  site*time_fct,
                correlation = corAR1(form = ~1|time),
                weights = ~ I(1/SE), 
                data = emmeans_phos)

phosmod2 <- gls(log(emmean) ~ site*splitplot + 
                  site*plot +
                  site*time_fct,
                correlation = corAR1(form = ~1|time),
                weights = ~ I(1/SE), 
                data = emmeans_phos)

AIC(phosmod1, phosmod2) # phosmod1 is better

plot(phosmod1) # 
anova(phosmod1, type = "marginal")
# site:time***, F12,52=7.0; time***, plot* F1,52=4.2, splitplot* F1,52=5.5, site*

emmeans_results <- emmeans(phosmod1, ~ splitplot)
pairs(emmeans_results)
# site:time: high>all others at T1, low>ctrl, mod and mod>ctrl at T2, ctrl>low at T4
# plot: S>E
# splitplot: C>T

enzP <- phos_stat[,c(1,7:13)]
######### Figures:
phos.plot <- phos_stat %>%
  group_by(time, TREAT) %>%
  summarise(relabund = sum(sum_abund),
            SE = std.error(sum_abund)) %>%
  mutate(Day = case_when(time == 0 ~ 0,
                         time == 1 ~ 5,
                         time == 2 ~ 10,
                         time == 3 ~ 15,
                         time == 4 ~ 30))
library(ggsci)
npg_colors <- (pal_npg("nrc", alpha =1)(10))
npg_colors2 <- c(npg_colors, "black")

ggplot(data = phos.plot, aes(x=Day, y=relabund, color = TREAT, fill = TREAT, group = TREAT, shape = TREAT))+
  geom_errorbar(aes(ymax=relabund+SE, ymin=relabund-SE), width = 0.5)+
  geom_line()+
  geom_point(size = 3, color = "black")+
  theme_minimal()+
  ylab("Relative Abundance (%)")+
  scale_fill_manual(values = npg_colors2[c(12,1,6)], name = "Carcass Access")+
  scale_color_manual(values = npg_colors2[c(12,1,6)], name = "Carcass Access")+
  scale_shape_manual(values = c(21, 22, 24), name = "Carcass Access")+
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical")

ggsave("outputs/SI acid phosphatase fig.png", height = 4, width = 5, units = "in",
       dpi = 350, bg = "white")

################################################################################

#### NITROGEN ENZYMES #####

nit <- interest[6:7,4]

######

# Step 2: Extract only genes of interest:
nit_reads <- kegg %>%
  filter(keggfxn %in% nit) %>%  
  pivot_longer(-keggfxn, names_to = "sample", values_to = "reads") %>%
  mutate(total_reads = total_reads_per_sample[sample],
         rel_abundance = (reads / total_reads)*100) %>%
  group_by(sample) %>%
  summarise(totalabund = sum(rel_abundance)) 

nit2 <- left_join(nit_reads, meta, by = "sample")

# are they significantly different? (if yes, do real stats)
nitmod <- aov(totalabund ~ trt, data = nit2)
anova(nitmod)
# trt***

emmeans(nitmod, pairwise ~ trt)
# n enzymes lower under E than EC and S than SC, but not diff at ctrls or btwn S, E


############ Do the stats right: stage 1
colnames(nit2)
# change rep to "block"
# need plot: E/S--first part of trt
# need splitplot: T/C
# need sample_site: block_plot_splitplot
# need time_fct

nit_stat <- nit2 %>%
  mutate(block = as.factor(rep),
         plot = str_sub(trt, 1, 1),
         splitplot = case_when(trt == "E" ~ "T",
                               trt == 'EC' ~ "C",
                               trt == "S" ~ "T",
                               trt == "SC" ~ "C"),
         sample_site = paste(block, plot, splitplot, sep = "_"),
         time_fct = as.factor(time),
         rel_abundance = as.numeric(totalabund),
         TREAT = case_when(trt == "EC" ~ "Control",
                           trt == "E" ~ "Excluded",
                           trt == "SC" ~ "Control",
                           trt == "S" ~ "Open"))

enzN <- nit_stat[,c(1,13)]
# for quantifying:
nit.agg <- nit_stat %>%
  group_by(trt) %>%
  summarise(mean = mean(rel_abundance, na.rm = T),
            SE = std.error(rel_abundance))
# percent change = -9% for E, -5% for S

# split into sites for stage 1:
nitA <- nit_stat %>% filter(site == "A")
nitB <- nit_stat %>% filter(site == "B")
nitG <- nit_stat %>% filter(site == "G")
nitT <- nit_stat %>% filter(site == "T")

corr_str1 = corAR1(form = ~ time|block/plot/sample_site, value = 0.2, fixed = FALSE)

hist(nitT$rel_abundance) 
# A - pretty normal
# B - left skew
# G - normal with a high outlier
# T - left skew

# remove any outliers:
nitG2 <- nitG %>% filter(sample != "1039" & sample != "1012")
nitT2 <- nitT %>% filter(sample != "1104")

# Fit the model using lme
m0 <- lme(rel_abundance ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = nitB, na.action = na.exclude)

m1 <- lme(log(rel_abundance) ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = nitT, na.action = na.exclude)

AIC(m0, m1)
# A: m0 
# B: m0
# G: m0
# T: m0

# Print diagnostic plots to console to assess fit
qqnorm(m0) 

nitT$residuals <- residuals(m0)  # Assign residuals directly
res_df <- nitT[order(-abs(nitT$residuals)), ]  

# # View top residuals
head(res_df, 10)
# G: drop 1039 and 1012
# T: drop 1104


shapiro.test(m0$residuals) #

anova(m0, type = "marginal") 

# Calculate emmeans
eres <- emmeans(m0, ~ splitplot | time_fct, type = 'link')
pairs(eres)  

## get emmeans to write to table
e = emmeans(m0, ~ splitplot*plot | time_fct, type = 'response') 
edf <- as.data.frame(summary(e))

#hiemm <- edf
#hiemm$site <- "High"
# 
#modemm <- edf
#modemm$site <- "Mod"
# 
#lowemm <- edf
#lowemm$site <- "Low"
# 
#ctrlemm <- edf
#ctrlemm$site <- "Ctrl"

emmeans_nit <- rbind(hiemm, modemm, lowemm, ctrlemm)
#write.xlsx(emmeans_phos, "outputs/emmeans_n enzymes.xlsx")

####### Stage 2: #########
# only want interactings including site, since this is for ACROSS SITE comparisons
hist(emmeans_nit$emmean) # some left skew
# sampling_loc = splitplot, carcsas_loc = plot
emmeans_nit$time <- as.numeric(emmeans_nit$time_fct)

nitmod1 <- gls(emmean ~ site*splitplot + 
                 site*plot +
                 site*time_fct,
               correlation = corAR1(form = ~1|time),
               weights = ~ I(1/SE), 
               data = emmeans_nit)

nitmod2 <- gls(log(emmean) ~ site*splitplot + 
                 site*plot +
                 site*time_fct,
               correlation = corAR1(form = ~1|time),
               weights = ~ I(1/SE), 
               data = emmeans_nit)

AIC(nitmod1, nitmod2) # mod1 is better

plot(nitmod1) # 
anova(nitmod1, type = "marginal")
# time*** F4,52=10.9; split** F1,52=9.5

emmeans_results <- emmeans(nitmod1, ~ time_fct)
pairs(emmeans_results)
# decrease from T0-1, 1-2, not afterwards (fast change, then stable)
# split: C>T


######### Figures:
nit.plot <- nit_stat %>%
  group_by(time, TREAT) %>%
  summarise(relabund = sum(rel_abundance),
            SE = std.error(rel_abundance)) %>%
  mutate(Day = case_when(time == 0 ~ 0,
                         time == 1 ~ 5,
                         time == 2 ~ 10,
                         time == 3 ~ 15,
                         time == 4 ~ 30))

ggplot(data = nit.plot, aes(x=Day, y=relabund, color = TREAT, fill = TREAT, group = TREAT, shape = TREAT))+
  geom_errorbar(aes(ymax=relabund+SE, ymin=relabund-SE), width = 0.5, color = "black")+
  geom_line()+
  geom_point(size = 3, color = "black")+
  theme_minimal()+
  ylab("Relative Abundance (%)")+
  scale_fill_manual(values = npg_colors2[c(12,1,6)], name = "Carcass Access")+
  scale_color_manual(values = npg_colors2[c(12,1,6)], name = "Carcass Access")+
  scale_shape_manual(values = c(21, 22, 24), name = "Carcass Access")+
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical")

ggsave("outputs/SI N enzymes fig.png", height = 4, width = 5, units = "in",
       dpi = 350, bg = "white")

################################################################################

#### CARBON ENZYMES #####

carb <- interest[1:5,4]
c2 <- unlist(strsplit(carb, ", "))


######

# Step 2: Extract only genes of interest:
carb_reads <- kegg %>%
  filter(keggfxn %in% c2) %>%  
  pivot_longer(-keggfxn, names_to = "sample", values_to = "reads") %>%
  mutate(total_reads = total_reads_per_sample[sample],
         rel_abundance = (reads / total_reads)*100) %>%
  group_by(sample) %>%
  summarise(totalabund = sum(rel_abundance)) 

carb2 <- left_join(carb_reads, meta, by = "sample")

# are they significantly different? (if yes, do real stats)
carbmod <- aov(totalabund ~ trt, data = carb2)
anova(carbmod)
# trt***

emmeans(carbmod, pairwise ~ trt)
# C enzymes lower under E than EC and S than SC, but not diff at ctrls or btwn S, E


############ Do the stats right: stage 1
colnames(carb2)
# change rep to "block"
# need plot: E/S--first part of trt
# need splitplot: T/C
# need sample_site: block_plot_splitplot
# need time_fct

carb_stat <- carb2 %>%
  mutate(block = as.factor(rep),
         plot = str_sub(trt, 1, 1),
         splitplot = case_when(trt == "E" ~ "T",
                               trt == 'EC' ~ "C",
                               trt == "S" ~ "T",
                               trt == "SC" ~ "C"),
         sample_site = paste(block, plot, splitplot, sep = "_"),
         time_fct = as.factor(time),
         rel_abundance = as.numeric(totalabund),
         TREAT = case_when(trt == "EC" ~ "Control",
                           trt == "E" ~ "Excluded",
                           trt == "SC" ~ "Control",
                           trt == "S" ~ "Open"))

enzC <- carb_stat[,c(1,13)]
# for quantifying:
carb.agg <- carb_stat %>%
  group_by(trt) %>%
  summarise(mean = mean(rel_abundance, na.rm = T),
            SE = std.error(rel_abundance))
# percent change = -16% for E, -10% for S

# split into sites for stage 1:
carbA <- carb_stat %>% filter(site == "A")
carbB <- carb_stat %>% filter(site == "B")
carbG <- carb_stat %>% filter(site == "G")
carbT <- carb_stat %>% filter(site == "T")

corr_str1 = corAR1(form = ~ time|block/plot/sample_site, value = 0.2, fixed = FALSE)

hist(carbT$rel_abundance) 

# remove any outliers:
carbA2 <- carbA %>% filter(sample != "1010")

# Fit the model using lme
m0 <- lme(rel_abundance ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = carbT, na.action = na.exclude)

m1 <- lme(log(rel_abundance) ~ time_fct*plot*splitplot,
          random = ~ 1 | block/plot/sample_site,
          corr = corAR1(form = ~ time | block/plot/sample_site, value = 0.2, fixed = FALSE),
          data = carbG, na.action = na.exclude)

AIC(m0, m1)


# Print diagnostic plots to console to assess fit
qqnorm(m0) 

carbA$residuals <- residuals(m0)  # Assign residuals directly
res_df <- carbA[order(-abs(carbA$residuals)), ]  

# # View top residuals
head(res_df, 10)
# A: drop 1010, 1200
# T: drop 1104


shapiro.test(m0$residuals) #

anova(m0, type = "marginal") 
# A/low: time:plot:split*** F4,74=9.3
# B: time:split* F4,22=3.1
# G: time:split*** F4,79 = 5.2
# T: time:split*** F4,78=8.5

# Calculate emmeans
eres <- emmeans(m0, ~ splitplot | time_fct, type = 'link')
pairs(eres)  
# A (low): no diffs at T0, EC>ET at T1, no diffs T2, EC>ET and SC>ST at T3, no diffs T4
# B (no): C>T at T2
# G (high): C>T at T0,2,3
# T (mod): C>T at T1

## get emmeans to write to table
e = emmeans(m0, ~ splitplot*plot | time_fct, type = 'response') 
edf <- as.data.frame(summary(e))

#hiemm <- edf
#hiemm$site <- "High"
# 
#modemm <- edf
#modemm$site <- "Mod"
# 
#lowemm <- edf
#lowemm$site <- "Low"
# 
#ctrlemm <- edf
#ctrlemm$site <- "Ctrl"

emmeans_carb <- rbind(hiemm, modemm, lowemm, ctrlemm)
write.xlsx(emmeans_carb, "outputs/emmeans_C enzymes.xlsx")

####### Stage 2: #########
# only want interactings including site, since this is for ACROSS SITE comparisons
hist(emmeans_carb$emmean) # some left skew
# sampling_loc = splitplot, carcsas_loc = plot
emmeans_carb$time <- as.numeric(emmeans_nit$time_fct)

carbmod1 <- gls(emmean ~ site*splitplot + 
                  site*plot +
                  site*time_fct,
                correlation = corAR1(form = ~1|time),
                weights = ~ I(1/SE), 
                data = emmeans_carb)

carbmod2 <- gls(log(emmean) ~ site*splitplot + 
                  site*plot +
                  site*time_fct,
                correlation = corAR1(form = ~1|time),
                weights = ~ I(1/SE), 
                data = emmeans_carb)

AIC(carbmod1, carbmod2) # mod1 is better

plot(carbmod1) # 
anova(carbmod1, type = "marginal")
# time*** F4,52=15.5, split* F1,52=6.7, site* F3,52=3.9

emmeans_results <- emmeans(carbmod1, ~ site)
pairs(emmeans_results)
# time: decrease T0-T1 and increase T3-T4
# split: C>T
# site: diffs between Low and all other sites


######### Figures:
carb.plot <- carb_stat %>%
  group_by(time, TREAT) %>%
  summarise(relabund = sum(rel_abundance),
            SE = std.error(rel_abundance)) %>%
  mutate(Day = case_when(time == 0 ~ 0,
                         time == 1 ~ 5,
                         time == 2 ~ 10,
                         time == 3 ~ 15,
                         time == 4 ~ 30))

ggplot(data = carb.plot, aes(x=Day, y=relabund, color = TREAT, fill = TREAT, group = TREAT, shape = TREAT))+
  geom_errorbar(aes(ymax=relabund+SE, ymin=relabund-SE), width = 0.5)+
  geom_line()+
  geom_point(size = 3, color = "black")+
  theme_minimal()+
  ylab("Relative Abundance (%)")+
  scale_fill_manual(values = npg_colors2[c(12,1,6)], name = "Carcass Access")+
  scale_color_manual(values = npg_colors2[c(12,1,6)], name = "Carcass Access")+
  scale_shape_manual(values = c(21, 22, 24), name = "Carcass Access")+
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical")

ggsave("outputs/SI C enzymes fig.png", height = 4, width = 5, units = "in",
       dpi = 350, bg = "white")


#### data for all enzymes:
colnames(enzP)[7] <- "P_relabund"
colnames(enzN)[2] <- "N_relabund"
colnames(enzC)[2] <- "C_relabund"

enz <- full_join(enzP, enzN, by = "sample")
enz2 <- full_join(enz, enzC, by = "sample")
write.xlsx(enz2, "outputs/enzyme relabund.xlsx")

