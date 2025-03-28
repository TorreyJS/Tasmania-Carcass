################################################################################
# PICRUSt2.0 16S copy # results
# 1/23/25
#
################################################################################

## the goal is to estimate the 16S copy number for each sequence
# this is what picrust does. so, I need to find the output that has just the info 
# for sequence and 16S_rRNA_Count

### code from Ernie's github (https://github.com/eosburn/Coweeta-Microbes), edited:

library(biomformat);library(phyloseq);library(ggplot2);library(dplyr)
library(vegan);library(car);library(lme4);library(reshape);library(emmeans)
library(cowplot)

Bac_otus2 <- "Data/summer_16s_filtered_OTU_table.biom"

x1 <- read_biom(Bac_otus2)
x2 <- as.matrix(biom_data(x1))

head(x2)
colSums(x2) # all just about 10k

copies2 <- nsti_data <- read.table("Data/summer_marker_predicted_and_nsti.tsv.gz", header = TRUE, sep = "\t")
x3.1 <- x2[rownames(x2) %in% copies2$sequence,]
colSums(x3.1)

min(colSums(x3.1)) #9869
x4.1 <- t(x3.1)
x5.1 <- rrarefy(x4.1, 9869)

rowSums(x5.1) # now all are 9869

library(funrar)

x6.1 <- make_relative(x5.1)
rowSums(x6.1) # 1
copies2 <- copies2[order(copies2$sequence),]
x7.1 <- t(x6.1)
x8.1 <- cbind(rownames(x7.1), data.frame(x7.1, row.names=NULL,check.names=FALSE))
colnames(x8.1)[1] <- c("sequence")
x8.1 <- x8.1[order(x8.1$sequence),]
x9.1 = sapply(x8.1[,c(2:401)], '*', copies2$X16S_rRNA_Count) 
colSums(x9.1) # about 1.3 to 1.5
table.1 <- data.frame(t(x9.1))
table.1$copynumber <- rowSums(table.1)
table2.1 <- cbind(rownames(table.1), data.frame(table.1, row.names=NULL,check.names=FALSE))
colnames(table2.1)[1] <- c("id")
table2.1 <- table2.1[order(table2.1$id),]

meta <- read.csv("Data/summer_metadata.csv")
colnames(meta)[1] <- "id"
meta <- meta[order(meta$id),]

table2.1 <- cbind(table2.1, meta)
table3.1 <- table2.1[,c(1,29417:29425)] # table with just copy numbers and metadata

write.xlsx(table3.1, "outputs/16S gene copy numbers.xlsx")
table3.1 <- read.xlsx("outputs/16S gene copy numbers.xlsx")

table3.1 <- table3.1[,c(2,8,9,10)]
sumtab <- table3.1 %>%
  group_by(site, time, trt) %>%
  summarise(copies = mean(copynumber),
            SE = std.error(copynumber))

table3.1_time2 <- table3.1 %>% filter(time == 2)

sumtab2 <- table3.1_time2 %>%
  mutate(TREAT = case_when(trt == "EC" ~ "Control",
                           trt == "E" ~ "Carcass",
                           trt == "SC" ~ "Control",
                           trt == "S" ~ "Carcass")) %>%
  group_by(TREAT) %>%
  summarise(copies = mean(copynumber),
            SE = std.error(copynumber))
# day 10: 1.51 +/- 0.04 (T) vs 1.32 +/- 0.01 (C)


# all times except 0: 1.4858 (T) vs 1.3587 (C) or 9.35% higher


##### STATISTICS: ####

# restructure 
table3.1 <- table3.1 %>%
  mutate(plot = substr(trt,1,1),
         splt = substr(trt,2,2),
         splitplot = case_when(splt == "C" ~ "C",
                               splt == "" | is.na(splt) ~ "T"),
         time_fct = as.factor(time),
         time_num = as.numeric(time),
         block = substr(ID_2, nchar(ID_2), nchar(ID_2)),
         block = as.factor(block),
         sample_site = paste(block, plot, splitplot, sep = "_"))


#split into sites:
# stage 1: separate by site
ctrltab <- table3.1 %>% filter(site == "B")
hitab <- table3.1 %>% filter(site == "G")
modtab <- table3.1 %>% filter(site == "T")
lowtab <- table3.1 %>% filter(site == "A")

## fit the model -- do this 4x
library(nlme)

hist(lowtab$copynumber) # ctrl: log, hi, mod, low: sort of log?

mod <- lme(log(copynumber) ~ time_fct*plot*splitplot, 
           random = ~1|block/plot/sample_site, 
           corr = corAR1(form = ~time_num|block/plot/sample_site,
                         value = 0.2, fixed = F), 
           data = hitab2, na.action = na.exclude)

# check fit
model_resids = residuals(mod)
shapiro.test(model_resids) # p 
qqnorm(mod$residuals, main = "Normal Q-Q Plot", xlab = "Residuals", ylab = "Expected values")
qqline(mod$residuals) # looks pretty good, but can try dropping high outliers to get normal distribution

residuals <- resid(mod)
highs <- order(abs(residuals), decreasing = T)
highs # hi: 4, 35, 16, 116; mod: 3, 7, 113, 9 ; low: 100, 12
hitab2 <- hitab %>% filter(!(id %in% c(1012, 1070, 1039, 1387)))
modtab2 <- modtab %>% filter(!(id %in% c(1075, 1082, 1400, 1086)))
lowtab2 <- lowtab %>% filter(!(id %in% c(1340, 1057, 1011, 1351)))

# Anova
Anova(mod, type = "3")
# ctrl: time:split** (p<0.01, chisq=13.4). log transform, no drops
# hi: time:split*** (p<0.001, chisq=30.1). log transform, drop 4 outliers
# mod: time:split* (p<0.05, chisq = 10.2). log transform, drop 4 outliers
# low: time:split* (p<0.05, chisq = 10.6). log transform, drop 4 outliers

# Calculate emmeans
emmeans(mod, pairwise ~ splitplot| time_fct, type = 'link')
# ctrl: T>C at time 2*
# hi: T>C at time 2**, C>T at time 3*
# moderate: T>C at time 1,2***
# low: T>C at time 2*
# copy numbers were higher at treatment sites than control
# r-strategists tend to have higher copy #s

e = emmeans(mod, ~ splitplot*plot | time_fct, type = 'response') 
edf <- as.data.frame(summary(e))

#hiemm <- edf
hiemm$site <- "High"
#modemm <- edf
modemm$site <- "Mod"
#lowemm <- edf
lowemm$site <- "Low"
#ctrlemm <- edf
ctrlemm$site <- "Ctrl"

emmeans_16scopy <- rbind(hiemm, modemm, lowemm, ctrlemm)
write.xlsx(emmeans_16scopy, "outputs/emmeans_16S copy numbers.csv")

emmeans_16scopy <- read.xlsx("outputs/emmeans_16S copy numbers.xlsx")

table3.1 <- table3.1 %>% 
  mutate(Treatment = case_when(trt == "EC" ~ "Control",
                               trt == "E" ~ "Excluded",
                               trt == "SC" ~ "Control",
                               trt == "S" ~ "Open"))

sum_operons <- aggregate(copynumber~ site+time+Treatment, data=table3.1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_operons <- do.call(data.frame, sum_operons)
sum_operons

sum_operons$se <- sum_operons$copynumber.sd / sqrt(sum_operons$copynumber.n)
head(sum_operons)

#Now for some plotting

library(ggplot2)
library(ggsci)

npg_colors <- (pal_npg("nrc", alpha =1)(10))

sum_operons$Treatment <- factor(sum_operons$Treatment, levels=c("Control","Open","Excluded"))
sum_operons <- sum_operons %>%
  mutate(site = case_when(site == "A" ~ "Medium Devils",
                          site == "B" ~ "No Devils",
                          site == "G" ~ "Low Devils",
                          site == "T" ~ "High Devils"))
sum_operons <- sum_operons %>%
  mutate(Day = case_when(time == 0 ~ 0,
                         time == 1 ~ 5,
                         time == 2 ~ 10,
                         time == 3 ~ 15,
                         time == 4 ~ 30))

sum_operons$Site <- factor(sum_operons$site, levels =c("Control","Low Devils","Medium Devils","High Devils"))

ggplot(sum_operons, aes(x=Day, y=copynumber.mean, color=Treatment, group = Treatment)) + 
  geom_point(size=3) +
  geom_line()+
  geom_errorbar(aes(ymin=copynumber.mean-se, ymax=copynumber.mean+se), width=.25, linewidth=.75) +
  labs(list(x ="")) +
  scale_color_manual(values = c(npg_colors[c(10,5,3)]))+
  theme_classic() +
  theme(axis.title=element_text(size=20),
        text=element_text(size=20),
        legend.text = element_text(color = "black", size = 15),
        axis.line.x=element_line(colour="black", size=1),
        axis.line.y=element_line(colour="black", size=1),
        legend.title=element_blank())+
  ylab(bquote('Average 16S operon number')) +
  facet_wrap(~ site)


ggsave("outputs/average 16S operon number.png", height = 6, width = 8,
       units = "in", dpi = 350, bg = "white")

#### stage 2: all sites together ####
# model structure: gls(emmean ~ site*carcass_loc+site*time_fct, correlation=corAR1(form=~1|time_num, weights=~I(1/SE), data=allsites))

emmeans_16scopy <- emmeans_16scopy %>%
  mutate(time_num = as.numeric(time_fct),
         TRT = case_when(plot == "E" & splitplot == "C" ~ "Control",
                         plot == "E" & splitplot == "T" ~ "Excluded",
                         plot == "S" & splitplot == "C" ~ "Control",
                         plot == "S" & splitplot == "T" ~ "Open"))

mod2 <- gls(response ~ site*plot + site*time_fct + time_fct*splitplot,
            correlation = corAR1(form=~1|time_num), weights=~I(1/SE), 
            data = emmeans_16scopy)

plot(mod2)

car::Anova(mod2, type = 3) 
# site:time* (p<0.05, X2=21.8), time* (p<0.05, X2=12.8), time:TRT*** (p<0.001, X2=27.8)

emmeans_results <- emmeans(mod2, ~ splitplot |time_fct , type = "response")
pairs(emmeans_results, adjust = "fdr")
# s:t - ctrl site > L,M,H at day 15
# time - increase at day 15 from day 0, then decrease at day 30 from day 15
# plot - E>S but not sig overall
# TRT: Excl>Ctrl overall
# for plot: trt by time: T0: no diffs, T1: E>ctrl and open, T2: E>ctrl, O>ctrl, T3: E>ctrl, T4: none
# trt>ctrl for 1*,2**,3*


#### stage 3: seasonal ####

### FIRST, bring in the 2022 data:

Bac_otus <- "Data/winter_biom_for_PICRUST.biom"

x1 <- read_biom(Bac_otus)
x2 <- as.matrix(biom_data(x1))

head(x2)
colSums(x2) # most 10950

copies2 <- nsti_data <- read.table("Data/winter_marker_predicted_and_nsti.tsv.gz", header = TRUE, sep = "\t")
x3.1 <- x2[rownames(x2) %in% copies2$sequence,]
colSums(x3.1)

min(colSums(x3.1)) # 10644
x4.1 <- t(x3.1)
x5.1 <- rrarefy(x4.1, 10644)

rowSums(x5.1) # now all are 10644

library(funrar)

x6.1 <- make_relative(x5.1)
rowSums(x6.1) # 1
copies2 <- copies2[order(copies2$sequence),]
x7.1 <- t(x6.1)
x8.1 <- cbind(rownames(x7.1), data.frame(x7.1, row.names=NULL,check.names=FALSE))
colnames(x8.1)[1] <- c("sequence")
x8.1 <- x8.1[order(x8.1$sequence),]
x9.1 = sapply(x8.1[,c(2:276)], '*', copies2$X16S_rRNA_Count) # to number of cols
colSums(x9.1) # about 1.3 to 1.5
table.1 <- data.frame(t(x9.1))
table.1$copynumber <- rowSums(table.1)
table2.1 <- cbind(rownames(table.1), data.frame(table.1, row.names=NULL,check.names=FALSE))
colnames(table2.1)[1] <- c("id")
table2.1 <- table2.1[order(as.numeric(table2.1$id)),]
table2.1$id <- as.numeric(table2.1$id)

meta <- read.csv("Data/16S_devils_metadata_winter.csv")
colnames(meta)[1] <- "id"
meta <- meta[order(as.numeric(meta$id)),]

table2.1 <- left_join(table2.1, meta, by = "id")

table3.1 <- table2.1[,c(1,11917:11925)] # table with just copy numbers and metadata
table3.2 <- table3.1 %>% filter(trt != "L" & trt != "LC")

library(openxlsx)
write.xlsx(table3.2, "outputs/2022_16S gene copy numbers.xlsx")

table3.2 <- read.xlsx("outputs/2022_16S gene copy numbers.xlsx")

library(plotrix)
sumtab <- table3.2 %>%
  group_by(site, time, trt) %>%
  summarise(copies = mean(copynumber),
            SE = std.error(copynumber))

table3.3 <- table3.2 %>% filter(time != 0)

sumtab2 <- table3.3 %>%
  mutate(TREAT = case_when(trt == "EC" ~ "Control",
                           trt == "E" ~ "Carcass",
                           trt == "SC" ~ "Control",
                           trt == "S" ~ "Carcass")) %>%
  group_by(TREAT) %>%
  summarise(copies = mean(copynumber),
            SE = std.error(copynumber))
