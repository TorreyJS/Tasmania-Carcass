#################################################################################
# FunGuild assignment, plotting, statistics (both summer AND both seasons)
# 
# 1/15/25
#################################################################################
devtools::install_github("brendanf/FUNGuildR")

library(FUNGuildR);library(biomformat);library(phyloseq);library(funrar)
library(dplyr);library(stringr);library(tidyr);library(openxlsx)
library(lme4);library(car);library(emmeans);library(vegan)

#First, there is some data processing that has to be done before FUNGuild can run. 
# I'm using the un-rarefied data since we won't be doing diversity, just using 
# relative abundance data

#############################
location <- "Data/ITS_ps_soilonly.rds"
summerps <- readRDS(location) # unrarefied

# Check the updated sample data
sample_data(summerps)
unique(sample_data(summerps)$time) #0, 4, and 'l'--blank?

#############################
fung_otus <- data.frame(t(otu_table(summerps))) #27291 obs of 164 vars
taxa2 <- data.frame(tax_table(summerps))
its_asvs <- cbind(rownames(fung_otus), data.frame(fung_otus, row.names=NULL))

colnames(its_asvs)[1] <- "OTUID"

# pull off all the prefixes
taxa3 <- taxa2 %>%
  mutate_all(~ str_remove(., "^[a-z]__"))

# make a column with the full taxonomy, and make a column for ASV ID
taxa3 <- taxa3 %>% 
  mutate(Taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";"),
         id = rownames(taxa3))

# add taxonomy column to asvs dataframe too
its_asvs$Taxonomy <- taxa3$Taxonomy

#Now fun FUNGuild!
guilds <- funguild_assign(its_asvs) #27291 obs

#Now lets take out everything that did not assign to a guild or assigned with low confidence.
guilds2<-guilds[!(guilds$confidenceRanking=="Possible"),] 
guilds2 <- guilds2[!is.na(guilds2$guild),] #9478 obs

write.xlsx(guilds2, "outputs/raw funguild_hi confidence.xlsx")

# Sadly, we lose most of the data this way. This is inevitable, unfortunately. 
# This information simply isn't know about most fungi
guilds2 <- read.xlsx("outputs/raw funguild_hi confidence.xlsx")


### GUILD -- one option for analyzing FunGuild output

#
#


#Now we need to check and see how many sequences we have left for each sample and rarefy again
guilds3 <- guilds2 %>%
  group_by(guild) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

guilds.matrix <- t(guilds3[,c(2:165)]) #columns 2:end
rowSums(guilds.matrix)
min(rowSums(guilds.matrix)) # 567 is the minimum

### IF we choose to rarefy, use the below code:
#guilds.matrix <- rrarefy(guilds.matrix, 189)
#Now check to make sure all of the depths are equal
#rowSums(guilds.matrix2)

#Now for some data processing to get it ready for analysis
guilds.rel <- data.frame(make_relative(guilds.matrix))
rowSums(guilds.rel) # should be 1, because relative abundance sums to 1

colnames(guilds.rel) <- guilds3$guild # assign column names to guilds.rel

c <- data.frame(sample_data(summerps))
c <- c[c$time != "l",]
guilds4 <- cbind(guilds.rel, c)

# now, pivot longer
guilds5 <- guilds4 %>%
  pivot_longer(
    cols = 1:120, #1:end of names (not including metadata)
    names_to = "Guild",
    values_to = "Abundance",
    names_transform = list(Guild = as.factor)
  )

## summarise
guilds6 <- guilds5 %>%
  group_by(Guild, site, block, time, trt1, trt2) %>% # metadata columns
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')

write.xlsx(guilds6, "outputs/FunGuild-guild.xlsx")

guilds6 <- read.xlsx("outputs/FunGuild.xlsx")

# change to percents rather than decimals (*100) and aggregate
guilds_avg <- guilds6 %>%
  mutate(Abundance_scaled = Abundance * 100) %>%
  group_by(Guild) %>%
  summarise(Mean_Abundance = mean(Abundance_scaled, na.rm = TRUE), .groups = 'drop')

# reshape data to wide format
guilds7 <- guilds6 %>%
  pivot_wider(names_from = Guild, values_from = Abundance)

guild_bray <- vegdist(guilds7[,c(6:120)]) #exclude metadata

#set perm to account for block
perm <- how(nperm = 1000) 
setBlocks(perm) <- with(guilds7, block)

# make a single treatment column
guilds7$trt <- paste0(guilds7$trt1, guilds7$trt2)


############## looking at ectomycorrhizae:
myco_cols <- grep("mycorrhizal", colnames(guilds7), ignore.case = TRUE, value = TRUE)

myco <- guilds7 %>%
  select(1:5, all_of(myco_cols))

write.xlsx(myco, "outputs/mycorrhizal fungi.xlsx")
myco <- read.xlsx("outputs/mycorrhizal fungi.xlsx")

myco <- myco %>%
  mutate(sum_myco = rowSums(select(., 6:22), na.rm = TRUE),
         trt = paste0(trt1, trt2))

### is the SUM of mycorrhizae different between treatments??
hist(myco$sum_myco) # pretty normal, nice
myco$time <- as.factor(myco$time)
myco_final <- myco %>% filter(time == "4")

# how many samples have EM? AM?
AM <- myco_final %>% filter(`Arbuscular.Mycorrhizal` != 0) #38 of 82
mean(AM$`Arbuscular.Mycorrhizal`) # average 0.3% relative abundance
meanAM <- AM %>%
  group_by(site, trt1, trt2) %>%
  summarise(mean = mean(`Arbuscular.Mycorrhizal`)) # no apparent site or trt diffs


EM <- myco_final %>% filter(`Ectomycorrhizal` != 0) #82 of 82

EM <- EM %>%
  mutate(trt3 = case_when(trt == "EC" ~ "C",
                          trt == "SC" ~ "C",
                          trt == "ST" ~ "S",
                          trt == "ET" ~ "E"))
mean(EM$Ectomycorrhizal) # average 22.5% relative abundance
meanEM <- EM %>%
  group_by(trt3) %>%
  summarise(mean = mean(Ectomycorrhizal)*100,
            SE = std.error(Ectomycorrhizal)*100) #MUCH higher at control vs trt! (~2.4-4x higher). lowest at E


## do stats on this!!!
# for bounded variables, beta regression is best
myco_final <- myco_final %>%
  mutate(Treatment = paste0(trt1, trt2),
         block = as.factor(block))

hist(myco_final$Ectomycorrhizal)

EMmod <- aov(log(`Ectomycorrhizal`) ~ trt, 
               data = myco_final) #trt p<0.001, F3,78=6.5
hist(resid(EMmod))
shapiro.test(resid(EMmod))


summary(EMmod)
car::Anova(EMmod, type = 3)

# emmeans
emmeans(EMmod, pairwise ~ trt)
# EC>ET p<0.05, t=2.8
# SC>ST p<0.01, t=3.4
# no diffs at time 0, only at later times = these are not underlying diffs!!


an <- aov(sum_myco ~ trt, data = myco_final)
summary(an)

emmeans(an, pairwise ~ trt)
# total mycorrhizae ARE lower at ET than EC and ST than SC!!!! 
# p<0.001, F3,78 = 5.1

## what about SPECIFICALLY just AM/EM?
amem <- c("Ectomycorrhizal","Arbuscular.Mycorrhizal")

myco_amem <- myco_final %>% 
  select(1:5, all_of(amem)) %>%
  mutate(trt = paste0(trt1, trt2))

aa <- aov(`Arbuscular.Mycorrhizal` ~ trt, data = myco_amem)
summary(aa)

emmeans(aa, pairwise ~ trt) # EM p<0.001, F3,78=6.9. EC>ET (0.002), SC>ST (0.054)
# AM trt p<0.05, F3,78=2.9. no diffs btwn treatments
# probably not reportable--our ITS primers are biased against AM amplification

#Run model (test site, treatment, time)
adonis2(guild_bray~site, data=guilds7, permutation=perm) #site***
adonis2(guild_bray~trt, data=guilds7, permutation=perm) #trt*
adonis2(guild_bray~time, data=guilds7, permutation=perm) # time***

adonis2(guild_bray~site*time*trt, data=guilds7, permutation=perm)
# time:trt**, site***

# Interaction of time and trt, and effect of site--can parse out further later

##### Figures ####

########## Barchart ##############
# Step 1: Sum total abundance for each guild
total_abundance_by_guild <- aggregate(Abundance ~ site+time+trt1+trt2+Guild, guilds6, FUN = sum)

# Step 2: Sum total abundance across all guilds
total_abundance <- sum(total_abundance_by_guild$Abundance) #164

# Step 3: Calculate relative abundance for each guild
total_abundance_by_guild$RelativeAbundance <- (total_abundance_by_guild$Abundance / total_abundance) * 100

## can plot by a single guild if one stands out -- return to this later 

#g <- total_abundance_by_guild[total_abundance_by_guild$id != "78",]
#g <- g[g$id != "1",]
#g <- g[g$Guild == "Plant Pathogen",]
g <- total_abundance_by_guild[total_abundance_by_guild$Guild == "Plant Pathogen",]
library(plotrix)
#agg <- aggregate(g$RelativeAbundance, by = list(g$Time, g$Crop), FUN = mean) 
agg <- aggregate(g$RelativeAbundance, by = list(g$Crop), FUN = mean) 
#se <- aggregate(g$RelativeAbundance, by = list(g$Time, g$Crop), FUN = std.error)
se <- aggregate(g$RelativeAbundance, by = list(g$Crop), FUN = std.error)
agg$se <- se$x #add se maomc to sagg df
#colnames(agg)[1:4] <- c("Time","Crop","mean", "se")
colnames(agg)[1:3] <- c("Crop","mean", "se")



ggplot(data=agg, aes(x=Crop, y=mean, fill=Crop)) +
  geom_col(data=agg, width=0.8, position=position_dodge(),color="black", size= 0.7,aes(fill = Crop)) +
  scale_color_manual(values= Bpalc4)+
  scale_fill_manual(values= Bpalc4) +
  geom_errorbar(width = 0.8,position=position_dodge(),aes(ymin = (mean-se), ymax = (mean+se))) +
  #scale_y_continuous(limits = c(0,2.2), expand = c(0, 0)) +
  ylab("Plant Pathogen Relative Abundance") +
  xlab("") +
  theme_lildino() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
#facet_wrap(~Time)

#You can repeat this with as many guilds as you are interested in. 


#################################################################################
#### looking at trophic MODE rather than guild ####
library(plyr)

## these are all different ways to group the data--trophic mode is the broadest

####### TROPHIC MODE

##
##

unique(guilds2$trophicMode)
# remove spaces--there should only be 8
guilds2 <- guilds2 %>%
  mutate(trophicMode = gsub(" ", "", trophicMode))
unique(guilds2$trophicMode) # goooood

guilds3 <- guilds2 %>%
  group_by(trophicMode) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

# sample depth
guilds.matrix <- t(guilds3[,c(2:165)]) # get rid of blank
rowSums(guilds.matrix)
min(rowSums(guilds.matrix)) #567


#Now for some data processing to get it ready for analysis
guilds.rel <- data.frame(make_relative(guilds.matrix))
rowSums(guilds.rel)# all 1 -- make_relative worked
colnames(guilds.rel) <- as.character(guilds3[[1]])
guilds4 <- guilds.rel

guilds5 <- cbind(guilds4, c)

write.xlsx(guilds5, "outputs/trophic mode.xlsx")

guilds5 <- guilds5 %>%
  mutate(Treatment = case_when(trt1 == "E" & trt2 == "C" ~ "Control",
                               trt1 == "E" & trt2 == "T" ~ "Excluded",
                               trt1 == "S" & trt2 == "C" ~ "Control",
                               trt1 == "S" & trt2 == "T" ~ "Open"),
         Site = case_when(site == "A" ~ "Low-Impact",
                          site == "B" ~ "Control",
                          site == "G" ~ "High-Impact",
                          site == "T" ~ "Moderate-Impact"))
guilds5$Site <- factor(guilds5$Site, levels = c("Control","High-Impact",
                                                "Moderate-Impact","Low-Impact"))

# pivot longer
guilds6 <- gather(guilds5, trophicmode, Abundance, 1:8, factor_key=TRUE) #columns 1:8 if trophic mode
guilds6$trophicmode <- as.character(guilds6$trophicmode)


# group dataframe by phylum, calculate mean rel. abundance
means <- guilds6 %>%
  group_by(trophicmode) %>%
  summarise(mean = mean(Abundance, na.rm = TRUE))

means # 115 for guild, 8 for trophic mode

# find Phyla whose mean rel. abund. is less than 1%
remainder <- means[means$mean <= 0.01,]$trophicmode # 98 of 115 OR 1 of 8

# change their name to "Other" -- OR just leave it if only 1 (trophic mode)
#guilds6[guilds6$trophicmode %in% remainder,]$trophicmode <- 'Other'



# aggregate
guilds7 <- aggregate(Abundance~trophicmode+Site+time+Treatment, guilds6, FUN=sum)

aggregate(Abundance*100~trophicmode, guilds7, FUN=mean)

## cut to just time final (control acts as "time 0")

guilds8 <- guilds7 %>%
  filter(time == "4")

library(ggsci)
npg_colors <- (pal_npg("nrc", alpha =1)(10))

ggplot(guilds8, aes(fill=trophicmode, y=Abundance, x=Treatment)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Abundance") +
  scale_fill_manual(values = npg_colors) +
  theme(axis.text.x=element_text(angle=45, hjust= 0.98, vjust=1)) +
  theme(legend.box.margin=margin(-10,-15,-10,-10))+
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  #xlab("Diversity") +
  labs(fill = "Trophic Mode")+
  facet_wrap(.~Site, ncol = 2)



##### Statistics ####
# Now lets look at individual trophic modes
# could use guilds5 for this, but it doesn't have the "Other" assignments yet
# we coul pivot guilds6 (w/other) wider, but it's a huge pain (doesn't work)
# OR, we can use guild5, but only use the categories from guild6
library(glmmTMB)

guilds5$block <- as.factor(guilds5$block)
gstats <- guilds5 %>% filter(time == 4) # stats only on time final

unique(guilds6$trophicmode)
# [1] "Other"                             "Pathotroph"                       
# [3] "Pathotroph-Saprotroph"             "Pathotroph-Saprotroph-Symbiotroph"
# [5] "Pathotroph-Symbiotroph"            "Saprotroph"                       
# [7] "Saprotroph-Symbiotroph"            "Symbiotroph"                      


# go one mode at a time
hist(gstats$`Symbiotroph`) # looks log
# for Pathotroph, get rid of one observation = to 0
#g2 <- gstats %>% filter(Pathotroph > 0)
#g3 <- gstats %>% filter(`Pathotroph-Saprotroph-Symbiotroph`>0)
#g4 <- gstats %>% filter(`Pathotroph-Symbiotroph`>0)

# for bounded variables, beta regression is best
library(glmmTMB)
mod <- glmmTMB(`Symbiotroph` ~ Treatment * Site + (1|Site/block), 
               data = gstats, 
               family = beta_family(link = "logit"))

# bacterial pathogens
ps3 <- ps3 %>%
  mutate(Treatment = case_when(trt=="E" ~ "Excluded",
                               trt=="EC" ~ "Control",
                               trt=="S" ~ "Open",
                               trt=="SC" ~ "Control"))

hist(ps3$abundance)
mod <- lme(log(abundance) ~ Treatment * site, 
           random = ~1|site/block, 
           data = ps3)

car::Anova(mod, type = 3)

# emmeans
emmeans(mod, pairwise ~ Treatment) # E***, S***


#################################################################################

##### Seasonal analysis #######
#############################
location <- "Data/both seasons_ITS_phyloseq_NOTrarefied.rds" 
seasonps <- readRDS(location) # not rarefied

# Check the updated sample data
sample_data(seasonps)

#############################
fung_otus <- data.frame(t(otu_table(seasonps))) #38136 obs of 174 vars
taxa2 <- data.frame(tax_table(seasonps))
its_asvs <- cbind(rownames(fung_otus), data.frame(fung_otus, row.names=NULL))

colnames(its_asvs)[1] <- "OTUID"

# pull off all the prefixes
taxa3 <- taxa2 %>%
  mutate_all(~ str_remove(., "^[a-z]__"))

# make a column with the full taxonomy, and make a column for ASV ID
taxa3 <- taxa3 %>% 
  mutate(Taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";"),
         id = rownames(taxa3))

# add taxonomy column to asvs dataframe too
its_asvs$Taxonomy <- taxa3$Taxonomy

#Now fun FUNGuild!
guilds <- funguild_assign(its_asvs) #38136 obs

#Now lets take out everything that did not assign to a guild or assigned with low confidence.
guilds2<-guilds[!(guilds$confidenceRanking=="Possible"),] 
guilds2 <- guilds2[!is.na(guilds2$guild),] #14045 obs

# Sadly, we lose most of the data this way. This is inevitable, unfortunately. 
# This information simply isn't know about most fungi

#Now we need to check and see how many sequences we have left for each sample and rarefy again
guilds3 <- guilds2 %>%
  group_by(guild) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

guilds.matrix <- t(guilds3[,c(2:175)]) #columns 2:end
rowSums(guilds.matrix)
min(rowSums(guilds.matrix)) # 574 is the minimum

### IF we choose to rarefy, use the below code:
#guilds.matrix <- rrarefy(guilds.matrix, 189)
#Now check to make sure all of the depths are equal
#rowSums(guilds.matrix2)

#Now for some data processing to get it ready for analysis
guilds.rel <- data.frame(make_relative(guilds.matrix))
rowSums(guilds.rel) # should be 1, because relative abundance sums to 1

colnames(guilds.rel) <- guilds3$guild # assign column names to guilds.rel

c <- data.frame(sample_data(seasonps))
guilds4 <- cbind(guilds.rel, c)

# now, pivot longer
guilds5 <- guilds4 %>%
  pivot_longer(
    cols = 1:135, #1:end of names (not including metadata)
    names_to = "Guild",
    values_to = "Abundance",
    names_transform = list(Guild = as.factor)
  )

# time = 4 = summer; time = 2 = winter
guilds5 <- guilds5 %>%
  mutate(Season = case_when(time == "4" ~ "Summer",
                            time == "2" ~ "Winter"))

## summarise
guilds6 <- guilds5 %>%
  group_by(Guild, site, trt, Season, rep) %>% # metadata columns
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')

write.xlsx(guilds6, "outputs/FunGuild_seasonal.xlsx")
guilds6 <- read.xlsx("outputs/FunGuild_seasonal.xlsx")

# change to percents rather than decimals (*100) and aggregate
guilds_avg <- guilds6 %>%
  mutate(Abundance_scaled = Abundance * 100) %>%
  group_by(Guild) %>%
  summarise(Mean_Abundance = mean(Abundance_scaled, na.rm = TRUE), .groups = 'drop')


#Now we need to check and see how many sequences we have left for each sample 
library(plyr);library(dplyr)

## these are all different ways to group the data--trophic mode is the broadest
guilds3 <- ddply(guilds2,"trophicMode",numcolwise(sum))
#guilds3 <- ddply(guilds2,"taxon",numcolwise(sum))
#guilds3 <- ddply(guilds2,"guild",numcolwise(sum))

# sample depth
guilds.matrix <- t(guilds3[,c(2:175)])
rowSums(guilds.matrix)
min(rowSums(guilds.matrix)) #574


#Now for some data processing to get it ready for analysis
guilds.rel <- data.frame(make_relative(guilds.matrix))
rowSums(guilds.rel)# all 1 -- make_relative worked
colnames(guilds.rel) <- guilds3[,1]
guilds4 <- guilds.rel

guilds5 <- cbind(guilds4, c)
unique(guilds5$trt)

guilds5 <- guilds5 %>%
  mutate(Treatment = case_when(trt == "EC" ~ "Control",
                               trt == "E" ~ "Excluded",
                               trt == "S" ~ "Control",
                               trt == "SC" ~ "Open",
                               trt1 == "E" & trt2 == "C" ~ "Control",
                               trt1 == "E" & trt2 == "T" ~ "Excluded",
                               trt1 == "S" & trt2 == "C" ~ "Control",
                               trt1 == "S" & trt2 == "T" ~ "Open"),
         Site = case_when(site == "A" ~ "Low-Impact",
                          site == "B" ~ "Control",
                          site == "G" ~ "High-Impact",
                          site == "T" ~ "Moderate-Impact"),
         Season = case_when(time == 4 ~ "Summer",
                            time == 2 ~ "Winter"))

guilds5$Site <- factor(guilds5$Site, levels = c("Control","High-Impact",
                                                "Moderate-Impact","Low-Impact"))

# pivot longer
guilds6 <- gather(guilds5, trophicmode, Abundance, 1:12, factor_key=TRUE) #columns 1:10 if trophic mode
guilds6$trophicmode <- as.character(guilds6$trophicmode)


# group dataframe by phylum, calculate mean rel. abundance
means <- ddply(guilds6, ~trophicmode, function(x) c(mean=mean(x$Abundance)))
means # 12 for trophic mode

# find Phyla whose mean rel. abund. is less than 1%
remainder <- means[means$mean <= 0.01,]$trophicmode # 6 of 12

# change their name to "Other"
guilds6[guilds6$trophicmode %in% remainder,]$trophicmode <- 'Other'



# aggregate
str(guilds6)
guilds7 <- aggregate(Abundance~trophicmode+Season+Site+Treatment, guilds6, FUN=sum)

aggregate(Abundance*100~trophicmode, guilds7, FUN=mean)

library(ggsci);library(ggplot2)
npg_colors <- (pal_npg("nrc", alpha =1)(10))
customcolors <- c(
    "Other" = "#B09C85FF", 
    "Pathotroph" = "#91D1C2FF",
    "Pathotroph-Saprotroph" = "#F39B7FFF",
    "Pathotroph-Saprotroph-Symbiotroph" = "#3C5488FF",
    "Pathotroph-Symbiotroph" = "#8491B4FF",
    "Saprotroph" = "#E64B35FF",
    "Saprotroph-Symbiotroph" = "#4DBBD5FF",
    "Symbiotroph" = "#00A087FF")

ggplot(guilds7, aes(fill = trophicmode, y = Abundance, x = Treatment)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = customcolors) +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "italic"),  
        strip.placement = "inside",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + 
  labs(fill = "Trophic Mode") +
  facet_grid(Site ~ Season)  # Hide Season labels


ggsave("outputs/seasonal ITS Funguild.png", height = 8, width = 12,
       units = "in", dpi = 350, bg = "white")


##### Statistics ####
# Now lets look at individual trophic modes
# could use guilds5 for this, but it doesn't have the "Other" assignments yet
# we could pivot guilds6 (w/other) wider, but it's a huge pain (doesn't work)
# OR, we can use guild5, but only use the categories from guild6
library(glmmTMB)
guilds5$rep <- as.character(guilds5$rep)
guilds5 <- guilds5 %>%
  mutate(block = coalesce(block, rep))

guilds5$block <- as.factor(guilds5$block)

unique(guilds6$trophicmode)
# [1] "Other"                                              
# [3] "Pathotroph-Saprotroph"             "Pathotroph-Saprotroph-Symbiotroph"
# [5] "Pathotroph-Symbiotroph"            "Saprotroph"                       
# [7] "Saprotroph-Symbiotroph"            "Symbiotroph"                      


# go one mode at a time
hist(guilds5$`Saprotroph-Symbiotroph`) # looks log
# for Pathotroph, get rid of one observation = to 0
#g2 <- guilds5 %>% filter(Pathotroph > 0)
#g3 <- guilds5 %>% filter(`Pathotroph-Saprotroph-Symbiotroph` >0)
#g4 <- guilds5 %>% filter(`Pathotroph-Symbiotroph` >0)

# for bounded variables, beta regression is best
mod <- glmmTMB(`Saprotroph-Symbiotroph` ~ Treatment * Season + (1|Site/block), 
               data = guilds5, 
               family = beta_family(link = "logit"))

car::Anova(mod, type = 3)

# emmeans
emmeans(mod, pairwise ~ Treatment|Season)

#################################################################################
