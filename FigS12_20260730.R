#####
#Suppp figures

rm(list=ls())

#create a function that fixes the display of p-values
disp_pval <- function(p){
  if(p<0.001){
    p <- scales::scientific(p, digits= 2)
  } else {
    p <- round(p, 2)
  }
  return(p)
}

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\FigS12\\")

# Library
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(reshape2)
library(readxl)
library(lmerTest)
library(emmeans)
library(survival)
library(ggsurvfit)
library(ggpubr)

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#load in database
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260717.xlsx'))

#filter out antibiotics samples
mrg1 <- mrg <- mrg %>% 
  filter(!grepl("Yes", Antibiotics))   #remove anitbiotics samples


###################################################
#                                                 #
#                 FIGURE S12A                     #
#                                                 #
###################################################

#extract the kids who are intermittent carriers of staph.
mrg1.sa <- mrg1[mrg1$Assay == 'Staphylococcus_aureus',] #first select only SA data
#in the data frame mrg1.sa, get the frequency of rows per StudyID for which the column InfectionLongerTrue is 1
mrg1.sa_freq <- mrg1.sa %>%
  group_by(StudyID) %>%
  summarise(freq = sum(InfectionLongerTrue == 1) / n()) #get the frequency of positive samples per child
mrg1.sa_freq$type <- '0.neverStaph'
mrg1.sa_freq$type[mrg1.sa_freq$freq>0] <- '1.intermittentStaph'
mrg1.sa_freq$type[mrg1.sa_freq$freq>0.5] <- '2.frequentStaph'

#add ages:
age_map <- mrg %>%
  distinct(StudyID, AgeYears)
mrg1.sa_freq2 <- mrg1.sa_freq %>%
  left_join(age_map, by = "StudyID")

m <- glm(freq ~ type + AgeYears,
         family = quasibinomial(link = "logit"),
         data = mrg1.sa_freq2)
summary(m)

FigS12A <- ggplot(mrg1.sa_freq, aes(x = type, y = freq, colour= type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
FigS12A

ggsave('FigS12A_freqStaphGroups.pdf', width = 5, height = 4)


###################################################
#                                                 #
#                 FIGURE S12B                     #
#                                                 #
###################################################

mrg_clear3 <- merge(
  mrg_clear2,
  mrg1.sa_freq[, c("StudyID", "type")],
  by = "StudyID",
  all.x = TRUE
)

mrg_clear3$type <- factor(mrg_clear3$type)
mrg_clear3$Assay <- factor(mrg_clear3$Assay)
mrg_clear3$StudyID <- factor(mrg_clear3$StudyID)

mrg_clear3$type <- relevel(factor(mrg_clear3$type), ref = "1.intermittentStaph") #use intermittent as ref

m <- lmer(log10_avgconc ~ type + (1 | StudyID) + (1 | Assay), data = mrg_clear3)

coefs <- coef(summary(m))

stats_type <- data.frame(
  term = rownames(coefs),
  Estimate = coefs[, "Estimate"],
  Pval = coefs[, "Pr(>|t|)"],
  row.names = NULL,
  stringsAsFactors = FALSE
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    group1 = "1.intermittentStaph",
    group2 = case_when(
      grepl("0\\.neverStaph", term) ~ "0.neverStaph",
      grepl("2\\.frequentStaph", term) ~ "2.frequentStaph",
      TRUE ~ NA_character_
    ),
    label = paste0(
      "Est.: ", round(Estimate, 3),
      ", P=", scales::scientific(Pval, 2)
    )
  )

stats_type$y.position <- seq(
  from = max(mrg_clear3$log10_avgconc, na.rm = TRUE) + 0.5,
  by = 0.35,
  length.out = nrow(stats_type)
)

mrg_clear3$type <- factor(
  mrg_clear3$type,
  levels = c("0.neverStaph", "1.intermittentStaph", "2.frequentStaph")
)

ggplot(mrg_clear3, aes(x = type, y = log10_avgconc, fill = type)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_pvalue_manual(
    stats_type,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.6,
    size = 3
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) 
ggsave("FigExtra_StaphViralLoad.pdf", width = 4.5, height = 5)



###################################################
#                                                 #
#                 FIGURE S12C                     #
#                                                 #
###################################################

#stratify by staph aureus carriage
mrg2.sa <- mrg2b[mrg2b$Assay == 'Staphylococcus_aureus',] #first select only SA data
lengthInfections2Longer$StaphCarriagetype <- lengthInfections2Longer$StaphGranular <- lengthInfections2Longer$StaphAnytime <- '0.never' #make a new column on the survival curve graph

#make a loop that cycles through the lengthinfectionslonger and see what to put
for(i in 1:nrow(lengthInfections2Longer)){
  kid <- lengthInfections2Longer[i,]$StudyID #extract tehe study id
  
  if(kid %in% mrg2.sa$StudyID){ #check if the kid has Staph at any timepoint otherwise the loop finishes
    kid.sa <- mrg2.sa[mrg2.sa$StudyID == kid,] #subset on the studyid
    kid.virusdays <- as.numeric(lengthInfections2Longer[i,]$Onset) : as.numeric(lengthInfections2Longer[i,]$Clearance)#get the virus period
    kid.sa_list <- kid.sa[kid.sa$Day2 %in% kid.virusdays,]$InfectionLongerTrue #extract the list of SA carriage data
    
    if(max(kid.sa_list) == 1){lengthInfections2Longer$StaphAnytime[i] <- '1.StaphSometime'} #for children who have staph at any time set the column to pos
    
    if(min(kid.sa_list) == 1){lengthInfections2Longer$StaphGranular[i] <- '1.StaphThroughout'} #for children who have staph throughout indicate this to pos
    if(max(kid.sa_list) == 1 & kid.sa_list[length(kid.sa_list)] == 0){lengthInfections2Longer$StaphGranular[i] <- '2.StaphCleared'} #for children who have staph throughout indicate this to pos
    if(max(kid.sa_list) == 1 & kid.sa_list[length(kid.sa_list)] == 1){lengthInfections2Longer$StaphGranular[i] <- '3.StaphatEnd'} #for children who have staph throughout indicate this to pos
    
    if(sum(kid.sa$InfectionLongerTrue == 1) / nrow(kid.sa)>0){lengthInfections2Longer$StaphCarriagetype[i] <- '1.StaphIntermittant'}
    if(sum(kid.sa$InfectionLongerTrue == 1) / nrow(kid.sa)>0.5){lengthInfections2Longer$StaphCarriagetype[i] <- '2.StaphFrequent'}
  }
}

FigS12C <- survfit2(Surv(Duration, status) ~ StaphCarriagetype, data = lengthInfections2Longer) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +   xlim(c(0,28)) 
FigS12C 

ggsave('FigS12C_Staphgroup_viralsurvival.pdf', width = 4, height = 3.5)

#Stats:
summary(survfit2(Surv(Duration, status) ~ StaphCarriagetype, data = lengthInfections2Longer))
survdiff(Surv(Duration, status) ~ StaphCarriagetype, data = lengthInfections2Longer)



###################################################
#                                                 #
#                 FIGURE S12D                     #
#                                                 #
###################################################

#remove the never staphs
mrg_clear2_onlystaph <- mrg_clear2[mrg_clear2$StudyID %in% mrg1.sa_freq[mrg1.sa_freq$type != "0.neverStaph",]$StudyID,]
unique(mrg_clear2_onlystaph$StudyID) #check that it works

#rerun the model without the never staphs
m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = mrg_clear2_onlystaph, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - OR_adjusted$Std..Error)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + OR_adjusted$Std..Error)
names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))

ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10() 
ggsave('Fig4Extra_ORAdjustedORClearanceNextDay_virus_RemoveNeverStaph.pdf', width = 3, height = 2)

#rerun with CI and bootstrap for fixed effects
set.seed(123)
ci <- confint(m, parm = "beta_", method = "boot", nsim = 1000, quiet = TRUE)
or <- exp(fixef(m))
or
or_ci <- exp(ci)
or_ci



### only take the intermittant staphs
mrg_clear2_intermittantstaph <- mrg_clear2[mrg_clear2$StudyID %in% mrg1.sa_freq[mrg1.sa_freq$type == "1.intermittentStaph",]$StudyID,]
unique(mrg_clear2_intermittantstaph$StudyID) #check that it works

#rerun the model with only the intermittant staphs
m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = mrg_clear2_intermittantstaph, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - OR_adjusted$Std..Error)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + OR_adjusted$Std..Error)
names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))

ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10() 
ggsave('Fig4Extra_ORAdjustedORClearanceNextDay_virus_Onlyintermittent.pdf', width = 3, height = 2)

#rerun with CI and bootstrap for fixed effects
set.seed(123)
ci <- confint(m, parm = "beta_", method = "boot", nsim = 1000, quiet = TRUE)
or <- exp(fixef(m))
or
or_ci <- exp(ci)
or_ci



#remove the frequent staphs
mrg_clear2_nofrequentstaph <- mrg_clear2[mrg_clear2$StudyID %in% mrg1.sa_freq[mrg1.sa_freq$type != "2.frequentStaph",]$StudyID,]
unique(mrg_clear2_nofrequentstaph$StudyID) #check that it works

#rerun the model without the never staphs
m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = mrg_clear2_nofrequentstaph, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - OR_adjusted$Std..Error)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + OR_adjusted$Std..Error)
names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))

ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10() 
ggsave('Fig4Extra_ORAdjustedORClearanceNextDay_virus_RemoveFreqStaph.pdf', width = 3, height = 2)

#rerun with CI and bootstrap for fixed effects
set.seed(123)
ci <- confint(m, parm = "beta_", method = "boot", nsim = 1000, quiet = TRUE)
or <- exp(fixef(m))
or
or_ci <- exp(ci)
or_ci



###################################################
#                                                 #
#                 FIGURE S12E                     #
#                                                 #
###################################################

df_season <- mrg_clear2 %>%
  dplyr::mutate(month = as.integer(str_split_fixed(DateCollection_EditDate, "_", 3)[, 2])) %>%
  mutate(
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Autumn",
      TRUE ~ NA_character_
    )
  )%>%
  dplyr::mutate(winter = if_else(season == "Winter", "Winter", "Not winter")) %>%
  dplyr::arrange(month) %>%
  dplyr::mutate(month2 = factor(month,
                                levels = 1:12,
                                labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))) %>%
  dplyr::mutate(month2 = factor(month2, levels = month.abb))

df_season$DateCollection_EditDate #check what the date looks like
df_season$DaysSinceJan01 <- as.integer(format(as.Date(df_season$DateCollection_EditDate, format = "%Y_%m_%d"), "%j")) - 1 #set it to an integer
df_season$DaysCosinus <-  cos(2*pi*df_season$DaysSinceJan01/365) #change it to a cosinus

#coninuous seasonality model:
m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + DaysCosinus + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = df_season, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - OR_adjusted$Std..Error)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + OR_adjusted$Std..Error)
names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))

FigS12D <- ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10() 
FigS12D

ggsave('FigS12D_ORAdjustedORClearanceNextDay_virus_IncludingSeasonality.pdf', width = 3, height = 2)

#rerun with CI and bootstrap for fixed effects
set.seed(123)
ci <- confint(m, parm = "beta_", method = "boot", nsim = 1000, quiet = TRUE)
or <- exp(fixef(m))
or
or_ci <- exp(ci)
or_ci
