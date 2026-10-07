#####
#This makes plots of exposure
#update 9-6. Fixed the forest plot for exposure vs activities to be of standard error and not CI
#####

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

#create a function that counts samples per group
count_per_group <- function(data, group_by_cols) {
  data <- data[complete.cases(data[, group_by_cols]), ] # Drop rows with NA in any of the grouping columns
  data %>%
    group_by(across(all_of(group_by_cols))) %>%
    summarise(
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(across(all_of(group_by_cols)))
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
library(ggmosaic)
library(patchwork)

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#load in database
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))

#filter out antibiotics samples
mrg1 <- mrg <- mrg %>% 
  filter(!grepl("Yes", Antibiotics))   #remove antibiotics samples

mrg2 <- mrg1[!mrg1$Assay %in% c('16S', 'PIV1','PIV2','Spn20',
                              'Moraxella_catarrhalis', 'Haemophilus_influenzae',
                              'Spn(lytA)', 'Spn(piaB)'),] #remove 16S and the none or always detected assays

mrg2b <- mrg1 %>% #remove pathogens/child combinations where there is no infection
   group_by(StudyID, Assay) %>%
   filter(
     max(InfectionLongerTrue) > 0
   ) %>%
   ungroup()
 
 
###################################################
#                                                 #
#                 FIGURE S12A                     #
#                                                 #
###################################################

mrg_exposure <- mrg2[mrg2$ExposureEvent == 1,]

mrg_exposure$Virus <- NA
mrg_exposure$VirusNumber <- NA
mrg_exposure$SpnNumber  <- NA
mrg_exposure$Staph <- NA
mrg_exposure$Haemophilus <- NA
mrg_exposure$LytA <- NA

for(i in 1:nrow(mrg_exposure)){
  
  sample_check <- mrg_exposure$SampleID_Day[i]
  sample_check2 <- mrg[mrg$SampleID_Day == sample_check,]
  
  #fill the values
  mrg_exposure$LytA[i] <- sample_check2[sample_check2$Assay == "Spn(lytA)",'log10_avgconc']
  mrg_exposure$Staph[i] <- sample_check2[sample_check2$Assay == "Staphylococcus_aureus",'InfectionLongerTrue']
  mrg_exposure$Haemophilus[i] <- sample_check2[sample_check2$Assay == "Haemophilus_influenzae",'InfectionLongerTrue']
  mrg_exposure$VirusNumber [i] <-sum(sample_check2[sample_check2$AssayType == "Viral", 'InfectionLongerTrue'])
  mrg_exposure$SpnNumber[i] <- sum(sample_check2[grepl('Spn',sample_check2$Assay) & 
                                                   !grepl('lytA',sample_check2$Assay) & 
                                                   !grepl('piaB',sample_check2$Assay) & 
                                                   !grepl('Spn6C/D',sample_check2$Assay), 'InfectionLongerTrue'])
}
mrg_exposure$Spn <- 0 
mrg_exposure$Spn[mrg_exposure$SpnNumber >0] <- 1 #make a Spn yes.no for viral comparison
mrg_exposure$Virus <- 0 
mrg_exposure$Virus[mrg_exposure$VirusNumber >0] <- 1 #make a Spn yes.no for viral comparison


#remove from the ones where we see an acquisition for Spn one serotype to prevent it being counted against itself
mrg_exposure_Spn <- mrg_exposure_Spn0 <- mrg_exposure[grepl('Spn',mrg_exposure$Assay) & !grepl('Spn6C/D', mrg_exposure$Assay),]
mrg_exposure_Spn$SpnNumber[mrg_exposure_Spn$AcquisitionLonger == 1] <- mrg_exposure_Spn0$SpnNumber[mrg_exposure_Spn0$AcquisitionLonger == 1] - 1
mrg_exposure_Spn$Spn <- 0 
mrg_exposure_Spn$Spn[mrg_exposure_Spn$SpnNumber >0] <- 1 #make a Spn yes.no for viral comparison
table(mrg_exposure_Spn$AcquisitionLonger)
#do the same for virus
mrg_exposure_virus <- mrg_exposure_virus0 <- mrg_exposure[mrg_exposure$AssayType == "Viral",]
mrg_exposure_virus$VirusNumber[mrg_exposure_virus$AcquisitionLonger == 1] <- mrg_exposure_virus0$VirusNumber[mrg_exposure_virus0$AcquisitionLonger == 1] - 1
mrg_exposure_virus$Virus <- 0 
mrg_exposure_virus$Virus[mrg_exposure_virus$VirusNumber >0] <- 1 #make a Spn yes.no for viral comparison
mrg_exposure_bocavirus <- mrg_exposure[mrg_exposure$Assay == "bocavirus",]
mrg_exposure_HRV <- mrg_exposure[mrg_exposure$Assay == "HRV",]
mrg_exposure_staph <- mrg_exposure[mrg_exposure$Assay == "Staphylococcus_aureus",]





table1 <- table(mrg_exposure_staph$AcquisitionLonger, mrg_exposure_staph$Virus)
capture.output(table1, file = "FigS12A_cases.txt")
table2 <- prop.table(table(mrg_exposure_staph$AcquisitionLonger, mrg_exposure_staph$Virus), 2)
capture.output(table2, file = "FigS12A_proportions.txt")

m <- glmer(as.factor(AcquisitionLonger) ~ Virus + AgeYears + (1|StudyID), data = mrg_exposure_staph, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
capture.output(summary(m), file = "Fig12A_statssummary.txt")

timestamp <- format(Sys.time(), "%Y%m%d")
ggplot(mrg_exposure_staph) + 
  geom_mosaic(aes(x=product(AcquisitionLonger, Virus), fill = AcquisitionLonger)) +
  facet_grid(~Virus, space = 'free') + 
  scale_fill_grey(start = 0.8, end = 0.2) 
pdf_filename <- paste0("FigS12A_linkvirus_staphacquisition_", timestamp, ".pdf")
ggsave(pdf_filename, width = 4, height = 3)

#nr samples:
fig_data      <- mrg_exposure_staph
group_cols    <- c("AcquisitionLonger", "Virus")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig1c_samplesize.txt")



###################################################
#                                                 #
#                 FIGURE S12B                     #
#                                                 #
###################################################

#check if viral infection is linked to any bacteria or number of viruses
mrg3 <- mrg2 %>%
  group_by(StudyID, Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLonger) == 1
  ) %>%
  ungroup() #make a new dataframe and only select assays that are present

mrg2_virus <- mrg2[mrg2$AssayType == 'Viral',]

mrg_clear <- mrg3[mrg3$Assay == 'Staphylococcus_aureus',]
mrg_clear$ClearNextDay <- 0
mrg_clear$VirusNumber  <- NA

mrg_clear.split <- split(mrg_clear,  mrg_clear$StudyPath) 

for(h in 1:length(mrg_clear.split)){
  each <- data.frame(mrg_clear.split[[h]])  #get each child vs assay out
  each <- each[order(each$Day),] #make sure it is in the right order
  each2 <- each[-1,] #remove the first sample as there cannot be clearance
  each3 <- each2[each2$InfectionLongerTrue == 1 | each2$ClearanceLonger == 1,] #remove all samples for which there is no infection or clearance as there is then no possible clearance
  each4 <- each3[each3$AcquisitionLonger == 0,] #remove acquisition samples that are acquired as there is then no clearance possible
  
  for(i in 1:nrow(each4)){ #cycle through the possible samples where there could be a clearance
    
    #check if it is not the last sample and then say it will clear next day set the counter to yes
    if(i < nrow(each4)){
      if(each4$ClearanceLonger[i+1]){
        each4$ClearNextDay[i] <- 1 
      }
    }
    
    sample_check <- each4$SampleID_Day[i] #check which are the relavant sample IDs
    sample_check2 <- data.frame(mrg2_virus[mrg2_virus$SampleID_Day == sample_check,]) #get the total dataset for this for only viral assays
    each4$VirusNumber[i] <- sum(sample_check2[, 'InfectionLongerTrue'])
    #fill in the covariates of interest
    
  }
  mrg_clear.split[[h]] <- each4 #save the fixed file back in the list
}

mrg_clear2 <- do.call(rbind, mrg_clear.split)
mrg_clear2 <- mrg_clear2[mrg_clear2$ClearanceLonger == 0,] #remove the clearancelonger point
table(mrg_clear2$Assay)
table(mrg_clear2$ClearNextDay, mrg_clear2$VirusNumber)
mrg_clear2$Virus <- 0
mrg_clear2$Virus[mrg_clear2$VirusNumber>0] <- 1
table(mrg_clear2$ClearNextDay, mrg_clear2$Virus)


#add potential confounders like bacterial carriage
mrg_clear2$SpnNumber  <- NA
mrg_clear2$Haemophilus <- NA
mrg_clear2$LytA <- NA

for(i in 1:nrow(mrg_clear2)){
  
  sample_check <- mrg_clear2$SampleID_Day[i]
  sample_check2 <- mrg[mrg$SampleID_Day == sample_check,]
  
  #fill the values
  mrg_clear2$LytA[i] <- sample_check2[sample_check2$Assay == "Spn(lytA)",'log10_avgconc']
  mrg_clear2$Haemophilus[i] <- sample_check2[sample_check2$Assay == "Haemophilus_influenzae",'InfectionLongerTrue']
  mrg_clear2$SpnNumber[i] <- sum(sample_check2[grepl('Spn',sample_check2$Assay) & 
                                                 !grepl('lytA',sample_check2$Assay) & 
                                                 !grepl('piaB',sample_check2$Assay) & 
                                                 !grepl('Spn6C/D',sample_check2$Assay), 'InfectionLongerTrue'])
}
mrg_clear2$Spn <- 0 
mrg_clear2$Spn[mrg_clear2$SpnNumber >0] <- 1 #make a Spn yes.no for viral comparison


m <- glmer(as.factor(ClearNextDay) ~ Virus + AgeYears +  (1|StudyID), data = mrg_clear2, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
capture.output(summary(m), file = "FigS12B_statssummary.txt")

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)

names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))
capture.output(OR_adjusted, file = "FigS12B_statssummary_CIs.txt")

ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10() 
ggsave('FigS12B_AdjustedORClearanceNextDay_staph.pdf', width = 3, height = 2.5)



###################################################
#                                                 #
#                 FIGURE S12C                     #
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
capture.output(summary(m), file = "FigS12C_statssummary.txt")

FigS12C <- ggplot(mrg1.sa_freq, aes(x = type, y = freq, colour= type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
FigS12C

ggsave('FigS12C_freqStaphGroups.pdf', width = 5, height = 4)

#check nr per group:
mrg1.sa_freq %>%
  filter(!is.na(type), !is.na(StudyID)) %>%
  group_by(type) %>%
  summarise(
    n_unique_StudyID = n_distinct(StudyID),
    .groups = "drop"
  ) %>%
  arrange(type)



###################################################
#                                                 #
#                 FIGURE S12D                     #
#                                                 #
###################################################

mrg_clear3 <- merge(
  mrg,
  mrg1.sa_freq[, c("StudyID", "type")],
  by = "StudyID",
  all.x = TRUE
) %>% filter(AssayType == "Viral") %>%
  filter(log10_avgconc > 0)


mrg_clear3$type <- factor(mrg_clear3$type)
mrg_clear3$Assay <- factor(mrg_clear3$Assay)
mrg_clear3$StudyID <- factor(mrg_clear3$StudyID)

mrg_clear3$type <- relevel(factor(mrg_clear3$type), ref = "1.intermittentStaph") #use intermittent as ref

m <- lmer(log10_avgconc ~ type + (1 | StudyID) + (1 | Assay), data = mrg_clear3)
coefs <- coef(summary(m))
capture.output(summary(m), file = "FigS12D_statssummary.txt")

disp_pval <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < 0.001,
      scales::scientific(p, digits = 3),
      format(round(p, 3), nsmall = 3)
    )
  )
}

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
      ", P=", disp_pval(Pval)
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
ggsave("FigS12D_StaphViralLoad.pdf", width = 4.5, height = 5)



###################################################
#                                                 #
#                 FIGURE S12E                     #
#                                                 #
###################################################

lengthInfectionsLonger <- data.frame()

mrg8.split <- split(mrg2, mrg2$StudyPath)
for(i in 1:length(mrg8.split)){
  each <- mrg8.split[[i]]
  if(max(each$ClearanceLonger == 1)){ #check f there is any clearance
    clear <- which(each$ClearanceLonger == 1) -1  #get all the clearances out. Substract 1 to make sure that the last positive sample is taken
    for(clears in clear){
      for(j in clears:1){ #count back to start until there is no more infection
        if(each$InfectionLonger[j] == 1){
          start <- j
        }
        if(each$InfectionLonger[j] == 0){break} #breaks out of the loop when a 0 is seen
      }
      lengthInfectionLonger <- each$Day2[clears] - each$Day2[start] #calculate the length of infections and add
      onset = ifelse(start==1, 'UnknownOnset', 'Known') #if the first sample is positive or not
      lengthInfectionsLonger <- rbind(lengthInfectionsLonger,  #add the infection length to the total one
                                      c(each$Assay[1], each$StudyID[1], lengthInfectionLonger, 
                                        each$Day2[start], each$Day2[clears], onset))
    }
  } 
  if(each$InfectionLonger[nrow(each)] == 1){ #check if at last timepoint there is infection and then look at this one too
    clear <- nrow(each)
    for(j in nrow(each):1){ #count back to start until there is no more infection
      if(each$InfectionLonger[j] == 1){
        start <- j
      }
      if(each$InfectionLonger[j] == 0){break} #breaks out of the loop when a 0 is seen
    }
    lengthInfectionLonger <- each$Day2[clear] - each$Day2[start] #calculate the lengt of infections and add
    onset = 'UnknownEnd' #if the first sample is positive or not
    lengthInfectionsLonger <- rbind(lengthInfectionsLonger,  #add the infection length to the total one
                                    c(each$Assay[1], each$StudyID[1], lengthInfectionLonger, 
                                      each$Day2[start], each$Day2[clear], onset))
  }
}

colnames(lengthInfectionsLonger) <- c('Assay', 'StudyID','Duration', 'Onset', 'Clearance', 'DurationKnown')

#make a new table for survival analysis and add status where 0 means censored data and 1 the event occurring
lengthInfections2Longer <- lengthInfectionsLonger
lengthInfections2Longer$status <- 0
lengthInfections2Longer$status[lengthInfections2Longer$DurationKnown == 'Known'] <- 1
lengthInfections2Longer$Duration <- as.numeric(lengthInfections2Longer$Duration)

#stratify by staph aureus carriage
mrg2.sa <- mrg2b[mrg2b$Assay == 'Staphylococcus_aureus',] #first select only SA data
lengthInfections2Longer$StaphCarriagetype <- lengthInfections2Longer$StaphGranular <- lengthInfections2Longer$StaphAnytime <- '0.never' #make a new column on the survival curve graph

#make a loop that cycles through the lengthinfectionslonger and see what to put
for(i in 1:nrow(lengthInfections2Longer)){
  kid <- lengthInfections2Longer[i,]$StudyID #extract the study id
  
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

#select only viral assays:
unique_viral_assays <- unique(mrg$Assay[mrg$AssayType == "Viral"])
unique_viral_assays

lengthInfections3Longer <- lengthInfections2Longer %>%
  filter (Assay %in% unique_viral_assays)



FigS12E <- survfit2(Surv(Duration, status) ~ StaphCarriagetype, data = lengthInfections3Longer) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +   xlim(c(0,28)) 
FigS12E 



#Stats:
s <- summary(survfit2(Surv(Duration, status) ~ StaphCarriagetype, data = lengthInfections2Longer))
s
txt <- capture.output(print(s))
writeLines(txt, "Fig12E_statssummary.txt")

ss <- survdiff(Surv(Duration, status) ~ StaphCarriagetype, data = lengthInfections2Longer)
ss
txt2 <- capture.output(print(ss))
writeLines(txt2, "Fig12E_survdiff.txt")

chisq_val <- as.numeric(ss$chisq)
df_val    <- length(ss$n) - 1
p_val     <- 1 - pchisq(chisq_val, df_val)

label_text <- sprintf("χ²=%.1f; P=%.2f",
                      chisq_val, p_val)

FigS12E <- ggsurvfit(FigS12E_fit) +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) +
  add_confidence_interval() +
  xlim(c(0, 28)) +
  annotate(
    "text",
    x = 28, y = 1,
    hjust = 1, vjust = 1,
    label = label_text,
    size = 3.5
  )

FigS12E

ggsave('FigS12E_Staphgroup_viralsurvival.pdf', width = 4, height = 3.5)


###################################################
#                                                 #
#                 FIGURE S12F                     #
#                                                 #
###################################################

mrg_clear <- mrg %>%
  group_by(StudyID, Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLonger) == 1
  ) %>%
  ungroup() #make a new dataframe and only select assays that are present

#initialize the things relevant for the clearance changes
mrg_clear$Virus <- NA
mrg_clear$VirusNumber <- NA
mrg_clear$SpnNumber  <- NA
mrg_clear$Staph <- NA
mrg_clear$Haemophilus <- NA
mrg_clear$LytA <- NA
mrg_clear$ClearNextDay <- 0

mrg_clear.split <- split(mrg_clear,  mrg_clear$StudyPath) 

for(h in 1:length(mrg_clear.split)){
  each <- data.frame(mrg_clear.split[[h]])  #get each child vs assay out
  each <- each[order(each$Day),] #make sure it is in the right order
  each2 <- each[-1,] #remove the first sample as there cannot be clearance by the next day yet
  each3 <- each2[each2$InfectionLongerTrue == 1 | each2$ClearanceLonger == 1,] #remove all samples for which there is no infection or clearance as there is then no possible clearance
  each4 <- each3[each3$AcquisitionLonger == 0,] #remove acquisition samples that are acquired as there is then no clearance possible
  
  for(i in 1:nrow(each4)){ #cycle through the possible samples where there could be a clearance
    
    #check if it is not the last sample and then say it will clear next day set the counter to yes
    if(i < nrow(each4)){
      if(each4$ClearanceLonger[i+1]){
        each4$ClearNextDay[i] <- 1 
      }
    }
    
    sample_check <- each4$SampleID_Day[i] #check which are the relavant sample IDs
    sample_check2 <- data.frame(mrg[mrg$SampleID_Day == sample_check,]) #get the total dataset for this 
    
    #fill in the covariates of interest
    each4$LytA[i] <- sample_check2[sample_check2$Assay == "Spn(lytA)",'log10_avgconc']
    each4$Staph[i] <- sample_check2[sample_check2$Assay == "Staphylococcus_aureus",'InfectionLongerTrue']
    each4$Haemophilus[i] <- sample_check2[sample_check2$Assay == "Haemophilus_influenzae",'InfectionLongerTrue']
    each4$VirusNumber [i] <-sum(sample_check2[sample_check2$AssayType == "Viral", 'InfectionLongerTrue'])
    each4$SpnNumber[i] <- sum(sample_check2[grepl('Spn',sample_check2$Assay) & 
                                              !grepl('lytA',sample_check2$Assay) & 
                                              !grepl('piaB',sample_check2$Assay) & 
                                              !grepl('Spn6C/D',sample_check2$Assay), 'InfectionLongerTrue'])
  }
  mrg_clear.split[[h]] <- each4 #save the fixed file back in the list
}

mrg_clear2 <- do.call(rbind, mrg_clear.split)
mrg_clear2 <- mrg_clear2[mrg_clear2$ClearanceLonger == 0,] #remove the clearancelonger point

mrg_clear2$Spn <- 0 
mrg_clear2$Spn[mrg_clear2$SpnNumber >1] <- 1 #make a Spn yes.no for when there is more than 1 Spn at the same time
mrg_clear2$Virus <- 0 
mrg_clear2$Virus[mrg_clear2$VirusNumber >0] <- 1 #make a Spn yes.no for viral comparison
mrg_clear2$SpnNumber2 <- mrg_clear2$SpnNumber
mrg_clear2$SpnNumber2[mrg_clear2$SpnNumber>3] <- '4+'

#select only viral assays:
unique(mrg_clear2$Assay)
unique_viral_assays <- unique(mrg$Assay[mrg$AssayType == "Viral"])
unique_viral_assays

mrg_clear3 <- mrg_clear2 %>%
  filter (Assay %in% unique_viral_assays)
unique(mrg_clear3$Assay)

#
##
### 1 - remove the never staphs

mrg_clear3_onlystaph <- mrg_clear3[mrg_clear3$StudyID %in% mrg1.sa_freq[mrg1.sa_freq$type != "0.neverStaph",]$StudyID,]
unique(mrg_clear3_onlystaph$StudyID) #check that it works

#rerun the model without the never staphs
m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = mrg_clear3_onlystaph, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)

names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))
write_tsv(OR_adjusted, "Fig12F_1_statssummary.txt")

p1 <- ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10() 
p1


#
##
### 2 - only take the intermittant staphs

mrg_clear3_intermittantstaph <- mrg_clear3[mrg_clear3$StudyID %in% mrg1.sa_freq[mrg1.sa_freq$type == "1.intermittentStaph",]$StudyID,]
unique(mrg_clear3_intermittantstaph$StudyID) #check that it works

#rerun the model with only the intermittant staphs
m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = mrg_clear3_intermittantstaph, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))

summary(m)

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)

names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))
write_tsv(OR_adjusted, "Fig12F_2_statssummary.txt")

p2 <- ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10() 
p2


#
##
### 3 - remove the frequent staphs

mrg_clear3_nofrequentstaph <- mrg_clear3[mrg_clear3$StudyID %in% mrg1.sa_freq[mrg1.sa_freq$type != "2.frequentStaph",]$StudyID,]
unique(mrg_clear3_nofrequentstaph$StudyID) #check that it works

#rerun the model without the never staphs
m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = mrg_clear3_nofrequentstaph, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)

names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))
write_tsv(OR_adjusted, "Fig12F_3_statssummary.txt")

p3 <- ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10() 
p3


##Combine all

combined <- p1 / p1 / p3   # '/' stacks vertically; '|' would place side by side
combined

#save:
filename <- 'FigS12F_BacterialStability.pdf'
pdf(filename, width = 3, height = 6)
combined
dev.off()




###################################################
#                                                 #
#                 FIGURE S12G                     #
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

#select only viral assays:
unique(df_season$Assay)
unique_viral_assays <- unique(mrg$Assay[mrg$AssayType == "Viral"])
unique_viral_assays

df_season2 <- df_season %>%
  filter (Assay %in% unique_viral_assays)
unique(df_season2$Assay)


#coninuous seasonality model:
m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + DaysCosinus + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = df_season2, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)

names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))
write_tsv(OR_adjusted, "FigS12G_statssummary.txt")

FigS12G <- ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10() 
FigS12G

ggsave('FigS12G_ORAdjustedORClearanceNextDay_virus_IncludingSeasonality.pdf', width = 3, height = 2)




###############################################
save.image('figureS10_staphvirus.Rdata')



