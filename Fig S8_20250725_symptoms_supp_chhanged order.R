#####
#This script looks at the effect of virus on bacteiral densities
#the order of figures was changed
###################################
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


setwd("Z:\\Projects\\SAMSAM\\Manuscripts\\Figures\\Supplemental figure 8_symptoms/")


# Library
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(reshape2)
library(readxl)
library(lmerTest)
library(emmeans)
library(UpSetR)
library(mediation)
library(pheatmap)
library(mlVAR)
library(lmms)



#save the session info
sink("sessionInfo.txt")
sessionInfo()
sink()




#load in database
mrg <- data.frame(read_xlsx('Z:\\Projects\\SAMSAM\\IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20250722.xlsx'))


#######################################################
#
#Figure S8A
#
##############################################################
#fix the symptoms
mrg$SnottyNose[is.na(mrg$SnottyNose)] <- 'NA'
mrg$SnottyNose_Level[is.na(mrg$SnottyNose_Level)] <- 'NA'
mrg$RespiratoryComplaints[is.na(mrg$RespiratoryComplaints)] <- 'NA'
table(mrg$SnottyNose, mrg$SnottyNose_Level)
mrg$SnottyNose_Level[mrg$SnottyNose == 'No'] <- 0
mrg$SnottyNose[mrg$SnottyNose_Level == '0'] <-  'No'
table(mrg$SnottyNose, mrg$SnottyNose_Level)


#make a total symptom score which is just the sum of the individual ones
mrg$SneezeYN <- NA
mrg$SneezeYN[mrg$Sneeze == 'No'] <- 0
mrg$SneezeYN[mrg$Sneeze == 'Yes'] <- 1
mrg$SnottyNoseYN <- NA
mrg$SnottyNoseYN[mrg$SnottyNose == 'No'] <- 0
mrg$SnottyNoseYN[mrg$SnottyNose == 'Yes'] <- 1
mrg$CoughYN <- NA
mrg$CoughYN[mrg$Cough == 'No'] <- 0
mrg$CoughYN[mrg$Cough == 'Yes'] <- 1
mrg$FeverYN <- NA
mrg$FeverYN[mrg$Fever == 'No'] <- 0
mrg$FeverYN[mrg$Fever == 'Yes'] <- 1
mrg$ThroatYN <- NA
mrg$ThroatYN[mrg$Throat == 'No'] <- 0
mrg$ThroatYN[mrg$Throat == 'Yes'] <- 1
mrg$EarYN <- NA
mrg$EarYN[mrg$Ear == 'No'] <- 0
mrg$EarYN[mrg$Ear == 'Yes'] <- 1
mrg$SnottyNose_Level <- as.numeric(mrg$SnottyNose_Level)
mrg$TSS <- mrg$SnottyNose_Level + mrg$SneezeYN +mrg$CoughYN +
  mrg$FeverYN+mrg$EarYN
mrg$RespiratoryComplaintsYN <- NA
mrg$RespiratoryComplaintsYN[mrg$RespiratoryComplaints == 'No'] <- 0
mrg$RespiratoryComplaintsYN[mrg$RespiratoryComplaints == 'Yes'] <- 1

table(mrg$RespiratoryComplaintsYN, mrg$RespiratoryComplaints)
table(mrg$SnottyNoseYN)
table(mrg$CoughYN)
table(mrg$FeverYN)
table(mrg$EarYN)





###################################
#make a plot for HRV over time select only HRV positive children
###########################################
mrg_HRV <- mrg %>% 
  filter(Assay == 'HRV')%>% 
  group_by(StudyID) %>%
  filter(
    max(AcquisitionLonger) > 0
  ) %>%
  #filter(                                  #these lines of code remove the children who have multiple viruses
  #  max(InfectionNumberLonger) == 1
  #) %>%
  ungroup()

mrg_HRV.split <- split(mrg_HRV, mrg_HRV$StudyID)
for(i in 1:length(mrg_HRV.split)){
  each <- mrg_HRV.split[[i]] #cycle through the childre
  dates <- each[each$AcquisitionLonger == 1,]$Day2 #select the day of viral acquisition
  each$DayToVirus <- each$Day2 - dates[1] #creat a new time column that normalizes to the FIRST time a virus is acquired
  mrg_HRV.split[[i]] <- each #save back in the list
}
mrg_HRV2 <- do.call(rbind, mrg_HRV.split)


mrg_HRV2.split <- split(mrg_HRV2, mrg_HRV2$AgeYears)
for(i in 1:length(mrg_HRV2.split)){
  mrg_HRV2.age <- mrg_HRV2.split[[i]]
  mrg_HRV2.age <- mrg_HRV2.age[!is.na(mrg_HRV2.age$RespiratoryComplaintsYN) &
                                 !is.na(mrg_HRV2.age$CoughYN) &
                                 !is.na(mrg_HRV2.age$SnottyNoseYN) ,]
HRV_complaints_age <- data.frame('Any'=prop.table(table(mrg_HRV2.age$DayToVirus, mrg_HRV2.age$RespiratoryComplaints), margin =1 )[,2],
                                'Cough' = prop.table(table(mrg_HRV2.age$DayToVirus, mrg_HRV2.age$Cough), margin =1 )[,2],
                             'SnottyNose' = prop.table(table(mrg_HRV2.age$DayToVirus, mrg_HRV2.age$SnottyNose), margin =1 )[,2]

)
HRV_complaints_age$Day <- as.numeric(rownames(HRV_complaints_age))
HRV_complaints_age$Age <- i

if(i == 1){
  HRV_complaints_age2 <- HRV_complaints_age
} else {
  HRV_complaints_age2 <- rbind(HRV_complaints_age2, HRV_complaints_age)
}
}

HRV_complaints2_age <- melt(HRV_complaints_age2, id = c('Day', 'Age'))


ggplot(HRV_complaints2_age, aes(x=Day, y = value*100, colour = as.factor(Age))) + geom_line(linewidth = 1) + 
  xlim(c(-3,10)) + 
  ylab('Percentage of children with symptoms') + 
  xlab('Days since HRV acquisition') + 
  geom_vline(xintercept = -1, linetype = 'dashed') + 
  facet_grid(.~variable)
timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("FigS8A_symptoms_HRV_over_time_byage",
                       timestamp, ".pdf")
ggsave(pdf_filename, width = 8,height = 3)




#####################################################
#
#Figure S8C
#
#################################################

mrg2 <- mrg[mrg$Assay == 'Spn(lytA)' | 
               grepl("Haemophilus_influenzae", mrg$Assay)| 
               grepl("aureus", mrg$Assay)|
               grepl("Moraxella_catarrhalis", mrg$Assay) |
              mrg$Assay == 'Streptococcus_pyogenes' |
              mrg$Assay == 'Mycoplasma_pneumoniae' | 
              mrg$AssayType == 'Viral' ,] #select only assays of interest and remove the abx samples
 
mrg3 <- mrg2[mrg2$Antibiotics != 'Yes' & 
              mrg2$Assay != 'PIV1' & 
               mrg2$Assay != 'PIV2',]


mrg_all <- mrg3 %>% 
  filter(AssayType == 'Viral')


mrg_all2 <- data.frame(mrg_all[,c(1,69,97)])
mrg_all3 <- mrg_all2 %>%
  pivot_wider(names_from = Assay, values_from = InfectionLongerTrue, values_fill = 0)

#make the name that says whether there is infection or not and for which virus
mrg_all3$InfectionType = '0.None'
mrg_all3$InfectionType[mrg_all3$enterovirus == 1] <- 'enterovirus'
mrg_all3$InfectionType[mrg_all3$HRV == 1] <- 'HRV'
mrg_all3$InfectionType[mrg_all3$influenza_A == 1] <- 'influenza_A'
mrg_all3$InfectionType[mrg_all3$bocavirus == 1] <- 'bocavirus'
mrg_all3$InfectionType[mrg_all3$PIV4 == 1] <- 'PIV4'
mrg_all3$InfectionType[mrg_all3$influenza_B == 1] <- 'influenza_B'
mrg_all3$InfectionType[mrg_all3$'sars-cov-2' == 1] <- 'sars-cov-2'
mrg_all3$InfectionType[mrg_all3$coronavirus_NL63 == 1] <- 'coronavirus_NL63'
mrg_all3$InfectionType[mrg_all3$coronavirus_OC43 == 1] <- 'coronavirus_OC43'
mrg_all3$InfectionType[mrg_all3$adenovirus == 1] <- 'adenovirus'
mrg_all3$InfectionType[mrg_all3$PIV3 == 1] <- 'PIV3'
mrg_all3$InfectionType[mrg_all3$coronavirus_HKU == 1] <- 'coronavirus_HKU'
mrg_all3$InfectionType[mrg_all3$coronavirus_229E == 1] <- 'coronavirus_229E'
mrg_all3$InfectionType[rowSums(mrg_all3[,2:15]) > 1] <- 'z.multiple_viruses'
table(mrg_all3$InfectionType)

#make also one for viral family to increase sample size of subtypes
mrg_all3$InfectionTypeBroad <- mrg_all3$InfectionType
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'sars-cov-2'] <- 'coronavirus'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'coronavirus_OC43'] <- 'coronavirus'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'coronavirus_HKU'] <- 'coronavirus'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'coronavirus_229E'] <- 'coronavirus'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'coronavirus_NL63'] <- 'coronavirus'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'PIV3'] <- 'PIV'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'PIV4'] <- 'PIV'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'influenza_A'] <- 'influenza'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'influenza_B'] <- 'influenza'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'HRV'] <- 'entero_rhino_virus'
mrg_all3$InfectionTypeBroad[mrg_all3$InfectionType == 'enterovirus'] <- 'entero_rhino_virus'
table(mrg_all3$InfectionTypeBroad)

#creat also with an infection by any virus that can be used for analysis
mrg_all3$InfectionAnyVirus <- '0.None'
mrg_all3$InfectionAnyVirus[mrg_all3$InfectionType != '0.None' ] <- 'AnyVirus' 
table(mrg_all3$InfectionAnyVirus)  

#merge with the dataset and make sure for hte kid that has no virus it is set to none
mrg_all4 <- merge(mrg3, mrg_all3, by= 'SampleID_Day', all.x=T)
mrg_all4$InfectionType[is.na(mrg_all4$InfectionType)] <- '0.None'
mrg_all4$InfectionTypeBroad[is.na(mrg_all4$InfectionTypeBroad)] <- '0.None'
mrg_all4$InfectionAnyVirus[is.na(mrg_all4$InfectionAnyVirus)] <- '0.None'


mrg_all4_Example <- data.frame(mrg_all4[mrg_all4$Assay == 'HRV',]) #just pick one to not have duplication of samples


#cycle through the ones to have also the bacterial info in there
mrg_all4_Example$Spn  <- NA
mrg_all4_Example$Staph <- NA
mrg_all4_Example$Haemophilus <- NA
mrg_all4_Example$Moraxella <- NA
mrg_all4_Example$GAS <- NA
mrg_all4_Example$LytADens <- NA
mrg_all4_Example$StaphDens <- NA
mrg_all4_Example$HaemophilusDens <- NA
mrg_all4_Example$MoraxellaDens <- NA
mrg_all4_Example$GASDens <- NA


for(i in 1:nrow(mrg_all4_Example)){
  
  sample_check <- mrg_all4_Example$SampleID_Day[i]
  sample_check2 <- mrg[mrg$SampleID_Day == sample_check,]
  
  #fill the values
  mrg_all4_Example$Spn[i] <- sample_check2[sample_check2$Assay == "Spn(lytA)",'InfectionLongerTrue']
  mrg_all4_Example$Staph[i] <- sample_check2[sample_check2$Assay == "Staphylococcus_aureus",'InfectionLongerTrue']
  mrg_all4_Example$Haemophilus[i] <- sample_check2[sample_check2$Assay == "Haemophilus_influenzae",'InfectionLongerTrue']
  mrg_all4_Example$Moraxella[i] <- sample_check2[sample_check2$Assay == "Moraxella_catarrhalis",'InfectionLongerTrue']
  mrg_all4_Example$GAS[i] <- sample_check2[sample_check2$Assay == "Streptococcus_pyogenes",'InfectionLongerTrue']
  
  
  mrg_all4_Example$LytADens[i] <- sample_check2[sample_check2$Assay == "Spn(lytA)",'log10_avgconc']
  mrg_all4_Example$StaphDens[i] <- sample_check2[sample_check2$Assay == "Staphylococcus_aureus",'log10_avgconc']
  mrg_all4_Example$HaemophilusDens[i] <- sample_check2[sample_check2$Assay == "Haemophilus_influenzae",'log10_avgconc']
  mrg_all4_Example$MoraxellaDens[i] <- sample_check2[sample_check2$Assay == "Moraxella_catarrhalis",'log10_avgconc']
  mrg_all4_Example$GASDens[i] <- sample_check2[sample_check2$Assay == "Streptococcus_pyogenes",'log10_avgconc']
  
}


#check if any of the viruses is linked to the symptoms, PIV and adeno are removed as there are too few observations
mrg_all5 <- mrg_all4_Example[mrg_all4_Example$InfectionTypeBroad != 'adenovirus' &
                               mrg_all4_Example$InfectionTypeBroad != 'PIV',]   




mrg_all4_Example.split <- split(mrg_all4_Example, mrg_all4_Example$StudyID)

for(i in 1:length(mrg_all4_Example.split)){
  split <- mrg_all4_Example.split[[i]] 
  split$OnsetAny <- 'No'
  split$OnsetSnot <- 'No'
  split$OnsetCough <- 'No'
  split2 <- split[!is.na(split$RespiratoryComplaintsYN),]
  split2 <- split2[order(split2$Day),]
  
  for(j in 1:nrow(split2)){ #cycle through all days for any symptom
    try(if(split2$RespiratoryComplaints[j] == 'Yes'){
      if(split2$RespiratoryComplaints[j - 1] =='No' &
         split2$RespiratoryComplaints[j - 2] =='No' &
         split2$RespiratoryComplaints[j - 3] =='No' &
         split2$RespiratoryComplaints[j + 1] =='Yes'){
        split2$OnsetAny[j] <- 'Yes'
        break
      }
    })
  }
  
  for(j in 1:nrow(split2)){ #cycle through all days for cough
    try(if(split2$Cough[j] == 'Yes'){
      if(split2$Cough[j - 1] =='No' &
         split2$Cough[j - 2] =='No' &
         split2$Cough[j - 3] =='No' &
         split2$Cough[j + 1] =='Yes'){
        split2$OnsetCough[j] <- 'Yes'
        break
      }
    })
  }
  
  for(j in 1:nrow(split2)){ #cycle through all days for snotty nose
    try(if(split2$SnottyNose[j] == 'Yes'){
      if(split2$SnottyNose[j - 1] =='No' &
         split2$SnottyNose[j - 2] =='No' &
         split2$SnottyNose[j - 3] =='No' &
         split2$SnottyNose[j + 1] =='Yes'){
        split2$OnsetSnot[j] <- 'Yes'
        break
      }
    })
  }
  
  mrg_all4_Example.split[[i]] <- split2
}
mrg_symptomonset <- do.call(rbind, mrg_all4_Example.split)
table(mrg_symptomonset$OnsetAny)
table(mrg_symptomonset$OnsetCough)
table(mrg_symptomonset$OnsetSnot)



#merge with the full dataset
mrg_all50 <- merge(mrg_all4, mrg_symptomonset[,c(1, 133:135)], by = 'SampleID_Day')



###################
#filter out kids with a snotty nose onset only
mrg_all_snotonset <- mrg_all50 %>% 
  group_by(StudyID) %>%
  filter(
    max(OnsetSnot == 'Yes') > 0
  ) %>%
  ungroup()

table(mrg_all_snotonset$OnsetSnot)


#se;ect for the target and normalize for day -1 to -3 relative to snotty nose onse
#first for spn
mrg7.split <- split(mrg_all_snotonset, mrg_all_snotonset$StudyID)
target = "Spn(lytA)"
for(i in 1:length(mrg7.split)){
  each <- mrg7.split[[i]] #cycle through the childre
  dates <- each[each$OnsetSnot == 'Yes',]$Day2 #select the day of viral acquisition
  each$DayToSNOTonset <- each$Day2 - dates[1] #creat a new time column that normalizes to the FIRST time a virus is acquired
  each2 <-  each[each$Assay == target,] #select only bacterial pathogen
  base <- mean(each2$log10_avgconc[which(each2$DayToSNOTonset == -3 | each2$DayToSNOTonset == -2 | each2$DayToSNOTonset == -1)], na.rm=T) #make a baseline average conc
  #base <-each2$log10_avgconc[which(each2$DayToSNOTonset == 0)] #make a baseline average conc
  each2$FC_dens <- each2$log10_avgconc - base #make a fold change relative to the baseline
  mrg7.split[[i]] <- each2 #save back in the list
}
mrg8a <- do.call(rbind, mrg7.split)

#then for haemophilus
mrg7.split <- split(mrg_all_snotonset, mrg_all_snotonset$StudyID)
target = "Haemophilus_influenzae"
for(i in 1:length(mrg7.split)){
  each <- mrg7.split[[i]] #cycle through the childre
  dates <- each[each$OnsetSnot == 'Yes',]$Day2 #select the day of viral acquisition
  each$DayToSNOTonset <- each$Day2 - dates[1] #creat a new time column that normalizes to the FIRST time a virus is acquired
  each2 <-  each[each$Assay == target,] #select only bacterial pathogen
  base <- mean(each2$log10_avgconc[which(each2$DayToSNOTonset == -3 | each2$DayToSNOTonset == -2 | each2$DayToSNOTonset == -1)], na.rm=T) #make a baseline average conc
  #base <-each2$log10_avgconc[which(each2$DayToSNOTonset == 0)] #make a baseline average conc
  each2$FC_dens <- each2$log10_avgconc - base #make a fold change relative to the baseline
  mrg7.split[[i]] <- each2 #save back in the list
}
mrg8b <- do.call(rbind, mrg7.split)

#do stats
mrg8a1 <- mrg8a[mrg8a$DayToSNOTonset > -8 & mrg8a$DayToSNOTonset < 11, ] #filter out smaples between x days before or after viral onset
mrg8a1$DayToSNOTonset2 <- mrg8a1$DayToSNOTonset #make a new column and change the 5 days prior to onset to the baseline
mrg8a1$DayToSNOTonset2[mrg8a1$DayToSNOTonset2 %in% c(-3:-1)] <- '-1_baseline'
lm <- lmer(log10_avgconc ~ as.factor(DayToSNOTonset2) +(1|StudyID), data = mrg8a1) #do stats for all days relative to baseline
acq_res <- data.frame(lmerTest:::get_coefmat(lm)[2:nrow(lmerTest:::get_coefmat(lm)),c(1,5)]) #extract the results from the linear model

mrg8b1 <- mrg8b[mrg8b$DayToSNOTonset > -8 & mrg8b$DayToSNOTonset < 11, ] #filter out smaples between x days before or after viral onset
mrg8b1$DayToSNOTonset2 <- mrg8b1$DayToSNOTonset #make a new column and change the 5 days prior to onset to the baseline
mrg8b1$DayToSNOTonset2[mrg8b1$DayToSNOTonset2 %in% c(-3:-1)] <- '-1_baseline'
lm <- lmer(log10_avgconc ~ as.factor(DayToSNOTonset2) +(1|StudyID), data = mrg8b1) #do stats for all days relative to baseline
acq_res2 <- data.frame(lmerTest:::get_coefmat(lm)[2:nrow(lmerTest:::get_coefmat(lm)),c(1,5)]) #extract the results from the linear model

#make the dataframe for stats
acq_res$Assay <- 'Spn(lytA)'
acq_res2$Assay <- 'Haemophilus_influenzae'
acq_res$x <- gsub('as.factor[(]DayToSNOTonset2[)]', '', rownames(acq_res))
acq_res2$x <- gsub('as.factor[(]DayToSNOTonset2[)]', '', rownames(acq_res2))
acq_res3 <- rbind(acq_res, acq_res2)
acq_res3$P2 <- ''
acq_res3$P2[acq_res3$Pr...t.. *2 < 0.05] <- '*'   #correct x4 for doing 2 symptoms and 2 bacteria
acq_res3$P2[acq_res3$Pr...t.. *2< 0.01] <- '**'
acq_res3$P2[acq_res3$Pr...t.. *2 < 0.001] <- '***'



#make the same graph for cough onset
#se;ect for the target and normalize for day -1 to -3 relative to cough onset
mrg_all_coughonset <- mrg_all50 %>% 
  group_by(StudyID) %>%
  filter(
    max(OnsetCough == 'Yes') > 0
  ) %>%
  ungroup()

table(mrg_all_coughonset$OnsetSnot)


mrg7.split <- split(mrg_all_coughonset, mrg_all_coughonset$StudyID)
target = "Spn(lytA)"
for(i in 1:length(mrg7.split)){
  each <- mrg7.split[[i]] #cycle through the childre
  dates <- each[each$OnsetCough == 'Yes',]$Day2 #select the day of viral acquisition
  each$DayToCOUGHonset <- each$Day2 - dates[1] #creat a new time column that normalizes to the FIRST time a virus is acquired
  each2 <-  each[each$Assay == target,] #select only bacterial pathogen
  base <- mean(each2$log10_avgconc[which(each2$DayToCOUGHonset == -3 | each2$DayToCOUGHonset == -2 | each2$DayToCOUGHonset == -1)], na.rm=T) #make a baseline average conc
  #base <-each2$log10_avgconc[which(each2$DayToCOUGHonset == 0)] #make a baseline average conc
  each2$FC_dens <- each2$log10_avgconc - base #make a fold change relative to the baseline
  mrg7.split[[i]] <- each2 #save back in the list
}
mrg9a <- do.call(rbind, mrg7.split)



mrg7.split <- split(mrg_all_coughonset, mrg_all_coughonset$StudyID)
target = "Haemophilus_influenzae"
for(i in 1:length(mrg7.split)){
  each <- mrg7.split[[i]] #cycle through the childre
  dates <- each[each$OnsetCough == 'Yes',]$Day2 #select the day of viral acquisition
  each$DayToCOUGHonset <- each$Day2 - dates[1] #creat a new time column that normalizes to the FIRST time a virus is acquired
  each2 <-  each[each$Assay == target,] #select only bacterial pathogen
  base <- mean(each2$log10_avgconc[which(each2$DayToCOUGHonset == -3 | each2$DayToCOUGHonset == -2 | each2$DayToCOUGHonset == -1)], na.rm=T) #make a baseline average conc
  #base <-each2$log10_avgconc[which(each2$DayToCOUGHonset == 0)] #make a baseline average conc
  each2$FC_dens <- each2$log10_avgconc - base #make a fold change relative to the baseline
  mrg7.split[[i]] <- each2 #save back in the list
}
mrg9b <- do.call(rbind, mrg7.split)


#do stats
mrg9a1 <- mrg9a[mrg9a$DayToCOUGHonset > -8 & mrg9a$DayToCOUGHonset < 11, ] #filter out smaples between x days before or after viral onset
mrg9a1$DayToCOUGHonset2 <- mrg9a1$DayToCOUGHonset #make a new column and change the 5 days prior to onset to the baseline
mrg9a1$DayToCOUGHonset2[mrg9a1$DayToCOUGHonset2 %in% c(-3:-1)] <- '-1_baseline'
lm <- lmer(log10_avgconc ~ as.factor(DayToCOUGHonset2) +(1|StudyID), data = mrg9a1) #do stats for all days relative to baseline
acq_res <- data.frame(lmerTest:::get_coefmat(lm)[2:nrow(lmerTest:::get_coefmat(lm)),c(1,5)]) #extract the results from the linear model

mrg9b1 <- mrg9b[mrg9b$DayToCOUGHonset > -8 & mrg9b$DayToCOUGHonset < 11, ] #filter out smaples between x days before or after viral onset
mrg9b1$DayToCOUGHonset2 <- mrg9b1$DayToCOUGHonset #make a new column and change the 5 days prior to onset to the baseline
mrg9b1$DayToCOUGHonset2[mrg9b1$DayToCOUGHonset2 %in% c(-3:-1)] <- '-1_baseline'
lm <- lmer(log10_avgconc ~ as.factor(DayToCOUGHonset2) +(1|StudyID), data = mrg9b1) #do stats for all days relative to baseline
acq_res2 <- data.frame(lmerTest:::get_coefmat(lm)[2:nrow(lmerTest:::get_coefmat(lm)),c(1,5)]) #extract the results from the linear model

#make the dataframe for stats
acq_res$Assay <- 'Spn(lytA)'
acq_res2$Assay <- 'Haemophilus_influenzae'
acq_res$x <- gsub('as.factor[(]DayToCOUGHonset2[)]', '', rownames(acq_res))
acq_res2$x <- gsub('as.factor[(]DayToCOUGHonset2[)]', '', rownames(acq_res2))
acq_res4 <- rbind(acq_res, acq_res2)
acq_res4$P2 <- ''
acq_res4$P2[acq_res4$Pr...t.. *2 < 0.05] <- '*'   #correct x4 for doing 2 symptoms and 2 bacteria
acq_res4$P2[acq_res4$Pr...t.. *2< 0.01] <- '**'
acq_res4$P2[acq_res4$Pr...t.. *2 < 0.001] <- '***'


mrg9 <- rbind(mrg9a, mrg9b)


###########################
#associate symptom onset with frequency of virus
table(mrg9$DayToCOUGHonset, mrg9$InfectionAnyVirus)
coughonset_virus <- data.frame(prop.table(table(mrg9a$DayToCOUGHonset, mrg9a$InfectionAnyVirus), margin = 1)[,2])
snotonset_virus <- data.frame(prop.table(table(mrg8a$DayToSNOTonset, mrg8a$InfectionAnyVirus), margin = 1)[,2])
coughonset_virus$Symptom <- 'Cough'
snotonset_virus$Symptom <- 'Snot'
coughonset_virus$Day <- as.numeric(rownames(coughonset_virus))
snotonset_virus$Day <- as.numeric(rownames(snotonset_virus))
colnames(snotonset_virus)[1] <- colnames(coughonset_virus)[1] <- 'Freq'
onset_virus <- rbind(coughonset_virus, snotonset_virus)

ggplot(onset_virus, aes(x=Day, y = Freq*100, colour = Symptom)) + geom_line() + 
  geom_point()+
  xlim(c(-10,10))  +
  ylab('Percentage of children with virus infection') + 
  xlab('Days since symptom onset') + 
  geom_vline(xintercept = 0, linetype = 'dashed')

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("FigS8C_virus_relative_to_symptoms_time",
                       timestamp, ".pdf")
ggsave(pdf_filename, width = 3,height = 2)





###########################################################################
#
#Figure S8E
#
###########################################################################
#do VAR to see how lyta,virus and symptoms interact
mycols <- c("StudyID",
            "Day",
            "LytADens",
            'HaemophilusDens',
            'InfectionAnyVirus',
            'SnottyNose_Level',
            'SnottyNoseYN', 
            'TSS'
)


mrg_mlvar <- mrg_all4_Example %>% dplyr::select(all_of(mycols)) 
mrg_mlvar$Day <- as.numeric(gsub('D', '', mrg_mlvar$Day)) 
mrg_mlvar <- mrg_mlvar[!is.na(mrg_mlvar$Day) & !is.na(mrg_mlvar$SnottyNose_Level),]
mrg_mlvar$InfectionAnyVirus <- as.numeric(as.factor(mrg_mlvar$InfectionAnyVirus))-1




VAR2 <- mlVAR(mrg_mlvar, vars = c('HaemophilusDens', 'SnottyNose_Level', 'InfectionAnyVirus'), 
              idvar = 'StudyID', lags = 1, contemporaneous = 'correlated',
              estimator = 'lmer',)


pdf('FigS8E_VAR_Haemophilus_snot_virus_lag1.pdf', width = 5, height = 3)
par(mfrow=c(1,2))
plot(VAR2, posCol = 'Red', negCol ='Blue', title = 'Temporal HI - lag 1', label.cex = 1, edge.width=10,
     edge.labels =T, label.scale = F, edge.label.color = 'black', mode = 'direct')
dev.off()





#######################################################################
#
#Figure S8B
#
#######################################################################
bacteria <- mrg_all4[mrg_all4$Assay %in% c('Haemophilus_influenzae',
                                           'Spn(lytA)'),]

bacteria_pos <- bacteria[bacteria$InfectionLongerTrue ==1     , ]

#make graph for cough
bacteria_cough <- bacteria_pos[!is.na(bacteria_pos$CoughYN),]
pl6b <- ggplot(bacteria_cough, aes(x=as.factor(CoughYN), y = log10_avgconc))  + 
  geom_violin(aes(fill=as.factor(CoughYN))) + 
  geom_boxplot(width = 0.2) + 
  facet_grid(Assay~InfectionAnyVirus) 


pdf('FigS8B_density and virus vs symptoms.pdf', width = 5, height = 4)
print(pl6b)
dev.off()



Spn_cough <- bacteria_cough[bacteria_cough$Assay == 'Spn(lytA)',]
HI_cough <- bacteria_cough[bacteria_cough$Assay != 'Spn(lytA)',]
m <- glmer(as.factor(CoughYN) ~ log10_avgconc   + InfectionAnyVirus+ AgeYears +(1|StudyID), 
           data = Spn_cough, family ='binomial', control = glmerControl(optimizer = "bobyqa"))
summary(m)
m <- glmer(as.factor(CoughYN) ~ log10_avgconc   + InfectionAnyVirus+ AgeYears +  (1|StudyID), 
           data = HI_cough, family ='binomial', control = glmerControl(optimizer = "bobyqa"))
summary(m)



#######################################################################
#
#Figure S8D
#
#######################################################################
pl7a <- ggplot(mrg9, aes(x=DayToCOUGHonset, y=FC_dens, colour = Assay)) + 
  xlim(c(-7,10)) + 
  geom_hline(yintercept=0, linetype='dashed')  + 
  geom_vline(xintercept=0, linetype='dashed')  +
  #stat_summary(fun.data = mean_se,  geom = "errorbar", width=0.2) +
  stat_summary(fun.data = mean_cl_normal,  
               geom = "errorbar", width=0.2)+
  stat_summary(fun.y = mean, geom = "point",  size =4) +
  
  ylab('Density normalized to onset (log10)') + xlab('Day to snotty nose onset') +
  facet_grid(.~Assay) + 
  geom_text(data= acq_res4, aes(label=P2, x=as.numeric(x)), y=2.8, colour = 'black')

print(pl7a)
ggsave('FigS8D_bacterial_density_relative_to_cough.pdf', width = 8, height = 4)





###########################################################################33
save.image('FigS8_symptoms.Rdata')

