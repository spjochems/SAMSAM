#####
#This script looks at the effect of virus on bacteiral densities

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

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\Fig5\\")

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
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#load in database
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))
mrg <- mrg[mrg$Antibiotics != 'Yes',] #remove the antibiotics samples as it may affect results

mrg$Rhinorrhea <- mrg$SnottyNose
mrg$Rhinorrhea_Level <- mrg$SnottyNose_Level
mrg$RhinorrheaYN <- mrg$SnottyNoseYN 

#add seasonality:
mrg <- mrg %>%
  #filter(Day == "D01") %>%
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
                                           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))) 
#dplyr::mutate(month2 = factor(month2, levels = month.abb)) 

#enocde the season as a cosinus
mrg$DateCollection_EditDate #check what the date looks like
mrg$DaysSinceJan01 <- as.integer(format(as.Date(mrg$DateCollection_EditDate, format = "%Y_%m_%d"), "%j")) - 1 #set it to an integer
mrg$DaysCosinus <-  cos(2*pi*mrg$DaysSinceJan01/365) #change it to a cosinus


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


###################################################
#                                                 #
#                 FIGURE 5A                       #
#                                                 #
###################################################

#fix the symptoms
mrg$Rhinorrhea[is.na(mrg$Rhinorrhea)] <- 'NA'
mrg$Rhinorrhea_Level[is.na(mrg$Rhinorrhea_Level)] <- 'NA'
mrg$RespiratoryComplaints[is.na(mrg$RespiratoryComplaints)] <- 'NA'
table(mrg$Rhinorrhea, mrg$Rhinorrhea_Level)
mrg$Rhinorrhea_Level[mrg$Rhinorrhea == 'No'] <- 0
mrg$Rhinorrhea[mrg$Rhinorrhea_Level == '0'] <-  'No'
table(mrg$Rhinorrhea, mrg$Rhinorrhea_Level)


example <- mrg[mrg$Assay == '16S',]
upsetList <- list(Rhinorrhea = which(example$Rhinorrhea == 'Yes'),
                  Sneeze = which(example$Sneeze == 'Yes'),
                  Cough = which(example$Cough == 'Yes'),
                  Fever = which(example$Fever == 'Yes'),
                  Ear = which(example$Ear == 'Yes'))
Fig5A <- upset(fromList(upsetList), order.by = "freq", nsets=5, point.size = 2.5, line.size = 0)
Fig5A

table(example$RespiratoryComplaints)
table(example$Rhinorrhea)
table(example$Cough)

sum(table(example$Sneeze, example$StudyID)[3,]>1) #how many children sneeze

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig5A_upsetplot_symptoms_",
                       timestamp, ".pdf")
pdf(pdf_filename, width = 5,height = 5)
Fig5A
dev.off()


#make a total symptom score which is just the sum of the individual ones
mrg$SneezeYN <- NA
mrg$SneezeYN[mrg$Sneeze == 'No'] <- 0
mrg$SneezeYN[mrg$Sneeze == 'Yes'] <- 1
mrg$RhinorrheaYN <- NA
mrg$RhinorrheaYN[mrg$Rhinorrhea == 'No'] <- 0
mrg$RhinorrheaYN[mrg$Rhinorrhea == 'Yes'] <- 1
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
mrg$Rhinorrhea_Level <- as.numeric(mrg$Rhinorrhea_Level)
mrg$TSS <- mrg$Rhinorrhea_Level + mrg$SneezeYN +mrg$CoughYN +
  mrg$FeverYN+mrg$EarYN
mrg$RespiratoryComplaintsYN <- NA
mrg$RespiratoryComplaintsYN[mrg$RespiratoryComplaints == 'No'] <- 0
mrg$RespiratoryComplaintsYN[mrg$RespiratoryComplaints == 'Yes'] <- 1

table(mrg$RespiratoryComplaintsYN, mrg$RespiratoryComplaints)
table(mrg$RhinorrheaYN)
table(mrg$CoughYN)
table(mrg$FeverYN)
table(mrg$EarYN)




###################################################
#                                                 #
#                 FIGURE 5B                       #
#                                                 #
###################################################

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

#check for age relationship to symptoms
age_complaints <- data.frame('Any'=prop.table(table(mrg_all5$AgeYears, mrg_all5$RespiratoryComplaints), margin =1 )[,3],
                             'Sneeze' = prop.table(table(mrg_all5$AgeYears, mrg_all5$Sneeze), margin =1 )[,3],
                             'Cough' = prop.table(table(mrg_all5$AgeYears, mrg_all5$Cough), margin =1 )[,3],
                             'Rhinorrhea' = prop.table(table(mrg_all5$AgeYears, mrg_all5$Rhinorrhea), margin =1 )[,3],
                             'Ear' = prop.table(table(mrg_all5$AgeYears, mrg_all5$Ear), margin =1 )[,3],
                             'Fever' = prop.table(table(mrg_all5$AgeYears, mrg_all5$Fever), margin =1 )[,3]
)

age_complaints$AgeYears <- as.numeric(rownames(age_complaints))
age_complaints2 <- melt(age_complaints, id = 'AgeYears')

#age_complaints2$variable <- factor(age_complaints2$variable,
#  levels = c("Any", "Rhinorrhea","Cough", "Sneeze",  "Ear", "Fever"))

Fig5B <- ggplot(age_complaints2, aes(x=AgeYears, y = value*100, colour = variable)) + geom_line(linewidth = 2) + 
  ylab('Percentage of questionnaires with symptoms') + 
  xlab('Age in Years') 
Fig5B

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig5B_symptoms_by_age",
                       timestamp, ".pdf")
ggsave(pdf_filename, width = 5,height = 4)


m <- glmer(as.factor(RhinorrheaYN) ~  AgeYears +  (1|StudyID), data = mrg_all5, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
capture.output(summary(m), file = "Fig5B_statssummary.txt")

#nr samples:
samplesize <- mrg_all5 %>%
  filter(!is.na(AgeYears)) %>%
  group_by(AgeYears) %>%
  summarise(
    Any        = sum(RespiratoryComplaints == "Yes", na.rm = TRUE),
    Sneeze     = sum(Sneeze == "Yes", na.rm = TRUE),
    Cough      = sum(Cough == "Yes", na.rm = TRUE),
    Rhinorrhea = sum(Rhinorrhea == "Yes", na.rm = TRUE),
    Ear        = sum(Ear == "Yes", na.rm = TRUE),
    Fever      = sum(Fever == "Yes", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(AgeYears)
samplesize
capture.output(samplesize, file = "Fig5B_samplesize.txt")




###################################################
#                                                 #
#                 FIGURE 5C                       #
#                                                 #
###################################################

m <- glmer(as.factor(RespiratoryComplaintsYN) ~ InfectionTypeBroad + AgeYears +  (1|StudyID), data = mrg_all5, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

m1 <- glmer(as.factor(RespiratoryComplaintsYN) ~ InfectionAnyVirus + AgeYears +  (1|StudyID), data = mrg_all5, 
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1)
Any_virus <- rbind(summary(m)[[10]], summary(m1)[[10]])


#link to snotty nose
m <- glmer(as.factor(RhinorrheaYN) ~ InfectionTypeBroad + AgeYears +  (1|StudyID), data = mrg_all5, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
m1 <- glmer(as.factor(RhinorrheaYN) ~ InfectionAnyVirus + AgeYears +  (1|StudyID), data = mrg_all5, 
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1)
Snotty_virus <- rbind(summary(m)[[10]], summary(m1)[[10]])


#link to cough
m <- glmer(as.factor(CoughYN) ~ InfectionTypeBroad + AgeYears +  (1|StudyID), data = mrg_all5, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
m1 <- glmer(as.factor(CoughYN) ~ InfectionAnyVirus + AgeYears +  (1|StudyID), data = mrg_all5, 
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1)
Cough_virus <- rbind(summary(m)[[10]], summary(m1)[[10]])


#link to sneeze
m <- glmer(as.factor(SneezeYN) ~ InfectionTypeBroad + AgeYears +  (1|StudyID), data = mrg_all5, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
m1 <- glmer(as.factor(SneezeYN) ~ InfectionAnyVirus + AgeYears +  (1|StudyID), data = mrg_all5, 
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1)
Sneeze_virus <- rbind(summary(m)[[10]], summary(m1)[[10]])


#link to ear pain
m <- glmer(as.factor(EarYN) ~ InfectionTypeBroad + AgeYears +  (1|StudyID), data = mrg_all5, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
m1 <- glmer(as.factor(EarYN) ~ InfectionAnyVirus + AgeYears +  (1|StudyID), data = mrg_all5, 
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1)
Ear_virus <- rbind(summary(m)[[10]], summary(m1)[[10]])


all_virus_p <- data.frame(Any = Any_virus[c(1:6,9),4], 
                          Rhinorrhea = Snotty_virus[c(1:6,9),4],
                          Cough = Cough_virus[c(1:6,9),4],
                          Sneeze = Sneeze_virus[c(1:6,9),4],
                          Ear = Ear_virus[c(1:6,9),4])
all_virus_p2 <- all_virus_p

all_virus_p2[1,] <- 1
rownames(all_virus_p2)[1] <- '0.None'
rownames(all_virus_p2) <- gsub('InfectionTypeBroad', '', rownames(all_virus_p2))
rownames(all_virus_p2) <- gsub('InfectionAnyVirus', '', rownames(all_virus_p2))
all_virus_p3 <- all_virus_p2
all_virus_p3[all_virus_p2*5 > 0.05] <- '' 
all_virus_p3[all_virus_p2*5 < 0.05] <- '*' 
all_virus_p3[all_virus_p2*5 < 0.01] <- '**' 
all_virus_p3[all_virus_p2*5 < 0.001] <- '***' 



virus_complaints <- data.frame('Any'=prop.table(table(mrg_all5$InfectionTypeBroad, mrg_all5$RespiratoryComplaints), margin =1 )[,3],
                               'Sneeze' = prop.table(table(mrg_all5$InfectionTypeBroad, mrg_all5$Sneeze), margin =1 )[,3],
                               'Cough' = prop.table(table(mrg_all5$InfectionTypeBroad, mrg_all5$Cough), margin =1 )[,3],
                               'Rhinorrhea' = prop.table(table(mrg_all5$InfectionTypeBroad, mrg_all5$Rhinorrhea), margin =1 )[,3],
                               'Ear' = prop.table(table(mrg_all5$InfectionTypeBroad, mrg_all5$Ear), margin =1 )[,3]
)


virus_complaints2 <- data.frame('Any'=prop.table(table(mrg_all5$InfectionAnyVirus, mrg_all5$RespiratoryComplaints), margin =1 )[,3],
                                'Sneeze' = prop.table(table(mrg_all5$InfectionAnyVirus, mrg_all5$Sneeze), margin =1 )[,3],
                                'Cough' = prop.table(table(mrg_all5$InfectionAnyVirus, mrg_all5$Cough), margin =1 )[,3],
                                'Rhinorrhea' = prop.table(table(mrg_all5$InfectionAnyVirus, mrg_all5$Rhinorrhea), margin =1 )[,3],
                                'Ear' = prop.table(table(mrg_all5$InfectionAnyVirus, mrg_all5$Ear), margin =1 )[,3]
)


virus_complaints3 <- rbind(virus_complaints, virus_complaints2[2,])
virus_complaints3 <- virus_complaints3[,match(colnames(all_virus_p3), colnames(virus_complaints3))]
Fig5C <- pheatmap::pheatmap(virus_complaints3, display_numbers = all_virus_p3, cluster_rows = F, cluster_cols = F,
                            fontsize_number = 20)
Fig5C

pdf('Fig5C_heatmap_virus_vs_symptoms.pdf', width = 4, height = 4,onefile=T)
Fig5C
dev.off()


#sample size:
virus_symptom_counts_type <- mrg_all5 %>%
  filter(!is.na(InfectionTypeBroad)) %>%
  group_by(InfectionTypeBroad) %>%
  summarise(
    n_total    = n(),
    Any        = sum(RespiratoryComplaints == "Yes", na.rm = TRUE),
    Sneeze     = sum(Sneeze == "Yes", na.rm = TRUE),
    Cough      = sum(Cough == "Yes", na.rm = TRUE),
    Rhinorrhea = sum(Rhinorrhea == "Yes", na.rm = TRUE),
    Ear        = sum(Ear == "Yes", na.rm = TRUE),
    Fever      = sum(Fever == "Yes", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(InfectionTypeBroad)
virus_symptom_counts_type





###################################################
#                                                 #
#                 FIGURE 5D                       #
#                                                 #
###################################################

mrg_infected_virus <- mrg_all[mrg_all$InfectionLongerTrue == 1,]

m <- lmer(log10_avgconc ~as.factor(RespiratoryComplaintsYN) + AgeYears +  (1|Assay) +  (1|StudyID), data = mrg_infected_virus)
summary(m)
capture.output(table2, file = "Fig5D_statssummary.txt")

mrg_infected_virus2 <- mrg_infected_virus[!is.na(mrg_infected_virus$RespiratoryComplaintsYN),]
mrg_infected_virus2$RespiratoryComplaintsYN2 <- ifelse(
  mrg_infected_virus2$RespiratoryComplaintsYN == 0, "No",
  ifelse(mrg_infected_virus2$RespiratoryComplaintsYN == 1, "Yes", NA)
)

Fig5D <- ggplot(mrg_infected_virus2, aes(x=as.factor(RespiratoryComplaintsYN2), y=log10_avgconc)) + 
  geom_jitter(aes(colour = Assay), height = 0, width = 0.3)+
  geom_boxplot(outlier.colour = NA,  alpha= 0.5) +
  ylab('Log10 (genome copies/mL +1)') + 
  xlab('Any symptom') 
Fig5D

pdf_filename <- paste0("Fig5D_symptoms_related_to_density",
                       timestamp, ".pdf")
ggsave(pdf_filename, width = 5,height = 4)




###################################################
#                                                 #
#                 FIGURE 5E                       #
#                                                 #
###################################################

mrg_HRV <- mrg %>% 
  filter(Assay == 'HRV')%>% 
  group_by(StudyID) %>%
  filter(
    max(AcquisitionLonger) > 0
  ) %>%
  ungroup()

mrg_HRV.split <- split(mrg_HRV, mrg_HRV$StudyID)
for(i in 1:length(mrg_HRV.split)){
  each <- mrg_HRV.split[[i]] #cycle through the childre
  dates <- each[each$AcquisitionLonger == 1,]$Day2 #select the day of viral acquisition
  each$DayToVirus <- each$Day2 - dates[1] #creat a new time column that normalizes to the FIRST time a virus is acquired
  mrg_HRV.split[[i]] <- each #save back in the list
}
mrg_HRV2 <- do.call(rbind, mrg_HRV.split)


HRV_complaints <- data.frame('Any'=prop.table(table(mrg_HRV2$DayToVirus, mrg_HRV2$RespiratoryComplaints), margin =1 )[,3],
                             'Sneeze' = prop.table(table(mrg_HRV2$DayToVirus, mrg_HRV2$Sneeze), margin =1 )[,3],
                             'Cough' = prop.table(table(mrg_HRV2$DayToVirus, mrg_HRV2$Cough), margin =1 )[,3],
                             'Rhinorrhea' = prop.table(table(mrg_HRV2$DayToVirus, mrg_HRV2$Rhinorrhea), margin =1 )[,3],
                             'Ear' = prop.table(table(mrg_HRV2$DayToVirus, mrg_HRV2$Ear), margin =1 )[,3]
)
HRV_complaints$Day <- as.numeric(rownames(HRV_complaints))
HRV_complaints2 <- melt(HRV_complaints, id = 'Day')


#do stats per time
mrg_HRV2$DayToVirus2 <- mrg_HRV2$DayToVirus #make a new column and change the 5 days prior to onset to the baseline
mrg_HRV2$DayToVirus2[mrg_HRV2$DayToVirus2 %in% c(-3:-1)] <- '-1_baseline'
mrg_HRV3 <- mrg_HRV2[mrg_HRV2$DayToVirus > -4 & mrg_HRV2$DayToVirus < 11,] #only select the good timepoints


#make the models based on time for any symptoms, snotty nose and couching
m <- glmer(as.factor(RhinorrheaYN) ~ as.factor(DayToVirus2) + AgeYears + DaysCosinus + (1|StudyID), data = mrg_HRV3, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
model1 <- summary(m)[[10]]
model1a <- model1[2:nrow(model1), 1]
model1b <- model1[2:nrow(model1), 4]


m <- glmer(as.factor(CoughYN) ~ as.factor(DayToVirus2) + AgeYears + DaysCosinus + (1|StudyID), data = mrg_HRV3, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
model2 <- summary(m)[[10]]
model2a <- model2[2:nrow(model2), 1]
model2b <- model2[2:nrow(model2), 4]


m <- glmer(as.factor(RespiratoryComplaintsYN) ~ as.factor(DayToVirus2) + AgeYears + DaysCosinus + (1|StudyID), data = mrg_HRV3, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
model3 <- summary(m)[[10]]
model3a <- model3[2:nrow(model3), 1]
model3b <- model3[2:nrow(model3), 4]

#combine the datasets into model estimates
model_estimates <- rbind(model3a, model1a, model2a)
model_pvals <- rbind(model3b, model1b, model2b)
rownames(model_pvals) <-  rownames(model_estimates) <- c('Any', 'Rhinorrhea', 'Coughing')
colnames(model_pvals) <-  colnames(model_estimates) <- gsub('as.factor[(]DayToVirus2[)]', '', colnames(model_estimates))
model_estimates <- model_estimates[,order(as.numeric(colnames(model_estimates)))]
model_pvals <- model_pvals[,order(as.numeric(colnames(model_pvals)))]
model_pvals2 <- model_pvals
model_pvals2[model_pvals*3 > 0.05] <- ' ' #multiply by 3 for testing 3 different symptoms
model_pvals2[model_pvals*3 < 0.05] <- '*' #multiply by 3 for testing 3 different symptoms
model_pvals2[model_pvals*3 < 0.01] <- '**' 
model_pvals2[model_pvals*3 < 0.001] <- '***' 


model_pvals3 <- melt(model_pvals2)
model_pvals3$Var1 <- gsub('Coughing', 'Cough', model_pvals3$Var1)
HRV_complaints2$ID <- paste(HRV_complaints2$Day, HRV_complaints2$variable)
model_pvals3$ID <- paste(model_pvals3$Var2, model_pvals3$Var1)
HRV_complaints3 <- merge(HRV_complaints2, model_pvals3, by= 'ID', all.x=T)


#fix colours to be the same as in B
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(6)

Fig5E <- ggplot(HRV_complaints3, aes(x=Day, y = value.x*100, colour = variable)) +  
  geom_line(linewidth = 2.5) + 
  xlim(c(-3,10)) + 
  ylab('Percentage of children \nwith symptoms') + 
  xlab('Days since HRV acquisition') + 
  geom_vline(xintercept = -1, linetype = 'dashed')  +
  geom_text(aes(label=value.y), nudge_y = 10, size = 7, fontface = "bold") +  
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  scale_colour_manual(values = cols) 
Fig5E

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig5E_symptoms_HRV_over_time_incstats",
                       timestamp, ".pdf")
ggsave(pdf_filename, width = 5,height = 4)

#sample size:
HRV_time_symptom_N <- mrg_HRV3 %>%
  filter(!is.na(DayToVirus)) %>%
  #group_by(DayToVirus) %>%
  summarise(
    n_total    = n(),
    Any        = sum(RespiratoryComplaints == "Yes", na.rm = TRUE),
    Sneeze     = sum(Sneeze == "Yes", na.rm = TRUE),
    Cough      = sum(Cough == "Yes", na.rm = TRUE),
    Rhinorrhea = sum(Rhinorrhea == "Yes", na.rm = TRUE),
    Ear        = sum(Ear == "Yes", na.rm = TRUE),
    Fever      = sum(Fever == "Yes", na.rm = TRUE),
    .groups = "drop"
  )# %>%
  #arrange(DayToVirus)
HRV_time_symptom_N




###################################################
#                                                 #
#                 FIGURE 5F                       #
#                                                 #
###################################################

#there is no filter on samples that are absent and then negative bacteria are put as 0 density so presence/absence is included in the model.
m <- glmer(as.factor(RespiratoryComplaintsYN) ~ LytADens + StaphDens + HaemophilusDens +
             MoraxellaDens +  DaysCosinus + AgeYears +  (1|StudyID), data = mrg_all4_Example, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
Any_bac <- summary(m)[[10]]


m <- glmer(as.factor(RhinorrheaYN) ~ LytADens + StaphDens + HaemophilusDens +
             MoraxellaDens  +  DaysCosinus +  AgeYears +   (1|StudyID), data = mrg_all4_Example, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
Snotty_bac <- summary(m)[[10]]


m <- glmer(as.factor(SneezeYN) ~ LytADens + StaphDens + HaemophilusDens +
             MoraxellaDens  + DaysCosinus +  AgeYears +   (1|StudyID), data = mrg_all4_Example, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
Sneeze_bac <- summary(m)[[10]]


m <- glmer(as.factor(CoughYN) ~ LytADens + StaphDens + HaemophilusDens +
             MoraxellaDens  + DaysCosinus +   AgeYears +   (1|StudyID), data = mrg_all4_Example, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
Cough_bac <- summary(m)[[10]]


m <- glmer(as.factor(EarYN) ~ LytADens + StaphDens + HaemophilusDens +
             MoraxellaDens  + DaysCosinus +  AgeYears +   (1|StudyID), data = mrg_all4_Example, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
Ear_bac <- summary(m)[[10]]


#combine the different symptoms
all_bacteria_p <- data.frame(Any = Any_bac[2:5,4], 
                             Rhinorrhea = Snotty_bac[2:5,4], 
                             Cough = Cough_bac[2:5,4], 
                             Sneeze = Sneeze_bac[2:5,4], 
                             Ear = Ear_bac[2:5,4])
all_bacteria_p2 <- all_bacteria_p

#export p values:
all_bacteria_p4 <- (all_bacteria_p2*5)
capture.output(all_bacteria_p4, file = "Fig5F_statssummary.txt")

rownames(all_bacteria_p2) <- gsub('Dens', '', rownames(all_bacteria_p2))
all_bacteria_p3 <- all_bacteria_p2
all_bacteria_p3[all_bacteria_p2*5 > 0.05] <- '' 
all_bacteria_p3[all_bacteria_p2*5 < 0.05] <- '*' 
all_bacteria_p3[all_bacteria_p2*5 < 0.01] <- '**' 
all_bacteria_p3[all_bacteria_p2*5 < 0.001] <- '***' 


all_bacteria_est <- data.frame(Any = Any_bac[2:5,1], 
                               Rhinorrhea = Snotty_bac[2:5,1], 
                               Cough = Cough_bac[2:5,1], 
                               Sneeze = Sneeze_bac[2:5,1], 
                               Ear = Ear_bac[2:5,1])
rownames(all_bacteria_est) <- gsub('Dens', '', rownames(all_bacteria_est))

# Convert to matrices
mat_est  <- as.matrix(all_bacteria_est)
mat_annot <- as.matrix(all_bacteria_p3)

# Ensure same dimnames
all.equal(rownames(mat_est), rownames(mat_annot))
all.equal(colnames(mat_est), colnames(mat_annot))

Fig5F <- pheatmap(mat_est, display_numbers = mat_annot, cluster_cols = F,
                  fontsize_number = 20)
Fig5F

pdf('Fig5F_heatmap_bac_density_vs_symptoms.pdf', width = 4.5, height = 4,onefile=T)
Fig5F
dev.off()

#sample size:
symptom_bac_table <- mrg_all5 %>%
  mutate(
    LytA_pos     = as.integer(LytADens > 0),
    Staph_pos    = as.integer(StaphDens > 0),
    Haemophilus_pos = as.integer(HaemophilusDens > 0),
    Moraxella_pos   = as.integer(MoraxellaDens > 0)
  )

# Put symptoms in long format
symptom_bac_long <- symptom_bac_table %>%
  pivot_longer(
    cols = c(RespiratoryComplaints, Sneeze, Cough, Rhinorrhea, Ear, Fever),
    names_to = "Symptom",
    values_to = "SymYes"
  )

sym_bac_counts <- symptom_bac_long %>%
  filter(SymYes == "Yes") %>%   # only symptom-positive rows
  group_by(Symptom) %>%
  summarise(
    N_symptom_pos   = n(),
    LytA_n          = sum(LytA_pos == 1, na.rm = TRUE),
    Staph_n         = sum(Staph_pos == 1, na.rm = TRUE),
    Haemophilus_n   = sum(Haemophilus_pos == 1, na.rm = TRUE),
    Moraxella_n     = sum(Moraxella_pos == 1, na.rm = TRUE),
    .groups = "drop"
  )

sym_bac_counts

write.table(sym_bac_counts, "Fig5f_samplesize.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


###################################################
#                                                 #
#                 FIGURE 5G                       #
#                                                 #
###################################################

bacteria <- mrg_all4[mrg_all4$Assay %in% c('Haemophilus_influenzae',
                                           'Spn(lytA)'),]

bacteria_pos <- bacteria[bacteria$InfectionLongerTrue ==1     , ]

bacteria_snot <- bacteria_pos[!is.na(bacteria_pos$Rhinorrhea_Level),]

Fig5G <- ggplot(bacteria_snot, aes(x = as.factor(Rhinorrhea_Level), y = log10_avgconc)) + 
  geom_violin(aes(fill = Rhinorrhea_Level, alpha = Rhinorrhea_Level)) + 
  geom_boxplot(width = 0.2) + 
  ylab("Log10 (genome copies/mL +1)") + 
  xlab("Rhinorrhea level") + 
  facet_grid(Assay ~ InfectionAnyVirus) +
  theme(
    legend.position = "none",
    strip.background = element_blank()
  )
Fig5G

pdf('Fig5G_density and virus vs Snotty nose.pdf', width = 4, height = 4)
Fig5G
dev.off()


#look at whether density is better explained by virus or by symptoms
Spn_snot <- bacteria_snot[bacteria_snot$Assay == 'Spn(lytA)',]
HI_snot <- bacteria_snot[bacteria_snot$Assay != 'Spn(lytA)',]
m <- lmer(Rhinorrhea_Level ~ log10_avgconc  + DaysCosinus + AgeYears + InfectionAnyVirus+  (1|StudyID), data = Spn_snot)
summary(m)
capture.output(summary(m), file = "Fig5G_statssummary_RhinoSpn.txt")

m <- lmer(Rhinorrhea_Level ~ log10_avgconc + DaysCosinus  + AgeYears + InfectionAnyVirus+  (1|StudyID), data = HI_snot)
summary(m)
capture.output(summary(m), file = "Fig5G_statssummary_RhinoHI.txt")

#stratify stats per virus or not
Spn_snot_virus <- bacteria_snot[bacteria_snot$Assay == 'Spn(lytA)' & bacteria_snot$InfectionAnyVirus == 'AnyVirus',]
m <- lmer(Rhinorrhea_Level ~ log10_avgconc + DaysCosinus  + AgeYears + (1|StudyID), data = Spn_snot_virus)
summary(m)
capture.output(summary(m), file = "Fig5G_statssummary_SpnAnyVirus.txt")

Spn_snot_novirus <- bacteria_snot[bacteria_snot$Assay == 'Spn(lytA)' & bacteria_snot$InfectionAnyVirus != 'AnyVirus',]
m <- lmer(Rhinorrhea_Level ~ log10_avgconc + DaysCosinus  + AgeYears + (1|StudyID), data = Spn_snot_novirus)
summary(m)
capture.output(summary(m), file = "Fig5G_statssummary_SpnNoVirus.txt")

HI_snot_virus <- bacteria_snot[bacteria_snot$Assay != 'Spn(lytA)' & bacteria_snot$InfectionAnyVirus == 'AnyVirus',]
m <- lmer(Rhinorrhea_Level ~ log10_avgconc + DaysCosinus  + AgeYears + (1|StudyID), data = HI_snot_virus)
summary(m)
capture.output(summary(m), file = "Fig5G_statssummary_HinfAnyVirus.txt")

HI_snot_novirus <- bacteria_snot[bacteria_snot$Assay != 'Spn(lytA)' & bacteria_snot$InfectionAnyVirus != 'AnyVirus',]
m <- lmer(Rhinorrhea_Level ~ log10_avgconc + DaysCosinus  + AgeYears + (1|StudyID), data = HI_snot_novirus)
summary(m)
capture.output(summary(m), file = "Fig5G_statssummary_HinfNoVirus.txt")


#nr samples:
fig_data      <- bacteria_snot
group_cols    <- c("Rhinorrhea_Level", "Assay", "InfectionAnyVirus")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig5G_samplesize.txt")




###################################################
#                                                 #
#                 FIGURE 5H                       #
#                                                 #
###################################################

#make a symptom onset, snotty nose onset and cough onset
#we need at least 3 negative samples before and 2 positive samples afterwards, only look for first time
mrg_all4_Example.split <- split(mrg_all4_Example, mrg_all4_Example$StudyID)

################this loop will throw many errors, which are caught with the try
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
colnames(mrg_symptomonset) #take the onset ones
colnames(mrg_symptomonset)[141:143]
mrg_all50 <- merge(mrg_all4, mrg_symptomonset[,c(1, 141:143)], by = 'SampleID_Day')


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
lm <- lmer(log10_avgconc ~ as.factor(DayToSNOTonset2) + DaysCosinus + AgeYears + (1|StudyID), data = mrg8a1) #do stats for all days relative to baseline
acq_res <- data.frame(lmerTest:::get_coefmat(lm)[2:nrow(lmerTest:::get_coefmat(lm)),c(1,5)]) #extract the results from the linear model

mrg8b1 <- mrg8b[mrg8b$DayToSNOTonset > -8 & mrg8b$DayToSNOTonset < 11, ] #filter out smaples between x days before or after viral onset
mrg8b1$DayToSNOTonset2 <- mrg8b1$DayToSNOTonset #make a new column and change the 5 days prior to onset to the baseline
mrg8b1$DayToSNOTonset2[mrg8b1$DayToSNOTonset2 %in% c(-3:-1)] <- '-1_baseline'
lm <- lmer(log10_avgconc ~ as.factor(DayToSNOTonset2) + DaysCosinus + AgeYears + (1|StudyID), data = mrg8b1) #do stats for all days relative to baseline
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

#combine spn and haemophilus and plot
mrg8 <- rbind(mrg8a, mrg8b)
Fig5H <- ggplot(mrg8, aes(x=DayToSNOTonset, y=FC_dens, colour=Assay)) + 
  xlim(c(-7,10)) + 
  geom_hline(yintercept=0, linetype='dashed')  + 
  geom_vline(xintercept=0, linetype='dashed')  +
  #stat_summary(fun.data = mean_se,  geom = "errorbar", width=0.2) +
  stat_summary(fun.data = mean_cl_normal,  
               geom = "errorbar", width=0.2)+
  stat_summary(fun.y = mean, geom = "point",  size =4) +
  facet_grid(.~Assay) + 
  ylab("Fold-change density (Log10))") + 
  xlab("Day to rhinorrhea onset") + 
  theme(legend.position = "none",
        strip.background = element_blank()) +
  geom_text(data= acq_res3, aes(label=P2, x=as.numeric(x)), y=3.5, colour = 'black')
Fig5H

ggsave('Fig5H_bacterial_density_relative_to_symptoms.pdf', width = 5, height = 4)


#nr samples:
fig_data      <- mrg8a
group_cols    <- c("DayToSNOTonset", "Assay")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig5H_samplesize.txt")



###################################################
#                                                 #
#                 FIGURE 5I                       #
#                                                 #
###################################################

#do VAR to see how lyta,virus and symptoms interact
mycols <- c("StudyID",
            "Day",
            "LytADens",
            'HaemophilusDens',
            'InfectionAnyVirus',
            'Rhinorrhea_Level',
            'RhinorrheaYN', 
            'TSS'
)

mrg_mlvar <- mrg_all4_Example %>% dplyr::select(all_of(mycols)) 
mrg_mlvar$Day <- as.numeric(gsub('D', '', mrg_mlvar$Day)) 
mrg_mlvar <- mrg_mlvar[!is.na(mrg_mlvar$Day) & !is.na(mrg_mlvar$Rhinorrhea_Level),]
mrg_mlvar$InfectionAnyVirus <- as.numeric(as.factor(mrg_mlvar$InfectionAnyVirus))-1

VAR1 <- mlVAR(mrg_mlvar, vars = c('LytADens', 'Rhinorrhea_Level', 'InfectionAnyVirus'), 
              idvar = 'StudyID', lags = 1, contemporaneous = 'correlated',
              estimator = 'lmer',)

print(VAR1)
summary(VAR1)

pdf('Fig5I_VAR_Spn_snot_infection_lag1.pdf', width = 5, height = 3)
par(mfrow=c(1,2))
plot(VAR1, posCol = 'Red', negCol ='Blue', title = 'Temporal Spn - lag 1', label.cex = 1, edge.width=5,
     edge.labels =T, label.scale = F, edge.label.color = 'black', mode = 'direct')
dev.off()



###########################################################################33
save.image('Fig5_symptoms.Rdata')
