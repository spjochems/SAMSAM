#Pneumolines

rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\FigS11\\")

# Library
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(lmerTest)
library(reshape2)
library(readxl)
library(lmms)

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#load in database
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))
mrg1 <- mrg[mrg$Antibiotics != 'Yes',] #remove the antibiotics samples as it may affect results




###################################################
#                                                 #
#                 FIGURE S9A                      #
#                                                 #
###################################################
 
#plots of only ones with an infection
#filter out targets where we have at least 1 infection
mrg9 <- mrg %>%
  group_by(StudyID, Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLongerTrue) ==1 
  ) %>%
  ungroup()

spn2 <- mrg9 %>% filter(grepl("^Spn", Assay))  %>% 
  filter(!grepl("18A", Assay)) %>% filter(!grepl("19B", Assay))

# Create a PDF file to save all plots
ggplot(spn2, aes(x = Day2, y = log10_avgconc)) +  # Convert Day to numeric for regression
  geom_point(aes(color = Assay), size = 1) +  # Points colored by Assay
  #geom_line(aes(color = Assay), size = 0.7) +  # Line connecting the points
  geom_smooth(aes(color = Assay), method = "loess", se = FALSE, size = 0.7, span = .6) +    # Add separate smoothing line for each Assay
  theme_minimal() +
  facet_wrap(.~StudyID, ncol = 7) + 
  theme(
    axis.text.x = element_text(), 
    legend.position = "right"
  ) +  # Rotate x-axis labels
  labs(
    title = 'S. pneumoniae carriage',  # Title for each StudyID
    x = NULL, y = "log10(average_conc)"  # y-axis label with log scale
  )

ggsave("FigS11A_spn_atleast2.pdf", width = 12, height = 10)



###################################################
#                                                 #
#                 FIGURE S11B                      #
#                                                 #
###################################################

#check what is the chance of clearance on the next day of Spn and if it is related to other factors

mrg_clear <- mrg1 %>%
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
    
    sample_check <- each4$SampleID_Day[i] #check which are the relevant sample IDs
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

rm(mrg_clear.split)

#select only Spn to look at clearance effects
mrg_clear_spn <- mrg_clear2[grepl('Spn', mrg_clear2$Assay) &
                              !grepl('lytA', mrg_clear2$Assay) &
                              !grepl('piaB', mrg_clear2$Assay) &
                              !grepl('6C/D', mrg_clear2$Assay),]

tab <- table(mrg_clear_spn$SpnNumber2, mrg_clear_spn$ClearNextDay)
tab <- cbind(tab, `3` = round(tab[, "1"] / rowSums(tab) * 100, 1))
tab

colnames(tab) <- c("Does not clear next day", "Clears next day", "%cleared next day")
tab

capture.output(tab, file = "FigS11B_Spnclearance.txt")




###################################################
#                                                 #
#                 FIGURE S11C                      #
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


rm(mrg_clear.split)

table(mrg_clear2$ClearNextDay)

#select only Spn to look at clearance effects
mrg_clear_spn <- mrg_clear2[grepl('Spn', mrg_clear2$Assay) &
                              !grepl('lytA', mrg_clear2$Assay) &
                              !grepl('piaB', mrg_clear2$Assay) &
                              !grepl('6C/D', mrg_clear2$Assay),]

#try stats without correcting for age and serotype for reviewers:
m <- glm(as.factor(ClearNextDay) ~ SpnNumber2, data = mrg_clear_spn, family = binomial)
coef(summary(m))
coefs <- coef(summary(m))

OR_adjusted <- as.data.frame(coefs[-1, , drop = FALSE])
colnames(OR_adjusted) <- c("Estimate", "SE", "z", "p")

OR_adjusted$SpnNumber <- gsub("SpnNumber2", "", rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`SE`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`SE`)

OR_adjusted$stars <- ifelse(OR_adjusted$p < 0.001, "***",
                            ifelse(OR_adjusted$p < 0.01, "**",
                                   ifelse(OR_adjusted$p < 0.05, "*", "")))
write_tsv(OR_adjusted, "FigS11C_statssummary.txt")

ggplot(OR_adjusted, aes(x = OR, y = SpnNumber)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower),
                 size = 0.5, height = 0.2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") + 
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  xlab("Adjusted OR of clearing \na serotype next day") + 
  geom_vline(xintercept = 1, linetype = "dashed")

ggsave('FigS11C_AdjustedOR_Spn_clearanceNextDay_simple.pdf', width = 3, height = 2.5)
