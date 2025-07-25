#####
#This script uses the mlvar package to calculate correlations between bacteria
#https://arxiv.org/pdf/1609.04156
#https://cran.r-project.org/web/packages/mlVAR/mlVAR.pdf
#Spn (lytA) / Haemophilus / Moraxella / Staph 
#it uses the RAW (log-transformed) data or with data normalized to 16S 


rm(list=ls())

setwd("Z:\\Projects\\SAMSAM\\Manuscripts\\Figures\\Supplemental figure 7_bacterialinteractions/")


# Library
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(lmerTest)
library(reshape2)
library(readxl)
library(mlVAR)
library(lmms)
library(pheatmap)

#save the session info
sink("sessionInfo.txt")
sessionInfo()
sink()


#load in database
mrg <- read_xlsx('Z:\\Projects\\SAMSAM\\IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20250722.xlsx')
mrg0 <- mrg[mrg$Antibiotics != 'Yes',] #remove the antibiotics samples as it may affect results


mrg1 <- mrg[mrg$StudyID == 'SAM006',] #select the child with the antibiotics

#################################################################
#
#Figure S6A
#
#############################################################################
mrg1a <- mrg1 %>%
  group_by(Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLonger) == 1
  ) %>%
  ungroup() #make a new dataframe and only select assays that are present


xmin <- min(mrg1a$Day2[mrg1a$Antibiotics == 'Yes'])
xmax <- max(mrg1a$Day2[mrg1a$Antibiotics == 'Yes'])
ggplot(mrg1a, aes(x=Day2, y=log10_avgconc, colour = Assay)) +  
  geom_rect(xmin = xmin, xmax=xmax, ymin=0, ymax=7, colour = 'grey', fill = 'grey', alpha = 0.2)+
  geom_point(aes(shape = AssayType))  + 
  geom_line(aes(group=Assay, linetype = AssayType)) + 
  scale_color_brewer(palette = 'Paired') + theme_bw() + 
  facet_grid(.~AssayType) + 
  annotate('text', x= 10, y=7.2, label = 'Antibiotics (Phenethicillin)', colour = 'black')
ggsave('FigS7A_SAM006_lineplots.pdf', width =6, height =2.5)






###############################################
#
#Fig S7B-C
#loess regression VAR
###############################################
mycols <- c("StudyID",
            "Day",
            "Assay",
            "log10_avgconc")



mrg2 <- mrg0 


#remove 18PS,  EAV and 16S as assays
mrg2 <- mrg2[mrg2$Assay == 'Spn(lytA)' | 
               grepl("Haemophilus_influenzae", mrg2$Assay)| 
               grepl("aureus", mrg2$Assay)|
               grepl("Moraxella_catarrhalis", mrg2$Assay),]

mrg3 <- mrg2 %>% dplyr::select(all_of(mycols)) 
#mrg3$Day <- factor(mrg3$Day, levels = unique(mrg3$Day)) 
mrg3$Day <- as.numeric(gsub('D', '', mrg3$Day)) 
mrg3 <- mrg3[!is.na(mrg3$Day),]


mrg3$StudyID_Assay <- paste0(mrg3$StudyID, '_', mrg3$Assay)

mrg3.split <- split(mrg3, f= mrg3$StudyID_Assay)

for(i in 1:length(mrg3.split)){
  mrg5 <- data.frame(mrg3.split[[i]])
  names <- mrg5$X2[!is.na(mrg5$value)]
  test <- loess(log10_avgconc ~ Day, data=mrg5, span=0.3)
  mrg3.split[[i]]$LoessFit <- data.frame(fitted(test))[,1]
}
mrg3.loess <- do.call(rbind, mrg3.split)

ggplot(mrg3.loess, aes(x=as.numeric(Day), y=log10_avgconc, colour = Assay)) +
  geom_point()+
  facet_wrap(.~StudyID, scales = 'free_x') +
  geom_line(aes(colour = Assay, y=LoessFit)) + 
  ggtitle('Loess fit')


mycols2 <- c("StudyID",
             "Day",
             "Assay",
             "LoessFit")
mrg4.loess <- mrg3.loess %>% dplyr::select(all_of(mycols2)) 
#mrg3$Day <- factor(mrg3$Day, levels = unique(mrg3$Day)) 
mrg4.loess$Day <- as.numeric(gsub('D', '', mrg4.loess$Day)) 
mrg4.loess <- mrg4.loess[!is.na(mrg4.loess$Day),]

mrg5.loess <- mrg4.loess %>%
  pivot_wider(names_from = Assay, 
              values_from = c(LoessFit))
colnames(mrg5.loess)[colnames(mrg5.loess) == "Spn(lytA)"] <- 'Spn'      


VAR.splines <- mlVAR(mrg5.loess, vars = c('Spn', 'Haemophilus_influenzae',  'Staphylococcus_aureus', 'Moraxella_catarrhalis'), 
                     idvar = 'StudyID', lags = 1, contemporaneous = 'correlated',
                     estimator = 'lmer')


print(VAR.splines)
summary(VAR.splines)

dev.off()

pdf('FigS7BC_VAR_plots_loess_fits.pdf', width = 8, height = 3)
par(mfrow=c(1,3))
plot(VAR.splines, posCol = 'Red', negCol ='Blue', title = 'Temporal - loess', label.cex = 1, label.scale = F)
plot(VAR.splines, type = 'contemporaneous', posCol = 'Red', negCol ='Blue', title = 'Contemporaneous - loess', label.cex = 1, label.scale = F)
plot(VAR.splines, type = 'between', posCol = 'Red', negCol ='Blue', title = 'Between Subject - loess', label.cex = 1, label.scale = F)
dev.off()




###############################################
#
#Fig S7D-E
#16S-normalized VAR
###############################################
mrg2.norm <- mrg0[mrg0$Assay == 'Spn(lytA)' | 
                   grepl("Haemophilus_influenzae", mrg0$Assay)| 
                   grepl("aureus", mrg0$Assay)|
                   grepl("16S", mrg0$Assay)|
                   grepl("Moraxella_catarrhalis", mrg0$Assay),]

mrg3.norm <- mrg2.norm %>% dplyr::select(all_of(mycols)) 
#mrg3.norm$Day <- factor(mrg3.norm$Day, levels = unique(mrg3.norm$Day)) 
mrg3.norm$Day <- as.numeric(gsub('D', '', mrg3.norm$Day)) 
mrg3.norm <- mrg3.norm[!is.na(mrg3.norm$Day),]


mrg4.norm <- mrg3.norm %>%
  pivot_wider(names_from = Assay, 
              values_from = c(log10_avgconc))
colnames(mrg4.norm)[colnames(mrg4.norm) == "Spn(lytA)"] <- 'Spn'                    


mrg5.norm <- mrg4.norm #do a normalization based on concentration of 16S
for(i in 3:ncol(mrg5.norm)){
  mrg5.norm[,i] <- mrg5.norm[,i] / mrg4.norm[,colnames(mrg4.norm) =='16S']
}
mrg6.norm <- mrg5.norm[!is.na(mrg5.norm$'16S'),] #remove samples without 16S data


VAR.norm <- mlVAR(mrg6.norm, vars = c('Spn', 'Haemophilus_influenzae', 'Staphylococcus_aureus', 'Moraxella_catarrhalis'), 
                  idvar = 'StudyID', lags = 1, contemporaneous = 'correlated',
                  estimator = 'lmer')


print(VAR.norm)
summary(VAR.norm)
pdf('FigS7DE_VAR_plots_bacteria_16Snormalized_lag1.pdf', width = 8, height = 3)
par(mfrow=c(1,3))
plot(VAR.norm, posCol = 'Red', negCol ='Blue', title = 'Temporal - 16S normalized', label.cex = 1, label.scale = F)
plot(VAR.norm, type = 'contemporaneous', posCol = 'Red', negCol ='Blue', 
     title = 'Contemporaneous - 16S normalized', label.cex = 1, label.scale = F)
plot(VAR.norm, type = 'between', posCol = 'Red', negCol ='Blue', title = 'Between Subject - 16S normalized', label.cex = 1, label.scale = F)
dev.off()








############################
#
#Fig S7F-G
#set the counter of max number of times in a row we see a pathogen
##################################################################
virus_pos <- unique(mrg0[mrg0$AssayType == 'Viral' & mrg0$PositiveYesNo == 1,]$SampleID_Day) #remove any sample with a virus

mrg.a <- mrg0[!mrg0$SampleID_Day %in% virus_pos,] #remove all samples that are virus positive



mycols <- c("StudyID",
            "Day",
            "Assay",
            "log10_avgconc")



mrg2a <- mrg.a 


#remove 18PS,  EAV and 16S as assays
mrg2a <- mrg2a[mrg2a$Assay == 'Spn(lytA)' | 
                          grepl("Haemophilus_influenzae", mrg2a$Assay)| 
               grepl("aureus", mrg2a$Assay)|
                          grepl("Moraxella_catarrhalis", mrg2a$Assay),]

mrg3a <- mrg2a %>% dplyr::select(all_of(mycols)) 
mrg3a$Day <- as.numeric(gsub('D', '', mrg3a$Day)) 
mrg3a <- mrg3a[!is.na(mrg3a$Day),]

mrg4a <- mrg3a %>%
  pivot_wider(names_from = Assay, 
              values_from = c(log10_avgconc))


colnames(mrg4a)[colnames(mrg4a) == "Spn(lytA)"] <- 'Spn'


VARv <- mlVAR(mrg4a, vars = c('Spn', 'Haemophilus_influenzae',  'Staphylococcus_aureus', 'Moraxella_catarrhalis'), 
      idvar = 'StudyID', lags = 1, contemporaneous = 'correlated',
      estimator = 'lmer',)

print(VARv)
summary(VARv)

pdf('FigS7FG_VAR_plots_bacteria_rawdata_lag1_wovirus.pdf', width = 8, height = 3)
par(mfrow=c(1,3))
plot(VARv, posCol = 'Red', negCol ='Blue', title = 'Temporal - raw', label.cex = 1, label.scale = F)
plot(VARv, type = 'contemporaneous', posCol = 'Red', negCol ='Blue', title = 'Contemporaneous - raw', label.cex = 1, label.scale = F)
plot(VARv, type = 'between', posCol = 'Red', negCol ='Blue', title = 'Between Subject - raw', label.cex = 1, label.scale = F)
dev.off()






#####################################################################
#
#Table S7H
#
#####################################################################
######################################################
#check what is the chance of clearance on the next day of Spn and if it is related to other factors
###########################################################
mrg_clear <- mrg0 %>%
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


#select only Spn to look at clearance effects
mrg_clear_spn <- mrg_clear2[grepl('Spn', mrg_clear2$Assay) &
                              !grepl('lytA', mrg_clear2$Assay) &
                              !grepl('piaB', mrg_clear2$Assay) &
                              !grepl('6C/D', mrg_clear2$Assay),]

tab <- table(mrg_clear_spn$SpnNumber2, mrg_clear_spn$ClearNextDay)
tab
tab / rowSums(tab)



##############################################
save.image('FigureS7_bacteria.Rdata')
