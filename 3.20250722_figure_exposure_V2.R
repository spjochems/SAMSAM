#####
#This makes plots of exposure
#update 9-6. Fixed the forest plot for exposure vs activities to be of standard error and not CI

rm(list=ls())

setwd("Z:\\Projects\\SAMSAM\\IntegrativeAnalysis\\3.NA_detection\\1.20250126_linktodatabase/output/")

# Library
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(reshape2)
library(readxl)
library(survival)
library(ggsurvfit)
library(lmerTest)
library(ggmosaic)



#create a function that fixes the display of p-values
disp_pval <- function(p){
  if(p<0.001){
    p <- scales::scientific(p, digits= 3)
  } else {
    p <- round(p, 3)
  }
  return(p)
}


#load in database
mrg <- data.frame(read_xlsx('Z:/Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20250722.xlsx'))

setwd('Z:\\Projects\\SAMSAM\\Manuscripts\\Figures\\Figure2_exposure')


#save the session info
sink("sessionInfo.txt")
sessionInfo()
sink()


mrg2 <- mrg[!mrg$Assay %in% c('16S', 'PIV1','PIV2','Spn20',
                              'Moraxella_catarrhalis', 'Haemophilus_influenzae',
                              'Spn(lytA)', 'Spn(piaB)'),] #remove 16S and the none or always detected assays





######################
#
#Fig. 2B
#
###############################
#get association with school going. Select only the timepoints of exposure or not
#only select the ones where we have a potential exposure,
mrg4 <- mrg2[mrg2$Repeats < 2 & #remove any with ongoing infection
               mrg2$Day2 >4 & #remove any in first 4 days as there is no exposiure possible according to definition
               !(mrg2$Repeats ==0 & mrg2$InfectionLonger == 1),]#remove the ones where there is an ongoing infection but intermittant negative 

mrg4$Assay2 <- mrg4$Assay
mrg4$Assay2[mrg4$AssayType == 'Viral'] <- 'AnyVirus' 
mrg4$Assay2[grepl('Spn',mrg4$Assay) & 
              !grepl('lytA',mrg4$Assay) & 
              !grepl('piaB',mrg4$Assay) & 
              !grepl('Spn6C/D',mrg4$Assay)] <- 'AnySerotype'


mrg5a <- mrg4 %>% 
  filter(!grepl("NA", Swimmingpool)) %>% 
  filter(grepl("Evening", PeriodCollection) | grepl("Afternoon", PeriodCollection)  )

m <- glmer(as.factor(ExposureEvent) ~ SchoolDaycare + OtherActivity + Swimmingpool  + (1|Assay) + (1|StudyID), data = mrg5a, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)
#CIs <- confint(m)

summary(m)[[10]][,1:2]
OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Activity <- gsub('Yes', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - OR_adjusted$Std..Error)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + OR_adjusted$Std..Error)
OR_adjusted$Pval <- paste0('p=', sapply(OR_adjusted$Pr...z.., FUN = "disp_pval"))
ggplot(OR_adjusted, aes(x=OR, y=Activity)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  geom_vline(xintercept = 1, linetype = 'dashed') + 
  xlab('Adjusted Odds Ratio LMM') 
ggsave('Fig2B_AdjustedOddsRatioExposure.pdf', width = 3, height = 2)






####################################
#
#Fig. 2C
#
######################################
mrg_exposure <- mrg2[mrg2$ExposureEvent == 1,]
exposed_paths <- data.frame(table(table(mrg_exposure$SampleID_Day)))

pl3 <- ggplot(exposed_paths, 
              aes(y=Var1, x=Freq)) + 
  geom_bar(stat='identity', 
           fill="steelblue", 
           colour="black") + 
  scale_x_continuous(expand = c(0,0)) +
  geom_text(aes(label = Freq), hjust = 1, colour = 'white')+ 
  ylab('Number of co-detected exposure events')
pl3
timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("fig2C_numberofcoexposure_", timestamp, ".pdf")
ggsave(pdf_filename, width = 3, height = 2)

sum(as.numeric(exposed_paths$Var1) * exposed_paths$Freq ) #check the total number of exposures






####################################
#
#Fig. 2D
#
######################################
mrg.split <- split(mrg2, f=mrg2$Assay)
for(i in 1:length(mrg.split)){
  pres <- table(mrg.split[[i]]$ExposureEvent, mrg.split[[i]]$StudyID)
  acq <- table(mrg.split[[i]]$AcquisitionLonger, mrg.split[[i]]$StudyID)
  abort <- table(mrg.split[[i]]$AbortiveExposureLonger, mrg.split[[i]]$StudyID)
  
  #add the presence table
  if(nrow(pres)>1){
    if(i == 1){ 
      presence <- pres[2,]
    }else{
      presence <- cbind(presence, pres[2,])  
      colnames(presence)[ncol(presence)] <- names(mrg.split)[i]
    }
  }
  
  #add the acquisition table
  if(nrow(acq)>1){
    if(i == 1){ 
      acquisition <- acq[2,]
    }else{
      acquisition <- cbind(acquisition, acq[2,])  
      colnames(acquisition)[ncol(acquisition)] <- names(mrg.split)[i]
    }
  }
  
  
  #add the abortive table
  if(nrow(abort)>1){
    if(i == 1){ 
      abortive <- abort[2,]
    }else{
      abortive <- cbind(abortive, abort[2,])  
      colnames(abortive)[ncol(abortive)] <- names(mrg.split)[i]
    }
  }
  
}
colnames(abortive)[1] <- colnames(acquisition)[1] <- colnames(presence)[1] <- names(mrg.split)[1]


length(which(grepl('Spn', colnames(presence)) & !grepl('Spn6C/D', colnames(presence))))
colnames(presence)[which(colnames(presence) %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay))]



exposure_perpathogen <- data.frame(colSums(presence))
abortive_perpathogen <- data.frame(colSums(abortive))
acquisition_perpathogen <- data.frame(colSums(acquisition))
expab <- merge(exposure_perpathogen, abortive_perpathogen, by= 'row.names', all.x=T)
all_perpath <- merge(expab, acquisition_perpathogen, by.x= 'Row.names', by.y= 'row.names', all.x=T)
colnames(all_perpath) <- c('Assay', 'Exposures', 'Abortive', 'Acquisition')
all_perpath[is.na(all_perpath)] <- 0

#calculate the number of observations per spn or total virus
spna <- sum(all_perpath[grepl('Spn', all_perpath$Assay) & !grepl('Spn6C/D', all_perpath$Assay ),2])
spnb <- sum(all_perpath[grepl('Spn', all_perpath$Assay) & !grepl('Spn6C/D', all_perpath$Assay ),3])
spnc <- sum(all_perpath[grepl('Spn', all_perpath$Assay) & !grepl('Spn6C/D', all_perpath$Assay ),4])

vira <- sum(all_perpath[all_perpath$Assay %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay),2])
virb <- sum(all_perpath[all_perpath$Assay %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay),3])
virc <- sum(all_perpath[all_perpath$Assay %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay),4])


all_perpath2 <- rbind(all_perpath, 
                      c('AnySpnSerotype', spna, spnb, spnc),
                      c('AnyVirus', vira, virb, virc))
all_perpath2[,2] <- as.numeric(all_perpath2[,2] )
all_perpath2[,3] <- as.numeric(all_perpath2[,3] )
all_perpath2[,4] <- as.numeric(all_perpath2[,4] )

all_perpath2$ExposureRate <- all_perpath2[,2] / 45
all_perpath2$AcquisitionPercentage <-  all_perpath2[,4] / all_perpath2[,2] * 100
all_perpath2$AcquisitionPercentage2 <- round(all_perpath2$AcquisitionPercentage, 0) 
all_perpath2$AcquisitionPercentage2[all_perpath2$Exposures <20] <- NA


all_perpath2$Assay <- factor(all_perpath2$Assay, 
                             levels = all_perpath2$Assay[order(all_perpath2$ExposureRate)]) 

pl1 <- ggplot(all_perpath2, 
              aes(y=Assay, x=ExposureRate)) + 
  geom_bar(stat='identity', 
           fill="steelblue", 
           colour="black")+
  xlab('Exposure Rate per child follow up')

pl1

pl2 <- ggplot(all_perpath2, aes(x=1, y = Assay, fill = AcquisitionPercentage2)) +
  geom_tile(color = "black") +
  geom_text(aes(label = AcquisitionPercentage2), color = "white", size = 4) +
  coord_fixed()

plot_grid(pl1, pl2)
timestamp <- format(Sys.time(), "%Y%m%d")

pdf_filename <- paste0("fig2D_exposurerate_", timestamp, ".pdf")
ggsave(pdf_filename, width = 10, height = 6)







######################################################
#
#Fig 2E
#
###########################################################
mrg4_exposure_virus <- mrg4[mrg4$ExposureEvent == 1 & mrg4$Assay2 == 'AnyVirus',]
exposure_table_virus <- table(mrg4_exposure_virus$AgeYears, mrg4_exposure_virus$Assay)
mrg4_exposure_virus2 <- cbind( exposure_table_virus, table(mrg[which(!duplicated(mrg$StudyID)),]$AgeYears))
for(i in 1:ncol(mrg4_exposure_virus2)){mrg4_exposure_virus2[,i] <- mrg4_exposure_virus2[,i]/mrg4_exposure_virus2[,ncol(mrg4_exposure_virus2)]}
mrg4_exposure_virus3 <- t(mrg4_exposure_virus2[,-13])
mrg4_exposure_virus3 <- rbind(mrg4_exposure_virus3, colSums(mrg4_exposure_virus3 ))
rownames(mrg4_exposure_virus3)[nrow(mrg4_exposure_virus3)] <- 'AnyVirus'
pdf_filename <- paste0("fig2E_heatmap_virusexposure_", timestamp, ".pdf")
pdf(pdf_filename, width = 4, height = 4)
pheatmap::pheatmap(mrg4_exposure_virus3, cluster_cols = F, treeheight_row = 0)
dev.off()

#stats to go with them
mrg4_select <- mrg4[mrg4$Assay2 %in% c('AnyVirus'),]
m <- glmer(as.factor(ExposureEvent) ~ AgeYears + (1|Assay)  + (1|StudyID), data = mrg4_select, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m) #check if exposure is realted to age


mrg4_select <- mrg4[mrg4$Assay %in% c('bocavirus'),]
m <- glmer(as.factor(ExposureEvent) ~ AgeYears   + (1|StudyID), data = mrg4_select, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)




################################################
#
#Figure 2F
#
############################################################
PerChild <- data.frame(Spn_exposure = rowSums(presence[,which(grepl('Spn', colnames(presence)) & !grepl('Spn6C/D', colnames(presence)))]),
                       Spn_abortive = rowSums(abortive[,which(grepl('Spn', colnames(abortive)) & !grepl('Spn6C/D', colnames(abortive)))]),
                       Spn_acquisition = rowSums(acquisition[,which(grepl('Spn', colnames(acquisition)) & !grepl('Spn6C/D', colnames(acquisition)))]),
                       Virus_exposure = rowSums(presence[,which(colnames(presence) %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay))]),
                       Virus_abotive = rowSums(abortive[,which(colnames(abortive) %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay))]),
                       Virus_acquisition = rowSums(acquisition[,which(colnames(acquisition) %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay))]),
                       HRV_exposure = presence[,'HRV'],
                       HRV_abortive = abortive[,'HRV'],
                       HRV_acquisition = acquisition[,'HRV'],
                       bocavirus_exposure = presence[,'bocavirus'],
                       bocavirus_abortive = abortive[,'bocavirus'],
                       bocavirus_acquisition = acquisition[,'bocavirus'],
                       Staphylococcus_aureus_exposure = presence[,'Staphylococcus_aureus'],
                       Staphylococcus_aureus_abortive = abortive[,'Staphylococcus_aureus'],
                       Staphylococcus_aureus_acquisition = acquisition[,'Staphylococcus_aureus']
                       )

demo <- mrg[!duplicated(mrg$StudyID),c(2,6,7)]
PerChild2 <- merge(PerChild, demo, by.x = 'row.names', by.y = 'StudyID')
PerChild2$Spn_acquisitionSucces <- PerChild2$Spn_acquisition / PerChild2$Spn_exposure * 100
PerChild2$VirusacquisitionSucces <- PerChild2$Virus_acquisition / PerChild2$Virus_exposure * 100
PerChild2$HRVacquisitionSucces <- PerChild2$HRV_acquisition / PerChild2$HRV_exposure * 100
PerChild2$bocavirusacquisitionSucces <- PerChild2$bocavirus_acquisition / PerChild2$bocavirus_exposure * 100
PerChild2$StaphaureusacquisitionSucces <- PerChild2$Staphylococcus_aureus_acquisition / PerChild2$Staphylococcus_aureus_exposure * 100



ggplot(PerChild2, aes(x=as.factor(AgeYears), y=VirusacquisitionSucces)) + 
  geom_boxplot(outlier.color = NA, aes(fill = AgeYears)) + geom_jitter(width=0.2, height = 0) + 
  scale_fill_gradient2()
timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("fig2F_virus_acquisition_succes_", timestamp, ".pdf")
ggsave(pdf_filename, width = 5, height = 4)






################################################
#
#Figure 2G
#
############################################################
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
mrg_exposure_pyogenes <- mrg_exposure[mrg_exposure$Assay == "Streptococcus_pyogenes",]



ggplot(PerChild2, aes(x=as.factor(AgeYears), y=bocavirusacquisitionSucces)) + 
  geom_boxplot(outlier.color = NA, aes(fill = AgeYears)) + geom_jitter(width=0.2, height = 0) + 
  scale_fill_gradient2()

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("fig2G_bocavirus_acquisition_succes_", timestamp, ".pdf")
ggsave(pdf_filename, width = 5, height = 4)


m <- glmer(as.factor(AcquisitionLonger) ~ AgeYears + (1|StudyID), data = mrg_exposure_bocavirus, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)








###################################################################
#
#Figure 2H
#
##########################################################
table(mrg_exposure_virus$AcquisitionLonger, mrg_exposure_virus$Staph)
48/(48+29)
10/(10+15)



m <- glmer(as.factor(AcquisitionLonger) ~ Staph + AgeYears + (1|Assay) + (1|StudyID), data = mrg_exposure_virus, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)


ggplot(mrg_exposure_virus) + 
  geom_mosaic(aes(x=product(AcquisitionLonger, Staph), fill = Assay)) +
  facet_grid(~Staph, space = 'free') + 
  scale_fill_discrete() 
pdf_filename <- paste0("fig2H_linkstaph_virusacquisition_", timestamp, ".pdf")
ggsave(pdf_filename, width = 5, height = 4)



###################################################################
#
#Figure 2I
#
##########################################################
m <- glmer(as.factor(AcquisitionLonger) ~  Virus  + AgeYears +  (1|StudyID), data = mrg_exposure_pyogenes, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

table(mrg_exposure_pyogenes$AcquisitionLonger, mrg_exposure_pyogenes$Virus)
table(mrg_exposure_pyogenes$Virus)
prop.table(table(mrg_exposure_pyogenes$AcquisitionLonger, mrg_exposure_pyogenes$Virus), 2)

ggplot(mrg_exposure_pyogenes) + 
  geom_mosaic(aes(x=product(AcquisitionLonger, Virus), fill = AcquisitionLonger)) +
  facet_grid(~Virus, space = 'free') + 
  scale_fill_grey(start = 0.8, end = 0.2) 
pdf_filename <- paste0("fig2I_linkvirus_pyogenes_", timestamp, ".pdf")
ggsave(pdf_filename, width = 3, height = 4)











################################################
#
#Figure 2J
#
############################################################
ggplot(PerChild2, aes(x=as.factor(AgeYears), y=Spn_exposure)) + 
  geom_boxplot(outlier.color = NA, aes(fill = AgeYears)) + geom_jitter(width=0.2, height = 0) + 
  scale_fill_gradient2(high = 'orange')
timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("fig_2J_Spn_exposure_", timestamp, ".pdf")
ggsave(pdf_filename, width = 5, height = 4)





################################################
#
#Figure 2K
#
############################################################
mrg_exposure$Virus <- NA
mrg_exposure$VirusNumber <- NA
mrg_exposure$SpnNumber  <- NA
mrg_exposure$Staph <- NA
mrg_exposure$Haemophilus <- NA
mrg_exposure$LytA <- NA
mrg_exposure$Spnexposures <- NA
mrg_exposure$Totalexposures <- NA




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
  mrg_exposure$Spnexposures[i] <- sum(sample_check2[grepl('Spn',sample_check2$Assay) & 
                                                      !grepl('lytA',sample_check2$Assay) & 
                                                      !grepl('piaB',sample_check2$Assay) & 
                                                      !grepl('Spn6C/D',sample_check2$Assay), 'ExposureEvent'])
  mrg_exposure$Totalexposures[i] <- sum(sample_check2[ !grepl('piaB',sample_check2$Assay) & 
                                                         !grepl('Spn6C/D',sample_check2$Assay), 'ExposureEvent'])
  
  
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

mrg_exposure_Spn$Spn
mrg_exposure_Spn$Spnexposures2 <- mrg_exposure_Spn$Spnexposures
mrg_exposure_Spn$Spnexposures2[mrg_exposure_Spn$Spnexposures>1] <- '2+'
mrg_exposure_Spn$LytA2 <- 0
mrg_exposure_Spn$LytA2[mrg_exposure_Spn$LytA >4] <- 1


m <- glmer(as.factor(AcquisitionLonger) ~ SpnNumber + AgeYears + as.factor(Virus) +  as.factor(Spnexposures2) +(1|Assay)+ (1|StudyID), data = mrg_exposure_Spn, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

summary(m)[[10]][,1:2]
OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Variable <- gsub('as.factor', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - OR_adjusted$Std..Error)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + OR_adjusted$Std..Error)
OR_adjusted$Pval <- paste0('p=', sapply(OR_adjusted$Pr...z.., FUN = "disp_pval"))
ggplot(OR_adjusted, aes(x=OR, y=Variable)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  geom_vline(xintercept = 1, linetype = 'dashed') + 
  xlab('Adjusted Odds Ratio LMM') 
ggsave('Fig2K_AdjustedOR_Spn_exposure.pdf', width = 3, height = 2)




save.image('figure2_exposures.Rdata')



