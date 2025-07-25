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


setwd("Z:\\Projects\\SAMSAM\\Manuscripts\\Figures\\Figure4_bacteria_viruses")


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
library(ggmosaic)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

#save the session info
sink("sessionInfo.txt")
sessionInfo()
sink()


#load in database
mrg <- data.frame(read_xlsx('Z:\\Projects\\SAMSAM\\IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20250722.xlsx'))

#filter out antibiotics samples
mrg1 <- mrg %>% 
  filter(!grepl("Yes", Antibiotics))   #remove anitbiotics samples

mrg2b <- mrg1 %>% #remove pathogens/child combinations where there is no infection
  group_by(StudyID, Assay) %>%
  filter(
    max(InfectionLongerTrue) > 0
  ) %>%
  ungroup()

##############################################
#
#Figure 4A
#
###########################################################
mrg2 <- mrg1[mrg1$AssayType == 'Viral', ] #select only the viral assays
mrg2_infected <- mrg2[mrg2$InfectionLongerTrue == 1,] #take only infected samples to look at the assay and kids


infections_num <- infections2 <- infections <- table(mrg2_infected$Assay, mrg2_infected$StudyID) #make a table with number of infections
infections2[infections>0] <- 1 #make a matrix with kids per sample

#add the one child with no viral infections
infections3 <- cbind(infections2, rep(0,13)) 
colnames(infections3)[45] <- 'SAM-K1-25'

infections_num3 <- cbind(infections_num, rep(0,13)) 
colnames(infections_num3)[45] <- 'SAM-K1-25'

#get the age of the kids
demo <- mrg[!duplicated(mrg$StudyID),c('StudyID',  'AgeYears', 'Gender')]
demo2 <- demo[match(colnames(infections_num3), demo$StudyID),] #order them in same way as heatmap

#make a heatmap of viral infections, make annotation for age and number of infections
columnha = HeatmapAnnotation(age = demo2$AgeYears, child = anno_barplot(colSums(infections3)),
                             col = list(age =  c("1" = "lightblue", "2" = "steelblue", 
                                                 "3" = "blue",  "4" = "darkblue", '5' = 'black')))
rowha = rowAnnotation(virus = anno_barplot(rowSums(infections3)))

pl1 <- Heatmap(infections_num3, 
               rect_gp = gpar(col = "white", lwd = .5),
               col = colorRamp2(c(0, 1, 14, 30), c("blue", "white", 'orange', "red")),
               top_annotation = columnha,
               right_annotation = rowha)

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig4A_heatmap_virus_days_", timestamp, ".pdf")
pdf(pdf_filename, width = 10, height = 6)
print(pl1)
dev.off()



###########################################################
#
#Fig 4B
#
############################################################
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
      lengthInfectionLonger <- each$Day2[clears] - each$Day2[start] #calculate the lengt of infections and add
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

Surv(lengthInfections2Longer$Duration, lengthInfections2Longer$status)
s1 <- survfit(Surv(Duration, status) ~ 1, data = lengthInfections2Longer)
str(s1)


pl3 <- survfit2(Surv(Duration, status) ~ 1, data = lengthInfections2Longer) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +   xlim(c(0,27)) 

summary(survfit2(Surv(Duration, status) ~ 1, data = lengthInfections2Longer))


pdf_filename <- paste0("Fig4B_survival_virus_days_", timestamp, ".pdf")
pdf(pdf_filename, width = 3, height = 3)
print(pl3)
dev.off()






###########################################################
#
#Fig 4C
#
############################################################
#check if viral infection is linked to any bacteria or number of viruses
mrg3 <- mrg2 %>%
  group_by(StudyID, Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLonger) == 1
  ) %>%
  ungroup() #make a new dataframe and only select assays that are present



mrg_clear <- mrg3
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
    sample_check2 <- data.frame(mrg2[mrg2$SampleID_Day == sample_check,]) #get the total dataset for this for only viral assays
    each4$VirusNumber[i] <- sum(sample_check2[, 'InfectionLongerTrue'])
    #fill in the covariates of interest
    
  }
  mrg_clear.split[[h]] <- each4 #save the fixed file back in the list
}

mrg_clear2 <- do.call(rbind, mrg_clear.split)
mrg_clear2 <- mrg_clear2[mrg_clear2$ClearanceLonger == 0,] #remove the clearancelonger point
table(mrg_clear2$Assay)
table(mrg_clear2$ClearNextDay, mrg_clear2$VirusNumber)

#add potential confounders like bacterial carriage
mrg_clear2$SpnNumber  <- NA
mrg_clear2$Staph <- NA
mrg_clear2$Haemophilus <- NA
mrg_clear2$LytA <- NA

for(i in 1:nrow(mrg_clear2)){
  
  sample_check <- mrg_clear2$SampleID_Day[i]
  sample_check2 <- mrg1[mrg1$SampleID_Day == sample_check,]
  
  #fill the values
  mrg_clear2$LytA[i] <- sample_check2[sample_check2$Assay == "Spn(lytA)",'log10_avgconc']
  mrg_clear2$Staph[i] <- sample_check2[sample_check2$Assay == "Staphylococcus_aureus",'InfectionLongerTrue']
  mrg_clear2$Haemophilus[i] <- sample_check2[sample_check2$Assay == "Haemophilus_influenzae",'InfectionLongerTrue']
  mrg_clear2$SpnNumber[i] <- sum(sample_check2[grepl('Spn',sample_check2$Assay) & 
                                                 !grepl('lytA',sample_check2$Assay) & 
                                                 !grepl('piaB',sample_check2$Assay) & 
                                                 !grepl('Spn6C/D',sample_check2$Assay), 'InfectionLongerTrue'])
}
mrg_clear2$Spn <- 0 
mrg_clear2$Spn[mrg_clear2$SpnNumber >0] <- 1 #make a Spn yes.no for viral comparison


m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = mrg_clear2, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)


OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - OR_adjusted$Std..Error)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + OR_adjusted$Std..Error)
ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_x_log10() 
ggsave('Fig4C_AdjustedORClearanceNextDay_virus.pdf', width = 3, height = 2)







###########################################################
#
#Fig 4D
#
############################################################
mrg_all <- mrg2 
#check for infection by any virus
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
mrg_all3$InfectionType[rowSums(mrg_all3[,2:17]) > 1] <- 'z.multiple_viruses'
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
mrg_all4 <- merge(mrg2b, mrg_all3, by= 'SampleID_Day', all.x=T)
mrg_all4$InfectionType[is.na(mrg_all4$InfectionType)] <- '0.None'
mrg_all4$InfectionTypeBroad[is.na(mrg_all4$InfectionTypeBroad)] <- '0.None'
mrg_all4$InfectionAnyVirus[is.na(mrg_all4$InfectionAnyVirus)] <- '0.None'

mrg_bacteria <- mrg_all4[ mrg_all4$Assay == 'Staphylococcus_aureus',]
m <- glmer(as.factor(InfectionLongerTrue) ~ as.factor(InfectionAnyVirus) + AgeYears + (1|StudyID), 
           data = mrg_bacteria,    family = binomial, control = glmerControl(optimizer = "bobyqa"))
anyvirus <- summary(m)[[10]][2,]
#anyvirus_conf <- confint(m)

mrg_bacteria2 <- mrg_bacteria[ mrg_bacteria$InfectionTypeBroad != 'PIV' &
                                 mrg_bacteria$InfectionTypeBroad != 'adenovirus'  ,]
m1<- glmer(as.factor(InfectionLongerTrue) ~ as.factor(InfectionTypeBroad) + AgeYears+ (1|StudyID), 
           data = mrg_bacteria2,    family = binomial, control = glmerControl(optimizer = "bobyqa"))
spec <- summary(m1)[[10]][2:6,]

#########these tables are part of the figure
table(mrg_bacteria$InfectionTypeBroad, mrg_bacteria$InfectionLongerTrue)
table(mrg_bacteria$InfectionAnyVirus, mrg_bacteria$InfectionLongerTrue)


prop.table(table(mrg_bacteria$InfectionTypeBroad, mrg_bacteria$InfectionLongerTrue), margin = 1)
prop.table(table(mrg_bacteria$InfectionAnyVirus, mrg_bacteria$InfectionLongerTrue), margin = 1)


OR_adjusted <- data.frame(rbind(anyvirus, spec))
OR_adjusted$Virus <- gsub('as.factor[(]InfectionTypeBroad[)]', '', rownames(OR_adjusted))
OR_adjusted$Perc <- exp(OR_adjusted$Estimate)
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - OR_adjusted$Std..Error)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + OR_adjusted$Std..Error)
OR_adjusted$Virus <- factor(OR_adjusted$Virus, levels = rev(c('anyvirus', 'entero_rhino_virus',
                                                              'bocavirus', 'z.multiple_viruses',
                                                              'coronavirus', 'influenza',
                                                              'adenovirus', 'PIV')))
ggplot(OR_adjusted, aes(x=Perc, y=Virus)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('Adjusted OR of staph aureus detection') + 
  geom_vline(xintercept = 1, linetype = 'dashed')
ggsave('Fig4D_OR of staph aureus detection.pdf', width = 3, height = 2)





###########################################################
#
#Fig 4E
#
############################################################
targets0 <- c('Haemophilus_influenzae', 'Spn(lytA)', 'Staphylococcus_aureus',
              'Moraxella_catarrhalis', 'Streptococcus_pyogenes')


counter  = 0 #set a counter that is used to merge
for(target in targets0){ #cycle through the bacteria
  
  counter  = counter+1
  
  #select the bacterial targets
  #mrg_all5 <- mrg_all4[grepl(target, mrg_all4$Assay),]
  
  #select the bacterial target but only when detected, so remove the non-detected ones
  mrg_all5 <- mrg_all4[ mrg_all4$Assay == target &
                         mrg_all4$Average_conc > 0,]
  
  #creat the tally to display
  mrg_all5_summary <- data.frame(mrg_all5 %>% group_by(InfectionTypeBroad) %>% tally())
  mrg_all5_summary$Assay <- target
    
  #perform stats using LMM, correct for age as it matters  for density and use child id as random effect
  lmm <- lmer(log10_avgconc ~ as.factor(InfectionTypeBroad) + AgeYears + (1|StudyID) , data = mrg_all5)
  mrg_all5_inf <- data.frame(lmerTest:::get_coefmat(lmm)[2:(nrow(lmerTest:::get_coefmat(lmm))-1),c(1,5)])
  

  #fix the names etc
  mrg_all5_inf$InfectionTypeBroad <- gsub('as.factor[(]InfectionTypeBroad[)]', '', rownames(mrg_all5_inf))
  colnames(mrg_all5_inf)[2] <- 'Pval'
  mrg_all5_inf$Pdisp <- mrg_all5_inf$Pval #creat a column with a p-value for displating in graphs
  for(r in 1:nrow(mrg_all5_inf)){mrg_all5_inf$Pdisp[r] <- disp_pval(mrg_all5_inf$Pval[r])} #make it correct for disdplay
  mrg_all5_inf$Assay <- target
  
  
  #combine into one large dataframe for all. If counter is 1 the new frames are created, otherwise they are bound
  if(counter == 1){
    mrg_all6 <- mrg_all5
    mrg_all6_summary <- mrg_all5_summary
    mrg_all6_inf <- mrg_all5_inf
      } else {
    mrg_all6 <- rbind(mrg_all6, mrg_all5)
    mrg_all6_summary <- rbind(mrg_all6_summary, mrg_all5_summary)
    mrg_all6_inf <- rbind(mrg_all6_inf, mrg_all5_inf)
  }
}


#create a merger that is used to remove the assays with not enough observations   
mrg_all6_summary$MergeID <- paste(mrg_all6_summary$Assay, mrg_all6_summary$InfectionTypeBroad, sep='_')
mrg_all6_inf$MergeID <- paste(mrg_all6_inf$Assay, mrg_all6_inf$InfectionTypeBroad, sep='_')
mrg_all6$MergeID <- paste(mrg_all6$Assay, mrg_all6$InfectionTypeBroad, sep='_')

#keep only combinations with at least 10 observations
mrg_all7_summary <- mrg_all6_summary[mrg_all6_summary$n > 9,]
mrg_all7_inf <- mrg_all6_inf[mrg_all6_inf$MergeID %in% mrg_all7_summary$MergeID,]
mrg_all7 <- mrg_all6[mrg_all6$MergeID %in% mrg_all7_summary$MergeID,]

#make an adjusted p-value as we are looking at multiple bacteria and viruses combination
mrg_all7_inf$Padj <- p.adjust(mrg_all7_inf$Pval, method = 'BH')
mrg_all7_inf$Pdisp_adj <- mrg_all7_inf$Padj #creat a column with a p-value for displatying it in graphs
for(r in 1:nrow(mrg_all7_inf)){mrg_all7_inf$Pdisp_adj[r] <- disp_pval(mrg_all7_inf$Padj[r])}
mrg_all7_inf$Pdisp_adj2 <- ''
mrg_all7_inf$Pdisp_adj2[mrg_all7_inf$Padj<0.05] <- '*'
mrg_all7_inf$Pdisp_adj2[mrg_all7_inf$Padj<0.01] <- '**'
mrg_all7_inf$Pdisp_adj2[mrg_all7_inf$Padj<0.001] <- '***'
  
#set the factor level to have them in the right order
mrg_all7$Assay <- factor(mrg_all7$Assay, levels= c( 'Spn(lytA)', 'Haemophilus_influenzae', 
                         'Moraxella_catarrhalis', 'Staphylococcus_aureus', 'Streptococcus_pyogenes'))


#make the plot
    pl33 <- ggplot(mrg_all7, aes(x=InfectionTypeBroad, y=log10_avgconc)) + 
    geom_violin(aes(fill = InfectionTypeBroad)) + geom_boxplot(width = 0.2)+
    xlab('') +
    facet_grid(.~Assay, scale = 'free_x', space = 'free', ) + 
    geom_text(data= mrg_all7_summary, aes(InfectionTypeBroad, Inf, label=n), vjust = 1) +
      ggtitle(paste0('Infection virus effect on bacterial density')) + 
      geom_text(data= mrg_all7_inf, aes(InfectionTypeBroad, Inf, label=Pdisp_adj2), vjust = 2.5)
    
#save it
timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig4E_virus_bacterial_density_",
                       timestamp, ".pdf")
pdf(pdf_filename, width = 10,height = 3)
print(pl33)
dev.off()





########################################################################
#
#Figure 4F
#check relative to virus onset 
#######################################################################
viruses2 = c('HRV')
targets1 <- c('Haemophilus_influenzae', 'Spn(lytA)', 'Moraxella_catarrhalis', 'Staphylococcus_aureus')

#initialize a heatmap for timestats, the ncol needs to match the number of timepoints tested
timestat_estClear <- timestat_estAcquire <- matrix(ncol=15, nrow=length(targets1))
timestat_pClear <-timestat_pAcquire <- matrix(ncol=15, nrow=length(targets1))



for(virus in viruses2){
r2 <- 0

for(target in targets1){

r2=r2+1 #add to counter

#se;ect for the target and normalize for day -1 to -3 relative to viral infection onset
mrg7.split <- split(mrg2b, mrg2b$StudyID)
for(i in 1:length(mrg7.split)){
each <- mrg7.split[[i]] #cycle through the childre
dates <- each[grepl(virus, each$Assay) & each$AcquisitionLonger == 1,]$Day2 #select the day of viral acquisition
each$DayToVirus <- each$Day2 - dates[1] #creat a new time column that normalizes to the FIRST time a virus is acquired
each2 <-  each[each$Assay == target,] #select only bacterial pathogen
base <- mean(each2$log10_avgconc[which(each2$DayToVirus == -3 | each2$DayToVirus == -2 | each2$DayToVirus == -1)], na.rm=T) #make a baseline average conc
each2$FC_dens <- each2$log10_avgconc - base #make a fold change relative to the baseline
mrg7.split[[i]] <- each2 #save back in the list
}
mrg8 <- do.call(rbind, mrg7.split)

mrg9 <- mrg8[mrg8$DayToVirus > -6 & mrg8$DayToVirus < 15, ] #filter out smaples between x days before or after viral onset
mrg9$DayToVirus2 <- mrg9$DayToVirus #make a new column and change the 5 days prior to onset to the baseline
mrg9$DayToVirus2[mrg9$DayToVirus2 %in% c(-5:-1)] <- '-1_baseline'


lm <- lmer(log10_avgconc ~ as.factor(DayToVirus2) +(1|StudyID), data = mrg9) #do stats for all days relative to baseline
acq_res <- lmerTest:::get_coefmat(lm)[2:nrow(lmerTest:::get_coefmat(lm)),c(1,5)] #extract the results from the linear model



#save the dataframe
assign(paste0(target,'_', virus,'_acqdataset'), mrg9)

#do the same as above but relative to viral clearance
mrg2c <- mrg2b[grepl(virus, mrg2b$Assay),] 
keep <- unique(mrg2c$StudyID[mrg2c$ClearanceLonger == 1])
mrg2d <- mrg2b[mrg2b$StudyID %in% keep,]
mrg7.split <- split(mrg2d, mrg2d$StudyID)
for(i in 1:length(mrg7.split)){
  each <- mrg7.split[[i]]
  dates <- each[grepl(virus, each$Assay) & each$ClearanceLonger == 1,]$Day2
  each$DayToVirus <- each$Day2 - dates[length(dates)] #after final clearance
  each2 <-  each[each$Assay == target,] #select only bacterial pathogen
  base <- mean(each2$log10_avgconc[which(each2$DayToVirus == -3 | each2$DayToVirus == -2 | each2$DayToVirus == -1)], na.rm=T)
  each2$FC_dens <- each2$log10_avgconc - base
  mrg7.split[[i]] <- each2
}
mrg10 <- do.call(rbind, mrg7.split)

mrg11 <- mrg10[mrg10$DayToVirus > -6 & mrg10$DayToVirus < 15, ]
mrg11$DayToVirus2 <- mrg11$DayToVirus
mrg11$DayToVirus2[mrg11$DayToVirus2 %in% c(-5:-1)] <- '-1_baseline'


lm <- lmer(log10_avgconc ~ as.factor(DayToVirus2) +(1|StudyID), data = mrg11)
lm_clear <- lmerTest:::get_coefmat(lm)[2:nrow(lmerTest:::get_coefmat(lm)),c(1,5)]

#add the stats to the time heatmaps per bacteria
timestat_estClear[r2,] <- lm_clear[,1]
timestat_pClear[r2,] <- lm_clear[,2]
timestat_estAcquire[r2,] <- acq_res[,1]
timestat_pAcquire[r2,] <- acq_res[,2]


}
rownames(timestat_estClear) <- rownames(timestat_pClear) <- rownames(timestat_estAcquire) <- rownames(timestat_pAcquire) <-targets1
colnames(timestat_estClear) <- colnames(timestat_pClear) <- colnames(timestat_estAcquire) <- colnames(timestat_pAcquire) <-gsub('as.factor(DayToVirus2)', '', rownames(lm_clear), fixed = T)


#initialize a stats for time effect, the ncol needs to match the number of timepoints tested
assign(x= paste0(virus, '_timestats_clear_est'), value = timestat_estClear) 
assign(x= paste0(virus, '_timestats_clear_pval'), value = timestat_pClear)
assign(x= paste0(virus, '_timestats_acquire_est'), value = timestat_estAcquire) 
assign(x= paste0(virus, '_timestats_acquire_pval'), value = timestat_pAcquire)


}

#heatmaps for time per virus
HRV_timestats_acquire_est2 <- HRV_timestats_acquire_est[,order(as.numeric(colnames(HRV_timestats_acquire_est)))]
HRV_timestats_acquire_pval3 <-  HRV_timestats_acquire_pval2 <- HRV_timestats_acquire_pval[,order(as.numeric(colnames(HRV_timestats_acquire_pval)))]*4 #multiply by 4 to correct for 4 bacteria
HRV_timestats_acquire_pval3[HRV_timestats_acquire_pval2 >0.05] <- ''
HRV_timestats_acquire_pval3[HRV_timestats_acquire_pval2 <0.05] <- '*'
HRV_timestats_acquire_pval3[HRV_timestats_acquire_pval2 <0.01] <- '**'
HRV_timestats_acquire_pval3[HRV_timestats_acquire_pval2 <0.001] <- '***'
pl6 <- pheatmap::pheatmap(HRV_timestats_acquire_est2, cluster_cols = F, main = paste0('HRV: days since acquisition'),
                          display_numbers = HRV_timestats_acquire_pval3, 
                          fontsize_number = 15 )

pdf('Fig4F_heatmap_time_virus_bacteria_effect.pdf', width = 10, height = 8,onefile=T)
print(pl6)
dev.off()




########################################################
#
#Figure 4G
#
########################################################33
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
ggplot(mrg_HRV2, aes(x = DayToVirus, y = log10_avgconc, colour = StudyID)) +  # make a plot.
  geom_line(aes(group = StudyID), alpha=0.5)+
  theme(
    legend.position = "right"
  ) +  # Rotate x-axis labels
  labs(
    title = 'HRV over time',  # Title for each StudyID
    x = 'Days since viral onset', y = "(log10) density"  # y-axis label with log scale
  )  + xlim(c(-3,14)) +
  stat_summary(fun = "mean", colour = "red", size = 2, geom = "line", na.rm=T) + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  theme(legend.position = 'none')
filename <- paste0('Fig4G_HRV_density_over_time.pdf')
ggsave(filename, width = 4, height = 4)

table(mrg_HRV2[mrg_HRV2$DayToVirus ==14,]$InfectionLongerTrue)




###########################################
#
#Fig 4H
#
##############################################
#combine all datasets for HRV
All_HRV <- rbind(Haemophilus_influenzae_HRV_acqdataset, `Spn(lytA)_HRV_acqdataset`, 
                Moraxella_catarrhalis_HRV_acqdataset, Staphylococcus_aureus_HRV_acqdataset)
All_HRV <- All_HRV[!is.na(All_HRV$Assay),]
ggplot(All_HRV, aes(x = DayToVirus, y = FC_dens, colour = Assay)) +  # make a plot.
  geom_vline(xintercept=0, linetype='dashed')+
  geom_hline(yintercept=0, linetype='dashed')+
  geom_line(aes(group = StudyID), alpha=0.5)+
  geom_smooth(aes(color = Assay, fill=Assay), method = "loess", se = T, size = 0.7, span = .6) +    # Add separate smoothing line for each Assay
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "right"
  ) +  # Rotate x-axis labels
  labs(
    title = 'HRV time effect',  # Title for each StudyID
    x = 'Days since viral onset', y = "FC (log10)"  # y-axis label with log scale
  )  + xlim(c(-5,14)) + 
  facet_grid(.~Assay)
filename <- paste0('Fig4H_HRV_bacteria_timeefeect.pdf')
ggsave(filename, width = 10, height = 4)




##############################################################################
save.image('Figure4_bacteria_virus.Rdata')
