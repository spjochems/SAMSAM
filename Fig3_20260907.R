#####
#This script uses the mlvar package to calculate correlations between bacteria
#https://arxiv.org/pdf/1609.04156
#https://cran.r-project.org/web/packages/mlVAR/mlVAR.pdf
#Spn (lytA) / Haemophilus / Moraxella / Staph 
#it uses the RAW (log-transformed) data or with data normalized to 16S 


rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\Fig3\\")

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
library(ComplexHeatmap)
library(survival)
library(ggsurvfit)
library(ppcor)
library(broom.mixed)
library(parallel)

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#load in database
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))
mrg <- mrg[mrg$Antibiotics != 'Yes',] #remove the antibiotics samples as it may affect results

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



###################################################
#                                                 #
#                 FIGURE 3A                       #
#                                                 #
###################################################

mrg9 <- mrg %>%
  group_by(StudyID, Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLongerTrue) ==1 
  ) %>%
  ungroup()

notspn2 <- mrg9[grepl("^Spn\\(lytA\\)$", mrg9$Assay) | !grepl("^Spn", mrg9$Assay), ]
notspn2 <- notspn2 %>% filter(Assay != "EAV")  %>% filter(Assay != "16S")
notspn3 <- notspn2 
notspn3$Assay <- factor(notspn3$Assay,
                        levels = c("Haemophilus_influenzae", "Moraxella_catarrhalis", "Spn(lytA)", "Staphylococcus_aureus", 
                                   "Streptococcus_pyogenes", "HRV", "bocavirus", "coronavirus_NL63"))
notspn3 <- notspn3 %>% 
  filter(StudyID %in% c("SAM-K1-17", "SAM-K1-21", "SAM-K1-41"))

ggplot(notspn3, aes(x = Day2, y = log10_avgconc)) +  # Convert Day to numeric for regression
  geom_point(aes(color = Assay, shape = AssayType), size = 1) +  # Points colored by Assay
  #geom_line(aes(color = Assay), size = 0.7) +  # Line connecting the points
  geom_smooth(aes(color = Assay, linetype = AssayType), method = "loess", se = FALSE, size = 0.7, span = .6) +    # Add separate smoothing line for each Assay
  theme_minimal() +
  facet_wrap(.~StudyID) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "right"
  ) +  # Rotate x-axis labels
  labs(
    title = 'Bacterial/Viral carriage',  # Title for each StudyID
    x = NULL, y = "log10(average_conc)"  # y-axis label with log scale
  )
ggsave("FIG3A.pdf", width = 4.8, height = 1.7)



###################################################
#                                                 #
#                 FIGURE 3B                       #
#                                                 #
###################################################

mrg1 <- mrg

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
pl1

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig3B_heatmap_virus_days_", timestamp, ".pdf")
pdf(pdf_filename, width = 10, height = 6)
print(pl1)
dev.off()

mean(colSums(infections3)) #average number of infections per child


###################################################
#                                                 #
#                 FIGURE 3C                       #
#                                                 #
###################################################

mrg1 <- mrg %>% 
  filter(!grepl("Yes", Antibiotics))   #remove anitbiotics samples

mrg2b <- mrg1 %>% #remove pathogens/child combinations where there is no infection
  group_by(StudyID, Assay) %>%
  filter(
    max(InfectionLongerTrue) > 0
  ) %>%
  ungroup()

mrg2 <- mrg1[mrg1$AssayType == 'Viral', ] #select only the viral assays
mrg2_infected <- mrg2[mrg2$InfectionLongerTrue == 1,] #take only infected samples to look at the assay and kids

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

targets0 <- c('Haemophilus_influenzae', 'Spn(lytA)', 'Moraxella_catarrhalis')

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
Fig3C <- ggplot(mrg_all7, aes(x=InfectionTypeBroad, y=log10_avgconc)) + 
  geom_violin(aes(fill = InfectionTypeBroad)) + geom_boxplot(width = 0.2)+
  xlab('') +
  facet_grid(.~Assay, scale = 'free_x', space = 'free', ) + 
  geom_text(data= mrg_all7_summary, aes(InfectionTypeBroad, Inf, label=n), vjust = 1) +
  ggtitle(paste0('Infection virus effect on bacterial density')) + 
  geom_text(data= mrg_all7_inf, aes(InfectionTypeBroad, Inf, label=Pdisp_adj2), vjust = 2.5)
Fig3C

#save it
timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig3C_virus_bacterial_density_",
                       timestamp, ".pdf")
pdf(pdf_filename, width = 8,height = 3)
Fig3C
dev.off()


###################################################
#                                                 #
#                 FIGURE 3D                       #
#                                                 #
###################################################

viruses2 = c('HRV')
targets1 <- c('Haemophilus_influenzae', 'Spn(lytA)', 'Moraxella_catarrhalis', 'Staphylococcus_aureus')

#initialize a heatmap for timestats, the ncol needs to match the number of timepoints tested
timestat_estClear <- timestat_estAcquire <- matrix(ncol=15, nrow=length(targets1))
timestat_pClear <-timestat_pAcquire <- matrix(ncol=15, nrow=length(targets1))

for(virus in viruses2){
  r2 <- 0
  
  for(target in targets1){
    
    r2=r2+1 #add to counter
    
    #select for the target and normalize for day -1 to -3 relative to viral infection onset
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

All_HRV <- rbind(Haemophilus_influenzae_HRV_acqdataset, `Spn(lytA)_HRV_acqdataset`, 
                 Moraxella_catarrhalis_HRV_acqdataset, Staphylococcus_aureus_HRV_acqdataset)
All_HRV <- All_HRV[!is.na(All_HRV$Assay),]

#reorder:
All_HRV <- All_HRV %>%
  mutate(
    Assay = factor(
      Assay,
      levels = c(
        "Haemophilus_influenzae",
        "Spn(lytA)",          # check exact spelling in your data
        "Staphylococcus_aureus",
        "Moraxella_catarrhalis"
      )
    )
  )

ggplot(All_HRV, aes(x = DayToVirus, y = FC_dens, colour = Assay)) +  # make a plot.
  geom_vline(xintercept=0, linetype='dashed')+
  geom_hline(yintercept=0, linetype='dashed')+
  geom_line(aes(group = StudyID), alpha=0.5)+
  geom_smooth(aes(color = Assay, fill=Assay), method = "loess", se = T, size = 0.7, span = .6) +    # Add separate smoothing line for each Assay
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0), # Rotate x-axis labels
    legend.position = "none", strip.background = element_rect(fill = NA, colour = NA)) +  
  labs(
    title = 'HRV time effect',  # Title for each StudyID
    x = 'Days since HRV onset', y = "Fold-change density (log10)" ) + # y-axis label with log scale
    xlim(c(-5,14)) + 
  facet_wrap(.~Assay, ncol = 2)

filename <- paste0('Fig3D_HRV_bacteria_timeefeect.pdf')
ggsave(filename, width = 3.5, height = 5)


#nr samples:
fig_data      <- All_HRV
group_cols    <- c("Assay", "DayToVirus")  
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig3D_samplesize.txt")




###################################################
#                                                 #
#                 FIGURE 3E                       #
#                                                 #
###################################################

#heatmaps for time per virus
HRV_timestats_acquire_est2 <- HRV_timestats_acquire_est[,order(as.numeric(colnames(HRV_timestats_acquire_est)))]
HRV_timestats_acquire_pval3 <-  HRV_timestats_acquire_pval2 <- HRV_timestats_acquire_pval[,order(as.numeric(colnames(HRV_timestats_acquire_pval)))]*4 #multiply by 4 to correct for 4 bacteria
HRV_timestats_acquire_pval3[HRV_timestats_acquire_pval2 >0.05] <- ''
HRV_timestats_acquire_pval3[HRV_timestats_acquire_pval2 <0.05] <- '*'
HRV_timestats_acquire_pval3[HRV_timestats_acquire_pval2 <0.01] <- '**'
HRV_timestats_acquire_pval3[HRV_timestats_acquire_pval2 <0.001] <- '***'

Fig3E <- pheatmap::pheatmap(HRV_timestats_acquire_est2, cluster_cols = F, main = paste0('HRV: days since acquisition'),
                          display_numbers = HRV_timestats_acquire_pval3, 
                          fontsize_number = 15 )
Fig3E

pdf('Fig4E_heatmap_time_virus_bacteria_effect.pdf', width = 8, height = 3,onefile=T)
Fig3E
dev.off()


###################################################
#                                                 #
#                 FIGURE 3F                       #
#                                                 #
###################################################

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

filename <- paste0('Fig3F_HRV_density_over_time.pdf')
ggsave(filename, width = 2, height = 2)

table(mrg_HRV2[mrg_HRV2$DayToVirus ==14,]$InfectionLongerTrue)

#nr samples:
fig_data      <- mrg_HRV2
group_cols    <- c("DayToVirus")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig3F_samplesize.txt")






###################################################
#                                                 #
#                 FIGURE 3G                       #
#                                                 #
###################################################

mycols <- c("StudyID",
            "Day",
            "Assay",
            "log10_avgconc")

mrg2 <- mrg 

#remove 18PS,  EAV and 16S as assays
mrg2 <- mrg2[mrg2$Assay == 'Spn(lytA)' | 
               grepl("Haemophilus_influenzae", mrg2$Assay)| 
               grepl("aureus", mrg2$Assay)|
               grepl("Moraxella_catarrhalis", mrg2$Assay),]

mrg3 <- mrg2 %>% dplyr::select(all_of(mycols)) 
#mrg3$Day <- factor(mrg3$Day, levels = unique(mrg3$Day)) 
mrg3$Day <- as.numeric(gsub('D', '', mrg3$Day)) 
mrg3 <- mrg3[!is.na(mrg3$Day),]

mrg4 <- mrg3 %>%
  pivot_wider(names_from = Assay, 
              values_from = c(log10_avgconc))


colnames(mrg4)[colnames(mrg4) == "Spn(lytA)"] <- 'Spn'

cor.test(mrg4$Spn, mrg4$Haemophilus_influenzae)  
cor.test(mrg4$Spn, mrg4$Staphylococcus_aureus)  
cor.test(mrg4$Spn, mrg4$Moraxella_catarrhalis) 

###########do a partial correlation analysis
partialcor <- pcor(mrg4[,3:6])
rounded <- round(partialcor$estimate, 2)

Fig3G <- ComplexHeatmap::pheatmap(
  partialcor$estimate,
  display_numbers = rounded,
  fontsize_number = 12,
  heatmap_legend_param = list(
    title = "Partial\ncorrelation\nanalysis",
    title_position = "topcenter"
  )
)
Fig3G

#save:
pdf('Fig3G_partialcorrelation_raw.pdf', width = 5, height = 4)
Fig3G
dev.off()



###################################################
#                                                 #
#                 FIGURE 3H&I                     #
#                                                 #
###################################################

VAR <- mlVAR(mrg4, vars = c('Spn', 'Haemophilus_influenzae',  'Staphylococcus_aureus', 'Moraxella_catarrhalis'), 
             idvar = 'StudyID', lags = 1, contemporaneous = 'correlated',
             estimator = 'lmer',)

print(VAR)
summary(VAR)

par(mfrow=c(1,3))
plot(VAR, posCol = 'Red', negCol ='Blue', title = 'Temporal - raw', label.cex = 1, label.scale = F)
plot(VAR, type = 'contemporaneous', posCol = 'Red', negCol ='Blue', title = 'Contemporaneous - raw', label.cex = 1, label.scale = F)
plot(VAR, type = 'between', posCol = 'Red', negCol ='Blue', title = 'Between Subject - raw', label.cex = 1, label.scale = F)

#Save pcor values:
par(mfrow=c(1,1))
Fig3H <- plot(VAR, type = "contemporaneous",
     posCol = "red", negCol = "blue",
     title = "Contemporaneous - raw",
     label.cex = 1, label.scale = FALSE,
     edge.labels = TRUE)
Fig3H

pdf('Fig3H_VAR_plots_bacteria_rawdata_lag1.pdf', width = 3, height = 3)
Fig3H
dev.off()

Fig3I <- plot(VAR, type = "between",
     posCol = "red", negCol = "blue",
     title = "Between Subject - raw",
     label.cex = 1, label.scale = FALSE,
     edge.labels = TRUE)
Fig3I

pdf('Fig3I_VAR_plots_bacteria_rawdata_lag1.pdf', width = 3, height = 3)
Fig3I
dev.off()



###################################################
#                                                 #
#                 FIGURE 3J                       #
#                                                 #
###################################################

mrg0 <- mrg

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

mrg3 <- mrg6.norm %>%
  pivot_longer(
    -c(StudyID, Day, '16S'),  # keep these fixed columns
    names_to = "Assay",
    values_to = "log10_avgconc"
  )

for(select_assay in c('Spn', 'Staphylococcus_aureus', 
                      'Haemophilus_influenzae', 'Moraxella_catarrhalis')){
  
  #select the assay of interest
  
  lyta_tab2 <- data.frame(mrg3[mrg3$Assay == select_assay,])
  lyta_tab2 <- lyta_tab2[,-which(colnames(lyta_tab2) == 'Assay')]
  
  #make the dataframe into a format to wrangle the autocorrleations
  lytawid <- data.frame(pivot_wider(lyta_tab2, 
                                    names_from = Day, 
                                    values_from = log10_avgconc))
  rownames(lytawid) <- lytawid[,1]
  lytawid2 <- lytawid[,-1]
  colnames(lytawid2) <- gsub( 'X', '',  colnames(lytawid2))
  lytawid3 <- lytawid2[,order(as.numeric(colnames(lytawid2)))]
  lytawid4 <- lytawid3[,1:29]
  
  #make a lagging list up to 21 days
  lags <- list()
  for(j in 1:21){
    lag1 <- lytawid4[,c(1,(1+j))]
    colnames(lag1) <- c('t0', 't1')
    lag1$child <- rownames(lytawid4)
    for(i in 2:28){
      if((i+j)<30){
        lag1.0 <-  lytawid4[,c(i,(i+j))]
        colnames(lag1.0) <- c('t0', 't1')
        lag1.0$child <- rownames(lytawid4)
        lag1 <- rbind(lag1, lag1.0)
      }
    }
    lags[[j]] <- lag1
  }
  
  
  #mixed effect model per childr
  all_lmer <- lapply(lags, function(i){
    lmm <- lmer(t1 ~ t0 + (1 | child), data=i)
    summary(lmm)[[10]][2,]
  })
  all_lmer2 <- data.frame(do.call(rbind, all_lmer))
  all_lmer2$lag <- 1:21
  all_lmer2$Assay <- select_assay
  all_lmer2 <- all_lmer2[1:10,]
  assign(paste0(select_assay, '_totalLMER'), all_lmer2)
  
  
  t1 <- lags[[1]]
  t3 <- lags[[3]]
  t9 <- lags[[9]]
  p1 <- ggplot(t1, aes(x=t0, y=t1)) + geom_point(alpha=0.8)  + ggtitle('Lag 1 day') + 
    theme(aspect.ratio=1)
  p2 <- ggplot(t3, aes(x=t0, y=t1)) + geom_point(alpha=0.8) + ylab('t3')  + ggtitle('Lag 3 days')+ 
    theme(aspect.ratio=1)
  p3 <- ggplot(t9, aes(x=t0, y=t1)) + geom_point(alpha=0.8) + ylab('t9') + ggtitle('Lag 9 days')+ 
    theme(aspect.ratio=1)
  #filename <- paste0('Fig3J_', select_assay,'autocorrelation_examples_16Normalised.pdf')
  #pdf(filename, width = 9, height = 3)
  print(plot_grid(p1,p2,p3, ncol=3))
  #dev.off()
  
  
  
  ###############look for stability per age
  ages <- data.frame(mrg2[!duplicated(mrg2$StudyID),c(2,7)])
  ages_list <- split(ages, ages$AgeYears)
  
  
  #########do this fper age
  k=0
  age_lmer <- list()
  for(age in ages_list){
    k=k+1
    
    lytawid5 <- lytawid4[rownames(lytawid4) %in% age$StudyID,]
    
    lags <- list()
    for(j in 1:21){
      lag1 <- lytawid5[,c(1,(1+j))]
      colnames(lag1) <- c('t0', 't1')
      lag1$child <- rownames(lytawid5)
      for(i in 2:28){
        if((i+j)<30){
          lag1.0 <-  lytawid5[,c(i,(i+j))]
          colnames(lag1.0) <- c('t0', 't1')
          lag1.0$child <- rownames(lytawid5)
          lag1 <- rbind(lag1, lag1.0)
        }
      }
      lags[[j]] <- lag1
    }
    
    #mixed effect model per childr
    all_lmer <- lapply(lags, function(i){
      lmm <- lmer(t1 ~ t0 + (1 | child), data=i)
      summary(lmm)[[10]][2,]
    })
    all_lmer2 <- data.frame(do.call(rbind, all_lmer))
    all_lmer2$lag <- 1:21
    all_lmer2 <- all_lmer2[1:14,]
    all_lmer2$Age <- k
    age_lmer[[k]] <- all_lmer2
  }
  
  age_lmers <- do.call(rbind, age_lmer)
  
  age_lmers <- age_lmers[age_lmers$lag < 11,] 
  age_lmers$Assay <- select_assay
  colnames(age_lmers)[5] <- 'Pval'
  age_lmers$Padj <- p.adjust(age_lmers$Pval, method='BH')
  assign(paste0(select_assay, '_agedLMER'), age_lmers)
  
  
  mycols3 <- c("Estimate",
               "lag",
               "Age")
  mycols4 <- c("Padj",
               "lag",
               "Age")
  age_lmers2 <- age_lmers %>% dplyr::select(all_of(mycols3)) 
  age_lmers3 <- age_lmers2 %>%
    pivot_wider(names_from = lag, 
                values_from = c(Estimate))
  age_lmers4 <- data.frame(age_lmers3[,2:ncol(age_lmers3)])
  colnames(age_lmers4)<- gsub('X', '', colnames(age_lmers4))
  
  age_lmersp <- age_lmers %>% dplyr::select(all_of(mycols4)) 
  age_lmersp2 <- age_lmersp %>%
    pivot_wider(names_from = lag, 
                values_from = c(Padj))
  age_lmersp4 <- age_lmersp3 <- data.frame(age_lmersp2[,2:ncol(age_lmersp2)])
  age_lmersp4[age_lmersp3>0.05] <- ''
  age_lmersp4[age_lmersp3<0.05] <- '*'
  age_lmersp4[age_lmersp3<0.01] <- '**'
  age_lmersp4[age_lmersp3<0.001] <- '***'
  
  #filename <- paste0('Fig3J_', select_assay,'_autocorrelation_lmm_perage_HM_16normalised.pdf')
  #pdf(filename, width = 6, height =3)
  pl2 <- pheatmap(age_lmers4, cluster_rows = F, cluster_cols=F, 
                  display_numbers = age_lmersp4, main = select_assay)
  print(pl2)
  
  #dev.off()
}

all_LMER <- rbind(Moraxella_catarrhalis_totalLMER, 
                  Haemophilus_influenzae_totalLMER,
                  Staphylococcus_aureus_totalLMER,
                  Spn_totalLMER)

colnames(all_LMER)[5] <- 'Pval'
all_LMER$Padj <- p.adjust(all_LMER$Pval, method='BH')
mycols3 <- c("Estimate",
             "lag",
             "Assay")
mycols4 <- c("Padj",
             "lag",
             "Assay")
all_LMER2 <- all_LMER %>% dplyr::select(all_of(mycols3)) 
all_LMER3 <- all_LMER2 %>%
  pivot_wider(names_from = lag, 
              values_from = c(Estimate))
all_LMER4 <- data.frame(all_LMER3[,2:ncol(all_LMER3)])
rownames(all_LMER4) <- data.frame(all_LMER3)[,1]
colnames(all_LMER4)<- gsub('X', '', colnames(all_LMER4))


all_LMERp <- all_LMER %>% dplyr::select(all_of(mycols4)) 
all_LMERp2 <- all_LMERp %>%
  pivot_wider(names_from = lag, 
              values_from = c(Padj))
all_LMERp4 <- all_LMERp3 <- data.frame(all_LMERp2[,2:ncol(all_LMERp2)])
all_LMERp4[all_LMERp3>0.05] <- ''
all_LMERp4[all_LMERp3<0.05] <- '*'
all_LMERp4[all_LMERp3<0.01] <- '**'
all_LMERp4[all_LMERp3<0.001] <- '***'

#plot:
Fig3J <- pheatmap(
  all_LMER4,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  display_numbers = all_LMERp4,
  show_rownames = TRUE,
  fontsize_number = 20, 
  angle_col = 0
)
Fig3J

#save:
filename <- 'Fig3J_all_assays_lmm_perassay_HM.pdf'
pdf(filename, width = 6, height =3)
Fig3J
dev.off()

#nr samples:
sum(mrg$Assay == "Staphylococcus_aureus" &
    mrg$log10_avgconc > 0,
  na.rm = TRUE)





###################################################
#                                                 #
#                 FIGURE 3K                       #
#                                                 #
###################################################

#take raw concentration, not 16S normalised
mycols <- c("StudyID",
            "Day",
            "Assay",
            "log10_avgconc")



mrg2 <- mrg 


#remove 18PS,  EAV and 16S as assays
mrg2 <- mrg2[grepl("Haemophilus_influenzae", mrg2$Assay),]

mrg3 <- mrg2 %>% dplyr::select(all_of(mycols)) 
#mrg3$Day <- factor(mrg3$Day, levels = unique(mrg3$Day)) 
mrg3$Day <- as.numeric(gsub('D', '', mrg3$Day)) 
mrg3 <- mrg3[!is.na(mrg3$Day),]

mrg3$Assay[mrg3$Assay == 'Spn(lytA)'] <- 'Spn'

for(select_assay in c('Haemophilus_influenzae')){
  
  #select the assay of interest
  
  lyta_tab2 <- data.frame(mrg3[mrg3$Assay == select_assay,])
  lyta_tab2 <- lyta_tab2[,-which(colnames(lyta_tab2) == 'Assay')]
  
  #make the dataframe into a format to wrangle the autocorrleations
  lytawid <- data.frame(pivot_wider(lyta_tab2, 
                                    names_from = Day, 
                                    values_from = log10_avgconc))
  rownames(lytawid) <- lytawid[,1]
  lytawid2 <- lytawid[,-1]
  colnames(lytawid2) <- gsub( 'X', '',  colnames(lytawid2))
  lytawid3 <- lytawid2[,order(as.numeric(colnames(lytawid2)))]
  lytawid4 <- lytawid3[,1:29]
  
  #make a lagging list up to 21 days
  lags <- list()
  for(j in 1:21){
    lag1 <- lytawid4[,c(1,(1+j))]
    colnames(lag1) <- c('t0', 't1')
    lag1$child <- rownames(lytawid4)
    for(i in 2:28){
      if((i+j)<30){
        lag1.0 <-  lytawid4[,c(i,(i+j))]
        colnames(lag1.0) <- c('t0', 't1')
        lag1.0$child <- rownames(lytawid4)
        lag1 <- rbind(lag1, lag1.0)
      }
    }
    lags[[j]] <- lag1
  }
  
  
  #mixed effect model per childr
  all_lmer <- lapply(lags, function(i){
    lmm <- lmer(t1 ~ t0 + (1 | child), data=i)
    summary(lmm)[[10]][2,]
  })
  all_lmer2 <- data.frame(do.call(rbind, all_lmer))
  all_lmer2$lag <- 1:21
  all_lmer2$Assay <- select_assay
  all_lmer2 <- all_lmer2[1:10,]
  assign(paste0(select_assay, '_totalLMER'), all_lmer2)
  
  
  t1 <- lags[[1]]
  t3 <- lags[[3]]
  t9 <- lags[[9]]
  p1 <- ggplot(t1, aes(x=t0, y=t1)) + geom_point(alpha=0.8)  + ggtitle('Lag 1 day') + 
    theme(aspect.ratio=1)
  p2 <- ggplot(t3, aes(x=t0, y=t1)) + geom_point(alpha=0.8) + ylab('t3')  + ggtitle('Lag 3 days')+ 
    theme(aspect.ratio=1)
  p3 <- ggplot(t9, aes(x=t0, y=t1)) + geom_point(alpha=0.8) + ylab('t9') + ggtitle('Lag 9 days')+ 
    theme(aspect.ratio=1)
  filename <- paste0('Fig3K_', select_assay,'autocorrelation_examples.pdf')
  pdf(filename, width = 9, height = 3)
  print(plot_grid(p1,p2,p3, ncol=3))
  dev.off()
  
  
  
  ###############look for stability per age
  ages <- data.frame(mrg2[!duplicated(mrg2$StudyID),c(2,7)])
  ages_list <- split(ages, ages$AgeYears)
  
  
  #########do this fper age
  k=0
  age_lmer <- list()
  for(age in ages_list){
    k=k+1
    
    lytawid5 <- lytawid4[rownames(lytawid4) %in% age$StudyID,]
    
    lags <- list()
    for(j in 1:21){
      lag1 <- lytawid5[,c(1,(1+j))]
      colnames(lag1) <- c('t0', 't1')
      lag1$child <- rownames(lytawid5)
      for(i in 2:28){
        if((i+j)<30){
          lag1.0 <-  lytawid5[,c(i,(i+j))]
          colnames(lag1.0) <- c('t0', 't1')
          lag1.0$child <- rownames(lytawid5)
          lag1 <- rbind(lag1, lag1.0)
        }
      }
      lags[[j]] <- lag1
    }
    
    #mixed effect model per childr
    all_lmer <- lapply(lags, function(i){
      lmm <- lmer(t1 ~ t0 + (1 | child), data=i)
      summary(lmm)[[10]][2,]
    })
    all_lmer2 <- data.frame(do.call(rbind, all_lmer))
    all_lmer2$lag <- 1:21
    all_lmer2 <- all_lmer2[1:14,]
    all_lmer2$Age <- k
    age_lmer[[k]] <- all_lmer2
  }
  
  age_lmers <- do.call(rbind, age_lmer)
  
  age_lmers <- age_lmers[age_lmers$lag < 11,] 
  age_lmers$Assay <- select_assay
  colnames(age_lmers)[5] <- 'Pval'
  age_lmers$Padj <- p.adjust(age_lmers$Pval, method='BH')
  assign(paste0(select_assay, '_agedLMER'), age_lmers)
  
  
  mycols3 <- c("Estimate",
               "lag",
               "Age")
  mycols4 <- c("Padj",
               "lag",
               "Age")
  age_lmers2 <- age_lmers %>% dplyr::select(all_of(mycols3)) 
  age_lmers3 <- age_lmers2 %>%
    pivot_wider(names_from = lag, 
                values_from = c(Estimate))
  age_lmers4 <- data.frame(age_lmers3[,2:ncol(age_lmers3)])
  colnames(age_lmers4)<- gsub('X', '', colnames(age_lmers4))
  
  age_lmersp <- age_lmers %>% dplyr::select(all_of(mycols4)) 
  age_lmersp2 <- age_lmersp %>%
    pivot_wider(names_from = lag, 
                values_from = c(Padj))
  age_lmersp4 <- age_lmersp3 <- data.frame(age_lmersp2[,2:ncol(age_lmersp2)])
  age_lmersp4[age_lmersp3>0.05] <- ''
  age_lmersp4[age_lmersp3<0.05] <- '*'
  age_lmersp4[age_lmersp3<0.01] <- '**'
  age_lmersp4[age_lmersp3<0.001] <- '***'
  
  filename <- paste0('Fig3K_', select_assay,'Fig3K_autocorrelation_lmm_perage_HM.pdf')
  pdf(filename, width = 6, height =3)
  pl2 <- pheatmap(age_lmers4, cluster_rows = F, cluster_cols=F, 
                  display_numbers = age_lmersp4, main = select_assay)
  print(pl2)
  #dev.off()
}

