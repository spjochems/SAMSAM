#####
#This script uses the mlvar package to calculate correlations between bacteria
#https://arxiv.org/pdf/1609.04156
#https://cran.r-project.org/web/packages/mlVAR/mlVAR.pdf
#Spn (lytA) / Haemophilus / Moraxella / Staph 
#it uses the RAW (log-transformed) data or with data normalized to 16S 

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

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\FigS10\\")

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
library(RColorBrewer)
library(patchwork)

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#load in database
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))

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

mrg$DateCollection_EditDate #check what the date looks like
mrg$DaysSinceJan01 <- as.integer(format(as.Date(mrg$DateCollection_EditDate, format = "%Y_%m_%d"), "%j")) - 1 #set it to an integer
mrg$DaysCosinus <-  cos(2*pi*mrg$DaysSinceJan01/365) #change it to a cosinus


#filter out antibiotics samples
mrg1 <- mrg <- mrg %>% 
  filter(!grepl("Yes", Antibiotics))   #remove anitbiotics samples

mrg2b <- mrg1 %>% #remove pathogens/child combinations where there is no infection
  group_by(StudyID, Assay) %>%
  filter(
    max(InfectionLongerTrue) > 0
  ) %>%
  ungroup()

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
#                 FIGURE S10A                      #
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

targets0 <- c('Staphylococcus_aureus', 'Streptococcus_pyogenes')

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

#colour scheme:
infection_colours <- c(
  "0.None"            = "#F8766D",
  "bocavirus"         = "#C2A500",
  "coronavirus"       = "#3DAF3D",
  "entero_rhino_virus" = "#2CB1B8",
  "influenza"         = "#7196C9",
  "z.multiple_viruses" = "#BE73A8"
)

#make the plot
FigS10A <- ggplot(mrg_all7, aes(x=InfectionTypeBroad, y=log10_avgconc)) + 
  geom_violin(aes(fill = InfectionTypeBroad)) + geom_boxplot(width = 0.2)+
  xlab('') +
  facet_grid(.~Assay, scale = 'free_x', space = 'free', ) + 
  geom_text(data= mrg_all7_summary, aes(InfectionTypeBroad, Inf, label=n), vjust = 1) +
  ggtitle(paste0('Infection virus effect on bacterial density')) + 
  geom_text(data= mrg_all7_inf, aes(InfectionTypeBroad, Inf, label=Pdisp_adj2), vjust = 2.5)  +
  scale_fill_manual(values = infection_colours, drop = FALSE)
FigS10A

#save it
timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("FigS10A_virus_bacterial_density_",
                       timestamp, ".pdf")
pdf(pdf_filename, width = 4,height = 3)
FigS10A
dev.off()


###################################################
#                                                 #
#                 FIGURE S10B                     #
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
  lmm <- lmer(log10_avgconc ~ as.factor(InfectionTypeBroad) + AgeYears + DaysCosinus + (1|StudyID) , data = mrg_all5)
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
FigS10B <- ggplot(mrg_all7, aes(x=InfectionTypeBroad, y=log10_avgconc)) + 
  geom_violin(aes(fill = InfectionTypeBroad)) + geom_boxplot(width = 0.2)+
  xlab('') +
  facet_grid(.~Assay, scale = 'free_x', space = 'free', ) + 
  geom_text(data= mrg_all7_summary, aes(InfectionTypeBroad, Inf, label=n), vjust = 1) +
  ggtitle(paste0('Infection virus effect on bacterial density')) + 
  geom_text(data= mrg_all7_inf, aes(InfectionTypeBroad, Inf, label=Pdisp_adj2), vjust = 2.5)
FigS10B

#save it
timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("FigS10B_virus_bacterial_density_",
                       timestamp, ".pdf")
pdf(pdf_filename, width = 8,height = 3)
FigS10B
dev.off()




###################################################
#                                                 #
#                 FIGURE S10C                     #
#                                                 #
###################################################

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
mrg_no_viral <- mrg_all7 %>%
  group_by(StudyID, Assay) %>%
  arrange(Day2, .by_group = TRUE) %>%
  mutate(
    start_ok = InfectionTypeBroad[Day2 == 1][1] == "0.None",
    keep = if_else(start_ok, cumall(InfectionTypeBroad == "0.None"), FALSE)
  ) %>%
  filter(keep) %>%
  filter(Assay != "Streptococcus_pyogenes") %>%
  filter(Day2 <= 14) %>%
  ungroup()

#check:
unique(mrg_no_viral$InfectionTypeBroad)
unique(mrg_no_viral$StudyID)
unique(mrg_no_viral$Assay)
length(unique(mrg_no_viral$SampleID_Day))
check <- mrg_no_viral %>%
  dplyr::select(StudyID, Day2, Assay) %>%
  filter(Assay == "Spn(lytA)")
mrg_all7 %>%
  dplyr::filter(StudyID == "SAM001", Day2 == 12) %>%
  dplyr::pull(InfectionTypeBroad)

#reorder:
mrg_no_viral <- mrg_no_viral %>%
  filter(Assay %in% c("Haemophilus_influenzae",
                      "Spn(lytA)",
                      "Staphylococcus_aureus",
                      "Moraxella_catarrhalis")) %>%
  mutate(Assay = factor(Assay,
                        levels = c("Haemophilus_influenzae",
                                   "Spn(lytA)",
                                   "Staphylococcus_aureus",
                                   "Moraxella_catarrhalis")))

#plot:
FigS10C1 <- ggplot(mrg_no_viral, aes(x = Day2, y = log10_avgconc, colour = Assay)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_line(aes(group = StudyID), alpha = 0.5) +
  geom_smooth(aes(fill = Assay), method = "loess", se = TRUE, size = 0.7, span = 0.6) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0),
    legend.position = "none"
  ) +
  labs(
    title = "Donors with consecutive days without viral infection",
    x = "Consecutive days without viral infection",
    y = "Log10 (genome copies/mL +1)"
  ) +
  facet_wrap(. ~ Assay, ncol = 4)
FigS10C1

#reorder:
mrg_all7 <- mrg_all7 %>%
  filter(Assay %in% c("Haemophilus_influenzae",
                      "Spn(lytA)",
                      "Staphylococcus_aureus",
                      "Moraxella_catarrhalis")) %>%
  mutate(Assay = factor(Assay,
                        levels = c("Haemophilus_influenzae",
                                   "Spn(lytA)",
                                   "Staphylococcus_aureus",
                                   "Moraxella_catarrhalis")))

#plot all samples:
FigS10C2 <- ggplot(mrg_all7[mrg_all7$Day2 <= 14 & mrg_all7$Assay != "Streptococcus_pyogenes",], aes(x = Day2, y = log10_avgconc, colour = Assay)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_line(aes(group = StudyID), alpha = 0.5) +
  geom_smooth(aes(fill = Assay), method = "loess", se = TRUE, size = 0.7, span = 0.6) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0),
    legend.position = "none"
  ) +
  labs(
    title = "Bacterial densities during first 14 dats, all participants",
    x = "Days",
    y = "Log10 (genome copies/mL +1)"
  ) +
  facet_wrap(. ~ Assay, ncol = 4)
FigS10C2


combined <- FigS10C2 / FigS10C1   # '/' stacks vertically; '|' would place side by side
combined

#save:
filename <- 'FigS10C_BacterialStability.pdf'
pdf(filename, width = 8, height = 5.5)
combined
dev.off()


#nr samples:
fig_data      <- mrg_no_viral
group_cols    <- c("Assay", "Day2")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "FigS10C1_samplesize.txt")


fig_data      <- mrg_all7
group_cols    <- c("Assay", "Day2")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "FigS10C2_samplesize.txt")





###################################################
#                                                 #
#                 FIGURE S10D                     #
#                                                 #
###################################################

viruses2 = c('HRV')

mrg2c <- mrg2b[grepl('Spn',mrg2b$Assay) & 
                 !grepl('lytA',mrg2b$Assay) & 
                 !grepl('piaB',mrg2b$Assay) & 
                 !grepl('Spn6C/D',mrg2b$Assay),]

targets1 <- unique(mrg2c$Assay)

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
      
      #if there is no base then put FC_dens at NA
      if(base == 0 | is.na(base)){each2$FC_dens <- NA}
      
      mrg7.split[[i]] <- each2 #save back in the list
    }
    mrg8 <- do.call(rbind, mrg7.split)
    
    mrg9 <- mrg8[mrg8$DayToVirus > -4 & mrg8$DayToVirus < 15, ] #filter out smaples between x days before or after viral onset
    mrg9$DayToVirus2 <- mrg9$DayToVirus #make a new column and change the 5 days prior to onset to the baseline
    mrg9$DayToVirus2[mrg9$DayToVirus2 %in% c(-3:-1)] <- '-1_baseline'
    
    
    #save the dataframe
    assign(paste0(target,'_', virus,'_acqdataset'), mrg9)
    
    
    if(r2==1){
      mrg10 <- mrg9
    } else{
      mrg10 <- rbind(mrg10, mrg9)
    }
  }
  
}

mrg11 <- mrg10[!is.na(mrg10$DayToVirus),]
mrg12 <- mrg11[mrg11$log10_avgconc != 0,]

colnames(mrg11)

mrg12$Assay <- factor(
  mrg12$Assay,
  levels = c(
    "Spn3",
    "Spn6A/B/C/D",
    "Spn7C",
    "Spn9A/L/N/V",
    "Spn10A",
    "Spn11A/B/C/D/F",
    "Spn12A/12B/12F/44/46",
    "Spn15A/B/C/F",
    "Spn16F",
    "Spn19A",
    "Spn21",
    "Spn22A/F",
    "Spn23A/B/F",
    "Spn25A/25F/38",
    "Spn33A/F/37",
    "Spn34/37/17A",
    "Spn35B"
  )
)

ggplot(mrg12, aes(x = DayToVirus, y = FC_dens, colour = Assay)) +  # make a plot.
  geom_vline(xintercept=0, linetype='dashed')+
  geom_hline(yintercept=0, linetype='dashed')+
  geom_line(aes(group = StudyPath), alpha=0.5)+
  geom_smooth(color = 'black', method = "loess", se = T, size = 0.7, span = .6) +    # Add separate smoothing line for each Assay
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "right"
  ) +  # Rotate x-axis labels
  labs(
    title = 'HRV time effect',  # Title for each StudyID
    x = 'Days since viral onset', y = "FC (log10)"  # y-axis label with log scale
  )  + xlim(c(-5,14)) 

filename <- paste0('FigS10D_Spn_serotype_timeefeect_preexistingonlyFC.pdf')
ggsave(filename, width = 6, height = 4)



###################################################
#                                                 #
#                 FIGURE S10E                      #
#                                                 #
###################################################

mrg12$New <- 'No'
mrg12$New[is.na(mrg12$FC_dens)] <- 'Yes'

ggplot(mrg12, aes(x = DayToVirus, y = log10_avgconc, colour = New, fill= New)) +  # make a plot.
  geom_vline(xintercept=0, linetype='dashed')+
  geom_hline(yintercept=0, linetype='dashed')+
  geom_line(aes(group = StudyPath), alpha=0.5)+
  geom_smooth(aes(color = New), method = "loess", se = T, size = 0.7, span = .6) +    # Add separate smoothing line for each Assay
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "right"
  ) +  # Rotate x-axis labels
  labs(
    title = 'HRV time effect',  # Title for each StudyID
    x = 'Days since viral onset', y = "log10 genome copies"  # y-axis label with log scale
  )  + xlim(c(-3,14)) 

filename <- paste0('Figs10E_SpnNewvsOld.pdf')
ggsave(filename, width = 4, height = 3)


#nr samples:
fig_data      <- mrg12
group_cols    <- c("New", "DayToVirus")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig10E_samplesize.txt")


###################################################
#                                                 #
#                 FIGURE S10F,G                    #
#                                                 #
###################################################

mycols <- c("StudyID",
            "Day",
            "Assay",
            "log10_avgconc")

mrg2 <- mrg1 


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

plot(VAR.splines, posCol = 'Red', negCol ='Blue', title = 'Temporal - loess', label.cex = 1, label.scale = F)
FigS10F <- plot(VAR.splines, type = 'contemporaneous', posCol = 'Red', negCol ='Blue', title = 'Contemporaneous - loess', 
               label.cex = 1, label.scale = F, edge.labels = TRUE)
FigS10F
FigS10G <- plot(VAR.splines, type = 'between', posCol = 'Red', negCol ='Blue', title = 'Between Subject - loess', 
               label.cex = 1, label.scale = F, edge.labels = TRUE)
FigS10G

pdf('FigS10FG_VAR_plots_loess_fits.pdf', width = 8, height = 3)
par(mfrow=c(1,2))
FigS10F
FigS10G
dev.off()


###################################################
#                                                 #
#      FIGURE S8H,I - 16S-normalized VAR          #
#                                                 #
###################################################

mrg2.norm <- mrg1[mrg1$Assay == 'Spn(lytA)' | 
                    grepl("Haemophilus_influenzae", mrg1$Assay)| 
                    grepl("aureus", mrg1$Assay)|
                    grepl("16S", mrg1$Assay)|
                    grepl("Moraxella_catarrhalis", mrg1$Assay),]

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

plot(VAR.norm, posCol = 'Red', negCol ='Blue', title = 'Temporal - 16S normalized', 
     label.cex = 1, label.scale = F, edge.labels = TRUE)
FigS10H <- plot(VAR.norm, type = 'contemporaneous', posCol = 'Red', negCol ='Blue', 
     title = 'Contemporaneous - 16S normalized', label.cex = 1, label.scale = F, edge.labels = TRUE)
FigS10H
FigS10I <- plot(VAR.norm, type = 'between', posCol = 'Red', negCol ='Blue', 
               title = 'Between Subject - 16S normalized', label.cex = 1, label.scale = F,
               edge.labels = TRUE)
FigS10I

pdf('FigS10HI_VAR_plots_bacteria_16Snormalized_lag1.pdf', width = 8, height = 3)
par(mfrow=c(1,2))
FigS10H
FigS10I
dev.off()



#################################################################################################
#                                                                                               #
#      FIGURE S8J,K - set the counter of max number of times in a row we see a pathogen         #
#                                                                                               #
#################################################################################################

virus_pos <- unique(mrg1[mrg1$AssayType == 'Viral' & mrg1$PositiveYesNo == 1,]$SampleID_Day) #remove any sample with a virus

mrg.a <- mrg1[!mrg1$SampleID_Day %in% virus_pos,] #remove all samples that are virus positive


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

plot(VARv, posCol = 'Red', negCol ='Blue', title = 'Temporal - raw', label.cex = 1, label.scale = F)
FigS10J <- plot(VARv, type = 'contemporaneous', posCol = 'Red', negCol ='Blue', 
               title = 'Contemporaneous - raw', label.cex = 1, label.scale = F, edge.labels = TRUE)
FigS10J
FigS10K <- plot(VARv, type = 'between', posCol = 'Red', negCol ='Blue', title = 'Between Subject - raw', 
               label.cex = 1, label.scale = F, edge.labels = TRUE)
FigS10K

pdf('FigS10JK_VAR_plots_bacteria_rawdata_lag1_wovirus.pdf', width = 8, height = 3)
par(mfrow=c(1,2))
FigS10J
FigS10K
dev.off()


###################################################
#                                                 #
#                 FIGURE S10L                     #
#                                                 #
###################################################

mrg0 <- mrg

mrg2.norm <- mrg0[mrg0$Assay == 'Spn(lytA)' | 
                    grepl("Haemophilus_influenzae", mrg0$Assay)| 
                    grepl("aureus", mrg0$Assay)|
                    grepl("16S", mrg0$Assay)|
                    grepl("Moraxella_catarrhalis", mrg0$Assay),]

mycols <- c("StudyID",
            "Day",
            "Assay",
            "log10_avgconc")

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
  filename <- paste0('Fig3XX_', select_assay,'autocorrelation_examples_16Normalised.pdf')
  pdf(filename, width = 9, height = 3)
  print(plot_grid(p1,p2,p3, ncol=3))
  dev.off()
  
  
  
  ###############look for stability per age
  ages <- data.frame(mrg1[!duplicated(mrg1$StudyID),c(2,7)])
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
  
  filename <- paste0('Fig3XX_', select_assay,'_autocorrelation_lmm_perage_HM_16normalised.pdf')
  pdf(filename, width = 6, height =3)
  pl2 <- pheatmap(age_lmers4, cluster_rows = F, cluster_cols=F, 
                  display_numbers = age_lmersp4, main = select_assay)
  print(pl2)
  dev.off()
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

mat4 <- as.matrix(all_LMER4)
storage.mode(mat4) <- "numeric"
star_mat <- as.matrix(all_LMERp4)
storage.mode(star_mat) <- "character"

p1 <- pheatmap(mat4, cluster_rows = T, cluster_cols=F, 
               display_numbers = star_mat, main = "All assays",show_rownames = T, 
               fontsize_number = 14 )
p1

filename <- 'FigS10L_all_assays_lmm_perassay_HM_16SNormalised.pdf'
pdf(filename, width = 6, height =3)
p1
dev.off()



###################################################
#                                                 #
#                 FIGURE S10M                      #
#                                                 #
###################################################

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

mrg3$Assay[mrg3$Assay == 'Spn(lytA)'] <- 'Spn'

#autocorrelation:
autocorrlag_plots <- list()
autocorr_plots <- list()

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
    theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5, size = 14),
          axis.title = element_text(size = 13), axis.text = element_text(size = 12))
  p2 <- ggplot(t3, aes(x=t0, y=t1)) + geom_point(alpha=0.8) + ylab('t3')  + ggtitle('Lag 3 days')+ 
    theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5, size = 14),
          axis.title = element_text(size = 13), axis.text = element_text(size = 12))
  p3 <- ggplot(t9, aes(x=t0, y=t1)) + geom_point(alpha=0.8) + ylab('t9') + ggtitle('Lag 9 days')+ 
    theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5, size = 14),
          axis.title = element_text(size = 13), axis.text = element_text(size = 12))
  
  autocorrlag_plots[[select_assay]] <- list(
    lag1 = p1,
    lag3 = p2,
    lag9 = p3
  )
  
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
  
  pl2 <- pheatmap(age_lmers4, cluster_rows = F, cluster_cols=F, 
                  fontsize = 12, fontsize_number = 18,
                  display_numbers = age_lmersp4, main = select_assay)
  
  autocorr_plots[[select_assay]] <- pl2
  dev.off()
}

#plot:
microbe <- c("Moraxella_catarrhalis")
FigS10M <- autocorr_plots[[microbe]]
FigS10M
dev.off()

#save:
pdf('FigS10M_autocorrelationsLag.pdf', width = 4, height = 3)
FigS10L
dev.off()


#nr samples:
fig_data      <- mrg %>% filter(Assay == "Moraxella_catarrhalis") %>% 
  filter(log10_avgconc > 0) %>% group_by(StudyID) %>% 
  slice_sample(n = 1) %>% ungroup()
group_cols    <- c("AgeYears")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig10M_samplesize.txt")







##############################################
save.image('FigureS10_bacteria.Rdata')

