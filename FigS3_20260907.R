#####
#This script creates UMAPS combining all the pathogen data and see if they associate
#with parameters
#look at effect of covariates on denstity

rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\FigS3\\")

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

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#load in database and select key bacteria
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))
colnames(mrg)

#add the cohort
mrg$Cohort <- 'Pilot'
mrg$Cohort[grepl('K1', mrg$StudyID)] <-'Main'
table(mrg$Cohort)

mrg2 <- mrg[mrg$Antibiotics != 'Yes',] #remove the antibiotics samples as it may affect results

#create a function that fixes the display of p-values
disp_pval <- function(p){
  if(p<0.001){
    p <- scales::scientific(p, digits= 3)
  } else {
    p <- round(p, 3)
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
#                 FIGURE S3A                      #
#                                                 #
###################################################

#make associations with age
mrg2a <- mrg2[grepl("^Spn\\(lytA\\)$", mrg2$Assay) | 
               grepl("Haemophilus", mrg2$Assay) |
               grepl("pyogenes", mrg2$Assay) |
               grepl("Moraxella", mrg2$Assay)|
               grepl("aureus", mrg2$Assay), ]
mrg3 <- mrg2a[mrg2a$log10_avgconc>0,]

#check for age
stats <- data.frame(Estimate = NA, Pval = NA)
mrg3.split <- split(mrg3, f=mrg3$Assay)
for(i in 1:length(mrg3.split)){
  test <- mrg3.split[[i]]
  m <- lmer(log10_avgconc ~ Cohort + AgeYears +(1|StudyID), data = test)
  res <- lmerTest:::get_coefmat(m)[2,c(1,5)]
  stats <- rbind(stats, res)
}
stats_age <- stats[-1,]
rownames(stats_age) <- stats_age$Assay <- names(mrg3.split)
stats_age$Padj <- 5*stats_age$Pval
stats_age$Padj[stats_age$Padj >1] <- 1
stats_age$Label = paste0('Est: ', round(stats_age$Estimate, 2), ', P=', sapply(stats_age$Padj, disp_pval) )

mrg3.cohort <- merge(mrg3, stats_age, by= 'Assay')
mrg3.cohort$Label2 <- mrg3.cohort$Label
mrg3.cohort$Label2[duplicated(mrg3.cohort$Label)] <- NA
mrg3.cohort$Cohort <- factor(mrg3.cohort$Cohort, levels = c('Pilot', 'Main'))
pl3 <- ggplot(mrg3.cohort, aes(x = Cohort, y = log10_avgconc)) +  
  facet_wrap(. ~ Assay,  ncol=5) + 
  geom_violin(aes(group = Cohort, fill=Assay)) +   
  geom_boxplot(aes(group = Cohort), width = 0.2) +  
  ggtitle("Effect of Cohort, all kids/samples") + 
  geom_text(aes(label = Label2), y = Inf, x = Inf, hjust = 1.2, vjust = 1, colour = "red")
pl3

pdf_filename <- paste0("FigS3A_density_cohort.pdf")
ggsave(pdf_filename, width = 10, height = 3)


#nr samples:
fig_data      <- mrg3
group_cols    <- c("Assay", "Cohort")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "FigS3A_samplesize.txt")



###################################################
#                                                 #
#                 FIGURE S3B                      #
#                                                 #
################################################### 

mrg_virus <- mrg[mrg$AssayType == 'Viral',]
mrg2_infected <- mrg_virus[mrg_virus$InfectionLongerTrue == 1,] #take only infected samples to look at the assay and kids


infections_num <- infections2 <- infections <- table(mrg2_infected$Assay, mrg2_infected$StudyID) #make a table with number of infections
infections2[infections>0] <- 1 #make a matrix with kids per sample

#add the one child with no viral infections
infections3 <- cbind(infections2, rep(0,13)) 
colnames(infections3)[45] <- 'SAM-K1-25'

infections_num3 <- cbind(infections_num, rep(0,13)) 
colnames(infections_num3)[45] <- 'SAM-K1-25'

#get the data of the kids
demo <- mrg[!duplicated(mrg$StudyID),c('StudyID', 'Cohort', 'AgeYears', 'Gender')]
demo2 <- demo[match(colnames(infections_num3), demo$StudyID),] #order them in same way as heatmap

demo3 <- cbind(demo2, data.frame(Number = colSums(infections3)))

round(prop.table(table(demo3$Number, demo3$Cohort), 2)*100, 1)
table(demo3$Number, demo3$Cohort)
fisher.test(table(demo3$Number, demo3$Cohort))


###################################################
#                                                 #
#                 FIGURE S3C                      #
#                                                 #
###################################################

D1 <- mrg2[mrg2$Day == 'D01',]
D2 <- mrg2[mrg2$Day == 'D02',]
D1$PersonCollect <- 'Clinician'
D2$PersonCollect <- 'Parent'
mrg2_main <- mrg2[mrg2$Cohort == 'Main',]
D3 <- mrg2_main[mrg2_main$Day == 'D28',]
D3$PersonCollect <- 'Parent'
D4 <- mrg2_main[mrg2_main$Day == 'D29' | mrg2_main$Day == 'D30' |
                  mrg2_main$Day == 'D31' | mrg2_main$Day == 'D32',]
D4$PersonCollect <- 'Clinician'

D12 <- rbind(D1, D2, D3, D4)

D12a <- D12[grepl("^Spn\\(lytA\\)$", D12$Assay) | 
                grepl("Haemophilus", D12$Assay) |
                grepl("pyogenes", D12$Assay) |
                grepl("Moraxella", D12$Assay)|
                grepl("aureus", D12$Assay), ]
D12b <- D12a[D12a$log10_avgconc>0,]



#check for age
stats <- data.frame(Estimate = NA, Pval = NA)
D12b.split <- split(D12b, f=D12b$Assay)
for(i in 1:length(D12b.split)){
  test <- D12b.split[[i]]
  m <- lmer(log10_avgconc ~ PersonCollect + (1|StudyID), data = test)
  res <- lmerTest:::get_coefmat(m)[2,c(1,5)]
  stats <- rbind(stats, res)
}
stats_sampler <- stats[-1,]
rownames(stats_sampler) <- stats_sampler$Assay <- names(D12b.split)
stats_sampler$Padj <- 5*stats_sampler$Pval
stats_sampler$Padj[stats_sampler$Padj >1] <- 1
stats_sampler$Label = paste0('Est: ', round(stats_sampler$Estimate, 2), ', P=', sapply(stats_sampler$Padj, disp_pval))

D12b.cohort <- merge(D12b, stats_sampler, by= 'Assay')
D12b.cohort$Label2 <- D12b.cohort$Label
D12b.cohort$Label2[duplicated(D12b.cohort$Label)] <- NA

pl4 <- ggplot(D12b.cohort, aes(x = PersonCollect, y = log10_avgconc)) +  
  facet_wrap(. ~ Assay,  ncol=5) + 
  geom_violin(aes(group = PersonCollect, fill=Assay)) +   
  geom_boxplot(aes(group = PersonCollect), width = 0.2) +  
  ggtitle("Effect of PersonCollect, all kids/samples") + 
  geom_text(aes(label = Label2), y = Inf, x = Inf, hjust = 1, vjust = 1, colour = "red")
pl4

pdf_filename <- paste0("FigS3C_density_parent.pdf")
ggsave(pdf_filename, width = 10, height = 3)


#nr samples:
fig_data      <- D12b.cohort
group_cols    <- c("Assay", "PersonCollect")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "FigS3C_samplesize.txt")






###############################################
save.image('FigureS3_QC_data.Rdata')


