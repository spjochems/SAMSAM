#####
#This script creates UMAPS combining all the pathogen data and see if they associate
#with parameters
#look at effect of covariates on denstity

rm(list=ls())

setwd("Z:\\Projects\\SAMSAM\\Manuscripts\\Figures\\Supplemental figure 6_timeofday")

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
sink("sessionInfo.txt")
sessionInfo()
sink()


#load in database and select key bacteria
mrg <- read_xlsx('Z:\\Projects\\SAMSAM\\IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20250722.xlsx')

#get out numbers for tables
table(mrg$Assay)

mrg_ex <- mrg[!duplicated(mrg$StudyID),]
table(mrg_ex$AgeYears, mrg_ex$Gender)

mrg_ex2 <- mrg[mrg$Assay == '16S',]
table(table(mrg_ex2$StudyID))


mrg_ex3 <- mrg[!duplicated(mrg$SampleID_Day),]
ques <- mrg_ex3[mrg_ex3$QuestionereFill == 'YesQFill',]
data.frame(table(ques$SchoolDaycare, ques$StudyID)[2,])




#make associations with age
mrg2a <- mrg[grepl("^Spn\\(lytA\\)$", mrg$Assay) | 
              grepl("Haemophilus", mrg$Assay) |
               grepl("pyogenes", mrg$Assay) |
              grepl("Moraxella", mrg$Assay)|
              grepl("aureus", mrg$Assay), ]
mrg2 <- mrg2a[mrg2a$log10_avgconc>0,]




###########################################################
#
#Fig S6A
#
#######################################################################
#effect of period of collection
mrg3 <- mrg2[mrg2$PeriodCollection != 'NA', ]
mrg3$PeriodCollection <- factor(mrg3$PeriodCollection,
                                levels = c('Morning', 'Afternoon', 'Evening'))

stats <- data.frame(EstimatePM = NA,
                    PvalPM = NA,
                    EstimateEve = NA,
                    PvalEve = NA)

mrg3.split <- split(mrg3, f=mrg3$Assay)
for(i in 1:length(mrg3.split)){
  test <- mrg3.split[[i]]
  m <- lmer(log10_avgconc ~ PeriodCollection +(1|StudyID), data = test)
  resAft <- lmerTest:::get_coefmat(m)[2,c(1,5)]
  resEve <- lmerTest:::get_coefmat(m)[3,c(1,5)]
  stats <- rbind(stats, c(resAft, resEve))
}
stats_period <- stats[-1,]
rownames(stats_period) <- stats_period$Assay <- names(mrg3.split)
stats_period$Label = paste0('EstPM: ', round(stats_period$EstimatePM, 3), ', PvalPM: ', scales::scientific(stats_period$PvalPM, 2), '\n',
                            'EstEve: ', round(stats_period$EstimateEve, 3), ', PvalEve: ', scales::scientific(stats_period$PvalEve, 2))

mrg3.period <- merge(mrg3, stats_period, by= 'Assay')
mrg3.period$Label2 <- mrg3.period$Label
mrg3.period$Label2[duplicated(mrg3.period$Label)] <- NA

pl3<- ggplot(mrg3.period, aes(x=as.factor(PeriodCollection), y=log10_avgconc)) + 
  facet_wrap(.~Assay, scale='free', ncol=5) + 
  geom_violin(aes(fill=PeriodCollection)) +   geom_boxplot(width=0.2) + 
  ggtitle('Period of collection, all kids/samples') + 
  geom_text(aes(label = Label2), y=Inf, x=Inf, hjust = 1, vjust = 1, colour = 'red')
pl3
pdf_filename <- paste0("fig_S6A_density_periodcollection_3periods_no_swithcing_times.pdf")
ggsave(pdf_filename, width = 15, height = 3)





###########################################################
#
#Fig S6B
#
#######################################################################
mrg3$Period2 <- 'PM'
mrg3$Period2[mrg3$PeriodCollection == 'Morning'] <- 'AM'
stats <- data.frame(Estimate = NA, Pval = NA)
mrg3.split <- split(mrg3, f=mrg3$Assay)
for(i in 1:length(mrg3.split)){
  test <- mrg3.split[[i]]
  m <- lmer(log10_avgconc ~ Period2 +(1|StudyID), data = test)
  res <- lmerTest:::get_coefmat(m)[2,c(1,5)]
  stats <- rbind(stats, res)
}
stats_period2 <- stats[-1,]
rownames(stats_period2) <- stats_period2$Assay <- names(mrg3.split)
stats_period2$Label = paste0('Estimate: ', round(stats_period2$Estimate, 3), ', Pval: ', scales::scientific(stats_period2$Pval, 2))

mrg3.period2 <- merge(mrg3, stats_period2, by= 'Assay')
mrg3.period2$Label2 <- mrg3.period2$Label
mrg3.period2$Label2[duplicated(mrg3.period2$Label)] <- NA

pl4<-ggplot(mrg3.period2, aes(x=Period2, y=log10_avgconc)) + 
  facet_wrap(.~Assay, scale='free', ncol=5) + 
  geom_violin(aes(fill=Period2)) +  
  geom_boxplot(width=0.2) + 
  ggtitle('Period of collection AM/PM, all times') + 
  geom_text(aes(label = Label2), y=Inf, x=Inf, hjust = 1, vjust = 1, colour = 'red')
pl4
pdf_filename <- paste0("fig_S6B_density_periodcollection_no_swithcing_times.pdf")
ggsave(pdf_filename, width = 10, height = 3)



save.image('Figure_S6_switching times.Rdata')


