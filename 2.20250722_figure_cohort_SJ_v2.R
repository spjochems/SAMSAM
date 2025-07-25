#####
#This script creates UMAPS combining all the pathogen data and see if they associate
#with parameters
#look at effect of covariates on denstity

rm(list=ls())

setwd("Z:\\Projects\\SAMSAM\\Manuscripts\\Figures\\Figure1_cohort\\")

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
table(mrg_ex$AgeYears, mrg_ex$StudyID)


mrg_ex2 <- mrg[mrg$Assay == '16S',]
table(table(mrg_ex2$StudyID))

(table(mrg_ex2$QuestionereFill))
table(mrg_ex2$SchoolDaycare)

table(mrg_ex2$PeriodCollection)


###################################################
#
#FIGURE 1B
#
################################################
#check per child how many positives
mrg.split <- split(mrg, f=mrg$StudyID)
for(i in 1:length(mrg.split)){
  pres <- table(mrg.split[[i]]$Assay, mrg.split[[i]]$Positive>0)
  if(i == 1){ 
    presence0 <- pres[,2]
  }else{
    presence0 <- cbind(presence0, pres[,2])  
  }
}
colnames(presence0) <- names(mrg.split)
presence <- presence0[-1,]

#how many detected per child
presence_child <- data.frame(rowSums(presence>0) / ncol(presence))
presence_child$Assay <- rownames(presence_child)
colnames(presence_child)[1] <- 'PresenceAtLeastOnce'

presence_child$Assay <- factor(presence_child$Assay, 
                               levels = presence_child$Assay[order(presence_child$PresenceAtLeastOnce)]) 

pl1 <- ggplot(presence_child, 
              aes(y=Assay, x=PresenceAtLeastOnce*100)) + 
  geom_bar(stat='identity', 
           fill="steelblue", 
           colour="black")+
  xlab('%children positive at least once')

pl1

presence_melt <- melt(apply(presence, 2, function(u) u/max(u)), 
                      by=row.names)



presence_melt$Var1 <- factor(presence_melt$Var1, 
                             levels = presence_child$Assay[order(presence_child$PresenceAtLeastOnce)]) 

presence_melt2 <- presence_melt[presence_melt$value>0,]


pl2 <- ggplot(presence_melt2, aes(x=Var1, y=value*100)) + 
  geom_boxplot(fill="steelblue") + 
  coord_flip() + 
  ylab('% of positive samples per child')

pl2

plot_grid(pl1, pl2)
timestamp <- format(Sys.time(), "%Y%m%d")

pdf_filename <- paste0("fig1B_presandfreq_", timestamp, ".pdf")
ggsave(pdf_filename, width = 20, height = 10)



mrg_ex3 <- mrg[!duplicated(mrg$SampleID_Day),]
ques <- mrg_ex3[mrg_ex3$QuestionereFill == 'YesQFill',]
data.frame(table(ques$SchoolDaycare, ques$StudyID)[2,])



###################################################
#
#FIGURE 1C
#
################################################
#make associations with age
mrg2a <- mrg[grepl("^Spn\\(lytA\\)$", mrg$Assay) | 
              grepl("Haemophilus", mrg$Assay) |
               grepl("pyogenes", mrg$Assay) |
              grepl("Moraxella", mrg$Assay)|
              grepl("aureus", mrg$Assay), ]
mrg2 <- mrg2a[mrg2a$log10_avgconc>0,]

#check for age
stats <- data.frame(Estimate = NA, Pval = NA)
mrg2.split <- split(mrg2, f=mrg2$Assay)
for(i in 1:length(mrg2.split)){
  test <- mrg2.split[[i]]
  m <- lmer(log10_avgconc ~ AgeYears +(1|StudyID), data = test)
  res <- lmerTest:::get_coefmat(m)[2,c(1,5)]
  stats <- rbind(stats, res)
}
stats_age <- stats[-1,]
rownames(stats_age) <- stats_age$Assay <- names(mrg2.split)
stats_age$Label = paste0('Linear Est: ', round(stats_age$Estimate, 3), ', Pval: ', scales::scientific(stats_age$Pval, 2))

mrg2.age <- merge(mrg2, stats_age, by= 'Assay')
mrg2.age$Label2 <- mrg2.age$Label
mrg2.age$Label2[duplicated(mrg2.age$Label)] <- NA

pl3 <- ggplot(mrg2.age, aes(x = AgeYears, y = log10_avgconc)) +  
  facet_wrap(. ~ Assay,  ncol=5) + 
  geom_violin(aes(group = AgeYears, fill=Assay, alpha=AgeYears)) +   
  geom_boxplot(aes(group = AgeYears), width = 0.2) +  
  ggtitle("Effect of age, all kids/samples") + 
  geom_text(aes(label = Label2), y = Inf, x = Inf, hjust = 1, vjust = 1, colour = "red")
pl3
pdf_filename <- paste0("fig1C_density_age_", timestamp, ".pdf")
ggsave(pdf_filename, width = 10, height = 3)



###################################################
#
#FIGURE 1D
#
################################################
#effect of period of collection
mrg3 <- mrg2a[mrg2a$PeriodCollection != 'NA', ]
mrg3$PeriodCollection <- factor(mrg3$PeriodCollection,
                                levels = c('Morning', 'Afternoon', 'Evening'))
mrg3$Period2 <- 'PM'
mrg3$Period2[mrg3$PeriodCollection == 'Morning'] <- 'AM'

ex <- mrg3[mrg3$Assay == 'Haemophilus_influenzae',]
ex.split <- split(ex, f=ex$StudyID)
for(i in 1:length(ex.split)){
  ex2 <- ex.split[[i]]
  for(j in 2:nrow(ex2)){
    if(ex2$Period2[j] != ex2$Period2[j-1]){
      keep <- append(keep, ex2$SampleID_Day[j-1])
      keep <- append(keep, ex2$SampleID_Day[j])
      keep <- append(keep, ex2$SampleID_Day[j+1])
    }
  }
}
mrg3.selecteda <- mrg3[mrg3$SampleID_Day %in% unique(keep) ,]
mrg3.selected <- mrg3.selecteda[mrg3.selecteda$log10_avgconc>0,]



#########only selected AM vs PM
stats <- data.frame(Estimate = NA, Pval = NA)
mrg3.split <- split(mrg3.selected, f=mrg3.selected$Assay)
for(i in 1:length(mrg3.split)){
  test <- mrg3.split[[i]]
  m <- lmer(log10_avgconc ~ Period2 +(1|StudyID), data = test)
  res <- lmerTest:::get_coefmat(m)[2,c(1,5)]
  stats <- rbind(stats, res)
}
stats_period2 <- stats[-1,]
rownames(stats_period2) <- stats_period2$Assay <- names(mrg3.split)
stats_period2$Label = paste0('Estimate: ', round(stats_period2$Estimate, 3), ', Pval: ', scales::scientific(stats_period2$Pval, 2))

mrg3.period2 <- merge(mrg3.selected, stats_period2, by= 'Assay')
mrg3.period2$Label2 <- mrg3.period2$Label
mrg3.period2$Label2[duplicated(mrg3.period2$Label)] <- NA

pl6<-ggplot(mrg3.period2, aes(x=Period2, y=log10_avgconc)) + 
  facet_wrap(.~Assay, scale='free', ncol=5) + 
  geom_violin(aes(fill=Period2)) +  
  geom_boxplot(width=0.2) + 
  ggtitle('Period of collection AM/PM, switching times') + 
  geom_text(aes(label = Label2), y=Inf, x=Inf, hjust = 1, vjust = 1, colour = 'red')
pl6
pdf_filename <- paste0("fig1D_density_periodcollection_", timestamp, ".pdf")
ggsave(pdf_filename, width = 10, height = 3)




###################################################
#
#FIGURE 1E
#
################################################
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
  facet_wrap(.~StudyID) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "right"
  ) +  # Rotate x-axis labels
  labs(
    title = 'S. pneumoniae carriage',  # Title for each StudyID
    x = NULL, y = "log10(average_conc)"  # y-axis label with log scale
  )
ggsave("spn_atleast2.pdf", width = 12, height = 10)



#plots for the not spn
notspn2 <- mrg9[grepl("^Spn\\(lytA\\)$", mrg9$Assay) | !grepl("^Spn", mrg9$Assay), ]
notspn2 <- notspn2 %>% filter(Assay != "EAV")  %>% filter(Assay != "16S")
ggplot(notspn2, aes(x = Day2, y = log10_avgconc)) +  # Convert Day to numeric for regression
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
ggsave("bac_vir_atleast2.pdf", width = 12, height = 10)



save.image('Figure1_cohort.Rdata')


