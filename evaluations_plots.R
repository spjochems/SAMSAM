rm(list = ls())

library(ggplot2)
library(reshape2)
library(cowplot)
library(ggbeeswarm)

setwd("R:/Para-CIH/CIH Group member folders/Simon/Results/SamSam/evaluation forms")
data <- read.delim('2024_01_31_evaluationforms_all.txt')

data_complete <- data[data$Dropout == 'No',]
data_complete$Study <- data_complete$Study.group
data_complete$Study[data_complete$Study.group == 'Pilot'] <- ' Pilot (n=10)'
data_complete$Study[data_complete$Study.group == 'Main'] <- 'Main (n=35)'

ggplot(data_complete, aes(x=AgeYears, fill = sex)) + 
    geom_bar() + 
    theme_bw() +
    xlab('Age (years)') +
    ylab('Number of children')
ggsave2('Age_distribution.pdf', width = 5, height = 5, units = 'cm')


ggplot(data_complete, aes(x=AgeYears, fill = sex)) + 
  geom_bar() + 
  theme_bw() +
  xlab('Age (years)') +
  ylab('Number of children') + 
  facet_grid(.~Study)

ggsave2('Age_distribution_percohort.pdf', width = 10, height = 5, units = 'cm')




data_questionnaire <- data[data$FilledQuestionnaire == 'Yes',]
data2 <- melt(data_questionnaire, id =  colnames(data_questionnaire)[c(1:10,21)])
data2$Study <- data2$Study.group
data2$Study[data2$Study.group == 'Pilot'] <- ' Pilot (n=10)'
data2$Study[data2$Study.group == 'Main'] <- 'Main (n=29, 6 did not fill Q)'
ggplot(data2, aes(x=value, y= variable)) + 
    geom_jitter(aes(colour = as.factor(AgeYears), shape = sex), width = 0.1, height = 0.1) + 
    stat_summary(fun = 'median', geom = 'point', pch = 18, size = 5)+
    theme_bw() +
    xlab('Agree') +
    ylab('') + 
  facet_grid(.~Study)
    
ggsave2('Questionnaire_results_percohort.pdf', width = 20, height = 10, units = 'cm')



data3 <- data2[data2$variable == 'X10_Collection_NOTtaxing' | 
                 data2$variable == 'X9_SelfSamplingVsNurse_Prefer' |
                 data2$variable == 'X7_NasalCollect_Well' |
                 data2$variable == 'X5_Instr_Clear',]
ggplot(data3, aes(x=value, y= variable)) + 
  geom_jitter(aes(colour = as.factor(AgeYears), shape = sex), width = 0.1, height = 0.1) + 
  stat_summary(fun = 'median', geom = 'point', pch = 18, size = 5)+
  theme_bw() +
  xlab('Agree') +
  ylab('') + 
  facet_grid(.~Study)
ggsave2('Questionnaire_results_tolerated_percohort.pdf', width = 20, height = 6, units = 'cm')

