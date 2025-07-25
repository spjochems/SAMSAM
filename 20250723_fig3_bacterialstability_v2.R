#####
#This script uses the mlvar package to calculate correlations between bacteria
#https://arxiv.org/pdf/1609.04156
#https://cran.r-project.org/web/packages/mlVAR/mlVAR.pdf
#Spn (lytA) / Haemophilus / Moraxella / Staph 
#it uses the RAW (log-transformed) data or with data normalized to 16S 


rm(list=ls())

setwd("Z:\\Projects\\SAMSAM\\Manuscripts\\Figures\\Figure3_bacterialstability")


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
library(survival)
library(ggsurvfit)
library(ppcor)


#save the session info
sink("sessionInfo.txt")
sessionInfo()
sink()


#load in database
mrg <- read_xlsx('Z:\\Projects\\SAMSAM\\IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20250722.xlsx')
mrg <- mrg[mrg$Antibiotics != 'Yes',] #remove the antibiotics samples as it may affect results



############################
#
#Fig 3A
#
##################################################################
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
pdf('Fig3A_partialcorrelation_raw.pdf', width = 5, height = 4)
pheatmap(partialcor$estimate, display_numbers = rounded)
dev.off()


###############################################
#
#Fig3B/C
#
###############################################
VAR <- mlVAR(mrg4, vars = c('Spn', 'Haemophilus_influenzae',  'Staphylococcus_aureus', 'Moraxella_catarrhalis'), 
      idvar = 'StudyID', lags = 1, contemporaneous = 'correlated',
      estimator = 'lmer',)

print(VAR)
summary(VAR)

pdf('Fig3BC_VAR_plots_bacteria_rawdata_lag1.pdf', width = 10, height = 3)
par(mfrow=c(1,3))
plot(VAR, posCol = 'Red', negCol ='Blue', title = 'Temporal - raw', label.cex = 1, label.scale = F)
plot(VAR, type = 'contemporaneous', posCol = 'Red', negCol ='Blue', title = 'Contemporaneous - raw', label.cex = 1, label.scale = F)
plot(VAR, type = 'between', posCol = 'Red', negCol ='Blue', title = 'Between Subject - raw', label.cex = 1, label.scale = F)
dev.off()



##########################################
#
#Fig. 3E-F
#
##############################################
mrg3$Assay[mrg3$Assay == 'Spn(lytA)'] <- 'Spn'

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
filename <- paste0('Fig3E_', select_assay,'autocorrelation_examples.pdf')
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

filename <- paste0('Fig3F_', select_assay,'_autocorrelation_lmm_perage_HM.pdf')
pdf(filename, width = 6, height =3)
pl2 <- pheatmap(age_lmers4, cluster_rows = F, cluster_cols=F, 
         display_numbers = age_lmersp4, main = select_assay)
print(pl2)
dev.off()
}


###########################
#
#Fig 3D
#
#############################
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

filename <- 'Fig3D_all_assays_lmm_perassay_HM.pdf'
pdf(filename, width = 6, height =3)
pheatmap(all_LMER4, cluster_rows = T, cluster_cols=F, 
                display_numbers = all_LMERp4, main = "All assays",show_rownames = T)
dev.off()




##############################################################
#
#Fig 3G
#
#####################################################################
mrg_Spn <- mrg[grepl('Spn', mrg$Assay),]
lytapiaBPos <- mrg_Spn[(mrg_Spn$Assay == "Spn(lytA)" & mrg_Spn$InfectionLongerTrue == 1) | 
                     (mrg_Spn$Assay == "Spn(piaB)" & mrg_Spn$InfectionLongerTrue == 1),1]
mrg_Spn2 <- mrg_Spn[mrg_Spn$SampleID_Day %in% lytapiaBPos$SampleID_Day & 
                        mrg_Spn$Assay != "Spn(lytA)" &
                         mrg_Spn$Assay != "Spn(piaB)" &
                        mrg_Spn$Assay != "Spn6C/D",]

table(mrg_Spn2$Assay)
table(mrg_Spn2$Assay, mrg_Spn2$InfectionLongerTrue)
Total_serotype <- data.frame(table(table(mrg_Spn2$SampleID_Day, mrg_Spn2$InfectionLongerTrue)[,2]))
Total_serotype$Perc <- round(Total_serotype$Freq / sum(Total_serotype$Freq) * 100, 1)
pl3 <- ggplot(Total_serotype, 
              aes(y=Var1, x=Freq)) + 
  geom_bar(stat='identity', 
           fill="steelblue", 
           colour="black") + 
  scale_x_continuous(expand = c(0,0)) +
  geom_text(aes(label = Freq), hjust = 1, colour = 'white')+ 
  ylab('Number of co-detected Spn serotypes')
pl3
timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig3G_numberofcoserotypes_", timestamp, ".pdf")
ggsave(pdf_filename, width = 3, height = 2)

sum(Total_serotype$Freq[3:7])

#####################################3
#
#Fig. 3H 
#
########################################
lengthInfectionsLonger <- data.frame()

mrg8.split <- split(mrg, mrg$StudyPath)
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


lengthSpnLonger <- lengthInfections2Longer[grepl('Spn',  lengthInfections2Longer$Assay) & 
                                             !grepl('piaB',  lengthInfections2Longer$Assay)& 
                                             !grepl('Spn6C/D',  lengthInfections2Longer$Assay)   ,]
lengthSpnLonger$Assay2 <- lengthSpnLonger$Assay
lengthSpnLonger$Assay2[lengthSpnLonger$Assay != 'Spn(lytA)' & lengthSpnLonger$Assay !="Spn(piaB)"] <- 'Serotype'


pl7a <- survfit2(Surv(Duration, status) ~ Assay2, data = lengthSpnLonger) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability Spn"
  ) + 
  add_confidence_interval() +
  xlim(c(0,27)) + ggtitle('Spn survival')
pl7a


survdiff(Surv(Duration, status) ~ Assay2, data = lengthSpnLonger)
summary(survfit(Surv(Duration, status) ~ Assay2, data = lengthSpnLonger), times = 27)
summary(survfit(Surv(Duration, status) ~ Assay2, data = lengthSpnLonger))

filename <- 'Fig3H_survival_curves_Spn.pdf'
pdf(filename, width = 2.5, height =3)
print(pl7a)
dev.off()



######################################################
#
#Fig 3I
#
###########################################################
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




m <- glmer(as.factor(ClearNextDay) ~ AgeYears + SpnNumber2 + (1|Assay) + (1|StudyID), data = mrg_clear_spn, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

summary(m)[[10]][,1:2]
OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$SpnNumber <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - OR_adjusted$Std..Error)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + OR_adjusted$Std..Error)
ggplot(OR_adjusted, aes(x=OR, y=SpnNumber)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') 
ggsave('Fig3I_AdjustedORClearance.pdf', width = 2, height = 1.5)




##############################################
#
#Fig 3J
#
#######################################################
#remove any assays that are never found
mrg5 <- mrg %>%
  group_by(StudyID, Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLonger) == 1
  ) %>%
  ungroup() #make a new dataframe and only select assays that are present

#select only serotypes
mrg6 <- mrg5[grepl('Spn', mrg5$Assay) &
               mrg5$Assay != "Spn(lytA)" &
               mrg5$Assay != "Spn(piaB)" &
               mrg5$Assay != "Spn6C/D",  ]

mrg6$StudyID_Assay <- paste0(mrg6$StudyID, '_', mrg6$Assay)

#make a time since acquisition point
mrg6.split <- split(mrg6, mrg6$StudyID_Assay)
for(i in 1:length(mrg6.split)){
  each <- data.frame(mrg6.split[[i]])
  each$DaysSinceAcquisition <- NA
  if(sum(each$AcquisitionLonger == 1)){ #only do if there is one acquisition day otherwise it is hard to link to which one
    day <- each[which(each$AcquisitionLonger == 1),]$Day2[1]
    each$DaysSinceAcquisition <- each$Day2 - day
  }
  mrg6.split[[i]] <- each
}
mrg7 <- do.call(rbind, mrg6.split )
mrg7$DaysSinceAcquisition2 <- mrg7$DaysSinceAcquisition
mrg7$DaysSinceAcquisition2[mrg7$DaysSinceAcquisition2 < 0] <- NA

#divide the values by lyta and indicate whether it is dominant or not
mrg7$LytANorm <- NA
mrg7$Dominant <- 0
for(i in 1:nrow(mrg7)){
  p <- mrg7[i,]$SampleID_Day
  lyt <- mrg5[mrg5$SampleID_Day == p & mrg5$Assay == 'Spn(lytA)',]$Average_conc
  mrg7[i,]$LytANorm <- mrg7[i,]$Average_conc / lyt
}
mrg7$Dominant[mrg7$LytANorm>0.1] <- 1 #indicate if the serotype is >10% of the lytA then it is dominant
mrg8 <- mrg7[is.finite(mrg7$LytANorm) & !is.na(mrg7$LytANorm),] #remove any ones that are lytA negative

table(mrg8$Dominant, mrg8$Assay)


mrg4a <- mrg8[mrg8$InfectionLongerTrue == 1,]
keep <- names(which(table(mrg4a$SampleID_Day, mrg4a$InfectionLongerTrue)[,1] == 2))
mrg4.doublecarriage <- mrg4a[mrg4a$SampleID_Day %in% keep,]

mrg4.doublecarriage

mycols5 <- c("SampleID_Day", 'StudyID', 'Day', 'LytANorm',
             "log10_avgconc", "Assay")
mrg10 <- mrg4.doublecarriage %>% dplyr::select(all_of(mycols5)) 

mrg10a <- mrg10[!duplicated(mrg10$SampleID_Day),]
mrg10b <- mrg10[duplicated(mrg10$SampleID_Day),]
mrg10c <- merge(mrg10a, mrg10b[,-c(2,3)], by= 'SampleID_Day', suffixes = c('_1', '_2'))
mrg10d <- mrg10c[!is.na(mrg10c$LytANorm_1) & is.finite(mrg10c$LytANorm_1),]
cor.test(mrg10c$LytANorm_1, mrg10c$LytANorm_2)

mrg10d$Dominant <- 'none'
mrg10d$Dominant[mrg10d$LytANorm_1 >= 0.1 & mrg10d$LytANorm_2 < 0.1] <- mrg10d[mrg10d$LytANorm_1 >= 0.1 & mrg10d$LytANorm_2 < 0.1,]$Assay_1
mrg10d$Dominant[mrg10d$LytANorm_2 >= 0.1 & mrg10d$LytANorm_1 < 0.1] <- mrg10d[mrg10d$LytANorm_2 >= 0.1 & mrg10d$LytANorm_1 < 0.1,]$Assay_2
mrg10d$Dominant[mrg10d$LytANorm_2 >= 0.1 & mrg10d$LytANorm_1 >= 0.1] <- 'both'
table(mrg10d$Dominant, mrg10d$StudyID)

mrg10d$Dominant2 <- mrg10d$Dominant 
mrg10d$Dominant2[mrg10d$LytANorm_1 >= 0.1 & mrg10d$LytANorm_2 < 0.1] <- 'Serotype_1'
mrg10d$Dominant2[mrg10d$LytANorm_2 >= 0.1 & mrg10d$LytANorm_1 < 0.1] <- 'Serotype_2'

ggplot(mrg10d, aes(x=log10(LytANorm_1+.0000001), y = log10(LytANorm_2+.0000001))) + geom_point(aes(colour = Dominant2)) + 
  theme(aspect.ratio = 1) + 
  geom_hline(yintercept = -1, linetype= 'dashed', colour = 'black')  + 
  geom_vline(xintercept = -1, linetype= 'dashed', colour = 'black') + 
  xlab('Serotype_1 density normalized to lytA (log10)') + 
  ylab('Serotype_2 density normalized to lytA (log10)')
filename <- paste0('Fig3J_double_carriage.pdf')
ggsave(filename, width = 4, height = 4)



##############################################
#
#Fig 3K
#
#######################################################
#get the number of Spn serotypes found at the same time
Spn_number <- data.frame(table(mrg8$SampleID_Day, mrg8$InfectionLongerTrue)[,2])
colnames(Spn_number) <- 'SpnNumber'
mrg9 <- merge(mrg8, Spn_number, by.x= 'SampleID_Day', by.y=0)

table(mrg9$Dominant, mrg9$SpnNumber) 
prop.table(table(mrg9$Dominant, mrg9$SpnNumber), margin = 2 )


#only look at the ones with single-infection
mrg50 <- mrg9[mrg9$InfectionLongerTrue == 1 & 
                mrg9$SpnNumber  == 1  ,] 
table(mrg50$Dominant, mrg50$AcquisitionLonger)
table(mrg50$Dominant, mrg50$Assay)



# look at the ones with at least one-infection to see what is the proportion of dominant serotypes
mrg501 <- mrg9[mrg9$InfectionLongerTrue == 1 & 
                 mrg9$SpnNumber  > 0  ,] 
table(mrg501$Dominant, mrg501$Assay)
table(mrg501$Dominant, mrg501$SpnNumber)
dominant_prop <- data.frame(prop.table(table(mrg501$Dominant, mrg501$SpnNumber), margin = 2 )[2,])
colnames(dominant_prop) <- 'Proportion'
dominant_prop$Var2 <- 1:6
dominant_prop$Percentage <- round(dominant_prop$Proportion, 2)
dominant_prop$Y <- as.numeric(table(mrg501$SpnNumber))+10
dominant_number <- data.frame(table(mrg501$Dominant, mrg501$SpnNumber))


#only look at the ones with co-infection
mrg11 <- mrg9[mrg9$InfectionLongerTrue == 1 & 
                mrg9$SpnNumber  >1  ,] 

#look at whether some serotypes are more often dominant
table(mrg11$Dominant, mrg11$Assay)
orders <- names(table(mrg501$Assay))[order(table(mrg501$Assay), decreasing = F)]
dominant_types <- data.frame(table(mrg501$Dominant, mrg501$Assay, mrg501$SpnNumber))
dominant_types$Var2 <- factor(dominant_types$Var2, levels= orders)
ggplot(dominant_types, aes(y=Var2, x=Freq, fill=Var1)) + geom_bar(stat='identity') + 
  facet_grid(.~Var3) + scale_fill_manual(values = c('steelblue', 'black'))
filename <- paste0('Fig3K_dominance_per_spnnumber.pdf')
ggsave(filename, width = 8, height = 4)




dominant_types <- data.frame(prop.table(table(mrg11$Dominant, mrg11$Assay, mrg11$SpnNumber), margin = 2)) 


keep50<- names(which(table(mrg11$Assay) > 50))
mrg11_keep50 <- mrg11[mrg11$Assay %in% keep50,]

m <- glmer(as.factor(Dominant) ~ 0 + Assay + SpnNumber + (1|StudyID), data = mrg11_keep50, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)



##############################################
#
#Fig 3L
#
#######################################################
#look at time since acquisition
table(mrg11$Dominant, mrg11$DaysSinceAcquisition2, mrg11$SpnNumber)
table(mrg11$Dominant, mrg11$DaysSinceAcquisition2)
prop.table(table(mrg11$Dominant, mrg11$DaysSinceAcquisition2), margin = 2 )


mrg62 <- mrg11[mrg11$SpnNumber == 2,]
m <- glmer(as.factor(Dominant) ~ as.factor(DaysSinceAcquisition2)  + (1|Assay) +(1|StudyID), data = mrg62, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

duration_prop <- data.frame(prop.table(table(mrg62$Dominant, mrg62$DaysSinceAcquisition2), margin = 2))
duration_prop2 <- duration_prop[as.numeric(duration_prop$Var2) <4,]

ggplot(duration_prop2, aes(x=Var2, y=Freq, fill=Var1)) + geom_bar(stat='identity') + 
  scale_fill_manual(values = c('steelblue', 'black'))
filename <- paste0('Fig3L_dominance_per_timesinceacquisitiondoublecolonization.pdf')
ggsave(filename, width = 3, height = 2)



####################
#
#Figure 3M
#check for clearance and if it is linked to being dominant or not
##################################################################
mrg_clear7 <- mrg7
mrg_clear7$ClearNextDay <- 0
mrg_clear7$SpnNumber  <- NA

mrg_clear7.split <- split(mrg_clear7,  mrg_clear7$StudyPath) 

for(h in 1:length(mrg_clear7.split)){
  each <- data.frame(mrg_clear7.split[[h]])  #get each child vs assay out
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
    each4$SpnNumber[i] <- sum(sample_check2[grepl('Spn',sample_check2$Assay) & 
                                              !grepl('lytA',sample_check2$Assay) & 
                                              !grepl('piaB',sample_check2$Assay) & 
                                              !grepl('Spn6C/D',sample_check2$Assay), 'InfectionLongerTrue'])
    #fill in the covariates of interest
    
  }
  mrg_clear7.split[[h]] <- each4 #save the fixed file back in the list
}

mrg_clear8 <- do.call(rbind, mrg_clear7.split)
mrg_clear8$SpnNumber2 <- mrg_clear8$SpnNumber
mrg_clear8$SpnNumber2[mrg_clear8$SpnNumber>3] <- '4+'
mrg_clear_spn80 <- mrg_clear8[mrg_clear8$ClearanceLonger == 0,] #remove the clearancelonger point

mrg_clear_spn8 <- mrg_clear_spn80[mrg_clear_spn80$SpnNumber > 0,]
table(mrg_clear_spn8$Assay)
table(mrg_clear_spn8$ClearNextDay, mrg_clear_spn8$Dominant)


mrg_clear_spn_fin <- mrg_clear_spn8[is.finite(mrg_clear_spn8$LytANorm),]

#####make a model that incorporates density and dominance and the number of serotypes
m <- glmer(as.factor(ClearNextDay) ~ Dominant + SpnNumber2 + log10_avgconc + (1|StudyID), data = mrg_clear_spn_fin, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)


#make a plot of the model
OR_adjusted2 <- data.frame(summary(m)[[10]][-1,])
OR_adjusted2$Parameter <- gsub('Number2', '', rownames(OR_adjusted2))
OR_adjusted2$Ratio <- exp(OR_adjusted2$Estimate)
OR_adjusted2$Lower <- exp(OR_adjusted2$Estimate - OR_adjusted2$Std..Error)
OR_adjusted2$Upper <- exp(OR_adjusted2$Estimate + OR_adjusted2$Std..Error)
ggplot(OR_adjusted2, aes(x=Ratio, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing serotype next day') + 
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 'dashed')
filename <- paste0('Fig3N_Dominant_serotype_number_riskClearance_incDens.pdf')
ggsave(filename, width = 3, height = 2.3)




####
save.image('Z:/Projects/SAMSAM/Manuscripts\\Figures\\Figure3_bacterialstability/Figure3_bacteria.Rdata')
