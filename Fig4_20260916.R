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

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\Fig4\\")

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
library(broom.mixed)
library(parallel)

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#date:
timestamp <- format(Sys.time(), "%Y%m%d")

#load in database
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))

#enocde the season as a cosinus
mrg$DateCollection_EditDate #check what the date looks like
mrg$DaysSinceJan01 <- as.integer(format(as.Date(mrg$DateCollection_EditDate, format = "%Y_%m_%d"), "%j")) - 1 #set it to an integer
mrg$DaysCosinus <-  cos(2*pi*mrg$DaysSinceJan01/365) #change it to a cosinus

#filter out antibiotics samples
mrg1 <- mrg <- mrg %>% 
  filter(!grepl("Yes", Antibiotics))   #remove anitbiotics samples

mrg2b <- mrg2 <- mrg1 %>% #remove pathogens/child combinations where there is no infection
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
#                 FIGURE 4A                       #
#                                                 #
###################################################

mrg9 <- mrg %>%
  group_by(StudyID, Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLongerTrue) ==1 
  ) %>%
  ungroup()

spn2 <- mrg9 %>% filter(grepl("^Spn", Assay))  %>% 
  filter(!grepl("18A", Assay)) %>% filter(!grepl("19B", Assay))

#Fig 4A_interaction  Spn serotypes
ggplot(spn2[spn2$StudyID %in% c("SAM-K1-22", "SAM-K1-36",  "SAM-K1-38","SAM010") &
              spn2$Assay != "Spn(piaB)", ], 
       aes(x = Day2, y = log10_avgconc)) +  # Convert Day to numeric for regression
  geom_point(aes(color = Assay), size = 1) +  # Points colored by Assay
  #geom_line(aes(color = Assay), size = 0.7) +  # Line connecting the points
  geom_smooth(aes(color = Assay), method = "loess", se = FALSE, size = 0.7, span = .6) +    # Add separate smoothing line for each Assay
  theme_minimal() +
  facet_wrap(.~StudyID, ncol = 4) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "right"
  ) +  # Rotate x-axis labels
  labs(
    title = 'S. pneumoniae carriage',  # Title for each StudyID
    x = NULL, y = "log10(average_conc)"  # y-axis label with log scale
  )
ggsave("Fig4A_spn_atleast2_SpnComp.pdf", width = 6, height = 1.7)



###################################################
#                                                 #
#                 FIGURE 4B                       #
#                                                 #
###################################################

mrg_Spn <- mrg[grepl('Spn', mrg$Assay),]
lytapiaBPos <- mrg_Spn[
  (mrg_Spn$Assay == "Spn(lytA)" & mrg_Spn$InfectionLongerTrue == 1) |
    (mrg_Spn$Assay == "Spn(piaB)" & mrg_Spn$InfectionLongerTrue == 1),
]
mrg_Spn2 <- mrg_Spn[
  mrg_Spn$SampleID_Day %in% lytapiaBPos$SampleID_Day &
    mrg_Spn$Assay != "Spn(lytA)" &
    mrg_Spn$Assay != "Spn(piaB)" &
    mrg_Spn$Assay != "Spn6C/D",
]

table(mrg_Spn2$Assay)
table(mrg_Spn2$Assay, mrg_Spn2$InfectionLongerTrue)
Total_serotype <- data.frame(table(table(mrg_Spn2$SampleID_Day, mrg_Spn2$InfectionLongerTrue)[,2]))
Total_serotype$Perc <- round(Total_serotype$Freq / sum(Total_serotype$Freq) * 100, 1)

#plot:
Fig4B <- ggplot(Total_serotype, aes(y=Var1, x=Freq)) + 
  geom_bar(stat = 'identity', fill="steelblue", colour="black") + 
  scale_x_continuous(expand = c(0,0)) +
  geom_text(
    aes(
      label = Freq,
      x = ifelse(Freq < 26, Freq + 5, Freq - 5),
      hjust = ifelse(Freq < 26, 0, 1),
      colour = ifelse(Freq < 26, "small", "large")
    )) +
  scale_colour_manual(values = c("small" = "black", "large" = "white"), guide = "none") +
  ylab('Colonizing serotypes') +  
  xlab("Frequency") +
  coord_cartesian(clip = "off")
Fig4B

#save:
timestamp <- format(Sys.time(), "%Y%m%d")
Fig4B
filename <- paste0("Fig4B_", timestamp, ".pdf")
ggsave(filename , width = 3, height = 2)

sum(Total_serotype$Freq[3:7])



###################################################
#                                                 #
#                 FIGURE 4C                       #
#                                                 #
###################################################

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

#stats:
survdiff(Surv(Duration, status) ~ Assay2, data = lengthSpnLonger)
summary(survfit(Surv(Duration, status) ~ Assay2, data = lengthSpnLonger), times = 27)
summary(survfit(Surv(Duration, status) ~ Assay2, data = lengthSpnLonger))

sd <- survdiff(Surv(Duration, status) ~ Assay2, data = lengthSpnLonger)
p_val <- pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
p_txt <- formatC(p_val, format = "e", digits = 2)

#plot:
Fig4C <- survfit2(Surv(Duration, status) ~ Assay2, data = lengthSpnLonger) %>% 
  ggsurvfit() +
  labs(
    title = paste0("P = ", p_txt),
    x = "Days",
    y = "Overall survival probability Spn"
  ) + 
  add_confidence_interval() +
  xlim(c(0, 27))
Fig4C

#save:
filename <- 'Fig4C_survival_curves_Spn.pdf'
pdf(filename, width = 2.5, height =3)
Fig4C
dev.off()


#nr samples:
fig_data      <- lengthSpnLonger
group_cols    <- c("Assay2", "Duration")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig4c_samplesize.txt")







###################################################
#                                                 #
#                 FIGURE 4D                       #
#                                                 #
###################################################

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

#stats:
m <- glmer(as.factor(ClearNextDay) ~ AgeYears + SpnNumber2 + (1|Assay) + (1|StudyID), data = mrg_clear_spn, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

summary(m)[[10]][,1:2]
OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$SpnNumber <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)

names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))
write_tsv(OR_adjusted, "Fig4D_statssummary.txt")

#plot:
ggplot(OR_adjusted, aes(x = OR, y = SpnNumber)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower),
                 size = 0.5, height = 0.2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") + 
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  xlab("Adjusted OR of clearing \na serotype next day") + 
  geom_vline(xintercept = 1, linetype = "dashed")
ggsave('Fig4D_AdjustedOR_Spn_clearanceNextDat.pdf', width = 3, height = 2.5)



###################################################
#                                                 #
#                 FIGURE 4E                       #
#                                                 #
###################################################

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

mrg10d$Dominant2 <- dplyr::recode(
  mrg10d$Dominant2,
  `both` = "Both",
  `none` = "None",
  `Serotype_1` = "Serotype A",
  `Serotype_2` = "Serotype B"
)

ggplot(mrg10d, aes(x=log10(LytANorm_1+.0000001), y = log10(LytANorm_2+.0000001))) + geom_point(aes(colour = Dominant2)) + 
  theme(aspect.ratio = 1) + 
  geom_hline(yintercept = -1, linetype= 'dashed', colour = 'black')  + 
  geom_vline(xintercept = -1, linetype= 'dashed', colour = 'black') + 
  xlab('Serotype A density \nnormalized to lytA (log10)') + 
  ylab('Serotype B density \nnormalized to lytA (log10)')

filename <- paste0('Fig4E_double_carriage.pdf')
ggsave(filename, width = 4, height = 4)



###################################################
#                                                 #
#                 FIGURE 4F                       #
#                                                 #
###################################################

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

ggplot(dominant_types[dominant_types$Var3 != 1,], aes(y=Var2, x=Freq, fill=Var1)) + geom_bar(stat='identity') + 
  facet_grid(.~Var3) + scale_fill_manual(values = c('steelblue', 'black'))

filename <- paste0('Fig4F_dominance_per_spnnumber.pdf')
ggsave(filename, width = 8, height = 4)

dominant_types <- data.frame(prop.table(table(mrg11$Dominant, mrg11$Assay, mrg11$SpnNumber), margin = 2)) 

keep50<- names(which(table(mrg11$Assay) > 50))
mrg11_keep50 <- mrg11[mrg11$Assay %in% keep50,]

m <- glmer(as.factor(Dominant) ~ 0 + Assay + SpnNumber + (1|StudyID), data = mrg11_keep50, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

table2 <- tidy(m) %>%
  filter(effect == "fixed") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value)
table2

write_tsv(table2, "Fig4F_statssummary.txt")


###################################################
#                                                 #
#                 FIGURE 4G                       #
#                                                 #
###################################################

#look at time since acquisition
table(mrg11$Dominant, mrg11$DaysSinceAcquisition2, mrg11$SpnNumber)
table(mrg11$Dominant, mrg11$DaysSinceAcquisition2)
prop.table(table(mrg11$Dominant, mrg11$DaysSinceAcquisition2), margin = 2 )

mrg62 <- mrg11[mrg11$SpnNumber == 2,]
m <- glmer(as.factor(Dominant) ~ as.factor(DaysSinceAcquisition2)  + (1|Assay) +(1|StudyID), data = mrg62, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

table2 <- tidy(m) %>%
  filter(effect == "fixed") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value)
table2
write_tsv(table2, "Fig4G_statssummary.txt")

duration_prop <- data.frame(prop.table(table(mrg62$Dominant, mrg62$DaysSinceAcquisition2), margin = 2))
duration_prop2 <- duration_prop[as.numeric(duration_prop$Var2) <4,]

ggplot(duration_prop2, aes(x=Var2, y=Freq, fill=Var1)) + geom_bar(stat='identity') + 
  scale_fill_manual(values = c('transparent', 'black')) +
  ylab("Chance of becoming the \ndominant S. pneumoniae serotype") +
  xlab("Days since acquisition") +
  theme(legend.position = "none")

filename <- paste0('Fig4G_dominance_per_timesinceacquisitiondoublecolonization.pdf')
ggsave(filename, width = 3, height = 2)



###################################################
#                                                 #
#                 FIGURE 4H                       #
#                                                 #
###################################################

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
unique(mrg_clear_spn_fin$Assay)

#####make a model that incorporates density and dominance and the number of serotypes
m <- glmer(as.factor(ClearNextDay) ~ Dominant + SpnNumber2 + log10_avgconc + (1|StudyID), data = mrg_clear_spn_fin, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

#make a plot of the model
OR_adjusted2 <- data.frame(summary(m)[[10]][-1,])
OR_adjusted2$Parameter <- gsub('Number2', '', rownames(OR_adjusted2))
OR_adjusted2$Ratio <- exp(OR_adjusted2$Estimate)

#95% Wald CIs:
OR_adjusted2$Lower <- exp(OR_adjusted2$Estimate - 1.96 * OR_adjusted2$`Std..Error`)
OR_adjusted2$Upper <- exp(OR_adjusted2$Estimate + 1.96 * OR_adjusted2$`Std..Error`)

names(OR_adjusted2)
OR_adjusted2$stars <- with(OR_adjusted2,
                           ifelse(Pr...z.. < 0.001, "***",
                                  ifelse(Pr...z.. < 0.01,  "**",
                                         ifelse(Pr...z.. < 0.05,  "*", ""))))
write_tsv(OR_adjusted2, "Fig4H_statssummary.txt")


#plot:
ggplot(OR_adjusted2, aes(x=Ratio, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  xlab('OR of clearing serotype next day') + 
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 'dashed')

#save:
filename <- paste0('Fig4H_Dominant_serotype_number_riskClearance_incDens.pdf')
ggsave(filename, width = 3, height = 2.3)



###################################################
#                                                 #
#                 FIGURE 4I                       #
#                                                 #
###################################################

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
unique(lengthInfections2Longer$Assay)

#select only viral assays:
unique_viral_assays <- unique(mrg$Assay[mrg$AssayType == "Viral"])
unique_viral_assays
lengthInfections3Longer <- lengthInfections2Longer %>%
  filter (Assay %in% unique_viral_assays)
unique(lengthInfections3Longer$Assay)

Surv(lengthInfections3Longer$Duration, lengthInfections3Longer$status)
s1 <- survfit(Surv(Duration, status) ~ 1, data = lengthInfections3Longer)
str(s1)

pl3 <- survfit2(Surv(Duration, status) ~ 1, data = lengthInfections3Longer) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +   xlim(c(0,27)) 
pl3

s <- summary(survfit2(Surv(Duration, status) ~ 1, data = lengthInfections2Longer))
s
txt <- capture.output(print(s))
writeLines(txt, "Fig4I_statssummary.txt")

pdf_filename <- paste0("Fig4I_survival_virus_days_", timestamp, ".pdf")
pdf(pdf_filename, width = 3, height = 3)
print(pl3)
dev.off()


###################################################
#                                                 #
#                 FIGURE 4J                       #
#                                                 #
###################################################

mrg2 <- mrg

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

unique(mrg_clear2$Assay)

#select only viral:
unique_viral_assays <- unique(mrg$Assay[mrg$AssayType == "Viral"])
unique_viral_assays

mrg_clear3 <- mrg_clear2 %>%
  filter (Assay %in% unique_viral_assays)
unique(mrg_clear3$Assay)

#stats:
m <- glmer(as.factor(ClearNextDay) ~ VirusNumber + Staph  + AgeYears + (1|Assay) + (1|StudyID), data = mrg_clear3, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Parameter <- gsub('SpnNumber2', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)
names(OR_adjusted)

OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))
write_tsv(OR_adjusted, "Fig4J_statssummary.txt")

ggplot(OR_adjusted, aes(x=OR, y=Parameter)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  xlab('OR of clearing next day') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") 
  #scale_x_log10() 
ggsave('Fig4J_AdjustedORClearanceNextDay_virus.pdf', width = 3, height = 2)




###################################################
#                                                 #
#                 FIGURE 4K                       #
#                                                 #
###################################################

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

virus_cols <- c(
  "enterovirus", "HRV", "influenza_A", "bocavirus", "PIV4", "influenza_B",
  "sars-cov-2", "coronavirus_NL63", "coronavirus_OC43",
  "adenovirus", "PIV3", "coronavirus_HKU", "coronavirus_229E"
)
virus_cols <- intersect(virus_cols, names(mrg_all3))
mrg_all3$InfectionType[rowSums(mrg_all3[, virus_cols]) > 1] <- "z.multiple_viruses"

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
unique(mrg_bacteria$Assay)

m <- glmer(as.factor(InfectionLongerTrue) ~ as.factor(InfectionAnyVirus) + AgeYears + (1|StudyID), 
           data = mrg_bacteria,    family = binomial, control = glmerControl(optimizer = "bobyqa"))
anyvirus <- summary(m)[[10]][2,]
summary(m)
capture.output(summary(m), file = "Fig4K_statssummary_anyvirus.txt")
#anyvirus_conf <- confint(m)

mrg_bacteria2 <- mrg_bacteria[ mrg_bacteria$InfectionTypeBroad != 'PIV' &
                                 mrg_bacteria$InfectionTypeBroad != 'adenovirus'  ,]
unique(mrg_bacteria2$Assay)

m1<- glmer(as.factor(InfectionLongerTrue) ~ as.factor(InfectionTypeBroad) + AgeYears+ (1|StudyID), 
           data = mrg_bacteria2,    family = binomial, control = glmerControl(optimizer = "bobyqa"))
spec <- summary(m1)[[10]][2:6,]

summary(m1)
capture.output(summary(m1), file = "Fig4K_statssummary_all.txt")

#########these tables are part of the figure
table1 <- table(mrg_bacteria$InfectionTypeBroad, mrg_bacteria$InfectionLongerTrue)
table1
capture.output(table1, file = "Fig4K_tableall.txt")
table2 <- table(mrg_bacteria$InfectionAnyVirus, mrg_bacteria$InfectionLongerTrue)
table2
capture.output(table2, file = "Fig4K_tableanyvirus.txt")

prop.table(table(mrg_bacteria$InfectionTypeBroad, mrg_bacteria$InfectionLongerTrue), margin = 1)
prop.table(table(mrg_bacteria$InfectionAnyVirus, mrg_bacteria$InfectionLongerTrue), margin = 1)

OR_adjusted <- data.frame(rbind(anyvirus, spec))
OR_adjusted$Virus <- gsub('as.factor[(]InfectionTypeBroad[)]', '', rownames(OR_adjusted))
OR_adjusted$Perc <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)

OR_adjusted$Virus <- factor(OR_adjusted$Virus, levels = rev(c('anyvirus', 'entero_rhino_virus',
                                                              'bocavirus', 'z.multiple_viruses',
                                                              'coronavirus', 'influenza',
                                                              'adenovirus', 'PIV')))
names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))
write_tsv(OR_adjusted, "Fig4K_statssummary.txt")


Fig4K <- ggplot(OR_adjusted, aes(x=Perc, y=Virus)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  xlab('Adjusted OR of staph aureus detection') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_x_log10() 
Fig4K

ggsave('Fig4K_OR of staph aureus detection.pdf', width = 3, height = 2)



###################################################
save.image('Z:/Projects/SAMSAM/Manuscripts\\Figures\\Figure3_bacterialstability/version_newflu/Figure3_bacteria.Rdata')

