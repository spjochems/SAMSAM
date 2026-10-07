#####
#This script creates UMAPS combining all the pathogen data and see if they associate
#with parameters
#look at effect of covariates on denstity

rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\Fig1\\")

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
library(lme4)  
library(broom)
library(broom.mixed)

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#load in database and select key bacteria
mrg <- read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx')

#get out numbers for tables
table(mrg$Assay)

mrg_ex <- mrg[!duplicated(mrg$StudyID),]
table(mrg_ex$AgeYears, mrg_ex$Gender)
table(mrg_ex$AgeYears, mrg_ex$StudyID)

mrg_ex2 <- mrg[mrg$Assay == '16S',]
table(table(mrg_ex2$StudyID))

table(mrg_ex2$QuestionereFill)
table(mrg_ex2$QuestionereFill, mrg_ex2$StudyID)
table(mrg_ex2$SchoolDaycare)

table(mrg_ex2$PeriodCollection)

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
#                 FIGURE 1B                       #
#                                                 #
###################################################

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

#median(presence_melt2$value[presence_melt2$Var1 == "Streptococcus_pyogenes"], na.rm = TRUE)*100
#median(presence_melt2$value[presence_melt2$Var1 == "Staphylococcus_aureus"], na.rm = TRUE)*100
#median(presence_melt2$value[presence_melt2$Var1 == "Spn(lytA)"], na.rm = TRUE)*100

Fig1B <- plot_grid(pl1, pl2)
Fig1B
timestamp <- format(Sys.time(), "%Y%m%d")

pdf_filename <- paste0("Fig1B_presandfreq_", timestamp, ".pdf")
ggsave(pdf_filename, width = 20, height = 10)

mrg_ex3 <- mrg[!duplicated(mrg$SampleID_Day),]
ques <- mrg_ex3[mrg_ex3$QuestionereFill == 'YesQFill',]
data.frame(table(ques$SchoolDaycare, ques$StudyID)[2,])

###################################################
#                                                 #
#                 FIGURE 1C                       #
#                                                 #
###################################################

#make associations with age
mrg2a <- mrg[  grepl("Haemophilus", mrg$Assay) |
               grepl("Moraxella", mrg$Assay)|
               grepl("^Spn\\(lytA\\)$", mrg$Assay) |
               grepl("aureus", mrg$Assay),]
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

summary(m)

table2 <- tidy(m) %>%
  filter(effect == "fixed") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value)
table2

capture.output(table2, file = "Fig1c_statssummary.txt")


stats_age <- stats[-1,]
rownames(stats_age) <- stats_age$Assay <- names(mrg2.split)
stats_age$Label = paste0('Linear Est: ', round(stats_age$Estimate, 3), ', Pval: ', scales::scientific(stats_age$Pval, 2))

#BH correction
mrg2.age <- merge(mrg2, stats_age, by = "Assay")
mrg2.age$Pval_adj <- p.adjust(mrg2.age$Pval, method = "BH")
mrg2.age$Label <- paste0(
  "Est:", round(mrg2.age$Estimate, 3),
  ", P=", sapply(mrg2.age$Pval_adj, disp_pval)
)
mrg2.age$Label2 <- mrg2.age$Label
mrg2.age$Label2[duplicated(mrg2.age$Label2)] <- NA

mrg2.age$Assay <- factor(mrg2.age$Assay,
                         levels = c("Haemophilus_influenzae", "Moraxella_catarrhalis","Spn(lytA)", "Staphylococcus_aureus"))

Fig1C <- ggplot(mrg2.age, aes(x = AgeYears, y = log10_avgconc)) +  
  facet_wrap(. ~ Assay,  ncol=4) + 
  geom_violin(aes(group = AgeYears, fill=Assay, alpha=AgeYears)) +   
  geom_boxplot(aes(group = AgeYears), width = 0.2) +  
  ggtitle("Effect of age, all kids/samples") + 
  geom_text(aes(label = Label2), y = Inf, x = Inf, hjust = 1.1, vjust = 1, colour = "red") +
  theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA))
Fig1C

pdf_filename <- paste0("Fig1C_density_age_new", ".pdf")
ggsave(pdf_filename, width= 8, height = 3)

#nr samples:
fig_data      <- mrg2.age
group_cols    <- c("Assay", "AgeYears")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig1c_samplesize.txt")


###################################################
#                                                 #
#                 FIGURE 1D                       #
#                                                 #
###################################################

#select only M. cattharalis and Spn:
mrg2b <- mrg2a %>%
  filter(Assay %in% c("Moraxella_catarrhalis", "Spn(lytA)"))

#effect of period of collection
mrg3 <- mrg2b[mrg2b$PeriodCollection != 'NA', ]
mrg3$PeriodCollection <- factor(mrg3$PeriodCollection,
                                levels = c('Morning', 'Afternoon', 'Evening'))
mrg3$Period2 <- 'PM'
mrg3$Period2[mrg3$PeriodCollection == 'Morning'] <- 'AM'

ex <- mrg3[mrg3$Assay == 'Moraxella_catarrhalis',]
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
mrg3.episodes <- mrg3.selected %>%
  mutate(Day_num = as.numeric(Day2)) %>%
  arrange(StudyID, Assay, Day_num) %>%
  group_by(StudyID, Assay) %>%
  mutate(
    day_gap = Day_num - lag(Day_num),
    new_episode = is.na(day_gap) | day_gap > 1,
    episode = cumsum(new_episode)
  ) %>%
  ungroup()

stats <- data.frame(Estimate = NA, Pval = NA)

mrg3.split <- split(mrg3.episodes, f = mrg3.episodes$Assay)

stats_list <- vector("list", length(mrg3.split))

for(i in seq_along(mrg3.split)){
  test <- mrg3.split[[i]]
  
  m <- lmer(log10_avgconc ~ Period2  + (1|StudyID:episode), data = test)
  coefs <- lmerTest:::get_coefmat(m)
  
  stats_list[[i]] <- data.frame(
    Assay = names(mrg3.split)[i],
    Estimate = coefs[2, 1],
    Pval = coefs[2, 5],
    stringsAsFactors = FALSE
  )
}

summary(m)

table2 <- tidy(m) %>%
  filter(effect == "fixed") %>%
 dplyr:: select(term, estimate, std.error, statistic, p.value)
table2

write_tsv(table2, "Fig1D_statssummary.txt")


stats_period2 <- bind_rows(stats_list) %>%
  mutate(
    adj.Pval = p.adjust(Pval, method = "BH"),
    Label = paste0(
      "Est.:", round(Estimate, 3),
      ", P=", sapply(Pval, disp_pval)
    )
  )

mrg3.period2 <- mrg3.episodes %>%
  left_join(stats_period2, by = "Assay") %>%
  group_by(Assay) %>%
  mutate(Label2 = ifelse(row_number() == 1, Label, NA_character_)) %>%
  ungroup()

pl6 <- ggplot(mrg3.period2, aes(x = Period2, y = log10_avgconc)) + 
  facet_wrap(. ~ Assay, scales = "free", ncol = 4) + 
  geom_violin(aes(fill = Period2)) +  
  geom_boxplot(width = 0.2) + 
  ggtitle("Period of collection AM/PM, switching times") + 
  geom_text(aes(label = Label2), y = Inf, x = Inf, hjust = 1.2, vjust = 1, colour = "red") +
  theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA))
pl6

pdf_filename <- paste0("Fig1D_density_AmvsPM", ".pdf")
ggsave(pdf_filename, width = 5, height = 3.5)

#nr samples:
fig_data      <- mrg3.period2
group_cols    <- c("Assay", "Period2")   # 2–4 columns as needed
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig1D_samplesize.txt")



# Export table with samples used:
class(mrg3.period2$log10_avgconc)

table1 <- mrg3.period2

table2 <- table1 %>%
  dplyr::filter(Assay == "Moraxella_catarrhalis") %>%
  dplyr::select(StudyID, Day, Period2, episode, log10_avgconc)

write_tsv(table2, "Fig1D_AmvsPM_sampletable.txt")



#Samples per pathogen now
dplyr::filter(mrg3.period2, Assay == "Moraxella_catarrhalis") %>% 
  dplyr::tally() %>% 
  dplyr::pull()

#Samples per pathogen start
dplyr::filter(mrg, Assay == "Moraxella_catarrhalis") %>% 
  dplyr::tally() %>% 
  dplyr::pull()


