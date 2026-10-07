#####
#This script creates UMAPS combining all the pathogen data and see if they associate
#with parameters
#look at effect of covariates on denstity

rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\FigS8\\")

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
library(knitr)
library(broom.mixed)

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()


#load in database and select key bacteria
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))

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

#create a function that fixes the display of p-values
disp_pval <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < 0.001,
      scales::scientific(p, digits = 3),
      format(round(p, 3), nsmall = 3)
    )
  )
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
#                 FIGURE S8A                      #
#                                                 #
###################################################

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
stats_period$Label <- paste0(
  "Afternoon: Est: ", round(stats_period$EstimatePM, 3),
  ", P = ", disp_pval(stats_period$PvalPM),
  "\n",
  "Evening: Est: ", round(stats_period$EstimateEve, 3),
  ", P = ", disp_pval(stats_period$PvalEve)
)

mrg3.period <- merge(mrg3, stats_period, by= 'Assay')
mrg3.period$Label2 <- mrg3.period$Label
mrg3.period$Label2[duplicated(mrg3.period$Label)] <- NA

stats_by_assay <- mrg3.period %>%
  group_by(Assay) %>%
  summarise(
    Label2 = first(Label2),   # or your own summary string
    .groups = "drop"
  )

# Merge back so each row has the per-assay label
mrg3.period <- mrg3.period %>%
  left_join(stats_by_assay, by = "Assay", suffix = c("", ".stat")) %>%
  mutate(
    AssayLabel = paste0(Assay, "\n", Label2.stat)  # use the joined column
  )

pl3 <- ggplot(mrg3.period, aes(x = as.factor(PeriodCollection), y = log10_avgconc)) + 
  facet_wrap(. ~ AssayLabel, scale = "free", ncol = 5) + 
  geom_violin(aes(fill = PeriodCollection)) +   
  geom_boxplot(width = 0.2) + 
  ggtitle("Period of collection, all kids/samples") +
  theme(
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text = element_text(colour = "black", size = 9),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "red")
  )

pl3

pdf_filename <- paste0("FigS8A_density_periodcollection_3periods_no_swithcing_times.pdf")
ggsave(pdf_filename, width = 12, height = 3)


#nr samples:
fig_data      <- mrg3.period
group_cols    <- c("Assay", "PeriodCollection")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "FigS8A_samplesize.txt")



###################################################
#                                                 #
#                 FIGURE S8B                      #
#                                                 #
###################################################


#select only M. cattharalis and Spn:
mrg2b <- mrg2a %>%
  filter(Assay %in% c("Haemophilus_influenzae", "Moraxella_catarrhalis", "Spn(lytA)", "Staphylococcus_aureus", "Streptococcus_pyogenes"))

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
  select(term, estimate, std.error, statistic, p.value)
table2

write_tsv(table2, "FigS8B_statssummary.txt")

disp_pval <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < 0.001,
      scales::scientific(p, digits = 3),
      format(round(p, 3), nsmall = 3)
    )
  )
}

stats_period2 <- bind_rows(stats_list) %>%
  mutate(
    adj.Pval = p.adjust(Pval, method = "BH"),
    Label = paste0(
      "Est.: ", round(Estimate, 3),
      ", P = ", disp_pval(Pval)
    )
  )

mrg3.period2 <- mrg3.episodes %>%
  left_join(stats_period2, by = "Assay") %>%
  group_by(Assay) %>%
  mutate(Label2 = ifelse(row_number() == 1, Label, NA_character_)) %>%
  ungroup()

pl6 <- ggplot(mrg3.period2, aes(x = Period2, y = log10_avgconc)) + 
  facet_wrap(. ~ Assay, scales = "free", ncol = 5) + 
  geom_violin(aes(fill = Period2)) +  
  geom_boxplot(width = 0.2) + 
  ggtitle("Period of collection AM/PM, switching times") +
  scale_y_continuous(
    expand = expansion(mult = c(0.0, 0.15))  # more space at top
  ) +
  geom_text(aes(label = Label2), y = Inf, x = 1.5, hjust =0.5, vjust = 1.2, colour = "red") +
  theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA))
pl6

pdf_filename <- paste0("FigS8B_density_AmvsPM", ".pdf")
ggsave(pdf_filename, width = 11, height = 3)


#nr samples:
fig_data      <- mrg3.period2
group_cols    <- c("Assay", "Period2")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "FigS8B_samplesize.txt")


# check samples used:
class(mrg3.period2$log10_avgconc)
table1 <- mrg3.period2
table2 <- table1 %>%
  dplyr::filter(Assay == "Moraxella_catarrhalis") %>%
  dplyr::select(StudyID, Day, Period2, episode)
table2

# Export nr table with samples used:
table3 <- mrg3.period %>%
  filter(Assay == "Moraxella_catarrhalis") %>%
  count(StudyID, name = "SamplesTotal") %>%
  left_join(
    mrg3.period2 %>%
      filter(Assay == "Moraxella_catarrhalis") %>%
      count(StudyID, name = "SamplesSelected"),
    by = "StudyID"
  ) %>%
  mutate(
    SamplesSelected = coalesce(SamplesSelected, 0L)
  ) %>%
  arrange(StudyID)
table3

write_tsv(table3, "FigS8B_AmvsPM_sampletable.txt")

