# %pos samples, Age and gender

rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\FigS6\\")

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
library(lme4)
library(scales)
library(glmmTMB)
library(broom.mixed)

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
#                 FIGURE S6A                      #
#                                                 #
###################################################

mrg_d01 <- subset(mrg, Day == "D01")
mrg_d01 <- mrg_d01 %>% filter(Assay != "16S")

presence_d01 <- aggregate(
  Positive > 0 ~ Assay,
  data = mrg_d01,
  FUN = mean
)

colnames(presence_d01)[2] <- "PctD01"

presence_d01$Assay <- factor(
  presence_d01$Assay,
  levels = presence_d01$Assay[order(presence_d01$PctD01)]
)

pl_d01 <- ggplot(presence_d01, aes(y = Assay, x = PctD01*100)) +
  geom_bar(stat = "identity", fill = "steelblue", colour = "black") +
  xlab("% samples positive on D01")
pl_d01

#now together:
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
presence_child <- data.frame(rowSums(presence>0) / ncol(presence))
presence_child$Assay <- rownames(presence_child)
colnames(presence_child)[1] <- 'PresenceAtLeastOnce'

presence_child$Assay <- factor(presence_child$Assay, 
                               levels = presence_child$Assay[order(presence_child$PresenceAtLeastOnce)]) 

plot_df <- merge(
  presence_child,
  presence_d01,
  by = "Assay",
  all = TRUE
)

colnames(plot_df)[colnames(plot_df) == "PresenceAtLeastOnce"] <- "PctAll"

FigS6A <- ggplot(plot_df, aes(y = Assay)) +
  geom_col(aes(x = PctAll*100), fill = "steelblue", colour = "black",
           position = "identity") +
  geom_col(aes(x = PctD01*100), fill = "red", colour = "black",
           position = "identity") +
  xlab("% positive") +
  ylab("Assay")
FigS6A

pdf_filename <- paste0("SupFig6A_Percd01vsall_", ".pdf")
ggsave(pdf_filename, width = 4, height = 6)



###################################################
#                                                 #
#         FIGURE S6B - %pos >1 per sex            #
#                                                 #
###################################################

#Presence child (check per child how many positives)
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

presence_child <- data.frame(rowSums(presence>0) / ncol(presence)) #how many detected per child
presence_child$Assay <- rownames(presence_child)
colnames(presence_child)[1] <- 'PresenceAtLeastOnce'

#per gender:
child_assay <- mrg %>%
  filter(Assay != "16S") %>%
  group_by(StudyID, Gender, Assay) %>%
  summarise(Positive = max(Positive > 0, na.rm = TRUE), .groups = "drop")

n_gender <- mrg %>%   # total number of children per gender
  distinct(StudyID, Gender) %>%
  count(Gender, name = "n_gender")

presence_gender <- child_assay %>%   # % positive at least once within each gender, per assay
  group_by(Gender, Assay) %>%
  summarise(
    n_pos = sum(Positive),
    .groups = "drop"
  ) %>%
  left_join(n_gender, by = "Gender") %>%
  mutate(
    PresenceAtLeastOnce = n_pos / n_gender * 100
  )

#reorder:
assay_order <- presence_child %>%
  arrange(PresenceAtLeastOnce) %>%
  pull(Assay)

presence_gender$Assay <- factor(presence_gender$Assay, levels = assay_order)


ggplot(presence_gender, aes(y = Assay, x = PresenceAtLeastOnce, fill = Gender)) +
  geom_col(position = position_dodge(width = 0.8), colour = "black") +
  xlab("% children positive at least once") +
  ylab("Assay") 

pdf_filename <- paste0("SupFig6B_ExposurePerAge_", ".pdf")
ggsave(pdf_filename, width = 4, height = 6)




###################################################
#                                                 #
#                 FIGURE S6C                      #
#                                                 #
###################################################

#make associations with sex
assays_keep <- c(
  "Spn(lytA)",
  "Haemophilus_influenzae",
  "Moraxella_catarrhalis",
  "Staphylococcus_aureus")

mrg2a <- mrg[
    grepl("^Spn\\(lytA\\)$", mrg$Assay) |
    grepl("Haemophilus", mrg$Assay) |
    grepl("Moraxella", mrg$Assay) |
    grepl("aureus", mrg$Assay),
]

mrg2 <- mrg2a %>%
  filter(log10_avgconc > 0) %>%
  mutate(
    Assay = factor(Assay, levels = assays_keep),
    Gender = factor(Gender),
    Gender = relevel(Gender, ref = "F")
  ) %>%
  droplevels()

# Check the reference level and group sizes
levels(mrg2$Gender)

mrg2 %>%
  count(Assay, Gender)

#stats:
stats <- data.frame(
  Estimate = numeric(0),
  Pval = numeric(0),
  Lower = numeric(0),
  Upper = numeric(0),
  n = integer(0),
  n_study = integer(0),
  note = character(0)
)

mrg2.split <- split(mrg2, mrg2$Assay)

for (i in seq_along(mrg2.split)) {
  
  test <- droplevels(mrg2.split[[i]])
  assay_name <- names(mrg2.split)[i]
  
  # Cannot estimate Gender effect if only one level is present
  if (length(unique(na.omit(test$Gender))) < 2) {
    
    res <- data.frame(
      Estimate = NA_real_,
      Pval = NA_real_,
      Lower = NA_real_,
      Upper = NA_real_,
      n = nrow(test),
      n_study = length(unique(test$StudyID)),
      note = "Only one Gender level"
    )
    
  } else {
    
    m <- lmer(
      log10_avgconc ~ Gender + AgeYears + (1 | StudyID),
      data = test,
      REML = TRUE
    )
    
    coef_table <- as.data.frame(summary(m)$coefficients)
    
    # Find the Gender coefficient by name, rather than assuming row 2
    gender_row <- grep("^Gender", rownames(coef_table))
    
    if (length(gender_row) != 1) {
      
      res <- data.frame(
        Estimate = NA_real_,
        Pval = NA_real_,
        Lower = NA_real_,
        Upper = NA_real_,
        n = nrow(test),
        n_study = length(unique(test$StudyID)),
        note = "Gender coefficient not uniquely identified"
      )
      
    } else {
      
      estimate <- coef_table$Estimate[gender_row]
      se <- coef_table$`Std. Error`[gender_row]
      pval <- coef_table$`Pr(>|t|)`[gender_row]
      df <- coef_table$df[gender_row]
      
      critical_value <- qt(0.975, df = df)
      
      res <- data.frame(
        Estimate = estimate,
        Pval = pval,
        Lower = estimate - critical_value * se,
        Upper = estimate + critical_value * se,
        n = nrow(test),
        n_study = length(unique(test$StudyID)),
        note = NA_character_
      )
    }
  }
  
  res$Assay <- assay_name
  stats <- rbind(stats, res)
}

stats_sex <- stats %>%
  dplyr::select(Assay, everything()) %>%
  mutate(
    Pval_adj = p.adjust(Pval, method = "BH"),
    Label = paste0(
      "Est:", round(Estimate, 3),
      ", P=", sapply(Pval_adj, disp_pval)
    )
  )

stats_sex
capture.output(stats_sex, file = "FigS6C_statssummary_new.txt")


#plot:
mrg2.sex <- mrg2 %>%
  left_join(
    stats_sex %>%
      dplyr::select(Assay, Estimate, Pval, Pval_adj, Lower, Upper, Label),
    by = "Assay"
  ) %>%
  mutate(
    Assay = factor(Assay, levels = assays_keep)
  )

annotation_data <- mrg2.sex %>%
  group_by(Assay) %>%
  summarise(
    Label = first(Label),
    y = max(log10_avgconc, na.rm = TRUE),
    .groups = "drop"
  )


FigS6C <- ggplot(
  mrg2.sex,
  aes(x = Gender, y = log10_avgconc)
) +
  facet_wrap(
    . ~ Assay,
    ncol = 2, scales = "free"
  ) +
  geom_violin(
    aes(
      group = Gender,
      fill = Gender,
      alpha = Gender
    )
  ) +
  geom_boxplot(
    aes(group = Gender),
    width = 0.2
  ) +
  ggtitle("Effect of sex, all kids/samples") +
  geom_text(
    data = annotation_data,
    aes(
      x = 1.5,
      y = Inf,
      label = Label
    ),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 1.5,
    colour = "red"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = NA,
      colour = NA
    ),
    strip.text = element_text(
      colour = "black"
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.35)
    )
  )

FigS6C

ggsave('FigS6C_Persex.pdf', width = 5, height = 5)

#nr samples:
fig_data      <- mrg2.sex
group_cols    <- c("Assay", "Gender")   
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "FigS6C_samplesize.txt")






###################################################
#                                                 #
#         FIGURE S6D - %pos >1 per age            #
#                                                 #
###################################################

# Make one row per Assay × child × age
Assay_age <- mrg %>%
  mutate(Assay = if_else(
      Assay == "Spn(lytA)",
      "S. pneumoniae",
      Assay)) %>%
      filter(!Assay %in% c("16S", "PIV1", "PIV2", "RSVA_B" )) %>%
  filter(!grepl("^Spn", Assay)) %>%
  group_by(Assay, StudyID, AgeYears) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(Positive > 0, na.rm = TRUE),
    n_negative = n_samples - n_positive,
    percent_positive = n_positive / n_samples,
    .groups = "drop"
  ) %>%
  mutate(
    AgeYears = as.numeric(as.character(AgeYears))
  )

#check sample size:
n_sample_per_assay <- Assay_age %>%
  dplyr::select(Assay, AgeYears, n_samples) 
n_sample_per_assay

# Assays to test
assays_to_run <- unique(Assay_age$Assay)

# One age-trend binomial GLMM per assay
results_list <- lapply(assays_to_run, function(a) {
  
  d <- Assay_age %>%
    filter(Assay == a)
  
  # Need variation in outcome and age for a meaningful model
  if (
    nrow(d) < 5 ||
    dplyr::n_distinct(d$AgeYears) < 2 ||
    sum(d$n_positive, na.rm = TRUE) == 0 ||
    sum(d$n_negative, na.rm = TRUE) == 0
  ) {
    return(tibble(
      Assay = a,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_
    ))
  }
  
  fit <- tryCatch(
    glmmTMB(
      cbind(n_positive, n_negative) ~ AgeYears + (1 | StudyID),
      data = d,
      family = binomial()
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      Assay = a,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_
    ))
  }
  
  broom.mixed::tidy(fit, exponentiate = FALSE) %>%
    dplyr::filter(
      effect == "fixed",
      term == "AgeYears"
    ) %>%
    dplyr::select(
      term,
      estimate,
      std.error,
      statistic,
      p.value
    ) %>%
    dplyr::mutate(
      Assay = a,
      .before = 1
    )
})

Assay_age_results <- bind_rows(results_list) 

disp_pval <- function(p){
  out <- ifelse(
    p < 0.001,
    scales::scientific(p, digits = 3),
    round(p, 3)
  )
  return(out)
}

Assay_age_results2 <- Assay_age_results %>%
       mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  mutate(across(c(p.value, p.adj),
      ~ disp_pval(.x),
      .names = "{.col}_fmt")) %>% 
      arrange(p.adj)
Assay_age_results2


Assay_age_results3 <- Assay_age_results2 %>%
  mutate(
    p.value = p.adj_fmt,
    stars = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    estimate  = sprintf("%.3f", estimate),
    std.error = sprintf("%.3f", std.error)
  ) %>%
  dplyr::select(
    Assay,
    estimate,
    std.error,
    p.value,
    stars
  ) 
Assay_age_results3

table_wide <- Assay_age_results3 %>%
  mutate(
    p.value = as.character(p.value),
    estimate    = as.character(estimate),
    std.error   = as.character(std.error),
    stars       = as.character(stars)
  ) %>%
  dplyr::select(Assay, estimate, std.error, p.value, stars) %>%
  pivot_longer(
    cols = c(estimate, std.error, p.value, stars),
    names_to = "stat",
    values_to = "value"
  ) %>%
  mutate(
    stat = factor(
      stat,
      levels = c("estimate", "std.error", "p.value", "stars")
    )
  ) %>%
  pivot_wider(
    names_from = Assay,
    values_from = value
  )
table_wide

write_tsv(table_wide, "FigS6D_statstable.txt")






