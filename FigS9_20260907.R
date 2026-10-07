#####
#This makes plots of exposure
#update 9-6. Fixed the forest plot for exposure vs activities to be of standard error and not CI

rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\FigS9\\")

# Library
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(reshape2)
library(readxl)
library(survival)
library(ggsurvfit)
library(lme4)
library(lmerTest)
library(ggmosaic)
library(reshape2)
library(circlize)
library(ComplexHeatmap)


#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#create a function that fixes the display of p-values
disp_pval <- function(p){
  if(p<0.001){
    p <- scales::scientific(p, digits= 3)
  } else {
    p <- round(p, 3)
  }
  return(p)
}

#load in database
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))


###################################################
#                                                 #
#                 FIGURE S9A                      #
#                                                 #
###################################################

mrg9 <- mrg %>%
  group_by(StudyID, Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLongerTrue) ==1 
  ) %>%
  ungroup()

#plots for the non-spn
notspn2 <- mrg9[grepl("^Spn\\(lytA\\)$", mrg9$Assay) | !grepl("^Spn", mrg9$Assay), ]
notspn2 <- notspn2 %>% filter(Assay != "EAV")  %>% filter(Assay != "16S")

ggplot(notspn2, aes(x = Day2, y = log10_avgconc)) +  # Convert Day to numeric for regression
  geom_point(aes(color = Assay, shape = AssayType), size = 1) +  # Points colored by Assay
  #geom_line(aes(color = Assay), size = 0.7) +  # Line connecting the points
  geom_smooth(aes(color = Assay, linetype = AssayType), method = "loess", se = FALSE, size = 0.7, span = .6) +    # Add separate smoothing line for each Assay
  theme_minimal() +
  facet_wrap(.~StudyID, ncol = 7) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "right"
  ) +  # Rotate x-axis labels
  labs(
    title = 'Bacterial/Viral carriage',  # Title for each StudyID
    x = NULL, y = "log10(average_conc)"  # y-axis label with log scale
  )
ggsave("FigS9A_bac_vir_atleast2.pdf", width = 12, height = 10)


save.image('FigS9A_cohort.Rdata')


###################################################
#                                                 #
#                 FIGURE S9B                      #
#                                                 #
###################################################

mrg1b <- mrg[mrg$Antibiotics == 'Yes',] #keep antibiotics users
unique(mrg1b$StudyID)
mrg1 <- mrg[mrg$StudyID == 'SAM006',] #select the child with the antibiotics

mrg1a <- mrg1[mrg1$Assay == 'Spn(lytA)' | 
                    mrg1$Assay == 'Spn(piaB)'| 
                    grepl("Haemophilus_influenzae", mrg1$Assay)| 
                    grepl("aureus", mrg1$Assay)|
                    grepl("Spn23A/B/F", mrg1$Assay)|grepl("16S", mrg1$Assay)|
                    grepl("bocavirus", mrg1$Assay)|
                    grepl("influenza_A", mrg1$Assay)|
                    grepl("aureus", mrg1$Assay)|
                    grepl("Moraxella_catarrhalis", mrg1$Assay),]


mrg1a <- mrg1a %>%
  group_by(Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    max(InfectionLonger) == 1
  ) %>%
  ungroup() #make a new dataframe and only select assays that are present

xmin <- min(mrg1a$Day2[mrg1a$Antibiotics == 'Yes'])
#xmax <- max(mrg1a$Day2[mrg1a$Antibiotics == 'Yes'])
xmax <- 12

ggplot(mrg1a, aes(x = Day2, y = log10_avgconc, colour = Assay)) +  
  geom_rect(
    xmin = xmin, xmax = xmax, ymin = 0, ymax = 7,
    colour = 'grey', fill = 'grey', alpha = 0.2
  ) +
  geom_point(aes(shape = AssayType)) + 
  geom_line(aes(group = Assay, linetype = AssayType)) + 
  scale_color_brewer(palette = 'Paired') +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.15))  # more space at top
  ) +
  facet_grid(. ~ AssayType) + 
  annotate(
    'text',
    x = mean(range(mrg1a$Day2, na.rm = TRUE)),  # center in x
    y = Inf,
    label = 'Antibiotics (Pheneticillin)',
    colour = 'black',
    hjust = 0.5, vjust = 1.2
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text = element_text(colour = "black"),
    legend.position = "right"
  )

ggsave('FigS9B_SAM006_lineplots.pdf', width =6, height =2.5)


