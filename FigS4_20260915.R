# Author:Youvika Singh

# Script details:
# script merges database lyta conventional and biomark results 
# and comapres Lyta conventional concentration and biomark 
# concentration

# Last edit: Youvika Singh
# edited for ct value replacement for unidentified (999.000) to 40
# NEW CHIP 8  and subsequent chips, and updated database


rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\FigS4\\")
output_dir <- getwd()

# Library
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(reshape2)
library(scales)
library(psych)
library(irr)
library(patchwork)
library(scales)

#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#chip <- read.delim("analysis/R/output/all_biomark_20250402.tsv")
#colnames(chip)[1] <- "SampleID_Day"

# Select appropriate file to work
chip <- read.delim("\\\\vf-lucid-r-i.lumcnet.prod.intern/lucid-r-i$//Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/0.20250128_main_Biomark_SAMSaliva_Allchips_YS/analysis/R/output/all_biomark_20250722.tsv")
mrg <- data.frame(read_xlsx("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx"))


#chip <- chip %>% filter(!grepl("Chip8", chip))

#create a function that fixes the display of p-values
disp_pval <- function(p){
  if(p<0.001){
    p <- scales::scientific(p, digits= 2)
  } else {
    p <- round(p, 2)
  }
  return(p)
}



###################################################
#                                                 #
#                 FIGURE S4A                      #
#                                                 #
###################################################

spn <- chip %>% filter(grepl("^Spn", Assay))

spn <- spn %>%
  filter(grepl("^SAM", Sample))  %>%
  filter(!grepl("NPS$", Sample)) %>%
  filter(!grepl("_S", Sample))  %>%
  filter(!grepl("MOCK", Sample))  

spn <- spn %>%
  separate(Sample, into = c("Sample", "Day"), sep = "_")

unique(spn$Day)
spn %>%
  filter(is.na(Day)) %>%
  distinct(Sample)

spn <- spn %>%
  filter(!grepl("N15\\+SAM-K1-13", Day)) %>%  # escape '+'
  filter(!is.na(Day))
unique(spn$Day)
unique(spn$Sample)

spn <- spn %>%
  mutate(Is_Positive = ifelse(Average_conc > 0, 1, 0))

str(spn)

spn_filt <-  spn %>%
  group_by(Sample, Day) %>%
  filter(
    all(Is_Positive[Assay %in% c("Spn(lytA)", "Spn(piaB)")] == 0) &  # Both lytA and piaB are 0
      any(Is_Positive[!Assay %in% c("Spn(lytA)", "Spn(piaB)")] == 1)   # Any other assay is 1
  ) %>%
  ungroup()
str(spn)


#remove the saliva samples from the serotype specificity
spn_filt <- spn_filt[!grepl('S', spn_filt$Day),]


# Summarize data by remaining serotypes
serotype_summary <- spn_filt %>%
  group_by(Assay) %>%
  summarize(
    Total_Positive = sum(Is_Positive),
    Total_Samples = n(),
    Total_Flag = sum(Flag)  # Sum of Flag column
  ) %>%
  mutate(
    Specificity = Total_Positive / Total_Samples * 100
  )

# Filter rows with Specificity != 0
serotype_summary <- serotype_summary %>% filter(Specificity != 0.0)

# Barplot with Flag counts added as labels
p <- ggplot(serotype_summary, aes(x = reorder(Assay, Specificity), 
                                  y = Specificity, fill = Assay)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = paste0("(N=", Total_Positive, ")")),          # or whatever your count column is called
    position = position_identity(),
    hjust = -0.1,                # just outside the bar end
    size = 3.5
  ) +
  theme_minimal() +
  labs(
    title = "Serotypes specificity",
    x = "Assay",
    y = "% times observed in absence of lytA/piaB"
  ) +
  theme(
    axis.text.x = element_text(),
    legend.position = "none"
  ) +  
  coord_flip()
print(p)

timestamp <- format(Sys.time(), "%Y%m%d")

pdf_filename <- paste0(output_dir,"/FigS4A_bar_specificity_spn_serotypes_", timestamp, ".pdf")
pdf(pdf_filename, width = 4, height = 3)
p
dev.off()



###################################################
#                                                 #
#                 FIGURE S4B                      #
#                                                 #
###################################################

#Figure to see how often concentration is above the lyta
#for each assay in serotype_summary$Assay, go to the spn object and then for each sample/day subtract the Average_Ct value from 
#the Spn(lytA) assay from that assay and return the object 

spn$Conc_log <- log10(spn$Average_conc+1)

spn_filt2 <-  spn %>%
  group_by(Sample, Day) %>%
  filter(
      any(Is_Positive[!Assay %in% c("Spn(lytA)", "Spn(piaB)")] == 1)   # Any other assay is 1
  ) %>%
  filter(
    any(Is_Positive[Assay %in% c("Spn(lytA)")] == 1)   # Any other assay is 1
  ) %>%
  ungroup()
str(spn)

unique(
  spn_filt2$Is_Positive[
    spn_filt2$Assay == "Spn(lytA)"
  ]
)

too_high <- data.frame(ncol=2, nrow = length(serotype_summary$Assay))
i=1

for(assay in serotype_summary$Assay){
  print(assay)
  
  #filter from spn_filt2 the rows where Assay is "Spn(lytA)" or the assay
  spn_assay <- spn_filt2 %>%
    filter(Assay %in% c("Spn(lytA)", assay)) %>%
    dplyr::select(Sample, Day, Assay, Conc_log) %>%
    pivot_wider(names_from = Assay, values_from = Conc_log)
  
  too_high[i,] <- c(assay, sum(spn_assay[,4] - spn_assay[,3] > 1))   #Times Spn log_conc exceeds the lytA > 1
  i=i+1
}

too_high <- too_high %>%
  rename(Specificity = ncol,
         nr_obs = nrow) %>%
  mutate(
    Specificity = factor(
      Specificity,
      levels = Specificity[order(nr_obs, decreasing = FALSE)]
    )
  ) %>%
  mutate(nr_obs = as.integer(nr_obs))

#make %
pos_per_assay <- spn_filt2 %>%
  filter(Is_Positive == 1) %>%
  count(Assay, name = "n_pos")

too_high2 <- too_high %>%
  rename(too_high_count = nr_obs) %>%
  mutate(too_high_count = as.integer(too_high_count)) %>%
  inner_join(pos_per_assay, by = c("Specificity" = "Assay")) %>%
  mutate(
    n_pos = as.integer(n_pos),
    too_high_pct = ifelse(n_pos > 0, 100 * too_high_count / n_pos, NA_real_),
    Specificity = factor(
      Specificity,
      levels = Specificity[order(too_high_count, decreasing = FALSE)]
    )
  ) 

# Levels and colors from the first plot
assay_levels <- unique(serotype_summary$Assay)
n_colors <- length(assay_levels)
pal_colors <- hue_pal()(n_colors)
assay_cols <- setNames(pal_colors, assay_levels)


p2 <- ggplot(too_high2, aes(x = Specificity, y = too_high_count, fill = Specificity)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Serotypes specificity",
    x = "Assay",
    y = "Times serotype/group concentration is 10-fold higher then lytA"
  ) +
  theme(
    axis.text.x = element_text(),
    legend.position = "none"
  ) +  
  coord_flip() +
  scale_fill_manual(values = assay_cols) 
print(p2)

timestamp <- format(Sys.time(), "%Y%m%d")

pdf_filename <- paste0(output_dir,"/FigS4B_bar_SpnVsLytA_", timestamp, ".pdf")
pdf(pdf_filename, width = 4, height = 3)
p2
dev.off()

