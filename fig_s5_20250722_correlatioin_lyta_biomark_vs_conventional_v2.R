# Author:Youvika Singh

# Script details:
# script merges database lyta conventional and biomark results 
# and comapres Lyta conventional concentration and biomark 
# concentration

# Last edit: Youvika Singh
# edited for ct value replacement for unidentified (999.000) to 40
# NEW CHIP 8  and subsequent chips, and updated database


rm(list = ls())


# Library
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

#chip <- read.delim("analysis/R/output/all_biomark_20250402.tsv")

#colnames(chip)[1] <- "SampleID_Day"

setwd("Z:/Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/0.20250128_main_Biomark_SAMSaliva_Allchips_YS/")


# Select appropriate file to work

chip <- read.delim("analysis/R/output/all_biomark_20250722.tsv")



setwd("Z:\\Projects\\SAMSAM\\Manuscripts\\Figures\\Supplemental figure 5_specificity")



#save the session info
sink("sessionInfo.txt")
sessionInfo()
sink()




#chip <- chip %>% filter(!grepl("Chip8", chip))
##############################################################
#
#Figure S5A
#
################################################################
spn <- chip %>% filter(grepl("^Spn", Assay))

spn <- spn %>%
  separate(Sample, into = c("Sample", "Day"), sep = "_")

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
  theme_minimal() +
  labs(
    title = "Serotypes specificity",
    x = "Assay",
    y = "Number of times observed in absence of lytA/piaB"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Remove the legend if it's unnecessary
  ) +  
  coord_flip()


print(p)

timestamp <- format(Sys.time(), "%Y%m%d")
output_dir <- "Z:/Projects/SAMSAM/Manuscripts\\Figures\\Supplemental figure 5_specificity/"

pdf_filename <- paste0(output_dir,"S5A_bar_specificity_spn_serotypes_", timestamp, ".pdf")

pdf(pdf_filename, width = 15, height = 12)
p
dev.off()



##############################################################
#
#Figure S5B
#
################################################################



conven <- read.delim("Z:/Projects/SAMSAM/MetaData/Files/Raw_files/exported_averages_LytA_all.txt")

db <- readxl::read_excel("Z:/Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDatabaseFinalVersionBiomark_Flu24_Correct_20250722.xlsx")

#mrg <- merge(db, conven, by = "SampleID_Day")

mrg2 <- merge(db, conven, 
              by.x = "SampleID_Day", 
              by.y = "SampleID")

#colnames(mrg2)

# Rows in db but not in chip
#db_not_in_chip <- anti_join(db, chip, by = "SampleID_Day")

# Rows in chip but not in db
#chip_not_in_db <- anti_join(chip, db, by = "SampleID_Day")

mycols <- c("Study",
            "StudyID", 
            "Day", 
            "SampleID_Day",
            "Type" ,
            "Assay",
            "Ct_Value_1",
            "Ct_Value_2",
            "Conc_1",
            "Conc_2",
            "Average_conc",
            "Flag",
            "Positive",
            "CV_percenatge",
            "LytA_Copies_average",
            "chip")

mrg_filt <- mrg2 %>% dplyr::select(all_of(mycols)) 
mrg_filt$log10_avgconc_biomark <- log10(mrg_filt$Average_conc+1)
mrg_filt$log10_avgconc_conventional <- log10(mrg_filt$LytA_Copies_average+1)

# mrg_filt <- mrg_filt %>%
#   group_by(StudyID, Assay) %>%
#   filter(
#     # Keep groups where at least one value is valid
#     any(
#       !(grepl("^999", as.character(Ct_Value_1)) | Ct_Value_1 == 40),
#       !(grepl("^999", as.character(Ct_Value_2)) | Ct_Value_2 == 40),
#       !(grepl("^999", as.character(Conc_1)) | Conc_1 == 0),
#       !(grepl("^999", as.character(Conc_2)) | Conc_2 == 0)
#     )
#   ) %>%
#   ungroup() 

mrg_filt$Day2 <- as.numeric(gsub('D', '', mrg_filt$Day)) #fix teh column as the order was not always right....
mrg_filt$Day <- factor(mrg_filt$Day, levels = unique(mrg_filt$Day)) 

mrg_filt <- mrg_filt%>% filter(Assay == "Spn(lytA)")

mrg_filt_long <- mrg_filt %>% 
  pivot_longer(cols = c("log10_avgconc_biomark", 
                        "log10_avgconc_conventional"),
               names_to = "variable", 
               values_to = "value")

mrg_filt_long <- mrg_filt_long %>% 
  dplyr::select(StudyID, Day2, variable, value)

mrg_filt_long <- mrg_filt_long %>%
  dplyr::mutate(
    type = dplyr::case_when(
      grepl("conven", variable) ~ "conventional",
      grepl("biomark", variable) ~ "Biomark",
      TRUE ~ "unknown" # Optional: catch-all for unexpected cases
    )
  )


timestamp <- format(Sys.time(), "%Y%m%d")



db <- db %>%
  mutate(across(starts_with("Ct"), ~ ifelse(. == 40, 30, .)))
db$Ct_Average_biomark <- (db$Ct_Value_1 + db$Ct_Value_2)/2


cor.test(db_lytpia$lytA_Ct, db_lytpia$piaB_Ct)

db_lytpia$Concordant  <- db_lytpia$StudyID
db_lytpia$Concordant[!db_lytpia$StudyID %in% c('SAM-K1-25', 'SAM-K1-26',
                                               'SAM-K1-30','SAM-K1-34',
                                               'SAM-K1-37', 'SAM-K1-38')] <-  NA

p1 <- ggplot(db_lytpia, aes(x = lytA_Ct, 
                      y = piaB_Ct, 
                      colour = Concordant)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  ggtitle("Lyta Ct average biomark vs conventional") 


pdf('correlation_lytApiaB.pdf', width =5, height = 5)
print(p1)
dev.off()



#---------------------------------
# Correlation plots for Ct values
#--------------------------------

mycols_forct <- c("Study",
                  "StudyID",
                  "Day",
                  "SampleID_Day",
                  "Type" ,
                  "Assay",
                  "Ct_Value_1",
                  "Ct_Value_2",
                  "Conc_1",
                  "Conc_2",
                  "Average_conc",
                  "Flag",
                  "Positive",
                  "CV_percenatge",
                  "LytA_Ct_average",
                  "LytA_Copies_average",
                  "chip")

mrg_filt2 <- mrg2 %>% 
  dplyr::select(all_of(mycols_forct)) 

mrg_filt2 <- mrg_filt2 %>%
  group_by(StudyID, Assay) %>%
  filter(
    # Keep groups where at least one value is valid
    any(
      !(grepl("^999", as.character(Ct_Value_1)) | Ct_Value_1 == 40),
      !(grepl("^999", as.character(Ct_Value_2)) | Ct_Value_2 == 40),
      !(grepl("^999", as.character(Conc_1)) | Conc_1 == 0),
      !(grepl("^999", as.character(Conc_2)) | Conc_2 == 0)
    )
  ) %>%
  ungroup() 

mrg_filt2$Day2 <- as.numeric(gsub('D', '', mrg_filt2$Day)) #fix teh column as the order was not always right....
mrg_filt2$Day <- factor(mrg_filt2$Day, 
                        levels = unique(mrg_filt2$Day)) 

# Replace 999 in columns starting with "Ct" with 0

# mrg_filt2 <- mrg_filt2 %>%
#   mutate(across(starts_with("Ct"), ~ ifelse(. == 999, 0, .)))
# 

mrg_filt2 <- mrg_filt2 %>%
  mutate(across(starts_with("Ct"), ~ ifelse(. == 40, 30, .)))


mrg_filt2 <- mrg_filt2%>% filter(Assay == "Spn(lytA)")

#mrg_filt2 <- mrg_filt2 %>%
 # mutate(Ct_Average_biomark = rowMeans(select(., starts_with("Ct_")), na.rm = TRUE))

mrg_filt2$Ct_Average_biomark <- (mrg_filt2$Ct_Value_1 + mrg_filt2$Ct_Value_2)/2


mrg_filt2$log10_avgconc_biomark <- log10(mrg_filt2$Average_conc+1)
mrg_filt2$log10_avgconc_conventional <- log10(mrg_filt2$LytA_Copies_average+1)


mrg_filt2 <- mrg_filt2 %>%
  mutate(Ct_Average_biomark = ifelse(Ct_Average_biomark == 0.000000,
                                     30.00000, 
                                     Ct_Average_biomark))

mrg_filt2$outlier <- 'No'
mrg_filt2$outlier[(mrg_filt2$LytA_Ct_average - (mrg_filt2$Ct_Average_biomark+14)>4 |
                 mrg_filt2$LytA_Ct_average - (mrg_filt2$Ct_Average_biomark+14)< -5)& 
                    mrg_filt2$Ct_Average_biomark<30] <- 'Yes'    
mrg_filt2$Name <- mrg_filt2$SampleID_Day
mrg_filt2$Name[mrg_filt2$outlier == 'No'] <- NA

pdf_filename <- paste0(output_dir,
                       "S5B_line_lyta_ct_biomark_vs_conventional_",
                       timestamp,
                       ".pdf")

pdf(pdf_filename, width =5, height = 5)

p1 <- ggplot(mrg_filt2, aes(x = LytA_Ct_average, 
                            y = Ct_Average_biomark, 
                            colour = chip)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  ggtitle("Lyta Ct average biomark vs conventional") 

print(p1)
dev.off()

cor.test(mrg_filt2$LytA_Ct_average, mrg_filt2$Ct_Average_biomark)


##############################################################
#
#Table S5C
#
################################################################
table(mrg_filt2$LytA_Ct_average < 40, mrg_filt2$Ct_Average_biomark < 30)
table(mrg_filt2$LytA_Ct_average < 40)
19/(19+1067)
table(mrg_filt2$LytA_Ct_average < 40, mrg_filt2$Ct_Average_biomark < 30) / 1231





##############################################################
#
#Figure S5D
#
################################################################


FN <- mrg_filt2[mrg_filt2$LytA_Ct_average < 40 & mrg_filt2$Ct_Average_biomark == 30,]$SampleID_Day

FNs <- db[db$SampleID_Day %in% FN & grepl('Spn', db$Assay), ]


heatmap_data <- FNs %>%
  select(SampleID_Day, Assay, log10_avgconc) %>%
  pivot_wider(names_from = SampleID_Day, 
              values_from = log10_avgconc,
              values_fill = 0)
heatmap_data <- data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data[,1]
heatmap_data2 <- heatmap_data[,-1]
heatmap_data3 <- heatmap_data2[rowSums(heatmap_data2)>0,]


breaksList = seq(0, 6, by = 0.2)


FP <- mrg_filt2[mrg_filt2$LytA_Ct_average == 40 & mrg_filt2$Ct_Average_biomark < 30,]$SampleID_Day
FPs <- db[db$SampleID_Day %in% FP & grepl('Spn', db$Assay), ]
heatmap_FP <- FPs %>%
  select(SampleID_Day, Assay, log10_avgconc) %>%
  pivot_wider(names_from = SampleID_Day, 
              values_from = log10_avgconc,
              values_fill = 0)
heatmap_FP <- data.frame(heatmap_FP)
rownames(heatmap_FP) <- heatmap_FP[,1]
heatmap_FP2 <- heatmap_FP[,-1]
heatmap_FP3 <- heatmap_FP2[rowSums(heatmap_FP2)>0,]



table(colSums(heatmap_data3)>0)


pdf_filename <- paste0(output_dir,
                       "S5D_heatmap_FP_biomark_",
                       timestamp,
                       ".pdf")

pdf(pdf_filename, width =5, height = 5)
pheatmap(heatmap_FP3, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()




##############################################################
#
#Figure S5E
#
################################################################


pdf_filename <- paste0(output_dir,
                       "S5E_heatmap_FN_biomark_S5D",
                       timestamp,
                       ".pdf")

pdf(pdf_filename, width =10, height = 5)
pheatmap(heatmap_data3, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()








##############################################################
#
#Figure S5F
#
################################################################
db_pia <- db[db$Assay == 'Spn(piaB)',]
db_lyt <- db[db$Assay == 'Spn(lytA)',]

identical(db_pia$SampleID_Day, db_lyt$SampleID_Day)
colnames(db_pia)

db_lytpia <- cbind(db_lyt[,c(1:7, 78,81,80)], db_pia[,c(81,80)])
colnames(db_lytpia)[9:12] <- c('lytA_Ct', 'lytA_Conc', 'piaB_Ct', 'piaB_Conc')



cor.test(db_lytpia$lytA_Ct, db_lytpia$piaB_Ct)

db_lytpia$Concordant  <- db_lytpia$StudyID
db_lytpia$Concordant[!db_lytpia$StudyID %in% c('SAM-K1-25', 'SAM-K1-26',
                                               'SAM-K1-30','SAM-K1-34',
                                               'SAM-K1-37', 'SAM-K1-38')] <-  NA

p1 <- ggplot(db_lytpia, aes(x = lytA_Ct, 
                            y = piaB_Ct, 
                            colour = Concordant)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  ggtitle("Lyta Ct average biomark vs conventional") 


pdf('S5F_correlation_lytApiaB.pdf', width =5, height = 5)
print(p1)
dev.off()












save.image('Figure_s5_specificity.Rdata')



