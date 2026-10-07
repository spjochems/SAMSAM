#####
#This makes plots of exposure
#update 9-6. Fixed the forest plot for exposure vs activities to be of standard error and not CI

rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\Fig2\\")

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
library(ggmosaic)
library(reshape2)
library(circlize)
library(ComplexHeatmap)
library(broom.mixed)
library(parallel)
library(glmmTMB)

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

#enocde the season as a cosinus
mrg$DateCollection_EditDate #check what the date looks like
mrg$DaysSinceJan01 <- as.integer(format(as.Date(mrg$DateCollection_EditDate, format = "%Y_%m_%d"), "%j")) - 1 #set it to an integer
mrg$DaysCosinus <-  cos(2*pi*mrg$DaysSinceJan01/365) #change it to a cosinus


mrg2 <- mrg[!mrg$Assay %in% c('16S', 'PIV1','PIV2','Spn20',
                              'Moraxella_catarrhalis', 'Haemophilus_influenzae',
                              'Spn(lytA)', 'Spn(piaB)'),] #remove 16S and the none or always detected assays

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
#                 FIGURE 2A                       #
#                                                 #
###################################################

mrg3 <- mrg %>%
  mutate(
    Category = case_when(
      str_detect(Assay, "16S|Mycoplasma|Staphylococcus|Moraxella|Streptococcus|Haemophilus") ~ "Bacteria",
      str_detect(Assay, "Spn") ~ "Bacteria_spn",
      str_detect(Assay, "sars-cov-2|coronavirus") ~ "Corona_virus",
      TRUE ~ "Other_virus"
    )
  )%>%
  arrange(SampleID_Day, Assay) %>%  
  distinct(SampleID_Day, Assay, .keep_all = TRUE)

mrgBackup <- mrg
removed_rows <- anti_join(mrgBackup, mrg3, by = colnames(mrgBackup))

allkids <- "SAM-K1-33"

for (everykid in allkids) {
  
  kid <- mrg3 %>% filter(StudyID == everykid)
  
  # Convert 'Day' to a character to avoid factor issues
  kid$Day <- as.character(kid$Day)
  
  # Remove "D" and convert to numeric for sorting
  kid$Day_num <- as.numeric(sub("D", "", kid$Day))
  
  # Order by the numeric value of 'Day_num'
  kid <- kid[order(kid$Day_num), ]
  
  # Reassign 'Day' as a factor with ordered levels based on sorted 'Day_num'
  kid$Day <- factor(kid$Day, levels = unique(kid$Day[order(kid$Day_num)]))
  
  # Categorize log10_avgconc
  kid <- kid %>%
    mutate(
      Detection = case_when(
        log10_avgconc == 0 ~ "Undetected",
        log10_avgconc < 2 ~ "<10^2",
        log10_avgconc >= 2 & log10_avgconc < 4 ~ "10^2-10^4",
        log10_avgconc >= 4 ~ ">10^4"
      )
    )
  
  # Convert Detection categories to numeric for heatmap
  detection_levels <- c("Undetected" = 0, 
                        "<10^2" = 1, 
                        "10^2-10^4" = 2, 
                        ">10^4" = 3)
  
  kid <- kid %>%
    mutate(DetectionNumeric = detection_levels[Detection])
  
  
  # annotation data
  col_annotation <- kid %>%
    dplyr::select(Day, 
                  #SnottyNose,
                  #SnottyNose_Level,
                  #Sneeze,
                  #Cough,
                  #Throat, 
                  #Fever, 
                  #Ear,
                  #OtherComplaints,
                  #NoActivity,
                  SchoolDaycare,
                  #Swimmingpool,
                  OtherActivity)%>%
                  #Antibiotics
                  distinct() %>%
    column_to_rownames("Day")
  
  col_annotation[] <- lapply(col_annotation, function(x) ifelse(is.na(x), "NA", x))
  #col_annotation <- col_annotation %>%
   # mutate(across(!SnottyNose_Level, ~ ifelse(. == "NA", "No", .)))
  
  
  # remove annotation rows with No
  
  col_annotation <- col_annotation %>%
    dplyr::select(dplyr::where(~ any(. != "No")))
  
  # # Pivot data for heatmap
  #   heatmap_data <- kid %>%
  #   select(Assay, Day, DetectionNumeric) %>%
  #   pivot_wider(names_from = Day, values_from = DetectionNumeric) %>%
  #   column_to_rownames(var = "Assay") %>%
  #   as.matrix()
  
  heatmap_data <- kid %>%
    dplyr::select(Day, Assay, log10_avgconc) %>%
    pivot_wider(names_from = Day, 
                values_from = log10_avgconc,
                values_fill = 0)
  
  # remove EAV, 19b/f and 18a/b/c
  heatmap_data <- heatmap_data[heatmap_data$Assay != "EAV" & 
                                 heatmap_data$Assay != "Spn19B/F" &
                                 heatmap_data$Assay != "Spn18A/B/C", ]
  
  # Prepare the heatmap matrix
  heatmap_matrix <- as.matrix(heatmap_data[,-1])
  rownames(heatmap_matrix) <- heatmap_data$Assay
  heatmap_matrix <- heatmap_matrix[rowSums(heatmap_matrix != 0) > 0, ]
  
  # Define color scale
  breaks <- c(-1, 0, 2, 4, max(heatmap_matrix, na.rm = TRUE))
  colors <- c("grey", "yellow", "orange", "red")
  
  
  # Filter `kid` dataframe to match the rows remaining in `heatmap_matrix`
  row_split <- kid %>%
    filter(Assay %in% rownames(heatmap_matrix)) %>%
    distinct(Assay, Category)
  
  # Define row splits based on the filtered `kid` dataframe
  row_split <- row_split %>%
    arrange(match(rownames(heatmap_matrix), Assay)) %>%  # Ensure the split order matches heatmap rows
    pull(Category)
  
  
  # # Prepare annotations for columns
  # col_annotation <- col_annotation %>%
  #   lapply(function(x) ifelse(x == "NA", "No", x)) %>%
  #   as.data.frame()
  # 
  # # Filter annotation columns with at least one "Yes"
  # col_annotation <- col_annotation %>%
  #   select(where(~ any(. != "No")))
  # 
  # Column annotation colors
  annotation_colors <- list(
    #SnottyNose = c("Yes" = "burlywood3", "No" = "white"),
    #SnottyNose_Level = c("0" = "white", "1" = "pink", "2"= "palevioletred","NA" = "grey90"),
    #Sneeze = c("Yes" = "burlywood3", "No" = "white"),
    #Cough = c("Yes" = "burlywood3", "No" = "white"),
    #Throat = c("Yes" = "burlywood3", "No" = "white"),
    #Fever = c("Yes" = "burlywood3", "No" = "white"),
    #Ear = c("Yes" = "burlywood3", "No" = "white"),
    #OtherComplaints = c("Yes" = "burlywood3", "No" = "white"),
    #NoActivity = c("Yes" = "skyblue", "No" = "white"),
    SchoolDaycare = c("Yes" = "skyblue", "No" = "white", "NA" = "grey90"),
    #Swimmingpool = c("Yes" = "skyblue", "No" = "white"),
    OtherActivity = c("Yes" = "skyblue", "No" = "white", "NA" = "grey90")
    #Antibiotics = c("Yes" = "purple", "No" = "white")  # Corrected closing
  )
  
    # Create column annotations for the heatmap
  column_annot <- HeatmapAnnotation(
    df = col_annotation,
    col = annotation_colors,
    show_legend = FALSE, # for removing the legend as plot is noisy
    annotation_legend_param = list(
      title = "Symptoms/Activities",
      legend_direction = "vertical"
    ),
    annotation_name_side = "left",  # Adjust annotation name side
    annotation_name_gp = gpar(fontsize = 10),  # Optional: control annotation name font size
    gap = unit(1, "mm"),  # Adding a gap between columns of the annotation (this adds the white line)
    border = TRUE  # Adds borders around the individual annotation columns
  )
  
  # prepare column for higher ct value (>2) difference and discordant duplicates
  # code when influenza cutoff if NOT used
  
  #kid$ct_value_diff <- abs(kid$Ct_Value_1 - kid$Ct_Value_2) 
  
  kid$finalcol_hp <- ifelse(kid$ExposureEvent == 1, "#", "")
  
  # prepare table for text
  text_matrix <- kid %>%
    dplyr::select(Day, Assay, finalcol_hp) %>%
    pivot_wider(names_from = Day, values_from = finalcol_hp, 
                values_fill = "") %>%
    column_to_rownames(var = "Assay") %>%
    as.matrix()
  
  # Ensure the text matrix matches the heatmap_matrix dimensions
  text_matrix <- text_matrix[rownames(heatmap_matrix), 
                             colnames(heatmap_matrix)]
  
  p <- Heatmap(
    heatmap_matrix,
    #name = "Microbial Detection",
    col = colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = column_annot,
    row_split = row_split,
    column_title = paste0(everykid, sep = "_", "Days"),
    column_title_side = "bottom",
    row_title = "Assays",
    rect_gp = gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(
      title = "log10(avg conc)",
      at = c(0, 2, 4, max(heatmap_matrix, na.rm = TRUE)),
      labels = c("Undetected", "<2", "2-4", ">4")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        text_matrix[i, j],  # Add corresponding value from text_matrix
        x = x, y = y,
        gp = gpar(fontsize = 12, col = "black")  # Customize font size and color
      )
    }
  )
  draw(p)
}
p 

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig2A_heatmap_SAM033_", timestamp, ".pdf")
pdf(pdf_filename, width = 15,height = 10)
p
dev.off()


###################################################
#                                                 #
#                 FIGURE 2B                       #
#                                                 #
###################################################

#get association with school going. Select only the timepoints of exposure or not
#only select the ones where we have a potential exposure,
mrg4 <- mrg2[mrg2$Repeats < 2 & #remove any with ongoing infection
               mrg2$Day2 >4 & #remove any in first 4 days as there is no exposiure possible according to definition
               !(mrg2$Repeats ==0 & mrg2$InfectionLonger == 1),]#remove the ones where there is an ongoing infection but intermittant negative 

mrg4$Assay2 <- mrg4$Assay
mrg4$Assay2[mrg4$AssayType == 'Viral'] <- 'AnyVirus' 
mrg4$Assay2[grepl('Spn',mrg4$Assay) & 
              !grepl('lytA',mrg4$Assay) & 
              !grepl('piaB',mrg4$Assay) & 
              !grepl('Spn6C/D',mrg4$Assay)] <- 'AnySerotype'

mrg5a <- mrg4 %>% 
  filter(!grepl("NA", Swimmingpool)) %>% 
  filter(grepl("Evening", PeriodCollection) | grepl("Afternoon", PeriodCollection)  )

m <- glmer(as.factor(ExposureEvent) ~ SchoolDaycare + OtherActivity + Swimmingpool + DaysCosinus + (1|Assay) + (1|StudyID), data = mrg5a, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)


summary(m)[[10]][,1:2]
OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Activity <- gsub('Yes', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)

OR_adjusted$Pval <- paste0('p=', sapply(OR_adjusted$`Pr(>|z|)`, FUN = "disp_pval"))
names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(`Pr...z..` < 0.001, "***",
                                 ifelse(`Pr...z..` < 0.01,  "**",
                                        ifelse(`Pr...z..` < 0.05,  "*", ""))))
write_tsv(OR_adjusted, "Fig2B_statssummary.txt")


ggplot(OR_adjusted, aes(x=OR, y=Activity)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  geom_vline(xintercept = 1, linetype = 'dashed') + 
  xlab('Adjusted Odds Ratio LMM') 

ggsave('Fig2B_AdjustedOddsRatioExposure.pdf', width = 3, height = 2)


#Count activities events:
activities <- c("Swimmingpool", "SchoolDaycare", "OtherActivity")

activity_summary <- lapply(activities, function(col) {
  mrg4 %>%
    filter(grepl("Yes", .data[[col]])) %>%
    summarise(
      Activity = col,
      total_unique_studyids = n_distinct(StudyID),
      total_unique_samples = n_distinct(SampleID_Day)
    )
}) %>%
  bind_rows()

activity_summary




###################################################
#                                                 #
#                 FIGURE 2C                       #
#                                                 #
###################################################

mrg_exposure <- mrg2[mrg2$ExposureEvent == 1,]
exposed_paths <- data.frame(table(table(mrg_exposure$SampleID_Day)))
sum(as.numeric(exposed_paths[,1]) * as.numeric(exposed_paths[,2])) #total number of exposires
sum(as.numeric(exposed_paths[,1]) * as.numeric(exposed_paths[,2])) / 45 #average number of exposires/child

exposed_paths[1,2] / sum(as.numeric(exposed_paths[,1]) * as.numeric(exposed_paths[,2]))

Fig2C <- ggplot(exposed_paths, aes(y = Var1, x = Freq)) + 
  geom_bar(stat = "identity", fill = "steelblue", colour = "black") + 
  scale_x_continuous(expand = c(0, 0)) +
  geom_text(
    aes(
      label = Freq,
      x = ifelse(Freq < 20, Freq + 5, Freq - 5),
      hjust = ifelse(Freq < 20, 0, 1),
      colour = ifelse(Freq < 20, "small", "large")
    )) +
  scale_colour_manual(values = c("small" = "black", "large" = "white"), guide = "none") +
  ylab("Number of co-detected exposure events") +
  xlab("Number of exposure events") +
  coord_cartesian(clip = "off")
Fig2C

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig2C_", timestamp, ".pdf")
ggsave(pdf_filename, width = 3, height = 2)

sum(as.numeric(exposed_paths$Var1) * exposed_paths$Freq ) #check the total number of exposures


#add multiplication:
plot_df <- exposed_paths %>%
  mutate(Freq2 = Freq * as.numeric(as.character(Var1))) 

plot_text <- plot_df %>%
  pivot_longer(c(Freq, Freq2), names_to = "Type", values_to = "Value") %>%
  mutate(
    label_x = ifelse(Value < 66, Value + 5, Value - 5),
    label_hjust = ifelse(Value < 66, 0, 1),
    label_col = ifelse(Value < 66, "small", "large")
  )

Fig2C <- ggplot(plot_df, aes(y = Var1)) + 
  geom_col(aes(x = Freq2), fill = "lightblue", colour = "black", position = "identity") +
  geom_col(aes(x = Freq), fill = "steelblue", colour = "black", position = "identity") +
  geom_text(
    data = plot_text,
    aes(x = label_x,
      y = Var1,
      label = Value,
      hjust = label_hjust,
      colour = label_col)) +
  scale_colour_manual(values = c("small" = "black", "large" = "white"), guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  ylab("Number of co-detected exposure events") +
  xlab("Number of exposure events") +
  coord_cartesian(clip = "off")
Fig2C

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig2C_", timestamp, ".pdf")
ggsave(pdf_filename, width = 3, height = 2)



###################################################
#                                                 #
#                 FIGURE 2D                       #
#                                                 #
###################################################

mrg.split <- split(mrg2, f=mrg2$Assay)
for(i in 1:length(mrg.split)){
  pres <- table(mrg.split[[i]]$ExposureEvent, mrg.split[[i]]$StudyID)
  acq <- table(mrg.split[[i]]$AcquisitionLonger, mrg.split[[i]]$StudyID)
  abort <- table(mrg.split[[i]]$AbortiveExposureLonger, mrg.split[[i]]$StudyID)
  
  #add the presence table
  if(nrow(pres)>1){
    if(i == 1){ 
      presence <- pres[2,]
    }else{
      presence <- cbind(presence, pres[2,])  
      colnames(presence)[ncol(presence)] <- names(mrg.split)[i]
    }
  }
  
  #add the acquisition table
  if(nrow(acq)>1){
    if(i == 1){ 
      acquisition <- acq[2,]
    }else{
      acquisition <- cbind(acquisition, acq[2,])  
      colnames(acquisition)[ncol(acquisition)] <- names(mrg.split)[i]
    }
  }
  
  
  #add the abortive table
  if(nrow(abort)>1){
    if(i == 1){ 
      abortive <- abort[2,]
    }else{
      abortive <- cbind(abortive, abort[2,])  
      colnames(abortive)[ncol(abortive)] <- names(mrg.split)[i]
    }
  }
  
}
colnames(abortive)[1] <- colnames(acquisition)[1] <- colnames(presence)[1] <- names(mrg.split)[1]

length(which(grepl('Spn', colnames(presence)) & !grepl('Spn6C/D', colnames(presence))))
colnames(presence)[which(colnames(presence) %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay))]

exposure_perpathogen <- data.frame(colSums(presence))
abortive_perpathogen <- data.frame(colSums(abortive))
acquisition_perpathogen <- data.frame(colSums(acquisition))
expab <- merge(exposure_perpathogen, abortive_perpathogen, by= 'row.names', all.x=T)
all_perpath <- merge(expab, acquisition_perpathogen, by.x= 'Row.names', by.y= 'row.names', all.x=T)
colnames(all_perpath) <- c('Assay', 'Exposures', 'Abortive', 'Acquisition')
all_perpath[is.na(all_perpath)] <- 0

#calculate the number of observations per spn or total virus
spna <- sum(all_perpath[grepl('Spn', all_perpath$Assay) & !grepl('Spn6C/D', all_perpath$Assay ),2])
spnb <- sum(all_perpath[grepl('Spn', all_perpath$Assay) & !grepl('Spn6C/D', all_perpath$Assay ),3])
spnc <- sum(all_perpath[grepl('Spn', all_perpath$Assay) & !grepl('Spn6C/D', all_perpath$Assay ),4])

vira <- sum(all_perpath[all_perpath$Assay %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay),2])
virb <- sum(all_perpath[all_perpath$Assay %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay),3])
virc <- sum(all_perpath[all_perpath$Assay %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay),4])

all_perpath2 <- rbind(all_perpath, 
                      c('AnySpnSerotype', spna, spnb, spnc),
                      c('AnyVirus', vira, virb, virc))
all_perpath2[,2] <- as.numeric(all_perpath2[,2] )
all_perpath2[,3] <- as.numeric(all_perpath2[,3] )
all_perpath2[,4] <- as.numeric(all_perpath2[,4] )

all_perpath2$ExposureRate <- all_perpath2[,2] / 45
all_perpath2$AcquisitionPercentage <-  all_perpath2[,4] / all_perpath2[,2] * 100
all_perpath2$AcquisitionPercentage2 <- round(all_perpath2$AcquisitionPercentage, 0) 
all_perpath2$AcquisitionPercentage2[all_perpath2$Exposures <20] <- NA

all_perpath2$Assay <- factor(all_perpath2$Assay, 
                             levels = all_perpath2$Assay[order(all_perpath2$ExposureRate)]) 

pl1 <- ggplot(all_perpath2, aes(y = Assay, x = ExposureRate)) + 
  geom_bar(stat = "identity", fill = "steelblue", colour = "black") +
  xlab("Exposure Rate per child follow up") +
  theme(
    plot.margin = margin(5.5, -10, 5.5, 5.5)
  )

pl2 <- ggplot(all_perpath2, aes(x = 1, y = Assay, fill = AcquisitionPercentage2)) +
  geom_tile(color = "transparent") +
  geom_text(aes(label = AcquisitionPercentage2), color = "white", size = 3.5) +
  scale_fill_gradient(
    na.value = "transparent",
    name = "Aquisition\npercentage"
  ) +
  coord_fixed() +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5.5, -25, 5.5, -35),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, -25),
    legend.margin = margin(0, 0, 0, 0)
  ) + theme(aspect.ratio = 13)
pl2

leg <- get_legend(
  pl2 +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "right",
      legend.justification = c(0, 0.5),
      legend.box.margin = margin(0, 0, 0, -35),
      legend.margin = margin(0, 0, 0, 0)
    )
)
pl2_noleg <- pl2 + theme(legend.position = "none")

plot_grid(
  pl1,
  pl2_noleg,
  leg,
  nrow = 1,
  rel_widths = c(1, 0.10, 0.06)
)

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig2D_", timestamp, ".pdf")
ggsave(pdf_filename, width = 5, height = 6)



###################################################
#                                                 #
#                 FIGURE 2E                       #
#                                                 #
###################################################

mrg4_exposure_virus <- mrg4[mrg4$ExposureEvent == 1 & mrg4$Assay2 == 'AnyVirus',]
exposure_table_virus <- table(mrg4_exposure_virus$AgeYears, mrg4_exposure_virus$Assay)
mrg4_exposure_virus2 <- cbind( exposure_table_virus, table(mrg[which(!duplicated(mrg$StudyID)),]$AgeYears))
for(i in 1:ncol(mrg4_exposure_virus2)){mrg4_exposure_virus2[,i] <- mrg4_exposure_virus2[,i]/mrg4_exposure_virus2[,ncol(mrg4_exposure_virus2)]}
mrg4_exposure_virus3 <- t(mrg4_exposure_virus2[,-13])
mrg4_exposure_virus3 <- rbind(mrg4_exposure_virus3, colSums(mrg4_exposure_virus3 ))
rownames(mrg4_exposure_virus3)[nrow(mrg4_exposure_virus3)] <- 'AnyVirus'

#stats 
assays <- rownames(mrg4_exposure_virus3)

res <- lapply(assays, function(a) {
  dat <- subset(mrg4, Assay == a)
  n_events <- sum(dat$ExposureEvent == 1, na.rm = TRUE)
  
  if (n_events < 5) {
    return(data.frame(
      Assay = a,
      n_events = n_events,
      p_age = NA_real_,
      stars = ""
    ))
  }
  
  fit <- tryCatch(
    glmer(as.factor(ExposureEvent) ~ AgeYears + (1 | StudyID),
          data = dat,
          family = binomial,
          control = glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      Assay = a,
      n_events = n_events,
      p_age = NA_real_,
      stars = ""
    ))
  }
  
  coefs <- coef(summary(fit))
  if (!"AgeYears" %in% rownames(coefs)) {
    return(data.frame(
      Assay = a,
      n_events = n_events,
      p_age = NA_real_,
      stars = ""
    ))
  }
  
  p <- coefs["AgeYears", "Pr(>|z|)"]
  s <- ifelse(p < 0.001, "***",
              ifelse(p < 0.01, "**",
                     ifelse(p < 0.05, "*", "")))
  
  data.frame(
    Assay = a,
    n_events = n_events,
    p_age = p,
    stars = s
  )
})

res_df <- do.call(rbind, res)
res_df

#check also Anyvirus
mrg4_select <- mrg4[mrg4$Assay2 %in% c('AnyVirus'),]
m_any <- glmer(as.factor(ExposureEvent) ~ AgeYears + (1|Assay)  + (1|StudyID), data = mrg4_select, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m_any) #check if exposure is realted to age

table2 <- tidy(m_any) %>%
  filter(effect == "fixed") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value)
table2

write_tsv(table2, "Fig2E_statssummary.txt")


coefs_any <- coef(summary(m_any))

p_any <- if ("AgeYears" %in% rownames(coefs_any)) coefs_any["AgeYears", "Pr(>|z|)"] else NA_real_

s_any <- if (is.na(p_any)) "" else
  if (p_any < 0.001) "***" else if (p_any < 0.01) "**" else if (p_any < 0.05) "*" else ""

any_row <- data.frame(
  Assay = "AnyVirus",
  n_events = sum(mrg4_select$ExposureEvent == 1, na.rm = TRUE),
  p_age = p_any,
  stars = s_any
)

res_df <- rbind(
  res_df,
  any_row
)
res_df <- res_df %>%
  filter(!(Assay == "AnyVirus" & n_events == 0))

write_tsv(res_df, "Fig2E_statssummary_PerVirus.txt")

star_map <- setNames(res_df$stars, res_df$Assay)

row_lab <- rownames(mrg4_exposure_virus3)
row_lab <- ifelse(
  row_lab %in% names(star_map) & star_map[row_lab] != "",
  paste0(row_lab, " ", star_map[row_lab]),
  row_lab
)
res_df

# AnyVirus separated, then all other rows in alphabetical order
rn <- rownames(mrg4_exposure_virus3)
new_order <- c("AnyVirus", setdiff(rn, "AnyVirus"))

mat2 <- mrg4_exposure_virus3[new_order, , drop = FALSE]
row_lab2 <- row_lab[new_order]

#plot:
Fig2E <- ComplexHeatmap::pheatmap(
  mat2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = 1,
  treeheight_row = 0,
  angle_col = "0",
  labels_row = row_lab2,
  fontsize_row = 10,
  fontsize_col = 10,
  heatmap_legend_param = list(
    title = "Average\nexposures\nper child"
  )
)
Fig2E

#save:
pdf_filename <- paste0("Fig2E_heatmap_virusexposure_", timestamp, ".pdf")
pdf(pdf_filename, width = 4, height = 4)
Fig2E
dev.off()

#table(mrg4_select$ExposureEvent[mrg4_select$Assay == "coronavirus_NL63"])



###################################################
#                                                 #
#                 FIGURE 2F                       #
#                                                 #
###################################################

PerChild <- data.frame(Spn_exposure = rowSums(presence[,which(grepl('Spn', colnames(presence)) & !grepl('Spn6C/D', colnames(presence)))]),
                       Spn_abortive = rowSums(abortive[,which(grepl('Spn', colnames(abortive)) & !grepl('Spn6C/D', colnames(abortive)))]),
                       Spn_acquisition = rowSums(acquisition[,which(grepl('Spn', colnames(acquisition)) & !grepl('Spn6C/D', colnames(acquisition)))]),
                       Virus_exposure = rowSums(presence[,which(colnames(presence) %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay))]),
                       Virus_abotive = rowSums(abortive[,which(colnames(abortive) %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay))]),
                       Virus_acquisition = rowSums(acquisition[,which(colnames(acquisition) %in% unique(mrg2[mrg2$AssayType == 'Viral',]$Assay))]),
                       HRV_exposure = presence[,'HRV'],
                       HRV_abortive = abortive[,'HRV'],
                       HRV_acquisition = acquisition[,'HRV'],
                       bocavirus_exposure = presence[,'bocavirus'],
                       bocavirus_abortive = abortive[,'bocavirus'],
                       bocavirus_acquisition = acquisition[,'bocavirus'],
                       Staphylococcus_aureus_exposure = presence[,'Staphylococcus_aureus'],
                       Staphylococcus_aureus_abortive = abortive[,'Staphylococcus_aureus'],
                       Staphylococcus_aureus_acquisition = acquisition[,'Staphylococcus_aureus']
)

demo <- mrg[!duplicated(mrg$StudyID),c(2,6,7)]
PerChild2 <- merge(PerChild, demo, by.x = 'row.names', by.y = 'StudyID')
PerChild2$Spn_acquisitionSucces <- PerChild2$Spn_acquisition / PerChild2$Spn_exposure * 100
PerChild2$VirusacquisitionSucces <- PerChild2$Virus_acquisition / PerChild2$Virus_exposure * 100
PerChild2$HRVacquisitionSucces <- PerChild2$HRV_acquisition / PerChild2$HRV_exposure * 100
PerChild2$bocavirusacquisitionSucces <- PerChild2$bocavirus_acquisition / PerChild2$bocavirus_exposure * 100
PerChild2$StaphaureusacquisitionSucces <- PerChild2$Staphylococcus_aureus_acquisition / PerChild2$Staphylococcus_aureus_exposure * 100

Fig2F <- ggplot(PerChild2, aes(x=as.factor(AgeYears), y=Spn_exposure)) + 
  geom_boxplot(outlier.color = NA, aes(fill = AgeYears)) + geom_jitter(width=0.2, height = 0) + 
  scale_fill_gradient2(high = 'orange')
Fig2F

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig2F_Spn_exposure_", timestamp, ".pdf")
ggsave(pdf_filename, width = 5, height = 4)


#nr samples:
fig_data      <- PerChild2
group_cols    <- c("AgeYears")   # 2–4 columns as needed
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig2F_samplesize.txt")



###################################################
#                                                 #
#                 FIGURE 2G                       #
#                                                 #
###################################################

mrg_exposure$Virus <- NA
mrg_exposure$VirusNumber <- NA
mrg_exposure$SpnNumber  <- NA
mrg_exposure$Staph <- NA
mrg_exposure$Haemophilus <- NA
mrg_exposure$LytA <- NA
mrg_exposure$Spnexposures <- NA
mrg_exposure$Totalexposures <- NA

for(i in 1:nrow(mrg_exposure)){
  
  sample_check <- mrg_exposure$SampleID_Day[i]
  sample_check2 <- mrg[mrg$SampleID_Day == sample_check,]
  
  #fill the values
  mrg_exposure$LytA[i] <- sample_check2[sample_check2$Assay == "Spn(lytA)",'log10_avgconc']
  mrg_exposure$Staph[i] <- sample_check2[sample_check2$Assay == "Staphylococcus_aureus",'InfectionLongerTrue']
  mrg_exposure$Haemophilus[i] <- sample_check2[sample_check2$Assay == "Haemophilus_influenzae",'InfectionLongerTrue']
  mrg_exposure$VirusNumber [i] <-sum(sample_check2[sample_check2$AssayType == "Viral", 'InfectionLongerTrue'])
  mrg_exposure$SpnNumber[i] <- sum(sample_check2[grepl('Spn',sample_check2$Assay) & 
                                                   !grepl('lytA',sample_check2$Assay) & 
                                                   !grepl('piaB',sample_check2$Assay) & 
                                                   !grepl('Spn6C/D',sample_check2$Assay), 'InfectionLongerTrue'])
  mrg_exposure$Spnexposures[i] <- sum(sample_check2[grepl('Spn',sample_check2$Assay) & 
                                                      !grepl('lytA',sample_check2$Assay) & 
                                                      !grepl('piaB',sample_check2$Assay) & 
                                                      !grepl('Spn6C/D',sample_check2$Assay), 'ExposureEvent'])
  mrg_exposure$Totalexposures[i] <- sum(sample_check2[ !grepl('piaB',sample_check2$Assay) & 
                                                         !grepl('Spn6C/D',sample_check2$Assay), 'ExposureEvent'])
  
  
}
mrg_exposure$Spn <- 0 
mrg_exposure$Spn[mrg_exposure$SpnNumber >0] <- 1 #make a Spn yes.no for viral comparison
mrg_exposure$Virus <- 0 
mrg_exposure$Virus[mrg_exposure$VirusNumber >0] <- 1 #make a Spn yes.no for viral comparison


#remove from the ones where we see an acquisition for Spn one serotype to prevent it being counted against itself
mrg_exposure_Spn <- mrg_exposure_Spn0 <- mrg_exposure[grepl('Spn',mrg_exposure$Assay) & !grepl('Spn6C/D', mrg_exposure$Assay),]
mrg_exposure_Spn$SpnNumber[mrg_exposure_Spn$AcquisitionLonger == 1] <- mrg_exposure_Spn0$SpnNumber[mrg_exposure_Spn0$AcquisitionLonger == 1] - 1
mrg_exposure_Spn$Spn <- 0 
mrg_exposure_Spn$Spn[mrg_exposure_Spn$SpnNumber >0] <- 1 #make a Spn yes.no for viral comparison
table(mrg_exposure_Spn$AcquisitionLonger)

mrg_exposure_Spn$Spn
mrg_exposure_Spn$Spnexposures2 <- mrg_exposure_Spn$Spnexposures
mrg_exposure_Spn$Spnexposures2[mrg_exposure_Spn$Spnexposures>1] <- '2+'
mrg_exposure_Spn$LytA2 <- 0
mrg_exposure_Spn$LytA2[mrg_exposure_Spn$LytA >4] <- 1


m <- glmer(as.factor(AcquisitionLonger) ~ SpnNumber + AgeYears + DaysCosinus + as.factor(Virus) +  as.factor(Spnexposures2) +(1|Assay)+ (1|StudyID), data = mrg_exposure_Spn, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)


summary(m)[[10]][,1:2]
OR_adjusted <- data.frame(summary(m)[[10]][-1,])
OR_adjusted$Variable <- gsub('as.factor', '', rownames(OR_adjusted))
OR_adjusted$OR <- exp(OR_adjusted$Estimate)

#95% Wald CIs:
OR_adjusted$Lower <- exp(OR_adjusted$Estimate - 1.96 * OR_adjusted$`Std..Error`)
OR_adjusted$Upper <- exp(OR_adjusted$Estimate + 1.96 * OR_adjusted$`Std..Error`)

OR_adjusted$Pval <- paste0('p=', sapply(OR_adjusted$Pr...z.., FUN = "disp_pval"))
names(OR_adjusted)
OR_adjusted$stars <- with(OR_adjusted,
                          ifelse(Pr...z.. < 0.001, "***",
                                 ifelse(Pr...z.. < 0.01,  "**",
                                        ifelse(Pr...z.. < 0.05,  "*", ""))))
write_tsv(OR_adjusted, "Fig2G_statssummary.txt")

ggplot(OR_adjusted, aes(x=OR, y=Variable)) + 
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")  + 
  geom_text(aes(label = stars, x = Upper * 1.08),
            hjust = 0, size = 5, fontface = "bold") +
  geom_vline(xintercept = 1, linetype = 'dashed') + 
  xlab('Adjusted Odds Ratio LMM') 
ggsave('Fig2G_AdjustedOR_Spn_exposure.pdf', width = 3, height = 2)


#sample size:
count_per_group(mrg_exposure_Spn, "SpnNumber")
count_per_group(mrg_exposure_Spn, "AgeYears")
count_per_group(mrg_exposure_Spn, "Virus")
count_per_group(mrg_exposure_Spn, "Spnexposures2")





###################################################
#                                                 #
#                 FIGURE 2H                       #
#                                                 #
###################################################

Fig2H <- ggplot(PerChild2, aes(x = as.factor(AgeYears), y = VirusacquisitionSucces)) + 
  geom_boxplot(outlier.color = NA, aes(fill = AgeYears)) +
  geom_jitter(width = 0.2, height = 0) + 
  scale_fill_gradient2() +
  labs(
    title = "Any virus",
    x = "Age (years)",
    y = "Acquisition (% of exposure)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
Fig2H

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig2H_virus_acquisition_succes_", timestamp, ".pdf")
ggsave(pdf_filename, width = 5, height = 4)


#nr samples:
fig_data      <- PerChild2
group_cols    <- c("AgeYears")   # 2–4 columns as needed
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig2H_samplesize.txt")



###################################################
#                                                 #
#                 FIGURE 2I                       #
#                                                 #
###################################################

mrg_exposure$Virus <- NA
mrg_exposure$VirusNumber <- NA
mrg_exposure$SpnNumber  <- NA
mrg_exposure$Staph <- NA
mrg_exposure$Haemophilus <- NA
mrg_exposure$LytA <- NA

for(i in 1:nrow(mrg_exposure)){
  
  sample_check <- mrg_exposure$SampleID_Day[i]
  sample_check2 <- mrg[mrg$SampleID_Day == sample_check,]
  
  #fill the values
  mrg_exposure$LytA[i] <- sample_check2[sample_check2$Assay == "Spn(lytA)",'log10_avgconc']
  mrg_exposure$Staph[i] <- sample_check2[sample_check2$Assay == "Staphylococcus_aureus",'InfectionLongerTrue']
  mrg_exposure$Haemophilus[i] <- sample_check2[sample_check2$Assay == "Haemophilus_influenzae",'InfectionLongerTrue']
  mrg_exposure$VirusNumber [i] <-sum(sample_check2[sample_check2$AssayType == "Viral", 'InfectionLongerTrue'])
  mrg_exposure$SpnNumber[i] <- sum(sample_check2[grepl('Spn',sample_check2$Assay) & 
                                                   !grepl('lytA',sample_check2$Assay) & 
                                                   !grepl('piaB',sample_check2$Assay) & 
                                                   !grepl('Spn6C/D',sample_check2$Assay), 'InfectionLongerTrue'])
}
mrg_exposure$Spn <- 0 
mrg_exposure$Spn[mrg_exposure$SpnNumber >0] <- 1 #make a Spn yes.no for viral comparison
mrg_exposure$Virus <- 0 
mrg_exposure$Virus[mrg_exposure$VirusNumber >0] <- 1 #make a Spn yes.no for viral comparison

#remove from the ones where we see an acquisition for Spn one serotype to prevent it being counted against itself
mrg_exposure_Spn <- mrg_exposure_Spn0 <- mrg_exposure[grepl('Spn',mrg_exposure$Assay) & !grepl('Spn6C/D', mrg_exposure$Assay),]
mrg_exposure_Spn$SpnNumber[mrg_exposure_Spn$AcquisitionLonger == 1] <- mrg_exposure_Spn0$SpnNumber[mrg_exposure_Spn0$AcquisitionLonger == 1] - 1
mrg_exposure_Spn$Spn <- 0 
mrg_exposure_Spn$Spn[mrg_exposure_Spn$SpnNumber >0] <- 1 #make a Spn yes.no for viral comparison
table(mrg_exposure_Spn$AcquisitionLonger)
#do the same for virus
mrg_exposure_virus <- mrg_exposure_virus0 <- mrg_exposure[mrg_exposure$AssayType == "Viral",]
mrg_exposure_virus$VirusNumber[mrg_exposure_virus$AcquisitionLonger == 1] <- mrg_exposure_virus0$VirusNumber[mrg_exposure_virus0$AcquisitionLonger == 1] - 1
mrg_exposure_virus$Virus <- 0 
mrg_exposure_virus$Virus[mrg_exposure_virus$VirusNumber >0] <- 1 #make a Spn yes.no for viral comparison
mrg_exposure_bocavirus <- mrg_exposure[mrg_exposure$Assay == "bocavirus",]
mrg_exposure_HRV <- mrg_exposure[mrg_exposure$Assay == "HRV",]
mrg_exposure_staph <- mrg_exposure[mrg_exposure$Assay == "Staphylococcus_aureus",]
mrg_exposure_pyogenes <- mrg_exposure[mrg_exposure$Assay == "Streptococcus_pyogenes",]

#stats:
m <- glmer(as.factor(AcquisitionLonger) ~ AgeYears + VirusNumber + (1|StudyID), data = mrg_exposure_bocavirus, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

table2 <- tidy(m) %>%
  filter(effect == "fixed") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value)
table2

write_tsv(table2, "Fig2I_statssummary.txt")

coefs <- summary(m)$coefficients

est <- coefs["AgeYears", "Estimate"]
p   <- coefs["AgeYears", "Pr(>|z|)"]

plot_title <- paste0(
  "Bocavirus\n",
  "Est: ", round(est, 3),
  ", P=", sapply(p, disp_pval)
)

Fig2I <- ggplot(PerChild2, aes(x = as.factor(AgeYears), y = bocavirusacquisitionSucces)) + 
  geom_boxplot(outlier.color = NA, aes(fill = AgeYears)) +
  geom_jitter(width = 0.2, height = 0) + 
  scale_fill_gradient2() +
  labs(
    title = plot_title,
    x = "Age (years)",
    y = "Acquisition (% of exposure)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
Fig2I

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig2I_bocavirus_acquisition_succes_", timestamp, ".pdf")
ggsave(pdf_filename, width = 2.5, height = 3)


#nr samples:
fig_data      <- mrg_exposure_bocavirus
group_cols    <- c("AgeYears")   # 2–4 columns as needed
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig2I_samplesize.txt")





###################################################
#                                                 #
#                 FIGURE 2J                       #
#                                                 #
###################################################


boca_child <- mrg %>%
  filter(Assay == "bocavirus") %>%
  group_by(StudyID, AgeYears) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(Positive > 0, na.rm = TRUE),
    n_negative = n_samples - n_positive,
    percent_positive = n_positive / n_samples,
    .groups = "drop"
  )  %>%
  mutate(AgeYears = as.numeric(AgeYears))

#stats:
m <- glmmTMB(
  cbind(n_positive, n_negative) ~ AgeYears + (1 | StudyID),
  data = boca_child,
  family = binomial()
)
summary(m)

table2 <- tidy(m) %>%
  filter(effect == "fixed") %>%
  select(term, estimate, std.error, statistic, p.value)
table2

write_tsv(table2, "Fig2J_statssummary.txt")

est <- table2$estimate[table2$term == "AgeYears"]
p   <- table2$p.value[table2$term == "AgeYears"]

plot_title <- paste0(
  "Bocavirus\n",
  "Est: ", round(est, 3),
  ", P=", sapply(p, disp_pval)
)


boca_child <- boca_child %>%
  mutate(
    AgeNum = as.numeric(as.character(AgeYears))  # numeric for gradient
  )

Fig2J <- ggplot(boca_child, aes(x = factor(AgeYears), y = percent_positive)) + 
  geom_boxplot(outlier.color = NA, aes(fill = AgeNum)) +
  geom_jitter(width = 0.2, height = 0) + 
  scale_x_discrete(drop = FALSE) +
  scale_fill_gradient2() +
  labs(
    title = plot_title,
    x = "Age (years)",
    y = "% positive samples per child"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
Fig2J

pdf_filename <- paste0("Fig2J_PerBocaviruspos_age_", timestamp, ".pdf")
ggsave(pdf_filename, width = 2.5, height = 3)


#nr samples:
fig_data      <- boca_child
group_cols    <- c("AgeYears")   # 2–4 columns as needed
samplesize <- count_per_group(fig_data, group_cols)
samplesize
capture.output(samplesize, file = "Fig2J_samplesize.txt")



###################################################
#                                                 #
#                 FIGURE 2K                       #
#                                                 #
###################################################

mrg_exposure_virus2 <- mrg_exposure_virus %>%
  mutate(AcquisitionLonger2 = if_else(AcquisitionLonger == 0, "No", "Yes")) %>%
  mutate(Staph2 = if_else(Staph == 0, "No", "Yes"))

m <- glmer(as.factor(AcquisitionLonger2) ~ Staph2 +  (1|Assay) + (1|StudyID), data = mrg_exposure_virus2, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

table2 <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(effect == "fixed") %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)
table2

write_tsv(table2, "Fig2K_statssummary.txt")

p_val <- summary(m)$coefficients["Staph2Yes", "Pr(>|z|)"]
est_logOR <- summary(m)$coefficients["Staph2Yes", "Estimate"]
se_logOR  <- summary(m)$coefficients["Staph2Yes", "Std. Error"]

OR   <- exp(est_logOR)
lower <- exp(est_logOR - 1.96 * se_logOR)
upper <- exp(est_logOR + 1.96 * se_logOR)

p_txt   <- sapply(p_val, disp_pval)
OR_txt  <- formatC(OR, format = "f", digits = 3)
lower_txt <- formatC(lower, format = "f", digits = 3)
upper_txt <- formatC(upper, format = "f", digits = 3)


#m <- glmer(as.factor(AcquisitionLonger2) ~ Staph2 + AgeYears + (1|Assay) + (1|StudyID), data = mrg_exposure_virus2, 
#           family = binomial, control = glmerControl(optimizer = "bobyqa"))
#summary(m)

table(mrg_exposure_virus2$AcquisitionLonger2, mrg_exposure_virus2$Staph2)
47/(47+35)
10/(10+18)

table(mrg_exposure_virus2$Staph2)

virus_counts <- table(mrg_exposure_virus2$Staph2)
virus_lab <- setNames(
  paste0(names(virus_counts), "\n(n=", as.integer(virus_counts), ")"),
  names(virus_counts)
)

prop.table(table(mrg_exposure_virus2$AcquisitionLonger2, mrg_exposure_virus2$Staph2), 2)

pt <- prop.table(
  table(mrg_exposure_virus2$AcquisitionLonger2,
        mrg_exposure_virus2$Staph2), 2
)

ann <- as.data.frame(pt)
names(ann)[1:2] <- c("AcquisitionLonger2", "Staph2")
ann <- ann %>%
  mutate(label = paste0(round(Freq * 100, 1), "%"))

p0 <- ggplot(mrg_exposure_virus2) +
  geom_mosaic(aes(x = product(AcquisitionLonger2, Staph2),
                  fill = AcquisitionLonger2),
              show.legend = FALSE) +
  facet_grid(~Staph2, space = "free",
             labeller = labeller(Staph2 = virus_lab),
             switch = "x") +
  scale_fill_grey(start = 0.8, end = 0.2) +
  labs(title = "test", x = NULL, y = "S. pyogenes acquisition") +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 8))
  )

gb <- ggplot_build(p0)

# adjust this if needed; inspect gb$data[[1]] to match the exact columns in your setup
mosaic_df <- gb$data[[1]]

mosaic_df <- mosaic_df %>%
  mutate(
    Staph2 = rep(levels(factor(mrg_exposure_virus2$Staph2)),
                 each = length(unique(mrg_exposure_virus2$AcquisitionLonger2)))
  )

Fig2K <- p0 +
  geom_text(
    data = mosaic_df,
    aes(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = ann$label
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    title = paste0("OR:", OR_txt, "(CI:", lower_txt, "–", upper_txt, "), P= ", p_txt),
    x = "Virus carriage",
    y = "S. pyogenes acquisition"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
Fig2K

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig2K_", timestamp, ".pdf")
ggsave(pdf_filename, width = 2, height = 3)



###################################################
#                                                 #
#                 FIGURE 2L                       #
#                                                 #
###################################################

mrg_exposure_pyogenes2 <- mrg_exposure_pyogenes %>%
  mutate(AcquisitionLonger2 = if_else(AcquisitionLonger == 0, "No", "Yes")) %>%
  mutate(Virus2 = if_else(Virus == 0, "No", "Yes"))

m <- glmer(as.factor(AcquisitionLonger2) ~  Virus2  + (1|AgeYears) +  (1|StudyID), data = mrg_exposure_pyogenes2, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m)

table2 <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(effect == "fixed") %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)
table2

write_tsv(table2, "Fig2L_statssummary.txt")

p_val <- summary(m)$coefficients["Virus2Yes", "Pr(>|z|)"]
est_logOR <- summary(m)$coefficients["Virus2Yes", "Estimate"]
se_logOR  <- summary(m)$coefficients["Virus2Yes", "Std. Error"]

OR   <- exp(est_logOR)
lower <- exp(est_logOR - 1.96 * se_logOR)
upper <- exp(est_logOR + 1.96 * se_logOR)

p_txt   <- sapply(p_val, disp_pval)
OR_txt  <- formatC(OR, format = "f", digits = 3)
lower_txt <- formatC(lower, format = "f", digits = 3)
upper_txt <- formatC(upper, format = "f", digits = 3)

#m <- glmer(as.factor(AcquisitionLonger2) ~  Virus2  +  (1|StudyID), data = mrg_exposure_pyogenes2, 
#           family = binomial, control = glmerControl(optimizer = "bobyqa"))
#summary(m)

table(mrg_exposure_pyogenes2$AcquisitionLonger2, mrg_exposure_pyogenes2$Virus2)
table(mrg_exposure_pyogenes2$Virus2)
virus_counts <- table(mrg_exposure_pyogenes2$Virus2)
virus_lab <- setNames(
  paste0(names(virus_counts), "\n(n=", as.integer(virus_counts), ")"),
  names(virus_counts)
)

prop.table(table(mrg_exposure_pyogenes2$AcquisitionLonger2, mrg_exposure_pyogenes2$Virus2), 2)

pt <- prop.table(
  table(mrg_exposure_pyogenes2$AcquisitionLonger2,
        mrg_exposure_pyogenes2$Virus2), 2
)

ann <- as.data.frame(pt)
names(ann)[1:2] <- c("AcquisitionLonger2", "Virus2")
ann <- ann %>%
  mutate(label = paste0(round(Freq * 100, 1), "%"))

p0 <- ggplot(mrg_exposure_pyogenes2) +
  geom_mosaic(aes(x = product(AcquisitionLonger2, Virus2),
                  fill = AcquisitionLonger2),
              show.legend = FALSE) +
  facet_grid(~Virus2, space = "free",
             labeller = labeller(Virus2 = virus_lab),
             switch = "x") +
  scale_fill_grey(start = 0.8, end = 0.2) +
  labs(title = "test", x = NULL, y = "S. pyogenes acquisition") +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 8))
  )

gb <- ggplot_build(p0)

# adjust this if needed; inspect gb$data[[1]] to match the exact columns in your setup
mosaic_df <- gb$data[[1]]

mosaic_df <- mosaic_df %>%
  mutate(
    Virus2 = rep(levels(factor(mrg_exposure_pyogenes2$Virus2)),
                 each = length(unique(mrg_exposure_pyogenes2$AcquisitionLonger2)))
  )

Fig2L <- p0 +
  geom_text(
    data = mosaic_df,
    aes(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = ann$label
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    title = paste0("OR:", OR_txt, "(CI:", lower_txt, "–", upper_txt, "), P= ", p_txt),
    x = "Virus carriage",
    y = "S. pyogenes acquisition"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
Fig2L

timestamp <- format(Sys.time(), "%Y%m%d")
pdf_filename <- paste0("Fig2L_", timestamp, ".pdf")
ggsave(pdf_filename, width = 2, height = 3)




############################################################

save.image('figure2_exposures.Rdata')
