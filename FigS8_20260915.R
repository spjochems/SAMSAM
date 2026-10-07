#seasonality:

rm(list=ls())

setwd("\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects\\SAMSAM\\Manuscripts\\FiguresLM2026\\FigS8\\")

# Library
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(lmerTest)
library(reshape2)
library(readxl)
library(lmms)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(tibble)
library(broom.mixed)
library(patchwork)
library(glmmTMB)


#save the session info
date_str <- format(Sys.Date(), "%Y-%m-%d")
fname <- paste0("sessionInfo_", date_str, ".txt")
sink(fname)
sessionInfo()
sink()

#load in database
mrg <- data.frame(read_xlsx('\\\\vf-lucid-r-i.lumcnet.prod.intern\\lucid-r-i$\\Projects/SAMSAM/IntegrativeAnalysis/3.NA_detection/Database/ClinicalDb_Biomark_incDefsSkip3days_20260816.xlsx'))
mrg <- mrg[mrg$Antibiotics != 'Yes',] #remove the antibiotics samples as it may affect results

df_season <- mrg %>%
  #filter(Day == "D01") %>%
  dplyr::mutate(month = as.integer(str_split_fixed(DateCollection_EditDate, "_", 3)[, 2])) %>%
  mutate(
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Autumn",
      TRUE ~ NA_character_
    )
  )%>%
  dplyr::mutate(winter = if_else(season == "Winter", "Winter", "Not winter")) %>%
  dplyr::arrange(month) %>%
  dplyr::mutate(month2 = factor(month,
                                levels = 1:12,
                                labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))) %>%
  dplyr::mutate(month2 = factor(month2, levels = month.abb))

#enocde the season as a cosinus
df_season$DateCollection_EditDate #check what the date looks like
df_season$DaysSinceJan01 <- as.integer(format(as.Date(df_season$DateCollection_EditDate, format = "%Y_%m_%d"), "%j")) - 1 #set it to an integer
df_season$DaysCosinus <-  cos(2*pi*df_season$DaysSinceJan01/365) #change it to a cosinus

#create a function that fixes the display of p-values
disp_pval <- function(p){
  if(p<0.001){
    p <- scales::scientific(p, digits= 3)
  } else {
    p <- round(p, 3)
  }
  return(p)
}



###################################################
#                                                 #
#                 FIGURE S8A                      #
#                                                 #
###################################################


#select pathogens
df_season2 <- df_season[grepl("16S", df_season$Assay) | 
                          grepl("^Spn\\(lytA\\)$", df_season$Assay) |
                          grepl("Haemophilus", df_season$Assay) |
                          grepl("pyogenes", df_season$Assay) |
                          grepl("Moraxella", df_season$Assay)|
                          grepl("aureus", df_season$Assay) |
                          grepl("HRV", df_season$Assay) |
                          grepl("bocavirus", df_season$Assay)
                        , ]
df_season2 <- df_season2[df_season2$log10_avgconc>0,]

#check if this makes sense
P1 <- ggplot(df_season2, aes(x=as.factor(month), y=DaysCosinus)) + geom_boxplot() #plot it to check
P2 <- ggplot(df_season2, aes(x=as.factor(month), y=DaysSinceJan01)) + geom_boxplot() #plot it to check
P3 <- ggplot(df_season2, aes(x=DaysSinceJan01, y=DaysCosinus)) + geom_point() #plot it to check

combined <- P1 | P2 | P3
combined

ggsave(
  "FigS8A_SeasonCosinusConversion3.pdf",
  plot = combined,
  width = 9,   # wider to fit three panels
  height = 3
)



###################################################
#                                                 #
#                 FIGURE S8B                      #
#                                                 #
###################################################

# stats seasonality:
stats <- data.frame(Estimate = NA, Pval = NA)

df_season3 <- df_season2 %>%
  filter(Assay %in% c("Haemophilus_influenzae", "Moraxella_catarrhalis","Spn(lytA)", "Staphylococcus_aureus"))

unique(df_season3$Assay)

df_season2.split <- split(df_season3, f=df_season3$Assay)

#stats <- NULL  # <-- reset if rerun

for(i in 1:length(df_season2.split)){
  test <- df_season2.split[[i]]
  m <- lmer(log10_avgconc ~ DaysCosinus + AgeYears +(1|StudyID), data = test)
  res <- lmerTest:::get_coefmat(m)[2,c(1,5)]
  stats <- rbind(stats, res)
}

stats_season <- stats[-1,]
rownames(stats_season) <- stats_season$Assay <- names(df_season2.split)
stats_season$Pval_adj <- p.adjust(stats_season$Pval, method = "BH")
stats_season$Label = paste0('Est:', round(stats_season$Estimate, 3), ', P=', sapply(stats_season$Pval_adj, disp_pval))
mrg2.season <- merge(df_season2, stats_season, by= 'Assay')
mrg2.season$Label2 <- mrg2.season$Label
mrg2.season$Label2[duplicated(mrg2.season$Label)] <- NA


mrg2.season$Assay <- factor(mrg2.season$Assay,
                            levels = c("Haemophilus_influenzae", "Moraxella_catarrhalis","Spn(lytA)", "Staphylococcus_aureus"))


FigS8B <- ggplot(mrg2.season, aes(x = DaysSinceJan01, y = log10_avgconc)) +
  facet_wrap(. ~ Assay,  ncol=4, scale = 'free_y') +
  geom_point(aes(colour = StudyID)) +
  stat_smooth(method = "lm", formula = y ~ splines::ns(x, 3), se = TRUE, colour = "blue") + #formula = y ~ splines::ns(x, 3) tells it to fit a natural spline with 3 basis degrees of freedom, which gives a curved line.
  ggtitle("Effect of seasonality, all kids/samples") +
  geom_text(aes(label = Label2), y = Inf, x = Inf, hjust = 1, vjust = 1, colour = "red") +
  theme(
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text = element_text(colour = "black")) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))  # more space at top
  )
FigS8B

ggsave('FigS8B_PerSeason.pdf', width = 14, height = 3)



###################################################
#                                                 #
#                 FIGURE S8C                      #
#                                                 #
###################################################


#Samples per month:
monthly_positive <- df_season %>%
  filter(
    Assay %in% c(
      "Haemophilus_influenzae",
      "Moraxella_catarrhalis",
      "Spn(lytA)",
      "Staphylococcus_aureus"
    ),
    !is.na(InfectionLongerTrue),
    !is.na(month)
  ) %>%
  mutate(
    month = as.character(month),
    positive = InfectionLongerTrue == 1
  ) %>%
  group_by(StudyID,Assay, month) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(positive),
    n_negative = n_samples - n_positive,
    percent_positive = 100* n_positive / n_samples,
    .groups = "drop"
  ) %>%
  arrange(Assay, month)

monthly_positive

monthly_positive2 <- monthly_positive %>%
  mutate(
    month = factor(
      as.character(month),
      levels = as.character(1:12)
    )
  )

p1 <- ggplot(monthly_positive2, aes(x = month, y = percent_positive)) + 
  facet_wrap(. ~ Assay, scales = "free", ncol = 5) + 
  geom_violin(aes(fill = month)) +  
  geom_boxplot(width = 0.2) + 
  ggtitle("Carriage positivity per month") +
  #geom_text(aes(label = Label2), y = Inf, x = 1.5, hjust =0.5, vjust = 1.2, colour = "red") +
  theme(legend.position = "none", strip.background = element_rect(fill = NA, colour = NA))+
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0.0, 0.15))
  )
p1

#stats:
carriage_season_data <- df_season %>%
  filter(
    Assay %in% c(
      "Haemophilus_influenzae",
      "Moraxella_catarrhalis",
      "Spn(lytA)",
      "Staphylococcus_aureus"
    ),
    !is.na(InfectionLongerTrue),
    !is.na(DaysCosinus),
    !is.na(StudyID)
  ) %>%
  mutate(
    positive = as.integer(InfectionLongerTrue == 1)
  )
carriage_season_data

HI_carriage <- carriage_season_data %>% filter(Assay == "Haemophilus_influenzae")
m_HI <- glmer(
  positive ~ DaysCosinus + AgeYears + (1 | StudyID),
  data = HI_carriage,
  family = binomial
)
summary(m_HI)
#capture.output(summary(m_HI), file = "FigS5C_HI_statssummary_Spn.txt")


MC_carriage <- carriage_season_data %>% filter(Assay == "Moraxella_catarrhalis")
m_MC<- glmer(
  positive ~ DaysCosinus +  (1 | StudyID),
  data = MC_carriage,
  family = binomial
)
summary(m_MC)
#capture.output(summary(m_MC), file = "FigS5C_MC_statssummary_Spn.txt")


Spn_carriage <- carriage_season_data %>% filter(Assay == "Spn(lytA)")
m_Spn <- glmer(
  positive ~ DaysCosinus +  (1 | StudyID),
  data = Spn_carriage,
  family = binomial
)
summary(m_Spn)
#capture.output(summary(m_Spn), file = "FigS5C_Spn_statssummary_Spn.txt")


SA_carriage <- carriage_season_data %>% filter(Assay == "Staphylococcus_aureus")
m_SA <- glmer(
  positive ~ DaysCosinus +  (1 | StudyID),
  data = SA_carriage,
  family = binomial
)
summary(m_SA)
#capture.output(summary(m_SA), file = "FigS5C_SA_statssummary_Spn.txt")


disp_pval2 <- function(p) {
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

#collect P values:
season_results <- dplyr::bind_rows(
  broom.mixed::tidy(m_HI, effects = "fixed") %>%
    filter(term == "DaysCosinus") %>%
    mutate(Assay = "Haemophilus_influenzae"),
  
  broom.mixed::tidy(m_MC, effects = "fixed") %>%
    filter(term == "DaysCosinus") %>%
    mutate(Assay = "Moraxella_catarrhalis"),
  
  broom.mixed::tidy(m_Spn, effects = "fixed") %>%
    filter(term == "DaysCosinus") %>%
    mutate(Assay = "Spn(lytA)"),
  
  broom.mixed::tidy(m_SA, effects = "fixed") %>%
    filter(term == "DaysCosinus") %>%
    mutate(Assay = "Staphylococcus_aureus")
) %>%
  dplyr::select(
    Assay,
    estimate,
    std.error,
    statistic,
    p.value
  ) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    estimate_fmt = sprintf("%.3f", estimate),
    p.value_fmt = disp_pval2(p.value),
    p.adj_fmt = disp_pval2(p.adj),
    stars = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Label = paste0(
      "Est:", estimate_fmt,
      ", P=", p.adj_fmt
    )
  )
season_results
capture.output(summary(season_results), file = "FigS8C_ALL_statssummary_Spn.txt")

season_labels <- season_results %>%
  mutate(
    month = factor(
      "1",
      levels = as.character(1:12)
    ),
    y_position = 98
  )

#plot:
p1 <- ggplot(
  monthly_positive2,
  aes(x = month, y = percent_positive)
) +
  facet_wrap(. ~ Assay, scales = "free", ncol = 5) +
  geom_violin(aes(fill = month)) +
  geom_boxplot(width = 0.2) +
  geom_text(
    data = season_labels,
    aes(
      x = month,
      y = y_position,
      label = Label
    ),
    inherit.aes = FALSE,
    hjust = -0.4,
    vjust = -1.2,
    size = 3.5
  ) +
  ggtitle("Carriage positivity per month") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = NA, colour = NA)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.15))
  )
p1

ggsave('FigS8C_PerSeasonCarriage.pdf', width = 10, height = 3)




###################################################
#                                                 #
#                 FIGURE S8D                      #
#                                                 #
###################################################
#season virus exposure

mrg2 <-  mrg1 <- df_season

# Keep all viral assay observations
mrg2 <- df_season %>%
  filter(AssayType == "Viral") %>%
  mutate(
    exposed = ExposureEvent == 1, 
    Month_num = parse_number(as.character(month))
  )


#make a model to see if the infection of viruses differs per season
m2 <- glmer(ExposureEvent ~ DaysCosinus + AgeYears + (1|Assay) +(1|StudyID), data = mrg2, family = 'binomial')
summary(m2)

table2 <- tidy(m2) %>%
  filter(effect == "fixed") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    
    std.error = sprintf("%.3f", std.error),
    estimate = sprintf("%.3f",estimate),
    statistic = disp_pval2(statistic),
    p.value    = disp_pval2(p.value),
    p.adj    = disp_pval2(p.adj),
    
    stars = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    ))
table2

capture.output(table2, file = "FigS8D_season_viralexposure_statssummary.txt")

# check per pathogen:
assays <- unique(mrg2$Assay)

res_list <- lapply(assays, function(a) {
  d <- mrg2 %>% filter(Assay == a)
  
  # Check variation in outcome
  y_tab <- table(d$ExposureEvent)
  
  # Try to fit model, with tryCatch to avoid stopping the loop
  mb <- tryCatch(
    glmer(
      ExposureEvent ~ DaysCosinus + AgeYears + (1 | StudyID),
      data = d,
      family = binomial,
      control = glmerControl(
        optimizer = "bobyqa",
        check.conv.grad = "ignore",      # relax gradient check
        check.conv.singular = "ignore"   # allow near-singular fits
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(mb)) {
    return(NULL)
  }
  
  # Optionally check convergence code
  if (!is.null(mb@optinfo$conv$lme4$messages) &&
      length(mb@optinfo$conv$lme4$messages) > 0) {
    # You can either skip or keep with a warning flag; here we keep.
  }
  
  t <- tidy(mb, effects = "fixed") %>%
    filter(term == "DaysCosinus") %>%
    mutate(
      Assay = a,
      converged = isTRUE(mb@optinfo$conv$lme4$conv == 0)
    )
  
  t
})

season_results_byAssay <- bind_rows(res_list) %>%
  select(
    Assay,
    estimate,
    std.error,
    statistic,
    p.value
  ) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    
    std.error = sprintf("%.3f", std.error),
    estimate = sprintf("%.3f", estimate),
    statistic = sprintf("%.3f", statistic),
    p.value    = disp_pval2(p.value),
    p.adj    = disp_pval2(p.adj),
    
    stars = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    
    Label = paste0(
      "Est:", estimate,
      ", P=", p.adj
    )
  )
season_results_byAssay

capture.output(season_results_byAssay, file = "FigS8D_season_viralexposure_PerAssay_statssummary.txt")


#Plot % exposures per virus:
sample_assay <- mrg2 %>%
  group_by(
    month,
    Month_num,
    Day2,
    StudyID,
    Assay
  ) %>%
  summarise(
    exposed = any(exposed, na.rm = TRUE), 
    .groups = "drop"
  )

percentage_individual <- sample_assay %>%
  group_by(
    month,
    Month_num,
    Assay
  ) %>%
  summarise(
    n_samples = n(),
    n_exposed = sum(exposed),
    percent_exposures = 100 * n_exposed / n_samples,
    .groups = "drop"
  )

# % exposures virus combined:
sample_combined <- sample_assay %>%
  group_by(
    month,
    Month_num,
    Day2,
    StudyID
  ) %>%
  summarise(
    exposed = any(exposed),
    .groups = "drop"
  ) %>%
  mutate(
    Assay = "All viral assays"
  )

percentage_combined <- sample_combined %>%
  group_by(
    month,
    Month_num,
    Assay
  ) %>%
  summarise(
    n_samples = n(),
    n_exposed = sum(exposed), 
    percent_exposures = 100 * n_exposed / n_samples,
    .groups = "drop"
  )

percentage_data <- bind_rows(
  percentage_individual,
  percentage_combined
)

#remove assays never positive and order assays:
assay_order_individual <- percentage_data %>%
  group_by(Assay) %>%
  summarise(
    total_exposed = sum(n_exposed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_exposed > 0) %>%
  arrange(desc(total_exposed)) %>%  
  pull(Assay)

assay_order <- c(
  assay_order_individual)

plot_data <- percentage_data %>%
  filter(Assay %in% assay_order) %>%
  mutate(
    Assay = factor(
      Assay,
      levels = assay_order
    )
  )

#make matrixes
percentage_matrix <- plot_data %>%
  dplyr::select(Assay, month, percent_exposures) %>%  
  pivot_wider(
    names_from = month,
    values_from = percent_exposures, 
    values_fill = 0
  ) %>%
  column_to_rownames("Assay") %>%
  as.matrix()

percentage_matrix <- percentage_matrix[assay_order, , drop = FALSE] 
assay_order <- intersect(assay_order, rownames(percentage_matrix)) #reorder assays
percentage_matrix <- percentage_matrix[assay_order, , drop = FALSE]

count_matrix <- plot_data %>%
  dplyr::select(Assay, month, n_exposed) %>%  
  pivot_wider(
    names_from = month,
    values_from = n_exposed,
    values_fill = 0
  ) %>%
  column_to_rownames("Assay") %>%
  as.matrix()

count_matrix <- count_matrix[
  rownames(percentage_matrix),
  colnames(percentage_matrix),
  drop = FALSE
]

# percentage positive across ALL viral samples per month
month_positive <- sample_combined %>%
  group_by(month) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(exposed, na.rm = TRUE),
    percent_positive = 100 * n_positive / n_samples,
    .groups = "drop"
  ) %>%
  slice(match(colnames(percentage_matrix), month))

# Match the annotation values exactly to heatmap-column order
month_positive <- month_positive %>%
  slice(match(colnames(percentage_matrix), month))


#percentage positive across ALL months per assay
assay_positive <- sample_assay %>%
  group_by(Assay) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(exposed, na.rm = TRUE),
    percent_positive = 100 * n_positive / n_samples,
    .groups = "drop"
  )

# Add combined viral-assay exposure row, if it is in the heatmap
combined_positive <- sample_combined %>%
  summarise(
    n_samples = n(),
    n_positive = sum(exposed, na.rm = TRUE),
    percent_positive = 100 * n_positive / n_samples
  ) %>%
  mutate(Assay = "All viral assays")

assay_positive <- bind_rows(assay_positive, combined_positive)

# Match annotation values exactly to heatmap-row order
assay_positive <- assay_positive %>%
  slice(match(rownames(percentage_matrix), Assay))

columnha <- HeatmapAnnotation(
  `% exposures/month` = anno_barplot(
    month_positive$percent_positive,
    gp = gpar(fill = "darkgrey"),
    border = FALSE,
    ylim = c(0, 15),
    height = unit(2, "cm")
  ),
  annotation_name_side = "left"
)

rowha <- rowAnnotation(
  `% exposures/pathogen` = anno_barplot(
    assay_positive$percent_positive,
    gp = gpar(fill = "darkgrey"),
    border = FALSE,
    ylim = c(0, 15),
    width = unit(3, "cm")
  )
)

#separate virus combined:
is_viral <- grepl("viral", rownames(percentage_matrix), ignore.case = TRUE)  
row_group <- factor(ifelse(is_viral, "viral", "other"),
                    levels = c("viral", "other"))

pl1 <- Heatmap(
  percentage_matrix,
  name = "% exposures",
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  col = circlize::colorRamp2(
    c(0, 2.5, 5, 7.5, 10),
    c("blue", "lightblue", "white", "orange", "red")
  ),
  
  top_annotation = columnha,
  right_annotation = rowha,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_split = row_group,
  row_gap = unit(2, "mm"),
  
  column_title = "Percentage of viral exposures per month",
  row_title = "Viral assay"
)
pl1


pdf("FigS8D_Heatmap_season_viralexposures.pdf", width = 10, height = 8)
draw(pl1)
dev.off()



###################################################
#                                                 #
#                 FIGURE S8E                      #
#                                                 #
###################################################
#season virus carriage

mrg1 <- df_season

# Keep all viral assay observations
mrg2 <- mrg1 %>%
  filter(AssayType == "Viral") %>%
  mutate(
    positive = InfectionLongerTrue == 1,
    Month_num = parse_number(as.character(month))
  )

#stats, make a model to see if the infection of viruses differs per season
m2 <- glmer(InfectionLonger ~ DaysCosinus + AgeYears + (1|Assay) +(1|StudyID), data = mrg2, family = 'binomial')
summary(m2)

table2 <- tidy(m2) %>%
  filter(effect == "fixed") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    
    std.error = sprintf("%.3f", std.error),
    estimate = sprintf("%.3f",estimate),
    statistic = disp_pval2(statistic),
    p.value    = disp_pval2(p.value),
    p.adj    = disp_pval2(p.adj),
    
    stars = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    ))
table2

capture.output(table2, file = "FigS8E_season_viralCarriage_statssummary.txt")

# check per pathogen:
assays <- unique(mrg2$Assay)

res_list <- lapply(assays, function(a) {
  d <- mrg2 %>% filter(Assay == a)
  
  # Check variation in outcome
  y_tab <- table(d$InfectionLongerTrue)
  
  # Try to fit model, with tryCatch to avoid stopping the loop
  mb <- tryCatch(
    glmer(
      InfectionLongerTrue ~ DaysCosinus + AgeYears + (1 | StudyID),
      data = d,
      family = binomial,
      control = glmerControl(
        optimizer = "bobyqa",
        check.conv.grad = "ignore",      # relax gradient check
        check.conv.singular = "ignore"   # allow near-singular fits
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(mb)) {
    return(NULL)
  }
  
  if (any(y_tab < 10)) {
    # Very few events or non-events -> skip (adjust threshold if you like)
   return(NULL)
  }
  
  # Optionally check convergence code
  if (!is.null(mb@optinfo$conv$lme4$messages) &&
      length(mb@optinfo$conv$lme4$messages) > 0) {
    # You can either skip or keep with a warning flag; here we keep.
  }
  
  t <- tidy(mb, effects = "fixed") %>%
    filter(term == "DaysCosinus") %>%
    mutate(
      Assay = a,
      converged = isTRUE(mb@optinfo$conv$lme4$conv == 0)
    )
  
  t
})

season_results_byAssay <- bind_rows(res_list) %>%
  select(
    Assay,
    estimate,
    std.error,
    statistic,
    p.value
  ) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    
    std.error = sprintf("%.3f", std.error),
    estimate = sprintf("%.3f", estimate),
    statistic = sprintf("%.3f", statistic),
    p.value    = disp_pval2(p.value),
    p.adj    = disp_pval2(p.adj),
    
    stars = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    
    Label = paste0(
      "Est:", estimate,
      ", P=", p.adj
    )
  )
season_results_byAssay

capture.output(season_results_byAssay, file = "FigS8E_season_viralCarriage_PerAssay_statssummary.txt")

#check sample size:
mrg2 %>%
  group_by(Assay) %>%
  summarise(
    n_total   = n(),
    n_pos     = sum(InfectionLongerTrue, na.rm = TRUE),
    n_donors  = n_distinct(StudyID[InfectionLongerTrue == 1]),
    .groups = "drop"
  ) %>%
  arrange(Assay)







# Plot % positive per virus:
sample_assay <- mrg2 %>%
  group_by(
    month,
    Month_num,
    Day2,
    StudyID,
    Assay
  ) %>%
  summarise(
    positive = any(positive, na.rm = TRUE),
    .groups = "drop"
  )

percentage_individual <- sample_assay %>%
  group_by(
    month,
    Month_num,
    Assay
  ) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(positive),
    percent_positive = 100 * n_positive / n_samples,
    .groups = "drop"
  )

# % positive virus combined:
sample_combined <- sample_assay %>%
  group_by(
    month,
    Month_num,
    Day2,
    StudyID
  ) %>%
  summarise(
    positive = any(positive),
    .groups = "drop"
  ) %>%
  mutate(
    Assay = "All viral assays"
  )

percentage_combined <- sample_combined %>%
  group_by(
    month,
    Month_num,
    Assay
  ) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(positive),
    percent_positive = 100 * n_positive / n_samples,
    .groups = "drop"
  )

percentage_data <- bind_rows(
  percentage_individual,
  percentage_combined
)

#remove assays never positive and order assays:
assay_order_individual <- percentage_data %>%
  group_by(Assay) %>%
  summarise(
    total_positive = sum(n_positive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_positive > 0) %>%
  arrange(desc(total_positive)) %>%
  pull(Assay)

assay_order <- c(
  assay_order_individual)

plot_data <- percentage_data %>%
  filter(Assay %in% assay_order) %>%
  mutate(
    Assay = factor(
      Assay,
      levels = assay_order
    )
  )

#make matrixes
percentage_matrix <- plot_data %>%
  dplyr::select(Assay, month, percent_positive) %>%
  pivot_wider(
    names_from = month,
    values_from = percent_positive,
    values_fill = 0
  ) %>%
  column_to_rownames("Assay") %>%
  as.matrix()

percentage_matrix <- percentage_matrix[assay_order, , drop = FALSE] 
assay_order <- intersect(assay_order, rownames(percentage_matrix)) #reorder assays
percentage_matrix <- percentage_matrix[assay_order, , drop = FALSE]

count_matrix <- plot_data %>%
  dplyr:: select(Assay, month, n_positive) %>%
  pivot_wider(
    names_from = month,
    values_from = n_positive,
    values_fill = 0
  ) %>%
  column_to_rownames("Assay") %>%
  as.matrix()

count_matrix <- count_matrix[
  rownames(percentage_matrix),
  colnames(percentage_matrix),
  drop = FALSE
]

# percentage positive across ALL viral samples per month
month_positive <- sample_combined %>%
  group_by(month) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(positive, na.rm = TRUE),
    percent_positive = 100 * n_positive / n_samples,
    .groups = "drop"
  ) %>%
  slice(match(colnames(percentage_matrix), month))

# Match the annotation values exactly to heatmap-column order
month_positive <- month_positive %>%
  slice(match(colnames(percentage_matrix), month))


#percentage positive across ALL months per assay
assay_positive <- sample_assay %>%
  group_by(Assay) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(positive, na.rm = TRUE),
    percent_positive = 100 * n_positive / n_samples,
    .groups = "drop"
  )

# Add combined viral-assay exposure row, if it is in the heatmap
combined_positive <- sample_combined %>%
  summarise(
    n_samples = n(),
    n_positive = sum(positive, na.rm = TRUE),
    percent_positive = 100 * n_positive / n_samples
  ) %>%
  mutate(Assay = "All viral assays")

assay_positive <- bind_rows(assay_positive, combined_positive)

# Match annotation values exactly to heatmap-row order
assay_positive <- assay_positive %>%
  slice(match(rownames(percentage_matrix), Assay))

columnha <- HeatmapAnnotation(
  `% positive/month` = anno_barplot(
    month_positive$percent_positive,
    gp = gpar(fill = "darkgrey"),
    border = FALSE,
    ylim = c(0, 100),
    height = unit(2, "cm")
  ),
  annotation_name_side = "left"
)

rowha <- rowAnnotation(
  `% positive/pathogen` = anno_barplot(
    assay_positive$percent_positive,
    gp = gpar(fill = "darkgrey"),
    border = FALSE,
    ylim = c(0, 75),
    width = unit(3, "cm")
  )
)

#separate virus combined:
is_viral <- grepl("viral", rownames(percentage_matrix), ignore.case = TRUE)  
row_group <- factor(ifelse(is_viral, "viral", "other"),
                    levels = c("viral", "other"))

pl1 <- Heatmap(
  percentage_matrix,
  name = "% positive",
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  col = circlize::colorRamp2(
    c(0, 15, 30, 45, 60),
    c("blue", "lightblue", "white", "orange", "red")
  ),
  
  top_annotation = columnha,
  right_annotation = rowha,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_split = row_group,
  row_gap = unit(2, "mm"),
  
  column_title = "Percentage of viral positive samples per month",
  row_title = "Viral assay"
)
pl1





