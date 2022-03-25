# Author: Marie-Madlen Pust
# Last updated: 25 March 2022

# empty global environment
rm(list = ls())

# specify required packages
req_packages <- c('readr', 'purrr', 'vegan','ggrepel', 'viridis', 'tidyr', 'Hmisc', 'dplyr', 
                  'igraph', 'matrixStats', 'ggpubr', 'ggthemes', 'hrbrthemes', 'plyr', 'readxl', 
                  'randomForest', 'Boruta', 'stringr', 'tidyverse', 'compositions', 'pheatmap', 
                  'ggplotify', 'factoextra', 'matrixStats', 'ggVennDiagram')

# define percentage for prevalence filtering
prev_filt_abund = 30

# define p-value for correlation analysis
p_corr = 0.05

# function for installing/importing packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)}

# load or install R packages
ipak(req_packages)

# import data file for the phage community
count_bacteria_phages <- read_delim("count_bacteria_phages_3.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
count_bacteria_phages <- data.frame(count_bacteria_phages)
count_bacteria_phages <- count_bacteria_phages[,colnames(count_bacteria_phages) != "Sp_282"]
rownames(count_bacteria_phages) <- count_bacteria_phages$organism
count_bacteria_phages$organism <- NULL
count_bacteria_phages$sum_rows <- rowSums(count_bacteria_phages)
count_bacteria_phages_filt <- count_bacteria_phages[count_bacteria_phages$sum_rows>0,]
count_bacteria_phages_filt$sum_rows <- NULL
count_bacteria_phages_filt_t  <- data.frame(t(count_bacteria_phages_filt))
count_bacteria_phages_filt_t <- count_bacteria_phages_filt_t[order(rownames(count_bacteria_phages_filt_t)),]

# import meta data
meta_data_sp <- read_delim("meta_data_sp.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
meta_data_sp <- data.frame(meta_data_sp)
remove_meta_data <- c("Sp_291", "Sp_104", "Sp_823")
meta_data_sp <- meta_data_sp[!meta_data_sp$sample_name %in% remove_meta_data,]
meta_data_sp <- meta_data_sp[order(meta_data_sp$sample_name),]
# filter meta data
meta_data_sp$age_group <- with(meta_data_sp, ifelse(age_years > 9 & age_years < 20, "10_19",
                                                ifelse(age_years >19 & age_years < 30, "20_29",
                                                    ifelse(age_years > 29 & age_years < 40, "30_39",
                                                        ifelse(age_years > 39 & age_years < 50, "40_49",
                                                            ifelse(age_years > 49, "50_60", NA))))))
meta_data_sp_noBlank <- meta_data_sp[!grepl("Blank",meta_data_sp$sample_name),]
rownames(meta_data_sp_noBlank) <- meta_data_sp_noBlank$sample_name
meta_data_sp_noBlank$sample_name <- NULL
count_bacteria_phages_filt_t$age_group <- meta_data_sp$age_group

# prevalence filtering (age group: 10-14)
count_bacteria_phages_filt_t_10_19 <- subset(count_bacteria_phages_filt_t, age_group == "10_19")
count_bacteria_phages_filt_t_10_19$age_group <- NULL
count_bacteria_phages_filt_10_19 <- data.frame(t(count_bacteria_phages_filt_t_10_19))
n_10_19=round(length(colnames(count_bacteria_phages_filt_10_19)) * prev_filt_abund / 100) 
count_bacteria_phages_filt_10_19[count_bacteria_phages_filt_10_19 == 0] <- NA
count_bacteria_phages_filt_10_19 <- count_bacteria_phages_filt_10_19[rowSums(is.na(count_bacteria_phages_filt_10_19)) < n_10_19, ]
count_bacteria_phages_filt_10_19[is.na(count_bacteria_phages_filt_10_19)] <- 0
count_bacteria_phages_filt_10_19$Species <- rownames(count_bacteria_phages_filt_10_19)

# prevalence filtering (age group: 20-24)
count_bacteria_phages_filt_t_20_29 <- subset(count_bacteria_phages_filt_t, age_group == "20_29")
count_bacteria_phages_filt_t_20_29$age_group <- NULL
count_bacteria_phages_filt_20_29 <- data.frame(t(count_bacteria_phages_filt_t_20_29))
n_20_29=round(length(colnames(count_bacteria_phages_filt_20_29)) * prev_filt_abund / 100) 
count_bacteria_phages_filt_20_29[count_bacteria_phages_filt_20_29 == 0] <- NA
count_bacteria_phages_filt_20_29 <- count_bacteria_phages_filt_20_29[rowSums(is.na(count_bacteria_phages_filt_20_29)) < n_20_29, ]
count_bacteria_phages_filt_20_29[is.na(count_bacteria_phages_filt_20_29)] <- 0
count_bacteria_phages_filt_20_29$Species <- rownames(count_bacteria_phages_filt_20_29)

# prevalence filtering (age group: 30_39)
count_bacteria_phages_filt_t_30_39 <- subset(count_bacteria_phages_filt_t, age_group == "30_39")
count_bacteria_phages_filt_t_30_39$age_group <- NULL
count_bacteria_phages_filt_30_39 <- data.frame(t(count_bacteria_phages_filt_t_30_39))
n_30_39=round(length(colnames(count_bacteria_phages_filt_30_39)) * prev_filt_abund / 100) 
count_bacteria_phages_filt_30_39[count_bacteria_phages_filt_30_39 == 0] <- NA
count_bacteria_phages_filt_30_39 <- count_bacteria_phages_filt_30_39[rowSums(is.na(count_bacteria_phages_filt_30_39)) < n_30_39, ]
count_bacteria_phages_filt_30_39[is.na(count_bacteria_phages_filt_30_39)] <- 0
count_bacteria_phages_filt_30_39$Species <- rownames(count_bacteria_phages_filt_30_39)

# prevalence filtering (age group: 40_49)
count_bacteria_phages_filt_t_40_49 <- subset(count_bacteria_phages_filt_t, age_group == "40_49")
count_bacteria_phages_filt_t_40_49$age_group <- NULL
count_bacteria_phages_filt_40_49 <- data.frame(t(count_bacteria_phages_filt_t_40_49))
n_40_49=round(length(colnames(count_bacteria_phages_filt_40_49)) * prev_filt_abund / 100) 
count_bacteria_phages_filt_40_49[count_bacteria_phages_filt_40_49 == 0] <- NA
count_bacteria_phages_filt_40_49 <- count_bacteria_phages_filt_40_49[rowSums(is.na(count_bacteria_phages_filt_40_49)) < n_40_49, ]
count_bacteria_phages_filt_40_49[is.na(count_bacteria_phages_filt_40_49)] <- 0
count_bacteria_phages_filt_40_49$Species <- rownames(count_bacteria_phages_filt_40_49)

# prevalence filtering (age group: 50_60)
count_bacteria_phages_filt_t_50_60 <- subset(count_bacteria_phages_filt_t, age_group == "50_60")
count_bacteria_phages_filt_t_50_60$age_group <- NULL
count_bacteria_phages_filt_50_60 <- data.frame(t(count_bacteria_phages_filt_t_50_60))
n_50_60=round(length(colnames(count_bacteria_phages_filt_50_60)) * prev_filt_abund / 100) 
count_bacteria_phages_filt_50_60[count_bacteria_phages_filt_50_60 == 0] <- NA
count_bacteria_phages_filt_50_60 <- count_bacteria_phages_filt_50_60[rowSums(is.na(count_bacteria_phages_filt_50_60)) < n_50_60, ]
count_bacteria_phages_filt_50_60[is.na(count_bacteria_phages_filt_50_60)] <- 0
count_bacteria_phages_filt_50_60$Species <- rownames(count_bacteria_phages_filt_50_60)

# prevalence filtering (blank controls)
count_bacteria_phages_filt_t_blank <- count_bacteria_phages_filt_t[grepl("Blank", rownames(count_bacteria_phages_filt_t)),]
count_bacteria_phages_filt_t_blank$age_group <- NULL
count_bacteria_phages_filt_blank <- data.frame(t(count_bacteria_phages_filt_t_blank))
n_blank=round(length(colnames(count_bacteria_phages_filt_blank)) * prev_filt_abund / 100) 
count_bacteria_phages_filt_blank[count_bacteria_phages_filt_blank == 0] <- NA
count_bacteria_phages_filt_blank <- count_bacteria_phages_filt_blank[rowSums(is.na(count_bacteria_phages_filt_blank)) < n_blank, ]
count_bacteria_phages_filt_blank[is.na(count_bacteria_phages_filt_blank)] <- 0
count_bacteria_phages_filt_blank$Species <- rownames(count_bacteria_phages_filt_blank)
blank_species_sp <- rownames(count_bacteria_phages_filt_blank)

# re-merge all data-frames
df_list <- list(count_bacteria_phages_filt_10_19, 
                count_bacteria_phages_filt_20_29,
                count_bacteria_phages_filt_30_39, 
                count_bacteria_phages_filt_40_49, 
                count_bacteria_phages_filt_50_60)  

new_df = df_list %>% reduce(full_join, by='Species')
new_df <- new_df[!new_df$Species %in% blank_species_sp,]

# clean species names
new_df$Species <- str_replace_all(new_df$Species, "chromosome", "")
new_df$Species <- str_replace_all(new_df$Species, "genome", "")
new_df$Species <- str_replace_all(new_df$Species, "complete", "")
new_df$Species <- str_replace_all(new_df$Species, "sequence", "")
new_df$Species <- str_replace_all(new_df$Species, "BAC", "")
new_df$Species <- str_replace_all(new_df$Species, "strain", "")
new_df$Species <- str_replace_all(new_df$Species, "subsp", "")
new_df$Species <- str_replace_all(new_df$Species, "ATCC", "")
new_df$Species <- str_replace_all(new_df$Species, "DSM", "")
new_df$Species <- str_replace_all(new_df$Species, "____", "")
new_df$Species0 <- unlist(map(str_split(new_df$Species, "_"),3))
new_df$Species1 <- unlist(map(str_split(new_df$Species, "_"),4))
new_df$Species2 <- unlist(map(str_split(new_df$Species, "_"),5))

new_df$Genus <- ifelse(new_df$Species0 == "1" | new_df$Species0 == "2" | new_df$Species0 == "3", new_df$Species1, new_df$Species0)
new_df$Species3 <- ifelse(grepl("phage", new_df$Species), "phage",
                          ifelse(grepl("virus", new_df$Species), "virus", new_df$Species2))
new_df$Species4 <- paste(new_df$Genus, new_df$Species3 )
new_df$Species5 <- ifelse(new_df$Species3=="", paste(new_df$Species4, new_df$Species1),new_df$Species4)

new_df$Species <- NULL
new_df$Species0 <- NULL
new_df$Species1 <- NULL
new_df$Species2 <- NULL
new_df$Species3 <- NULL
new_df$Species4 <- NULL
new_df$Genus <- NULL

new_df[is.na(new_df)] <- 0
new_df_filt <- ddply(new_df,"Species5",numcolwise(sum))
new_df_filt$Species5 <- str_replace(new_df_filt$Species5, "Human virus", "Human herpes virus")
rownames(new_df_filt) <- new_df_filt$Species5
new_df_filt$Species5 <- NULL

new_df_filt_phages <- new_df_filt[grepl("phage",rownames(new_df_filt)),]
new_df_filt_viruses <- new_df_filt[grepl("virus",rownames(new_df_filt)),]
new_df_filt_phages_viruses <- data.frame(rbind(new_df_filt_phages, new_df_filt_viruses))

new_df_filt_bacteria <- new_df_filt[!grepl("phage",rownames(new_df_filt)),]
new_df_filt_bacteria <- new_df_filt_bacteria[!grepl("virus",rownames(new_df_filt_bacteria)),]


# species abundance filtering (BACTERIA) #
new_df_filt_bacteria_t <- data.frame(t(new_df_filt_bacteria))
new_df_filt_bacteria_t$age_group <- meta_data_sp_noBlank$age_group

# species abundance filtering, 10_19
new_df_filt_bacteria_t_10_19_t <- subset(new_df_filt_bacteria_t, age_group == "10_19")
new_df_filt_bacteria_t_10_19_t$age_group <- NULL
new_df_filt_bacteria_t_10_19 <- data.frame(t(new_df_filt_bacteria_t_10_19_t))
new_df_filt_bacteria_t_10_19$sum_row <- rowSums(new_df_filt_bacteria_t_10_19) 
new_df_filt_bacteria_t_10_19$sum_row_2 <- new_df_filt_bacteria_t_10_19$sum_row / sum(new_df_filt_bacteria_t_10_19$sum_row) * 100
new_df_filt_bacteria_t_10_19 <- new_df_filt_bacteria_t_10_19[order(-new_df_filt_bacteria_t_10_19$sum_row_2),]
new_df_filt_bacteria_t_10_19$cumsum_col <- cumsum(new_df_filt_bacteria_t_10_19$sum_row_2)
count_bacteria_filt_10_19_abundant <- new_df_filt_bacteria_t_10_19[new_df_filt_bacteria_t_10_19$cumsum_col<95,]
count_bacteria_filt_10_19_abundant$Species <- rownames(count_bacteria_filt_10_19_abundant)
count_bacteria_filt_10_19_abundant$cumsum_col <- NULL
count_bacteria_filt_10_19_abundant$sum_row <- NULL
count_bacteria_filt_10_19_abundant$sum_row_2 <- NULL
count_bacteria_filt_10_19_rare <- new_df_filt_bacteria_t_10_19[new_df_filt_bacteria_t_10_19$cumsum_col>=95,]
count_bacteria_filt_10_19_rare$Species <- rownames(count_bacteria_filt_10_19_rare)
count_bacteria_filt_10_19_rare$cumsum_col <- NULL
count_bacteria_filt_10_19_rare$sum_row <- NULL
count_bacteria_filt_10_19_rare$sum_row_2 <- NULL

# species abundance filtering, 20_29
new_df_filt_bacteria_t_20_29_t <- subset(new_df_filt_bacteria_t, age_group == "20_29")
new_df_filt_bacteria_t_20_29_t$age_group <- NULL
new_df_filt_bacteria_t_20_29 <- data.frame(t(new_df_filt_bacteria_t_20_29_t))
new_df_filt_bacteria_t_20_29$sum_row <- rowSums(new_df_filt_bacteria_t_20_29) 
new_df_filt_bacteria_t_20_29$sum_row_2 <- new_df_filt_bacteria_t_20_29$sum_row / sum(new_df_filt_bacteria_t_20_29$sum_row) * 100
new_df_filt_bacteria_t_20_29 <- new_df_filt_bacteria_t_20_29[order(-new_df_filt_bacteria_t_20_29$sum_row_2),]
new_df_filt_bacteria_t_20_29$cumsum_col <- cumsum(new_df_filt_bacteria_t_20_29$sum_row_2)
count_bacteria_filt_20_29_abundant <- new_df_filt_bacteria_t_20_29[new_df_filt_bacteria_t_20_29$cumsum_col<95,]
count_bacteria_filt_20_29_abundant$Species <- rownames(count_bacteria_filt_20_29_abundant)
count_bacteria_filt_20_29_abundant$cumsum_col <- NULL
count_bacteria_filt_20_29_abundant$sum_row <- NULL
count_bacteria_filt_20_29_abundant$sum_row_2 <- NULL
count_bacteria_filt_20_29_rare <- new_df_filt_bacteria_t_20_29[new_df_filt_bacteria_t_20_29$cumsum_col>=95,]
count_bacteria_filt_20_29_rare$Species <- rownames(count_bacteria_filt_20_29_rare)
count_bacteria_filt_20_29_rare$cumsum_col <- NULL
count_bacteria_filt_20_29_rare$sum_row <- NULL
count_bacteria_filt_20_29_rare$sum_row_2 <- NULL

# species abundance filtering, 30_39
new_df_filt_bacteria_t_30_39_t <- subset(new_df_filt_bacteria_t, age_group == "30_39")
new_df_filt_bacteria_t_30_39_t$age_group <- NULL
new_df_filt_bacteria_t_30_39 <- data.frame(t(new_df_filt_bacteria_t_30_39_t))
new_df_filt_bacteria_t_30_39$sum_row <- rowSums(new_df_filt_bacteria_t_30_39) 
new_df_filt_bacteria_t_30_39$sum_row_2 <- new_df_filt_bacteria_t_30_39$sum_row / sum(new_df_filt_bacteria_t_30_39$sum_row) * 100
new_df_filt_bacteria_t_30_39 <- new_df_filt_bacteria_t_30_39[order(-new_df_filt_bacteria_t_30_39$sum_row_2),]
new_df_filt_bacteria_t_30_39$cumsum_col <- cumsum(new_df_filt_bacteria_t_30_39$sum_row_2)
count_bacteria_filt_30_39_abundant <- new_df_filt_bacteria_t_30_39[new_df_filt_bacteria_t_30_39$cumsum_col<95,]
count_bacteria_filt_30_39_abundant$Species <- rownames(count_bacteria_filt_30_39_abundant)
count_bacteria_filt_30_39_abundant$cumsum_col <- NULL
count_bacteria_filt_30_39_abundant$sum_row <- NULL
count_bacteria_filt_30_39_abundant$sum_row_2 <- NULL
count_bacteria_filt_30_39_rare <- new_df_filt_bacteria_t_30_39[new_df_filt_bacteria_t_30_39$cumsum_col>=95,]
count_bacteria_filt_30_39_rare$Species <- rownames(count_bacteria_filt_30_39_rare)
count_bacteria_filt_30_39_rare$cumsum_col <- NULL
count_bacteria_filt_30_39_rare$sum_row <- NULL
count_bacteria_filt_30_39_rare$sum_row_2 <- NULL

# species abundance filtering, 40_49
new_df_filt_bacteria_t_40_49_t <- subset(new_df_filt_bacteria_t, age_group == "40_49")
new_df_filt_bacteria_t_40_49_t$age_group <- NULL
new_df_filt_bacteria_t_40_49 <- data.frame(t(new_df_filt_bacteria_t_40_49_t))
new_df_filt_bacteria_t_40_49$sum_row <- rowSums(new_df_filt_bacteria_t_40_49) 
new_df_filt_bacteria_t_40_49$sum_row_2 <- new_df_filt_bacteria_t_40_49$sum_row / sum(new_df_filt_bacteria_t_40_49$sum_row) * 100
new_df_filt_bacteria_t_40_49 <- new_df_filt_bacteria_t_40_49[order(-new_df_filt_bacteria_t_40_49$sum_row_2),]
new_df_filt_bacteria_t_40_49$cumsum_col <- cumsum(new_df_filt_bacteria_t_40_49$sum_row_2)
count_bacteria_filt_40_49_abundant <- new_df_filt_bacteria_t_40_49[new_df_filt_bacteria_t_40_49$cumsum_col<95,]
count_bacteria_filt_40_49_abundant$Species <- rownames(count_bacteria_filt_40_49_abundant)
count_bacteria_filt_40_49_abundant$cumsum_col <- NULL
count_bacteria_filt_40_49_abundant$sum_row <- NULL
count_bacteria_filt_40_49_abundant$sum_row_2 <- NULL
count_bacteria_filt_40_49_rare <- new_df_filt_bacteria_t_40_49[new_df_filt_bacteria_t_40_49$cumsum_col>=95,]
count_bacteria_filt_40_49_rare$Species <- rownames(count_bacteria_filt_40_49_rare)
count_bacteria_filt_40_49_rare$cumsum_col <- NULL
count_bacteria_filt_40_49_rare$sum_row <- NULL
count_bacteria_filt_40_49_rare$sum_row_2 <- NULL

# species abundance filtering, 50_60
new_df_filt_bacteria_t_50_60_t <- subset(new_df_filt_bacteria_t, age_group == "50_60")
new_df_filt_bacteria_t_50_60_t$age_group <- NULL
new_df_filt_bacteria_t_50_60 <- data.frame(t(new_df_filt_bacteria_t_50_60_t))
new_df_filt_bacteria_t_50_60$sum_row <- rowSums(new_df_filt_bacteria_t_50_60) 
new_df_filt_bacteria_t_50_60$sum_row_2 <- new_df_filt_bacteria_t_50_60$sum_row / sum(new_df_filt_bacteria_t_50_60$sum_row) * 100
new_df_filt_bacteria_t_50_60 <- new_df_filt_bacteria_t_50_60[order(-new_df_filt_bacteria_t_50_60$sum_row_2),]
new_df_filt_bacteria_t_50_60$cumsum_col <- cumsum(new_df_filt_bacteria_t_50_60$sum_row_2)
count_bacteria_filt_50_60_abundant <- new_df_filt_bacteria_t_50_60[new_df_filt_bacteria_t_50_60$cumsum_col<95,]
count_bacteria_filt_50_60_abundant$Species <- rownames(count_bacteria_filt_50_60_abundant)
count_bacteria_filt_50_60_abundant$cumsum_col <- NULL
count_bacteria_filt_50_60_abundant$sum_row <- NULL
count_bacteria_filt_50_60_abundant$sum_row_2 <- NULL
count_bacteria_filt_50_60_rare <- new_df_filt_bacteria_t_50_60[new_df_filt_bacteria_t_50_60$cumsum_col>=95,]
count_bacteria_filt_50_60_rare$Species <- rownames(count_bacteria_filt_50_60_rare)
count_bacteria_filt_50_60_rare$cumsum_col <- NULL
count_bacteria_filt_50_60_rare$sum_row <- NULL
count_bacteria_filt_50_60_rare$sum_row_2 <- NULL

# re-merge all data-frames with abundant species
df_abundant_list <- list(count_bacteria_filt_10_19_abundant, count_bacteria_filt_20_29_abundant,
                         count_bacteria_filt_30_39_abundant, 
                count_bacteria_filt_40_49_abundant, count_bacteria_filt_50_60_abundant)  

new_df_abundant_bac = df_abundant_list %>% reduce(full_join, by='Species')
rownames(new_df_abundant_bac) <- new_df_abundant_bac$Species
new_df_abundant_bac$Species <- NULL
new_df_abundant_bac[is.na(new_df_abundant_bac)] <- 0

# re-merge all data-frames with rare species
df_rare_list <- list(count_bacteria_filt_10_19_rare, count_bacteria_filt_20_29_rare,
                         count_bacteria_filt_30_39_rare, 
                         count_bacteria_filt_40_49_rare, count_bacteria_filt_50_60_rare)  

new_df_rare_bac = df_rare_list %>% reduce(full_join, by='Species')
rownames(new_df_rare_bac) <- new_df_rare_bac$Species
new_df_rare_bac$Species <- NULL
new_df_rare_bac[is.na(new_df_rare_bac)] <- 0


# species abundance filtering (Viruses) #################
new_df_filt_phages_viruses_t <- data.frame(t(new_df_filt_phages_viruses))
new_df_filt_phages_viruses_t$age_group <- meta_data_sp_noBlank$age_group

# species abundance filtering, 10_19
new_df_filt_phages_t_10_19_t <- subset(new_df_filt_phages_viruses_t, age_group == "10_19")
new_df_filt_phages_t_10_19_t$age_group <- NULL
new_df_filt_phages_t_10_19 <- data.frame(t(new_df_filt_phages_t_10_19_t))
new_df_filt_phages_t_10_19$sum_row <- rowSums(new_df_filt_phages_t_10_19) 
new_df_filt_phages_t_10_19$sum_row_2 <- new_df_filt_phages_t_10_19$sum_row / sum(new_df_filt_phages_t_10_19$sum_row) * 100
new_df_filt_phages_t_10_19 <- new_df_filt_phages_t_10_19[order(-new_df_filt_phages_t_10_19$sum_row_2),]
new_df_filt_phages_t_10_19$cumsum_col <- cumsum(new_df_filt_phages_t_10_19$sum_row_2)
count_phages_filt_10_19_abundant <- new_df_filt_phages_t_10_19[new_df_filt_phages_t_10_19$cumsum_col<95,]
count_phages_filt_10_19_abundant$Species<- rownames(count_phages_filt_10_19_abundant)
count_phages_filt_10_19_abundant$cumsum_col <- NULL
count_phages_filt_10_19_abundant$sum_row <- NULL
count_phages_filt_10_19_abundant$sum_row_2 <- NULL
count_phages_filt_10_19_rare <- new_df_filt_phages_t_10_19[new_df_filt_phages_t_10_19$cumsum_col>=95,]
count_phages_filt_10_19_rare$Species<- rownames(count_phages_filt_10_19_rare)
count_phages_filt_10_19_rare$cumsum_col <- NULL
count_phages_filt_10_19_rare$sum_row <- NULL
count_phages_filt_10_19_rare$sum_row_2 <- NULL

# species abundance filtering, 20_29
new_df_filt_phages_t_20_29_t <- subset(new_df_filt_phages_viruses_t, age_group == "20_29")
new_df_filt_phages_t_20_29_t$age_group <- NULL
new_df_filt_phages_t_20_29 <- data.frame(t(new_df_filt_phages_t_20_29_t))
new_df_filt_phages_t_20_29$sum_row <- rowSums(new_df_filt_phages_t_20_29) 
new_df_filt_phages_t_20_29$sum_row_2 <- new_df_filt_phages_t_20_29$sum_row / sum(new_df_filt_phages_t_20_29$sum_row) * 100
new_df_filt_phages_t_20_29 <- new_df_filt_phages_t_20_29[order(-new_df_filt_phages_t_20_29$sum_row_2),]
new_df_filt_phages_t_20_29$cumsum_col <- cumsum(new_df_filt_phages_t_20_29$sum_row_2)
count_phages_filt_20_29_abundant <- new_df_filt_phages_t_20_29[new_df_filt_phages_t_20_29$cumsum_col<95,]
count_phages_filt_20_29_abundant$Species<- rownames(count_phages_filt_20_29_abundant)
count_phages_filt_20_29_abundant$cumsum_col <- NULL
count_phages_filt_20_29_abundant$sum_row <- NULL
count_phages_filt_20_29_abundant$sum_row_2 <- NULL
count_phages_filt_20_29_rare <- new_df_filt_phages_t_20_29[new_df_filt_phages_t_20_29$cumsum_col>=95,]
count_phages_filt_20_29_rare$Species<- rownames(count_phages_filt_20_29_rare)
count_phages_filt_20_29_rare$cumsum_col <- NULL
count_phages_filt_20_29_rare$sum_row <- NULL
count_phages_filt_20_29_rare$sum_row_2 <- NULL


# species abundance filtering, 30_39
new_df_filt_phages_t_30_39_t <- subset(new_df_filt_phages_viruses_t, age_group == "30_39")
new_df_filt_phages_t_30_39_t$age_group <- NULL
new_df_filt_phages_t_30_39 <- data.frame(t(new_df_filt_phages_t_30_39_t))
new_df_filt_phages_t_30_39$sum_row <- rowSums(new_df_filt_phages_t_30_39) 
new_df_filt_phages_t_30_39$sum_row_2 <- new_df_filt_phages_t_30_39$sum_row / sum(new_df_filt_phages_t_30_39$sum_row) * 100
new_df_filt_phages_t_30_39 <- new_df_filt_phages_t_30_39[order(-new_df_filt_phages_t_30_39$sum_row_2),]
new_df_filt_phages_t_30_39$cumsum_col <- cumsum(new_df_filt_phages_t_30_39$sum_row_2)
count_phages_filt_30_39_abundant <- new_df_filt_phages_t_30_39[new_df_filt_phages_t_30_39$cumsum_col<95,]
count_phages_filt_30_39_abundant$Species<- rownames(count_phages_filt_30_39_abundant)
count_phages_filt_30_39_abundant$cumsum_col <- NULL
count_phages_filt_30_39_abundant$sum_row <- NULL
count_phages_filt_30_39_abundant$sum_row_2 <- NULL
count_phages_filt_30_39_rare <- new_df_filt_phages_t_30_39[new_df_filt_phages_t_30_39$cumsum_col>=95,]
count_phages_filt_30_39_rare$Species<- rownames(count_phages_filt_30_39_rare)
count_phages_filt_30_39_rare$cumsum_col <- NULL
count_phages_filt_30_39_rare$sum_row <- NULL
count_phages_filt_30_39_rare$sum_row_2 <- NULL

# species abundance filtering, 40_49
new_df_filt_phages_t_40_49_t <- subset(new_df_filt_phages_viruses_t, age_group == "40_49")
new_df_filt_phages_t_40_49_t$age_group <- NULL
new_df_filt_phages_t_40_49 <- data.frame(t(new_df_filt_phages_t_40_49_t))
new_df_filt_phages_t_40_49$sum_row <- rowSums(new_df_filt_phages_t_40_49) 
new_df_filt_phages_t_40_49$sum_row_2 <- new_df_filt_phages_t_40_49$sum_row / sum(new_df_filt_phages_t_40_49$sum_row) * 100
new_df_filt_phages_t_40_49 <- new_df_filt_phages_t_40_49[order(-new_df_filt_phages_t_40_49$sum_row_2),]
new_df_filt_phages_t_40_49$cumsum_col <- cumsum(new_df_filt_phages_t_40_49$sum_row_2)
count_phages_filt_40_49_abundant <- new_df_filt_phages_t_40_49[new_df_filt_phages_t_40_49$cumsum_col<95,]
count_phages_filt_40_49_abundant$Species <- rownames(count_phages_filt_40_49_abundant)
count_phages_filt_40_49_abundant$cumsum_col <- NULL
count_phages_filt_40_49_abundant$sum_row <- NULL
count_phages_filt_40_49_abundant$sum_row_2 <- NULL
count_phages_filt_40_49_rare <- new_df_filt_phages_t_40_49[new_df_filt_phages_t_40_49$cumsum_col>=95,]
count_phages_filt_40_49_rare$Species<- rownames(count_phages_filt_40_49_rare)
count_phages_filt_40_49_rare$cumsum_col <- NULL
count_phages_filt_40_49_rare$sum_row <- NULL
count_phages_filt_40_49_rare$sum_row_2 <- NULL

# species abundance filtering, 50_60
new_df_filt_phages_t_50_60_t <- subset(new_df_filt_phages_viruses_t, age_group == "50_60")
new_df_filt_phages_t_50_60_t$age_group <- NULL
new_df_filt_phages_t_50_60 <- data.frame(t(new_df_filt_phages_t_50_60_t))
new_df_filt_phages_t_50_60$sum_row <- rowSums(new_df_filt_phages_t_50_60) 
new_df_filt_phages_t_50_60$sum_row_2 <- new_df_filt_phages_t_50_60$sum_row / sum(new_df_filt_phages_t_50_60$sum_row) * 100
new_df_filt_phages_t_50_60 <- new_df_filt_phages_t_50_60[order(-new_df_filt_phages_t_50_60$sum_row_2),]
new_df_filt_phages_t_50_60$cumsum_col <- cumsum(new_df_filt_phages_t_50_60$sum_row_2)
count_phages_filt_50_60_abundant <- new_df_filt_phages_t_50_60[new_df_filt_phages_t_50_60$cumsum_col<95,]
count_phages_filt_50_60_abundant$Species <- rownames(count_phages_filt_50_60_abundant)
count_phages_filt_50_60_abundant$cumsum_col <- NULL
count_phages_filt_50_60_abundant$sum_row <- NULL
count_phages_filt_50_60_abundant$sum_row_2 <- NULL
count_phages_filt_50_60_rare <- new_df_filt_phages_t_50_60[new_df_filt_phages_t_50_60$cumsum_col>=95,]
count_phages_filt_50_60_rare$Species<- rownames(count_phages_filt_50_60_rare)
count_phages_filt_50_60_rare$cumsum_col <- NULL
count_phages_filt_50_60_rare$sum_row <- NULL
count_phages_filt_50_60_rare$sum_row_2 <- NULL

# re-merge all data-frames with abundant species
df_abundant_phages_list <- list(count_phages_filt_10_19_abundant, count_phages_filt_20_29_abundant,
                         count_phages_filt_30_39_abundant, 
                         count_phages_filt_40_49_abundant, count_phages_filt_50_60_abundant)  

new_df_abundant_pha = df_abundant_phages_list %>% reduce(full_join, by='Species')
rownames(new_df_abundant_pha) <- new_df_abundant_pha$Species
new_df_abundant_pha$Species <- NULL
new_df_abundant_pha[is.na(new_df_abundant_pha)] <- 0

# re-merge all data-frames with rare species
df_rare_phageslist <- list(count_phages_filt_10_19_rare, count_phages_filt_20_29_rare,
                     count_phages_filt_30_39_rare, 
                     count_phages_filt_40_49_rare, count_phages_filt_50_60_rare)  

new_df_rare_pha = df_rare_phageslist %>% reduce(full_join, by='Species')
rownames(new_df_rare_pha) <- new_df_rare_pha$Species
new_df_rare_pha$Species <- NULL
new_df_rare_pha[is.na(new_df_rare_pha)] <- 0

############################################################################################
# Diversity estimates
data_bac <- data.frame(rbind(new_df_rare_bac, new_df_abundant_bac))
data_bac_t <- data.frame(t(data_bac))
data_bac_t <- data_bac_t[order(rownames(data_bac_t)),]
data_bac_rare <- data.frame(t(new_df_rare_bac))
data_bac_rare <- data_bac_rare[order(rownames(data_bac_rare)),]
data_bac_core <- data.frame(t(new_df_abundant_bac))
data_bac_core <- data_bac_core[order(rownames(data_bac_core)),]

# total bacteria
meta_data_sp_noBlank$shannon_bac <- vegan::diversity(data_bac_t, index = "shannon")
meta_data_sp_noBlank$specnumber_bac <- vegan::specnumber(data_bac_t)
meta_data_sp_noBlank$evenness_bac <- meta_data_sp_noBlank$shannon_bac / log(meta_data_sp_noBlank$specnumber_bac)
meta_data_sp_noBlank$specnumber_bac <- NULL

# low-abundance bacteria
meta_data_sp_noBlank$shannon_rare_bac <- vegan::diversity(data_bac_rare, index="shannon")
meta_data_sp_noBlank$simpson_rare_bac <- vegan::diversity(data_bac_rare, index="simpson")
meta_data_sp_noBlank$specnumber_rare_bac <- vegan::specnumber(data_bac_rare)
meta_data_sp_noBlank$evenness_rare_bac <- meta_data_sp_noBlank$shannon_rare_bac / log(meta_data_sp_noBlank$specnumber_rare_bac)
meta_data_sp_noBlank$richness_rare <- meta_data_sp_noBlank$specnumber_rare_bac / sqrt(length(meta_data_sp_noBlank$specnumber_rare_bac))

# high-abundance bacteria
meta_data_sp_noBlank$shannon_core_bac <- vegan::diversity(data_bac_core, index="shannon")
meta_data_sp_noBlank$simpson_core_bac <- vegan::diversity(data_bac_core, index="simpson")
meta_data_sp_noBlank$specnumber_core_bac <- vegan::specnumber(data_bac_core)
meta_data_sp_noBlank$evenness_core_bac <- meta_data_sp_noBlank$shannon_core_bac / log(meta_data_sp_noBlank$specnumber_core_bac)
meta_data_sp_noBlank$richness_core <- meta_data_sp_noBlank$specnumber_core_bac / sqrt(length(meta_data_sp_noBlank$specnumber_core_bac))
data_phage <- data.frame(rbind(new_df_rare_pha, new_df_abundant_pha))
data_phage_t <- data.frame(t(data_phage))
data_phage_t <- data_phage_t[order(rownames(data_phage_t)),]

# phage diversity (total)
meta_data_sp_noBlank$shannon_div_phage <- vegan::diversity(data_phage_t, index = "shannon")
meta_data_sp_noBlank$specnumber_phage <- vegan::specnumber(data_phage_t)
meta_data_sp_noBlank$evenness_phage <- meta_data_sp_noBlank$shannon_div_phage / log(meta_data_sp_noBlank$specnumber_phage)
meta_data_sp_noBlank$specnumber_phage <- NULL
meta_data_sp_noBlank$age_group_2 <- str_replace_all(meta_data_sp_noBlank$age_group, "_", "-")

# make Venn Diagrams 
# Age group (10-19), core
data_bac_core$Age_group <- meta_data_sp_noBlank$age_group_2
data_bac_core_group1 <- data_bac_core[data_bac_core$Age_group=="10-19",]
data_bac_core_group1$Age_group <- NULL
data_bac_core_group1_t <- data.frame(t(data_bac_core_group1))
data_bac_core_group1_t
data_bac_core_group1_t$med <- rowMedians(as.matrix(data_bac_core_group1_t))
data_bac_core_group1_t <- data_bac_core_group1_t[data_bac_core_group1_t$med>0,]
list_core_spec_10_19 <- rownames(data_bac_core_group1_t)

# Age group (10-19), rare
data_bac_rare$Age_group <- meta_data_sp_noBlank$age_group_2
data_bac_rare_group1 <- data_bac_rare[data_bac_rare$Age_group=="10-19",]
data_bac_rare_group1$Age_group <- NULL
data_bac_rare_group1_t <- data.frame(t(data_bac_rare_group1))
data_bac_rare_group1_t$med <- rowMedians(as.matrix(data_bac_rare_group1_t))
data_bac_rare_group1_t <- data_bac_rare_group1_t[data_bac_rare_group1_t$med>0,]
list_rare_spec_10_19 <- rownames(data_bac_rare_group1_t)

# Age group (20-29), core
data_bac_core_group2 <- data_bac_core[data_bac_core$Age_group=="20-29",]
data_bac_core_group2$Age_group <- NULL
data_bac_core_group2_t <- data.frame(t(data_bac_core_group2))
data_bac_core_group2_t$med <- rowMedians(as.matrix(data_bac_core_group2_t))
data_bac_core_group2_t <- data_bac_core_group2_t[data_bac_core_group2_t$med>0,]
list_core_spec_20_29 <- rownames(data_bac_core_group2_t)

# Age group (20-29), rare
data_bac_rare_group2 <- data_bac_rare[data_bac_rare$Age_group=="20-29",]
data_bac_rare_group2$Age_group <- NULL
data_bac_rare_group2_t <- data.frame(t(data_bac_rare_group2))
data_bac_rare_group2_t$med <- rowMedians(as.matrix(data_bac_rare_group2_t))
data_bac_rare_group2_t <- data_bac_rare_group2_t[data_bac_rare_group2_t$med>0,]
list_rare_spec_20_29 <- rownames(data_bac_rare_group2_t)

# Age group >30, core
data_bac_core_group3 <- data_bac_core[data_bac_core$Age_group=="30-39" | 
                                        data_bac_core$Age_group=="40-49" | 
                                        data_bac_core$Age_group=="50-60",]
data_bac_core_group3$Age_group <- NULL
data_bac_core_group3_t <- data.frame(t(data_bac_core_group3))
data_bac_core_group3_t$med <- rowMedians(as.matrix(data_bac_core_group3_t))
data_bac_core_group3_t <- data_bac_core_group3_t[data_bac_core_group3_t$med>0,]
list_core_spec_30_above <- rownames(data_bac_core_group3_t)

# Age group >30, rare
data_bac_rare_group3 <- data_bac_rare[data_bac_rare$Age_group=="30-39"| 
                                        data_bac_rare$Age_group=="40-49" | 
                                        data_bac_rare$Age_group=="50-60",]
data_bac_rare_group3$Age_group <- NULL
data_bac_rare_group3_t <- data.frame(t(data_bac_rare_group3))
data_bac_rare_group3_t$med <- rowMedians(as.matrix(data_bac_rare_group3_t))
data_bac_rare_group3_t <- data_bac_rare_group3_t[data_bac_rare_group3_t$med>0,]
list_rare_spec_30_above <- rownames(data_bac_rare_group3_t)


rare_age_venn <- list(A = list_rare_spec_10_19, B = list_rare_spec_20_29, C = list_rare_spec_30_above)
core_age_venn <- list(A = list_core_spec_10_19, B = list_core_spec_20_29, C = list_core_spec_30_above)

gg_venn_rare <-
  ggVennDiagram(rare_age_venn, size=1.5, 
              category.names = c('10-19', '20-29', '>30'), label = 'both', label_alpha=0, percent_digits = 0) +
  theme_bw(base_size = 6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_text(color="black", size=10), 
        legend.position = "right",
        plot.margin = unit(c(-0.6,-0.3,-0.3,-0.3), "cm")) +
   scale_fill_gradient(name="Species count", low="beige", high = "darkred", limits=c(0,40)) + 
   scale_colour_manual(values=c('seashell2', 'seashell2', 'seashell2')) 

gg_venn_core <-
  ggVennDiagram(core_age_venn, size=1.5, 
                category.names = c('10-19', '20-29', '>30'), label = 'both', label_alpha=0, percent_digits = 0) +
  theme_bw(base_size = 6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_text(color="black", size=10), 
        legend.position = "none",
        plot.margin = unit(c(-0.3,-0.3,-0.3,-0.3), "cm")) +
  scale_fill_gradient(name="Species count", low="beige", high = "darkred", limits=c(0,40)) + 
  scale_colour_manual(values=c('seashell2', 'seashell2', 'seashell2')) 

venn_plots_merged <- ggarrange(gg_venn_core, gg_venn_rare, labels=c("A", "B"), common.legend = TRUE, legend = "right")

pdf("venn_plots_merged.pdf", width=10, height=7)
venn_plots_merged
dev.off()


data_bac_rare$Age_group <- NULL
data_bac_core$Others <- rowSums(data_bac_rare)
data_bac_core$Samples <- rownames(data_bac_core)

# randomy sub-samble with 10 samples per age group
set.seed(1)
data_bac_core_sample <- data_bac_core %>% group_by(Age_group) %>% sample_n(10)
data_bac_core_sample <- data.frame(data_bac_core_sample)

data_bac_core_rel_plot <-gather(data_bac_core_sample, key="Species", value="Relative_abundance", -c(Samples, Age_group))

data_bac_core_rel_plot$Species <- str_replace(data_bac_core_rel_plot$Species, "\\.", " ")
data_bac_core_rel_plot$Species <- str_replace(data_bac_core_rel_plot$Species, " \\.", " ")

cols=c("Actinomyces meyeri"="cyan2",
       "Atopobium parvulum"="grey20",
       "Campylobacter concisus"="grey",
       "Capnocytophaga gingivalis"="blue",    
       "Capnocytophaga leadbetteri"="blue",
       "Eubacterium sulci"="green",
       "Fusobacterium hwasookii"="orange",
       "Fusobacterium periodonticum"="orange",   
       "Gemella haemolysans"="lightgreen",
       "Gemella morbillorum"="lightgreen",
       "Gemella sanguinis"="lightgreen",
       "Haemophilus parainfluenzae"="firebrick2",    
       "Leptotrichia buccalis"="firebrick",
       "Neisseria elongata"="khaki1",
       "Neisseria meningitidis"="khaki1",
       "Neisseria mucosa"="khaki1",   
       "Others"="white",
       "Porphyromonas gingivalis"="gold",
       "Prevotella denticola"="deeppink4",
       "Prevotella enoeca"="deeppink4",
       "Prevotella fusca"="deeppink4",              
       "Prevotella intermedia"="deeppink4",
       "Prevotella jejuni"="deeppink4",
       "Prevotella melaninogenica"="deeppink4",
       "Rothia aeria"="black",                  
       "Rothia dentocariosa"="black",
       "Rothia mucilaginosa"="black",
       "Schaalia odontolytica"="cyan2",
       "Streptococcus constellatus"="dodgerblue4",    
       "Streptococcus equinus"="dodgerblue4",
       "Streptococcus gordonii"="dodgerblue4",
       "Streptococcus mitis"="dodgerblue4",
       "Streptococcus oralis"="dodgerblue4",          
       "Streptococcus parasanguinis"="dodgerblue4",
       "Streptococcus pneumoniae"="dodgerblue4",
       "Streptococcus pseudopneumoniae"="dodgerblue4",
       "Streptococcus salivarius"="dodgerblue4",      
       "Streptococcus sanguinis"="dodgerblue4",
       "Tannerella forsythia"="forestgreen",
       "Veillonella atypica"="forestgreen",
       "Veillonella dispar"="forestgreen",           
       "Veillonella parvula"="forestgreen")

rel_abund <-
  ggplot(data_bac_core_rel_plot, aes(x=Samples, y=Relative_abundance, fill=Species)) +
  geom_bar(position = "fill", stat = "identity", colour="black",width = 0.9) +
  facet_wrap(~Age_group, nrow=1, scales = "free_x") +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=4)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(face="italic"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill=NA)) + 
  xlab("Healthy airways") + 
  ylab("Relative abundance")

core_shannon_plot <-
  ggplot(meta_data_sp_noBlank, aes(y=shannon_core_bac, x=age_group_2)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width=0.05) + 
  theme_bw() +
  scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4.5)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank()) +
  xlab(" ") + ylab("Shannon diversity")

core_simpson_plot <-
  ggplot(meta_data_sp_noBlank, aes(y=simpson_core_bac, x=age_group_2)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width=0.05) + 
  theme_bw() +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1.3)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank()) +
  xlab(" ") + ylab("Simpson diversity index")

core_evenness_plot <-
  ggplot(meta_data_sp_noBlank, aes(y=evenness_core_bac, x=age_group_2)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width=0.05) + 
  theme_bw() +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1.5)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=90)) +
  xlab("Age [years]") + ylab("Pielou's evenness index")

core_richness_plot <-
  ggplot(meta_data_sp_noBlank, aes(y=richness_core, x=age_group_2)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width=0.05) + 
  theme_bw() +
  scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4.5)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank()) +
  xlab(" ") + ylab("Menhinick's index")

div_plots <- ggarrange(core_shannon_plot, core_simpson_plot, core_richness_plot, core_evenness_plot, nrow = 4, labels = c("B", "C", "D", "E"))
core_plots <- ggarrange(rel_abund, div_plots, nrow=1, widths = c(1,0.3), labels = c("A", "B"))

pdf("core_plots.pdf", width = 10, height = 11)
core_plots
dev.off()


rare_shannon_plot <-
  ggplot(meta_data_sp_noBlank, aes(y=shannon_rare_bac, x=age_group_2)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width=0.05) + 
  theme_bw() +
  scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4.5)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank()) +
  xlab(" ") + ylab("Shannon diversity")

rare_simpson_plot <-
  ggplot(meta_data_sp_noBlank, aes(y=simpson_rare_bac, x=age_group_2)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width=0.05) + 
  theme_bw() +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1.3)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank()) +
  xlab(" ") + ylab("Simpson diversity index")

rare_evenness_plot <-
  ggplot(meta_data_sp_noBlank, aes(y=evenness_rare_bac, x=age_group_2)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width=0.05) + 
  theme_bw() +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1.5)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=90)) +
  xlab("Age [years]") + ylab("Pielou's evenness index")

compare_richness_rare <- list(c("10-19", "20-29"), c("20-29", "50-60"))

rare_richness_plot <-
  ggplot(meta_data_sp_noBlank, aes(y=richness_rare, x=age_group_2)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width=0.05) + 
  theme_bw() +
  stat_compare_means(comparisons = compare_richness_rare, label = "p.signif", label.y = c(4.8, 5.5)) +
  scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,6)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=90)) +
  xlab(" ") + ylab("Menhinick's index")

#kruskal.test(meta_data_sp_noBlank$richness_rare, g=meta_data_sp_noBlank$age_group_2)
#rcompanion::epsilonSquared(meta_data_sp_noBlank$richness_rare, g=meta_data_sp_noBlank$age_group_2, ci=TRUE)

rare_correlation_plot <-
  ggplot(meta_data_sp_noBlank, aes(y=richness_rare, x=age_years)) +
  geom_point() +
  geom_smooth(method = "loess", colour="darkred", fill="khaki", size=0.8) +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,5)) +
  theme(panel.grid = element_blank()) +
  xlab("Age [years]") + ylab("Menhinick's index")

div_plots_rare <- ggarrange(rare_shannon_plot, rare_simpson_plot, rare_richness_plot, rare_evenness_plot, 
                            nrow=2,ncol=2, widths = c(0.95,1), heights = c(0.87,1), labels = c("B", "C", "D", "E"))
rare_plots <- ggarrange(rare_correlation_plot, div_plots_rare, nrow=1, widths = c(1,1), labels = c("A", "B"))

pdf("rare_plots.pdf", width = 8, height = 6)
rare_plots
dev.off()


#############################################################################################
# Spearman's rank correlation (10_19)
# merge high and low abundant taxa
count_bacteria_filt_10_19_rare$Species <- NULL
list_bac_10_19_rare <- rownames(count_bacteria_filt_10_19_rare)
count_bacteria_filt_10_19_abundant$Species <- NULL
list_bac_10_19_abundant <- rownames(count_bacteria_filt_10_19_abundant)
count_bacteria_10_19_all <- data.frame(rbind(count_bacteria_filt_10_19_abundant, count_bacteria_filt_10_19_rare))
count_bacteria_10_19_all_clr <- clr(count_bacteria_10_19_all)
length(colnames(count_bacteria_10_19_all_clr))
# perform Spearman's rank correlation analysis 
list_10_19_bac_mean_abundance <- rowMedians(as.matrix(count_bacteria_10_19_all_clr))
names(list_10_19_bac_mean_abundance) <- rownames(count_bacteria_10_19_all_clr)
spearman_norm_10_19_bac <- rcorr(as.matrix(t(count_bacteria_10_19_all_clr)), type = 'spearman')
# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_norm_10_19_bac <- spearman_norm_10_19_bac$P
# extract and store correlation coefficients of analysis
spearman_r_norm_10_19_bac <- spearman_norm_10_19_bac$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
spearman_10_19_bac_edges <- reshape2::melt(spearman_p_norm_10_19_bac)
# take data in wide format and stack columns into a single column of data for correlation coefficients
spearman_10_19_bac_cor_edges <- reshape2::melt(spearman_r_norm_10_19_bac)
# merge tables
spearman_10_19_bac_edges$COR <- spearman_10_19_bac_cor_edges$value
# round p-values to five decimal places
spearman_10_19_bac_edges$value <- round(spearman_10_19_bac_edges$value, 5)
# store row names in own column
spearman_10_19_bac_edges$Label <- row.names(spearman_10_19_bac_cor_edges)
# rename columns
colnames(spearman_10_19_bac_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_10_19_bac_edges) <- NULL
# remove NAs
spearman_10_19_bac_edges <- spearman_10_19_bac_edges[complete.cases(spearman_10_19_bac_edges), ]
# extract significant correlations
spearman_10_19_bac_edges_sig <- spearman_10_19_bac_edges[spearman_10_19_bac_edges$pValue < p_corr,]
spearman_10_19_bac_edges_sig <- subset(spearman_10_19_bac_edges_sig, Weight <= -0.2 | Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_10_19_bac_edges_sig$Correlation <- ifelse(spearman_10_19_bac_edges_sig$Weight > 0, "pos", "neg")
spearman_10_19_bac_edges_sig$Type <- "Directed"
spearman_10_19_bac_edges_sig_final <- select(spearman_10_19_bac_edges_sig, c("Source", "Target", "Type", "Correlation"))
spearman_10_19_bac_edges_sig_final$Genus <- sapply(strsplit(as.character(spearman_10_19_bac_edges_sig_final$Target)," "), `[`, 1)
#select both columns to create node list
spearman_10_19_bac_nodes_sig <- select(spearman_10_19_bac_edges_sig, c("Target"))
spearman_10_19_bac_nodes_sig$Target2 <- spearman_10_19_bac_nodes_sig$Target
spearman_10_19_bac_nodes_sig$Correlation <- spearman_10_19_bac_edges_sig$Correlation

# remove duplicate entries
spearman_10_19_bac_nodes_sig = spearman_10_19_bac_nodes_sig[!duplicated(spearman_10_19_bac_nodes_sig$Target),]
# make data frame
spearman_10_19_bac_nodes_sig <- data.frame(spearman_10_19_bac_nodes_sig)
# re-index rows
rownames(spearman_10_19_bac_nodes_sig) <- NULL
# rename columns
colnames(spearman_10_19_bac_nodes_sig) <- c("Id", "Label", "Correlation")
# add genus information
spearman_10_19_bac_nodes_sig$Genus <- sapply(strsplit(as.character(spearman_10_19_bac_nodes_sig$Label)," "), `[`, 1)
mean_list_10_19_bac = c()
for (items in spearman_10_19_bac_nodes_sig$Id){
  mean_list_10_19_bac = append(mean_list_10_19_bac, list_10_19_bac_mean_abundance[items])}
spearman_10_19_bac_nodes_sig$mean_abundance <- mean_list_10_19_bac

# export node list
write.table(spearman_10_19_bac_nodes_sig, file="network_tables/nodes_10_19_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)
# export edge list
write.table(spearman_10_19_bac_edges_sig_final, file="network_tables/edges_10_19_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)


############################################################################################
# Spearman's rank correlation (20_29)
# merge high and low abundant taxa
count_bacteria_filt_20_29_rare$Species <- NULL
list_bac_20_29_rare <- rownames(count_bacteria_filt_20_29_rare)
count_bacteria_filt_20_29_abundant$Species <- NULL
list_bac_20_29_abundant <- rownames(count_bacteria_filt_20_29_abundant)
count_bacteria_20_29_all <- data.frame(rbind(count_bacteria_filt_20_29_abundant, count_bacteria_filt_20_29_rare))
count_bacteria_20_29_all_clr <- clr(count_bacteria_20_29_all)


# perform Spearman's rank correlation analysis 
list_20_29_bac_mean_abundance <- rowMedians(as.matrix(count_bacteria_20_29_all_clr))
names(list_20_29_bac_mean_abundance) <- rownames(count_bacteria_20_29_all_clr)
spearman_norm_20_29_bac <- rcorr(as.matrix(t(count_bacteria_20_29_all_clr)), type = 'spearman')
# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_norm_20_29_bac <- spearman_norm_20_29_bac$P
# extract and store correlation coefficients of analysis
spearman_r_norm_20_29_bac <- spearman_norm_20_29_bac$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
spearman_20_29_bac_edges <- reshape2::melt(spearman_p_norm_20_29_bac)
# take data in wide format and stack columns into a single column of data for correlation coefficients
spearman_20_29_bac_cor_edges <- reshape2::melt(spearman_r_norm_20_29_bac)
# merge tables
spearman_20_29_bac_edges$COR <- spearman_20_29_bac_cor_edges$value
# round p-values to five decimal places
spearman_20_29_bac_edges$value <- round(spearman_20_29_bac_edges$value, 5)
# store row names in own column
spearman_20_29_bac_edges$Label <- row.names(spearman_20_29_bac_cor_edges)
# rename columns
colnames(spearman_20_29_bac_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_20_29_bac_edges) <- NULL
# remove NAs
spearman_20_29_bac_edges <- spearman_20_29_bac_edges[complete.cases(spearman_20_29_bac_edges), ]
# extract significant correlations
spearman_20_29_bac_edges_sig <- spearman_20_29_bac_edges[spearman_20_29_bac_edges$pValue < p_corr,]
spearman_20_29_bac_edges_sig <- subset(spearman_20_29_bac_edges_sig, Weight <= -0.2 | Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_20_29_bac_edges_sig$Correlation <- ifelse(spearman_20_29_bac_edges_sig$Weight > 0, "pos", "neg")
spearman_20_29_bac_edges_sig$Type <- "Directed"
spearman_20_29_bac_edges_sig_final <- select(spearman_20_29_bac_edges_sig, c("Source", "Target", "Type", "Correlation"))
spearman_20_29_bac_edges_sig_final$Genus <- sapply(strsplit(as.character(spearman_20_29_bac_edges_sig_final$Target)," "), `[`, 1)
#select both columns to create node list
spearman_20_29_bac_nodes_sig <- select(spearman_20_29_bac_edges_sig, c("Target"))
spearman_20_29_bac_nodes_sig$Target2 <- spearman_20_29_bac_nodes_sig$Target
spearman_20_29_bac_nodes_sig$Correlation <- spearman_20_29_bac_edges_sig$Correlation

# remove duplicate entries
spearman_20_29_bac_nodes_sig = spearman_20_29_bac_nodes_sig[!duplicated(spearman_20_29_bac_nodes_sig$Target),]
# make data frame
spearman_20_29_bac_nodes_sig <- data.frame(spearman_20_29_bac_nodes_sig)
# re-index rows
rownames(spearman_20_29_bac_nodes_sig) <- NULL
# rename columns
colnames(spearman_20_29_bac_nodes_sig) <- c("Id", "Label", "Correlation")
# add genus information
spearman_20_29_bac_nodes_sig$Genus <- sapply(strsplit(as.character(spearman_20_29_bac_nodes_sig$Label)," "), `[`, 1)
mean_list_20_29_bac = c()
for (items in spearman_20_29_bac_nodes_sig$Id){
  mean_list_20_29_bac = append(mean_list_20_29_bac, list_20_29_bac_mean_abundance[items])}
spearman_20_29_bac_nodes_sig$mean_abundance <- mean_list_20_29_bac
# export node list
write.table(spearman_20_29_bac_nodes_sig, file="network_tables/nodes_20_29_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)
# export edge list
write.table(spearman_20_29_bac_edges_sig_final, file="network_tables/edges_20_29_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)



############################################################################################
# Spearman's rank correlation (30_39)
# merge high and low abundant taxa
count_bacteria_filt_30_39_rare$Species <- NULL
list_bac_30_39_rare <- rownames(count_bacteria_filt_30_39_rare)
count_bacteria_filt_30_39_abundant$Species <- NULL
list_bac_30_39_abundant <- rownames(count_bacteria_filt_30_39_abundant)
count_bacteria_30_39_all <- data.frame(rbind(count_bacteria_filt_30_39_abundant, count_bacteria_filt_30_39_rare))
count_bacteria_30_39_all_clr <- clr(count_bacteria_30_39_all)


# perform Spearman's rank correlation analysis 
list_30_39_bac_mean_abundance <- rowMedians(as.matrix(count_bacteria_30_39_all_clr))
names(list_30_39_bac_mean_abundance) <- rownames(count_bacteria_30_39_all_clr)
spearman_norm_30_39_bac <- rcorr(as.matrix(t(count_bacteria_30_39_all_clr)), type = 'spearman')
# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_norm_30_39_bac <- spearman_norm_30_39_bac$P
# extract and store correlation coefficients of analysis
spearman_r_norm_30_39_bac <- spearman_norm_30_39_bac$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
spearman_30_39_bac_edges <- reshape2::melt(spearman_p_norm_30_39_bac)
# take data in wide format and stack columns into a single column of data for correlation coefficients
spearman_30_39_bac_cor_edges <- reshape2::melt(spearman_r_norm_30_39_bac)
# merge tables
spearman_30_39_bac_edges$COR <- spearman_30_39_bac_cor_edges$value
# round p-values to five decimal places
spearman_30_39_bac_edges$value <- round(spearman_30_39_bac_edges$value, 5)
# store row names in own column
spearman_30_39_bac_edges$Label <- row.names(spearman_30_39_bac_cor_edges)
# rename columns
colnames(spearman_30_39_bac_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_30_39_bac_edges) <- NULL
# remove NAs
spearman_30_39_bac_edges <- spearman_30_39_bac_edges[complete.cases(spearman_30_39_bac_edges), ]
# extract significant correlations
spearman_30_39_bac_edges_sig <- spearman_30_39_bac_edges[spearman_30_39_bac_edges$pValue < p_corr,]
spearman_30_39_bac_edges_sig <- subset(spearman_30_39_bac_edges_sig, Weight <= -0.2 | Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_30_39_bac_edges_sig$Correlation <- ifelse(spearman_30_39_bac_edges_sig$Weight > 0, "pos", "neg")
spearman_30_39_bac_edges_sig$Type <- "Directed"
spearman_30_39_bac_edges_sig_final <- select(spearman_30_39_bac_edges_sig, c("Source", "Target", "Type", "Correlation"))
spearman_30_39_bac_edges_sig_final$Genus <- sapply(strsplit(as.character(spearman_30_39_bac_edges_sig_final$Target)," "), `[`, 1)
#select both columns to create node list
spearman_30_39_bac_nodes_sig <- select(spearman_30_39_bac_edges_sig, c("Target"))
spearman_30_39_bac_nodes_sig$Target2 <- spearman_30_39_bac_nodes_sig$Target
spearman_30_39_bac_nodes_sig$Correlation <- spearman_30_39_bac_edges_sig$Correlation

# remove duplicate entries
spearman_30_39_bac_nodes_sig = spearman_30_39_bac_nodes_sig[!duplicated(spearman_30_39_bac_nodes_sig$Target),]
# make data frame
spearman_30_39_bac_nodes_sig <- data.frame(spearman_30_39_bac_nodes_sig)
# re-index rows
rownames(spearman_30_39_bac_nodes_sig) <- NULL
# rename columns
colnames(spearman_30_39_bac_nodes_sig) <- c("Id", "Label", "Correlation")
# add genus information
spearman_30_39_bac_nodes_sig$Genus <- sapply(strsplit(as.character(spearman_30_39_bac_nodes_sig$Label)," "), `[`, 1)
mean_list_30_39_bac = c()
for (items in spearman_30_39_bac_nodes_sig$Id){
  mean_list_30_39_bac = append(mean_list_30_39_bac, list_30_39_bac_mean_abundance[items])}
spearman_30_39_bac_nodes_sig$mean_abundance <- mean_list_30_39_bac
# export node list
write.table(spearman_30_39_bac_nodes_sig, file="network_tables/nodes_30_39_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)
# export edge list
write.table(spearman_30_39_bac_edges_sig_final, file="network_tables/edges_30_39_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)


############################################################################################
# Spearman's rank correlation (40_49)
# merge high and low abundant taxa
count_bacteria_filt_40_49_rare$Species <- NULL
list_bac_40_49_rare <- rownames(count_bacteria_filt_40_49_rare)
count_bacteria_filt_40_49_abundant$Species <- NULL
list_bac_40_49_abundant <- rownames(count_bacteria_filt_40_49_abundant)
count_bacteria_40_49_all <- data.frame(rbind(count_bacteria_filt_40_49_abundant, count_bacteria_filt_40_49_rare))
count_bacteria_40_49_all_clr <- clr(count_bacteria_40_49_all)


# perform Spearman's rank correlation analysis 
list_40_49_bac_mean_abundance <- rowMedians(as.matrix(count_bacteria_40_49_all_clr))
names(list_40_49_bac_mean_abundance) <- rownames(count_bacteria_40_49_all_clr)
spearman_norm_40_49_bac <- rcorr(as.matrix(t(count_bacteria_40_49_all_clr)), type = 'spearman')
# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_norm_40_49_bac <- spearman_norm_40_49_bac$P
# extract and store correlation coefficients of analysis
spearman_r_norm_40_49_bac <- spearman_norm_40_49_bac$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
spearman_40_49_bac_edges <- reshape2::melt(spearman_p_norm_40_49_bac)
# take data in wide format and stack columns into a single column of data for correlation coefficients
spearman_40_49_bac_cor_edges <- reshape2::melt(spearman_r_norm_40_49_bac)
# merge tables
spearman_40_49_bac_edges$COR <- spearman_40_49_bac_cor_edges$value
# round p-values to five decimal places
spearman_40_49_bac_edges$value <- round(spearman_40_49_bac_edges$value, 5)
# store row names in own column
spearman_40_49_bac_edges$Label <- row.names(spearman_40_49_bac_cor_edges)
# rename columns
colnames(spearman_40_49_bac_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_40_49_bac_edges) <- NULL
# remove NAs
spearman_40_49_bac_edges <- spearman_40_49_bac_edges[complete.cases(spearman_40_49_bac_edges), ]
# extract significant correlations
spearman_40_49_bac_edges_sig <- spearman_40_49_bac_edges[spearman_40_49_bac_edges$pValue < p_corr,]
spearman_40_49_bac_edges_sig <- subset(spearman_40_49_bac_edges_sig, Weight <= -0.2 | Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_40_49_bac_edges_sig$Correlation <- ifelse(spearman_40_49_bac_edges_sig$Weight > 0, "pos", "neg")
spearman_40_49_bac_edges_sig$Type <- "Directed"
spearman_40_49_bac_edges_sig_final <- select(spearman_40_49_bac_edges_sig, c("Source", "Target", "Type", "Correlation"))
spearman_40_49_bac_edges_sig_final$Genus <- sapply(strsplit(as.character(spearman_40_49_bac_edges_sig_final$Target)," "), `[`, 1)
#select both columns to create node list
spearman_40_49_bac_nodes_sig <- select(spearman_40_49_bac_edges_sig, c("Target"))
spearman_40_49_bac_nodes_sig$Target2 <- spearman_40_49_bac_nodes_sig$Target
spearman_40_49_bac_nodes_sig$Correlation <- spearman_40_49_bac_edges_sig$Correlation

# remove duplicate entries
spearman_40_49_bac_nodes_sig = spearman_40_49_bac_nodes_sig[!duplicated(spearman_40_49_bac_nodes_sig$Target),]
# make data frame
spearman_40_49_bac_nodes_sig <- data.frame(spearman_40_49_bac_nodes_sig)
# re-index rows
rownames(spearman_40_49_bac_nodes_sig) <- NULL
# rename columns
colnames(spearman_40_49_bac_nodes_sig) <- c("Id", "Label", "Correlation")
# add genus information
spearman_40_49_bac_nodes_sig$Genus <- sapply(strsplit(as.character(spearman_40_49_bac_nodes_sig$Label)," "), `[`, 1)
mean_list_40_49_bac = c()
for (items in spearman_40_49_bac_nodes_sig$Id){
  mean_list_40_49_bac = append(mean_list_40_49_bac, list_40_49_bac_mean_abundance[items])}
spearman_40_49_bac_nodes_sig$mean_abundance <- mean_list_40_49_bac
# export node list
write.table(spearman_40_49_bac_nodes_sig, file="network_tables/nodes_40_49_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)
# export edge list
write.table(spearman_40_49_bac_edges_sig_final, file="network_tables/edges_40_49_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)


############################################################################################
# Spearman's rank correlation (50_60)
# merge high and low abundant taxa
count_bacteria_filt_50_60_rare$Species <- NULL
list_bac_50_60_rare <- rownames(count_bacteria_filt_50_60_rare)
count_bacteria_filt_50_60_abundant$Species <- NULL
list_bac_50_60_abundant <- rownames(count_bacteria_filt_50_60_abundant)
count_bacteria_50_60_all <- data.frame(rbind(count_bacteria_filt_50_60_abundant, count_bacteria_filt_50_60_rare))
count_bacteria_50_60_all_clr <- clr(count_bacteria_50_60_all)
length(colnames(count_bacteria_50_60_all_clr))

# perform Spearman's rank correlation analysis 
list_50_60_bac_mean_abundance <- rowMedians(as.matrix(count_bacteria_50_60_all_clr))
names(list_50_60_bac_mean_abundance) <- rownames(count_bacteria_50_60_all_clr)
spearman_norm_50_60_bac <- rcorr(as.matrix(t(count_bacteria_50_60_all_clr)), type = 'spearman')
# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_norm_50_60_bac <- spearman_norm_50_60_bac$P
# extract and store correlation coefficients of analysis
spearman_r_norm_50_60_bac <- spearman_norm_50_60_bac$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
spearman_50_60_bac_edges <- reshape2::melt(spearman_p_norm_50_60_bac)
# take data in wide format and stack columns into a single column of data for correlation coefficients
spearman_50_60_bac_cor_edges <- reshape2::melt(spearman_r_norm_50_60_bac)
# merge tables
spearman_50_60_bac_edges$COR <- spearman_50_60_bac_cor_edges$value
# round p-values to five decimal places
spearman_50_60_bac_edges$value <- round(spearman_50_60_bac_edges$value, 5)
# store row names in own column
spearman_50_60_bac_edges$Label <- row.names(spearman_50_60_bac_cor_edges)
# rename columns
colnames(spearman_50_60_bac_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_50_60_bac_edges) <- NULL
# remove NAs
spearman_50_60_bac_edges <- spearman_50_60_bac_edges[complete.cases(spearman_50_60_bac_edges), ]
# extract significant correlations
spearman_50_60_bac_edges_sig <- spearman_50_60_bac_edges[spearman_50_60_bac_edges$pValue < p_corr,]
spearman_50_60_bac_edges_sig <- subset(spearman_50_60_bac_edges_sig, Weight <= -0.2 | Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_50_60_bac_edges_sig$Correlation <- ifelse(spearman_50_60_bac_edges_sig$Weight > 0, "pos", "neg")
spearman_50_60_bac_edges_sig$Type <- "Directed"
spearman_50_60_bac_edges_sig_final <- select(spearman_50_60_bac_edges_sig, c("Source", "Target", "Type", "Correlation"))
spearman_50_60_bac_edges_sig_final$Genus <- sapply(strsplit(as.character(spearman_50_60_bac_edges_sig_final$Target)," "), `[`, 1)
#select both columns to create node list
spearman_50_60_bac_nodes_sig <- select(spearman_50_60_bac_edges_sig, c("Target"))
spearman_50_60_bac_nodes_sig$Target2 <- spearman_50_60_bac_nodes_sig$Target
spearman_50_60_bac_nodes_sig$Correlation <- spearman_50_60_bac_edges_sig$Correlation

# remove duplicate entries
spearman_50_60_bac_nodes_sig = spearman_50_60_bac_nodes_sig[!duplicated(spearman_50_60_bac_nodes_sig$Target),]
# make data frame
spearman_50_60_bac_nodes_sig <- data.frame(spearman_50_60_bac_nodes_sig)
# re-index rows
rownames(spearman_50_60_bac_nodes_sig) <- NULL
# rename columns
colnames(spearman_50_60_bac_nodes_sig) <- c("Id", "Label", "Correlation")
# add genus information
spearman_50_60_bac_nodes_sig$Genus <- sapply(strsplit(as.character(spearman_50_60_bac_nodes_sig$Label)," "), `[`, 1)
mean_list_50_60_bac = c()
for (items in spearman_50_60_bac_nodes_sig$Id){
  mean_list_50_60_bac = append(mean_list_50_60_bac, list_50_60_bac_mean_abundance[items])}
spearman_50_60_bac_nodes_sig$mean_abundance <- mean_list_50_60_bac
# export node list
write.table(spearman_50_60_bac_nodes_sig, file="network_tables/nodes_50_60_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)
# export edge list
write.table(spearman_50_60_bac_edges_sig_final, file="network_tables/edges_50_60_bac.csv", sep=";", col.names = TRUE, row.names = FALSE)



###########################################################################################################################

# Gephi analysis
# 10_19 age group 
gephi_out_10_19 <- read_csv("C:/Users/marie/Desktop/Ajith_Tavaran/Gephi_analysis/gephi_out_10er.csv", 
                            col_types = cols(correlation = col_factor(levels = c("pos", "neg")), indegree = col_number(), 
                                             outdegree = col_number(), Degree = col_number(), 
                                             `weighted indegree` = col_number(), 
                                             `weighted outdegree` = col_number(), 
                                             `Weighted Degree` = col_number(), 
                                             Eccentricity = col_number(), closnesscentrality = col_number(), 
                                             harmonicclosnesscentrality = col_number(), 
                                             betweenesscentrality = col_number(), 
                                             Authority = col_number(), Hub = col_number(), 
                                             modularity_class = col_number(), 
                                             pageranks = col_number(), componentnumber = col_number(), 
                                             strongcompnum = col_number(), clustering = col_number(), 
                                             eigencentrality = col_number()))
gephi_out_10_19 <- data.frame(gephi_out_10_19)
gephi_out_10_19$age_group <- "10_19"
gephi_out_10_19$species_type <- ifelse(gephi_out_10_19$Id %in% list_bac_10_19_rare, "rare", "core")
gephi_out_10_19$timeset <- NULL
gephi_out_10_19$modularity_class_cat <- ifelse(gephi_out_10_19$modularity_class == 0, "module_1",
                                               ifelse(gephi_out_10_19$modularity_class == 1, "module_2", "module_3"))

gephi_out_10_19$Id <- str_replace_all(gephi_out_10_19$Id, "Schaalia..odontolytica", "Schaalia.odontolytica")
gephi_out_10_19$Id <- str_replace_all(gephi_out_10_19$Id, "Veillonella..dispar", "Veillonella.dispar")
gephi_out_10_19$Id <- str_replace_all(gephi_out_10_19$Id, "Gemella..sanguinis", "Gemella.sanguinis")
gephi_out_10_19$Id <- str_replace_all(gephi_out_10_19$Id, "Gemella..haemolysans", "Gemella.haemolysans")
gephi_out_10_19$Id <- str_replace_all(gephi_out_10_19$Id, "Gemella..morbillorum", "Gemella.morbillorum")
gephi_out_10_19$Id <- str_replace_all(gephi_out_10_19$Id, "Lautropia..mirabilis", "Lautropia.mirabilis")

# 20_29 age group 
gephi_out_20_29 <- read_csv("C:/Users/marie/Desktop/Ajith_Tavaran/Gephi_analysis/gephi_out_20er.csv", 
                            col_types = cols(correlation = col_factor(levels = c("pos", "neg")), indegree = col_number(), 
                                             outdegree = col_number(), Degree = col_number(), 
                                             `weighted indegree` = col_number(), 
                                             `weighted outdegree` = col_number(), 
                                             `Weighted Degree` = col_number(), 
                                             Eccentricity = col_number(), closnesscentrality = col_number(), 
                                             harmonicclosnesscentrality = col_number(), 
                                             betweenesscentrality = col_number(), 
                                             Authority = col_number(), Hub = col_number(), 
                                             modularity_class = col_number(), 
                                             pageranks = col_number(), componentnumber = col_number(), 
                                             strongcompnum = col_number(), clustering = col_number(), 
                                             eigencentrality = col_number()))
gephi_out_20_29 <- data.frame(gephi_out_20_29)
gephi_out_20_29$age_group <- "20_29"
gephi_out_20_29$species_type <- ifelse(gephi_out_20_29$Id %in% list_bac_20_29_rare, "rare", "core")
gephi_out_20_29$timeset <- NULL
gephi_out_20_29$modularity_class_cat <- ifelse(gephi_out_20_29$modularity_class == 0, "module_1",
                                               ifelse(gephi_out_20_29$modularity_class == 1, "module_2", "module_3"))

gephi_out_20_29$Id <- str_replace_all(gephi_out_20_29$Id, "Schaalia..odontolytica", "Schaalia.odontolytica")
gephi_out_20_29$Id <- str_replace_all(gephi_out_20_29$Id, "Veillonella..dispar", "Veillonella.dispar")
gephi_out_20_29$Id <- str_replace_all(gephi_out_20_29$Id, "Gemella..sanguinis", "Gemella.sanguinis")
gephi_out_20_29$Id <- str_replace_all(gephi_out_20_29$Id, "Gemella..haemolysans", "Gemella.haemolysans")
gephi_out_20_29$Id <- str_replace_all(gephi_out_20_29$Id, "Gemella..morbillorum", "Gemella.morbillorum")
gephi_out_20_29$Id <- str_replace_all(gephi_out_20_29$Id, "Lautropia..mirabilis", "Lautropia.mirabilis")

gephi_out_30_39 <- read_csv("C:/Users/marie/Desktop/Ajith_Tavaran/Gephi_analysis/gephi_out_30er.csv", 
                            col_types = cols(correlation = col_factor(levels = c("pos", "neg")), indegree = col_number(), 
                                             outdegree = col_number(), Degree = col_number(), 
                                             `weighted indegree` = col_number(), 
                                             `weighted outdegree` = col_number(), 
                                             `Weighted Degree` = col_number(), 
                                             Eccentricity = col_number(), closnesscentrality = col_number(), 
                                             harmonicclosnesscentrality = col_number(), 
                                             betweenesscentrality = col_number(), 
                                             Authority = col_number(), Hub = col_number(), 
                                             modularity_class = col_number(), 
                                             pageranks = col_number(), componentnumber = col_number(), 
                                             strongcompnum = col_number(), clustering = col_number(), 
                                             eigencentrality = col_number()))
gephi_out_30_39 <- data.frame(gephi_out_30_39)
gephi_out_30_39$age_group <- "30_39"
gephi_out_30_39$species_type <- ifelse(gephi_out_30_39$Id %in% list_bac_30_39_rare, "rare", "core")
gephi_out_30_39$timeset <- NULL
gephi_out_30_39$modularity_class_cat <- ifelse(gephi_out_30_39$modularity_class == 1, "module_1",
                                               ifelse(gephi_out_30_39$modularity_class == 0, "module_2", "module_3"))

gephi_out_30_39$Id <- str_replace_all(gephi_out_30_39$Id, "Schaalia..odontolytica", "Schaalia.odontolytica")
gephi_out_30_39$Id <- str_replace_all(gephi_out_30_39$Id, "Veillonella..dispar", "Veillonella.dispar")
gephi_out_30_39$Id <- str_replace_all(gephi_out_30_39$Id, "Gemella..sanguinis", "Gemella.sanguinis")
gephi_out_30_39$Id <- str_replace_all(gephi_out_30_39$Id, "Gemella..haemolysans", "Gemella.haemolysans")
gephi_out_30_39$Id <- str_replace_all(gephi_out_30_39$Id, "Gemella..morbillorum", "Gemella.morbillorum")
gephi_out_30_39$Id <- str_replace_all(gephi_out_30_39$Id, "Lautropia..mirabilis", "Lautropia.mirabilis")

gephi_out_40_49 <- read_csv("C:/Users/marie/Desktop/Ajith_Tavaran/Gephi_analysis/gephi_out_40er.csv", 
                            col_types = cols(correlation = col_factor(levels = c("pos", "neg")), indegree = col_number(), 
                                             outdegree = col_number(), Degree = col_number(), 
                                             `weighted indegree` = col_number(), 
                                             `weighted outdegree` = col_number(), 
                                             `Weighted Degree` = col_number(), 
                                             Eccentricity = col_number(), closnesscentrality = col_number(), 
                                             harmonicclosnesscentrality = col_number(), 
                                             betweenesscentrality = col_number(), 
                                             Authority = col_number(), Hub = col_number(), 
                                             modularity_class = col_number(), 
                                             pageranks = col_number(), componentnumber = col_number(), 
                                             strongcompnum = col_number(), clustering = col_number(), 
                                             eigencentrality = col_number()))
gephi_out_40_49 <- data.frame(gephi_out_40_49)
gephi_out_40_49$age_group <- "40_49"
gephi_out_40_49$species_type <- ifelse(gephi_out_40_49$Id %in% list_bac_40_49_rare, "rare", "core")
gephi_out_40_49$timeset <- NULL
gephi_out_40_49$modularity_class_cat <- ifelse(gephi_out_40_49$modularity_class == 0, "module_1",
                                               ifelse(gephi_out_40_49$modularity_class == 2, "module_2", "module_3"))

gephi_out_40_49$Id <- str_replace_all(gephi_out_40_49$Id, "Schaalia..odontolytica", "Schaalia.odontolytica")
gephi_out_40_49$Id <- str_replace_all(gephi_out_40_49$Id, "Veillonella..dispar", "Veillonella.dispar")
gephi_out_40_49$Id <- str_replace_all(gephi_out_40_49$Id, "Gemella..sanguinis", "Gemella.sanguinis")
gephi_out_40_49$Id <- str_replace_all(gephi_out_40_49$Id, "Gemella..haemolysans", "Gemella.haemolysans")
gephi_out_40_49$Id <- str_replace_all(gephi_out_40_49$Id, "Gemella..morbillorum", "Gemella.morbillorum")
gephi_out_40_49$Id <- str_replace_all(gephi_out_40_49$Id, "Lautropia..mirabilis", "Lautropia.mirabilis")


gephi_out_50_59 <- read_csv("C:/Users/marie/Desktop/Ajith_Tavaran/Gephi_analysis/gephi_out_50er.csv", 
                            col_types = cols(correlation = col_factor(levels = c("pos", "neg")), indegree = col_number(), 
                                             outdegree = col_number(), Degree = col_number(), 
                                             `weighted indegree` = col_number(), 
                                             `weighted outdegree` = col_number(), 
                                             `Weighted Degree` = col_number(), 
                                             Eccentricity = col_number(), closnesscentrality = col_number(), 
                                             harmonicclosnesscentrality = col_number(), 
                                             betweenesscentrality = col_number(), 
                                             Authority = col_number(), Hub = col_number(), 
                                             modularity_class = col_number(), 
                                             pageranks = col_number(), componentnumber = col_number(), 
                                             strongcompnum = col_number(), clustering = col_number(), 
                                             eigencentrality = col_number()))
gephi_out_50_59 <- data.frame(gephi_out_50_59)
gephi_out_50_59$age_group <- "50_59"
gephi_out_50_59$species_type <- ifelse(gephi_out_50_59$Id %in% list_bac_50_60_rare, "rare", "core")
gephi_out_50_59$timeset <- NULL
gephi_out_50_59$modularity_class_cat <- ifelse(gephi_out_50_59$modularity_class == 0, "module_1",
                                               ifelse(gephi_out_50_59$modularity_class == 1, "module_2", 
                                                      ifelse(gephi_out_50_59$modularity_class == 3,"module_3", "module_4")))

gephi_out_50_59$Id <- str_replace_all(gephi_out_50_59$Id, "Schaalia..odontolytica", "Schaalia.odontolytica")
gephi_out_50_59$Id <- str_replace_all(gephi_out_50_59$Id, "Veillonella..dispar", "Veillonella.dispar")
gephi_out_50_59$Id <- str_replace_all(gephi_out_50_59$Id, "Gemella..sanguinis", "Gemella.sanguinis")
gephi_out_50_59$Id <- str_replace_all(gephi_out_50_59$Id, "Gemella..haemolysans", "Gemella.haemolysans")
gephi_out_50_59$Id <- str_replace_all(gephi_out_50_59$Id, "Gemella..morbillorum", "Gemella.morbillorum")
gephi_out_50_59$Id <- str_replace_all(gephi_out_50_59$Id, "Lautropia..mirabilis", "Lautropia.mirabilis")

pathway_recons_gapseq <- read_delim("pathway_reconstruction_gapseq.csv", delim = ";", escape_double = FALSE, 
                                    col_types = cols(Completeness = col_number(), 
                                                     VagueReactions = col_number(), 
                                                     KeyReactions = col_number(), 
                                                     ReactionsFound = col_number()), trim_ws = TRUE)
pathway_recons_gapseq <- data.frame(pathway_recons_gapseq)
pathway_dict <- pathway_recons_gapseq$Name
names(pathway_dict) <- pathway_recons_gapseq$Id

# clean species names
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "chromosome", "")
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "genome", "")
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "complete", "")
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "sequence", "")
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "BAC", "")
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "strain", "")
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "subsp", "")
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "ATCC", "")
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "DSM", "")
pathway_recons_gapseq$Species <- str_replace_all(pathway_recons_gapseq$Species, "____", "")
pathway_recons_gapseq$Species1 <- unlist(map(str_split(pathway_recons_gapseq$Species, "_"),4))
pathway_recons_gapseq$Species2 <- unlist(map(str_split(pathway_recons_gapseq$Species, "_"),5))
pathway_recons_gapseq$Species3 <- paste(pathway_recons_gapseq$Species1, pathway_recons_gapseq$Species2)

pathway_recons_gapseq$Species <- NULL
pathway_recons_gapseq$Species1 <- NULL
pathway_recons_gapseq$Species2 <- NULL
pathway_recons_gapseq$Prediction <- NULL
pathway_recons_gapseq$VagueReactions <- NULL
pathway_recons_gapseq$KeyReactions <- NULL
pathway_recons_gapseq$ReactionsFound <- NULL
pathway_recons_gapseq$Species3 <- str_replace_all(pathway_recons_gapseq$Species3, " ", ".")
pathway_recons_gapseq$Species3 <- str_replace_all(pathway_recons_gapseq$Species3, "acidifaciens.", "Bacteroides.acidifaciens")
pathway_recons_gapseq$Species3 <- str_replace_all(pathway_recons_gapseq$Species3, "sanguinis.", "Streptococcus.sanguinis")
pathway_recons_gapseq$Species3 <- str_replace_all(pathway_recons_gapseq$Species3, "odontolytica.", "Schaalia.odontolytica")
pathway_recons_gapseq$Species3 <- str_replace_all(pathway_recons_gapseq$Species3, "mirabilis.", "Lautropia.mirabilis")
pathway_recons_gapseq$Species3 <- str_replace_all(pathway_recons_gapseq$Species3, "morbillorum.", "Gemella.morbillorum")
pathway_recons_gapseq$Species3 <- str_replace_all(pathway_recons_gapseq$Species3, "dispar.", "Veillonella.dispar")
pathway_recons_gapseq$Species3 <- str_replace_all(pathway_recons_gapseq$Species3, "haemolysans.", "Gemella.haemolysans")


# Gephi, Age = 10-19, Module 1-2
gephi_out_10_19_module_01 <- gephi_out_10_19[gephi_out_10_19$modularity_class_cat=="module_1",]
gephi_10_19_module_01 <- c(gephi_out_10_19_module_01$Id)
pathway_recons_gapseq_10_19_module_01 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_10_19_module_01,]
pathway_recons_gapseq_10_19_module_01_short <- ddply(pathway_recons_gapseq_10_19_module_01, "Id", numcolwise(median))
pathway_recons_gapseq_10_19_module_01_short$age_group <- "10_19"
pathway_recons_gapseq_10_19_module_01_short$module <- "Module_01"
pathway_recons_gapseq_10_19_module_01_short$age_module <- paste(pathway_recons_gapseq_10_19_module_01_short$age_group,
                                                                pathway_recons_gapseq_10_19_module_01_short$module)

gephi_out_10_19_module_02 <- gephi_out_10_19[gephi_out_10_19$modularity_class_cat=="module_2",]
gephi_10_19_module_02 <- c(gephi_out_10_19_module_02$Id)
pathway_recons_gapseq_10_19_module_02 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_10_19_module_02,]
pathway_recons_gapseq_10_19_module_02_short <- ddply(pathway_recons_gapseq_10_19_module_02, "Id", numcolwise(median))
pathway_recons_gapseq_10_19_module_02_short$age_group <- "10_19"
pathway_recons_gapseq_10_19_module_02_short$module <- "Module_02"
pathway_recons_gapseq_10_19_module_02_short$age_module <- paste(pathway_recons_gapseq_10_19_module_02_short$age_group,
                                                                pathway_recons_gapseq_10_19_module_02_short$module)

gephi_out_10_19_module_03 <- gephi_out_10_19[gephi_out_10_19$modularity_class_cat=="module_3",]
gephi_10_19_module_03 <- c(gephi_out_10_19_module_03$Id)
pathway_recons_gapseq_10_19_module_03 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_10_19_module_03,]
pathway_recons_gapseq_10_19_module_03_short <- ddply(pathway_recons_gapseq_10_19_module_03, "Id", numcolwise(median))
pathway_recons_gapseq_10_19_module_03_short$age_group <- "10_19"
pathway_recons_gapseq_10_19_module_03_short$module <- "Module_03"
pathway_recons_gapseq_10_19_module_03_short$age_module <- paste(pathway_recons_gapseq_10_19_module_03_short$age_group,
                                                                pathway_recons_gapseq_10_19_module_03_short$module)


# Gephi, Age = 20-29, Module 1-2
gephi_out_20_29_module_01 <- gephi_out_20_29[gephi_out_20_29$modularity_class_cat=="module_1",]
gephi_20_29_module_01 <- c(gephi_out_20_29_module_01$Id)
pathway_recons_gapseq_20_29_module_01 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_20_29_module_01,]
pathway_recons_gapseq_20_29_module_01_short <- ddply(pathway_recons_gapseq_20_29_module_01, "Id", numcolwise(median))
pathway_recons_gapseq_20_29_module_01_short$age_group <- "20_29"
pathway_recons_gapseq_20_29_module_01_short$module <- "Module_01"
pathway_recons_gapseq_20_29_module_01_short$age_module <- paste(pathway_recons_gapseq_20_29_module_01_short$age_group,
                                                                pathway_recons_gapseq_20_29_module_01_short$module)

gephi_out_20_29_module_02 <- gephi_out_20_29[gephi_out_20_29$modularity_class_cat=="module_2",]
gephi_20_29_module_02 <- c(gephi_out_20_29_module_02$Id)
pathway_recons_gapseq_20_29_module_02 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_20_29_module_02,]
pathway_recons_gapseq_20_29_module_02_short <- ddply(pathway_recons_gapseq_20_29_module_02, "Id", numcolwise(median))
pathway_recons_gapseq_20_29_module_02_short$age_group <- "20_29"
pathway_recons_gapseq_20_29_module_02_short$module <- "Module_02"
pathway_recons_gapseq_20_29_module_02_short$age_module <- paste(pathway_recons_gapseq_20_29_module_02_short$age_group,
                                                                pathway_recons_gapseq_20_29_module_02_short$module)

gephi_out_20_29_module_03 <- gephi_out_20_29[gephi_out_20_29$modularity_class_cat=="module_3",]
gephi_20_29_module_03 <- c(gephi_out_20_29_module_03$Id)
pathway_recons_gapseq_20_29_module_03 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_20_29_module_03,]
pathway_recons_gapseq_20_29_module_03_short <- ddply(pathway_recons_gapseq_20_29_module_03, "Id", numcolwise(median))
pathway_recons_gapseq_20_29_module_03_short$age_group <- "20_29"
pathway_recons_gapseq_20_29_module_03_short$module <- "Module_03"
pathway_recons_gapseq_20_29_module_03_short$age_module <- paste(pathway_recons_gapseq_20_29_module_03_short$age_group,
                                                                pathway_recons_gapseq_20_29_module_03_short$module)


# Gephi, Age = 30-39, Module 1-2
gephi_out_30_39_module_01 <- gephi_out_30_39[gephi_out_30_39$modularity_class_cat=="module_1",]
gephi_30_39_module_01 <- c(gephi_out_30_39_module_01$Id)
pathway_recons_gapseq_30_39_module_01 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_30_39_module_01,]
pathway_recons_gapseq_30_39_module_01_short <- ddply(pathway_recons_gapseq_30_39_module_01, "Id", numcolwise(median))
pathway_recons_gapseq_30_39_module_01_short$age_group <- "30_39"
pathway_recons_gapseq_30_39_module_01_short$module <- "Module_01"
pathway_recons_gapseq_30_39_module_01_short$age_module <- paste(pathway_recons_gapseq_30_39_module_01_short$age_group,
                                                                pathway_recons_gapseq_30_39_module_01_short$module)

gephi_out_30_39_module_02 <- gephi_out_30_39[gephi_out_30_39$modularity_class_cat=="module_2",]
gephi_30_39_module_02 <- c(gephi_out_30_39_module_02$Id)
pathway_recons_gapseq_30_39_module_02 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_30_39_module_02,]
pathway_recons_gapseq_30_39_module_02_short <- ddply(pathway_recons_gapseq_30_39_module_02, "Id", numcolwise(median))
pathway_recons_gapseq_30_39_module_02_short$age_group <- "30_39"
pathway_recons_gapseq_30_39_module_02_short$module <- "Module_02"
pathway_recons_gapseq_30_39_module_02_short$age_module <- paste(pathway_recons_gapseq_30_39_module_02_short$age_group,
                                                                pathway_recons_gapseq_30_39_module_02_short$module)

gephi_out_30_39_module_03 <- gephi_out_30_39[gephi_out_30_39$modularity_class_cat=="module_3",]
gephi_30_39_module_03 <- c(gephi_out_30_39_module_03$Id)
pathway_recons_gapseq_30_39_module_03 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_30_39_module_03,]
pathway_recons_gapseq_30_39_module_03_short <- ddply(pathway_recons_gapseq_30_39_module_03, "Id", numcolwise(median))
pathway_recons_gapseq_30_39_module_03_short$age_group <- "30_39"
pathway_recons_gapseq_30_39_module_03_short$module <- "Module_03"
pathway_recons_gapseq_30_39_module_03_short$age_module <- paste(pathway_recons_gapseq_30_39_module_03_short$age_group,
                                                                pathway_recons_gapseq_30_39_module_03_short$module)

# Gephi, Age = 40-49, Module 1-2
gephi_out_40_49_module_01 <- gephi_out_40_49[gephi_out_40_49$modularity_class_cat=="module_1",]
gephi_40_49_module_01 <- c(gephi_out_40_49_module_01$Id)
pathway_recons_gapseq_40_49_module_01 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_40_49_module_01,]
pathway_recons_gapseq_40_49_module_01_short <- ddply(pathway_recons_gapseq_40_49_module_01, "Id", numcolwise(median))
pathway_recons_gapseq_40_49_module_01_short$age_group <- "40_49"
pathway_recons_gapseq_40_49_module_01_short$module <- "Module_01"
pathway_recons_gapseq_40_49_module_01_short$age_module <- paste(pathway_recons_gapseq_40_49_module_01_short$age_group,
                                                                pathway_recons_gapseq_40_49_module_01_short$module)

gephi_out_40_49_module_02 <- gephi_out_40_49[gephi_out_40_49$modularity_class_cat=="module_2",]
gephi_40_49_module_02 <- c(gephi_out_40_49_module_02$Id)
pathway_recons_gapseq_40_49_module_02 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_40_49_module_02,]
pathway_recons_gapseq_40_49_module_02_short <- ddply(pathway_recons_gapseq_40_49_module_02, "Id", numcolwise(median))
pathway_recons_gapseq_40_49_module_02_short$age_group <- "40_49"
pathway_recons_gapseq_40_49_module_02_short$module <- "Module_02"
pathway_recons_gapseq_40_49_module_02_short$age_module <- paste(pathway_recons_gapseq_40_49_module_02_short$age_group,
                                                                pathway_recons_gapseq_40_49_module_02_short$module)


gephi_out_40_49_module_03 <- gephi_out_40_49[gephi_out_40_49$modularity_class_cat=="module_3",]
gephi_40_49_module_03 <- c(gephi_out_40_49_module_02$Id)
pathway_recons_gapseq_40_49_module_03 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_40_49_module_03,]
pathway_recons_gapseq_40_49_module_03_short <- ddply(pathway_recons_gapseq_40_49_module_03, "Id", numcolwise(median))
pathway_recons_gapseq_40_49_module_03_short$age_group <- "40_49"
pathway_recons_gapseq_40_49_module_03_short$module <- "Module_03"
pathway_recons_gapseq_40_49_module_03_short$age_module <- paste(pathway_recons_gapseq_40_49_module_03_short$age_group,
                                                                pathway_recons_gapseq_40_49_module_03_short$module)


# Gephi, Age = 50-59, Module 1-2
gephi_out_50_59_module_01 <- gephi_out_50_59[gephi_out_50_59$modularity_class_cat=="module_1",]
gephi_50_59_module_01 <- c(gephi_out_50_59_module_01$Id)
pathway_recons_gapseq_50_59_module_01 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_50_59_module_01,]
pathway_recons_gapseq_50_59_module_01_short <- ddply(pathway_recons_gapseq_50_59_module_01, "Id", numcolwise(median))
pathway_recons_gapseq_50_59_module_01_short$age_group <- "50_59"
pathway_recons_gapseq_50_59_module_01_short$module <- "Module_01"
pathway_recons_gapseq_50_59_module_01_short$age_module <- paste(pathway_recons_gapseq_50_59_module_01_short$age_group,
                                                                pathway_recons_gapseq_50_59_module_01_short$module)

gephi_out_50_59_module_02 <- gephi_out_50_59[gephi_out_50_59$modularity_class_cat=="module_2",]
gephi_50_59_module_02 <- c(gephi_out_50_59_module_02$Id)
pathway_recons_gapseq_50_59_module_02 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_50_59_module_02,]
pathway_recons_gapseq_50_59_module_02_short <- ddply(pathway_recons_gapseq_50_59_module_02, "Id", numcolwise(median))
pathway_recons_gapseq_50_59_module_02_short$age_group <- "50_59"
pathway_recons_gapseq_50_59_module_02_short$module <- "Module_02"
pathway_recons_gapseq_50_59_module_02_short$age_module <- paste(pathway_recons_gapseq_50_59_module_02_short$age_group,
                                                                pathway_recons_gapseq_50_59_module_02_short$module)

gephi_out_50_59_module_03 <- gephi_out_50_59[gephi_out_50_59$modularity_class_cat=="module_3",]
gephi_50_59_module_03 <- c(gephi_out_50_59_module_03$Id)
pathway_recons_gapseq_50_59_module_03 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_50_59_module_03,]
pathway_recons_gapseq_50_59_module_03_short <- ddply(pathway_recons_gapseq_50_59_module_03, "Id", numcolwise(median))
pathway_recons_gapseq_50_59_module_03_short$age_group <- "50_59"
pathway_recons_gapseq_50_59_module_03_short$module <- "Module_03"
pathway_recons_gapseq_50_59_module_03_short$age_module <- paste(pathway_recons_gapseq_50_59_module_03_short$age_group,
                                                                pathway_recons_gapseq_50_59_module_03_short$module)

gephi_out_50_59_module_04 <- gephi_out_50_59[gephi_out_50_59$modularity_class_cat=="module_4",]
gephi_50_59_module_04 <- c(gephi_out_50_59_module_04$Id)
pathway_recons_gapseq_50_59_module_04 <- pathway_recons_gapseq[pathway_recons_gapseq$Species3 %in% gephi_50_59_module_04,]
pathway_recons_gapseq_50_59_module_04_short <- ddply(pathway_recons_gapseq_50_59_module_04, "Id", numcolwise(median))
pathway_recons_gapseq_50_59_module_04_short$age_group <- "50_59"
pathway_recons_gapseq_50_59_module_04_short$module <- "Module_04"
pathway_recons_gapseq_50_59_module_04_short$age_module <- paste(pathway_recons_gapseq_50_59_module_04_short$age_group,
                                                                pathway_recons_gapseq_50_59_module_04_short$module)

pathway_recons_gapseq_all <- data.frame(rbind(pathway_recons_gapseq_10_19_module_01_short,
                                              pathway_recons_gapseq_10_19_module_02_short,
                                              pathway_recons_gapseq_10_19_module_03_short,
                                              pathway_recons_gapseq_20_29_module_01_short,
                                              pathway_recons_gapseq_20_29_module_02_short,
                                              pathway_recons_gapseq_20_29_module_03_short,
                                              pathway_recons_gapseq_30_39_module_01_short,
                                              pathway_recons_gapseq_30_39_module_02_short,
                                              pathway_recons_gapseq_30_39_module_03_short,
                                              pathway_recons_gapseq_40_49_module_01_short,
                                              pathway_recons_gapseq_40_49_module_02_short,
                                              pathway_recons_gapseq_40_49_module_03_short,
                                              pathway_recons_gapseq_50_59_module_01_short,
                                              pathway_recons_gapseq_50_59_module_02_short,
                                              pathway_recons_gapseq_50_59_module_03_short,
                                              pathway_recons_gapseq_50_59_module_04_short))

pathway_recons_gapseq_all$module <- NULL
pathway_recons_gapseq_all$age_group <- NULL
pathway_recons_gapseq_W <- spread(pathway_recons_gapseq_all, key=Id, value=Completeness)
rownames(pathway_recons_gapseq_W) <- pathway_recons_gapseq_W$age_module
pathway_recons_gapseq_W$age_module <- NULL
pathway_recons_gapseq_W_t <- data.frame(t(pathway_recons_gapseq_W))


# make principal component analysis (center and scale FM scores)
keep_all <- apply(pathway_recons_gapseq_W_t, 1, function(x) length(unique(x[!is.na(x)])) != 1)
pathway_recons_gapseq_W_02 <- pathway_recons_gapseq_W_t[keep_all, ]

pca_test <- prcomp(t(pathway_recons_gapseq_W_02), scale = TRUE, center = TRUE)

pca_df <- data.frame(pca_test$x)
pca_df$age <- unlist(map(str_split(rownames(pca_df), "\\."),1))
pca_df$module <- unlist(map(str_split(rownames(pca_df), "\\."),2))
pca_df$age <- str_replace(pca_df$age, "X", "")
pca_df$age <- str_replace(pca_df$age, "_", "-")
pca_df$module <- str_replace(pca_df$module, "_", " ")

pca_plot <-
  ggplot(pca_df, aes(x=PC1, y=PC2, colour=age, shape=module)) + 
  geom_jitter(size=4) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.1) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.1) +
  scale_colour_manual("Age group [years]", 
                      values=c("10-19"="darkblue", 
                               "20-29"="lightblue", 
                               "30-39"="orange", 
                               "40-49"="red", 
                               "50-59"="darkred")) +
  scale_shape_manual("Modularity class", 
                     values=c("Module 01"=16, 
                              "Module 02"=17, 
                              "Module 03"=3, 
                              "Module 04"=12)) +
  theme_bw() + 
  theme(legend.position = "right", 
        panel.grid = element_blank(),
        strip.background = element_rect(fill=NA)) + 
  xlab("PC1 (28 %)") + ylab("PC2 (22 %)")


# extract PCA results
res.var <- get_pca_var(pca_test)
# Contributions to the PCs
variable_contributions <- data.frame(res.var$contrib)
variable_contributions$Pathway <- rownames(variable_contributions)
variable_coord <- data.frame(res.var$coord)
variable_coord$Pathway <- rownames(variable_coord)

# Dimension 1
variable_contributions_dim1 <- select(variable_contributions, c("Pathway", "Dim.1"))
rownames(variable_contributions_dim1) <- NULL
variable_contributions_dim1 <- variable_contributions_dim1[order(variable_contributions_dim1$Dim.1, decreasing = TRUE),]
variable_contributions_dim1$cumsum <- cumsum(variable_contributions_dim1$Dim.1)
variable_contributions_dim1_mostVariance <- variable_contributions_dim1[1:50,]
variable_contributions_dim1_mostVariance$cumsum <- NULL

variable_coord_dim1 <- select(variable_coord, c("Pathway", "Dim.1"))
colnames(variable_coord_dim1) <- c("Pathway", "Dim.1_coord")
rownames(variable_coord_dim1) <- NULL
variable_coord_dim1_mostVariance <- variable_coord_dim1[variable_coord_dim1$Pathway %in% variable_contributions_dim1_mostVariance$Pathway,]
variable_dim_1 <- inner_join(variable_contributions_dim1_mostVariance, variable_coord_dim1_mostVariance, by="Pathway")
variable_dim_1$value <- ifelse(variable_dim_1$Dim.1_coord < 0, variable_dim_1$Dim.1 * -1, variable_dim_1$Dim.1)
variable_dim_1$value <- round(variable_dim_1$value, 2)
variable_dim_1$charge <- ifelse(variable_dim_1$value > 0, "pos", "neg")
variable_dim_1$Pathway <- str_replace_all(variable_dim_1$Pathway, "\\|", "")
variable_dim_1$Dimension <- "Dim01"
colnames(variable_dim_1) <- c("Pathway", "Contribution", "Coordinate", "Value", "Charge", "Dimension")


# Dimension 2
variable_contributions_dim2 <- select(variable_contributions, c("Pathway", "Dim.2"))
rownames(variable_contributions_dim2) <- NULL
variable_contributions_dim2 <- variable_contributions_dim2[order(variable_contributions_dim2$Dim.2, decreasing = TRUE),]
variable_contributions_dim2$cumsum <- cumsum(variable_contributions_dim2$Dim.2)
variable_contributions_dim2_mostVariance <- variable_contributions_dim2[1:50,]
variable_contributions_dim2_mostVariance$cumsum <- NULL

variable_coord_dim2 <- select(variable_coord, c("Pathway", "Dim.2"))
colnames(variable_coord_dim2) <- c("Pathway", "Dim.2_coord")
rownames(variable_coord_dim2) <- NULL
variable_coord_dim2_mostVariance <- variable_coord_dim2[variable_coord_dim2$Pathway %in% variable_contributions_dim2_mostVariance$Pathway,]
variable_dim_2 <- inner_join(variable_contributions_dim2_mostVariance, variable_coord_dim2_mostVariance, by="Pathway")
variable_dim_2$value <- ifelse(variable_dim_2$Dim.2_coord < 0, variable_dim_2$Dim.2 * -1, variable_dim_2$Dim.2)
variable_dim_2$value <- round(variable_dim_2$value, 2)
variable_dim_2$charge <- ifelse(variable_dim_2$value > 0, "pos", "neg")
variable_dim_2$Pathway <- str_replace_all(variable_dim_2$Pathway, "\\|", "")
variable_dim_2$Dimension <- "Dim02"
colnames(variable_dim_2) <- c("Pathway", "Contribution", "Coordinate", "Value", "Charge", "Dimension")

rownames(pathway_recons_gapseq_W_t) <- str_replace_all(rownames(pathway_recons_gapseq_W_t), "\\|", "")
pathway_recons_gapseq_W_t_module_dim02 <- pathway_recons_gapseq_W_t[rownames(pathway_recons_gapseq_W_t)%in%variable_dim_2$Pathway,]
colnames(pathway_recons_gapseq_W_t_module_dim02) <- c("10-19: Mod_01", "10-19: Mod_02", "10-19: Mod_03",
                                                      "20-29: Mod_01", "20-29: Mod_02", "20-29: Mod_03",
                                                      "30-39: Mod_01", "30-39: Mod_02", "30-39: Mod_03",
                                                      "40-49: Mod_01", "40-49: Mod_02", "40-49: Mod_03",
                                                      "50-59: Mod_01", "50-59: Mod_02", "50-59: Mod_03", "50-59: Mod_04")
dim1_names_df <-data.frame(Pathway=rownames(pathway_recons_gapseq_W_t_module_dim01), Dimension="PC1")
dim2_names_df <-data.frame(Pathway=rownames(pathway_recons_gapseq_W_t_module_dim02), Dimension="PC2")
dim_names_df <- data.frame(rbind(dim1_names_df, dim2_names_df))
# write.table(dim_names_df, "MetaCyc_nomenclature_Figure06.csv", row.names = FALSE, col.names = TRUE, sep=";")

pheatmap1 <- pheatmap(pathway_recons_gapseq_W_t_module_dim02, 
                      scale = "none", 
                      cutree_cols = 3, 
                      cutree_rows = 3,
                      cellwidth = 16,
                      cellheight = 8,
                      angle_col = 90,
                      show_rownames = TRUE,
                      border_color = "black",
                      fontsize_row = 5,
                      fontsize_col =7,
                      clustering_distance_cols = "euclidean",
                      clustering_method = "ward.D", 
                      treeheight_row = 0, treeheight_col = 12)
pheatmap1_gg <- as.ggplot(pheatmap1)

pathway_recons_gapseq_W_t_module_dim01 <- pathway_recons_gapseq_W_t[rownames(pathway_recons_gapseq_W_t)%in%variable_dim_1$Pathway,]
colnames(pathway_recons_gapseq_W_t_module_dim01) <- c("10-19: Mod_01", "10-19: Mod_02", "10-19: Mod_03",
                                                      "20-29: Mod_01", "20-29: Mod_02", "20-29: Mod_03",
                                                      "30-39: Mod_01", "30-39: Mod_02", "30-39: Mod_03",
                                                      "40-49: Mod_01", "40-49: Mod_02", "40-49: Mod_03",
                                                      "50-59: Mod_01", "50-59: Mod_02", "50-59: Mod_03", "50-59: Mod_04")
pheatmap2 <- pheatmap(pathway_recons_gapseq_W_t_module_dim01, 
                      scale = "none", 
                      cutree_cols = 4, 
                      cutree_rows = 3,
                      cellwidth = 16,
                      cellheight = 8,
                      angle_col = 90,
                      show_rownames = TRUE,
                      border_color = "black",
                      fontsize_row = 5,
                      fontsize_col = 7,
                      legend = FALSE,
                      clustering_distance_cols = "euclidean",
                      clustering_method = "ward.D", 
                      treeheight_row = 0, treeheight_col = 12)
pheatmap2_gg <- as.ggplot(pheatmap2)

#####################################################################################################################
# Network statistics
gephi_out <- data.frame(rbind(gephi_out_10_19, gephi_out_20_29, gephi_out_30_39, gephi_out_40_49, gephi_out_50_59))
gephi_out$age_group_rare <- paste(gephi_out$age_group,gephi_out$species_type)

compare_it_list <- list(c("10_19", "20_29"),
                        c("10_19", "30_39"),
                        c("10_19", "40_49"),
                        c("10_19", "50_59"),
                        c("30_39", "50_59"))
closness_plot <-
  ggplot(gephi_out, aes(x=age_group, y=closnesscentrality)) +
  geom_violin(width=0.5) +
  geom_jitter(width=0.05, size=0.8) + 
  theme_bw() +
  stat_summary(
    mapping = aes(x = age_group, y = closnesscentrality),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median, colour="red") +
  scale_x_discrete(labels=c("10-19", "20-29", "30-39", "40-49", "50-59")) +
  stat_compare_means(comparisons = compare_it_list, label = "p.signif") + xlab(" ") + ylab("Closness centrality") +
  theme(panel.grid = element_blank())

betweeness_plot <-
  ggplot(gephi_out, aes(x=age_group, y=log10(betweenesscentrality+1))) +
  geom_violin(width=0.5) +
  geom_jitter(width=0.05, size=0.8) + 
  theme_bw() +
  stat_summary(
    mapping = aes(x = age_group, y = log10(betweenesscentrality+1)),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median, colour="red") +
  scale_x_discrete(labels=c("10-19", "20-29", "30-39", "40-49", "50-59")) +
  stat_compare_means(comparisons = compare_it_list, label = "p.signif") + xlab("Age [years]") + ylab("Betweeness centrality") +
  theme(panel.grid = element_blank())

degree_plot <-
  ggplot(gephi_out, aes(x=age_group, y=Weighted.Degree)) +
  geom_violin(width=0.5) +
  geom_jitter(width=0.05, size=0.8) + 
  theme_bw() +
  stat_summary(
    mapping = aes(x = age_group, y = (Weighted.Degree)),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median, colour="red") +
  scale_x_discrete(labels=c("10-19", "20-29", "30-39", "40-49", "50-59")) +
  stat_compare_means(comparisons = compare_it_list, label = "p.signif") + xlab(" ") + ylab("Weighted degree") +
  theme(panel.grid = element_blank())


authority_plot <-
  ggplot(gephi_out, aes(x=age_group, y=Authority)) +
  geom_violin(width=0.5) +
  geom_jitter(width=0.05, size=0.8) + 
  theme_bw() +
  stat_summary(
    mapping = aes(x = age_group, y = (Authority)),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median, colour="red") +
  scale_x_discrete(labels=c("10-19", "20-29", "30-39", "40-49", "50-59")) +
  stat_compare_means(comparisons = compare_it_list, label = "p.signif") + xlab("Age [years]") + ylab("Authority score") +
  theme(panel.grid = element_blank())

# rcompanion::epsilonSquared(gephi_out$closnesscentrality, g=gephi_out$age_group, ci=TRUE) # e2 = 0.52, CI = 0.46 - 0.58
# rcompanion::epsilonSquared(gephi_out$betweenesscentrality, g=gephi_out$age_group, ci=TRUE) # e2=0.12. CI=0.1-0.2
# rcompanion::epsilonSquared(gephi_out$Weighted.Degree, g=gephi_out$age_group, ci=TRUE) # e2=0.41, CI=0.33-0.5


net_work_plots <- ggarrange(closness_plot, betweeness_plot, degree_plot, nrow=1, labels = c("A", "B", "C"))

pdf("net_work_statistics_age.pdf", width=10, height = 4)
net_work_plots
dev.off()

table_it <- data.frame(table(gephi_out$modularity_class_cat, gephi_out$age_group, gephi_out$species_type))

beta_network <- select(gephi_out, c("Label", "mean_abundance", "species_type", "age_group", "modularity_class_cat"))
beta_network$age_type_mod <- paste(beta_network$age_group,"_",beta_network$modularity_class_cat)
beta_network$species_type <- NULL
beta_network$modularity_class_cat <- NULL
beta_network$age_group <- NULL
beta_network_df <- spread(beta_network, key=Label, value=mean_abundance)
rownames(beta_network_df) <- beta_network_df$age_type_mod 
beta_network_df[is.na(beta_network_df)] <- 0
beta_network_df$age_type_mod <- NULL

beta_network_df
pca_df <- prcomp(t(beta_network_df), center = TRUE)
pc_df_2 <- data.frame(pca_df$rotation)

PoV <- pca_df$sdev^2/sum(pca_df$sdev^2)
PoV_2 <- round(PoV * 100)
pc_df_2$age_group <- unlist(map(str_split(rownames(pc_df_2)," _ "),1))
pc_df_2$mod <- unlist(map(str_split(rownames(pc_df_2)," _ "),2))
pc_df_2$age_group <- str_replace(pc_df_2$age_group, "_", "-")
pc_df_2$mod <- str_replace(pc_df_2$mod, "module_", "Module ")
pc_df_2$mod <- str_replace(pc_df_2$mod, "1", "01")
pc_df_2$mod <- str_replace(pc_df_2$mod, "2", "02")
pc_df_2$mod <- str_replace(pc_df_2$mod, "3", "03")
pc_df_2$mod <- str_replace(pc_df_2$mod, "4", "04")


modularity_comp_plot <-
  ggplot(pc_df_2, aes(x=PC1, y=PC2, colour=age_group, shape=mod, group=age_group)) +
  geom_jitter(size=4) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.1) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.1) +
  scale_colour_manual("Age group [years]", 
                      values=c("10-19"="darkblue", 
                               "20-29"="lightblue", 
                               "30-39"="orange", 
                               "40-49"="red", 
                               "50-59"="darkred")) +
  scale_shape_manual("Modularity class", 
                     values=c("Module 01"=16, 
                              "Module 02"=17, 
                              "Module 03"=3, 
                              "Module 04"=12)) +
  xlab(paste("PC1 (",PoV_2[1],"%)")) +
  ylab(paste("PC2 (",PoV_2[2],"%)")) +
  theme_bw() + 
  theme(legend.position = "right", 
        panel.grid = element_blank(),
        strip.background = element_rect(fill=NA)) + xlim(-0.5,1) + ylim(-1,0.5)

pheatmap12 <- ggarrange(pheatmap1_gg, pheatmap2_gg, common.legend = TRUE, labels = c("B", "C"))
pheatmap_pca <- ggarrange(pca_plot, pheatmap12, nrow=2, heights = c(0.4,1), labels = c("A", "B"))

pca_plots <- ggarrange(modularity_comp_plot,pca_plot, common.legend = TRUE, legend = "right",labels = c("A", "B"))
pheatmap_plots <- ggarrange(pheatmap2_gg, pheatmap1_gg, common.legend = TRUE, labels = c("C", "D"))
pca_pheatmaps <- ggarrange(pca_plots, pheatmap_plots, nrow=2, heights = c(0.6,1))

pdf("pheatmap_pca_function.pdf", width = 12, height = 11)
pca_pheatmaps
dev.off()



##########################################################################################
# Corrplot analysis
pathway_abund <- read_delim("pathabundance_merged_norm.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
pathway_abund <- data.frame(pathway_abund)
colnames(pathway_abund) <- str_replace(colnames(pathway_abund), "Meta_", "")
colnames(pathway_abund) <- str_replace(colnames(pathway_abund), "_R1R2_kneaddata_Abundance", "")
colnames(pathway_abund) <- sub("\\_S[0-9]*", "", colnames(pathway_abund))
rownames(pathway_abund) <- pathway_abund$Pathway
pathway_abund$Pathway <- NULL
remove_items <- c("UNMAPPED", "UNINTEGRATED")
pathway_abund <- pathway_abund[!rownames(pathway_abund)%in%remove_items,]
list_pathways <- unlist(map(str_split(rownames(pathway_abund), ": "),2))
names(list_pathways) <- unlist(map(str_split(rownames(pathway_abund), ": "),1))

rownames(pathway_abund) <- unlist(map(str_split(rownames(pathway_abund), ": "),1))
pathway_abund_clr <- clr(pathway_abund)
pathway_abund_clr_t <- data.frame(t(pathway_abund_clr))
pathway_abund_noBlank_t <- pathway_abund_clr_t[!grepl("Blank",rownames(pathway_abund_clr_t)),]
pathway_abund_clr_t_noBlank <- pathway_abund_noBlank_t[order(rownames(pathway_abund_noBlank_t)),]

pathway_corr_t <- data.frame(t(pathway_abund))
pathway_corr_t_no_blank <- pathway_corr_t[!grepl("Blank",rownames(pathway_corr_t)),]
pathway_corr_t_no_blank <- pathway_corr_t_no_blank[order(rownames(pathway_corr_t_no_blank)),]

which_one_to_remove <- c("Sp_282", "Sp_823")
pathway_abund_clr_t_noBlank <- pathway_abund_clr_t_noBlank[!rownames(pathway_abund_clr_t_noBlank)%in%which_one_to_remove,]
pathway_corr_t_no_blank <- pathway_corr_t_no_blank[!rownames(pathway_corr_t_no_blank)%in%which_one_to_remove,]
metadata_remove <- c("Sp_799", "Sp_827", "Sp_84", "Sp_856", "Sp_860", "Sp_862")

meta_data_sp_noBlank <- meta_data_sp_noBlank[!rownames(meta_data_sp_noBlank)%in%metadata_remove,]
meta_data_sp_noBlank <- meta_data_sp_noBlank[order(rownames(meta_data_sp_noBlank)),]
meta_data_sp_noBlank$pathway_div <- vegan::diversity(pathway_corr_t_no_blank, index="shannon")
meta_data_sp_noBlank$pathway_num <- vegan::specnumber(pathway_corr_t_no_blank)
meta_data_sp_noBlank$pathway_eveness <- meta_data_sp_noBlank$pathway_div / log(meta_data_sp_noBlank$pathway_num)
meta_data_sp_noBlank$gender_num <- ifelse(meta_data_sp_noBlank$gender=="W", 0, 1)
meta_data_sp_noBlank_cor <- select(meta_data_sp_noBlank, 
                                   c(age_years, BMI, gender_num, 
                                     shannon_div_phage, shannon_rare_bac, shannon_core_bac,
                                     specnumber_rare_bac, specnumber_core_bac, pathway_div, pathway_num))

colnames(meta_data_sp_noBlank_cor) <- c("Age [years]", "BMI", "Gender", 
                                        "Shannon diversity (phages)",
                                        "Shannon diversity (rare bacteria)", 
                                        "Shannon diversity (core bacteria)",
                                        "Species number (rare bacteria)",
                                        "Species number (core bacteria)", 
                                        "Shannon diversity (metabolic pathways)", 
                                        "Metabolic pathway number")

cor_all <- cor(meta_data_sp_noBlank_cor, method = "spearman")

testRes = corrplot::cor.mtest(meta_data_sp_noBlank_cor, conf.level = 0.95, method="spearman")
pAdj <- p.adjust(testRes$p, method = "BH")
resAdj <- matrix(pAdj, ncol = dim(testRes[[1]])[1])
colnames(resAdj) <- colnames(testRes$p)
rownames(resAdj) <- colnames(testRes$p)
pdf("corrplot.pdf", width = 7, height = 7)
corrplot::corrplot(cor_all, 
                   order = 'hclust', 
                   method = 'square',
                   type = 'upper', 
                   p.mat = resAdj, sig.level = 0.05, insig='blank', diag = FALSE,
                   tl.cex = 0.8, tl.col = "black")
dev.off()

###########################################################################################
# Prapare random forest

rare_bacteria <- data.frame(t(new_df_rare_bac))
core_bacteria <- data.frame(t(new_df_abundant_bac))
rare_phages <- data.frame(t(new_df_rare_pha))
core_phages <- data.frame(t(new_df_abundant_pha))

all_bacteria <- data.frame(cbind(rare_bacteria, core_phages))
all_bacteria <- all_bacteria[!rownames(all_bacteria)%in%metadata_remove,]
all_bacteria <- all_bacteria[order(rownames(all_bacteria)),]
all_taxonomy_function <- data.frame(all_bacteria, pathway_abund_clr_t_noBlank)
all_taxonomy_function$age_group <- NULL
all_taxonomy_function$BMI <- meta_data_sp_noBlank$BMI
all_taxonomy_function$age <- meta_data_sp_noBlank$age_years
all_taxonomy_function$shannon_bac <- meta_data_sp_noBlank$shannon_bac
all_taxonomy_function$shannon_div_phage <- meta_data_sp_noBlank$shannon_div_phage
all_taxonomy_function$evenness_bac <- meta_data_sp_noBlank$evenness_bac
all_taxonomy_function$evenness_phage <- meta_data_sp_noBlank$evenness_phage
all_taxonomy_function$physical_activity <- factor(meta_data_sp_noBlank$physical_activity, levels= c("1", "2", "3"))

set.seed(1)
random_seeds <- sample(1:200, 20, replace = FALSE) # repeat random forest with 100 different seeds set
imp_final = NULL

for (i in random_seeds){
  set.seed(i)
  boruta_rf <- Boruta(BMI~., data = all_taxonomy_function, pValue = 0.05, mcAdj=TRUE, maxRuns=1000)
  boruta_rf <- as.data.frame(boruta_rf$finalDecision)
  boruta_rf$pathway <- rownames(boruta_rf)
  boruta_rf_subset <-boruta_rf[boruta_rf$`boruta_rf$finalDecision`!="Rejected",]
  
  pathway_list <- c(boruta_rf_subset$pathway, "BMI")
  all_taxonomy_function_2 <- all_taxonomy_function[,colnames(all_taxonomy_function)%in%pathway_list]
  final_rf <- randomForest(BMI ~ ., data=all_taxonomy_function_2, na.action = na.omit, ntree=15, mtry=4, importance=TRUE) # 84, 80
  
  oob_MSE <- sqrt(tail (final_rf$mse, 1)) # oob MSE
  importance_df_rf <- data.frame(final_rf$importance)
  importance_df_rf$Predictors <- rownames(importance_df_rf) 
  rownames(importance_df_rf) <- NULL
  importance_df_rf$oob_MSE <- oob_MSE
  importance_df_rf$seed <- i
  importance_df_rf$rsq <- tail(final_rf$mse, n=1)
  importance_df_rf$var_explained <- tail(final_rf$rsq, n=1)
  # store table globally
  imp_final <- rbind(imp_final, importance_df_rf)
}


median(imp_final$var_explained) # 0.31
median(imp_final$oob_MSE) # 3.2
mad(imp_final$oob_MSE) # 0.13

imp_final$type <- ifelse(imp_final$Predictors == "age", "Host variables & Microbial ecosystem",
                         ifelse(imp_final$Predictors == "shannon_div_phage", "Host variables & Microbial ecosystem",
                                ifelse(grepl("PWY", imp_final$Predictors), "Metabolic pathway abundance", 
                                       "Host variables & Microbial ecosystem")))
imp_final <- imp_final[imp_final$X.IncMSE>0,]

inc_node_pur_plot <-
  ggplot(imp_final, aes(x=IncNodePurity,y=reorder(Predictors, IncNodePurity))) +
  geom_point(size=0.6) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, colour="red") +
  ylab(" ") +
  facet_wrap(~type, nrow=3, scales = "free_y") + theme_bw() +
  scale_y_discrete(labels=c("PWY.7992"="Menaquinol-8 biosynthesis III",
                            "PWY.1269"="CMP-3-deoxy-D-manno-octulosonate biosynthesis",
                            "RIBOSYN2.PWY"="Flavin biosynthesis I",
                            "DTDPRHAMSYN.PWY"="dTDP-beta-L-rhamnose biosynthesis",
                            "PWY.4041"="Gamma-glutamyl cycle",
                            "NONMEVIPP.PWY"="Methylerythritol phosphate pathway I",
                            "PWY.5855"="Ubiquinol-7 biosynthesis",
                            "PWY.5189"="Tetrapyrrole biosynthesis II",
                            "PWY66.409"="Superpathway of purine nucleotide salvage",
                            "PANTOSYN.PWY"="Superpathway of coenzyme A biosynthesis I",
                            "PWY0.1297"="Superpathway of purine deoxyribonucleosides degradation",
                            "PYRIDNUCSYN.PWY"="NAD de novo biosynthesis I",
                            "PWY.6969"="TCA cycle V",
                            "PWY.7953"="UDP-N-acetylmuramoyl-pentapeptide biosynthesis III",
                            "COLANSYN.PWY"="Colanic acid building blocks biosynthesis",
                            "PWY0.1296"="Purine ribonucleosides degradation",
                            "PANTO.PWY"="Phosphopantothenate biosynthesis I",
                            "Bacteroides.cellulosilyticus"=expression(italic("Bacteroides cellulosilyticus")),
                            "Megasphaera.elsdenii"=expression(italic("Megasphaera elsdenii")),
                            "Treponema.denticola"=expression(italic("Treponema denticola")),
                            "Blautia.hansenii"=expression(italic("Blautia hansenii")),
                            "shannon_div_phage"="Shannon diversity (Bacteriophages)",
                            "Bacteroides.ovatus"=expression(italic("Bacteroides ovatus")),
                            "Achromobacter.phage"="Achromobacter phage",
                            "Actinomyces.gaoshouyii"=expression(italic("Actinomyces gaoshouyii")),
                            "Selenomonas.sputigena"=expression(italic("Selenomonas sputigena")),
                            "Bacteroides.vulgatus"=expression(italic("Bacteroides vulgatus")),
                            "Alistipes.finegoldii"=expression(italic("Alistipes finegoldii")),
                            "Odoribacter.splanchnicus"=expression(italic("Odoribacter splanchnicus")),
                            "Filifactor.alocis"=expression(italic("Filifactor alocis")),
                            "Roseburia.hominis"=expression(italic("Roseburia hominis")),
                            "Bacteroides.dorei"=expression(italic("Bacteroides dorei")),
                            "Bacteroides.caecimuris"=expression(italic("Bacteroides caecimuris")),
                            "age"="Age [years]"
  )) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white")) + xlab("Mean Decrease Gini (IncNodePurity)")

lncMSE_plot <-
  ggplot(imp_final, aes(x=X.IncMSE,y=reorder(Predictors, X.IncMSE))) +
  geom_point(size=0.6) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, colour="red") +
  ylab(" ") + xlab("Mean Decrease Accuracy (%lncMSE)") +
  facet_wrap(~type, nrow=3, scales = "free_y") + theme_bw() +
  scale_y_discrete(labels=c("PWY.7992"="Menaquinol-8 biosynthesis III",
                            "PWY.1269"="CMP-3-deoxy-D-manno-octulosonate biosynthesis",
                            "RIBOSYN2.PWY"="Flavin biosynthesis I",
                            "DTDPRHAMSYN.PWY"="dTDP-beta-L-rhamnose biosynthesis",
                            "PWY.4041"="Gamma-glutamyl cycle",
                            "NONMEVIPP.PWY"="Methylerythritol phosphate pathway I",
                            "PWY.5855"="Ubiquinol-7 biosynthesis",
                            "PWY.5189"="Tetrapyrrole biosynthesis II",
                            "PWY66.409"="Superpathway of purine nucleotide salvage",
                            "PANTOSYN.PWY"="Superpathway of coenzyme A biosynthesis I",
                            "PWY0.1297"="Superpathway of purine deoxyribonucleosides degradation",
                            "PYRIDNUCSYN.PWY"="NAD de novo biosynthesis I",
                            "PWY.6969"="TCA cycle V",
                            "PWY.7953"="UDP-N-acetylmuramoyl-pentapeptide biosynthesis III",
                            "COLANSYN.PWY"="Colanic acid building blocks biosynthesis",
                            "PWY0.1296"="Purine ribonucleosides degradation",
                            "PANTO.PWY"="Phosphopantothenate biosynthesis I",
                            "Bacteroides.cellulosilyticus"=expression(italic("Bacteroides cellulosilyticus")),
                            "Megasphaera.elsdenii"=expression(italic("Megasphaera elsdenii")),
                            "Treponema.denticola"=expression(italic("Treponema denticola")),
                            "Blautia.hansenii"=expression(italic("Blautia hansenii")),
                            "shannon_div_phage"="Shannon diversity (Bacteriophages)",
                            "Bacteroides.ovatus"=expression(italic("Bacteroides ovatus")),
                            "Achromobacter.phage"="Achromobacter phage",
                            "Actinomyces.gaoshouyii"=expression(italic("Actinomyces gaoshouyii")),
                            "Selenomonas.sputigena"=expression(italic("Selenomonas sputigena")),
                            "Bacteroides.vulgatus"=expression(italic("Bacteroides vulgatus")),
                            "Alistipes.finegoldii"=expression(italic("Alistipes finegoldii")),
                            "Odoribacter.splanchnicus"=expression(italic("Odoribacter splanchnicus")),
                            "Filifactor.alocis"=expression(italic("Filifactor alocis")),
                            "Roseburia.hominis"=expression(italic("Roseburia hominis")),
                            "Bacteroides.dorei"=expression(italic("Bacteroides dorei")),
                            "Bacteroides.caecimuris"=expression(italic("Bacteroides caecimuris")),
                            "age"="Age [years]"
  )) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"))


most_imp_random_forest_plots <- ggarrange(inc_node_pur_plot, lncMSE_plot, labels = c("A", "B"))

bmi_table <- names(table(imp_final$Predictors))
bmi_species_df <- data.frame(Predictors=bmi_table)
bmi_species_df$Predictors
bmi_species_df$Type <- c("Bacteriophage-associated", 
                         "Gastrointestinal tract-associated microbes", 
                         "Host and diversity indices", 
                         "Gastrointestinal tract-associated microbes",
                         "Gastrointestinal tract-associated microbes", 
                         "Gastrointestinal tract-associated microbes", 
                         "Gastrointestinal tract-associated microbes",
                         "Gastrointestinal tract-associated microbes", 
                         "Gastrointestinal tract-associated microbes", 
                         "Gastrointestinal tract-associated microbes",
                         "Carbohydrate Biosynthesis", 
                         "Carbohydrate Biosynthesis", 
                         "Respiratory tract-associated microbes",
                         "Gastrointestinal tract-associated microbes", 
                         "Terpenoid Biosynthesis", 
                         "Gastrointestinal tract-associated microbes", 
                         "Cofactor, Carrier, and Vitamin Biosynthesis", 
                         "Cofactor, Carrier, and Vitamin Biosynthesis",
                         "Carbohydrate Biosynthesis", 
                         "Cofactor, Carrier, and Vitamin Biosynthesis", 
                         "Tetrapyrrole Biosynthesis",
                         "Cofactor, Carrier, and Vitamin Biosynthesis", 
                         "Generation of Precursor Metabolites and Energy",
                         "Cell Structure Biosynthesis", 
                         "Cofactor, Carrier, and Vitamin Biosynthesis",
                         "Degradation/Utilization/Assimilation", 
                         "Degradation/Utilization/Assimilation",
                         "Purine Nucleotide Biosynthesis", 
                         "Cofactor, Carrier, and Vitamin Biosynthesis", 
                         "Cofactor, Carrier, and Vitamin Biosynthesis", 
                         "Gastrointestinal tract-associated microbes", 
                         "Respiratory tract-associated microbes", 
                         "Bacteriophage-associated", 
                         "Respiratory tract-associated microbes")
bmi_species_df$Species_type <- ifelse(bmi_species_df$Predictors %in% colnames(rare_bacteria), "low-abundance bacteria",
                                      ifelse(bmi_species_df$Predictors %in% colnames(rare_bacteria), "high-abundance bacteria",
                                             ifelse(bmi_species_df$Predictors %in% colnames(rare_phages), "low-abundance bacteriophages", NA)))
predictor_list <- c(bmi_species_df$Type)
names(predictor_list) <- bmi_species_df$Predictors

list1 = c()
for (items in imp_final$Predictors){
  list1 = append(list1,predictor_list[[items]])}

imp_final$Predictors_2 <- list1
imp_final$Predictors_3 <- ifelse(grepl("PWY",imp_final$Predictors), "Metabolic pathway abundances", imp_final$Predictors_2)
imp_final_2 <- ddply(imp_final, .(Predictors_3, seed), numcolwise(sum))
imp_final_3 <- imp_final
imp_final_3$Species_type <- ifelse(imp_final_3$Predictors %in% colnames(rare_bacteria), "Low-abundance bacteria",
                                   ifelse(imp_final_3$Predictors %in% colnames(core_bacteria), "High-abundance bacteria",
                                          ifelse(imp_final_3$Predictors %in% colnames(rare_phages), "Low-abundance bacteriophages", 
                                                 ifelse(imp_final_3$Predictors %in% colnames(core_phages), "High-abundance bacteriophages", 
                                                        ifelse(grepl("PWY", imp_final_3$Predictors), "Metabolic pathway abundances", "Host and diversity indices")))))

imp_final_3 <- ddply(imp_final_3, .(Species_type, seed), numcolwise(sum))


lncMSE_plot_merged <-
  ggplot(imp_final_2, aes(x=X.IncMSE,y=reorder(Predictors_3, X.IncMSE))) +
  geom_point(size=0.6) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, colour="red") +
  ylab(" ") + xlab("%lncMSE") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"))


IncNodePurity_plot_merged <-
  ggplot(imp_final_2, aes(x=IncNodePurity,y=reorder(Predictors_3, IncNodePurity))) +
  geom_point(size=0.6) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, colour="red") +
  ylab(" ") + xlab("IncNodePurity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"))


lncMSE_plot_merged_rare <-
  ggplot(imp_final_3, aes(x=X.IncMSE,y=reorder(Species_type, X.IncMSE))) +
  geom_point(size=0.6) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, colour="red") +
  ylab(" ") + xlab("%lncMSE") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"))

IncNodePurity_plot_merged_rare <-
  ggplot(imp_final_3, aes(x=IncNodePurity,y=reorder(Species_type, IncNodePurity))) +
  geom_point(size=0.6) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, colour="red") +
  ylab(" ") + xlab("IncNodePurity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"))


small_plots <- ggarrange(IncNodePurity_plot_merged, lncMSE_plot_merged,
                         IncNodePurity_plot_merged_rare, lncMSE_plot_merged_rare, nrow=1, ncol=4, labels=c("C", "D", "E", "F"))
final_rf_plot <- ggarrange(most_imp_random_forest_plots, small_plots, nrow=2, heights = c(1,0.3))

# Figure 07
pdf("random_forest.pdf", width = 16, height = 7)
final_rf_plot
dev.off()
