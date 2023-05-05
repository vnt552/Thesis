# Plots on the preliminary data
library(ggplot2)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(tidyr)
library(superheat)
library(ggforce)
library(viridis)
library(devtools)
library(ggpubr)
library(RColorBrewer)
library(cowplot)
library(upstartr)

# Load data
metaDMG_results <- read_csv("metaDMG/results.csv.gz")
metaDMG_results_50 <- metaDMG_results%>%  filter(N_reads>=50) # Filter reads, above 50

# Filter the reads by classification
metaDMG_superkingdom <- metaDMG_results_50 %>%  filter(tax_rank=="superkingdom")


######### Make abundance plots over superkingdom ###############
# First calculate the percentage of each

superkingdom_percentage <- metaDMG_superkingdom  %>%
  group_by(sample, tax_name) %>%
  summarise(n = sum(N_reads)) %>%
  mutate(percentage = n / sum(n))

# Place the samples in the rigth order
superkingdom_percentage$sample <- factor(superkingdom_percentage$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))

# make the plot
superkingdom_percentage_plot <- ggplot(superkingdom_percentage, aes(y=percentage, x = sample, fill=tax_name)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  xlab("Sample") +
  ylab("% Classified Reads") +
  labs(fill="Superkingdom") +
  ggtitle("Relative Abundance", subtitle = "Superkingdom") +
  coord_flip()
  
  
  
euka_kingdom <- metaDMG_results_50 %>%  filter(tax_rank == "kingdom", grepl("\\bEukaryote\\b", tax_path))

euka_kingdom_percentage <- euka_kingdom %>%
  group_by(sample, tax_name) %>%
  summarise(n = sum(N_reads)) %>%
  mutate(percentage = n / sum(n))

# Place the samples in the rigth order
euka_kingdom_percentage$sample <- factor(euka_kingdom_percentage$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))

# make the plot
euka_kingdom_percentage_plot <- ggplot(euka_kingdom_percentage, aes(y=percentage, x = sample, fill=tax_name)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  xlab("Sample") +
  ylab("% Classified Reads") +
  labs(fill="Kingdoms Eukaryote") +
  ggtitle("Relative Abundance", subtitle = "Eukaryote") +
  coord_flip()
 
