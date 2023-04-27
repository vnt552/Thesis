### metaDMG plots
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

metaDMG_results <- read_csv("metaDMG/results.csv.gz")


### Filter reads, above 50
metaDMG_results_50 <- metaDMG_results%>%  filter(N_reads>=50)




# Filter the reads by classification
metaDMG_superkingdom <- metaDMG_results_50 %>%  filter(tax_rank=="superkingdom")
metaDMG_kingdom <- metaDMG_results_50 %>%  filter(tax_rank=="kingdom")

### Make abundance plots over superkingdom
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
  ggtitle("Percentage Classified Reads", subtitle = "Superkingdom") +
  coord_flip()


### Make abundance plots over Eukaryotes kingdom
# First calculate the percentage of each
euka_results <- read.csv("metaDMG/all_samples_eukaryota.csv", header = FALSE)
newheader <- colnames(metaDMG_results)
colnames(euka_results) <- newheader
euka_kingdom <- euka_results%>%  filter(tax_rank == "kingdom")
euka_kingdom_50 <- euka_kingdom%>%  filter(N_reads>=50)

euka_kingdom_percentage <- euka_kingdom_50  %>%
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
  ggtitle("Percentage Classified Reads", subtitle = "Eukaryote") +
  coord_flip()

## Can't see vira and archea, so now make abundance plot only with those two. 
#archaea_results <- read.csv("metaDMG/all_samples_archaea.csv", header = FALSE)
#virae_results <- read.csv("metaDMG/all_samples_virae.csv", header = FALSE)
#colnames(archaea_results) <- newheader
#colnames(virae_results) <- newheader

virae_archaea_percentage <- superkingdom_percentage[superkingdom_percentage$tax_name %in% c('Viruses', 'Archaea'), ]
virae_archaea_percentage_plot <- ggplot(virae_archaea_percentage, aes(y=percentage, x = sample, fill=tax_name)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  xlab("Sample") +
  ylab("% Classified Reads") +
  labs(fill="Superkingdom") +
  ggtitle("Percentage Classified Reads", subtitle = "Virae and Archaea") +
  coord_flip()


################### Number of reads per sample
# to be used with the Pb figure
# Before filtering the data N_reads > 50
total_reads_root <- metaDMG_results_50 %>%  filter(tax_name=="root")
total_reads_root$sample <- factor(total_reads_root$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))


classified_reads <- ggplot(total_reads_root, aes(y = N_reads, x = sample, group=1)) +
  geom_line() +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of Classified Reads") +
  ggtitle("Number of Classified Reads") +
  coord_flip()
classified_reads




############################

read_results <- read_csv("metaDMG/Bat_reads_workflow.csv")
read_results$Sample <- factor(read_results$Sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))

ggplot(read_results, aes(x = Sample, y = aMAW, group=1)) +
  geom_line(color = "blue")+
  theme_minimal()+
  xlab("Sample") +
  ylab("aMAW") +
  coord_flip()

ggplot(read_results, aes(x = Sample, y = metaDMG, group=1)) +
  geom_line(color= "green")+
  theme_minimal()+
  xlab("Sample") +
  ylab("metaDMG") +
  coord_flip()

ggplot(read_results, aes(x = Sample, y = Raw_data, group=1)) +
  geom_line(color = "red")+
  theme_minimal()+
  xlab("Sample") +
  ylab("Raw_data") +
  coord_flip()

## Top organisms in bacteria and metazoa
#
bacteria <- read.csv("metaDMG/all_samples_bacteria.csv", header=FALSE)
colnames(bacteria) <- newheader
bacteria$sample <- factor(bacteria$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))
bacteria_50 <-bacteria%>%  filter(N_reads>=50)

## only reads classified to order level
bacteria_family <- bacteria_50 %>%  filter(tax_rank=="family")

## Find top 10 family in each sample
bacteria_top <- bacteria_family[order(bacteria_family$N_reads, decreasing = TRUE), ]
bacteria_top <- Reduce(rbind,                                 # Top N highest values by group
                    by(bacteria_top,
                       bacteria_top["sample"],
                       head,
                       n = 10))
bacteria_top <- bacteria_top[order(bacteria_top$N_reads, decreasing = TRUE), ]


## metazoa
metazoa <- read.csv("metaDMG/all_samples_metazoa.csv", header=FALSE)
colnames(metazoa) <- newheader

metazoa$sample <- factor(metazoa$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))
metazoa_50 <-metazoa%>%  filter(N_reads>=50)

## only reads classified to order level
metazoa_family <- metazoa_50 %>%  filter(tax_rank=="genus")

######################## heatmap for abundance try 2
metazoa_abundance <- metazoa_family  %>%
  group_by(sample, tax_name) %>%
  summarise(N_reads = sum(N_reads)) %>%
  mutate(RelativeAbundance = N_reads / sum(N_reads))

metazoa_top_abundance <- metazoa_abundance[order(metazoa_abundance$N_reads, decreasing = TRUE), ]
metazoa_top_abundance <- Reduce(rbind,                                 # Top N highest values by group
                      by(metazoa_top_abundance,
                         metazoa_top_abundance["sample"],
                         head,
                         n = 10))
metazoa_top_abundance<- metazoa_top_abundance %>% select(tax_name, sample, RelativeAbundance)
metazoa_top_abundance_wide <-metazoa_top_abundance %>% pivot_wider(names_from = sample, values_from = RelativeAbundance)
metazoa_top_abundance_header <- as.matrix(metazoa_top_abundance_wide[,-1])
rownames(metazoa_top_abundance_header) <-as.matrix(metazoa_top_abundance_wide[,1])
metazoa_top_abundance_header[is.na(metazoa_top_abundance_header)] <- 0

pheatmap(metazoa_top_abundance_header, color = hcl.colors(50, "BluYl"))

### this heatmap only shows the highest values, and many looks completly empty. Therefore scaling for all data
metazoa_matrix_scaled <- t(scale(t(metazoa_top_abundance_header)))
metazoa_scaled_heatmap_abundance <- pheatmap(metazoa_matrix_scaled, color = hcl.colors(50, "BluYl"))
# scaling per row
pheatmap(metazoa_top_abundance_header, scale = "row", color = hcl.colors(50, "BluYl"))

######### Bacteria Heatmap

bacteria_abundance <- bacteria_family  %>%
  group_by(sample, tax_name) %>%
  summarise(N_reads = sum(N_reads)) %>%
  mutate(RelativeAbundance = N_reads / sum(N_reads))

bacteria_top_abundance <- bacteria_abundance[order(bacteria_abundance$N_reads, decreasing = TRUE), ]
bacteria_top_abundance <- Reduce(rbind,                                 # Top N highest values by group
                                by(bacteria_top_abundance,
                                   bacteria_top_abundance["sample"],
                                   head,
                                   n = 10))
bacteria_top_abundance<- bacteria_top_abundance %>% select(tax_name, sample, RelativeAbundance)
bacteria_top_abundance_wide <-bacteria_top_abundance %>% pivot_wider(names_from = sample, values_from = RelativeAbundance)
bacteria_top_abundance_header <- as.matrix(bacteria_top_abundance_wide[,-1])
rownames(bacteria_top_abundance_header) <-as.matrix(bacteria_top_abundance_wide[,1])
bacteria_top_abundance_header[is.na(bacteria_top_abundance_header)] <- 0

pheatmap(bacteria_top_abundance_header, color = hcl.colors(50, "BluYl"))

### this heatmap only shows the highest values, and many looks completly empty. Therefore scaling for all data
bacteria_matrix_scaled <- t(scale(t(bacteria_top_abundance_header)))
bacteria_scaled_heatmap_abundance <- pheatmap(bacteria_matrix_scaled, color = hcl.colors(50, "BluYl"))
# scaling per row
#pheatmap(metazoa_top_abundance_header, scale = "row", color = hcl.colors(50, "BluYl"))



##############The top organisms in each sample, facet wrap plot
metazoa_genus <- metazoa_50 %>%  filter(tax_rank=="genus")
metazoa_top_genus <- metazoa_genus[order(metazoa_genus$N_reads, decreasing = TRUE), ]
metazoa_top_genus <- Reduce(rbind,                                 # Top N highest values by group
                      by(metazoa_top_genus,
                         metazoa_top_genus["sample"],
                         head,
                         n = 10))
metazoa_top_genus <- metazoa_top_genus[order(metazoa_top_genus$N_reads, decreasing = TRUE), ]

ggplot(metazoa_top) +
  geom_col(aes(reorder_within(tax_name, - N_reads, sample), N_reads)) +
  scale_x_reordered() +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  facet_wrap(~sample, scales = "free", ncol = 4) +
  xlab("") +
  ylab("Number of Classified Reads") +
  ggtitle("Top 10 Organisms from each sample, Family Level") +
  coord_flip()


##### Looking at how chiroptera changes
chiroptera_results <- read.csv("metaDMG/all_samples_chiroptera.csv", header=FALSE)
colnames(chiroptera_results) <- newheader

chiroptera_results_50 <- chiroptera_results %>%  filter(N_reads > 50)
chiroptera_species <- chiroptera_results_50 %>% filter(tax_rank == "species")
chiroptera_genus <- chiroptera_results_50 %>%  filter(tax_rank == "genus")
chiroptera_genus$sample <- factor(chiroptera_genus$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))


chiroptera_plot_genus <- ggplot(chiroptera_genus, aes(x= sample, y = N_reads, group = tax_name, color = tax_name)) + 
  geom_line() +
  geom_point() +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of Reads") +
  labs(color="Genus") +
  ggtitle("Number of Reads Classified to Chiroptera Across Samples", subtitle = "Genus Level") +
  coord_flip()

## Plots ##################################################################
superkingdom_percentage_plot
virae_archaea_percentage_plot
euka_kingdom_percentage_plot
classified_reads
metazoa_scaled_heatmap_abundance
bacteria_scaled_heatmap_abundance
chiroptera_plot_genus
chiroptera_plot_species

ggdraw() +
  draw_plot(classified_reads, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(superkingdom_percentage_plot, x = .5, y = 0.5, width = .5, height = .5) +
  draw_plot(virae_archaea_percentage_plot, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(euka_kingdom_percentage_plot, x = 0.5, y = 0, width = .5, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))


ggdraw() +
  draw_plot(superkingdom_percentage_plot, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(virae_archaea_percentage_plot, x = .5, y = 0.5, width = .5, height = .5) +
  draw_plot(euka_kingdom_percentage_plot, x = 0, y = 0, width = .5, height = .5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))

