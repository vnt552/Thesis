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

metaDMG_results <- read_csv("metaDMG/results.csv.gz")

#### Basic ####
#Number of OTU allocated to phylum
num_otu <- metaDMG_results %>%  filter(tax_rank=="species")
mean(num_otu$N_reads)

num_otu_list <- as.list(num_otu$tax_id)
length(unique(num_otu_list)) ### number of unique species the reads have been classified as. 
sum(num_otu$N_reads) # Number of reads classified to species level




### How many reads in total
# Use "roots" or "num_otu" depending on the total number of classified reads or only OTU(species level)
sum(roots$N_reads)
mean(roots$N_reads)
sd(roots$N_reads)


### Mean and SD per sample
sum(num_otu$N_reads)
mean(num_otu$N_reads)
sd(num_otu$N_reads)
by(num_otu$N_reads,num_otu$sample,mean)
by(num_otu$N_reads,num_otu$sample,sd)

### how many unique species are found in each sample?
length(unique(num_otu_list)) ### number of unique species the reads have been classified as. 
species_count <- count(num_otu$sample)
mean(species_count$freq)
sd(species_count$freq)

# Number singletons (organisms found with only one read)
singletons <- num_otu %>%  filter(N_reads==1)
singleton_list <- as.list(singletons$tax_id)
length(unique(singleton_list))

### Number OTU after filtering 50
num_otu_50 <- num_otu %>%  filter(N_reads>=50)

num_otu_list_50 <- as.list(num_otu_50$tax_id)
length(unique(num_otu_list_50))
sum(num_otu_50$N_reads)
mean(num_otu_50$N_reads)
sd(num_otu_50$N_reads)



### Filter data by level
metaDMG_genus <- metaDMG_results %>% filter(tax_rank=="genus")
metaDMG_species <- metaDMG_results %>%  filter(tax_rank=="species")
metaDMG_kingdom <- metaDMG_results %>%  filter(tax_rank=="kingdom")

### The order of the samples, from newest to oldest or the other way around
metaDMG_species$sample <- factor(metaDMG_species$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))
metaDMG_genus$sample <- factor(metaDMG_genus$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))
metaDMG_kingdom$sample <- factor(metaDMG_kingdom$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))
metaDMG_results$sample <- factor(metaDMG_results$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))
metaDMG_kingdom$sample <- factor(metaDMG_kingdom$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))
metaDMG_results$sample <- factor(metaDMG_results$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))
metaDMG_genus$sample <- factor(metaDMG_genus$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))

roots <- metaDMG_results %>% filter(tax_id==1)

### Plot 1 ###
##Længde plots
ggplot(roots, aes(x=sample, y=mean_L)) +
  geom_col() +
  theme_bw() +
  xlab("Sample") +
  ylab("Mean read length (bp)") +
  ggtitle("Mean read length per sample") 

classified_reads <- ggplot(roots, aes(x=sample, y=N_reads)) +
  geom_col() +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of Reads") +
  ggtitle("Number of Classified Reads") +
  coord_flip()

### Plot 2 ###
##Længde plots for alle samples facet wrap
length_of_all_with_trend <- ggplot(metaDMG_results, aes(x=mean_L)) +
  geom_histogram(binwidth=5) +
  geom_density(aes(y=5*..count..), colour="red", adjust=4) +
  facet_wrap(~sample, scales = "free") +
  theme_bw() +
  xlab("Mean length of reads (bp)") +
  ylab("Amount of reads") +
  ggtitle("Mean read length per sample") 




### Plot 3 ###
#Number of reads across samples, kingdom level
## First calculate percentage, then plot
data_1 <- all_kingdom %>%  filter(tax_rank =="order")
data <- data  %>%
  group_by(sample, kingdom) %>%
  summarise(n = sum(N_reads)) %>%
  mutate(percentage = n / sum(n))

ggplot(data) + geom_col(aes(x=percentage, y=sample))

percentage_kingdoms <- ggplot(data, aes(y=percentage, x = sample, fill=kingdom)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  xlab("Sample") +
  ylab("Percentage of Reads Classified to Order Level") +
  ggtitle("Percentage of Reads Classified to each Kingdom") +
  coord_flip()


### Plot 4 ###
## Mean read length of reads mapped to kingdom
ggplot(metaDMG_kingdom, aes(y=mean_L, x = sample, group=tax_name, col=tax_name)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  theme_bw() +
  xlab("Sample") +
  ylab("Mean read Length") +
  ggtitle("Mean Length for reads mapped to Kingdom") +
  coord_flip()


### Plot 5 ###
## How many reads are there per tax_name, species level
ggplot(metaDMG_species, aes(x=sample, y=N_reads)) +
  geom_violin() +
  geom_sina() +
  scale_y_log10() +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of Mapped Reads") +
  ggtitle("Number of mapped reads to each tax ID, species level")


### Plot 6 ###
## Number of taxa mapped in each sample
taxa_in_samples <- ggplot(metaDMG_species) +
  geom_bar(aes(x = sample)) +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of species mapped to") +
  ggtitle("Number of species mapped to per sample")


### Plot 7 ###
## plot damage vs. significance for species
damage_plot <- ggplot(metaDMG_species) +
  geom_point(aes(x = MAP_significance, y = MAP_damage)) +
  theme_bw() +
  xlab("Significance") +
  ylab("Damage") +
  ggtitle("Damage and significance for species")


##################### Del data op i mindre filer alt efter organisme ############
metazoa <- metaDMG_results %>%  filter(grepl("\\bMetazoa\\b", tax_path))
fungi <- metaDMG_results %>%  filter(grepl("\\bFungi\\b", tax_path))
viridiplants <- metaDMG_results %>%  filter(grepl("\\bViridiplantae\\b", tax_path))
bacteria <- metaDMG_results %>%  filter(grepl("\\bBacteria\\b", tax_path))
virae <- metaDMG_results %>%  filter(grepl("\\bViruses\\b", tax_path))
archaea <- metaDMG_results %>%  filter(grepl("\\bArchaea\\b", tax_path))

### Lav fil, hvor alle organismer har deres kingdom stående
metazoa$kingdom <- "Metazoa"
viridiplants$kingdom <- "Viridiplantae"
bacteria$kingdom <- "Bacteria"
fungi$kingdom <- "Fungi"
virae$kingdom <- "Virae"
archaea$kingdom <- "Archaea"
all_kingdom <- rbind(metazoa, viridiplants, bacteria, fungi, virae, archaea)
all_kingdom$sample <- factor(all_kingdom$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))
all_kingdom$sample <- factor(all_kingdom$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))

#### Test at all_kingdom fil indeholder alt data fra metaDMG_results
metaDMG_results%>%select(which(!(colnames(metaDMG_results) %in% colnames(all_kingdom))))
a1 <- metaDMG_results[,1:5]
a2 <- all_kingdom[,1:5]
in_results_not_allkingdom <- setdiff(a1,a2)



### Sæt rækkefølgen på samples 
metazoa$sample <- factor(metazoa$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))
metazoa$sample <- factor(metazoa$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))

viridiplants$sample <- factor(viridiplants$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))
bacteria$sample <- factor(bacteria$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))
fungi$sample <- factor(fungi$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))


### Plot 8 ###
## Species in metazoa
metazoa_genus <- metazoa %>% filter(tax_rank=="genus")
metazoa_order <- metazoa %>% filter(tax_rank=="order")
metazoa_order <- metazoa_order%>%  filter(N_reads>=100)

ggplot(metazoa_genus, aes(y=N_reads, x = sample)) +
  geom_col() +
  scale_y_log10() +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of Reads") +
  ggtitle("Number of Reads Mapped to Metazoa, Genus Level")



## Top 100 organisms found in each sample
all_kingdom_genus <- all_kingdom %>% filter(tax_rank=="genus")
all_kingdom$sample <- factor(all_kingdom$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))
threshold50 <- all_kingdom_genus %>%  filter(N_reads>=50)
ggplot(threshold50) + geom_bar(aes(x=sample))
top_genus <- all_kingdom_genus[order(all_kingdom_genus$N_reads, decreasing = TRUE), ]
top_genus <- Reduce(rbind,                                 # Top N highest values by group
                    by(top_genus,
                       top_genus["sample"],
                       head,
                       n = 100))
top_genus <- top_genus[order(top_genus$N_reads, decreasing = TRUE), ]


top_GA11 <- top_genus %>% filter(sample=="GA-11")
top_GA19 <- top_genus %>% filter(sample=="GA-19")
top_GA31 <- top_genus %>% filter(sample=="GA-31")
top_GA41 <- top_genus %>% filter(sample=="GA-41")
top_GA61 <- top_genus %>% filter(sample=="GA-61")
top_GA71 <- top_genus %>% filter(sample=="GA-71")
top_GA81 <- top_genus %>% filter(sample=="GA-81")
top_GA91 <- top_genus %>% filter(sample=="GA-91")
top_GA111 <- top_genus %>% filter(sample=="GA-111")
top_GA121 <- top_genus %>% filter(sample=="GA-121")
top_GA131 <- top_genus %>% filter(sample=="GA-131")
top_GA141 <- top_genus %>% filter(sample=="GA-141")
top_GA161 <- top_genus %>% filter(sample=="GA-161")
top_GA181 <- top_genus %>% filter(sample=="GA-181")
top_GA191 <- top_genus %>% filter(sample=="GA-191")
top_GACex <- top_genus %>% filter(sample=="GA-Cex")
top_GACli <- top_genus %>% filter(sample=="GA-Cli")


## Plot top reads enkeltvis per sample ##
ggplot(top_GA81) +
  geom_col(aes(reorder(tax_name, - N_reads), N_reads)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  xlab("Genus") +
  ylab("Number of Reads Mapped") +
  ggtitle("Top 100 in N_reads, GA-81")

### Plot 9 ###
top80 <- ggplot(top_genus) +
  geom_col(aes(reorder_within(tax_name, - N_reads, sample), N_reads, fill=kingdom)) +
  scale_fill_brewer(palette="Set2") +
  scale_x_reordered() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")+ 
  facet_wrap(~sample, scales = "free", ncol = 4) +
  xlab("") +
  ylab("Number of Classified Reads") +
  ggtitle("Top 80 Organisms from each sample, Genus Level")

top80 


### Plot 10 ###
## Top 100 organismer i en sample, navn på genus og linje på 50
#### Calcualting abundance
all_genus_50 <- all_kingdom %>% filter(tax_rank=="genus") %>%  filter(N_reads >= 50)
all_genus_50 <- summarise(all_genus_50$N_reads, groups=sample)
all_genus_50 %>%group

ggplot(top_GA19) +
  geom_col(aes(reorder(tax_name, - N_reads), N_reads)) +
  #geom_hline(yintercept=100, color = "red") 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  xlab("Genus") +
  ylab("Number of Reads Mapped") +
  ggtitle("Top 100 in N_reads, GA-191")




########### PCR plot ################
### PCA med alt data, Genus ###
superking <- metaDMG_results %>%  filter(tax_rank == "superkingdom")
metaDMG_pca <- select_if(superking, is.numeric) 
metaDMG_header <- superking[1]
metaDMG_header <- metaDMG_header %>%  mutate(group = case_when(grepl("GA-C", sample) ~ "Control",
                                             grepl("1", sample) ~"Sample"))

str(metaDMG_pca)
dim(metaDMG_pca)
metaDMG_pca <- prcomp(metaDMG_pca, scale=TRUE, center = TRUE)

percent_variance_metaDMG <- 
  summary(metaDMG_pca)$importance["Proportion of Variance",] * 100
summary(metaDMG_pca)

# add header
plot_metaDMG <- 
  as_tibble(metaDMG_pca$x) %>% 
  bind_cols(metaDMG_header)

plot_metaDMG

###Make PCA plot

pca_control_sample <- ggplot(plot_metaDMG, aes(x=PC1, y=PC2, col = group)) + geom_point() + 
  xlab(label = paste("PC1", percent_variance_metaDMG[1])) +
  ylab(label = paste("PC2", percent_variance_metaDMG[2])) +theme_bw() +
  ggtitle("PCA Genus Level")

ggplot(plot_metaDMG, aes(x=PC3, y=PC4, col = group)) + geom_point() + 
  xlab(label = paste("PC3", percent_variance_metaDMG[3])) +
  ylab(label = paste("PC4", percent_variance_metaDMG[4])) +theme_bw() +
  ggtitle("PCA Genus Level")

pca_control_sample
#### PCA med metazoa, viridiplantae og fungi ###

kingdom_pca <- select_if(all_kingdom_genus, is.numeric) %>% as.data.frame()
kingdom_pca <- kingdom_pca[!names(kingdom_pca) %in% c("tax_id")]
kingdom_header <- cbind(all_kingdom_genus[1], all_kingdom_genus[128])
kingdom_header <- kingdom_header %>%  mutate(group = case_when(grepl("GA-C", sample) ~ "Control",
                                                               grepl("1", sample) ~"Sample"))
kingdom_pca <- kingdom_pca %>%  as.matrix()
kingdom_pca <- prcomp(t(kingdom_pca), scale. = TRUE, center = TRUE)

str(kingdom_pca)
dim(kingdom_pca)

percent_variance_kingdom <- 
  summary(kingdom_pca)$importance["Proportion of Variance",] * 100
summary(kingdom_pca)

# add header
plot_kingdom_pca <- 
  as_tibble(kingdom_pca$x) %>% 
  bind_cols(kingdom_header)

plot_kingdom_pca

#Make PCA plot

ggplot(plot_kingdom_pca, aes(x=PC1, y=PC2, col=sample)) + geom_point() + 
  xlab(label = paste("PC1", percent_variance_kingdom[1])) +
  ylab(label = paste("PC2", percent_variance_kingdom[2])) +theme_bw() +
  ggtitle("PCA Genus Level")

ggplot(plot_kingdom_pca, aes(x=PC3, y=PC4, col=group)) + geom_point() + 
  xlab(label = paste("PC3", percent_variance_kingdom[3])) +
  ylab(label = paste("PC4", percent_variance_kingdom[4])) +theme_bw() +
  ggtitle("PCA Genus Level")


top_genus_kingdom <- all_kingdom_genus[order(all_kingdom_genus$N_reads, decreasing = TRUE), ]
top_genus_kingdom <- Reduce(rbind,                                 # Top N highest values by group
                    by(top_genus_kingdom,
                       top_genus_kingdom["sample"],
                       head,
                       n = 100))
top_genus_kingdom <- top_genus_kingdom[order(top_genus_kingdom$N_reads, decreasing = TRUE), ]


top_GA181 <- top_genus_kingdom %>% filter(sample=="GA-181")


### Uden navn, til stort plot
ggplot(top_GA181) +
  geom_col(aes(reorder(tax_name, - N_reads), N_reads, fill=kingdom)) +
  geom_hline(yintercept=50, color = "red") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  xlab("Genus") +
  ylab("Number of Mapped Reads") +
  ggtitle("Top 100 Genus, N_reads, GA-181")

### Reads left after each threshold
reads_left <- read.csv("bat_reads_filter2.csv")
reads_left$Sample <- factor(reads_left$Sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))

ggplot(reads_left) +
  geom_col(aes(x=Sample,y=Reads, fill= Filter_type), position="dodge") +
  scale_y_log10() +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of Reads") +
  ggtitle("Number of reads left after mapping")
  



################### LCA MAPS #######################
## Heatmap for the top genus, devided by their kingdom
## First, calculate the abundance in percentage. 
## Abundance is calculated by the sample group, as N_reads for taxa / Total N_reads per sample
top_20_genus <- all_kingdom_genus[order(all_kingdom_genus$N_reads, decreasing = TRUE), ]
top_20_genus <- Reduce(rbind,                                 # Top N highest values by group
                    by(top_20_genus,
                       top_20_genus["sample"],
                       head,
                       n = 20))
top_20_genus <- top_20_genus[order(top_20_genus$N_reads, decreasing = TRUE), ]
heatmap_top_metazoa <- top_20_genus %>%  filter(kingdom=="Metazoa")
heatmap_top_metazoa <-  transform(heatmap_top_metazoa,                             # Calculate percentage by group
                                  perc = ave(N_reads,
                                             sample,
                                             FUN = prop.table))
heatmap_metazoa <- heatmap_top_metazoa %>%  select(perc) %>% as.matrix()
row.names(heatmap_metazoa) <- heatmap_top_metazoa$sample
mytree <- hclust(dist(heatmap_metazoa))
plot(mytree)
split.tree <- cutree(mytree, k=17)
rect.hclust(mytree,k=17,border = "red")
pheatmap(heatmap_metazoa, show_rownames = FALSE)






##################Arranging plots with letteres #########################
### Plot with damage, mean length and number of reads ###
ggdraw() +
  draw_plot(length_of_all_with_trend, x = 0, y = .4, width = 1, height = .6) +
  draw_plot_label(label = "A", size = 15,
                  x = c(0), y = c(1))
ggdraw() +
  draw_plot(damage_plot, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(pca_control_sample, x = .5, y = 0.5, width = .5, height = .5) +
  draw_plot(classified_reads, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(percentage_kingdoms, x = 0.5, y = 0, width = .5, height = .5) +
  draw_plot_label(label = c("B", "C", "D", "E"), size = 15,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))

ggdraw() +
  draw_plot(top80, x = 0, y = 0, width = 1, height = 1) +
  draw_plot_label(label = "F", size = 15,
                  x = c(0), y = c(1))

