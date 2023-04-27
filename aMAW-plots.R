library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggforce)


##### Read damage file
tp_mdmg <- read_csv("aMAW/tp-mdmg.lca.weight-1.csv.gz")
tp_mdmg$sample <- factor(tp_mdmg$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))

### Fordeling af number of reads
ggplot(tp_mdmg) +
  geom_point(aes(x = sample, y =N_reads)) +
  scale_y_log10() +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of reads") +
  ggtitle("Number of reads per sample")

order_level <- tp_mdmg %>% filter(tax_rank=="genus")

ggplot(order_level) +
  geom_point(aes(x = sample, y =N_reads)) +
  scale_y_log10() +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of reads") +
  ggtitle("Number of reads mapped to Genus level")

## Fordeling af damage
ggplot(tp_mdmg) +
  geom_point(aes(x = sample, y =D_max)) +
  theme_bw() +
  xlab("Sample") +
  ylab("Damage") +
  ggtitle("Damage across all reads")

#### Damage vs. number of reads
ggplot(order_level) +
  geom_point(aes(x = mean_L, y =D_max, colour = sample)) +
  #scale_shape_manual(values=c(16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17))
  scale_color_manual(values=c("brown3", "brown3", "brown3", "brown3", "brown3", "brown3", "brown3", "brown3", "brown3", "brown3", "brown3", "brown3", "brown3", "brown3", "brown3" ,"black", "black")) +
  theme_bw() +
  scale_x_log10() +
  xlab("mean length of reads") +
  ylab("Damage Max") +
  ggtitle("Damage based on read length")

## Which taxonomic rank are the reads mapped to?
ggplot(above_10reads) +
  geom_bar(aes(x = tax_rank)) +
  theme_bw() +
  xlab("Taxonomic Rank") +
  ylab("Count") +
  ggtitle("Taxonomic Rank Across all samples")



### data with only number of reads > 200
above_200reads <- tp_mdmg %>% filter(N_reads>=10)

ggplot(above_200reads, aes(y = D_max, x = sample)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Damage on samples with more than 200 reads")

### filter data aligned to genus level
ggplot(order_level, aes(y = D_max, x = sample)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Damage on samples, genus level, reads >= 200")

### Coverage for reads >=200, genus level
ggplot(order_level, aes(x = sample, y = c)) +
  geom_point() +
  theme_bw() +
  xlab("Sample") +
  ylab("Coverage") +
  ggtitle("Coverage for reads > 200")



### Look into the stats (initial and extension files together)
tp_stats <- read_tsv("aMAW-stats.txt")
tp_stats$file <- factor(tp_stats$file , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))

ggplot(tp_stats) +
  geom_col(mapping = aes(x = file, y = avg_len, fill= Initial_extension)) +
  theme_bw() +
  xlab("Sample") +
  ylab("Average Length") +
  ggtitle("aMAW read length before and after extension")

change_avg_len <- as.data.frame(tp_stats[1:6]) 
long_change <- pivot_wider(data = change_avg_len, names_from = file)

#### Read tp-mapping.summery.tsv.gz file
tp_map_sum <- read_tsv("tp-mapping.summary.tsv.gz")
tp_map_sum$label <- factor(tp_map_sum$label , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))

ggplot(tp_map_sum) +
  geom_sina(mapping = aes(x = label, y=read_length_mean)) +
  theme_bw() +
  xlab("Sample") +
  ylab("Mean Read Length") +
  ggtitle("Mean Read Length across samples")


#### Read tp-mapping-filtered.sumery.tsv file
tp_map_filt_sum <- read_tsv("tp-mapping-filtered.summary.tsv.gz")
tp_map_filt_sum$label <- factor(tp_map_filt_sum$label , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))

#options(scipen=10000)
ggplot(tp_map_filt_sum) +
  geom_sina(mapping = aes(x = label, y=mean_covered_bases)) +
  theme_bw() +
  scale_y_log10() +
  xlab("Sample") +
  ylab("Mean Covered Bases") +
  ggtitle("Mean Covered Bases after mapping and filtering")


ggplot(tp_map_filt_sum, mapping = aes(x = label, y=coverage_mean)) +
  geom_violin() +
  theme_bw() +
  scale_y_log10() +
  xlab("Sample") +
  ylab("Coverage Mean") +
  ggtitle("Mean Coverage after mapping and filtering")



