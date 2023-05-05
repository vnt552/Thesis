######## SGA #################
sga <- read_csv("sga-trimming/sga.csv")
##plot
P1 <- ggplot(sga, aes(x=sga$Sample, y = sga$reads, fill = sga$threshold)) +
  geom_col(position="dodge") +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of reads") +
  ggtitle("SGA thresholds")
P1

## Loking into the length of sga files
sga1 <- read_table("sga-trimming/readlength-sga1.txt")
sga2 <- read_table("sga-trimming/readlength-sga2.txt")
sga3 <- read_table("sga-trimming/sga3_readlength.txt")
sga4 <- read_table("sga-trimming/sga4_readlength.txt")

total <- merge(sga2,sga3,by=c("Sample","Length"))
total2 <- merge(total,sga4,by=c("Sample","Length"))
total3 <- merge(total2,sga1,by=c("Sample","Length"))
long_total <- pivot_longer(total2, cols = starts_with("SGA"), names_to = "sga", values_to = "reads")

## plot
P2 <- ggplot(long_total, aes(x=Length, y=reads, fill=sga)) +
  geom_col(position="dodge") +
  theme_bw() +
  xlab("Length (bp)") +
  ylab("Number of Reads") +
  labs(fill="Dust Threshold") +
  ggtitle("Read Length for Dust 2, 3 and 4")
P2

### Percentage plots over SGA dust thresholds
sga_percent <- read_table("sga-trimming/sga_percent.txt")

sga_percent$Percent <- c(rep(NA,17),sga_percent$Reads[18:dim(sga_percent)[1]] / rep(sga_percent[1:17,]$Reads, 4))
sga_percent <- sga_percent[-c(0:17),]
sga_percent$Sample <- factor(sga_percent$Sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))


P3 <- ggplot(sga_percent, aes(x = Sample, y = Percent, fill=Filter)) +
  geom_col(position="dodge") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  xlab("Sample") +
  ylab("Percent Reads Left") +
  labs(fill="Dust Threshold") +
  ggtitle("Percentage Reads Left after Dust Threshold")
P3

### look into what the length of the reads removed in dust=3 are
length_diff <- merge(sga3,sga4, by=c("Sample","Length"))
length_diff <- length_diff %>% mutate(difference=SGA4-SGA3)
P4 <- ggplot(length_diff, aes(x=Length, y= difference, fill=Sample)) +
  geom_col(position="dodge") +
  theme_bw() +
  xlab("Length") +
  ylab("Difference in amount of reads") +
  labs(fill="Dust Threshold") +
  ggtitle("Dust=3 vs. Dust=4")
P4

#### Plotting of mapped reads for the dust thresholds
sga_thresholds <- read_table("sga-trimming/sga-thresholds.txt")
long_sga <- pivot_longer(sga_thresholds,cols = starts_with("sga"), names_to = "Dust", values_to = "mapped_reads")

P5 <- ggplot(long_sga, aes(x=Database, y=mapped_reads, fill=Dust)) +
  geom_col(position="dodge") +
  theme_bw() +
  xlab("Database Name") +
  ylab("Amount of Mapped Reads") +
  ggtitle("Mapped Reads per Dust Threshold") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
P5

## plots uden dust=1
sga_wo12 <- read_table("sga-trimming/sga-db_wo12.txt")
long_sga_wo12 <- pivot_longer(sga_wo12,cols = starts_with("sga"), names_to = "Dust_thresholds", values_to = "mapped_reads")


P6 <- ggplot(long_sga_wo12, aes(x=Database, y=mapped_reads, fill=Dust_thresholds)) +
  geom_col(position="dodge") +
  scale_y_sqrt() +
  theme_bw() +
  xlab("Database") +
  ylab("Amount of mapped reads") +
  labs(fill="Dust Threshold") +
  ggtitle("Mapped Reads pr. Dust Threshold, GA-11") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
P6


ggdraw() +
  draw_plot(P3, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(P2, x = .5, y = 0.5, width = .5, height = .5) +
  draw_plot(P6, x = 0, y = 0, width = 1, height = .5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))



################# ngsLCA #########################
library("maditr")
library("data.table")

# Looking into sample 161 and 61, to see if we get more reads by changing amount of mismatch from 95% to 90%
# load the data for the different threshold, filter for samples "GA-161" and "GA-61" and for key taxa
# normalize the data

#### metaDMG threshold 90,91,92,93,94,95
# download the files
results_90 <- read.csv("metaDMG/results_09.csv.gz") %>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_91 <- read.csv("metaDMG/results_091.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_92 <- read.csv("metaDMG/results_092.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_93 <- read.csv("metaDMG/results_093.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_94 <- read.csv("metaDMG/results_094.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_95 <- read.csv("metaDMG/results.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))


#Normaliseret matrix
Norm_key_91 <- sweep(norm_matrix_t_91, MARGIN=2, norm_matrix_euka_t_91, `/`)
Norm_key_92 <- sweep(norm_matrix_t_92, MARGIN=2, norm_matrix_euka_t_92, `/`)
Norm_key_93 <- sweep(norm_matrix_t_93, MARGIN=2, norm_matrix_euka_t_93, `/`)
Norm_key_94 <- sweep(norm_matrix_t_94, MARGIN=2, norm_matrix_euka_t_94, `/`)
Norm_key_90 <- sweep(norm_matrix_t_90, MARGIN=2, norm_matrix_euka_t_90, `/`)
Norm_key_95 <- sweep(norm_matrix_t_95, MARGIN=2, norm_matrix_euka_t_95, `/`)

Norm_key_90_df <- as.data.frame(as.table(Norm_key_90))
colnames(Norm_key_90_df) <- c("Taxa", "Sample","Reads")
Norm_key_90_df$threshold <- c("90")


#order the data
Norm_key_90_df$Taxa <- factor(Norm_key_90_df$Taxa , levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Chiroptera", "Insecta"))

norm_key_90_95 <- rbind(Norm_key_90_df, Norm_key_91_df, Norm_key_92_df, Norm_key_93_df, Norm_key_94_df, Norm_key_95_df)
norm_key_90_95_subset <- norm_key_90_95 %>% filter(Taxa %in% c("Insecta", "Chiroptera"))


## Plot
ggplot() +
  geom_point(Norm_key_90_df, mapping =aes(x=Taxa, y=Reads, color = threshold, shape = Sample)) +
  geom_line(Norm_key_90_df, mapping=aes(x=Taxa, y=Reads, group = Sample, color=threshold)) +
  geom_point(Norm_key_91_df, mapping =aes(x=Taxa, y=Reads, color = threshold, shape = Sample)) +
  geom_line(Norm_key_91_df, mapping=aes(x=Taxa, y=Reads, group = Sample, color=threshold)) +
  geom_point(Norm_key_92_df, mapping =aes(x=Taxa, y=Reads, color = threshold, shape = Sample)) +
  geom_line(Norm_key_92_df, mapping=aes(x=Taxa, y=Reads, group = Sample, color=threshold)) +
  geom_point(Norm_key_93_df, mapping =aes(x=Taxa, y=Reads, color = threshold, shape = Sample)) +
  geom_line(Norm_key_93_df, mapping=aes(x=Taxa, y=Reads, group = Sample, color=threshold)) +
  geom_point(Norm_key_94_df, mapping =aes(x=Taxa, y=Reads, color = threshold, shape = Sample)) +
  geom_line(Norm_key_94_df, mapping=aes(x=Taxa, y=Reads, group = Sample, color=threshold)) +
  geom_point(Norm_key_95_df, mapping =aes(x=Taxa, y=Reads, color = threshold, shape = Sample)) +
  geom_line(Norm_key_95_df, mapping=aes(x=Taxa, y=Reads, group = Sample, color=threshold)) +
  scale_color_brewer(palette="Set2") +
  xlab("Taxonomic Rank or Key Taxa") +
  ylab("Classified Reads") +
  ggtitle("Number of Classified Reads Per Threshold", subtitle = "Normalized to Eukaryote") +
  theme_bw()



ggplot(norm_key_90_95_subset, aes(x = Sample, y = Reads, fill = Taxa, group=interaction(threshold, Taxa))) +
  geom_col(position = "dodge", color= "Black") +
  geom_text(aes(label = threshold), vjust = -0.2, size = 3, position = position_dodge(.9)) +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  ggtitle("Key Taxa, 90 to 95 percentage similarity")

 
