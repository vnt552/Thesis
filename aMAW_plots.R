aMAW_mdmg <- read_csv("aMaw/tp-mdmg.lca.weight-1.csv.gz")

aMAW_mdmg_50 <- aMAW_mdmg %>% filter(N_reads > 50)

aMAW_kingdom <- aMAW_mdmg_50 %>%  filter(tax_rank=="kingdom")

kingdom_percentage <- aMAW_kingdom  %>%
  group_by(sample, tax_name) %>%
  summarise(n = sum(N_reads)) %>%
  mutate(percentage = n / sum(n))

# Place the samples in the rigth order
kingdom_percentage$sample <- factor(kingdom_percentage$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))


# make the plot
aMAW_kingdom_plot <- ggplot(kingdom_percentage, aes(y=percentage, x = sample, fill=tax_name)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  xlab("Sample") +
  ylab("% Classified Reads") +
  labs(fill="Kingdom") +
  ggtitle("Percentage Classified Reads, aMAW", subtitle = "Kingdom") +
  coord_flip()

######### Number of reads per sample, after filtering
amaw_total_reads_root <- aMAW_mdmg_50 %>%  filter(tax_name=="root")
amaw_total_reads_root$sample <- factor(amaw_total_reads_root$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))


classified_reads_amaw <- ggplot(amaw_total_reads_root, aes(y = N_reads, x = sample, group=1)) +
  geom_line() +
  theme_bw() +
  xlab("Sample") +
  ylab("Number of Classified Reads") +
  ggtitle("Number of Classified Reads", subtitle = "aMAW") +
  coord_flip()
classified_reads


######## Top 10 organisms per sample
## only reads classified to order level
amaw_family <- aMAW_mdmg_50 %>%  filter(tax_rank=="genus")

## Find top 10 family in each sample
amaw_abundance <- amaw_family  %>%
  group_by(sample, tax_name) %>%
  summarise(N_reads = sum(N_reads)) %>%
  mutate(RelativeAbundance = N_reads / sum(N_reads))

amaw_top_abundance <- amaw_abundance[order(amaw_abundance$N_reads, decreasing = TRUE), ]
amaw_top_abundance <- Reduce(rbind,                                 # Top N highest values by group
                                by(amaw_top_abundance,
                                   amaw_top_abundance["sample"],
                                   head,
                                   n = 10))
amaw_top_abundance<- amaw_top_abundance %>% select(tax_name, sample, RelativeAbundance)
amaw_top_abundance_wide <-amaw_top_abundance %>% pivot_wider(names_from = sample, values_from = RelativeAbundance)
amaw_top_abundance_header <- as.matrix(amaw_top_abundance_wide[,-1])
rownames(amaw_top_abundance_header) <-as.matrix(amaw_top_abundance_wide[,1])
amaw_top_abundance_header[is.na(amaw_top_abundance_header)] <- 0

pheatmap(amaw_top_abundance_header, color = hcl.colors(50, "BluYl"))

### this heatmap only shows the highest values, and many looks completly empty. Therefore scaling for all data
amaw_matrix_scaled <- t(scale(t(amaw_top_abundance_header)))
amaw_scaled_heatmap_abundance <- pheatmap(amaw_matrix_scaled, color = hcl.colors(50, "BluYl"))

################# PLOTS
aMAW_kingdom_plot
classified_reads_amaw
amaw_scaled_heatmap_abundance


ggdraw() +
  draw_plot(aMAW_kingdom_plot, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(classified_reads_amaw, x = .5, y = 0.5, width = .5, height = .5) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))
