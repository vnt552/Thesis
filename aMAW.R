### Relative Abundance

# load file
aMAW_mdmg <- read_csv("aMaw/tp-mdmg.lca.weight-1.csv.gz")

aMAW_kingdom <- filter(tax_rank=="kingdom")

kingdom_percentage <- aMAW_kingdom  %>%
  group_by(sample, tax_name) %>%
  summarise(n = sum(N_reads)) %>%
  mutate(percentage = n / sum(n))

# Place the samples in the rigth order
kingdom_percentage$sample <- factor(kingdom_percentage$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))

# make relative abundance  plot
aMAW_kingdom_plot <- ggplot(kingdom_percentage, aes(y=percentage, x = sample, fill=tax_name)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  xlab("Sample") +
  ylab("% Classified Reads") +
  labs(fill="Kingdom") +
  ggtitle("Relative Abundance", subtitle = "aMAW Kingdoms") +
  coord_flip()
  
  
  
### Heat map of Normalized abundance
bacti_dataframe <- t(Normalized_bacti) # Use normalized data matrix, tax_rank = species, read count >=200
bacti_dataframe <- as.data.frame(bacti_dataframe)

row_order <- c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-Cex")
bacti_dataframe <- tibble::rownames_to_column(bacti_dataframe, "sample")
bacti_dataframe <- bacti_dataframe %>% arrange(factor(sample, levels = row_order))

numeric_cols <- vapply(bacti_dataframe, is.numeric, logical(1))
bacti_dataframe[, numeric_cols] <- lapply(bacti_dataframe[, numeric_cols, drop = FALSE],
                                         function(x) x * 1000)

samp2 <- bacti_dataframe[,-1]
rownames(samp2) <- bacti_dataframe[,1]

nb.cols <- 15
my_colour2 <- colorRampPalette(brewer.pal(9, "GnBu"))(nb.cols)


pheatmap(
  mat               = samp2,
  color =  my_colour2,
  breaks = c(0, 20,40,60, 80, 100, 200, 300,400, 500, 600, 700,800, 900, 1000),
  fontsize          = 10,
  main              = "Bacterial reads, above N_reads>200, species leve",
  cluster_rows = F,
  cluster_cols = F
)

