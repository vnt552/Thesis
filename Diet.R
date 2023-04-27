#### Look at the amount of insects and chiroptera is classified in the samples

metaDMG_results <- read.csv("metaDMG/results.csv.gz")
metaDMG_results$sample[metaDMG_results$sample == "GA-Cex"] <- "Extraction_control"
metaDMG_results$sample[metaDMG_results$sample == "GA-Cli"] <- "Library_control"


Lepidoptera <- metaDMG_results %>% filter(N_reads >= 50, grepl("\\bLepidoptera\\b", tax_path), tax_rank=="family", grepl ("GA", sample))


## Make read count matrix
# each column is a sample, and each row is an organism.
Lepidoptera_samples <- Lepidoptera$sample
Lepidoptera_taxname <-Lepidoptera$tax_name
Lepidoptera_reads <- Lepidoptera$N_reads
raw_Lepidoptera <- data.frame(Lepidoptera_taxname, Lepidoptera_samples, Lepidoptera_reads)

data_wide_Lepidoptera <- dcast(raw_Lepidoptera, Lepidoptera_samples~ Lepidoptera_taxname, value.var = "Lepidoptera_reads", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide_Lepidoptera) #sum different reads per taxon

k_Lepidoptera <- ncol(data_wide_Lepidoptera)
dd1_Lepidoptera = as.data.frame(data_wide_Lepidoptera[,2:k_Lepidoptera])
rownames(dd1_Lepidoptera) <- data_wide_Lepidoptera$Lepidoptera_samples

norm_matrix_Lepidoptera <- as.matrix(dd1_Lepidoptera)
norm_matrix_t_Lepidoptera <- t(norm_matrix_Lepidoptera)


#### Now the reads need to be divided by the superkingdom from each sample
colnames(vector_euka)
colnames(norm_matrix_t_insect_order)

###Normaliseret matrix

Norm_Lepidoptera <- sweep(norm_matrix_t_Lepidoptera, MARGIN=2, vector_euka, `/`)



### NUmber of reads in insecta by order
#### insect order
Lepidoptera_dataframe <- t(Norm_Lepidoptera)
Lepidoptera_dataframe <- as.data.frame(Lepidoptera_dataframe)
Lepidoptera_dataframe <- tibble::rownames_to_column(Lepidoptera_dataframe, "sample")

Lepidoptera_long <- Lepidoptera_dataframe %>% pivot_longer(!sample, names_to = "tax_name", values_to = "N_reads")
Lepidoptera_long$sample <- factor(Lepidoptera_long$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))
Lepidoptera_long$N_reads<-1000*(Lepidoptera_long$N_reads)

Lepidoptera_long$sample <- factor(Lepidoptera_long$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))




###### DATA NEEDS TO BE NORMALIZED BEFORE PLOTTING ###########

lepi <- Lepidoptera_long
lepi$n_reads_cut <- cut(lepi$N_reads, breaks = c(0, 1, 50, 100, 200, 300, 400, 500, 600, 700),
                          include.lowest = TRUE, right = TRUE)
levels_count <- length(levels(lepi$n_reads_cut))
color_scales_fn <- colorRampPalette(c("#142D47", "#54AEF4"))
manual_color <- color_scales_fn(levels_count)
names(manual_color) <- levels(lepi$n_reads_cut)
manual_size <- seq(1, by = 1, length.out = levels_count)
names(manual_size) <- levels(lepi$n_reads_cut)




ggplot(lepi, aes(y= reorder(tax_name, N_reads), x = sample, color = n_reads_cut, size = n_reads_cut)) +
  geom_point() +
  scale_size_manual(values = manual_size, name="pr. 1000 Reads") +
  scale_color_manual(values = manual_color,name="pr. 1000 Reads") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust= -0.3)) +
  scale_x_discrete(position = "top") +
  ylab("Family") +
  xlab(" ") +
  ggtitle("Lepidopteran Families")


######
# is the diet pest arthropods?
insecta <- metaDMG_results %>% filter(N_reads >= 50, grepl("\\bInsecta\\b", tax_path), tax_rank=="species")

pest_arthropods <- select(insecta, sample, tax_name, N_reads)
wide_pest <- pivot_wider(pest_arthropods, names_from = sample, values_from = N_reads)





  