#### Basic ####
#Number of OTU allocated to phylum
num_otu <- metaDMG_results %>%  filter(tax_rank=="species")
mean(num_otu$N_reads)

num_otu_list <- as.list(num_otu$tax_id)
length(unique(num_otu_list)) ### number of unique species the reads have been classified as. 
sum(num_otu$N_reads) # Number of reads classified to species level
mean(num_otu$N_reads)
sd(num_otu$N_reads)

### Number OTU after filtering 50
num_otu_50 <- num_otu %>%  filter(N_reads>=50)

num_otu_list_50 <- as.list(num_otu_50$tax_id)
length(unique(num_otu_list_50))
sum(num_otu_50$N_reads)
mean(num_otu_50$N_reads)
sd(num_otu_50$N_reads)


# Key taxa Chiroptera
chiroptera_sp_count <- num_otu_50 %>%  filter(grepl("\\bChiroptera\\b", tax_path))
chiroptera_sp_list<- as.list(chiroptera_sp_count$tax_id)
length(unique(chiroptera_sp_list)) # antal species
sum(chiroptera_sp_count$N_reads) # antal reads


# Key taxa Insecta
insecta_sp_count <- num_otu_50 %>%  filter(grepl("\\bInsecta\\b", tax_path))
insecta_sp_list<- as.list(insecta_sp_count$tax_id)
length(unique(insecta_sp_list)) # antal species
sum(insecta_sp_count$N_reads) # antal reads


####################### DIET ###############################
metaDMG_results$sample[metaDMG_results$sample == "GA-Cex"] <- "Extraction_control"
metaDMG_results$sample[metaDMG_results$sample == "GA-Cli"] <- "Library_control"


Lepidoptera <- metaDMG_results %>% filter(N_reads >= 50, grepl("\\bLepidoptera\\b", tax_path), tax_rank=="family", grepl ("GA", sample))

###### DATA NEEDS TO BE NORMALIZED BEFORE PLOTTING

### NUmber of reads in insecta by order
#### insect order
Lepidoptera_dataframe <- t(Norm_Lepidoptera)
Lepidoptera_dataframe <- as.data.frame(Lepidoptera_dataframe)
Lepidoptera_dataframe <- tibble::rownames_to_column(Lepidoptera_dataframe, "sample")

Lepidoptera_long <- Lepidoptera_dataframe %>% pivot_longer(!sample, names_to = "tax_name", values_to = "N_reads")
Lepidoptera_long$sample <- factor(Lepidoptera_long$sample , levels=c("GA-Cli", "GA-Cex", "GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))
Lepidoptera_long$N_reads<-1000*(Lepidoptera_long$N_reads)

Lepidoptera_long$sample <- factor(Lepidoptera_long$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))
 

### plot Lepidopteran Families
lepi <- Lepidoptera_long
lepi$n_reads_cut <- cut(lepi$N_reads, breaks = c(0, 100,200, 300, 400, 500, 600, 700, 800),
                          include.lowest = TRUE, right = TRUE)
levels_count <- length(levels(lepi$n_reads_cut))
color_scales_fn <- colorRampPalette(c("#142D47", "#54AEF4"))
manual_color <- color_scales_fn(levels_count)
names(manual_color) <- levels(lepi$n_reads_cut)
manual_size <- seq(1, by = 1, length.out = levels_count)
names(manual_size) <- levels(lepi$n_reads_cut)


ggplot(lepi, aes(y= reorder(tax_name, N_reads), x = sample, color = n_reads_cut, size = n_reads_cut)) +
  geom_point() +
  scale_size_manual(values = manual_size, name="Normalized and thousandfold") +
  scale_color_manual(values = manual_color,name="Normalized and thousandfold") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust= -0.3)) +
  scale_x_discrete(position = "top") +
  ylab("Family") +
  xlab(" ") +
  ggtitle("Lepidopteran Families")

