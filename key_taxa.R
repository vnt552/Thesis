### Basics

metaDMG_results <- read_csv("metaDMG/results.csv.gz")

#### Basic ####
#Number of OTU allocated to phylum
num_otu <- metaDMG_results %>%  filter(tax_rank=="species")
mean(num_otu$N_reads)

num_otu_list <- as.list(num_otu$tax_id)
length(unique(num_otu_list)) ### number of unique species the reads have been classified as. 
sum(num_otu$N_reads) # Number of reads classified to species level


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


# Key taxa Chiroptera
insecta_sp_count <- num_otu_50 %>%  filter(grepl("\\bInsecta\\b", tax_path))
insecta_sp_list<- as.list(insecta_sp_count$tax_id)
length(unique(insecta_sp_list)) # antal species
sum(insecta_sp_count$N_reads) # antal reads
