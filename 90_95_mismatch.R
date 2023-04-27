library("maditr")
# Looking into sample 161 and 61, to see if we get more reads by changing amount of mismatch from 95% to 90%
results_90 <- read.csv("metaDMG/results_09.csv.gz")


# Først skal kun key species vælge, aka. chiroptera and insects
# normalize the data
chiroptera_90 <- results_90 %>% filter(N_reads >= 50, grepl("\\bChiroptera\\b", tax_name))
insecta_90 <- results_90 %>% filter(N_reads >= 50, grepl("\\bInsecta\\b", tax_name))
eukaryota_90 <- results_90 %>% filter(N_reads >= 50, tax_name == "Eukaryota")


# make the matrix for eukaryotes
################

euka_samples <- eukaryota_90$sample
euka_taxname <-eukaryota_90$tax_name
euka_reads <- eukaryota_90$N_reads
raw_euka <- data.frame(euka_taxname, euka_samples, euka_reads)

data_wide_euka <- dcast(raw_euka, euka_samples~ euka_taxname, value.var = "euka_reads", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide_euka) #sum different reads per taxon

k_euka <- ncol(data_wide_euka)
dd1_euka = data_wide_euka[,2:k_euka]


norm_matrix_euka <- as.matrix(dd1_euka)
rownames(norm_matrix_euka) <- data_wide_euka$euka_samples
norm_matrix_euka_t <- t(norm_matrix_euka)

#####
# make matrix for key taxa
################

key_taxa_90 <- rbind(chiroptera_90, insecta_90)

# make the matrix
key_samples <- key_taxa_90$sample
key_taxname <-key_taxa_90$tax_name
key_reads <- key_taxa_90$N_reads
raw_key <- data.frame(key_taxname, key_samples, key_reads)

data_wide_key <- dcast(raw_key, key_samples~ key_taxname, value.var = "key_reads", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide_key) #sum different reads per taxon

k <- ncol(data_wide_key)
dd1 = data_wide_key[,2:k]
row.names(dd1) <- data_wide_key$key_samples
norm_matrix <- as.matrix(dd1)
row.names(norm_matrix) <- data_wide_key$key_samples
norm_matrix_t <- t(norm_matrix)

#Normaliseret matrix
Norm_key <- sweep(norm_matrix_t, MARGIN=2, norm_matrix_euka_t, `/`)

##############

# Hent key data fra 95-100% og normaliser det
########
metaDMG_results <- read.csv("metaDMG/results.csv.gz")
chirop_61 <- metaDMG_results %>% filter(N_reads >= 50, sample == "GA-61", grepl("\\bChiroptera\\b", tax_name)) 
chirop_161 <- metaDMG_results %>% filter(N_reads >= 50, sample == "GA-161", grepl("\\bChiroptera\\b", tax_name)) 
insect_61 <- metaDMG_results %>% filter(N_reads >= 50, sample == "GA-61", grepl("\\bInsecta\\b", tax_name)) 
insect_161 <- metaDMG_results %>% filter(N_reads >= 50, sample == "GA-161", grepl("\\bInsecta\\b", tax_name)) 

key_taxa_95 <- rbind(chirop_161, chirop_61, insect_161, insect_61)


eukaryote_95_61 <- metaDMG_results %>% filter(N_reads >= 50, tax_name == "Eukaryota", sample == "GA-61")
eukaryote_95_161 <- metaDMG_results %>% filter(N_reads >= 50, tax_name == "Eukaryota", sample == "GA-161")
eukaryote_95 <- rbind(eukaryote_95_161, eukaryote_95_61)


eukaryote_95_samples <- eukaryote_95$sample
eukaryote_95_taxname <-eukaryote_95$tax_name
eukaryote_95_reads <- eukaryote_95$N_reads
raw_eukaryote_95 <- data.frame(eukaryote_95_taxname, eukaryote_95_samples, eukaryote_95_reads)

data_wide_eukaryote_95 <- dcast(raw_eukaryote_95, eukaryote_95_samples~ eukaryote_95_taxname, value.var = "eukaryote_95_reads", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide_eukaryote_95) #sum different reads per taxon

k_eukaryote_95 <- ncol(data_wide_eukaryote_95)
dd1_eukaryote_95 = data_wide_eukaryote_95[,2:k_eukaryote_95]

norm_matrix_eukaryote_95 <- as.matrix(dd1_eukaryote_95)
rownames(norm_matrix_eukaryote_95) <- data_wide_eukaryote_95$eukaryote_95_samples
norm_matrix_eukaryote_95_t <- t(norm_matrix_eukaryote_95)

#####
# make matrix for key taxa

key_95_samples <- key_taxa_95$sample
key_95_taxname <-key_taxa_95$tax_name
key_95_reads <- key_taxa_95$N_reads
raw_key_95 <- data.frame(key_95_taxname, key_95_samples, key_95_reads)

data_wide_key_95 <- dcast(raw_key_95, key_95_samples~ key_95_taxname, value.var = "key_95_reads", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide_key_95) #sum different reads per taxon

k_key_95 <- ncol(data_wide_key_95)
dd1_key_95 = data_wide_key[,2:k_key_95]

norm_matrix_key_95 <- as.matrix(dd1_key_95)
rownames(norm_matrix_key_95) <- data_wide_key_95$key_95_samples
norm_matrix_key_95_t <- t(norm_matrix_key_95)

#Normaliseret matrix
Norm_key_95 <- sweep(norm_matrix_key_95_t, MARGIN=2, norm_matrix_eukaryote_95_t, `/`)



#########

key_taxa_90_dataframe <- as.data.frame(as.table(Norm_key))
colnames(key_taxa_90_dataframe) <- c("Taxa", "Sample","Reads")
key_taxa_95_dataframe <- as.data.frame(as.table(Norm_key_95))
colnames(key_taxa_95_dataframe) <- c("Taxa", "Sample","Reads")

key_taxa_95_dataframe$mismatch <- c("95")
key_taxa_90_dataframe$mismatch <- c("90")

normalised_key_taxa <- rbind(key_taxa_90_dataframe, key_taxa_95_dataframe)
normalised_key_taxa$Sample <- factor(normalised_key_taxa$Sample , levels=c("GA-61", "GA-161"))



### Singletons in 90 and 95. THIS IS NOT NORMALISED DATA
singletons_90 <- results_90 %>% filter(N_reads==1, tax_rank == "species")
singleton_list_90 <- as.list(singletons_90$tax_id)
length(unique(singleton_list_90))

singletons_95_61 <- metaDMG_results %>% filter(N_reads==1, tax_rank == "species", grepl("\\bGA-61\\b", sample))
singletons_95_161 <- metaDMG_results %>% filter(N_reads==1, tax_rank == "species", grepl("\\bGA-161\\b", sample))
singletons_95 <- rbind(singletons_95_61, singletons_95_161)
singleton_list_95 <- as.list(singletons_95$tax_id)
length(unique(singleton_list_95))




## Plot
ggplot(normalised_key_taxa, aes(x = Sample, y = Reads, fill = Taxa, group=interaction(mismatch, Taxa))) +
  geom_col(position = "dodge", color= "Black") +
  geom_text(aes(label = mismatch), vjust = -0.2, size = 3, position = position_dodge(.9)) +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  ggtitle("Key Taxa, 90 and 95 percentage similarity")










### New taxa found in 90? NORMALISED DATA?
new_taxa_90 <- results_90 %>% filter(N_reads>=50, tax_rank=="species", grepl("\\bEukaryota\\b", tax_path))
new_taxa_95_61 <- metaDMG_results %>% filter(N_reads>=50, tax_rank=="species", sample=="GA-61", grepl("\\bEukaryota\\b", tax_path))
new_taxa_95_161 <- metaDMG_results %>% filter(N_reads>=50, tax_rank=="species", sample=="GA-161", grepl("\\bEukaryota\\b", tax_path))
new_taxa_95 <- rbind(new_taxa_95_161, new_taxa_95_61)

new_taxa_90_list <- as.list(new_taxa_90$tax_name)
new_taxa_95_list <- as.list(new_taxa_95$tax_name)


new_taxa <- setdiff(new_taxa_90_list, new_taxa_95_list)
setdiff(new_taxa_95$tax_name, new_taxa_90$tax_name)


## Antal insect og chiroptera reads i 95

#######
# ændringer mellem 90 and 95



#### How many reads are classified to key taxa in 90 and 95? 
chirop_90 <- results_90 %>% filter(tax_name=="Chiroptera")
insecta_90 <- results_90 %>% filter(tax_name=="Insecta")

chirop_95 <- metaDMG_results %>% filter(tax_name=="Chiroptera")
insecta_95<- metaDMG_results %>% filter(tax_name=="Insecta")

try1 <- metaDMG_results %>% filter(tax_name=="Chiroptera", sample=="GA-11")














#### metaDMG threshold 90,91,92,93,94,95
# download the files
results_90 <- read.csv("metaDMG/results_09.csv.gz") %>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_91 <- read.csv("metaDMG/results_091.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_92 <- read.csv("metaDMG/results_092.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_93 <- read.csv("metaDMG/results_093.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_94 <- read.csv("metaDMG/results_094.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))
results_95 <- read.csv("metaDMG/results.csv.gz")%>%  filter(grepl("\\bEukaryota\\b", tax_path))

res_95_61 <- results_95 %>%  filter(sample=="GA-61")
res_95_161 <- results_95 %>%  filter(sample=="GA-161")

results_95 <- rbind(res_95_161, res_95_61)

# eukaryote superkingdom for each file
eukaryota <- results_95 %>% filter(N_reads >= 50, tax_name == "Eukaryota")

# euka data for normalization
euka_samples <- eukaryota$sample
euka_taxname <-eukaryota$tax_name
euka_reads <- eukaryota$N_reads
raw_euka <- data.frame(euka_taxname, euka_samples, euka_reads)

data_wide_euka <- dcast(raw_euka, euka_samples~ euka_taxname, value.var = "euka_reads", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide_euka) #sum different reads per taxon

k_euka <- ncol(data_wide_euka)
dd1_euka = data_wide_euka[,2:k_euka]


norm_matrix_euka <- as.matrix(dd1_euka)
rownames(norm_matrix_euka) <- data_wide_euka$euka_samples
norm_matrix_euka_t <- t(norm_matrix_euka)

#####
# make matrix for key taxa, phyla, class, order, family, genus and species
################
chiroptera_norm <- results_95 %>% filter(N_reads >= 50, grepl("\\bChiroptera\\b", tax_name))
chiroptera_norm$key <- c("Chiroptera")

insecta_norm <- results_95 %>% filter(N_reads >= 50, grepl("\\bInsecta\\b", tax_name))
insecta_norm$key <- c("Insecta")

phyla_norm <- results_95 %>% filter(N_reads >= 50, grepl("\\bphylum\\b", tax_rank))
phyla_norm$key <- c("Phylum")

class_norm <- results_95 %>% filter(N_reads >= 50, grepl("\\bclass\\b", tax_rank))
class_norm$key <- c("Class")

order_norm <- results_95 %>% filter(N_reads >= 50, grepl("\\border\\b", tax_rank))
order_norm$key <- c("Order")

family_norm <- results_95 %>% filter(N_reads >= 50, grepl("\\bfamily\\b", tax_rank))
family_norm$key <- c("Family")

genus_norm <- results_95 %>% filter(N_reads >= 50, grepl("\\bgenus\\b", tax_rank))
genus_norm$key <- c("Genus")

species_norm <- results_95 %>% filter(N_reads >= 50, grepl("\\bspecies\\b", tax_rank))
species_norm$key <- c("Species")

key_data <- rbind(chiroptera_norm, insecta_norm, phyla_norm, class_norm, order_norm, family_norm, genus_norm, species_norm)

# make the matrix
key_samples <- key_data$sample
key_taxname <-key_data$tax_name
key_reads <- key_data$N_reads
key_rank <- key_data$key
raw_key <- data.frame(key_taxname, key_samples, key_reads, key_rank)

data_wide_key <- dcast(raw_key, key_samples~ key_rank, value.var = "key_reads", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide_key) #sum different reads per taxon

k <- ncol(data_wide_key)
dd1 = data_wide_key[,2:k]
row.names(dd1) <- data_wide_key$key_samples
norm_matrix <- as.matrix(dd1)
row.names(norm_matrix) <- data_wide_key$key_samples
norm_matrix_t <- t(norm_matrix)

#Normaliseret matrix
#Norm_key_91 <- sweep(norm_matrix_t, MARGIN=2, norm_matrix_euka_t, `/`)
#Norm_key_92 <- sweep(norm_matrix_t, MARGIN=2, norm_matrix_euka_t, `/`)
#Norm_key_93 <- sweep(norm_matrix_t, MARGIN=2, norm_matrix_euka_t, `/`)
#Norm_key_94 <- sweep(norm_matrix_t, MARGIN=2, norm_matrix_euka_t, `/`)
#Norm_key_90 <- sweep(norm_matrix_t, MARGIN=2, norm_matrix_euka_t, `/`)
#Norm_key_95 <- sweep(norm_matrix_t, MARGIN=2, norm_matrix_euka_t, `/`)

Norm_key_95_df <- as.data.frame(as.table(Norm_key_95))
colnames(Norm_key_95_df) <- c("Taxa", "Sample","Reads")
Norm_key_95_df$threshold <- c("95")


#order the data
Norm_key_90_df$Taxa <- factor(Norm_key_90_df$Taxa , levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Chiroptera", "Insecta"))
Norm_key_91_df$Taxa <- factor(Norm_key_91_df$Taxa , levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Chiroptera", "Insecta"))
Norm_key_92_df$Taxa <- factor(Norm_key_92_df$Taxa , levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Chiroptera", "Insecta"))
Norm_key_93_df$Taxa <- factor(Norm_key_93_df$Taxa , levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Chiroptera", "Insecta"))
Norm_key_94_df$Taxa <- factor(Norm_key_94_df$Taxa , levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Chiroptera", "Insecta"))
Norm_key_95_df$Taxa <- factor(Norm_key_95_df$Taxa , levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Chiroptera", "Insecta"))




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



 



