#### Normalizing the data by CPM (counts per million)

metaDMG_results <- read_csv("metaDMG/results.csv.gz")

metaDMG_results$sample[metaDMG_results$sample == "GA-Cex"] <- "Extraction_control"
metaDMG_results$sample[metaDMG_results$sample == "GA-Cli"] <- "Library_control"

#Subset dataframe to show only genus level, filtering reads (choose best for your data)
metaDMG_50_genus <- metaDMG_results %>% filter(N_reads > 50,
                              grepl("\\bgenus\\b", tax_rank), grepl ("GA", sample))


#### Make a matrix, where each column is a sample, and each row is an organism.
sample1 <- metaDMG_50_genus$sample
tax_id1 <-metaDMG_50_genus$tax_name
reads1 <- metaDMG_50_genus$N_reads

raw1 <- data.frame(tax_id1, sample1, reads1)

data_wide <- dcast(raw1, sample1~ tax_id1 , value.var = "reads1", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide) #sum different reads per taxon

k <- ncol(data_wide)
dd1 = data_wide[,2:k]
rownames(dd1) <- data_wide$sample1

######## renaming and converting to matrix 
genus1  <-  dd1
genus <- as.matrix(genus1)
genus_matrix <- t(genus)


tss_genus <- tss(genus_matrix)

##########################
#normalization on sequencing depth
##########################
#1 Standardize abundances to the sequencing depth, CPM
physeq_median = genus
total <- 1e6

for (i in 1:dim(physeq_median)[2]) {
  x = physeq_median[,i]
  physeq_median[,i] = round(total * (x / sum(x)))
}


cpm_genus <- apply(subset(genus), 2, 
             function(x) x/sum(as.numeric(x)) * 10^6)
colSums(cpm_genus)

#Standardize abundances to the median sequencing depth
physeq_median = genus
total <- median(colSums(genus))

for (i in 1:dim(physeq_median)[2]) {
  x = physeq_median[,i]
  physeq_median[,i] = round(total * (x / sum(x)))
}








#############################
metaDMG_ga11 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-11", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga19 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-19", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga31 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-31", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga41 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-41", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga61 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-61", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga71 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-71", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga81 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-81", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga91 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-91", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga111 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-111", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga121 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-121", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga131 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-131", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga141 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-141", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga161 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-161", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga181 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-181", grepl("\\bgenus\\b", tax_rank))
metaDMG_ga191 <- metaDMG_results %>% filter(N_reads > 50, sample == "GA-191", grepl("\\bgenus\\b", tax_rank))




#### Make a matrix, where each column is a sample, and each row is an organism.
sample_ga11 <- metaDMG_ga11$sample
tax_id_ga11 <-metaDMG_ga11$tax_name
reads_ga11 <- metaDMG_ga11$N_reads
rawga11 <- data.frame(tax_id_ga11, sample_ga11, reads_ga11)

data_wide_ga11 <- dcast(rawga11, sample_ga11~ tax_id_ga11, value.var = "reads_ga11", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide_ga11) #sum different reads per taxon

k <- ncol(data_wide_ga11)
dd1 = data_wide_ga11[,2:k]
rownames(dd1) <- data_wide_ga11$sample_ga11

norm_matrix <- as.matrix(dd1)
norm_matrix <- t(norm_matrix)


######## Scaling factors
scaling_11 <- 96631/1000000
scaling_19 <- 79939/1000000
scaling_31 <- 102261/1000000
scaling_41 <- 28090/1000000
scaling_61 <- 13788/1000000
scaling_71 <- 67814/1000000
scaling_81 <- 112484/1000000
scaling_91 <- 49824/1000000
scaling_111 <- 34363/1000000
scaling_121 <- 57796/1000000
scaling_131 <- 26041/1000000
scaling_141 <- 32411/1000000
scaling_161 <- 50342/1000000
scaling_181 <- 22744/1000000
scaling_191 <- 22552/1000000


Rpm_norm_ga11 <- apply(subset(norm_matrix), 2, 
                  function(x) x/scaling_11)

Rpm_norm_ga11 <- apply(subset(norm_matrix), 2, 
                       function(x) x/(sum(x)/1000000))

Rpm_norm_ga19 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_19))

Rpm_norm_ga31 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_31))

Rpm_norm_ga41 <- apply(subset(norm_matrix), 2, 
                       function(x) x/scaling_41)

Rpm_norm_ga61 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_61))

Rpm_norm_ga71 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_71))

Rpm_norm_ga81 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_81))

Rpm_norm_ga91 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_91))

Rpm_norm_ga111 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_111))

Rpm_norm_ga121 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_121))

Rpm_norm_ga131 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_131))

Rpm_norm_ga141 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_141))

Rpm_norm_ga161 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_161))

Rpm_norm_ga181 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_181))

Rpm_norm_ga191 <- apply(subset(norm_matrix), 2, 
                       function(x) (x/scaling_191))





df_11 <- cpm_norm_ga11 %>% as.data.frame
df_11$rn <- rownames(df_11)
df_19 <- cpm_norm_ga19 %>%  as.data.frame
df_19$rn <- rownames(df_19)
df_31 <- cpm_norm_ga31 %>%  as.data.frame
df_31$rn <- rownames(df_31)
df_41 <- cpm_norm_ga41 %>%  as.data.frame
df_41$rn <- rownames(df_41)
df_61 <- cpm_norm_ga61 %>%  as.data.frame
df_61$rn <- rownames(df_61)
df_71 <- cpm_norm_ga71 %>%  as.data.frame
df_71$rn <- rownames(df_71)
df_81 <- cpm_norm_ga81 %>%  as.data.frame
df_81$rn <- rownames(df_81)
df_91 <- cpm_norm_ga91 %>%  as.data.frame
df_91$rn <- rownames(df_91)
df_111 <- cpm_norm_ga111 %>%  as.data.frame
df_111$rn <- rownames(df_111)
df_121 <- cpm_norm_ga121 %>%  as.data.frame
df_121$rn <- rownames(df_121)
df_131 <- cpm_norm_ga131 %>%  as.data.frame
df_131$rn <- rownames(df_131)
df_141 <- cpm_norm_ga141 %>%  as.data.frame
df_141$rn <- rownames(df_141)
df_161 <- cpm_norm_ga161 %>%  as.data.frame
df_161$rn <- rownames(df_161)
df_181 <- cpm_norm_ga181 %>%  as.data.frame
df_181$rn <- rownames(df_181)
df_191 <- cpm_norm_ga191 %>%  as.data.frame
df_191$rn <- rownames(df_191)


data_frame_merge <- merge(df_11, df_19,df_31, by = rn, all = TRUE)
df <- join_all(list(df_11,df_19,df_31,df_41, df_61, df_71, df_81, df_91, df_111, df_121, df_131, df_141, df_161, df_181, df_191), by = 'rn', type = 'full')
df_genus_merged <- df[,-2]
rownames(df_genus_merged) <- df[,2]

df_genus_merged <- t(df_genus_merged)
df_genus_merged <- cbind(sample=rownames(df_genus_merged),df_genus_merged)
df_genus_merged <- as.data.frame(df_genus_merged)
df_long <- df_genus_merged %>%  pivot_longer(!sample, names_to = "tax_name", values_to = "N_reads")

