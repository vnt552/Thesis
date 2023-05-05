# Normalized of reads to their superkingdom

# first load the result file
metaDMG_results <- read.csv("metaDMG/results.csv.gz")

# If Controls are not used, filter them out by giving them a new name
metaDMG_results$sample[metaDMG_results$sample == "GA-Cex"] <- "Extraction_control"
metaDMG_results$sample[metaDMG_results$sample == "GA-Cli"] <- "Library_control"

# Eukaryote Superkingdom vector
eukaryote_superkingdom <- superkingdoms %>% filter(N_reads > 50, tax_name == "Eukaryota", grepl ("GA", sample))
euka_data <- eukaryote_superkingdom %>% select("N_reads")
rownames(euka_data) <- eukaryote_superkingdom$sample
matrix_euka_superkingdom <- as.matrix(euka_data)
matrix_euka_superkingdom <- t(matrix_euka_superkingdom)
vector_euka <- as_vector(matrix_euka_superkingdom)


##### CHIROPTERA GENUS
### Now make genus level files
chiroptera_genus <- metaDMG_results %>% filter(N_reads > 50, grepl("\\bChiroptera\\b", tax_path), tax_rank=="genus")


## Make read count matrix
# each column is a sample, and each row is an organism.
chiroptera_samples <- chiroptera_genus$sample
chiroptera_taxname <-chiroptera_genus$tax_name
chiroptera_reads <- chiroptera_genus$N_reads
raw_chiroptera <- data.frame(chiroptera_taxname, chiroptera_samples, chiroptera_reads)

data_wide_chiroptera <- dcast(raw_chiroptera, chiroptera_samples~ chiroptera_taxname, value.var = "chiroptera_reads", fun.aggregate = sum) #N_reads each Tax_name genus (=each sample)
head(data_wide_chiroptera) #sum different reads per taxon

k <- ncol(data_wide_chiroptera)
dd1 = data_wide_chiroptera[,2:k]
rownames(dd1) <- data_wide_chiroptera$chiroptera_samples

norm_matrix <- as.matrix(dd1)
norm_matrix_t <- t(norm_matrix)

# check to see if the superkingdom and the genus matrix contains the same samples
colnames(vector_euka)
colnames(norm_matrix_t)

# if they do not contain the same samples, drop them from the vector_euka matrix

matrix_euka_superkingdom_df <- as.data.frame(vector_euka
drop <- c("GA-181")
MyNewdataframe = matrix_euka_superkingdom_df[,!(names(matrix_euka_superkingdom_df) %in% drop)]
MyNewdataframe
MyNewMatrix <- as.matrix(MyNewdataframe)

MyNewMatrix <- MyNewMatrix[, order(as.integer(colnames(MyNewMatrix)))]

### Divide the Chiroptera Matrix by the Eukaryote Matrix.
Norm_chiroptera <- sweep(norm_matrix_t, MARGIN=2, MyNewMatrix, `/`)

