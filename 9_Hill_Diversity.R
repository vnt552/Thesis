## library
library(devtools)
library(ggplot2)
library(tidyverse)

# Filter for Key Taxa, Here Insecta
# Alpha diversity looks at rare species, don't filter data by read count or normalize first. 
alpha_matrix <- metaDMG_results %>% filter(grepl("\\bInsecta\\b", tax_path), grepl("\\bGA\\b", sample), tax_rank=="species") %>%  select(sample, tax_name, N_reads)


# make read count matrix
alpha_wide_matrix <- pivot_wider(alpha_matrix, names_from = sample , values_from = N_reads)
alpha_wide_matrix[is.na(alpha_wide_matrix)] <- 0
alpha_wide_matrix <- as.data.frame(alpha_wide_matrix)
alpha_hill_matrix <- alpha_wide_matrix[,-1]
alpha_hill_species_list <- alpha_wide_matrix[,1]
t_alpha <- t(alpha_hill_matrix)

# For Beta diversity, use normalized insect data
norm_insecta


# Diversity using hillR package
library(hillR)
set.seed(123)

# alpha diversity for different q values (can also make q values between 0,1 and 2 for diversity per order.q plot)
hill_q0 <-  as.data.frame(hill_taxa(t_alpha, q=0))
hill_q1 <- as.data.frame(hill_taxa(t_alpha, q=1))
hill_q2 <- as.data.frame(hill_taxa(t_alpha, q=2))

# Summarized alha and beta diversity
hill_parti_q0 <- hill_taxa_parti(t_alpha, q=0)
hill_parti_q1 <- hill_taxa_parti(t_alpha, q=1)
hill_parti_q2 <- hill_taxa_parti(t_alpha, q=2)

# alpha and beta diversity grouped
hill_parti_pairwise_q0 <- hill_taxa_parti_pairwise(norm_insecta, q=0, output = "matrix")
hill_parti_pairwise_q1 <- hill_taxa_parti_pairwise(norm_insecta, q=1, output = "matrix")
hill_parti_pairwise_q2 <- hill_taxa_parti_pairwise(norm_insecta, q=2, output = "matrix")


## Diversity per order.q plot

all_hill <- cbind(hill_q0, hill_q025, hill_q05, hill_q075, hill_q1,hill_q1.25, hill_q1.5, hill_q1.75, hill_q2, hill_q2.25, hill_q2.5, hill_q2.75)
all_hill <- cbind(rownames(all_hill), data.frame(all_hill, row.names=NULL))
colnames(all_hill) <- c("Sample", "0", "0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2", "2.25", "2.5", "2.75")
long_hill <- all_hill %>% pivot_longer(cols = !Sample,values_to = "Diversity", names_to = "Order.q")

### plot q-values
labels_q0 <- long_hill %>%  filter(Order.q==0)
ggplot(long_hill, aes(x=Order.q, y=Diversity, group = Sample, color=Sample)) +
  geom_line() +
  geom_label_repel(data=labels_q0, aes(x= Order.q, y=Diversity, label=Sample), size=2, hjust=0,nudge_x= -0.5, direction="y") +
  theme_bw() +
  xlab("Order of q") +
  ylab("Diversity") +
  ggtitle("Diversity by Order of q", subtitle = "Abundance")

################# Alpha Diversity and Pollution ###################
# Linaer regression for Shannon div and pb
shannon_pb <- cbind(pollution_matrix, hill_q1[1:15,], hill_q0[1:15,], hill_q2[1:15,])
colnames(shannon_pb) <- c('Sample','Depth','Pb', "Cd", "s1", "Pb_mean", "Shannon_div", "Richness", "Simpson_div")
shannon_pb <- shannon_pb %>% select(Sample, Depth, Pb_mean, Shannon_div, Richness, Simpson_div)
shannon_pb$Sample <- factor(shannon_pb$Sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191"))

reg1 <- lm(Pb_mean ~ Shannon_div + Depth, data = shannon_pb)
summary(reg1)

# plot Shannon vs. pb
Shan_pb <- ggplot(shannon_pb) +
  geom_point(aes(x=Pb_mean, y=Shannon_div, color=Sample)) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("Pb Mean Conc.") +
  ylab("Shannon Diversity (q=1)") +
  ggtitle("Shannon Diversity", subtitle = "Insecta")
  
Simp_pb <- ggplot(shannon_pb) +
  geom_point(aes(x=Pb_mean, y=Simpson_div, color=Sample)) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("Pb Mean Conc.") +
  ylab("Simpson Diversity (q=2)") +
  ggtitle("Simpson Diversity", subtitle = "Insecta")

Rich_pb <- ggplot(shannon_pb) +
  geom_point(aes(x=Pb_mean, y=Richness, color=Sample)) +
  theme_bw() +
  xlab("Pb Mean Conc.") +
  ylab("Richness (q=0)") +
  ggtitle("Richness", subtitle = "Insecta")

ggdraw() +
  draw_plot(Shan_pb, x = 0, y = 0, width = 0.3, height = 0.7) +
  draw_plot(Simp_pb, x = .3, y = 0, width = .3, height = .7) +
  draw_plot(Rich_pb, x = 0.6, y = 0, width = .4, height = .7) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.3, .6), y = c(0.7, 0.7, 0.7))




####################### Beta Diversity #####################
########### beta NMDS plot, q0
#jaccard
beta_matrix_jac <- hill_parti_pairwise_q0$region_similarity
beta_matrix_jac <- Matrix::forceSymmetric(beta_matrix_jac)
beta_matrix_jac <- as.matrix(beta_matrix_jac)
beta_matrix_jac[is.na(beta_matrix_jac)] = 0

#bray_curtis
beta_matrix_bc <- hill_parti_pairwise_q0$local_similarity
beta_matrix_bc <- Matrix::forceSymmetric(beta_matrix_bc)
beta_matrix_bc <- as.matrix(beta_matrix_bc)
beta_matrix_bc[is.na(beta_matrix_bc)] = 0


bray_curtis_nmds <- metaMDS(beta_matrix_bc, k=3, distance = "Bray-Curtis")
bray_curtis_nmds$stress

jaccard_nmds <- metaMDS(beta_matrix_jac, k=3, distance = "jaccard")
jaccard_nmds$stress

data.scores.jac = as.data.frame(scores(jaccard_nmds, "sites"))
data.scores.jac$Sample = beta_matrix_jac$Sample
data.scores.jac$site <- rownames(data.scores.jac)

data.scores.bc = as.data.frame(scores(bray_curtis_nmds, "sites"))
data.scores.bc$Sample = beta_matrix_bc$Sample
data.scores.bc$site <- rownames(data.scores.bc)


## make plot
nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot() + 
  geom_point(data=data.scores.jac,aes(x=NMDS1,y=NMDS2,colour=site),size=3) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=15), # remove x-axis labels
        axis.title.y = element_text(size=15), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()) +
  scale_color_manual(values = mycolors) +
  ggtitle("NMDS, Insecta, Jaccard Dissimilarity Matrix, Stress = 0.19", subtitle = "hillR, Normalized")


################### PCoA plots
pcoa.bray <- cmdscale(beta_matrix_bc, k = 3, eig = T)

# extract axis positions for each site from cmdscale object and create a dataframe for plotting
pcoa.bray.plotting <- as.data.frame(pcoa.bray$points)
colnames(pcoa.bray.plotting) <- c("axis_1", "axis_2","axis_3")
pcoa.bray.plotting$site <- rownames(pcoa.bray.plotting)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pcoa.bray$eig[1]/(sum(pcoa.bray$eig))
pcoa.bray$eig[2]/(sum(pcoa.bray$eig))
pcoa.bray$eig[3]/(sum(pcoa.bray$eig))

# create a PCoA plot
pcoa.bray.plot <- ggplot(pcoa.bray.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_color_manual(values = mycolors) +
  theme_bw() + 
  xlab("PCoA 1 (18.8%)") +
  ylab("PCoA 2 (18.6%)") +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1)

# repeat process with Jaccard dissimilarity matrix
pcoa.jac <- cmdscale(beta_matrix_jac, k = 3, eig = T)

pcoa.jac.plotting <- as.data.frame(pcoa.jac$points)
colnames(pcoa.jac.plotting) <- c("axis_1", "axis_2","axis_3")
pcoa.jac.plotting$site <- rownames(pcoa.jac.plotting)

pcoa.jac$eig[1]/(sum(pcoa.jac$eig))
pcoa.jac$eig[2]/(sum(pcoa.jac$eig))
pcoa.jac$eig[3]/(sum(pcoa.jac$eig))

pcoa.jac.plot <- ggplot(pcoa.jac.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_color_manual(values = mycolors) +
  theme_bw() + 
  xlab("PCoA 1 (24.5%)") +
  ylab("PCoA 2 (23.9%)") +
  annotate(geom = 'text', label = 'Jaccard', x = Inf, y = -Inf, hjust = 1.215, vjust = -1)

# extract plot legend
legend <- get_legend(pcoa.jac.plot)

# plot Bray-Curtis PCoA and Jaccard PCoA side by side
plot_grid(pcoa.bray.plot + theme(legend.position = 'none'), pcoa.jac.plot + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1,1,0.5))



