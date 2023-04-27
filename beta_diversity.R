### Diversity
metaDMG_results <- read.csv("metaDMG/results.csv.gz")



# Alpha diversity - non hill
#####
# alpha diversity
alpha_matrix <- metaDMG_results %>% filter(tax_rank=="species") %>% select(sample, tax_name, N_reads)

# random reads for more statistical evenness
random <- alpha_matrix %>% 
  uncount(N_reads) %>% 
  mutate(name=sample(tax_name)) %>% 
  count(sample, name, name="N_reads")


alpha_diversity <- random %>% group_by(sample) %>% 
  summarize(sobs = specnumber(N_reads), 
            shannon_index = diversity(N_reads, index = "shannon"),
            simpson_index = diversity(N_reads, index = "simpson"),
            inv_simpson = 1/simpson_index,
            evenness = shannon_index/log(sobs),
            n=sum(N_reads))

# pivot longer
alpha_long <- alpha_diversity %>%
  pivot_longer(cols = c(sobs, shannon_index, simpson_index, inv_simpson, evenness),
               names_to = "metric")

# change the names
alpha_long$metric <- factor(alpha_long$metric, levels = c("sobs", "shannon_index", "simpson_index", "inv_simpson", "evenness"),
                            labels = c("Richness", "Shannon Index", "Simpson Index", "Inverse Simpson Index", "Pielou's Evenness" )) 
alpha_long$sample <- factor(alpha_long$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))


#plot
alpha_random_diversity <- ggplot(alpha_long, aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow = 4, scales = "free_y") +
  xlab("Number of Reads") +
  ylab("Diversity Value") +
  ggtitle("Alpha Diversity with Random Sample Size") +
  theme_minimal()
alpha_random_diversity

alpha_not_random_diversity <- ggplot(alpha_long, aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow = 4, scales = "free_y") +
  xlab("Number of Reads") +
  ylab("Diversity Value") +
  ggtitle("Alpha Diversity without Standardization") +
  theme_minimal()


ggplot(alpha_long, aes(x = sample, y = value)) +
  geom_point() +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#### correlation with depth and pollution
alpha_wide <- alpha_long %>%  pivot_wider(names_from = metric, values_from = value)
alpha_wide$sample <- factor(alpha_wide$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))
alpha_wide$depth <- c(11,111,121,131,141,161,181,19,191,31,41,61,71,81,91,0,0)

summary(lm(shannon_index ~ depth, alpha_wide))

### shannon diversity vs. depth of the samples
depth.reg <- ggplot(alpha_wide, aes(x = depth, y = shannon_index, colour = sample)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Depth (m)") + 
  ylab("Shannon's H") +
  ggtitle("Sample Depth and Shannon's H")+
  theme_bw()

depth.reg
##### ADD POLLUTION DATA!
pollution_data <- read_csv("pollution2000.csv")

subset_pollution <- pollution_data[c(89,160,275,375,524,623,719,819,967,1064,1164,1264,1403,1572,1622),]
subset_pollution$sample <- c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191")
subset_alpha_div <- alpha_diversity[c(1:15),]

alpha_pollution_matrix <- merge(subset_pollution, subset_alpha_div, by="sample")  
alpha_pollution_matrix$sample <- factor(alpha_pollution_matrix$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191"))

summary(lm(shannon_index ~ Pb, alpha_pollution_matrix))



ggplot(alpha_pollution_matrix, aes(x = Pb, y = sobs, colour = sample)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Pb") + 
  ylab("Shannon's H'") +
  theme_bw()




















# Beta diversity
######
library(adespatial)

## make the matrix
no_cont <- alpha_matrix[alpha_matrix$sample != "GA-Cex" & alpha_matrix$sample != "GA-Cli", ] 

beta_matrix1 <- no_cont %>%  pivot_wider(names_from = tax_name, values_from = N_reads)
beta_matrix1 <- as.data.frame(beta_matrix1)
beta_matrix1$sample <- factor(beta_matrix1$sample , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191"))
beta_matrix <- beta_matrix1[,-1]
beta_matrix <- as.data.frame(beta_matrix)
sample_names <- c("GA-11","GA-111","GA-121","GA-131","GA-141","GA-161","GA-181","GA-19","GA-191","GA-31","GA-41","GA-61","GA-71","GA-81","GA-91")
rownames(beta_matrix) <- sample_names
beta_matrix[is.na(beta_matrix)] = 0


#####
#### beta diversity calculated with beta.div.comp
beta_div_j <- beta.div.comp(beta_matrix, coef = "J", quant=T)
beta_div_j$part

beta_div_s <- beta.div.comp(beta_matrix, coef = "S", quant=T)
beta_div_s$part

# local contribution to the beta diversity, richness 
local_beta_relp <- LCBD.comp(beta_div_j$repl, sqrt.D = T)
local_beta_relp

local_beta_rich <- LCBD.comp(beta_div_j$rich, sqrt.D = T)
local_beta_rich


# species contribution on beta diversity 
species_beta <- beta.div(beta_matrix, method = "hellinger")
species_beta

#####
######### BETA DIVERSITY USING VEGAN PACKAGE
library(vegan)

# find the smallest number of species found, use that number for the avgdist
no_cont %>%  group_by(sample) %>%  summarise(total = sum(N_reads)) %>% arrange(total)

vegan.bray <- vegdist(beta_matrix, method = "bray")
vegan.bray
vegan.jac <- vegdist(beta_matrix, method = "jaccard")
vegan.jac

# rarefying the data and calculating distance matrix

avg.bray <- avgdist(beta_matrix, sample = 13857, dmethod = "bray")
avg.bray
avg.jac <- avgdist(beta_matrix, sample = 559, dmethod = "jaccard")
avg.jac

# contribution of each species
simper(beta_matrix, permutations = 999)

#plotting PCA of Beta Diversity
#####
### Principal Coordinates Analysis (PCoA)
# calculate principal coordinates analysis (Bray-Curtis)
pcoa.avg.bray <- cmdscale(avg.bray, k = 2, eig = T)

# extract axis positions for each site from cmdscale object and create a dataframe for plotting
pcoa.avg.bray.plotting <- as.data.frame(pcoa.avg.bray$points)
colnames(pcoa.avg.bray.plotting) <- c("axis_1", "axis_2")
pcoa.avg.bray.plotting$site <- rownames(pcoa.avg.bray.plotting)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pcoa.avg.bray$eig[1]/(sum(pcoa.avg.bray$eig))
pcoa.avg.bray$eig[2]/(sum(pcoa.avg.bray$eig))

# create a PCoA plot
pcoa.avg.bray.plot <- ggplot(pcoa.avg.bray.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.01, end = 0.9) +
  theme_bw() + 
  xlab("PCoA 1 (23.6%)") +
  ylab("PCoA 2 (17.1%)") +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1)

# repeat process with Jaccard dissimilarity matrix
pcoa.avg.jac <- cmdscale(avg.jac, k = 2, eig = T)

pcoa.avg.jac.plotting <- as.data.frame(pcoa.avg.jac$points)
colnames(pcoa.avg.jac.plotting) <- c("axis_1", "axis_2")
pcoa.avg.jac.plotting$site <- rownames(pcoa.avg.jac.plotting)

pcoa.avg.jac$eig[1]/(sum(pcoa.avg.jac$eig))
pcoa.avg.jac$eig[2]/(sum(pcoa.avg.jac$eig))
pcoa.avg.jac.plot <- ggplot(pcoa.avg.jac.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.01, end = 0.9) +
  theme_bw() + 
  xlab("PCoA 1 (23.2%)") +
  ylab("PCoA 2 (17.8%)") +
  annotate(geom = 'text', label = 'Jaccard', x = Inf, y = -Inf, hjust = 1.215, vjust = -1)

# extract plot legend
legend <- get_legend(pcoa.avg.jac.plot)

# plot Bray-Curtis PCoA and Jaccard PCoA side by side
plot_grid(pcoa.avg.bray.plot + theme(legend.position = 'none'), pcoa.avg.jac.plot + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1,1,0.5))



#####
#### Anova test for beta diversity
pollution_matrix # environmental data

#####

#vegan distance
adonis2(vegan.bray ~ Pb_mean, data = pollution_matrix, permutations = 999)
adonis2(vegan.jac ~ Pb_mean, data = pollution_matrix, permutations = 999)
adonis2(vegan.jac ~ depth, data = pollution_matrix, permutations = 999)
adonis2(vegan.jac ~ Pb_mean + depth, data =pollution_matrix, permutations = 999)
adonis2(vegan.jac ~ Pb_mean *depth, data = pollution_matrix, permutations = 999)

#####
#avgerage distance
adonis2(avg.bray ~ Pb_mean, data = pollution_matrix, permutations = 999)
adonis2(avg.bray ~ depth, data = pollution_matrix, permutations = 999)
adonis2(avg.jac ~ Pb_mean, data = pollution_matrix, permutations = 999)
adonis2(avg.jac ~ depth, data = pollution_matrix, permutations = 999)
adonis2(avg.jac ~ Pb_mean + depth, data = pollution_matrix, permutations = 999)
adonis2(avg.jac ~ Pb_mean *depth, data = pollution_matrix, permutations = 999)





bd <- betadisper(avg.jac, pollution_matrix$Pb_mean)
anova(bd)
permutest(bd)

sim <- simper(beta_matrix, permutations = 999)
summary(sim)







# NMDS plots

set.seed(123)
nmds.jac = metaMDS(beta_matrix, distance = "jaccard")
nmds.bray = metaMDS(beta_matrix, distance = "bray")
#nmds.bray.bi = metaMDS(beta_matrix, distance = "bray", binary = TRUE)

# Show the beta diversity in basic NMDS plots
plot(nmds.jac, type='t')
plot(nmds.bray)
#plot(nmds.bray.bi)



# stress plotting 
stressplot(nmds.jac)
gof <- goodness(nmds.jac)
gof
plot(nmds.jac, display = "sites", type = "n")
points(nmds.jac, display = "sites", cex = 2*gof/mean(gof))

# Find the pollution data
depth_pb <- pollution_matrix %>%  select(depth, Pb_mean, Assemblage)
en = envfit(nmds.jac, depth_pb, permutations = 999, na.rm = TRUE)
plot(en)
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#get the NMDS data output
data.scores.bray = as.data.frame(scores(nmds.bray)$sites)
data.scores.bray$Sample = beta_matrix1$Sample

data.scores.bray.bi = as.data.frame(scores(nmds.bray.bi)$sites)
data.scores.bray.bi$Sample = beta_matrix1$Sample

data.scores.jac = as.data.frame(scores(nmds.jac)$sites)
data.scores.jac$Sample = beta_matrix1$Sample
data.scores.jac$site <- rownames(data.scores.jac)
#species
species.scores.jac <- as.data.frame(scores(nmds.jac, "species"))
species.scores.jac$species <- rownames(species.scores.jac)
species.scores.jac <- cbind(species.scores.jac, pval = en$vectors$pvals) 
signi.species.jac <- subset(species.scores.jac, pval<=0.05)
head(signi.species.jac)

# plot NMDS
nmds.jaccard.plot <- ggplot() + 
  geom_point(data=species.scores.jac,aes(x=NMDS1,y=NMDS2),alpha=0.5)  +
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
  ggtitle("NMDS for Jaccard Dissimilarity Matrix")



## the species driving the dissimiliarity, maybe only for key taxa

nmds.jaccard.plot + 
  geom_segment(data = signi.species.jac, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = signi.species.jac, aes(x=NMDS1, y=NMDS2, label = species), cex = 3, direction = "both", segment.size = 0.25)

## plots
nmds.jaccard.plot

  