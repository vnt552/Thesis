# diversity with iNEXT

## library
library(devtools)
library(iNEXT)
library(ggplot2)
library(tidyverse)

### Choose what you want to look at, everything or key taxa?
alpha_matrix <- metaDMG_results %>% filter(tax_rank=="species") %>% select(sample, tax_name, N_reads)
alpha_matrix <- metaDMG_results %>% filter(tax_rank=="species", grepl("\\bInsecta\\b", tax_path)) %>% select(sample, tax_name, N_reads)


###
alpha_wide_matrix <- pivot_wider(alpha_matrix, names_from = sample , values_from = N_reads)
alpha_wide_matrix[is.na(alpha_wide_matrix)] <- 0
alpha_wide_matrix <- as.data.frame(alpha_wide_matrix)
alpha_hill_matrix <- alpha_wide_matrix[,-1]
alpha_hill_species_list <- alpha_wide_matrix[,1]

hill_q0 <- iNEXT(alpha_hill_matrix, q=0, datatype="abundance", size=NULL, endpoint=NULL)
hill_q1 <- iNEXT(alpha_hill_matrix, q=1, datatype="abundance", size=NULL, endpoint=NULL)
hill_q2 <- iNEXT(alpha_hill_matrix, q=2, datatype="abundance", size=NULL, endpoint=NULL)

hill_q0_q1_q2_insect <- iNEXT(alpha_hill_matrix, q=c(0,1,2), datatype="abundance", size=NULL, endpoint=NULL)


#### Hill numbers on only a subset of data, so plots are easier to look at
sub11_31 <- alpha_hill_matrix %>%  select("GA-11", "GA-19", "GA-31")
sub41_71 <- alpha_hill_matrix %>%  select("GA-41", "GA-61", "GA-71")
sub81_111 <- alpha_hill_matrix %>%  select("GA-81", "GA-91", "GA-111")
sub121_141 <- alpha_hill_matrix %>%  select("GA-121", "GA-131", "GA-141")
sub161_191 <- alpha_hill_matrix %>%  select("GA-161", "GA-181", "GA-191")
sub_cont <- alpha_hill_matrix %>%  select("GA-Cex", "GA-Cli")

hill_sub11_31 <- iNEXT(sub11_31, q=c(0,1,2), datatype="abundance", size=NULL, endpoint=NULL)
hill_sub41_71 <- iNEXT(sub41_71, q=c(0,1,2), datatype="abundance", size=NULL, endpoint=NULL)
hill_sub81_111 <- iNEXT(sub81_111, q=c(0,1,2), datatype="abundance", size=NULL, endpoint=NULL)
hill_sub121_141 <- iNEXT(sub121_141, q=c(0,1,2), datatype="abundance", size=NULL, endpoint=NULL)
hill_sub161_191 <- iNEXT(sub161_191, q=c(0,1,2), datatype="abundance", size=NULL, endpoint=NULL)
hill_subcontr <- iNEXT(sub_cont, q=c(0,1,2), datatype="abundance", size=NULL, endpoint=NULL)
ggiNEXT(hill_sub11_31, type=3, facet.var="Assemblage")
ggiNEXT(hill_sub41_71, type=3, facet.var="Assemblage")
ggiNEXT(hill_sub81_111, type=3, facet.var="Assemblage")
ggiNEXT(hill_sub121_141, type=3, facet.var="Assemblage")
ggiNEXT(hill_sub161_191, type=3, facet.var="Assemblage")
ggiNEXT(hill_sub11_31, type=1, facet.var="Order.q", color.var="Assemblage")



### coverage plots
coverage_plot <- ggiNEXT(hill_subcontr, type=2, facet.var="None", color.var="Assemblage")
coverage_plot + theme_bw()



ggiNEXT(hill_q0_q1_q2_insect, type=1, facet.var="Assemblage")
ggiNEXT(hill_q0_q1_q2_insect, facet.var="Order.q")
ggiNEXT(hill_q0_q1_q2_insect, type=3, facet.var="Assemblage")
coverage_hill <- hill_q0_q1_q2$iNextEst$coverage_based

hill_all$DataInfo

diversity_hill <- hill_q0_q1_q2$AsyEst
diversity_hill$Assemblage <- factor(diversity_hill$Assemblage , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))

diversity_hill_insects <- hill_q0_q1_q2_insect$AsyEst
diversity_hill_insects$Assemblage <- factor(diversity_hill_insects$Assemblage , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191", "GA-Cex", "GA-Cli"))

## plot of diversity, observed and estimated

ggplot(diversity_hill, aes(x=Assemblage, y=Estimator, color="#9ebcda")) + 
  geom_pointrange(aes(ymin=LCL, ymax=UCL)) +
  geom_point(data= diversity_hill, aes(x=Assemblage, y=Observed, color="#fec44f")) +
  scale_color_manual(values = c("#9ebcda", "#fec44f"), labels=c("Asymptotic Diversity Estimate", "Observed")) +
  facet_wrap(~Diversity, ncol=1, scales="free_y") +
  ylab("Diversity") +
  xlab("Sample") +
  ggtitle("Hill Numbers for Insecta", subtitle = "Standardized by Coverage") +
  theme_bw()


### Estimate diversity

div_estimate <- estimateD(alpha_hill_matrix, datatype = "abundance", base = "coverage")
div_estimate



#####
## Plot hilnumbers vs. lead
subdata_diversity <- diversity_hill[c(1:45),]
subdata_diversity$Assemblage <- factor(subdata_diversity$Assemblage , levels=c("GA-11", "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81", "GA-91", "GA-111", "GA-121", "GA-131", "GA-141", "GA-161", "GA-181", "GA-191"))
colnames(subset_pollution) <- c("depth", "Pb", "Cd", "Assemblage")

colnames(mean_pb_pollution) <- c("Assemblage", "Pb_mean")
pollution_matrix <- merge(subset_pollution, mean_pb_pollution, by="Assemblage")  


alpha_pollution_matrix <- merge(pollution_matrix, subdata_diversity, by="Assemblage")  

ggplot(alpha_pollution_matrix, aes(x = Pb_mean, y = Estimator, colour = Assemblage)) +
  geom_point() +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Pb") + 
  theme_bw() +
  facet_wrap(~Diversity, scales = "free_y")


## anova test on Shannon diversity
shannon_diversity <- alpha_pollution_matrix %>%  filter(Diversity=="Shannon diversity")

anova_pb <- aov(Observed ~ Pb_mean, shannon_diversity)
anova_depth<- aov(Observed ~ depth, shannon_diversity)
two_way <- aov(Observed ~Pb_mean + depth, shannon_diversity)
two_interaction <- aov(Observed ~Pb_mean * depth, shannon_diversity)
summary(anova_pb)
summary(anova_depth)
summary(two_way)
summary(two_interaction)

## which test is best?
library(AICcmodavg)

model.set <- list(anova_hill,anova_depth, two_way, two_interaction)
model.names <- c("One.Way.Pb", "One.Way.Depth", "Two.Way", "Two.Interaction")

aictab(model.set, modnames = model.names)


######
#### Hill numbers as a function of q
hill_all <- iNEXT(alpha_hill_matrix, q=c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), datatype="abundance", size=NULL, endpoint=NULL)
hill_all_order.q <- hill_all$iNextEst$coverage_based

hill_order_observed <- hill_all_order.q %>% filter(Method=="Observed")
labels_q0 <- hill_order_observed %>%  filter(Order.q==0.00)

## plotting
ggplot(hill_order_observed, aes(x=Order.q, y=qD, group = Assemblage, color=Assemblage)) +
  geom_line() +
  geom_label_repel(data=labels_q0, aes(x= Order.q, y=qD, label=Assemblage), size=2, hjust=0,nudge_x= -0.5, direction="y") +
  theme_bw() +
  xlab("Order of q") +
  ylab("Diversity Estimate") +
  ggtitle("Diversity Estimate by Order of q", subtitle = "Observed Data")
  

