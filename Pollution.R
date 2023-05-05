### pollution plot
pollution_data <- read_csv("pollution2000.csv")


# only data from the samples
subset_pollution <- subset(pollution_data, depth %in% c(110, 195, 310, 410, 610, 710, 810, 910, 1110, 1210,1310, 1410, 1610, 1810, 1910))
subset_pollution$sample <- c("GA-11" , "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81","GA-91", "GA-111", "GA-121", "GA-131","GA-141", "GA-161", "GA-181", "GA-191")

large_subset_pollution <- subset(pollution_data, depth %in% c(110:115, 195:200, 310:315, 410:415, 610:615, 710:715, 810:815, 910:915, 1110:1115, 1210:1215,1310:1315, 1410:1415, 1610:1615, 1810:1815, 1910:1915))
large_subset_pollution$Assemblage <- rep(c("GA-11" , "GA-19", "GA-31", "GA-41", "GA-61", "GA-71", "GA-81","GA-91", "GA-111", "GA-121", "GA-131","GA-141", "GA-161", "GA-181", "GA-191"),each=6)
large_subset_pollution$Pb <- as.numeric(as.character(large_subset_pollution$Pb))
mean_pb_pollution <- aggregate(large_subset_pollution[, "Pb"], list(large_subset_pollution$Assemblage), mean)
mean_pb_pollution$Group.1 <-factor(mean_pb_pollution$Group.1 , levels= c("GA-191", "GA-181", "GA-161", "GA-141", "GA-131", "GA-121", "GA-111", "GA-91", "GA-81", "GA-71", "GA-61", "GA-41", "GA-31", "GA-19", "GA-11"))


# pollution plot
ggplot(pollution_data, aes(x=Pb, y=depth, group=1)) +
  geom_line(orientation = "y",size=0.05) +
  scale_x_continuous() +
  scale_y_reverse(breaks=c(0,110, 190, 310, 410, 610, 710, 810, 910, 1110, 1210,1310, 1410, 1610, 1810, 1910), 
                  labels = c("0","11", "19", "31", "41", "61", "71", "81","91", "111", "121", "131","141", "161", "181", "191")) +
  ylab("Depth (cm)") +
  theme_classic()



ggplot(mean_pb_pollution, aes(y=Group.1, x=Pb)) +
  geom_col() +
  ylab("Sample") +
  theme_bw() +
  ggtitle("Pb concentration throughout the sample", subtitle = "Mean of 6mm")
