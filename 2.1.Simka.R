# SIMKA Presence/Absence
 
mat_presenceAbsence_braycurtis <- read.csv("simka2/mat_presenceAbsence_braycurtis.csv", sep=";")
mydata_PRES <- mat_presenceAbsence_braycurtis

SimkaMetaData <- read.csv("simka2/simka_list.txt", sep=":", quote="\"", comment.char="", header = FALSE)
SimkaMetaData <- SimkaMetaData %>%  mutate(group = case_when(grepl("GA-C", V1) ~ "Control",
                                                             grepl("1", V1) ~"Sample"))
##################### replace NA with 0's ###########################

mydata_PRES[is.na(mydata_PRES)]=0 #remove NA

##################### create percentage table #######################
b1=as.matrix(mydata_PRES[,seq(2,18)])  # change number of coloumns to the total of mydata
rownames(b1)=mydata_PRES[,1]

data.matrix <- b1
pca <- prcomp(t(data.matrix), scale=TRUE, center=TRUE) 

## calculate the percentage of variation that each PC accounts for
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.var.per


## now make a fancy looking plot that shows the PCs and the variation:
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,3])
pca.data

pca1.data <- data.frame(Sample=rownames(pca1$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
                       
md1 <- factor(SimkaMetaData$V1)
shp <- factor(SimkaMetaData$group)
pca$X <- as.numeric(as.character(pca$X))
pca$Y <- as.numeric(as.character(pca$Y))

# PCA plot
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(colour = md1, shape = shp), size = 5) +
  xlab(paste("PC3 - ", pca.var.per[3], "%", sep="")) +
  ylab(paste("PC4 - ", pca.var.per[4], "%", sep="")) +
  theme_classic() +
  labs(col="Sample", shape="Type") +
  ggtitle("PCA Simka Present Absence")+
  geom_text()
  
  
# Heat map
heatmap(b1)
b2 <- b1
b2[b2 == 0] <- NA
