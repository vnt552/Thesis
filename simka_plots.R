mat_abundance_braycurtis <- read.csv("mat_abundance_braycurtis.csv", sep=";")
mydataALL <- mat_abundance_braycurtis 

mat_presenceAbsence_braycurtis <- read.csv("mat_presenceAbsence_braycurtis.csv", sep=";")
mydata_PRES <- mat_presenceAbsence_braycurtis

SimkaMetaData <- read.csv("simka_list.txt", sep=":", quote="\"", comment.char="", header = FALSE)
SimkaMetaData <- SimkaMetaData %>%  mutate(group = case_when(grepl("GA-C", V1) ~ "Control",
                                                             grepl("1", V1) ~"Sample"))
##################### replace NA with 0's ###########################

mydata_PRES[is.na(mydata_PRES)]=0 #remove NA
mydataALL[is.na(mydataALL)]=0

############################### renaming the colnames #######################

#colnames(mydata) <- c( "taxa","Layer_2B", "Layer_2P","Layer_3B", "Layer_3P","Layer_4B", "Layer_1B", "Layer_1P","P_blank", "B_blank","Lib_blank1","Lib_blank2" )

##################### create percentage table #######################
b1=as.matrix(mydata_PRES[,seq(2,18)])  # change number of coloumns to the total of mydata
rownames(b1)=mydata_PRES[,1]

a1=as.matrix(mydataALL[,seq(2,18)])  # change number of coloumns to the total of mydata
rownames(a1)=mydataALL[,1]
#colnames(b1)=lapply(colnames(b1),function(x){strsplit(x,'\\.')[[1]][1]})
#b2 <- prop.table(b1, margin=2)*100 # makes proportion table, needs 2 margins e.g. header and 1st row names
#colSums(prop.table(b1, margin=2)*100) # should give 100 for each coloumn

#####################Creating subset of table###################################################

###https://www.statmethods.net/management/subset.html 

######### removing only rows if no values are above x in any box in coloumn 1-4 (can be changed)
#tmp <- b2[apply(b2[,1:11], MARGIN = 1, function(x) any(x > 0.5)), ]



##############################

library(ggplot2)
#https://statquest.org/2017/12/18/statquest-mds-and-pcoa-in-r/ 

data.matrix <- b1
data.matrix2 <- a1
###################################################################
##
## 1) Just for reference, draw a PCA plot using this data...
##
###################################################################
pca <- prcomp(t(data.matrix), scale=TRUE, center=TRUE) 
pca1 <- prcomp(t(data.matrix2), scale=TRUE, center=TRUE) 
## calculate the percentage of variation that each PC accounts for...
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.var.per

pca1.var <- pca1$sdev^2
pca1.var.per <- round(pca1.var/sum(pca1.var)*100, 1)
pca1.var.per
## now make a fancy looking plot that shows the PCs and the variation:
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

pca1.data <- data.frame(Sample=rownames(pca1$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca1.data

md1 <- factor(SimkaMetaData$V1)
shp <- factor(SimkaMetaData$group)
pca$X <- as.numeric(as.character(pca$X))
pca$Y <- as.numeric(as.character(pca$Y))

pca1$X <- as.numeric(as.character(pca1$X))
pca1$Y <- as.numeric(as.character(pca1$Y))
require(ggplot2)
#group your points with shp and md1 above
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(colour = md1, shape = shp), size = 5) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_classic() +
  ggtitle("PCA Simka Present Absence")+
  geom_text()

#group your points with shp and md1 above
ggplot(data=pca1.data, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(colour = md1, shape = shp), size = 5) +
  xlab(paste("PC1 - ", pca1.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca1.var.per[2], "%", sep="")) +
  theme_classic() +
  ggtitle("PCA_Presence_Absence_BrayCurtis")+
  geom_text()


heatmap(b1)
b2 <- b1
b2[b2 == 0] <- NA

heatmap(b2)

heatmap(a1)
a2 <- a1
a2[a2 == 0] <- NA

heatmap(a2)

