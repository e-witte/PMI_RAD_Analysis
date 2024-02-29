####################
## ADEGENET 2.1.1 ##
##      DAPC      ##
####################
rm(list = ls())
install.packages("adegenet", dep=TRUE)
install.packages("ape", dep=TRUE)
install.packages("RColorBrewer", dep=TRUE)
install.packages("vcfR", dep=TRUE)

library("adegenet")
library("ape") 
library("RColorBrewer")

#show Adegenet version
packageDescription("adegenet", fields = "Version") 

#Load .vcf data output from populations (STACKS)
#install.packages("vcfR")
library(vcfR)
#load vcf file
vcf <- read.vcfR(file = "data/21324_populations.snps.vcf")
#vcf <- read.vcfR(file = "populations_allsamples_Indian_Car_PCADAPT.snps.vcf")

#vcf <- read.vcfR(file = "LGI_67_all_loci_W_PCADAPT_95_populations.snps.vcf")
#transform vcf to genlight object (genlight is better at storing large SNP data than genind)
data <- vcfR2genlight(vcf)
data

#__________________
#Label populations 
#------------------
#check what order your data is in
data$ind.names

#label your populations 
pop(data)<- c(rep("RED",10),rep("MED",56)) #stores population information in Adegenet



#not sure about these... but they don't affect
Dgen <- dist.gene(data@gen)
Dgeo <- dist(data)
ploidy(data) <- 2

##############
# Check data #
##############
#plot missing data
glPlot(data, posi="topleft")
glPlot(data, col=bluepal(6))

#plot allele frequencies
myFreq <- glMean(data)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies", main="Distribution of (second) allele frequencies")
temp <- density(myFreq) 
lines(temp$x, temp$y*1.8,lwd=3)

######################
#Principal Components#
######################
#Run PCA
data.pca <- glPca(data)
barplot(100*data.pca$eig/sum(data.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

scatter(data.pca, posi="bottomleft", cex=0.5, alpha=0.5)
title("PCA of data\n axes 1-2")

#Assess PCA using colors
myCol <- colorplot(data.pca$scores,data.pca$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey") 
add.scatter.eig(data.pca$eig[1:40],2,1,2, posi="bottomright", inset=.1, ratio=.2)
legend(10,-4, col=c("orange","red", "blue", "green"))

#View PCA by population
#To view the results of the PCA we can use the package ggplot2. 
#We need to convert the data frame that contains the principal components (data.pca$scores) 
#into the new object data.pca.scores. In addition, we will add the population values as a new 
#column in our data.pca.scores object, in order to be able to color samples by population.

data.pca.scores <- as.data.frame(data.pca$scores)
data.pca.scores$pop <- pop(data)

cols <- brewer.pal(n = nPop(data), name = "Dark2")
cols <- c(brewer.pal(n = 7, name="Dark2"), brewer.pal(name="Paired", n=1))

#check percent of explained variance per PCA eigenvalues
100*data.pca$eig/sum(data.pca$eig)
#first two values will be the two PC plotted
#add percentage in below labels

library(ggplot2)
set.seed(9)
plot_PMI <- ggplot(data.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=1.5, alpha=0.8) +
  stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()+
  theme(legend.title=element_blank())+
  labs(x = "PC1 (6.94%)", y = "PC2 (4.14%)")

#save this file

pdf(plot_PMI, file = "figs/PMI_ALL_LOCI_aes12.png")

save.image(file = "figs/PMI_all_loci.png")

library(ggplot2)
set.seed(9)
ggplot(data.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=1.5, alpha=0.8) +
  stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()+
  theme(legend.title=element_blank())+
  labs(x = "PC1 (13.28%)", y = "PC2 (2.98%)")

dev.off()

########
# DAPC #
########
#run DAPC
data.dapc <- dapc(data)#, n.pca = 5, n.da = 3)
#if you want to save the DAPC data to not have to re-run the analysis you can do
#save(data.dapc, file = "LGI_df_DAPC.Rdata")


#create color palette
cols <- brewer.pal(n = nPop(data), name = "Dark2")

#plot
scatter(data.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, cstar=1, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, ratio.pca = 0.2, ratio.da = 0.2)

#save this file

pdf(file = "DAPC_NCY_Indian_invasives_PCAadapt_allsamples.pdf")

scatter(data.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, cstar=1, scree.pca = FALSE, cleg = 0.75, ratio.pca = 0.2, ratio.da = 0.2)

dev.off()


#-------------------------------------------------------------------------
#Assess if there are alleles that differentiate your clusters per axis
library(ade4)
set.seed(4)
contrib <- loadingplot(data.dapc$var.contr, axis = 2, thres = 0.0035, lab.jitter = 1)
#change axis to correspond to the DAPC axis that drives more differences in your data
#change the threshold accordingly 

#take a closer look at individual loci per population
temp    <- seploc(data)       # seploc {adegenet} creates a list of individual loci.
snp11925  <- tab(temp[["11925"]]) # tab {adegenet} returns a matrix of genotypes
snp5793  <- tab(temp[["5793"]])

# The following two commands find the average allele frequencies per population
(freq906 <- apply(snp501, 2, function(e) tapply(e, pop(data), mean, na.rm = TRUE)))
#------------------------------------------------------------------------


# find the number of clusters in your dataset
grp <- find.clusters(data, max.n.clust=8)
table(pop(data), grp$grp)
table.value(table(pop(data), grp$grp), col.lab=paste("inf", 1:6), row.lab=paste("ori", 1:6))

#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(data.dapc,col = cols, posi = 'top')

#make a graph

pdf(file = "complot_6_5.pdf")

compoplot(data.dapc,col = cols)

dev.off()

#transform data to comoplot in ggplot
dapc.results <- as.data.frame(data.dapc$posterior)
dapc.results$pop <- pop(data)
dapc.results$indNames <- rownames(dapc.results)


library(reshape2)
dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = cols) +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))+
  labs(fill="Assigned Population")

###################################
#Cross Validation: DAPC analysis ##
###################################
# 
set.seed(4)
contrib <- loadingplot(data.dapc$var.contr, axis = 2, thres = 0.002, lab.jitter = 1)
#---------
data.gind <- vcfR2genind(vcf)
#SRE_8x_qc
pop(data.gind)<- c(rep("Corredor",5), rep("Los Cabos",1), rep("Bahia Magdalena",8), rep("Los Cabos",5), rep("Loreto",9), rep("Cabo Pulmo",8), rep("Isla Ángel de la Guarda",8), rep("Isla San Pedro Mártir",10)) #stores population information in Adegenet

freq6600_24 <- tab(genind2genpop(data.gind[loc=c("6600_24")]),freq=TRUE)
freq120092_59 <- tab(genind2genpop(data.gind[loc=c("120092_59")]),freq=TRUE)

par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,.1),las=3)
matplot(freq6600_24, pch=c("a","c"), type="b",
        xlab="population",ylab="allele frequency", xaxt="n",
        cex=1.5, main="SNP # 6600_24")
axis(side=1, at=1:7)#, lab=2001:2006)
matplot(freq120092_59, pch=c("c","t"), type="b", xlab="year",
        ylab="allele frequency", xaxt="n", cex=1.5,
        main="SNP # 120092_59")
axis(side=1, at=1:7)#, lab=2001:2006)

#-------------
#few PC's
temp <- summary(dapc(data, n.da=3, n.pca=30))$assign.per.pop*100
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual breed",
        horiz=TRUE, las=1)

#too many PC's
temp <- summary(dapc(data, n.da=3, n.pca=50))$assign.per.pop*100
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual breed",
        horiz=TRUE, las=1)

#subsample
x <- data
pop(x) <- sample(pop(x))
temp <- summary(dapc(x, n.da=3, n.pca=50))$assign.per.pop*100
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual breed", horiz=TRUE, las=1)
#populations have been randomized yet they are all 100% assigned to the same populations

#Check a-score
dapc2 <- dapc(data, n.da=3, n.pca=20)
temp <- a.score(dapc2)
names(temp)

temp$tab[1:5,1:5]
temp

dapc2 <- dapc(data, n.da=10, n.pca=50)
#obtain optimal 
temp <- optim.a.score(dapc2)

#Cross-validation
x <- data
mat <- tab(x, NA.method="mean")
grp <- pop(x)
xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)

xval[2:6]


