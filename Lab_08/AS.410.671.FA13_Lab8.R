# AS.410.671.81.FA13 – Lab #8
# Cluster Analysis
# Author: Shan Sabri

#######################
#          1          #
#######################

library(fibroEset)
data(fibroEset)

dat <- fibroEset
dat$species


#######################
#          2          #
#######################

dat.sub <- dat[sample(featureNames(dat), 50), ]

colnames(exprs(dat.sub)) <- as.character(dat.sub$species)

#######################
#          3          #
#######################

d.manhat <- dist(t(exprs(dat.sub)), method = "manhattan")
hc       <- hclust(d.manhat, method = "median")

plot(
	hc, labels = colnames(exprs(dat.sub)), sub="", xlab="", 
	main = "Hierarchical Clustering of 50 random genes from the fibroEset data"
	) 

 
#######################
#          4          #
#######################

hm.rg <- c(
	"#FF0000", "#CC0000", "#990000", "#660000", "#330000", "#000000",
	"#000000", "#0A3300", "#146600", "#1F9900", "#29CC00", "#33FF00"
	)
	
heatmap(
	as.matrix(exprs(dat.sub)),
	col = hm.rg,
	main = "Heatmap of 50 random genes from the fibroEset data",
	xlab = "Sample", 
	ylab = "Probeset"
	)

 
#######################
#          5          #
#######################

library(kernlab)

dat.kpca <- kpca(t(exprs(dat.sub)), kernel = "rbfdot", kpar = list(sigma = 0.002), features = 2)
pcv      <- pcv(dat.kpca)
rot      <- rotated(dat.kpca)
pcv.k    <- kmeans(pcv, centers = 3, iter.max = 20)
rot.k    <- kmeans(rot, centers = 3, iter.max = 20)


#######################
#          6          #
#######################

par(mfrow = c(2, 1))
plot(
	pcv, col = pcv.k$cluster, cex = 1,
	main = "PCA Scatter Plot with PC vectors (column-wise)\nk = 3",
	xlab = "P1", 
	ylab = "P2"
	)
points(pcv.k$centers, col = 1:3, pch = "*", cex = 2.5)

plot(
	rot, col = rot.k$cluster, cex = 1,
	main = "PCA Scatter Plot with PC Projected vectors\nk = 3",
	xlab = "P1", 
	ylab = "P2"
	)
points(rot.k$centers, col = 1:3, pch = "*", cex = 2.5)