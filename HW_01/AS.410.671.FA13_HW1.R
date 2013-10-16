# AS.410.671.81.FA13 – HW #2
# Author: Shan Sabri

#######################
#          1          #
####################### 

setwd("C:\\Users\\Shan\\Desktop\\JHU_Fall2013\\410.671_Microarrays&Analysis\\DataSets\\") 
dat <- read.table("renal_cell_carcinoma.txt", header = T, row.names = 1) 
 
dim(dat)

	
#######################
#          2          #
####################### 

ann        <- read.table("renal_carcinoma_annotation.txt", header = F) 
identity   <- as.character(ann[, 9]) 
names(dat) <- paste(names(dat), identity, sep = "_") 
 
dimnames(dat)[[2]] 

	
#######################
#          3          #
####################### 

# Correlation Plot

library(gplots) 
dat.cor <- cor(dat) 
layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2), 5, 2, byrow = TRUE)) 
par(oma=c(5, 7, 1, 1)) 
cx  <- rev(colorpanel(25, "yellow", "black", "blue")) 
leg <- seq(min(dat.cor, na.rm = T), max(dat.cor, na.rm = T), length = 10) 
image(dat.cor, main= "Correlation Plot (Heat Map) – Normal v. Tumor", axes = F, col = cx) 
axis(1, at = seq(0, 1, length=ncol(dat.cor)), label = dimnames(dat.cor)[[2]], cex.axis = 
0.9, las = 2) 
axis(2, at = seq(0, 1, length=ncol(dat.cor)), label = dimnames(dat.cor)[[2]], cex.axis = 
0.9, las = 2) 
 
par(mar = rep(2, 4)) 
image(as.matrix(leg), col = cx, axes = T) 
tmp <- round(leg, 2) 
axis(1, at = seq(0, 1, length = length(leg)), labels = tmp, cex.axis = 1) 

# Hierarchical clustering dendrogram

dat.trans       <- t(dat) 
dat.trans.dist  <- dist(dat.trans, method = "Euclidean") 
dat.trans.clust <- hclust(dat.trans.dist, method = "Single") 
 
plot(dat.trans.clust, labels = names(dat.trans), cex = 0.75) 

# CV vs. mean plot 

dat.mean <- apply(log2(dat), 2, mean) 
dat.sd   <- sqrt(apply(log2(dat), 2, var)) 
dat.cv   <- dat.sd/dat.mean 
 
plot(dat.mean, dat.cv, main = "cRCC Dataset\nCV v. Mean", xlab = "Mean", ylab = "CV", col = "Blue", cex = 1.5, type = "n") 
points(dat.mean, dat.cv, bg = "Red", col = 1, pch = 21) 
text(dat.mean, dat.cv, label = dimnames(dat)[[2]], pos = 1, cex = 0.5) 

# Average correlation plot

dat.cor <- cor(dat) 
dat.avg <- apply(dat.cor, 1, mean) 
 
par(oma = c(3, 0.1, 0.1, 0.1)) 
plot(c(1, length(dat.avg)), range(dat.avg), type = "n", xlab = "", ylab = "Average Correlation (r)", main = "Average Correlation Plot of Tumor v. Normal samples", axes = F) 
points(dat.avg, bg = "Blue", col = 1, pch = 21, cex = 1.25) 
axis(1, at = c(1:length(dat.avg)), labels = dimnames(dat)[[2]], las = 2, cex.lab = 0.4, 
cex.axis = 0.6) 
axis(2) 
abline(v = seq(0.5, 62.5, 1), col = "Grey") 


#######################
#          4          #
####################### 

source("http://bioconductor.org/biocLite.R") 
biocLite("impute") 
 
library(impute) 


#######################
#          5          #
####################### 

dat <- dat[, -c(10, 19)] # Removing outliers 'GSM146798_Normal' & 'GSM146799_Tumor' 
dim(dat) 
dimnames(dat)[[2]] 


#######################
#          6          #
#######################

KNG1 <- dat["206054_at", ] 
KNG2 <- dat["217512_at", ] 
AQP  <- dat["206672_at", ] 
all.probesets <- rbind(KNG1, KNG2, AQP) 
KNG.probesets <- rbind(KNG1, KNG2) 
 
par(mfrow=c(3, 1)) 
 
plot(c(1, ncol(dat)), range(all.probesets), type = 'n', main="Profile plot of KNG1 and AQP2", xlab = "Samples", ylab = "Expression Intensity", axes = F) 
axis(side=1, at=c(1:20), labels = dimnames(probesets)[[2]], cex.axis = 0.4, las = 2) 
axis(side=2) 
for(i in 1:length(all.probesets)) { 
	dat.all <- as.numeric(all.probesets[i, ]) 
	lines(c(1:ncol(all.probesets)), dat.all, col = i, lwd = 2) 
} 
 
plot(c(1, ncol(dat)), range(KNG.probesets), type = 'n', main="Profile plot of KNG1", xlab = "Samples", ylab = "Expression Intensity", axes = F) 
axis(side=1, at=c(1:20), labels = dimnames(probesets)[[2]], cex.axis = 0.4, las = 2) 
axis(side=2) 
for(i in 1:length(KNG.probesets)) { 
	dat.KNG <- as.numeric(KNG.probesets[i, ]) 
	lines(c(1:ncol(KNG.probesets)), dat.KNG, col = i, lwd = 2) 
} 
 
plot(c(1, ncol(dat)), range(AQP), type = 'n', main="Profile plot of AQP2", xlab = "Samples", ylab = "Expression Intensity", axes = F) 
axis(side=1, at=c(1:20), labels = dimnames(probesets)[[2]], cex.axis = 0.4, las = 2) 
axis(side=2) 
for(i in 1:length(AQP)) { 
	dat.AQP <- as.numeric(AQP[i, ]) 
	lines(c(1:ncol(AQP)), dat.AQP, col = i, lwd = 2) 
} 


#######################
#          7          #
#######################

dat.KNG <- as.matrix(dat) 
actual  <- dat.KNG["206054_at", "GSM146784_Normal"] 
 
dat.KNG["206054_at", "GSM146784_Normal"] <- NA 
dat.KNG["206054_at", "GSM146784_Normal"] 


#######################
#          8          #
#######################

dat.knn <- impute.knn(dat.KNG, 6) 
dat.knn$data["206054_at", "GSM146784_Normal"] 
dat.knn <- dat.knn$data 
 
expected <- dat.knn["206054_at", "GSM146784_Normal"] 


#######################
#          9          #
#######################

RR <- abs(expected-actual)/actual 
RR 
percentRR <- RR * 100 
percentRR 


#######################
#          10         #
#######################

source("http://bioconductor.org/biocLite.R") 
biocLite("pcaMethods") 
library(pcaMethods) 
 
dat.svd <- as.matrix(dat) 
dat.svd["206054_at", "GSM146784_Normal"] <- NA 
dat.pca <- pca(dat.svd, method = "svdImpute", nPcs = 9) 
dat.pca.result <- completeObs(dat.pca) 
dat.pca.result["206054_at", "GSM146784_Normal"] 


#######################
#          11         #
#######################

actual <- dat["206054_at", ] 
pca <- dat.pca.result["206054_at", ] 
knn <- dat.knn["206054_at", ] 
dat.pca.knn <- rbind(actual, knn, pca) 
 
plot(c(1, ncol(dat)), range(dat.pca.knn), type = "n", main = "Profile plot of probeset 206054_at\nwith actual and imputed values of GSM146784_Normal", xlab = "Samples", ylab = "Expression Intensity", axes = F) 
axis(side = 1, at = 1:length(dat), labels = dimnames(dat.pca.knn)[[2]], cex.axis = 0.4, las = 2) 
axis(side = 2) 
 
for(i in 1:length(dat.pca.knn)) { 
	dat.y <- as.numeric(dat.pca.knn[i, ]) 
	lines(c(1:ncol(dat.pca.knn)), dat.y, col = i, lwd = 2) 
} 
