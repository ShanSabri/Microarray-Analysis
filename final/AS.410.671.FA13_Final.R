# AS.410.671.81.FA13 – Final Exam
# Author: Shan Sabri

# *** Primary Reference *** : http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3295

# Title   :	
#	 - Age-associated aneuploidy model	Cluster Analysis
# Summary :	
#	- Analysis of germinal vesicle (GV)-intact oocytes and metaphase II (MII) eggs obtained from young 
#	  and old females. An increased incidence of aneuploidy is observed with increasing maternal age. 
#	  Results provide insight into the molecular basis underlying the age-associated increase in aneuploidy.
# Organism:	
#	- Mus musculus
# Platform:	 
#	- GPL1261: [Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array
# Citation:	
#	- Pan et al. "Age-associated increase in aneuploidy and changes in gene expression in mouse eggs."
#	- Dev Biol 2008 Apr 15;316(2):397-407. 
#	- PMID: 18342300
#	- Reference Series: GSE11667

library(GEOquery)
library(stats)
library(outliers)
library(EMA)
library(e1071)
library(MASS)
library(multtest)
library(kernlab)
library(scatterplot3d)

# Must load SVM-RFE.r from source to perform SVM feature ranking.
# Source can be downloaded to local machine via:
# https://github.com/johncolby/SVM-RFE/blob/master/msvmRFE.R
setwd("C:\\Users\\Shan\\Desktop\\Final\\")
source(file = "code\\msvmRFE.R")

# Reading in GEO Data 
gds      <- getGEO("GDS3295")
exp.data <- Table(gds)

# Reading in the corresponding GEO annotation file
ann.file  <-  "GPL1261.annot"
ann.data  <- read.delim(ann.file, header = T, row.names = 1, skip = 27, sep = "\t")
ann.data  <- ann.data[1:nrow(exp.data), ]
gene.info <- data.frame(Description = ann.data$Gene.title, Symbol = ann.data$Gene.symbol)
rownames(gene.info) <- rownames(ann.data)

# Appending development stage (GV/MII) to corresponding samples
for (i in 3:10)  { names(exp.data)[i] <- paste("GV",  names(exp.data[i]), sep = "_") }
for (j in 11:18) { names(exp.data)[j] <- paste("MII", names(exp.data[j]), sep = "_") }

samp.matrix           <- data.matrix(exp.data[, (3:ncol(exp.data))])
rownames(samp.matrix) <- rownames(ann.data)

samp.count     <- ncol(samp.matrix)
profile.count  <- nrow(samp.matrix)


# Descriptive Satistics
data.stdev    <- apply(samp.matrix, 1, sd, na.rm = TRUE)
data.rowMeans <- rowMeans(samp.matrix, na.rm = TRUE)
# mean(data.stdev) ; mean(data.rowMeans)
# IQR(data.stdev) ; IQR(data.rowMeans)

# Histograms illistrating initial spread
png("Histogram_data.rowMeans.png")
hist(
	data.rowMeans, 
	col  = "Red",
	xlab = "Mean expression value for GV/MII samples",
	ylab = "Frequency",
	main = paste("Histogram of mean expression values for",profile.count,"profiles")
	)
dev.off()

png("Histogram_data.stdev.png")
hist(
	data.stdev, 
	col  = "Blue",
	xlab = "Standard Deviation expression value for GV/MII samples",
	ylab = "Frequency",
	main = paste("Histogram of standard deviation expression values for",profile.count,"profiles")
	)
dev.off()

########################
#   Outlier Analysis   #
########################

# Correlation Matrix
cor.matrix  <- cor(samp.matrix, method = "pearson", use = "pairwise.complete.obs")

color <- c("#FF0000","#CC0000","#990000","#660000","#330000","#000000",
		   "#000000","#0A3300","#146600","#1F9900","#29CC00","#33FF00")

png("Heatmap_cor.matrix.png")		   
heatmap(
	cor.matrix, 
	col = color, 
	scale = "column", 
	xlab = "Sample Data", 
	ylab="Sample Data", 
	main = "Heatmap illistrating expression values of GV/MII samples"
	)
dev.off()

# Column Means v. Column Variances 
col.mean <- apply(log2(samp.matrix), 2, mean) 
col.var  <- apply(log2(samp.matrix), 2, var)  
cv       <- col.var / col.mean

png("Scatterplot_ColMeansCV.png")
plot(
	col.mean, 
	cv, 
	xlab = "log2(ColMean)",
	ylab = "log2(CV)",
	main = "Plot of Column Mean v. Column Variance for GV/MII samples",
	col  = c(rep("Green", samp.count/2), rep("Blue", samp.count/2)),
	pch  = c(rep(17, samp.count/2), rep(19, samp.count/2))
	)
legend("topright", c("GV oocyte", "MII oocyte"), pch = c(17, 19), col = c("Green", "Blue"))
text(col.mean, cv, labels = names(col.mean), cex = 0.5, offset = 10)
dev.off()

# Row Means v. Column Variances 
row.mean <- apply(log2(samp.matrix), 1, mean) 
row.var  <- apply(log2(samp.matrix), 1, var)  
r.cv     <- row.var / row.mean

png("Scatterplot_RowMeansCV.png")
plot(
	row.mean, 
	r.cv, 
	xlab = "log2(RowMeans)",
	ylab = "log2(CV)",
	main = "Plot of Row Mean v. Row Variance for GV/MII samples",
	col  = c(rep("Green", samp.count/2), rep("Blue", samp.count/2)),
	pch  = c(rep(17, samp.count/2), rep(19, samp.count/2))
	)
legend("topright", c("GV oocyte", "MII oocyte"), pch = c(17, 19), col = c("Green", "Blue"))
abline(v = 0, col = 2, lwd = 2)
dev.off()


# Correlation plot of Row Averages
cor.means <- apply(cor.matrix, 1, mean)

png("Scatterplot_CorMeans_2groups.png")
plot(
	c(1,length(cor.means)), 
	range(cor.means), 
	type = "n", 
	xlab = "",
	ylab = "Average correlation",
	main = "Avg correlation for GV/MII samples",
	axes = FALSE
	)
points(
	cor.means,
	col = c(rep("Green", samp.count/2), rep("Blue", samp.count/2)),
	pch = c(rep(17, samp.count/2), rep(19, samp.count/2))
	)
axis(1, at=c(1:length(cor.means)), labels = colnames(samp.matrix), las = 2, cex.lab = 0.4, cex.axis = 0.6)
axis(2)
grid(nx = 16, col = "grey")
legend(
	"topright", 
	c("GV oocyte", "MII oocyte"), 
	pch = c(17, 19), col = c("Green", "Blue"), bg = "white"
	)
dev.off()

png("Scatterplot_CorMeans_4groups.png")
plot(
	c(1,length(cor.means)), 
	range(cor.means), 
	type = "n", 
	xlab = "",
	ylab = "Average correlation",
	main = "Avg correlation for GV/MII samples with age",
	axes = FALSE
	)
points(
	cor.means,
	col = c(rep("Green", samp.count/4), rep("Blue", samp.count/4), rep("Red", samp.count/4), rep("Black", samp.count/4)),
	pch = c(rep(16, samp.count/4), rep(17, samp.count/4), rep(18, samp.count/4), rep(19, samp.count/4))
	)
axis(1, at=c(1:length(cor.means)), labels = colnames(samp.matrix), las = 2, cex.lab = 0.4, cex.axis = 0.6)
axis(2)
grid(nx = 16, col = "grey")
legend(
	"bottomleft", 
	c("GV oocyte - 66 weeks", "GV oocyte - 6 weeks", "MII oocyte - 66 weeks", "MII oocyte - 6 weeks"), 
	pch = c(16:19), col = c("Green", "Blue", "Red", "Black"), bg = "white"
	)
dev.off()

# Identifying Outlier(s) via outlier()
o       <- cor.means <=  outlier(cor.means)
outlier <- cor.means[o]
cat(sprintf("%s Outlier(s) identified!\n", length(outlier)))
outlier

# Remove Outlier(s) -- Not needed given the high correlation value (0.8955532)
data.no.outliers <- samp.matrix[, -(grep(names(outlier), colnames(samp.matrix)))]
# Note: this can also be accomplished using rm.outlier() 
# data.no.outliers <- rm.outlier(cor.means, fill = FALSE, median = FALSE, opposite = FALSE)

########################
#     Filter Genes     #
########################

quantile(log2(rowMeans(samp.matrix)))
# Output:
#        0%       25%       50%       75%      100% 
# -1.590414  3.508984  5.021335  6.880584 13.003654 

#################
# [ Stage 1/2 ] #
#################

# Eliminating probes with rowMeans less than 0 on a log2 scale
samp.matrix.filtered <- subset(samp.matrix, log2(rowMeans(samp.matrix)) > 0)
removed              <- nrow(samp.matrix) - nrow(samp.matrix.filtered)
cat(sprintf("%s probes removed with rowMeans < 0 on a log2 scale\n", removed))

# Use expFilter() to fine filter genes with low expression values 
# This step is essentially a fail-safe and not necessarily needed
# A gene is kept if at least 0.01*ncol(samp.matrix) of its values is higher than threshold.
dat.fil <- expFilter(log2(samp.matrix.filtered), graph = TRUE)
dat.fil <- subset(dat.fil, rowMeans(dat.fil) > 0)
num.lowexp <- nrow(samp.matrix.filtered) - nrow(dat.fil)
cat(sprintf("%s gene(s) identified and removed for low expression\n", num.lowexp))

# Row Means v. Column Variances on filtered data
fil.mean <- apply(dat.fil, 1, mean) 
fil.var  <- apply(dat.fil, 1, var)  
f.cv     <- fil.var / fil.mean

# Plotting filtered genes (Stage 1/2)
png("Scatterplot_RowMeansCV_Filtered_Stage1.png")
plot(
	fil.mean, 
	f.cv, 
	xlab = "log2(RowMeans)",
	ylab = "log2(CV)",
	main = "Plot of Row Mean v. Row Variance for Filtered GV/MII samples",
	col  = c(rep("Green", samp.count/2), rep("Blue", samp.count/2)),
	pch  = c(rep(17, samp.count/2), rep(19, samp.count/2))
	)
legend("topright", c("GV oocyte", "MII oocyte"), pch = c(17, 19), col = c("Green", "Blue"))
abline(v = 3, col = 2, lwd = 2) # Threshold determined for stage 2/2 of the filtering process
dev.off()

#################
# [ Stage 2/2 ] #
#################

# Eliminating probes with rowMeans less than 3 on a log2 scale
dat.filtered <- subset(dat.fil, rowMeans(dat.fil) > 3)
removed.2    <- nrow(dat.fil) - nrow(dat.filtered)
cat(sprintf("%s probes removed with rowMeans < 0 on a log2 scale\n", removed.2))

fil.mean.2 <- apply(dat.filtered, 1, mean) 
fil.var.2  <- apply(dat.filtered, 1, var)  
fil.cv.2   <- fil.var.2 / fil.mean.2

png("Scatterplot_RowMeansCV_Filtered_Stage2.png")
plot(
	fil.mean.2, 
	fil.cv.2,
	xlim = c(0, 12), ylim = c(0, 50), # Consistant window Size w/ previous plot
	xlab = "log2(RowMeans)",
	ylab = "log2(CV)",
	main = "Stage 2/2 Filtered Plot of Row Mean v. Row Variance for GV/MII samples",
	col  = c(rep("Green", samp.count/2), rep("Blue", samp.count/2)),
	pch  = c(rep(17, samp.count/2), rep(19, samp.count/2))
	)
legend("topright", c("GV oocyte", "MII oocyte"), pch = c(17, 19), col = c("Green", "Blue"))
abline(v = 3, col = 2, lwd = 2)
dev.off()

# Update gene.info data.frame
gene.info <- subset(gene.info, rownames(gene.info) %in% rownames(dat.filtered))

########################
#  Feature Selection   #
########################   

type <- lapply(colnames(dat.filtered), 
	function(x) {
		if(regexpr("GV", x) < 1) {"MII"} else {"GV"}
	})

# Determine individual P-values in the original expression scale (not log2)
pval <- c()
ptm <- proc.time()
for (i in seq(nrow(dat.filtered))){
	t <- t.test(2^dat.filtered[i, type == "GV"], 2^dat.filtered[i, type == "MII"], alternative = "two.sided")
	pval <- c(pval, t$p.value)
}
proc.time() - ptm

# Output
cat(sprintf("Total number of genes with p-value < 0.05 is %s\n", sum(pval < 0.05))) # 19604
cat(sprintf("Total number of genes with p-value < 0.01 is %s\n", sum(pval < 0.01))) # 15797


# via Lab_05
Bonferroni <- 0.05/length(pval)
cat(sprintf("Bonferroni correction is %s\n", Bonferroni))
cat(sprintf("Total number of genes with p-value < the Bonferroni correction is %s\n", sum(pval < Bonferroni))) # 5953

# Adjust P-values (for sake of comparasion)
# Note: the analyses that follow will utilize the raw p-value
# Method 1
p.holm <- p.adjust(pval, method = "holm")
p.raw  <- sort(pval)
p.holm <- sort(p.holm)
p1     <- as.matrix(p.raw)
p2     <- as.matrix(p.holm)
all.p  <- cbind(p1, p2)
colnames(all.p) <- c("Raw P-value", "Adjusted P-Value")

# # Method 2
# p.holm <- mt.rawp2adjp(pval, proc = c("holmes"))
# p.raw  <- sort(p.holm$adjp[, 1])
# p.holm <- sort(p.holm$adjp[, 2])
# p1     <- as.matrix(p.raw)
# p2     <- as.matrix(p.holm)
# all.p  <- cbind(p1, p2)
# colnames(all.p) <- c("Raw P-value", "Adjusted P-Value")

# Output
cat(sprintf("Total number of genes with p-value < 0.05 is %s\n", sum(p.holm < 0.05))) # 6111
cat(sprintf("Total number of genes with p-value < 0.01 is %s\n", sum(p.holm < 0.01))) # 4878

Bonferroni.2 <- 0.05/length(p.holm)
cat(sprintf("Bonferroni correction is %s\n", Bonferroni.2))
cat(sprintf("Total number of genes with p-value < the Bonferroni correction is %s\n", sum(pval < Bonferroni.2))) # 5953

# Plot the sorted raw P-values
png("Raw_PvaluePlot.png")
plot(
	p.raw, type = "b", pch = 1, col = "Pink",
	xlab = "Genes",
	ylab = "P-values",
	main = "Raw P-values for All Genes\n"
	)
dev.off()
	
# Plot Adjusted vs. Raw P-values
png("PvaluePlot.png")
matplot(
	all.p, type = "b", pch = 1, col = 1:2,
	xlab = "Genes",
	ylab = "P-values",
	main = "Adjusted Vs. Raw P-values for All Genes\n"
	)
legend("bottomright", legend = colnames(all.p), pch = 1, col = 1:2)
dev.off()

# Update gene.info data.frame with P-value information
gene.info$pvalue <- pval

# Hypothesis Pass/Fail Dataframe
thresh  <- 0.05
rnames  <- rownames(dat.filtered)
p.test  <- lapply(as.list(pval), function(f){ if (f < thresh) TRUE else FALSE })
pval.df <- as.data.frame(do.call(rbind, p.test), rname = rnames)
names(pval.df) <- c("Pval.Test")
rownames(pval.df) <- rownames(dat.filtered) 

# Update gene.info data.frame with P-value test results (pass/fail)
gene.info$Pval.Test <- pval.df$Pval.Test

cat(sprintf("Threshold: %s\n", thresh))
table(pval.df$Pval.Test)
#   Number of True        : 19604
#   Number of False       : 15266
#   Total (nrow(pval.df)) : 34870

# Histogram illistrating the distribution of P-values 
png("Histogram_pval.png")
hist(
	pval, 
	col  = "Pink",
	xlab = "P-Value",
	ylab = "Frequency",
	main = "Histogram of T-Test P-Values values for GDS3295 profiles"
	)
hist(
	p.holm, 
	col  = "Red",
	add  = T 
	)
legend("top", legend = colnames(all.p), pch = 16, col = c("Pink", "Red"))
dev.off()

# Determining fold change between treatments in log2
GV.matrix  <- dat.filtered[, type=="GV"]
GV.mean    <- apply(GV.matrix, 1, mean, na.rm = TRUE)
MII.matrix <- dat.filtered[, type=="MII"]
MII.mean   <- apply(MII.matrix, 1, mean, na.rm = TRUE)

fold <- GV.mean - MII.mean
# Note, we are working on a log2() scale
# 2^max(fold) # 180.2804
# 2^min(fold) # 0.0642937

# Update gene.info data.frame with fold information
gene.info$fold <- fold	

# Via Lab_05, determind probesets that demonstrate a 2x fold change
fold.test <- lapply(as.list(fold), function(x){ if (abs(x) > log2(2)) TRUE else FALSE })
fold.df   <- as.data.frame(do.call(rbind, fold.test))
names(fold.df) <- c("Fold.Test")

# Update gene.info data.frame with fold information
gene.info$Fold.Test <- fold.df$Fold.Test

table(fold.df$Fold.Test)
# Number of True        : 9556 
# Number of False       : 25314
# Total (nrow(fold.df)) : 34870

# Histogram illistrating the distribution of log2(fold change)
png("Histogram_fold.png")
hist(
	fold, 
	col  = "Lightblue",
	xlab = "Log2 Fold Change",
	ylab = "Frequency",
	main = paste("Histogram of Fold Change values for GDS3295 profiles")
	)
abline(v = log2(2), col = 2, lwd = 2)
abline(v = -log2(2), col = 2, lwd = 2)
dev.off()

# Overall Volcano plot showing cutoffs and differences within the dataset
p.transformed <- (-1 * log10(pval))
png("VolcanoPlot.png")
plot(
	range(p.transformed),
	range(fold),
	type = "n", xlab = "-1 * log10(P-Value)", ylab = "Fold Change",
	main = "Volcano Plot Illistrating  GV Oocyte and MII Oocyte Differences"
	)
points(
	p.transformed, fold,
	col = 1, bg = 1, pch = 21
	)
points(
	p.transformed[(p.transformed > -log10(0.05) & fold > log2(2))],
	fold[(p.transformed > -log10(0.05) & fold > log2(2))],
	col = 1, bg = 2, pch = 21
	)
points(
	p.transformed[(p.transformed > -log10(0.05) & fold < -log2(2))],
	fold[(p.transformed > -log10(0.05) & fold < -log2(2))],
	col = 1, bg = 3, pch = 21
	)
abline(v = -log10(0.05))
abline(h = -log2(2))
abline(h = log2(2))
dev.off()

# New data.frame with genes that have passed both (Fold and Pval) tests
true.genes <- subset(gene.info, (Fold.Test & Pval.Test) == TRUE)
cat(sprintf("Total number of genes that pass both (Pval and Fold) tests: %s\n", nrow(true.genes))) # 9230

# Write these genes with their corresponding values to an output .txt file
write.table(true.genes, file = "TrueGenes.csv", sep = ",", col.names = NA, qmethod = "double")

# Ordering the highest genes (by P-value) in the form of a data.frame
# Note: dat.filtered is still in log2 scale
best.genes    <- order(pval)[1:length(pval)]
best.genes.df <- data.frame(
					index = best.genes,
					exp   = 2^dat.filtered[best.genes, ],
					pval  = pval[best.genes]
				)

# Expression matrix with the 'best.genes' in the original scale (based on P-value ranking)
top.genes.matrix <- 2^dat.filtered[best.genes, ]

# Feature Selection via svmRFE which utilizes the library e1071
# Majority of code was adapted from Ref1 listed below.
# Ref1.: http://www.colbyimaging.com/wiki/statistics/msvm-rfe
# Ref2.: https://github.com/johncolby/SVM-RFE/blob/master/msvmRFE.R
# Ref3.: http://cran.r-project.org/web/packages/e1071/index.html
t.dat  <- t(top.genes.matrix)
label  <- as.vector(unlist(type))
svm.df <- data.frame(label, t.dat)

# Note: svmRFE will automatically normalize/scale the dataframe
# See msvmRFE.R code in Ref2 above for more information
ptm <- proc.time()
ranked.list <- svmRFE(svm.df, k = 10, halve.above = 1000)
proc.time() - ptm
# Run-time:
#   user  system elapsed 
# 700.65    8.30  755.24 

# Write the rankings to an output .txt file so that it can be read in later if needed
output <- data.frame(RankedOrder = ranked.list)
write.table(output, file = "RankedList.txt")

top.ranked.genes           <- top.genes.matrix[ranked.list, ]
rownames(top.ranked.genes) <- rownames(top.genes.matrix[ranked.list, ])

# Create a new genes.info data.frame for the ranked genes
top.genes.info             <- gene.info[rownames(top.ranked.genes), ]

tg <- top.genes.info$pvalue[top.genes.info$pvalue < thresh]

png("Histogram_ranked_pval.png")
hist(
	tg,
	col  = "Pink",
	xlab = "P-Value",
	ylab = "Frequency",
	main = "Histogram of T-Test P-Values values for ranked profiles"
	)
abline(v = thresh, col = 2, lwd = 2)
dev.off()

top5    <- head(top.genes.info, n = 5L)
bottom5 <- tail(top.genes.info, n = 5L)

###########################
#  Sample Classification  #
###########################  

# PCA analysis on the first 3 component vectors
pca       <- prcomp(t(top.ranked.genes))
pca.loads <- pca$x[, 1:3] 

# As per Lecture 10:
svp <- ksvm(pca.loads,label,type="C-svc")
svp

# Get fitted values
fit <- fitted(svp)

# error rates (incorrect classifications)
er1 <- sum(fit[label == "GV"]== "MII")	# 0 
er2 <- sum(fit[label == "MII"] == "GV") # 0  


# Plot code adapted from Lab_07 and Lab_08
# Plotting: Component 1 vs. Component 2
#			Component 1 vs. Component 3
#			Component 2 vs. Component 3

png("Scatterplot_PCA_1v2.png")
plot(
	range(pca.loads[, 1]), 
	range(pca.loads[, 2]), 
	type = "n",
	xlab = "Principal Component 1",
	ylab = "Principal Component 2",
	main = "PCA Plot for GDS3295 data\n PC1 vs. PC2"
	)
points(
	pca.loads[, 1][as.numeric(type == "GV") == 1], 
	pca.loads[, 2][as.numeric(type == "GV") == 1],
	col = "Red", pch = 15
	)
points(
	pca.loads[, 1][as.numeric(type == "MII") == 1], 
	pca.loads[, 2][as.numeric(type == "MII") == 1],
	col = "Blue", pch = 19
	)
legend(
	"bottomleft", 
	c("MII Oocyte", "GV Oocyte"), 
	col = c("Blue", "Red"), pch = c(19,15)
	)
dev.off()
	
png("Scatterplot_PCA_1v3.png")
plot(
	range(pca.loads[, 1]), 
	range(pca.loads[, 3]), 
	type = "n",
	xlab = "Principal Component 1",
	ylab = "Principal Component 3",
	main = "PCA Plot for GDS3295 data\n PC1 vs. PC3"
	)
points(
	pca.loads[, 1][as.numeric(type == "GV") == 1], 
	pca.loads[, 3][as.numeric(type == "GV") == 1],
	col = "Red", pch = 15
	)
points(
	pca.loads[, 1][as.numeric(type == "MII") == 1], 
	pca.loads[, 3][as.numeric(type == "MII") == 1],
	col = "Blue", pch = 19
	)
legend(
	"bottomright", 
	c("MII Oocyte", "GV Oocyte"), 
	col = c("Blue", "Red"), pch = c(19,15)
	)
dev.off()

png("Scatterplot_PCA_2v3.png")
plot(
	range(pca.loads[, 2]), 
	range(pca.loads[, 3]), 
	type = "n",
	xlab = "Principal Component 2",
	ylab = "Principal Component 3",
	main = "PCA Plot for GDS3295 data\n PC2 vs. PC3"
	)
points(
	pca.loads[, 2][as.numeric(type == "GV") == 1], 
	pca.loads[, 3][as.numeric(type == "GV") == 1],
	col = "Red", pch = 15
	)
points(
	pca.loads[, 2][as.numeric(type == "MII") == 1], 
	pca.loads[, 3][as.numeric(type == "MII") == 1],
	col = "Blue", pch = 19
	)
legend(
	"bottomleft", 
	c("MII Oocyte", "GV Oocyte"), 
	col = c("Blue", "Red"), pch = c(19,15)
	)
dev.off()
	
pca.var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
# pca.var
cat(sprintf("Note: Approximately %s variability is explained using the only the first two eigenvalues.\n", 
	sum(pca.var[1:2])))

# Scree Plot adapted from Lab_07 Illistrating Level of Variance
png("Screeplot.png")
plot(
	c(1:length(pca.var)), 
	pca.var, 
	type = "b", 
	xlab = "Components",
	ylab = "Percent Variance", 
	bg = "Blue", pch = 21
	)
title("Scree Plot Illistrating %-Variability Explained By Each Eigenvalue\n GV/MII Dataset")
dev.off()

# MDS Analysis via Kruskal’s Non-metric Approach
dat.dist <- dist(t(top.ranked.genes))
dat.mds  <- isoMDS(dat.dist)

png("MDSplot_Kruskal.png")
plot(dat.mds$points, type = "n")
points(
	dat.mds$points[, 1][as.numeric(type == "GV") == 1],
	dat.mds$points[, 2][as.numeric(type == "GV") == 1], 
	col = "Red", pch = 16, cex = 1.5
	)
points(
	dat.mds$points[, 1][as.numeric(type == "MII") == 1], 
	dat.mds$points[, 2][as.numeric(type == "MII") == 1], 
	col = "Blue", pch = 16, cex = 1.5)
title(main = "Kruskal’s Non-metric MDS plot for GDS3295 Dataset\nGV vs. MII")
legend("bottomleft", c("GV Oocyte", "MII Oocyte"), col = c("Red", "Blue"), fill = c("Red", "Blue"))
dev.off()

# MDS Analysis via the Classical Metric Approach
# Note: No stress value is provided
dat.loc  <- cmdscale(dat.dist) 

png("MDSplot_Classical.png")
plot(dat.loc, type = "n")
points(
	dat.loc[, 1][as.numeric(type == "GV") == 1],
	dat.loc[, 2][as.numeric(type == "GV") == 1], 
	col = "Red", pch = 16, cex = 1.5
	)
points(
	dat.loc[, 1][as.numeric(type == "MII") == 1], 
	dat.loc[, 2][as.numeric(type == "MII") == 1], 
	col = "Blue", pch = 16, cex = 1.5)
title(main = "Classical Metric MDS plot for GDS3295 Dataset\nGV vs. MII")
legend("bottomleft", c("GV Oocyte", "MII Oocyte"), col = c("Red", "Blue"), fill = c("Red", "Blue"))
dev.off()

# Weighted Laplacian Graph
# Note: The function GeneratePhi was taken from Lab_07 and Lecture_08
GeneratePhi <- function (X, qnt = NULL) {
	dist2full <- function(dis) {
		n     <- attr(dis, "Size")
		full  <- matrix(0, n, n)
		full[lower.tri(full)] <- dis
		full + t(full)
	}
	dat.dis <- dist(t(X),"euc")^2
	if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
	if(is.null(qnt))  {eps <- min(dat.dis[dat.dis != 0])}
	kernal         <- exp(-1 * dat.dis / (eps))
	K1             <- dist2full(kernal)
	diag(K1)       <- 0
	D              <- matrix(0, ncol = ncol(K1), nrow = ncol(K1))
	tmpe           <- apply(K1, 1, sum)
	tmpe[tmpe > 0] <- 1/sqrt(tmpe[tmpe > 0])
	tmpe[tmpe < 0] <- 0
	diag(D)        <- tmpe
	L              <- D%*% K1 %*% D
	X              <- svd(L)$u
	Y              <- X / sqrt(apply(X^2, 1, sum))
}

temp <- t(top.ranked.genes)
temp <- scale(temp, center = T, scale = T)
phi <- GeneratePhi(t(temp), qnt = NULL)

png("LaplacianPlot.png")
plot(
	range(phi[, 1]), range(phi[, 2]),
	xlab = "Phi 1", ylab = "Phi 2",
	main = "Weighted Graph Laplacian Plot for GDS3295 Dataset\nGV vs. MII"
	)
points(
	phi[, 1][as.numeric(type == "GV") == 1],
	phi[, 2][as.numeric(type == "GV") == 1], 
	col = "Red", pch = 16, cex = 1.5
	)
points(
	phi[, 1][as.numeric(type == "MII") == 1], 
	phi[, 2][as.numeric(type == "MII") == 1], 
	col = "Blue", pch = 16, cex = 1.5
	)
legend("top", c("GV Oocyte", "MII Oocyte"), col = c("Red", "Blue"), fill = c("Red", "Blue"))
dev.off()

########################
#   Cluster Analysis   #
########################

# Hierarchical Clustering via Manhattan 
top.dist <- dist(t(top.ranked.genes), method = "manhattan")
top.clus <- hclust(top.dist, method = "median")

png("Dendogram_TopGenes.png")
plot(
	top.clus,
	labels = colnames(top.ranked.genes),  
	xlab   = "Clustered Samples",
	ylab   = "Distance",
	main   = "Hierarchical Clustering Dendrogram\nRanked Oocyte Classification"
	)
dev.off()

# Heatmap of top.ranked.genes
# Note: Due to memory allocation, the top 50 ranked genes were plotted
png("Heatmap_Top50Genes.png")
heatmap(
	top.ranked.genes[1:50, ],
	col  = color,
	xlab = "Samples",
	ylab = "Top Ranked Genes",
	main = "Heatmap for the top 50 genes"
	)
dev.off()

# Heatmap of top.ranked.genes
# Note: Due to memory allocation, 50 randomly ranked genes were plotted 
png("Heatmap_TopRandom50Genes.png")
heatmap(
	top.ranked.genes[sample(top.ranked.genes, 50), ],
	col  = color,
	xlab = "Samples",
	ylab = "Randomly Ranked Genes",
	main = "Heatmap for 50 randomly ranked genes"
	)
dev.off()

# K-means clustering via PCA Analysis 
# Done though the kernlab library
# Extract out the first 5 component vectors and compute K-means with two centers 
dat.kpca <- kpca(t(top.ranked.genes), kernel = "rbfdot", kpar = list(sigma = 0.002), features = 5)
pcv      <- pcv(dat.kpca)
rot      <- rotated(dat.kpca)
pcv.k    <- kmeans(pcv, centers = 2, iter.max = 20)
rot.k    <- kmeans(rot, centers = 2, iter.max = 20)

# 2D scatterplot of the first 5 eigenfeatures (from PCA)
png("Scatterplot_PCA.png")
par(mfrow = c(2, 1))
plot(
	pcv, col = pcv.k$cluster, cex = 1,
	main = "PCA Scatter Plot with PC vectors (column-wise)\nk = 2",
	xlab = "P1",
	ylab = "P2"
	)
points(pcv.k$centers, col = 1:2, pch = "*", cex = 2.5)

plot(
	rot, col = rot.k$cluster, cex = 1,
	main = "PCA Scatter Plot with PC Projected vectors\nk = 2",
	xlab = "P1",
	ylab = "P2"
	)
points(rot.k$centers, col = 1:2, pch = "*", cex = 2.5)
dev.off()

# LDA Analysis
# Note: Code is adapted from Lab_09

training <- as.data.frame(rbind(t.dat[c(1, 2, 5, 6, 9, 10, 13, 14), ]))
test     <- as.data.frame(t.dat[ !(rownames(t.dat) %in% rownames(training)), ])

te.names <- rownames(test)
te.names[c(1, 2, 5, 6)] <- paste("66", te.names[c(1, 2, 5, 6)], sep = "_")
te.names[c(3, 4, 7, 8)] <- paste("6", te.names[c(3, 4, 7, 8)], sep = "_")
test.names <- factor(gsub('_GSM[[:digit:]]+', '', te.names))

tr.names <- rownames(training)
tr.names[c(1, 2, 5, 6)] <- paste("66", tr.names[c(1, 2, 5, 6)], sep = "_")
tr.names[c(3, 4, 7, 8)] <- paste("6", tr.names[c(3, 4, 7, 8)], sep = "_")
train.names <- factor(gsub('_GSM[[:digit:]]+', '', tr.names))

# Due to memory allocations LDA was run on the first 5000 genes
ptm <- proc.time()
train.lda.2      <- lda(train.names ~ ., data = training[, c(1:5000)])
train.pred.2.out <- predict(train.lda.2, test[, c(1:5000)])
proc.time() - ptm
table(train.pred.2.out$class, test.names)

png("LDAPlot_2D.png")
plot(
	range(train.pred.2.out$x[, 1]),
	range(train.pred.2.out$x[, 2]),
	type = "n",
	xlab = "LD1",
	ylab = "LD2",
	main = "LDA Plot of 5000 Genes\nDisciminant Scores on 2 Discriminant Variables",
	)
points(
	train.pred.2.out$x[, 1],
	train.pred.2.out$x[, 2],
	col = c(rep("Green", 2), rep("Blue", 2), rep("Black", 2), rep("Red", 2)),
	pch = c(rep(16, 4), rep(17, 4))
	)
legend(
	"bottom",
	c("GV - 66 Weeks", "GV - 6 Weeks", "MII - 66 Weeks", "MII6 - 6 Weeks"),
	col = c("Green", "Blue", "Black", "Red"),
	pch = c(16, 17)
	)
dev.off()

col = c(rep("Green", 2), rep("Blue", 2), rep("Black", 2), rep("Red", 2))
png("LDAPlot_3D.png")
scatterplot3d(
	train.pred.2.out$x[, 1],
	train.pred.2.out$x[, 2],
	train.pred.2.out$x[, 3],
	xlab = "LD1",
	ylab = "LD2",
	zlab = "LD3",
	col,
	pch = c(rep(16, 4), rep(17, 4)),
	col.axis = "blue",
	col.grid = "lightblue", 
	main = "LDA Plot of 5000 Genes\nDisciminant Scores on 3 Discriminant Variables"
	)
dev.off()