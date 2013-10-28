# AS.410.671.81.FA13 – HW #2
# Author: Shan Sabri

#######################
#          1          #
####################### 

library(limma)
library(marray)

dir.path <- "C:\\Users\\Shan\\Desktop\\HW_02\\" 
dat <- read.GenePix(
	path = dir.path, 
	name.Gf = "F532 Median", name.Gb = "B532 Median", 
	name.Rf = "F635 Median", name.Rb = "B635 Median", 
	name.W ="Flags" 
	) 
	

#######################
#          2          #
####################### 

# Subject 1

subject.one        <- dat[ , 1] 
subject.one.Median <- maNorm(subject.one, norm = c("median"))
subject.one.Loess  <- maNorm(subject.one, norm = c("loess"))
subject.one.PTL    <- maNorm(subject.one, norm = c("printTipLoess")) 
par(mfrow = c(4, 1)) 
maPlot(subject.one, main = "Subject 1 - No Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.one.Median, main = "Subject 1 - Global Median Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.one.Loess, main = "Subject 1 - Loess Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.one.PTL, main = "Subject 1 - Print-tip-group Loess Normalization", lines.func = NULL, legend.func = NULL) 
 
# Subject 2

subject.two        <- dat[ , 2] 
subject.two.Median <- maNorm(subject.two, norm = c("median"))
subject.two.Loess  <- maNorm(subject.two, norm = c("loess"))
subject.two.PTL    <- maNorm(subject.two, norm = c("printTipLoess")) 
par(mfrow = c(4, 1)) 
maPlot(subject.two, main = "Subject 2 - No Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.two.Median, main = "Subject 2 - Global Median Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.two.Loess, main = "Subject 2 - Loess Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.two.PTL, main = "Subject 2 - Print-tip-group Loess Normalization", lines.func = NULL, legend.func = NULL) 
 
# Subject 3

subject.three        <- dat[ , 3] 
subject.three.Median <- maNorm(subject.three, norm = c("median"))
subject.three.Loess  <- maNorm(subject.three, norm = c("loess"))
subject.three.PTL    <- maNorm(subject.three, norm = c("printTipLoess")) 
par(mfrow = c(4, 1)) 
maPlot(subject.three, main = "Subject 3 - No Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.three.Median, main = "Subject 3 - Global Median Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.three.Loess, main = "Subject 3 - Loess Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.three.PTL, main = "Subject 3 - Print-tip-group Loess Normalization", lines.func = NULL, legend.func = NULL) 
 
# Subject 4

subject.four        <- dat[ , 4] 
subject.four.Median <- maNorm(subject.four, norm = c("median"))
subject.four.Loess  <- maNorm(subject.four, norm = c("loess"))
subject.four.PTL    <- maNorm(subject.four, norm = c("printTipLoess")) 
par(mfrow = c(4, 1)) 
maPlot(subject.four, main = "Subject 4 - No Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.four.Median, main = "Subject 4 - Global Median Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.four.Loess, main = "Subject 4 - Loess Normalization", lines.func = NULL, legend.func = NULL) 
maPlot(subject.four.PTL, main = "Subject 4 - Print-tip-group Loess Normalization", lines.func = NULL, legend.func = NULL) 
 

#######################
#          3          #
####################### 

plot(
	density(na.omit(maM(subject.four))), 
	main = "Density Plot for Array 4 (pre- and post-normalized)", 
	ylim=c(0, 0.9),xlim=c(-6.6, 11.7),
	col="black"
	)
lines(density(na.omit(maM(subject.four.Median))), col = "green") 
lines(density(na.omit(maM(subject.four.Loess))), col = "blue") 
lines(density(na.omit(maM(subject.four.PTL))), col = "red") 

leg.txt <- c(
	"No normalization", 
	"Global Median Normalization", 
	"Loess Normalization", 
	"Print-tip-group Loess Normalization"
	)
legend(1, 0.85, legend = leg.txt, lty = c(1, 1), lwd = c(2.5, 2.5),col = c("black", "green", "blue", "red")) 


#######################
#          5          #
####################### 

for(i in 1:4){
	name <- paste("sample", i, sep = ".")
	bg <- maRb(dat[ , i])
	fg <- maRf(dat[ , i])
	diff <- fg - bg
	assign(name, log2(diff))
} 

data.prenorm <- cbind(sample.1, sample.2, sample.3, sample.4)
data.median  <- apply(data.prenorm, 2, median, na.rm = T)
data.norm    <- sweep(data.prenorm, 2, data.median)

colnames(data.norm) <- c("Array 1", "Array 2", "Array 3", "Array 4")

median(data.norm[ , 1], na.rm = T) 
median(data.norm[ , 2], na.rm = T)
median(data.norm[ , 3], na.rm = T)
median(data.norm[ , 4], na.rm = T)


#######################
#          6          #
####################### 		
				
s1 <- maM(subject.one.Loess)
s2 <- maM(subject.two.Loess)
s3 <- maM(subject.three.Loess)
s4 <- maM(subject.four.Loess)

names <- c("Array 1", "Array 2", "Array 3", "Array 4")

data.loess               <- cbind(s1, s2, s3, s4)
colnames(data.loess)     <- make.names(names, unique = TRUE)
data.loess.cor           <- cor(data.loess, use = "complete.obs", method = "spearman") 
colnames(data.loess.cor) <- make.names(names, unique = TRUE)
rownames(data.loess.cor) <- make.names(names, unique = TRUE)
data.loess.cor

data.cor             <- cor(data.norm, use = "complete.obs", method = "spearman") 
rownames(data.cor)   <- make.names(names, unique = TRUE)
colnames(data.cor)   <- make.names(names, unique = TRUE)
data.cor

pairs(data.norm, main = "Scatterplot Matrix\nGlobal Median Normalization")
pairs(data.loess, main = "Scatterplot Matrix\nLoess Normalization")


#######################
#          7          #
####################### 

for(i in 1:4){
	name <- paste("samp", i, sep = "")
	bg <- maRb(dat[, i])
	fg <- maRf(dat[, i])
	assign(name, fg - bg)
} 

samp           <- cbind(samp1, samp2, samp3, samp4)
colnames(samp) <- c("samp1", "samp2", "samp3", "samp4")
samp.sorted    <- apply(samp, 2, sort) 
row.means      <- rowMeans(samp.sorted, na.rm = FALSE)

order1 <- order(samp1)  
order2 <- order(samp2)
order3 <- order(samp3)
order4 <- order(samp4)

norm1 <- rep(NA, nrow(samp))
norm2 <- rep(NA, nrow(samp))
norm3 <- rep(NA, nrow(samp))
norm4 <- rep(NA, nrow(samp))

norm1[order1] <- row.means
norm2[order2] <- row.means
norm3[order3] <- row.means
norm4[order4] <- row.means

quant.norm <- cbind(norm1, norm2, norm3, norm4)
head(samp)
head(quant.norm)
# median(quant.norm[ , 1])
# median(quant.norm[ , 2])
# median(quant.norm[ , 3])
# median(quant.norm[ , 4])

######################################################
# Verifying results quantitatively  
######################################################

	# library(preprocessCore)
	# head(normalize.quantiles(samp, copy = TRUE))

	# head(y <- normalizeBetweenArrays(samp))
	# head(y)
	# median(y[ , 1]); median(y[ , 2]); median(y[ , 3]); median(y[ , 4]);
		
######################################################
# Verifying results visually
######################################################

	# par(mfrow = c(4, 1))
	# hist(log2(quant.norm[ , 1]))
	# hist(log2(quant.norm[ , 2])) 
	# hist(log2(quant.norm[ , 3]))
	# hist(log2(quant.norm[ , 4]))
			
	# boxplot(log2(quant.norm), xlab="Sample", ylab="Log2 Signal", axes=F)
	# axis(1, labels = 1:ncol(quant.norm), at = 1:ncol(quant.norm))
	# axis(2)


#######################
#          8          #
####################### 

quant.norm.log <- log2(quant.norm)	
quant.norm.cor <- cor(quant.norm.log, use = "complete.obs", method = "spearman") 
quant.norm.cor

pairs(quant.norm.log, main = "Scatterplot Matrix\nQuantile Normalization")
