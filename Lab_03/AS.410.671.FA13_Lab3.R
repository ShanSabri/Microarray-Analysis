# AS.410.671.81.FA13 – Lab #3
# Power and Sample Size
# Author: Shan Sabri

#######################
#          2          #
####################### 

dat <- read.table(
		file = "C:\\Users\\Shan\\Desktop\\JHU_Fall2013\\410.671_Microarrays&Analysis\\DataSets\\eisen.txt", 
		header = T, 
		na.strings="NA",
		blank.lines.skip=F,
		row.names = 1
		)
		
dat <- as.data.frame(dat)


#######################
#          3          #
####################### 

classLabel <- read.table(
		file = "C:\\Users\\Shan\\Desktop\\JHU_Fall2013\\410.671_Microarrays&Analysis\\DataSets\\eisenClasses.txt", 
		header = T
		)


#######################
#          4          #
####################### 

dimnames(dat)[[2]] 

cl  <- as.character(classLabel[, 2])
dat <- dat[, cl]

germinalCentre <- cl[1:19]
activated      <- cl[20:39]

dimnames(dat)[[2]] 


#######################
#          5          #
####################### 

x <- as.numeric(dat[8008, germinalCentre])
y <- as.numeric(dat[8008, activated])

x <- x[!is.na(x)]	
y <- y[!is.na(y)]

xy.list <- list(x, y) 

boxplot(
	xy.list, 
	col = c("red", "blue"), 
	main = "Gene 8008 from the DLBCL cDNA 2-channel dataset", 
	axes = F, 
	ylab = "log2(ratio intensity)"
	)
axis(2)
axis(1, at=c(1, 2), c("GerminalCentre", "Activated"))


par(mfrow=c(2, 1))
hist(
	x,
	col = "red",
	labels = T,
	main = "Gene 8008 from the DLBCL cDNA 2-channel dataset (GerminalCentre)", 
	axes = T, 
	ylab = "Frequency",
	xlab = "log2(ratio intensity) ranges"
	)
hist(
	y,
	col = "blue",
	labels = T,
	main = "Gene 8008 from the DLBCL cDNA 2-channel dataset (Activated)", 
	axes = T, 
	ylab = "Frequency",
	xlab = "log2(ratio intensity) ranges"
	)

	
#######################
#          6          #
####################### 

nx <- length(x)
ny <- length(y)

pool.var <- (((nx - 1)*var(x)) + ((ny - 1)*var(y)))/(nx + ny - 2)
pool.var  # 1.006196

dif.fold <- log2(1.5)/sqrt(pool.var)
dif.fold  # 0.5831587

pl.ss3 <- power.t.test(
		d = dif.fold,
		sig.level = .01, 
		power = 0.8, 
		type = "two.sample"
		)
pl.ss3


#######################
#          7          #
####################### 

dif <- abs(mean(x)-mean(y))/sqrt(pool.var)
dif  # 0.5769141

pl.ss3 <- power.t.test(
		d = dif,
		sig.level = .01, 
		power = 0.8, 
		type = "two.sample"
		)
pl.ss3


#######################
#          8          #
####################### 

library(ssize)
library(gdata) 

sd = apply(dat, 1, sd, na.rm = TRUE)

hist(
	sd,
	col = "Blue",
	labels = T,
	main = "Histogram of S.D. of all genes within the DLBCL dataset", 
	axes = T, 
	ylab = "Frequency",
	xlab = "Standard Deviation Ranges"
	)

	
#######################
#          9          #
####################### 

fold.change = 3.0
sig.level = 0.05 
power = 0.8

all.size <- ssize(
 	 sd = sd, 
	 delta = log2(fold.change), 
	 sig.level = sig.level, 
	 power = power
	 ) 	 
ssize.plot(all.size, lwd = 2, col = "lime green", xlim = c(1,20)) 
xmax <- par("usr")[2] - 1 
ymin <- par("usr")[3] + 0.05 
title("Sample Size to Detect 3-Fold Change")  
legend(
	x = xmax, 
	y = ymin, 
	legend = strsplit(paste(
		"fold change = ", fold.change,", ",
		"alpha = ", sig.level, ", ", 
		"power = ",power,", ", 
		"# genes = ", length(sd), sep=''), ", ")[[1]], 
	xjust = 1, 
	yjust = 0, 
	cex = 1.0
	) 
	 

