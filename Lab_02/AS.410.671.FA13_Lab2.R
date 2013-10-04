# AS.410.671.81.FA13 – Lab #2
# Data Visualization 
# Author: Shan Sabri

#######################
#          2          #
####################### 

dat <- read.table(
		file = "C:\\Users\\Shan\\Desktop\\JHU_Fall2013\\410.671_Microarrays&Analysis\\DataSets\\spellman.txt", 
		header = T, 
		row.names = 1
		)
dat <- as.data.frame(dat);

		

#######################
#          3          #
#######################
 
dim(dat)
# dimnames(dat)[[1]] #Row names
# dimnames(dat)[[2]] #Col names


#######################
#          4          #
####################### 

cdc15 <- dat[, 23:46]
# dimnames(cdc15)[[2]] #Only CDC15 mutants 

#######################
#          5          #
####################### 

library(gplots)
dat.cor <- cor(cdc15, use = "pairwise.complete.obs", method = "pearson")

layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2), 5, 2, byrow = TRUE))
par(oma = c(5, 7, 1, 1))
cx <- rev(colorpanel(25, "yellow", "black", "blue"))
leg <- seq(min(dat.cor, na.rm = T), max(dat.cor, na.rm = T), length = 10)
image(dat.cor, main = "Correlation plot of CDC15 mutant between time points", axes = F, col = cx)
axis(1, at = seq(0, 1, length = ncol(dat.cor)), label = dimnames(dat.cor)[[2]], cex.axis = 0.9, las = 2)
axis(2, at = seq(0, 1, length = ncol(dat.cor)), label = dimnames(dat.cor)[[2]], cex.axis = 0.9, las = 2)

image(as.matrix(leg), col = cx, axes = F)
tmp <- round(leg, 2)
axis(1, at = seq(0, 1, length = length(leg)), labels = tmp, cex.axis = 1)


#######################
#          6          #
####################### 

VPS8      <- cdc15["YAL002W", ] 
VPS8.num  <- as.numeric(VPS8)
VPS8.mean <- mean(VPS8.num, na.rm = TRUE)

for (i in colnames(VPS8)) { 
	for (j in rownames(VPS8)) {  
		VPS8[j,i] <- replace(VPS8[j,i], is.na(VPS8[j,i]), VPS8.mean)
	} 
} 

VPS8.mean
VPS8


#######################
#          7          #
####################### 

VPS8 <- t(scale(t(VPS8)))

colnames(VPS8) <- c("10" , "30" , "50" , "70" , "80" , "90" , 
   		    "100", "110", "120", "130", "140", "150",
		    "160", "170", "180", "190", "200", "210", 
		    "220", "230", "240", "250", "270", "290")
					
plot(c(1, ncol(VPS8)), range(VPS8), type = 'n', main = "Profile plot of VPS8", xlab = "Timepoints", ylab = "Expression", axes = F)
axis(side = 1, at = c(1:ncol(VPS8)), labels = dimnames(VPS8)[[2]], cex.axis = 0.4, las = 2)
axis(side = 2)
for(i in 1:length(VPS8)) {
	lines(c(1:ncol(VPS8)), VPS8, col = i, lwd = 2)
}


