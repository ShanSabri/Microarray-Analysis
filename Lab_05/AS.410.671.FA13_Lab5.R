# AS.410.671.81.FA13 – Lab #5
# Differential Expression 
# Author: Shan Sabri

#######################
#          2          #
####################### 

setwd("C:\\Users\\Shan\\Desktop\\JHU_Fall2013\\410.671_Microarrays&Analysis\\DataSets\\")
kd <- read.table("rat_KD.txt", header = T, row.names = 1)


#######################
#          3          #
####################### 

kd.log <- log2(kd)
cl     <- as.character(names(kd.log))
kd.log <- kd.log[ ,cl]
cd     <- cl[1:6]
kd     <- cl[7:11]

t.test.all.genes <- function(x, s1, s2){
	x1 <- x[s1]
	x2 <- x[s2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	t.out <- t.test(x1, x2, alternative = "two.sided", var.equal = T)
	out <- as.numeric(t.out$p.value)
	return(out)
}

pv <- apply(kd.log, 1, t.test.all.genes, s1 = cd, s2 = kd)


#######################
#          4          #
####################### 

length(pv[pv < 0.05])
# 5160 probesets have p < 0.05

length(pv[pv < 0.01])
# 2414 probesets have p < 0.01

x <- 0.05/length(pv)
length(pv[pv < x])
# 12 probesets have p < 0.05/15923

hist(
	pv,
	col      = "lightblue",
	xlab     = "p-values",
	main     = "P-value distance between\n Ketogenic and Control diets",
	cex.main = 0.9
	)
abline(v = 0.05, col = 2, lwd = 2)


#######################
#          5          #
####################### 

cd.m <- apply(dat.log[ ,cd], 1, mean,na.rm = T)
kd.m <- apply(dat.log[ ,kd], 1, mean,na.rm = T)

fold <- cd.m - kd.m


#######################
#          6          #
####################### 

2^max(fold)
# 55.15521
2^min(fold)
# 0.08240443

fold.bon <- 2^fold
names(pv[pv < x & abs(fold.bon) > 2])


#######################
#          7          #
####################### 

# Probset	    | Gene Title via NetAffy
# --------------------------------------------------
# 1367553_x_at	| Hemoglobin beta
# 1370239_at 	| Hemoglobin alpha, adult chain 2
# 1370240_x_at 	| Hemoglobin alpha, adult chain 2
# 1371102_x_at 	| Hemoglobin beta globin minor gene
# 1371245_a_at	| Beta globin minor gene
# 1388608_x_at	| Hemoglobin alpha, adult chain 2

# [General association: Hemoglobin] 


#######################
#          8          #
####################### 

p.transformed <- (-1 * log10(pv))

plot(
	range(p.transformed),
	range(fold),
	type = "n", xlab = "-1 * log10 (P-Value)", ylab = "Fold Change", 
	main = "Volcano Plot of Ketogenic and Control Diets Differences"
	)

points(
	p.transformed, fold, 
	col = 1, bg = 1, pch = 21 
	)
points(
	p.transformed[(p.transformed > -log10(0.05) & fold >  log2(2))], fold[(p.transformed > -log10(0.05) & fold >  log2(2))],
	col = 1, bg = 2, pch = 21
	)
points(
	p.transformed[(p.transformed > -log10(0.05) & fold < -log2(2))], fold[(p.transformed > -log10(0.05) & fold < -log2(2))], 
	col = 1, bg = 3, pch = 21
	)

abline(v = -log10(0.05))
abline(h = -log2(2))
abline(h =  log2(2))









