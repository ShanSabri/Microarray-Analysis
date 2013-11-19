
# AS.410.671.81.FA13 – Lab #6
# Multiple Testing
# Author: Shan Sabri

#######################
#          2          #
####################### 

setwd("C:\\Users\\Shan\\Desktop\\JHU_Fall2013\\410.671_Microarrays&Analysis\\DataSets\\")
dat <- read.table("agingStudy11FCortexAffy.txt", header = T, row.names = 1)
ann <- read.table("agingStudy1FCortexAffyAnn.txt", header = T, row.names = 1)


#######################
#          3          #
####################### 

dat.gender <- dat
age.sort   <- ann[sort.list(ann[, 2]), ]
age        <- paste(dimnames(age.sort)[[1]], age.sort[, 2], age.sort[, 1], sep = ".") 	
dat.age    <- dat[, age] 


#######################
#          4          #
####################### 

# Gender comparison gene vector
g.g <- c(1394,  1474,  1917,  2099,  2367,  2428, 2625,  3168,  3181,  3641,  3832,  4526,
		 4731,  4863,  6062,  6356,  6684,  6787,  6900,  7223,  7244,  7299,  8086,  8652,
		 8959,  9073,  9145,  9389, 10219, 11238, 11669, 11674, 11793)

# Age comparison gene vector
g.a <- c(25, 302,  1847,  2324,  246,  2757, 3222, 3675,  4429,  4430,  4912,  5640, 5835,
		 5856,  6803,  7229,  7833,  8133, 8579,  8822,  8994, 10101, 11433, 12039, 12353,
		 12404, 12442, 67, 88, 100)

t.test.all.genes <- function(x,s1,s2) {
	x1 <- x[s1]
	x2 <- x[s2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	t.out <- t.test(x1, x2, alternative = "two.sided", var.equal = TRUE)
	out <- as.numeric(t.out$p.value)
	return(out)
}

rawp.gender <- apply(dat.gender[g.g, ], 1, t.test.all.genes, s1 = ann[, 1] == "M", s2 = ann[, 1] == "F") 
rawp.age    <- apply(dat.age[g.a, ], 1, t.test.all.genes, s1 = age.sort[, 2] < 50, s2 = age.sort[, 2] >= 50) 

# Method 1
p.holm.gender <- p.adjust(rawp.gender, method = "holm") 
p.holm.age    <- p.adjust(rawp.age, method = "holm") 

# Method 2 
library(multtest)
holm.gender <- mt.rawp2adjp(rawp.gender, proc = c("Holmes"))
holm.age    <- mt.rawp2adjp(rawp.age, proc = c("Holmes"))


#######################
#          5          #
####################### 

rawp.gender.sort   <- sort(rawp.gender)
p.holm.gender.sort <- sort(p.holm.gender) 
gen1 <- as.matrix(rawp.gender.sort)
gen2 <- as.matrix(p.holm.gender.sort)
gen  <- cbind(gen1, gen2) 
colnames(gen) <- c("Raw P-value", "Adjusted P-Value") 
matplot(gen, type = "b", pch = 1, col = 1:2, 
		main = "Gender P-value Plot\n(p.adjust holm method)", 
		ylab = "P-values") 
legend("topleft", legend = colnames(gen), pch = 1, col = 1:2) 

rawp.age.sort      <- sort(rawp.age)
p.holm.age.sort    <- sort(p.holm.age) 
age1 <- as.matrix(rawp.age.sort)
age2 <- as.matrix(p.holm.age.sort)
age  <- cbind(age1, age2) 
colnames(age) <- c("Raw P-value", "Adjusted P-Value") 
matplot(age, type = "b", pch = 1, col = 1:2, 
		main = "Age P-value Plot\n(p.adjust holm method)", 
		ylab = "P-values") 
legend("topleft", legend = colnames(gen), pch = 1, col = 1:2) 


#######################
#          6          #
####################### 

bon.gender <- p.adjust(rawp.gender, method = "bonferroni") 
bon.gender.sort <- sort(bon.gender) 
gen1.b <- as.matrix(rawp.gender.sort)
gen2.b <- as.matrix(bon.gender.sort)
gen.b  <- cbind(gen1.b, gen2.b) 
colnames(gen.b) <- c("Raw P-value", "Adjusted P-Value") 
matplot(gen.b, type = "b", pch = 1, col = 1:2, 
		main = "Gender P-value Plot\n(p.adjust bonferroni method)", 
		ylab = "P-values") 
legend("topleft", legend = colnames(gen.b), pch = 1, col = 1:2) 

bon.age    <- p.adjust(rawp.age, method = "bonferroni") 
bon.age.sort <- sort(bon.age) 
age1.b <- as.matrix(rawp.age.sort)
age2.b <- as.matrix(bon.age.sort)
age.b  <- cbind(age1.b, age2.b) 
colnames(age.b) <- c("Raw P-value", "Adjusted P-Value") 
matplot(age.b, type = "b", pch = 1, col = 1:2, 
		main = "Age P-value Plot\n(p.adjust bonferroni method)", 
		ylab = "P-values") 
legend("topleft", legend = colnames(age.b), pch = 1, col = 1:2) 


# Additional plots
res1     <- mt.maxT(dat.gender[g.g, ], ann[, 1])
p        <- res1$rawp[order(res1$index)]
teststat <- res1$teststat[order(res1$index)]
procs    <- c("Bonferroni", "Holm")
res2     <- mt.rawp2adjp(p, procs)
allp     <- cbind(res2$adjp[order(res2$index),],res1$adjp[order(res1$index)])
ltypes   <- c(3,rep(1,3))
cols     <- c(1:4)
mt.plot(allp, teststat, plottype = "pvsr", proc = procs, leg = c(80,0.4), 
		main = "GENDER", lty = ltypes, col = cols, lwd = 2)

res1.b     <- mt.maxT(dat.age[g.a, ], ann[, 1])
p.b        <- res1.b$rawp[order(res1.b$index)]
teststat.b <- res1.b$teststat[order(res1.b$index)]
procs.b    <- c("Bonferroni", "Holm")
res2.b     <- mt.rawp2adjp(p.b, procs.b)
allp.b     <- cbind(res2.b$adjp[order(res2.b$index),],res1.b$adjp[order(res1.b$index)])
mt.plot(allp.b, teststat.b, plottype = "pvsr", proc = procs, leg = c(80,0.4), 
		main = "AGE", lty = ltypes, col = cols, lwd = 2)



