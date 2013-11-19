# AS.410.671.81.FA13 – HW #3
# Author: Shan Sabri

#######################
#          1          #
####################### 

library(Biobase)
library(multtest)
library(annotate)
data(golub)


#######################
#          2          #
####################### 

dat                <- as.data.frame(golub)
dimnames(dat)[[1]] <- paste("g", dimnames(dat)[[1]], sep = "")


#######################
#          3          #
#######################

ann                <- golub.cl
dimnames(dat)[[2]] <- ann 


#######################
#          4          #
#######################

wilcox.test.all.genes <- function(x,s1,s2) {
	x1 <- x[s1]
	x2 <- x[s2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	w.out <- wilcox.test(x1, x2, exact = FALSE, alternative = "two.sided", correct = TRUE, var.equal = T)
	out <- as.numeric(w.out$statistic)
	return(out)
}

original.wmw.run <- apply(dat, 1, wilcox.test.all.genes, s1 = colnames(dat) == 0, s2 = colnames(dat) == 1)


#######################
#          5          #
#######################

wilcox.iterate <- function(x){
	list <- list(1:500)
	for (i in 1:500){
		colnames(x) <- sample(colnames(x))
		stat        <- apply(x, 1, wilcox.test.all.genes, s1 = colnames(x) == 0, s2 = colnames(x) == 1)
		list[i]     <- max(stat)
	}
	result <- list
	return(result)
}

ptm      <- proc.time()
stat.max <- wilcox.iterate(dat)
proc.time() - ptm
stat.max <- t(stat.max)
stat.max <- as.character(stat.max)


#######################
#          6          #
#######################

cutoff <- sort(stat.max)[0.95 * length(stat.max)]
subset <- original.wmw.run[original.wmw.run > cutoff]
subset


#######################
#          7          #
#######################

library(limma)

zero <- dat[colnames(dat) == 0]
one  <- dat[colnames(dat) == 1]
x    <- cbind(a = 1, b = c(rep(1, length(zero)), rep(0, length(one))))
fit  <- lmFit(dat, x)
fit  <- eBayes(fit)

topTable(fit)
attributes(fit)

pval <- fit$p.value[ , 2]


#######################
#          8          #
#######################

n     <- length(subset)
sort  <- sort(pval)
sort  <- sort[1:n]
inter <- intersect(names(sort), names(subset))
length(inter)
inter


#######################
#          9          #
#######################

t.test.all.genes <- function(x, d1, d2){
	x1 <- x[d1]
	x2 <- x[d2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	t.out <- t.test(x1, x2, alternative="two.sided", var.equal=T)
	out <- as.numeric(t.out$p.value)
	return(out)
}

t.run   <- apply(dat, 1, t.test.all.genes, d1 = colnames(dat) == 0, d2 = colnames(dat) == 1)
t.run.p <- t.run[t.run < 0.01]
t.run.p <- as.matrix(t.run.p)

bayes.names <- as.matrix(fit$p.value)
bayes.names <- bayes.names[ , -2]
bayes.names <- as.matrix(bayes.names)

overlap           <- merge(t.run.p, bayes.names, by = "row.names", all = FALSE)
overlap           <- as.matrix(overlap)
rownames(overlap) <- overlap[ , 1]
overlap           <- overlap[ , -1]
colnames(overlap) <- c("Student T-Test", "Empirical Bayes")
class(overlap)    <- "numeric"

plot(c(1, nrow(overlap)), range(overlap), type = "n", xaxt ="n", ylab = "P-Value", xlab="Genes")
points(1:nrow(overlap), col = "Green", overlap[ , 2])
points(1:nrow(overlap), col = "Blue",  overlap[ , 1])
title(main="P-Value Plot\nStudent T-Test Vs. Empirical Bayes")
legend(350, 1 , colnames(overlap), col = c("Green", "Blue"), pch = 15, cex = 0.9)