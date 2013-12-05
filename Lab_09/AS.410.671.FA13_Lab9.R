# AS.410.671.81.FA13 – Lab #9
# Classification
# Author: Shan Sabri

#######################
#          1          #
#######################

setwd("C:\\Users\\Shan\\Desktop\\JHU_Fall2013\\410.671_Microarrays&Analysis\\DataSets\\")
dat <- read.table("lung_cancer.txt", header = T, row.names = 1)


#######################
#          2          #
#######################

library(MASS)

names <- colnames(dat)
length(names)

dat  <- data.frame(names, t(dat))
dim(dat)


#######################
#          3          #
#######################

training <- rbind(dat[1:6, ], dat[11:16, ], dat[20:22, ])
test     <- dat[ !(rownames(dat) %in% rownames(training)), ]

te.names    <- test$names
test.names  <- factor(gsub('[[:digit:]]+', '', te.names))
test$names  <- NULL


#######################
#          4          #
#######################

tr.names       <- training$names
train.names    <- factor(gsub('[[:digit:]]+', '', tr.names))
training$names <- NULL

train.lda.2      <- lda(train.names ~ ., data = training[, c(1, 2)])
train.pred.2.out <- predict(train.lda.2, test[, c(1, 2)])
table(train.pred.2.out$class, test.names)


#######################
#          5          #
#######################

plot(
	range(train.pred.2.out$x[, 1]), 
	range(train.pred.2.out$x[, 2]), 
	type = "n", 
	xlab = "LD1",
	ylab = "LD2",
	main = "LDA Plot of 2 Genes from Training Set",
	)
points(
	train.pred.2.out$x,
	col = c(rep("Green", 4), rep("Blue", 3), rep("Red", 2)),
	pch = c(rep(16, 4), rep(17, 3), rep(18, 2))
	)
legend(
	"bottomright", 
	c("Adeno", "SCLC", "Normal"), 
	col = c("Green", "Blue", "Red"),
	pch = c(16:18)
	)

	
#######################
#          6          #
#######################

train.lda <- lda(train.names ~ ., data = training)
train.out <- predict(train.lda, test)
table(train.out$class, test.names)


#######################
#          7          #
#######################

plot(
	range(train.out$x[, 1]), 
	range(train.out$x[, 2]), 
	type = "n", 
	xlab = "LD1",
	ylab = "LD2",
	main = "LDA Plot of All Genes from Training Set",
	)
points(
	train.out$x,
	col = c(rep("Green", 4), rep("Blue", 3), rep("Red", 2)),
	pch = c(rep(16, 4), rep(17, 3), rep(18, 2))
	)
legend(
	"bottomright", 
	c("Adeno", "SCLC", "Normal"), 
	col = c("Green", "Blue", "Red"),
	pch = c(16:18)
	)