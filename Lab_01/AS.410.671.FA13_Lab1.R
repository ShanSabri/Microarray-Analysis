# AS.410.671.81.FA13 – Lab #1
# Basic R syntax/plots with data solutions
# Author: Shan Sabri

#######################
#          2          #
####################### 

dat <- read.table(
		file="C:\\Users\\Shan\\Desktop\\JHU_Fall2013\\410.671_Microarrays&Analysis\\DataSets\\sle_b_cell.txt", 
		header=T, 
		row.names=1
		)
		

#######################
#          3          #
#######################
 
dim(dat)


#######################
#          4          #
####################### 

dimnames(dat)[[2]]


#######################
#          5          #
####################### 

x <- dat[, 18]
y <- dat[, 2]

plot(
	x, y, 
	xlab = "Normal", ylab = "SLE", 
	main = "SLE B cell sample vs. Normal B cell sample – all probesets"
	)
grid(
	col = "grey", lty = "dotted"
	) 


#######################
#          6          #
####################### 

x <- dat[1:20, 18]
y <- dat[1:20, 2]

plot(
	x, y, pch = 15, col = "blue", 
	xlab = "Normal", ylab = "SLE", 
	main = "SLE B cell sample vs. Normal B cell sample – 20 probesets"
	)
grid(
	col = "grey", lty = "dotted"
	)  


#######################
#          7          #
####################### 

x <- range(1:26)
y <- range(dat["211881_x_at",])

x.line <- c(1:26) 
y.line <- dat["211881_x_at",]
# mode(y.line)
y.line.num <- as.numeric(y.line)
# mode(y.line.num)

plot(
	x, y, type = "n",
	xlab = "Sample", ylab = "Intensity", 
	main = "Gene Profile Plot of IGLJ3"
	)
lines(
	x.line, y.line.num, col = 'blue'
	)
grid(
	col = "grey", lty = "dotted"
	)  


#######################
#          8          #
####################### 

f <- c(rep("SLE",17),rep("Control",9))
y <- dat["211881_x_at",]
y.num <- as.numeric(y)

boxplot(
	y.num ~ f,
	xlab = "Sample",
	ylab = "Intensity",
	main = "Gene Profile Boxplot of IGLJ3"
	)