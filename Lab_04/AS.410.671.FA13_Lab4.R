# AS.410.671.81.FA13 – Lab #4
# Normalization and Bioconductor 
# Author: Shan Sabri

#######################
#          1          #
####################### 

library(marray)
dat <- data(swirl)


#######################
#          2          #
####################### 

swirl.three <- swirl[ ,3]
maPlot(swirl.three, main = "MvA Plot of Array 3", lines.func = NULL, legend.func = NULL)


#######################
#          3          #
####################### 

swirl.Norm <- maNorm(swirl.three, norm = c("median"))


#######################
#          4          #
####################### 

maPlot(
	swirl.Norm, 
	main = "MvA Plot of Array 3 (Median Normalized)", 
	lines.func = NULL, 
	legend.func = NULL
	)


#######################
#          5          #
####################### 

# The normalized plot has shifted up. 
# The normalized plot is not symmetric about M = 0. 


#######################
#          6          #
####################### 

swirl.Loess <- maNorm(swirl.three, norm = c("loess"))
maPlot(
	swirl.Loess, 
	main = "MvA Plot of Array 3 (Loess Normalized)", 
	lines.func = NULL, 
	legend.func = NULL
	)


#######################
#          7          #
####################### 

# Loess. 
# The Lowess plot curvature is nearly eliminated and linearity is strong. 


#######################
#          8          #
####################### 

dir.path <- "C:\\Users\\Shan\\Desktop\\JHU_Fall2013\\410.671_Microarrays&Analysis\\DataSets\\"
a.cdna <- read.GenePix(
	path    = dir.path,
	name.Gf = "F532 Median", name.Gb = "B532 Median", 
	name.Rf = "F635 Median", name.Rb = "B635 Median",
	name.W  ="Flags"
	)

	
#######################
#          9          #
####################### 

patient.one       <- a.cdna[ ,1]
patient.one.Loess <- maNorm(patient.one, norm = c("printTipLoess"))
patient.one.PTM   <- maNorm(patient.one, norm = c("scalePrintTipMAD"))
par(mfrow = c(3,1))
maPlot(patient.one, main = "Patient 1 - No Normalization", lines.func = NULL, legend.func = NULL)
maPlot(patient.one.Loess, main = "Patient 1 - Print-tip Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(patient.one.PTM, main = "Patient 1 - Scale Print-tip Normalization using MAD", lines.func = NULL, legend.func = NULL)

patient.two       <- a.cdna[ ,2]
patient.two.Loess <- maNorm(patient.two, norm = c("printTipLoess"))
patient.two.PTM   <- maNorm(patient.two, norm = c("scalePrintTipMAD"))
par(mfrow = c(3,1))
maPlot(patient.two, main = "Patient 2 - No Normalization", lines.func = NULL, legend.func = NULL)
maPlot(patient.two.Loess, main = "Patient 2 - Print-tip Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(patient.two.PTM, main = "Patient 2 - Scale Print-tip Normalization using MAD", lines.func = NULL, legend.func = NULL)


#######################
#          10         #
####################### 

data       <- cbind(patient.one, patient.two)
data.Loess <- cbind(patient.one.Loess, patient.two.Loess)
data.PTM   <- cbind(patient.one.PTM, patient.two.PTM)

probe.ids  <- data@maGnames@maLabels

Loess <- maM(data.Loess)
PTM   <- maM(data.PTM)
# dim(Loess)
# dim(PTM)

dimnames(Loess)[[1]] <- probe.ids
dimnames(PTM)[[1]]   <- probe.ids


#######################
#          11         #
####################### 

library(affy)
library(limma)
library(simpleaffy) 
library(affyPLM) 
library(fpc)


#######################
#          12         #
####################### 

fns       <- sort(list.celfiles(path = dir.path, full.names = TRUE))
data.affy <- ReadAffy(filenames = fns, phenoData = NULL)


#######################
#          13         #
####################### 

affy.MAS <- justMAS(data.affy)
affy.rma <- call.exprs(data.affy, "rma")

dim(affy.rma) 
dim(affy.MAS)


#######################
#          14         #
####################### 

r <- cor(exprs(affy.rma))
m <- cor(exprs(affy.MAS))

r
m