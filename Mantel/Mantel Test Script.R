# Script to create distance matrices of Gram et al Hummock data and to perform Mantel tests
# Created for Multivariate Data Analysis Seminar, Clemson University
# blb 9/11/2007
# Modified for subsequent classes (2008, 2010) at CU
# Modified for use at VT by blb 1/20/2012 and 3/29/2014

# sets the working directory, the directory that is accessed for data sets.  

setwd("C:/Users/Heyryanmoore_LT/Documents/01 VT/01 Fall 21/MultivariateStat/07_Mantel/LabMaterials")

###                   
# Section for reading in data from raw data files
###

# Reads in the Geographic dataset, places it in a data matrix (a "frame" in R)
# The "header" option automatically assigns columns names from the data file
# This statement also adds row names (i.e., species)
Geographic <- data.frame(read.table(file="Geographic.csv", sep = ',',  header=TRUE, row.names=c(1)))
 
# Reads in the Soils dataset, places it in a data matrix 
 Soils <- data.frame(read.table(file="Soils.csv", sep = ',',  header=TRUE, row.names="Hummock"))
 # Renames the variables something a bit more meaningful
   names(Soils) <- c('OM', 'Phos', 'NaHCO3', 'TotalN', 'K', 'Mg', 'Ca', 'Ca:Mg', 'Na', 'pH', 'CEC', 'Sand', 'Silt', 'Clay')
   
# Species data
 Comm <- read.csv(file='sppbysmp.csv', header=TRUE)
   # remove the ID column
   Comm <- Comm[,-1]
 
# Hummock measurements
 Humm <- read.csv(file='HummFactors.csv', header=TRUE)
   # remove the ID column
   Humm <- Humm[,-1]

   
   ###
   #  Standardization section
   #  script for standarizing the raw data prior to calculating distance matrices
   #  The 'center' option does the subtraction of mean
   #  The 'scale' option does the dividing by standard deviation
   #  Remmove the "#" to un-comment the lines to apply standardization as you see fit
   ###
   Dist.Stand <- scale(Geographic, center=TRUE, scale=TRUE)
    Soils.Stand <- scale(Soils, center=TRUE, scale=TRUE)
   Comm.Stand <- scale(Comm, center=TRUE, scale=TRUE)
   Humm.Stand <- scale(Humm, center=TRUE, scale=TRUE)

### 
# Distance metric section
###

library(vegan)
#
# 'vegdist' is the function for converting data matrix to a distance matrix
# 'vegdist' has many options; check out the 'vegan' documentation for how to chose distance metrics
# In this case 'Dist.Geo' is an output matrix with distance metrics for the Geographic data set

Dist.Geo <- vegdist(Geographic, method='euclidean')
Dist.Geo.st = vegdist(Dist.Stand, method='euclidean')
# With the Soils data....
Dist.Soils <- vegdist(Soils, method='euclidean')
Dist.Soils.st = vegdist(Soils.Stand, method='euclidean')
# Community Data
Dist.Comm <- vegdist(Comm, method='jaccard')
# Dist.Comm.st <- vegdist(Comm.Stand, method='jaccard') # standardized, will want to use a diff. distance method
# Hummock variables
Dist.Humm <- vegdist(Humm, method='euclidean')
Dist.Humm.st <- vegdist(Humm.Stand, method='euclidean')
#########################
##### MANTEL TEST ########
##########################

# Commands for the regular Mantel test
# We are saving the Mantel results in a variable 'A'
# For the code below, you would simply type 'A' at a command line to show output

A <- mantel(Dist.Geo, Dist.Soils, method='pearson', permutations=1000)  
A.st = mantel(Dist.Geo.st, Dist.Soils.st, method='pearson', permutations=1000)
A.soils.st = mantel(Dist.Geo, Dist.Soils.st, method='pearson', permutations=1000)

# geo dist v comm. comp ####
# mantel test 2 in lab
gdcc = mantel(Dist.Geo, Dist.Comm, method = 'pearson', permutations = 1000)

# soils v comm. comp
# mantel test 3 in lab
scc = mantel(Dist.Soils.st, Dist.Comm, method = 'pearson', permutations = 1000)


# And here is script for the partial Mantel test

# standardized matrix for soils used
B <- mantel.partial(Dist.Comm, Dist.Soils.st, Dist.Geo, method='pearson', permutations=1000)

#####################################
##### VARIATION PARTITIONING ########
####################################

# removing some variables from the Soils matrix; otherwise it's too big and we run out of degrees of freedom
Soils.red <- Soils[,-c(8,9,11,13)]

# variation partitioning 
vp.comm <- varpart(Dist.Comm, Soils.red, Geographic)

# Calculating PCNMs
 # first need to convert our lat-longs to actual space
 # use the function distm() from the package geosphere
library(geosphere)
geo.euc <- distm(Geographic)

 # Now we use the PCNM function in vegan to calculate PCNM vectors
pcnms <- pcnm(geo.euc)

# variation partitioning using PCNMs for spatial data
vp.comm.2 <- varpart(Dist.Comm, Soils.red, pcnms$vectors[,1:3])

vp.comm.3 <- varpart(Soils.red, Dist.Comm, pcnms$vectors[,1:3])
