# Script to load the invertebrate data from Brown and Swan 2010
#  blb 2/9/2014, modified from blb 1/6/2008

library(vegan)  # vegan is the package with the distance functions we'll be using.  
library(dplyr)

setwd("C:/Users/Heyryanmoore_LT/Documents/01 VT/01 Fall 21/MultivariateStat/02_MultivariateSimilarity/LabMaterials")

Set = data.frame(read.table(file='YSC Inverts.csv', 
                            sep=',', header=TRUE))
Set.sim = Set[,c(47,33)]

# using the vegdist function to examine similarity between two observations
dist.a <- vegdist(Set[1:2,-c(1:7)], method='bray', binary=TRUE)

# using the vegdist function to produce a similarity matrix
dist.b <- vegdist(Set[,-c(1:7)], method='bray', binary=TRUE)

dist.b <- as.matrix(vegdist(Set[,-c(1:7)], method='bray', binary=TRUE))

dist.b[4,5]

dist.J <- vegdist(Set[,-c(1:7)], method='jaccard', binary=TRUE)

dist.diff = dist.b - dist.J

# rm(list=ls()) # a line that clears the workspace when you're finished
              #  you can comment them out (using #) if you don't want them for the moment

Set.YG = Set %>% 
  filter(BASIN == "YG")
dist.YG = vegdist(Set.YG[,-c(1:7)], method = 'bray', binary = TRUE)




