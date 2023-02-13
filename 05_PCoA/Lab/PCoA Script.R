# Script for performing Principal Coordinates on Data from Wepking et al in Proceedings B
# blb 9/28/2019

# set working directory
setwd("C:/Users/Heyryanmoore_LT/Documents/01 VT/01 Fall 21/MultivariateStat/05_PCoA/Lab")

# Import Karl's data
Design <- read.csv(file='Design.csv', header=TRUE) # the design structure of the experiment
Micros <- read.csv(file='Karls Data.csv', header=TRUE) # the microbial OTU data
Micros <- Micros[Micros$Tag!=70,] # removing this row because it was used as an internal reference

# running the PCoA
library(ape) # load the ape package with a nifty PCoA function
library(vegan) # load the vegan package for vegdist
library(dplyr)
library(tidyr)

dim(Micros)
Micros_g = Micros %>%
   pivot_longer(c(1:13446)) # turning columns into rows
unique(Micros_g)
Micros_OTU = Micros_g %>%
   filter(name != "Tag") # removing 22 "Tag" observations

Micros_val = as.data.frame(as.numeric(gsub("OTU_", "", Micros_OTU$name)))
max(Micros_val)
unique(Micros_val) # checking that 13445 OTUs are present in the dataset, no skips

micro.dist <- vegdist(log(Micros[,-1]+1), method='bray', binary=FALSE) # producing a distance matrix using vegdist, notice I applied a log transformation in there like in Wepking et al. ####
micro.dist = as.data.frame(as.table(micro.dist))
histogram(micro.dist$Freq)

micro.dist1 <- vegdist(Micros, method='bray', binary=FALSE)
micro.dist1 = as.data.frame(as.table(micro.dist1))
histogram(micro.dist1$Freq)

micro.dist <- vegdist(log(Micros[,-1]+1), method='bray', binary=FALSE)
Micro.pco <- pcoa(micro.dist, correction='none') # running the PCoA ####
summary(Micro.pco)
Micro.pco
Micro.pco$values
Micro.pco$vectors


micro.dist1 <- vegdist(Micros, method='bray', binary=FALSE)
Micro.pco1 = pcoa(micro.dist1, correction='none')
Micro.pco1$values

# making a figure similar to that in Wepking et al
Z <- data.frame(Design, Micro.pco$vectors[,1:2]) # merging the Design information with the PCoA output for the 1st 2 axes


plot(Z$Axis.1, Z$Axis.2, type='n', xlab='PCoA 1', ylab='PCoA 2', xlim=c(-0.5, 0.35))
   points(Z$Axis.1[Z$Manure=='High'], Z$Axis.2[Z$Manure=='High'], pch=21, bg='black', cex=1.2)
   points(Z$Axis.1[Z$Manure=='Low'], Z$Axis.2[Z$Manure=='Low'], pch=21, bg='grey85', cex=1.2)
   text((Z$Axis.1 + 0.03), Z$Axis.2, labels=Z$Site)
   legend('bottomleft', legend=c('+manure', 'reference'), col=c('black'), pch=21, pt.bg=c('black', 'grey85'), cex=1.2)
   

      
Z1 <- data.frame(Design, Micro.pco1$vectors[,1:2])
    
plot(Z1$Axis.1, Z1$Axis.2, type='n', xlab='PCoA 1', ylab='PCoA 2', xlim=c(-0.5, 0.35))
   points(Z1$Axis.1[Z1$Manure=='High'], Z1$Axis.2[Z1$Manure=='High'], pch=21, bg='black', cex=1.2)
   points(Z1$Axis.1[Z1$Manure=='Low'], Z1$Axis.2[Z1$Manure=='Low'], pch=21, bg='grey85', cex=1.2)
   text((Z1$Axis.1 + 0.03), Z1$Axis.2, labels=Z1$Site)
   legend('bottomleft', legend=c('+manure', 'reference'), col=c('black'), pch=21, pt.bg=c('black', 'grey85'), cex=1.2)   
   
   
   