# Script to load the invertebrate data + distance data from Brown and Swan 2010
# Also Combines the physical distances with the site codes and incorporates them into the "Set"
#  blb 2/9/2014, modified from blb 1/6/2008

# Import
library(vegan)  # vegan is the package with the distance functions we'll be using.  
library(dplyr)

setwd("C:/Users/Heyryanmoore_LT/Documents/01 VT/01 Fall 21/MultivariateStat/02_MultivariateSimilarity/LabMaterials")

Set = data.frame(read.table(file='YSC Inverts.csv', 
                            sep=',', header=TRUE))

Codes <- data.frame(read.table(file='Codes.csv', sep=',', header=TRUE))
DATA <- merge(Codes, Set, by="SITE") # merges the bug data with the site codes 
# that match the 'Distance' dataset 

Distance <- data.frame(read.table(file='Distance.csv', sep=',', header=TRUE))
names(Distance)[4:6] <- c('Network', 'Euclidean', 'Difference')  

###################  
Metric <- 'bray'     # Choice of distance metric for the similarity from the 'vegdist' function
                     #  several possible choices;  use R help to investigate
################### 
##################
Stand <- TRUE         #  Chooses to standardize distances within each watershed                                    
###################    # Raw distances are on fairly different scales
##################
Bin <- FALSE         # Turns on the BINARY=TRUE option in vegdist
##################   # use primarily with Jaccard

                                                                                  
# Producing distances for the invert data

#  data.frame 'DATA' has site IDs in columns 1:12, and species counts in the remainder
a <- DATA[, -c(1:12)]

# Compute similarites by looping to produce independent ID variables and thus ensure that 
# we get correct correspondence between Distance and Similarity
Similar <- c()
Order.list <- unique(DATA$ORDER)
 for(k in 1:length(Order.list)){
     DATAz <- DATA[DATA$ORDER==Order.list[k],]
     Code.list <- unique(DATAz$Code)                            
     Simile <- c()
     for(i in 1:length(Code.list)){
        a <- DATAz[DATAz$Code==Code.list[i],-c(1:12)]
        Code.list2 <- Code.list[(i-1):length(Code.list)]
        for(j in 1:length(Code.list2)){
        b <- DATAz[DATAz$Code==Code.list2[j],-c(1:12)]
        ab <- rbind(a,b)
        sim <- vegdist(ab, method=Metric, binary=Bin)  # the similarity function
           SIM <- 1-(as.numeric(sim))    # converting from default dissimilarity to similarity
        out <- data.frame(START=Code.list[i], END=Code.list2[j], WS=DATAz$SHEDNAME.x[DATAz$Code==Code.list[i]], Order=Order.list[k], SIM)
        Simile <- rbind(Simile,out)  
          }
        }
       Similar <- rbind(Similar, Simile) 
     }

# Output datasets    
mrg <- merge(Distance, Similar, all.x=TRUE) # merging the distances with the similarities
ALL <- mrg[which(is.na(mrg$SIM)!=TRUE),] 

sit = TRUE
# standardizing distances between watersheds
if(sit==TRUE){
  all.stor <- c()
  ws.list <- unique(ALL$WS)
  for(i in 1:length(ws.list)) {
    a <- ALL[ALL$WS==ws.list[i],]
    mu.euc <- mean(a$Euclidean)
    mu.net <- mean(as.numeric(a$Network))
    a$Euclidean <- a$Euclidean/mu.euc
    a$Network <- a$Network/mu.net  
    all.stor <- rbind(all.stor, a)  
    }
}           
ALL <- all.stor   
   save(ALL, file='Invert Similarity.RData')

# a <- objects()
# rm(list=a)

Order.list <- unique(DATA$ORDER)
for(k in 1:length(Order.list)){
  DATAz <- DATA[DATA$ORDER==Order.list[k],]
  Code.list <- unique(DATAz$Code)                            
  Simile <- c()
  for(i in 1:length(Code.list)){
    a <- DATAz[DATAz$Code==Code.list[i],-c(1:12)]
    Code.list2 <- Code.list[(i-1):length(Code.list)]
    for(j in 1:length(Code.list2)){
      b <- DATAz[DATAz$Code==Code.list2[j],-c(1:12)]
      ab <- rbind(a,b)
      sim <- vegdist(ab, method='bray', binary=FALSE)  # the similarity function
      SIM <- 1-(as.numeric(sim))    # converting from default dissimilarity to similarity
      out <- data.frame(START=Code.list[i], END=Code.list2[j], WS=DATAz$SHEDNAME.x[DATAz$Code==Code.list[i]], Order=Order.list[k], SIM)
      Simile <- rbind(Simile,out)  
    }
  }
  Similar <- rbind(Similar, Simile) 
}
   

