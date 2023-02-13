#  Analysis of the lentic macroinvertebrate data from Hampton and Duggan 2003
#  blb 10/6/2008, modified 3/20/2012, 2/25/2014, 9/28/2015 

setwd("C:/Users/Heyryanmoore_LT/Documents/01 VT/01 Fall 21/MultivariateStat/04_NMDS/Lab")
# Read in the data from raw data provided by Hampton and Duggan

# loading the vegan library in order to calculate distance matrices
library(vegan)

DATA <- data.frame(read.table(file='Hampton-Lab.csv', sep=',', header=TRUE))

#  Rename the non-ID part of the matrix;  it will save you some typing later
dat <- DATA[,-c(23:25)]
    
# creating transformed data using the 4th root like in Hampton and Duggan
dat.trans <- dat^(1/4)

# Running the NMDS

#  metaMDS integrates functions from several packages to perform NMDS.....
#  ....including our old friend 'vegdist' frm the vegan package
X <- metaMDS(dat.trans, distance='bray', k=2, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)

names(X)
X$stress
X$points
spp = as.data.frame(X$species)
plot(spp$MDS1, spp$MDS2, type = 'p')

#  Basic plot of all of the points
plot(X, display=c('sites', 'species'), choices=c(1,2), type='p')

plot(X, display=c('sites', 'species'), choices=c(1,2), type='t')

#### https://ourcodingclub.github.io/tutorials/ordination/
env.data = DATA[,c(23:25)]
area = DATA[,25]

ef <- envfit(X, env.data, permu = 999)
ef
plot(X, type = "t", display = "sites")
plot(ef, p.max = 0.05)


plot(X, display=c('sites', 'species'), choices=c(1,2), type='n')

#  Points divided into day vs. night
plot(X, display=c('sites','species'),choices=c(1,2), type='n')
  points(X$points[DATA$time=='night',1], X$points[DATA$time=='night',2], pch=23,bg='black')
  points(X$points[DATA$time=='day',1], X$points[DATA$time=='day',2], pch=21,bg='grey75')
  legend(min(X$points[,1]), max(X$points[,2]), legend=c('Night', 'Day'), pch=23, pt.bg=c('black', 'grey75')) 
  
#  Basic outline for making your own Scree plot
k <- c(1:6)

k1 <- metaMDS(dat.trans, distance='bray', k=1, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
k2 <- metaMDS(dat.trans, distance='bray', k=2, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
k3 <- metaMDS(dat.trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
k4 <- metaMDS(dat.trans, distance='bray', k=4, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
k5 <- metaMDS(dat.trans, distance='bray', k=5, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
k6 <- metaMDS(dat.trans, distance='bray', k=6, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)
Stress <- c(k1$stress,k2$stress,k3$stress,k4$stress,k5$stress,k6$stress)
plot(k, Stress, type='b', xlab='Dimensions', ylab='Stress', pch=21, bg='grey75', cex=1.3)

# https://ourcodingclub.github.io/tutorials/ordination/
dist = vegdist(dat.trans, method = 'bray')

NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 6), replicate(6, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 6),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:6) {
    points(rep(i + 1,6),replicate(6, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

NMDS.scree(dist)


# metaMDS() plots ####
plot(k1, display=c('sites', 'species'), choices=c(1,2), type='n', main = "k1")
  points(X$points[DATA$time=='night',1], X$points[DATA$time=='night',2], pch=23,bg='black')
  points(X$points[DATA$time=='day',1], X$points[DATA$time=='day',2], pch=21,bg='grey75')
  
plot(k2, display=c('sites', 'species'), choices=c(1,2), type='n', main = "k2")
  points(X$points[DATA$time=='night',1], X$points[DATA$time=='night',2], pch=23,bg='black')
  points(X$points[DATA$time=='day',1], X$points[DATA$time=='day',2], pch=21,bg='grey75')

plot(k3, display=c('sites', 'species'), choices=c(1,2), type='n', main = "k3")
  points(X$points[DATA$time=='night',1], X$points[DATA$time=='night',2], pch=23,bg='black')
  points(X$points[DATA$time=='day',1], X$points[DATA$time=='day',2], pch=21,bg='grey75')

plot(k4, display=c('sites', 'species'), choices=c(1,2), type='n', main = "k4")
  points(X$points[DATA$time=='night',1], X$points[DATA$time=='night',2], pch=23,bg='black')
  points(X$points[DATA$time=='day',1], X$points[DATA$time=='day',2], pch=21,bg='grey75')

plot(k5, display=c('sites', 'species'), choices=c(1,2), type='n', main = "k5")
  points(X$points[DATA$time=='night',1], X$points[DATA$time=='night',2], pch=23,bg='black')
  points(X$points[DATA$time=='day',1], X$points[DATA$time=='day',2], pch=21,bg='grey75')

plot(k6, display=c('sites', 'species'), choices=c(1,2), type='n', main = "k6")
  points(X$points[DATA$time=='night',1], X$points[DATA$time=='night',2], pch=23,bg='black')
  points(X$points[DATA$time=='day',1], X$points[DATA$time=='day',2], pch=21,bg='grey75')

plot(k2, display=c('sites', 'species'), choices=c(1,2), type='n', main = "k2")
  points(X$points[DATA$area=='littoral',1], X$points[DATA$area=='littoral',2], pch=21,bg='black')
  points(X$points[DATA$area=='open',1], X$points[DATA$area=='open',2], pch=21,bg='grey75')
  points(X$points[DATA$area=='edge',1], X$points[DATA$area=='edge',2], pch=21,bg='blue')
  
  # 3 dimension plot
install.packages("vegan3d")
library(vegan3d)
X3 <- metaMDS(dat.trans, distance='bray', k=3, trymax=20, autotransform=FALSE, pc=FALSE, plot=FALSE)


ordiplot3d(X3, display = c("sites",'species'), choices = 1:3, type = 'n')
  points(X$points[DATA$area=='littoral',1], X$points[DATA$area=='littoral',2], pch=21,bg='black')
  points(X$points[DATA$area=='open',1], X$points[DATA$area=='open',2], pch=21,bg='grey75')
  points(X$points[DATA$area=='edge',1], X$points[DATA$area=='edge',2], pch=21,bg='blue')
  legend(5, 6, legend=c('Littoral','Open', 'Edge'), pch=23, pt.bg=c('black', 'grey75', 'blue'))
  
  ordiplot3d(X3, display = c("sites",'species'), choices = 1:3, col = "black",
           ax.col = "red", arr.len = 0.1, arr.col = "blue", envfit,
           xlab, ylab, zlab, ...)