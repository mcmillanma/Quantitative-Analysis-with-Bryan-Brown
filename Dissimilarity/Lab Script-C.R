# Script to load the results of Lab Script-B.R
# Provides you with a nice dataset to analyze
#  blb 2/9/2014, modified from blb 1/6/2008

setwd('C:/Users/Heyryanmoore_LT/Documents/01 VT/01 Fall 21/MultivariateStat/02_MultivariateSimilarity/LabMaterials')

load(file='Invert Similarity.RData') # loads the results from 'Lab Script-B.r'
                                     # the data.frame is named "ALL"


# Simple plotting script to examine some of the results
   #  You can modify this script for your own purposes

par(mfrow=c(1,2))
plot(ALL$Network[ALL$Order==1], ALL$SIM[ALL$Order==1], xlab='Distance', ylab='Community Similarity', main='First Order')
  reg.1 <- lm(SIM[ALL$Order==1]~Network[ALL$Order==1], data=ALL)  # calculating a simple linear regression
  abline(reg.1$coefficients)   # Using the regression to draw a slope on the figure
plot(ALL$Network[ALL$Order!=1], ALL$SIM[ALL$Order!=1], xlab='Distance', ylab='', main='Higher Order')
  reg.2 <- lm(SIM[ALL$Order!=1]~Network[ALL$Order!=1], data=ALL)
  abline(reg.2$coefficients)
  
# a <- objects()  #  Notice these lines are commented out right now.  If they weren't, you wouldn't get anything!
#  rm(list=a)