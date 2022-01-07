# Name: SVMMaj
# Author: Tomas Miskov
# Date: 04/01/22
# Version: 1.0
# Purpose: Minimize SVM loss function using majorization
#-------------------------------------------------------
rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(SVMMaj)

#------
# DATA
#------
load('bank.RData')
prepData <- function(data, vInd){
  vY <- ifelse(data$y == 'yes', 1, -1)
  mX <- model.matrix(y ~., data = data)
  
  vY <- vY[vInd]
  mX <- mX[vInd,]
  
  return(list('X' = mX, 'Y' = vY))
}

vInd <- sample.int(nrow(bank), 500)
lData <- prepData(bank, vInd)

#--------------
#SVMMaj package
#--------------
benchmark <- svmmaj(lData$X, lData$Y, hinge = 'absolute', 
                    convergence = 10^-4)
summary(benchmark)
plot(benchmark)

#-------------------
# Our implementation
#-------------------
source("SVMMaj_lib.R")
model <- MajSVM(lData$X, lData$Y, dEpsilon = 10^-4, hinge = 'absolute', silent = FALSE)
model$w
benchmark$beta
abs(sum(model$w - benchmark$beta[-1]))
