# Name: SVMMaj
# Author: Tomas Miskov & Loek van Montfort
# Date: 08/01/22
# Version: 3.0
# Purpose: Minimize SVM loss function using majorization, compare with SVMMaj
#----------------------------------------------------------------------------
#------
# Imports
#------
rm(list=ls())
if (!require("pacman")) install.packages("pacman")
if (!require("dplyr")) install.packages("dplyr")
if (!require("caret")) install.packages("caret")
pacman::p_load(SVMMaj)
source("SVMMaj_lib.R")
#------
# Functions
#------
#' Performance of the predictions
#' 
#' Compute the performance of the predictions
#'
#' @param vY Vector of outcome values (2 categories, -1 & 1) of the observed class
#' @param vYHat Vector of outcome values (2 categories, -1 & 1) of the predicted class

fPredPerf <- function(vY, vYHat){
  lConMatrix <- confusionMatrix(data = factor(vYHat), 
                                reference = factor(vY),
                                positive = "1")
  dMissClassRate <- sum(vY != vYHat)/length(vY)
  dHitRate <- 1-dMissClassRate
  
  return(list("ConMatrix" = lConMatrix, "MissClassRate" = dMissClassRate, "HitRate" = dHitRate))
}

#' Print the performance of the SVMMaj package and the own implementation
#' 
#' Print the performance of the SVMMaj package and the own implementation, and 
#' compare the results of both methods. 
#'
#' @param lPredPerfImplement List of performance indicators of the implementation
#' @param lPredPerfPackage List of performance indicators of the SVMMaj package
#' @param modelPackage Model estimated with the SVMMaj package
#' @param sHinge String indicating the type of hinge error
#' @param modelImplement Model results of the implementation

fPrintPredPerf <- function(lPredPerfImplement, lPredPerfPackage, modelPackage, sHinge, modelImplement){
  cat("===== Estimation properties SVMMaj package ===== \n")
  print(lPredPerfPackage)
  cat("=== Estimation properties own implementation === \n")
  cat("Confusion matrix:")
  print(lPredPerfImplement$ConMatrix)
  cat(".\n") 
  cat("Misclassification rate: " , lPredPerfImplement$MissClassRate, ".\n") 
  cat("Hit rate: " , lPredPerfImplement$HitRate, ".\n") 
  cat("============ Comparison methods (implementation vs package) ==============\n")
  cat("Coefficients: ")
  dfResults <- data.frame(implementation = round(modelImplement$v,5), package = round(modelPackage$beta,5))
  print(dfResults)
  cat("Iterations: ", modelImplement$iter, modelPackage$iteration,".\n")
  cat("Loss function: ", modelImplement$endLoss, modelPackage$loss,".\n")
  cat("Type of loss function: ", sHinge,".\n")
}

#' Prepare the data
#' 
#' Clean up the data, by the addition of dummies and scaling of the columns. 
#' Obtain a train and a test set.
#'
#' @param dfData Dataframe containing the data
#' @param vInd Vector of indices of the sample
#' @param iNTrain Int indicating the size of the training set

fPrepData <- function(dfData, vInd, iNTrain){
  dfData <- dfData[vInd,]
  dfData <- select(dfData, -c('emp.var.rate', 'euribor3m'))

  vY <- ifelse(dfData$y == 'yes', 1, -1)
  mX <- model.matrix(y ~., data = dfData)
  
  mX <- mX[,!colnames(mX) %in% c('educationilliterate', 'defaultyes', 'monthdec')]
  mXTrain <- mX[1:iNTrain,]
  mXTest <- mX[(iNTrain+1):length(vInd),]
  mXTrainStd <- scale(mXTrain, center = TRUE, scale = TRUE) 
  mXTestStd <- scale(mXTest, center = attr(mXTrainStd,"scaled:center"), 
                     scale = attr(mXTrainStd,"scaled:scale")) 
  mXTrainStd[,1] = 1
  mXTestStd[,1] = 1
  
  return(list('XTrain' = mXTrainStd, 'YTrain' = vY[1:iNTrain], 
              'XTest' = mXTestStd, 'YTest' = vY[(iNTrain+1):length(vInd)]))
}
#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------
#------
# Magic numbers
#------
dLambda <- 1
dEpsilon <- 10^-5
bSilent <- TRUE
iNTrain <- 2000
iNTest <- 1000
iSeed <- 235167
sHinge <- 'absolute'                                 # 'absolute' or 'quadratic'
#------
# Initialisation
#------
set.seed(iSeed)
load('bank.RData')
vInd <- sample(1:nrow(bank),size = iNTrain + iNTest, replace = FALSE)
lData <- fPrepData(bank, vInd, iNTrain)
#--------------
# Estimation
#--------------
#--------------
# SVMMaj package
#--------------
modelPackage <- svmmaj(lData$XTrain, lData$YTrain, hinge = sHinge, 
                    convergence = dEpsilon, lambda = dLambda, scale = "none")
lPredPerfPackage <- predict(modelPackage, lData$XTest, lData$YTest)
#-------------------
# Our implementation
#-------------------
modelImplement <- MajSVM(lData$XTrain, lData$YTrain, dEpsilon = dEpsilon, sHinge = sHinge, bSilent = bSilent, dLambda = dLambda)
vQPred <- lData$XTest %*% modelImplement$v
vYPred <- ifelse(vQPred > 0, 1, -1)
lPredPerfImplement <- fPredPerf(lData$YTest, as.numeric(vYPred))
#--------------
# Output
#--------------
fPrintPredPerf(lPredPerfImplement, lPredPerfPackage, modelPackage, sHinge, modelImplement)
