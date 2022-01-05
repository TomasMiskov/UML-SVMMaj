# Name: SVMMaj_lib
# Author: Tomas Miskov
# Date: 04/01/22
# Version: 1.0
# Purpose: Functions library for minimizing SVM loss using majorization
#----------------------------------------------------------------------

Lsvm <- function(X, y, dC, vW, dLambda){
  vQ <- dC + X %*% vW                                            #compute the Q vector
  vPosInd <- y == 1
  dLoss <- sum(pmax(0, 1 + vQ[-vPosInd])) + sum(pmax(0, 1 - vQ[vPosInd])) + dLambda * (t(vW)%*%vW)
  return(dLoss)
}

MajSVM <- function(X, y, dLambda = 1){
  iN <- length(y)
  iP <- ncol(X)
  vPosInd <- y == 1     #TRUE/FALSE vector of y == 'yes'/1
  dC <- 1                            #initial constant c
  vW <- rep(1, iP)                   #initial vector of weights w
  vV <- c(dC, vW)                    #initial vector v consisting of c & w
  X <- cbind(rep(1,iN), X)           #add column of 1s to X
  dL0 <- Lsvm(X, y, dC, vW, dLambda) #starting SVM loss value
  
  epsilon <- 1e^5
  k <- 1
  while(k = 1 | (dLk - dlK1)/dlk > epsilon){
    k <- k + 1
    
    vQ <- X %*% vV                                                 #compute vQ values
    vErrors <- ifelse(y == -1, pmax(0, vQ + 1), pmax(0, vQ - 1))   #corresponding hinge errors
    vA <- ifelse(y == -1, 0.25*abs(vErrors + 1)^-1, 0.25*abs(1 - vErrors)^-1)
    
    mA <- diag(vA)
    vB <- ifelse(y == -1, -vA - 0.25, vA + 0.25)
    vV <- solve()
    
  }
}