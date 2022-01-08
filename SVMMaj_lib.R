# Name: SVMMaj_lib
# Author: Tomas Miskov & Loek van Montfort
# Date: 08/01/22
# Version: 3.0
# Purpose: Functions librarY for minimizing SVM loss using majorization
#----------------------------------------------------------------------

#' SVM Loss Function
#' 
#' Compute the SVM loss function
#' 
#' @param mX Matrix of regressors including a column of 1s
#' @param vY Vector of outcome values (2 categories, -1 & 1)
#' @param vV Vector of weights (+ a constant)
#' @param dLambda Ridge penalization parameter lambda
LossSVM <- function(mX, vY, vV, dLambda){
  vQ <- mX %*% vV                                                         
  dLoss <- sum(pmax(0, 1 - vY * vQ)) + dLambda * (t(vV[-1]) %*% vV[-1])
  return(dLoss)
}

#' SVM Quadratic Loss Function
#' 
#' Compute the SVM loss function for quadratic hinge
#' 
#' @param mX Matrix of regressors including a column of 1s
#' @param vY Vector of outcome values (2 categories, -1 & 1)
#' @param vV Vector of weights (+ a constant)
#' @param dLambda Ridge penalization parameter lambda
LossSVMQuadr <- function(mX, vY, vV, dLambda){
  vQ <- mX %*% vV                                                         
  dLoss <- sum(pmax(0, 1 - vY * vQ)**2) + dLambda * (t(vV[-1]) %*% vV[-1])
  return(dLoss)
}
#' SVM majorization
#' 
#' Solve the primal SVM problem using majorization
#' 
#' @param mX Matrix of regressors including a column of 1s
#' @param vY Vector of outcome values (2 categories, -1 & 1)
#' @param dC Initial value of the constant
#' @param vW Vector of weights
#' @param dLambda Ridge penalization parameter lambda
#' @param dEpsilon Accuracy parameter
#' @param sHinge Type of hinge error (absolute or quadratic)
#' @param bSilent Boolean argument, if FALSE, function prints the iterations
#' @export
MajSVM <- function(mX, vY, dC = 1, vW = rep(1, ncol(mX) - 1), 
                   dLambda = 1, dEpsilon = 10^(-6),
                   sHinge = 'absolute', bSilent = TRUE){
  iN <- length(vY)
  iP <- ncol(mX)
  vV0 <- c(dC, vW)                      
  dStartL <- LossSVM(mX, vY, vV0, dLambda)          #starting SVM loss value
  mP <- diag(iP)
  mP[1,1] <- 0
  
  if(sHinge == 'quadratic'){
    mZ <- solve(t(mX) %*% mX + dLambda * mP) %*% t(mX)
    dStartL <- LossSVMQuadr(mX, vY, vV0, dLambda)          #starting SVM loss value
  }
  
  dL0 <- dStartL
  dDecrease <- 0
  iK <- 1
  while((iK == 1) | dDecrease > dEpsilon){
    vQ <- mX %*% vV0  
  
    if(sHinge == 'absolute'){
      vA <- ifelse(abs(1 - vY * vQ) > dEpsilon, 1/(4 * abs(1 - vY * vQ)), 1/(4*dEpsilon))   
      vB <- vY * (vA + 0.25)
      mA <- diag(as.vector(vA))
      vV1 <- solve(t(mX) %*% mA %*% mX  + dLambda * mP, t(mX) %*% vB)
      dL1 <- LossSVM(mX, vY, vV1, dLambda)
    }
    
    if(sHinge == 'quadratic'){
      vB <- ifelse((vY == -1 & vQ > -1) | (vY == 1 & vQ < 1), vY * 1, vQ)
      vV1 <- mZ %*% vB
      dL1 <- LossSVMQuadr(mX, vY, vV1, dLambda)
    }
    
    dDecrease <- (dL0 - dL1)/dL0
    
    if (bSilent == FALSE) {
      cat("Iteration: ", iK, "Loss SVM: ", dL1, 
          "Relative difference loss SVM: ", dDecrease,"\n")
    }
    dL0 <- dL1
    vV0 <- vV1
    iK <- iK + 1
  }
  return(list('v' = vV0, 'startLoss' = dStartL, 'endLoss' = dL0,
              'iter' = iK))
}
