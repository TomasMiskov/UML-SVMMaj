# Name: SVMMaj_lib
# Author: Tomas Miskov
# Date: 04/01/22
# Version: 1.0
# Purpose: Functions librarY for minimizing SVM loss using majorization
#----------------------------------------------------------------------

#' Absolute SVM Loss Function
#' 
#' Compute the absolute SVM loss function
#' 
#' @param mX Matrix of regressors including a column of 1s
#' @param vY Vector of outcome values (2 categories, -1 & 1)
#' @param vV Vector of weights (+ a constant)
#' @param dLambda Ridge penalization parameter lambda
AbsLossSVM <- function(mX, vY, vV, dLambda){
  vQ <- mX %*% vV                                                         
  dLoss <- sum(pmax(0, 1 - vY * vQ)) + dLambda * (t(vV[-1]) %*% vV[-1])
  return(dLoss)
}

#' Quadratic SVM Loss Function
#' 
#' Compute the quadratic SVM loss function
#' 
#' @param mX Matrix of regressors including a column of 1s
#' @param vY Vector of outcome values (2 categories, -1 & 1)
#' @param vV Vector of weights (+ a constant)
#' @param dLambda Ridge penalization parameter lambda
QuadLossSVM <- function(mX, vY, vV, dLambda){
  vQ <- mX %*% vV                                                         
  dLoss <- sum(pmax(0, 1 - vY * vQ)^2) + dLambda * (t(vV[-1]) %*% vV[-1])
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
#' @param hinge Type of hinge error (absolute or quadratic)
#' @param silent Boolean argument, if FALSE, function prints the iterations
#' @export
MajSVM <- function(mX, vY, dC = 1, vW = rep(1, ncol(mX) - 1), 
                   dLambda = 1, dEpsilon = 10^(-6),
                   hinge = c('absolute', 'quadratic'), silent = TRUE){
  iN <- length(vY)
  iP <- ncol(mX)
  vV0 <- c(dC, vW)                      
  dStartL <- AbsLossSVM(mX, vY, vV0, dLambda)          #starting SVM loss value
  dL0 <- dStartL
  dDecrease <- 0
  mP <- diag(iP)
  mP[1,1] <- 0
  
  if(hinge == 'quadratic'){
    mZ <- solve(t(mX) %*% mX + dLambda * mP) %*% t(mX)
  }
  
  iK <- 1
  while((iK == 1) | dDecrease > dEpsilon){
    
    vQ <- mX %*% vV0  
    if(hinge == 'absolute'){
      vA <- ifelse(abs(vQ) > 1.001, 1/(4 * abs(1 - vY * vQ)), dEpsilon/10)   
      vB <- vY * (vA + 0.25)
      mA <- diag(as.vector(vA))
      vV1 <- solve(t(mX) %*% mA %*% mX  + dLambda * mP, t(mX) %*% vB)
      dL1 <- AbsLossSVM(mX, vY, vV1, dLambda)
    }
    
    if(hinge == 'quadratic'){
      vB <- ifelse((vY == -1 & vQ > -1) | (vY == 1 & vQ > 1), vY * 1, vQ)
      vV1 <- mZ %*% vB
      dL1 <- QuadLossSVM(mX, vY, vV1, dLambda)
    }
    
    dDecrease <- (dL0 - dL1)/dL0
    
    dL0 <- dL1
    vV0 <- vV1
    
    if (silent == FALSE) {
      cat("Iteration: ", iK, "Loss SVM: ", dL1, 
          "Relative difference loss SVM: ", dDecrease,"\n")
    }
    
    iK <- iK + 1
  }
  return(list('c' = vV0[1], 'w' = vV0[-1], 'startLoss' = dStartL, 'endLoss' = dL0,
              'iter' = iK, 'q' = vQ))
}
