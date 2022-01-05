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
path <- "C:\\Users\\misko\\OneDrive\\Desktop\\BDS\\Block 3\\Unsupervised Machine Learning\\Week 1\\bank.RData"
load(path)
y <- bank$y
X <- model.matrix(y ~., data = bank)

#--------------
#SVMMaj package
#--------------
benchmark <- svmmaj(X, y)
summary(benchmark)
plot(benchmark)

#-------------------
# Our implementation
#-------------------
y <- factor(y, levels = c('yes', 'no'), labels = c(1, -1))       #re-factor y into factors 1, -1
model <- MajSVM()