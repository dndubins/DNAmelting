#SVD - NN(6).R
#This script fits a 6-parameter nearest-neighbour model, where 
#the identity of the bases is important but directionality is not considered, 
#for a Propery of DNA duplex. Singular Value Decomposition is used
#(SVD).
#Written by: David Dubins 
#Date: Jan 20, 2026
#Platform: R-Studio 2026.01.0 Build 392, on R version 4.5.2
#External Libraries used: corrplot, stringr from R package library
#This program expects two columns of data:
# Label Sequence  Yobs

library(corrplot) #install the corrplot package before loading for the first time
library(stringr) # install the stringr package before loading for the first time

#Read in the data
#For PC:
#Dataset <- read.table(file="clipboard", header=TRUE, sep="\t")
#For MacOS:
#Dataset <- read.table(pipe("pbpaste"), header=TRUE, sep="\t")

#Change this label to match the observed data you are fitting
Ylabel <- "ΔSASA"
#These are the labels of the nearest-neighbour parameters for this model:
NNlabels <-c(
  "TT,AA","TA,AT", "TG,GT,AC,CA", "AG,GA,TC,CT", "GG,CC", "GC,CG"
)
  
#Function to count the 6 NN types in each oligo:
count_nn6 <- function(seq) {
  # Clean sequence
  s <- toupper(gsub("/", "", seq))
  
  # Extract nearest neighbors
  pairs <- substring(s, 1:(nchar(s) - 1), 2:nchar(s))
  
  # Define NN classes
  class_map <- list(
    TT_AA = c("TT", "AA"),
    TA_AT = c("TA", "AT"),
    TG_GT_AC_CA = c("TG", "GT", "AC", "CA"),
    AG_GA_TC_CT = c("AG", "GA", "TC", "CT"),
    GG_CC = c("GG", "CC"),
    GC_CG = c("GC", "CG")
  )
  
  # Count each class
  counts <- sapply(class_map, function(nn) sum(pairs %in% nn))
  
  return(counts)
}

#Apply to data:
X <- t(sapply(Dataset$Sequence, count_nn6))
colnames(X) <- NNlabels
X <- as.matrix(X)

#Model is Y = X*beta; beta=(6 NN parameters)
#The fitting routine used here is SVD (Singular Value Decomposition) 
#Solve using SVD:
svd_X <- svd(X)

# Moore–Penrose pseudoinverse
X_pinv <- svd_X$v %*% diag(1 / svd_X$d) %*% t(svd_X$u)
beta <- X_pinv %*% Dataset$Yobs
colnames(beta) <- Ylabel
rownames(beta) <- NNlabels

# Construct the hat matrix:
Y_hat <- X %*% beta
# Plot predicted vs. observed values (hopefully a straight line)
plot(Dataset$Yobs, Y_hat,
     xlab = paste("Observed", Ylabel, sep = " "), #combine strings
     ylab = paste("Predicted", Ylabel, sep = " "),
     pch = 19)
abline(0, 1, col = "red")

# Examine the residuals
resid <- Dataset$Yobs - Y_hat
summary(resid)

#Calculate the Covariance Matrix
n <- nrow(X)
p <- ncol(X)
sigma2 <- sum(resid^2) / (n - p)
Dinv2 <- diag(1 / svd_X$d^2)
cov_beta <- sigma2 * svd_X$v %*% Dinv2 %*% t(svd_X$v)
rownames(cov_beta) <- colnames(cov_beta) <- NNlabels
cov_beta

#Concert Covariance Matrix to the Correlation Matrix.
#Look for off-diagonal numbers > 0.8 here - that's bad!
sd_beta <- sqrt(diag(cov_beta))
cor_beta <- cov_beta / (sd_beta %*% t(sd_beta))
cor_beta

#Here's a colour plot of the correlation matrix
corrplot::corrplot(
  cor_beta,
  method = "color",
  type = "upper",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45
)

# Calculate the standard errors of parameter estimates
se_beta <- sqrt(diag(cov_beta))
names(se_beta) <- NNlabels
# Report the parameter estimates and their standard errors:
results <- data.frame(
  Estimate  = round(as.numeric(beta),2), #Adjust rounding as needed
  StdError  = round(se_beta,2)
)
colnames(results) <- c("Estimate","±StdErr" )
results

#Calculate and report the rank of the model. 
#This should match the number of parameters, otherwise 
#the SVD algorithm reduced the effective number of independent fit parameters.
tol <- max(dim(X)) * max(svd_X$d) * .Machine$double.eps
rankX <- sum(svd_X$d > tol)
rankX

#Calculate the chi^2 for the model fit: (smaller is better)
chisq <- sum(resid^2) #chisq here is really RSS (Residual Sum of Squares)
#if you have sigma, import into R and use this formula instead, which is the true chisq:
#chisq <- sum((resid/sigma)^2)
chisq

#Calculate the r^2 for the model fit:
SS_res <- sum(resid^2)
SS_tot <- sum((Dataset$Yobs - mean(Dataset$Yobs))^2)
rsq <- 1 - SS_res / SS_tot
rsq

#Calculate the Adjusted r^2 for the model fit:
rsq_adj <- 1 - (1 - rsq) * (n - 1) / (n - p - 1)
rsq_adj
