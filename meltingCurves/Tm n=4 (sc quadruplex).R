#Tm n=4 (sc quadruplex).R
#This script fits TM data assuming Q -> SS1 (self-complementary)
#Written by: David Dubins 
#Date: Jan 17, 2026
#Platform: R-Studio 2026.01.0 Build 392, on R version 4.5.2
#External Libraries used: polynom 1.4 
#(Tools -> Install Packages -> select polynom)
#This program expects two columns of data:
# Temperature   Abs

library(polynom) #install the polynom package before loading for the first time

#Read in the data (Temperature, Abs)
#Uncomment the appropriate line below for your platform.
#Select the data in Excel, press Ctrl+C, then run the appropriate line.
#For PC:
#MeltingCurve <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#MeltingCurve <- read.table(pipe("pbpaste"), header=TRUE)

attach(MeltingCurve)

#Enter the lower and upper baselines for the spectroscopic data (Fit using MS Excel):
Uint <- 0.13844         #upper baseline intercept (SS DNA)
Uslope <- 4.8009e-5     #upper baseline slope (SS DNA)
Lint <- 0.110637311     #lower baseline intercept (Quadruplex)
Lslope <- 0.000147884   #lower baseline slope (Quadruplex)

#First calculate Alpha vs. Temperature curve
#Alpha = (Abs - Lower baseline) / (Upper baseline - Lower baseline)
Alpha <- (Abs-(Lint+Lslope*Temperature))/((Uint+Uslope*Temperature)-(Lint+Lslope*Temperature))

#Define model parameters
nr <- 4  #number of roots (also the order of the polynomial to fit)
points <- dim(MeltingCurve)[1] #number of data points to fit
Ct <- 100e-6  #Initial concentration of SS: 100 uM

#We need this function to solve for alpha by finding a real root between [0 and 1].
findReal <- function(roots){  # finds a real root between 0 and 1, or returns 0.
  validroot <- 0
  for(j in 1:nr){
    if(is.complex(roots[j])==FALSE){
      if(roots[j]>=0 && roots[j]<=1.001){
        #if(roots[j]>=0){
        validroot <- roots[j]
      }
    } 
    if(is.complex(roots[j])==TRUE){
      if(Im(roots[j])==0 && Re(roots[j])>=0 && Re(roots[j])<=1.001){
        #if(Im(roots[j])==0 && Re(roots[j])>=0){
        validroot <- Re(roots[j])
      }
    }
  }
  validroot
}

#Define the function to fit (quadruplex melting to SS1, self-complementary)
Alphafit <- function(Tm, DH) {
  #Lint: Lower baseline (intercept)
  #Lslope: Lower baseline (slope)
  #Uint: Upper baseline (intercept)
  #Uslope: Upper baseline (slope)
  #Kval: Quadruplex formation equilibrium constant * 1e10
  #DH: Enthalpy of Quadruplex Formation * 1000 (kJ/mol)
  #Ct #initial concentration of single-stranded oligonucleotide SS1 (Ct) in Mol
  #From solving the model:
  #a = 1
  #b = -4
  #c = 6
  #d = -4 - (1/(4*Kval*Ct^3))
  #e = 1
  #Fitting routines have a hard time fitting very small, or very large numbers.
  #The fitting routine has an easier time with DH input as kJ/mol.
  
  #We need to do solve the roots (alpha) one point at a time:
  for(i in 1:points){
    Kval <- (2/(Ct^3))*exp((-DH*1000/(8.314))*((1/(Temperature[i]+273.15))-(1/(Tm+273.15))))
    a <- 1
    b <- -4
    c <- 6
    d <- -4 - (1/(4*Kval*Ct^3))
    e <- 1

    rawroots <- polynomial (c(e,d,c,b,a)) #polynomial goes UP in order. Modify this if changing order of polynomial.
    if (!all(is.finite(c(a,b,c,d,e)))) {
      stop("Non-finite coefficients at i = ", i)
    }
    foundroots <- solve(rawroots)

    alpha_i <- findReal(foundroots)
    if(i==1){
      results <- alpha_i
    } else {
      results <- c(results,alpha_i)
    }
  }
  results # return the results from the function
}

#Define the chisq function as (Yobs-Ymodel)^2 (we are going to minimize this)
chisq <- function(p) sum((Alpha - Alphafit(p[1],p[2]))^2)

#We can have a look at our baselines here:
plot(Temperature,Abs,main="Absorbance vs. Temperature", xlab="Temperature (°C)", ylab = "OD (260 nm)")
LowerBLfit <- Temperature * Lslope + Lint
UpperBLfit <- Temperature * Uslope + Uint
lines(spline(Temperature, LowerBLfit),col="red")
lines(spline(Temperature, UpperBLfit),col="red")

#Now guess the initial guesses for the parameter estimates on a new plot.
plot(Temperature,Alpha,main="Fraction of Quadruplex vs. Temperature", xlab="Temperature (°C)", ylab = expression(alpha ~ "(fraction of quadruplex)"))
#Play with these initial guesses until the fit gets close.
g1 <- 70    #first guess for Tm
g2 <- -700  #first guess for DH in J/mol

Alphafity <- Alphafit(g1,g2) # This is the model fit using our first guesses.
lines(spline(Temperature, Alphafity)) # plot the model on the current graph.

#Do do the actual fit with the starting values we found. This is done
#by minimizing the function we defined using the starting values we just tried:
#(replace 1, 2 with your best guesses)

fit1 <- nlm(chisq, p = c(g1, g2), hessian = TRUE) # Minimize chisq starting with initial guesses

#To see the results of the fitting routine:
fit1

#Now superimpose our minimized fit on a brand new plot:
# Tm: fit1$estimate[1]
# DH: fit1$estimate[2]
plot(Temperature,Alpha,main="Fraction DNA Quadruplex vs. Temperature", xlab="Temperature (°C)", ylab = expression(alpha ~ "(fraction duplex)"))
Alphafity <- Alphafit(fit1$estimate[1],fit1$estimate[2])
lines(spline(Temperature, Alphafity)) #add the model fit here

#We can go back and look at what the fit would look like if we were to convert it back to absorbance:
plot(Temperature,Abs,main="DNA Absorbance vs. Temperature", xlab="Temperature (°C)", ylab = "OD (260nm)")
#Absfit = alpha*(Upper baseline - Lower baseline) + Lower baseline)
Absfity <- Alphafity*((Uint+Uslope*Temperature)-(Lint+Lslope*Temperature))+(Lint+Lslope*Temperature)
lines(spline(Temperature, Absfity))
#Plot the baselines
LowerBLfit <- Temperature * Lslope + Lint
UpperBLfit <- Temperature * Uslope + Uint
lines(spline(Temperature, Absfity))
lines(spline(Temperature, LowerBLfit),col="red")
lines(spline(Temperature, UpperBLfit),col="red")

#Calculate and report the correlation matrix: 
#(off-diagonal elements > 0.7 --> too many parameters?)
cov.mat <- 2 * fit1$minimum / (length(Abs) - 2) * solve(fit1$hessian)
cor.mat <- cov2cor(cov.mat)
cor.mat

#Calculate and report the approximate standard errors of the parameter estimates:
se <- sqrt(diag(cov.mat))
se

#Calculate the rsq:
sstot = sum((Abs-mean(Abs))^2)
sserr = sum((Abs - Absfity)^2)
rsq = 1-sserr/sstot
rsq
