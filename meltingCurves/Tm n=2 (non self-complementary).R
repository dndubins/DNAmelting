#Tm n=2 (not self-complementary).R
#This script fits TM data assuming D -> SS1 + SS2 (not self-complementary)
#Written by: David Dubins 
#Date: Jan 17, 2026
#Platform: R-Studio 2026.01.0 Build 392, on R version 4.5.2
#External Libraries used: polynom 1.4 
#(Tools -> Install Packages -> select polynom)
#This program expects two columns of data:
# Temperature   Abs

#Read in the data (Temperature, Abs)
#Uncomment the appropriate line below for your platform.
#Select the data in Excel, press Ctrl+C, then run the appropriate line.
#For PC:
#MeltingCurve <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#MeltingCurve <- read.table(pipe("pbpaste"), header=TRUE)

attach(MeltingCurve)

#Enter the lower and upper baselines for the spectroscopic data (Fit using MS Excel):
Lint <- 1.438055444     #first guess for Lint
Lslope <- 0.000979435   #first guess for Lslope
Uint <- 1.58965671      #first guess for Uint
Uslope <- 0.001000519   #first guess for Uslope

#First calculate Alpha vs. Temperature curve
#Alpha = (Upper baseline - Alpha) / (upper baseline - lower baseline)
Alpha <- ((Uint+Uslope*Temperature)-Abs)/((Uint+Uslope*Temperature)-(Lint+Lslope*Temperature))

#Define model parameters
nr <- 2  #number of roots (also the order of the polynomial to fit)
points <- dim(MeltingCurve)[1] #number of data points to fit
Ct <- 4.0e-6  #4 uM in molar

library(polynom) #install the polynom package before loading for the first time

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

#Define the function to fit (duplex melting to SS1 + SS2, non self-complementary)
Alphafit <- function(Tm,DH) {
  #Lint: Lower baseline (intercept)
  #Lslope: Lower baseline (slope)
  #Uint: Upper baseline (intercept)
  #Uslope: Upper baseline (slope)
  #Kval: Duplex formation equilibrium constant
  #Ct #initial concentration of single-stranded oligonucleotide SS1 (Ct) in Mol (assume SS1=SS2)
  #From solving the model:
  #a = 1
  #b = -2 - (2/Kval*Ct)
  #c = 1
  #Fitting routines have a hard time fitting very small, or very large numbers.
  #The fitting routine has an easier time with DH input as kJ/mol.
  
  #We need to do solve the roots (alpha) one point at a time:
  for(i in 1:points){
    Kval <- (4/Ct)*exp(((-DH*1000)/(8.314*(Temperature[i]+273.15)))*(1-(Temperature[i]+273.15)/(Tm+273.15)))
    a <- 1
    b <- -2 - (2/(Kval*Ct))
    c <- 1

    rawroots <- polynomial (c(c,b,a)) #polynomial goes UP in order. Modify this if changing order of polynomial.
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
lines(spline(Temperature, LowerBLfit))
lines(spline(Temperature, UpperBLfit))

#Now guess the initial guesses for the parameter estimates on a new plot.
plot(Temperature,Alpha,main="Fraction of Duplex vs. Temperature", xlab="Temperature (°C)", ylab = expression(alpha ~ "(fraction of duplex)"))
#Play with these initial guesses until the fit gets close.
g1 <- 54    #first guess for Tm.
g2 <- -1000  #first guess for DH in kJ/mol

Alphafity <- Alphafit(g1,g2) # This is the model fit using our first guesses.
lines(spline(Temperature, Alphafity)) # plot the model on the current graph.

#Do do the actual fit with the starting values we found. This is done
#by minimizing the function we defined using the starting values we just tried:
#(replace 1, 2, and 3 with your best guesses)

fit1 <- nlm(chisq, p = c(g1, g2), hessian = TRUE) # Minimize chisq starting with initial guesses

#To see the results of the fitting routine:
fit1

#Now superimpose our minimized fit on a brand new plot:
# Tm: fit1$estimate[1]
# DH: fit1$estimate[2]
plot(Temperature,Alpha,main="Fraction DNA Duplex vs. Temperature", xlab="Temperature (°C)", ylab = expression(alpha ~ "(fraction duplex)"))
Alphafity <- Alphafit(fit1$estimate[1],fit1$estimate[2])
lines(spline(Temperature, Alphafity)) #add the model fit here

#We can go back and look at what the fit would look like if we were to convert it back to absorbance:
plot(Temperature,Abs,main="DNA Absorbance vs. Temperature", xlab="Temperature (°C)", ylab = "OD (260nm)")
#Absfit = Upper baseline - alpha*(Upper baseline - Lower baseline)
Absfity <- (Uint+Uslope*Temperature) - Alphafity*((Uint+Uslope*Temperature)-(Lint+Lslope*Temperature))
lines(spline(Temperature, Absfity))
#Plot the baselines
LowerBLfit <- Temperature * Lslope + Lint
UpperBLfit <- Temperature * Uslope + Uint
lines(spline(Temperature, Absfity))
lines(spline(Temperature, LowerBLfit))
lines(spline(Temperature, UpperBLfit))

#Calculate and report the correlation matrix: 
#(off-diagonal elements > 0.7 --> too many parameters?)
cov.mat <- 2 * fit1$minimum / (length(Abs) - 2) * solve(fit1$hessian)
cor.mat <- cov2cor(cov.mat)
cor.mat

#Calcualte and report the approximate standard errors of the parameter estimates:
se <- sqrt(diag(cov.mat))
se

#Calculate the rsq:
sstot = sum((Abs-mean(Abs))^2)
sserr = sum((Abs - Absfity)^2)
rsq = 1-sserr/sstot
rsq
