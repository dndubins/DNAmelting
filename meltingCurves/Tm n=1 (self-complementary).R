#Tm n=1 (self-complementary).R
#This script fits TM data assuming D -> SS1 (monomolecular)
#Written by: David Dubins 
#Date: Jan 17, 2026
#Platform: R-Studio 2026.01.0 Build 392, on R version 4.5.2
#External Libraries used: None
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
Lint <- 1.438689416     #lower baseline intercept
Lslope <- 0.000958356   #lower baseline slope
Uint <- 1.600283171     #upper baseline intercept
Uslope <- 0.000839927   #upper baseline slope

#First calculate Alpha vs. Temperature curve
#Alpha = (Upper baseline - Abs) / (Upper baseline - Lower baseline)
Alpha <- ((Uint+Uslope*Temperature)-Abs)/((Uint+Uslope*Temperature)-(Lint+Lslope*Temperature))

#Define the function to fit (self-complementary duplex melting)
Alphafit <- function(Tm,DH) {
  #Fitting routines have a hard time fitting very small, or very large numbers.
  #The fitting routine has an easier time with DH input as kJ/mol.
  #This is the function of alpha vs. Temperature
  #Alphafit = 1 / (1 + exp((DH/R) * (1 / (T+273.15) - (1 / Tm + 273.15))))
  1 / (1+exp((DH*1000/8.314)*((1/(Temperature+273.15))-(1/(Tm+273.15)))))
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
plot(Temperature,Alpha,main="Fraction of Duplex vs. Temperature", xlab="Temperature (°C)", ylab = expression(alpha ~ "(fraction of duplex)"))
#Play with these initial guesses until the fit gets close.
g1 <- 48    #first guess for Tm.
g2 <- -500  #first guess for DH in kJ/mol

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
