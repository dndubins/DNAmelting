#Tm n=1 (self-complementary).R
#This script fits TM data assuming mono-molecular melting
#Written by: David Dubins
#Date: Jan 17, 2026
#Platform: R-Studio 2026.01.0 Build 392, on R version 4.5.2
#External Libraries used: None

#Read in the data (Temperature, OD)
#For PC:
MeltingCurve <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#MeltingCurve <- read.table(pipe("pbpaste"), header=TRUE)

attach(MeltingCurve)

#Here are the lower and upper baselines for the spectroscopic data:
Lint <- 1.438055444     #first guess for Lint
Lslope <- 0.000979435   #first guess for Lslope
Uint <- 1.58965671      #first guess for Uint
Uslope <- 0.001000519   #first guess for Uslope

#Define the function to fit (monomolecular melting)
Absfit <- function(Tm,DH) {
#Alphafit:  (Exp((-XP(6, 1) / 8.314) * (1 / (amtx(i, 1) + 273.15) - (1 / (XP(1, 1) + 273.15))))) / (Exp((-XP(6, 1) / 8.314) * (1 / (amtx(i, 1) + 273.15) - (1 / (XP(1, 1) + 273.15)))) + 1)
  #This is the function of alpha vs. Temperature
  alpha <- 1 / (1+exp((DH/8.314)*((1/(Temperature+273.15))-(1/(Tm+273.15)))))
  #This converts alpha back to absorbance for the fit:
  (Uint+Uslope*Temperature)-alpha*((Uint+Uslope*Temperature)-(Lint+Lslope*Temperature))
}

#Define the chisq function as (Yobs-Ymodel)^2 (we are going to minimize this)

chisq <- function(p) sum((Abs - Absfit(p[1],p[2]))^2)

#Now guess some values on a new plot.
#Set up variable names g1 to g2 --> easier to change 1st guesses.
#(replace 1, 2 with your starting guesses)

plot(Temperature,Abs,main="Absorbance vs. Temperature", xlab="Temperature (?C)", ylab = "Absorbance")
g1 <- 48   #first guess for Tm
g2 <- -607836   #first guess for DH

Absfitx <-Temperature
Absfity <- Absfit(g1,g2)
LowerBLfit <- Temperature * Lslope + Lint
UpperBLfit <- Temperature * Uslope + Uint

lines(spline(Absfitx, Absfity))
lines(spline(Absfitx, LowerBLfit))
lines(spline(Absfitx, UpperBLfit))

#Do do the actual fit with the starting values we found. This is done
#by minimizing the function we defined using the starting values we just tried:
#(replace 1, 2, and 3 with your best guesses)

fit1 <- nlm(chisq, p = c(g1, g2), hessian = TRUE)

#To see the results of the fitting routine:
fit1

#Now superimpose our minimized fit on a brand new plot:
plot(Temperature,Abs,main="DNA Absorbance vs. Temperature", xlab="Temperature (°C)", ylab = "OD (260nm)")
Absfitx <- MeltingCurve$Temperature
Absfity <- Absfit(fit1$estimate[1],fit1$estimate[2]) # This is our fitted model
lines(spline(Absfitx, Absfity))
#Plot the baselines
LowerBLfit <- Temperature * Lslope + Lint
UpperBLfit <- Temperature * Uslope + Uint

lines(spline(Absfitx, Absfity))
lines(spline(Absfitx, LowerBLfit))
lines(spline(Absfitx, UpperBLfit))

#Let's make an alpha vs. T plot for fun. Our solved estimates are:
# Tm: fit1$estimate[1]
# DH: fit1$estimate[2]
Alpha <- ((Uint+Uslope*Temperature)-Abs)/((Uint+Uslope*Temperature)-(Lint+Lslope*Temperature))
Alphafit <- 1 / (1+exp((fit1$estimate[2]/8.314)*((1/(Temperature+273.15))-(1/(fit1$estimate[1]+273.15)))))
plot(Temperature,Alpha,main="Fraction DNA Duplex vs. Temperature", xlab="Temperature (°C)", ylab = expression(alpha ~ "(fraction duplex)"))
lines(spline(Temperature, Alphafit)) #add the model fit here

#To obtain the approximate standard errors of the parameter estimates:
sqrt(diag(2*fit1$minimum/(length(Abs)-2)*solve(fit1$hessian)))

#Calculate the rsq:
sstot = sum((Abs-mean(Abs))^2)
sserr = sum((Abs - Absfity)^2)
rsq = 1-sserr/sstot
rsq

