#Kb n=2 (binding K+).R
#This script fits Kb data assuming SS + 2K+ -> Q (self-complementary)
#Since we have fixed values for the baselines, it also fits Au (unfolded absorbance)
#and Af (folded absorbance), which results in a better fit.
#Written by: David Dubins 
#Date: Jan 18, 2026
#Platform: R-Studio 2026.01.0 Build 392, on R version 4.5.2
#External Libraries used: polynom 1.4 
#(Tools -> Install Packages -> select polynom)
#This program expects two columns of data:
# Ktot   Abs

#Read in the data (Temperature, Abs)
#Uncomment the appropriate line below for your platform.
#Select the data in Excel, press Ctrl+C, then run the appropriate line.
#For PC:
#BindingCurve <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#BindingCurve <- read.table(pipe("pbpaste"), header=TRUE)

attach(BindingCurve)

#Define model parameters
nr <- 3  #number of roots (also the order of the polynomial to fit)
points <- dim(BindingCurve)[1] #number of data points to fit
Ct <- 0.000036  #Initial concentration of SS: 3.6e-5M

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

#Define the function to fit (quadruplex melting to SS1, self-complementary)
Alphafit <- function(Kb,Au,Af) {
  #Au: Lower baseline (unfolded state)
  #Af: Upper baseline (folded state)
  #Kb: Quadruplex formation binding constant * 1000 (helps out fitting routine)
  #Ct #initial concentration of single-stranded oligonucleotide SS1 (Ct) in Mol
  #From solving the model:
  #a = -4Ct^2
  #b = 4Ct^2 + 4 Ct*Ktot
  #c = -4Ct*Ktot - Ktot^2 - 1/Kb
  #d = Ktot^2
  #Fitting routines have a hard time fitting very small, or very large numbers.
  #The fitting routine has an easier time with DH input as kJ/mol.
  
  #We need to do solve the roots (alpha) one point at a time:
  for(i in 1:points){
    a <- -4*(Ct^2)
    b <- 4*(Ct^2) + 4*Ct*Ktot[i]
    c <- -4*Ct*Ktot[i] - (Ktot[i]^2) - (1/(Kb*1000))
    d <- Ktot[i]^2

    rawroots <- polynomial (c(d,c,b,a)) #polynomial goes UP in order. Modify this if changing order of polynomial.
    if (!all(is.finite(c(a,b,c,d)))) {
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
  #Now calculate OD (we are fitting the baselines as well, so the function needs to return Absfit)
  #Absfit = Alpha * (Upper baseline - Lower Baseline) + (Lower baseline)
  Absfit <- (results * (Af - Au)) + Au
  Absfit # return the results from the function
}

#Define the chisq function as (Yobs-Ymodel)^2 (we are going to minimize this)
chisq <- function(p) sum((Abs - Alphafit(p[1],p[2],p[3]))^2)

#We can have a look at our baselines here:
plot(Ktot,Abs,main="Elipticity vs. [K+]", xlab="[K+] (M)", ylab = "Elipticity (268 nm)")
LowerBLfit <- Ktot*0 + Au
UpperBLfit <- Ktot*0 + Af
lines(spline(Ktot, LowerBLfit))
lines(spline(Ktot, UpperBLfit))

#Now guess the initial guesses for the parameter estimates on a new plot.
#plot(Ktot,Abs,main="Fraction of Quadruplex vs. Potassium Ion Concentration", xlab="[K+]tot (M)", ylab = expression(alpha ~ "(fraction of quadruplex)"))

plot(Ktot,Abs,main="Elipticity vs. [K+]", xlab="[K+] (M)", ylab = "Elipticity (268 nm)")
#Play with these initial guesses until the fit gets close.
g1 <- 1.5  #first guess for Kb
g2 <- 910.15194816877      #Baseline for for Au (unfolded state baseline)
g3 <- 6100.3962325024      #Baseline for Af (folded state baseline)

Absfity <- Alphafit(g1,g2,g3) # This is the model fit using our first guesses.
lines(spline(Ktot, Absfity)) # plot the model on the current graph.
LowerBLfit <- Ktot*0 + g2    # plot the lower baseline
UpperBLfit <- Ktot*0 + g3    # plot the upper baseline
lines(spline(Ktot, LowerBLfit),col="red")
lines(spline(Ktot, UpperBLfit),col="red")

#Do do the actual fit with the starting values we found. This is done
#by minimizing the function we defined using the starting values we just tried:
#(replace g1 with your best guesses)

fit1 <- nlm(chisq, p = c(g1,g2,g3), hessian = TRUE) # Minimize chisq starting with initial guesses

#To see the results of the fitting routine:
fit1

#Now superimpose our minimized fit on a brand new plot:
# Kb: fit1$estimate[1]
# Au: fit1$estimate[2]
# Af: fit1$estimate[3]
plot(Ktot,Abs,main="Elipticity vs. [K+]", xlab="[K+] (M)", ylab = "Elipticity (268 nm)")
Absfity <- Alphafit(fit1$estimate[1],fit1$estimate[2],fit1$estimate[3])
lines(spline(Ktot, Absfity)) #add the model fit here
LowerBLfit <- Ktot*0 + fit1$estimate[2]    # plot the lower baseline
UpperBLfit <- Ktot*0 + fit1$estimate[3]    # plot the upper baseline
lines(spline(Ktot, LowerBLfit),col="red")
lines(spline(Ktot, UpperBLfit),col="red")

#Calculate Alpha vs. [K+] curve
#Alpha = (Abs - Lower Baseline) / (Upper baseline - Lower baseline)
Alpha <- (Abs-fit1$estimate[2])/(fit1$estimate[3]-fit1$estimate[2])
plot(Ktot,Alpha,main="Fraction DNA Quadruplex vs. Potassium Ion Concentration", xlab="[K+](M)", ylab = expression(alpha ~ "(fraction duplex)"))
Alphafity <- (Absfity-fit1$estimate[2])/(fit1$estimate[3]-fit1$estimate[2])
lines(spline(Ktot, Alphafity)) #add the model fit here

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
