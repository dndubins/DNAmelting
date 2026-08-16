#Tm n=2 (not self-complementary), with kinetic fit of cooling curve.R
#This script fits TM data assuming D -> SS1 + SS2 (not self-complementary)
#Original heating-curve section written by: David Dubins
#Kinetic cooling-curve section added: fits kon (koff derived from kon/Keq(T), Keq(T) fixed by the heating-curve fit)
#Date: Jan 17, 2026
#Last Updated: 14-Aug-26
#Platform: R-Studio 2026.01.0 Build 392, on R version 4.5.2
#Kinetics solved by: Claude.AI Sonnet 5
#External Libraries used: polynom 1.4, deSolve
#(Tools -> Install Packages -> select polynom, deSolve)
#This program expects three columns of data:
# Time  Temperature   Abs
# 0 ,    70,         1.39034 (or whatever)

# DNA melting analysis
# --------------------
# 1. Import data
# 2. Fit baselines
# 3. Calculate alpha
# 4. Fit equilibrium parameters
# 5. Fit cooling kinetics
# 6. Calculate dependent parameters
# 7. Plot results

#Read in the data (Temperature, Abs)
#Uncomment the appropriate line below for your platform.
#Select the data in Excel, press Ctrl+C, then run the appropriate line.

#Import the Heating Curve:
#=========================
#For PC:
HC <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#HC <- read.table(pipe("pbpaste"), header=TRUE)

#Import the Cooling Curve:
#=========================
#For PC:
CC <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#CC <- read.table(pipe("pbpaste"), header=TRUE)

#Heating Curve: Bring your own baselines
#Enter the lower and upper baselines for the heating curve data (Fit using MS Excel):
#HC_Lint <- 14759.2277210086   #first guess for Lint
#HC_Lslope <- -107.701947979195   #first guess for Lslope
#HC_Uint <- 5071.67248566532    #first guess for Uint
#HC_Uslope <- -56      #first guess for Uslope
#Cooling Curve: Bring your own baselines
#CC_Lint <- 13185.2479462185   #first guess for Lint
#CC_Lslope <- -72.8998223289307   #first guess for Lslope
#CC_Uint <- 2272.18116647355    #first guess for Uint
#CC_Uslope <- -16      #first guess for Uslope

#Fitting the lower and upper baselines in R-project:

# Heating Curve Baselines first:
#Plot the data to decide start and end temperatures for baseline fitting.
plot(HC$Temperature,HC$Abs,main="Absorbance vs. Temperature", xlab="Temperature (°C)", ylab = "OD", col="red")

Lstart <- 42 # start temperature of lower baseline fit
Lstop <- 50  # stop temperature of lower baseline fit
Ustart <- 66 # start temperature of upper baseline fit
Ustop <- 70  # stop temperature of upper baseline fit

#Isolate and fit the lower baseline:
HC_lowerBL <- HC[(HC[, "Temperature"] >= Lstart & HC[, "Temperature"] <= Lstop), ]
BLfit <- lm(HC_lowerBL$Abs ~ HC_lowerBL$Temperature)
HC_Lint <- unname(coef(BLfit)[1])   #fitted 
HC_Lslope <- unname(coef(BLfit)[2])   #first guess for Lslope

#Isolate and fit the upper baseline:
HC_upperBL <- HC[(HC[, "Temperature"] >= Ustart & HC[, "Temperature"] <= Ustop), ]
BLfit <- lm(HC_upperBL$Abs ~ HC_upperBL$Temperature)
HC_Uint <- unname(coef(BLfit)[1])   #fitted 
HC_Uslope <- unname(coef(BLfit)[2])   #first guess for Lslope

#Calculate and display the heating curve baselines:
HC_LowerBLfit <- HC$Temperature * HC_Lslope + HC_Lint
HC_UpperBLfit <- HC$Temperature * HC_Uslope + HC_Uint
lines(spline(HC$Temperature, HC_LowerBLfit), col="red")
lines(spline(HC$Temperature, HC_UpperBLfit), col="red")
#Add cooling curve to plot: (comment out or don't run if it's distracting)
points(CC$Temperature,CC$Abs, col="blue")

###### END OF HEATING CURVE BASELINE FITTING HERE ######

# Cooling Curve Baselines next:
#Plot the data to decide start and end temperatures for baseline fitting.
plot(CC$Temperature,CC$Abs,main="Absorbance vs. Temperature", xlab="Temperature (°C)", ylab = "OD", col="blue")

Lstart <- 35 # start temperature of lower baseline fit
Lstop <- 40  # stop temperature of lower baseline fit
Ustart <- 68 # start temperature of upper baseline fit
Ustop <- 70  # stop temperature of upper baseline fit

#Isolate and fit the lower baseline:
CC_lowerBL <- CC[(CC[, "Temperature"] >= Lstart & CC[, "Temperature"] <= Lstop), ]
BLfit <- lm(CC_lowerBL$Abs ~ CC_lowerBL$Temperature)
CC_Lint <- unname(coef(BLfit)[1])   #fitted 
CC_Lslope <- unname(coef(BLfit)[2])   #first guess for Lslope

#Isolate and fit the upper baseline:
CC_upperBL <- CC[(CC[, "Temperature"] >= Ustart & CC[, "Temperature"] <= Ustop), ]
BLfit <- lm(CC_upperBL$Abs ~ CC_upperBL$Temperature)
CC_Uint <- unname(coef(BLfit)[1])   #fitted 
CC_Uslope <- unname(coef(BLfit)[2])   #first guess for Lslope

#Calculate and display the cooling curve baselines:
CC_LowerBLfit <- CC$Temperature * CC_Lslope + CC_Lint
CC_UpperBLfit <- CC$Temperature * CC_Uslope + CC_Uint
lines(spline(CC$Temperature, CC_LowerBLfit), col="blue")
lines(spline(CC$Temperature, CC_UpperBLfit), col="blue")
#Add heating curve to plot: (comment out or don't run if it's distracting)
points(HC$Temperature,HC$Abs, col="red")

###### END OF COOLING CURVE BASELINE FITTING HERE ######

#NOTE: FINALIZE THE BASELINES before proceeding.
#Visually examine the baselines. Are they nice? If not, change limits above.
#Also considering trimming out garbage data at temperature extremes.
#This plot needs to be pretty!

#Now calculate Alpha vs. Temperature curve
#Alpha = (Upper baseline - Alpha) / (upper baseline - lower baseline)
HC_Alpha <- ((HC_Uint+HC_Uslope*HC$Temperature)-HC$Abs)/((HC_Uint+HC_Uslope*HC$Temperature)-(HC_Lint+HC_Lslope*HC$Temperature))

#Define model parameters
nr <- 2  #number of roots (also the order of the polynomial to fit)
points <- dim(HC)[1] #number of data points to fit
Ct <- 2.0e-6  #2 uM in molar
R <- 8.314    #gas constant, J mol^-1 K^-1

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

#Generalized equilibrium-alpha calculator. Takes any Temp vector instead of being hardwired 
#to HC$Temperature - used below for the ODE's initial condition, and to overlay the equilibrium
#curve on the cooling-curve plot):

#Define the function to fit (duplex melting to SS1 + SS2, non self-complementary)
EqAlphafit <- function(Temp, Tm, DH) {
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
    Kval <- (4/Ct)*exp(((-DH*1000)/(R*(Temp[i]+273.15)))*(1-(Temp[i]+273.15)/(Tm+273.15)))
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
HC_chisq <- function(p) sum((HC_Alpha - EqAlphafit(HC$Temperature,p[1],p[2]))^2)

#Now plot normalized curve (Alpha vs Temperature) to get initial guesses close.
plot(
  HC$Temperature,
  HC_Alpha,main="Fraction of Duplex vs. Temperature", 
  xlab="Temperature (°C)", 
  ylab = expression(alpha ~ "(fraction of duplex)"),
  col="red"
)

#Play with these initial guesses until the fit gets close.
g1 <- 59.0  #first guess for Tm.
g2 <- -510  #first guess for DH in kJ/mol

HC_fittedAlpha <- EqAlphafit(HC$Temperature,g1,g2) # This is the model fit using our first guesses.
lines(spline(HC$Temperature, HC_fittedAlpha),col="red") # plot the model on the current graph.

#Do do the actual fit with the starting values we found. This is done
#by minimizing the function we defined using the starting values we just tried:
HC_fit1 <- nlm(HC_chisq, p = c(g1, g2), hessian = TRUE) # Minimize chisq starting with initial guesses

#To see the results of the fitting routine:
HC_fit1

#Now superimpose our minimized fit on a brand new plot:
# Tm: HC_fit1$estimate[1]
# DH: HC_fit1$estimate[2]
plot(
  HC$Temperature,
  HC_Alpha,main="Fraction DNA Duplex vs. Temperature", 
  xlab="Temperature (°C)", 
  ylab = expression(alpha ~ "(fraction duplex)"),
  col="red"
)
#Calculate then plot fitted heating curve alphas:
HC_fittedAlpha <- EqAlphafit(HC$Temperature,HC_fit1$estimate[1],HC_fit1$estimate[2])
lines(spline(HC$Temperature, HC_fittedAlpha),col="red") #add the model fit here

#(Not needed:)
#We can go back and look at what the fit would look like if we were to convert it back to absorbance:
plot(HC$Temperature,HC$Abs,main="DNA Absorbance vs. Temperature", xlab="Temperature (°C)", ylab = "OD (260nm)", col="red")
#Absfit = Upper baseline - alpha*(Upper baseline - Lower baseline)
HC_Absfity <- (HC_Uint+HC_Uslope*HC$Temperature) - HC_fittedAlpha*((HC_Uint+HC_Uslope*HC$Temperature)-(HC_Lint+HC_Lslope*HC$Temperature))
lines(spline(HC$Temperature, HC_Absfity),col="red")
#Plot the baselines
HC_LowerBLfit <- HC$Temperature * HC_Lslope + HC_Lint
HC_UpperBLfit <- HC$Temperature * HC_Uslope + HC_Uint
lines(spline(HC$Temperature, HC_Absfity),col="red")
lines(spline(HC$Temperature, HC_LowerBLfit),col="red")
lines(spline(HC$Temperature, HC_UpperBLfit),col="red")

#Calculate and report the correlation matrix: 
#(off-diagonal elements > 0.7 --> too many parameters?)
HC_cov.mat <- 2 * HC_fit1$minimum / (length(HC_Alpha) - 2) * solve(HC_fit1$hessian)
HC_cor.mat <- cov2cor(HC_cov.mat)
HC_cor.mat

#Calcualte and report the approximate standard errors of the parameter estimates:
HC_se <- sqrt(diag(HC_cov.mat))
HC_se

#Calculate the rsq:
HC_sstot = sum((HC_Alpha - mean(HC_Alpha))^2)
HC_sserr = sum((HC_Alpha - HC_fittedAlpha)^2)
HC_rsq = 1-HC_sserr/HC_sstot
HC_rsq

#############################################
#           HALF WAY POINT!!!!              #
#############################################
library(deSolve)

#First calculate Alpha vs. Temperature curve (cooling curve)
#Alpha = (Upper baseline - Alpha) / (upper baseline - lower baseline)
CC_Alpha <- ((CC_Uint+CC_Uslope*CC$Temperature)-CC$Abs)/((CC_Uint+CC_Uslope*CC$Temperature)-(CC_Lint+CC_Lslope*CC$Temperature))

#Now plot normalized alpha vs. Temp for both curves:
plot(
  CC$Temperature,
  CC_Alpha,
  main="Fraction DNA Duplex vs. Temperature",
  xlab="Temperature (°C)",
  ylab = expression(alpha ~ "(fraction duplex)"),
  col="blue"
)
#plot the normalized heating curve to see hysteresis:
points(HC$Temperature,HC_Alpha, col="red")

#Time to prepare the ODE fit!
#First define some useful constants from our heating curve fit:
Tm <- HC_fit1$estimate[1]    #Tm from our first fit
DH <- HC_fit1$estimate[2]    #DH in kJ/mol from our first fit
#Ct <- 2.0e-6  #2 uM in molar (already defined above)
#R <- 8.314    #gas constant, J mol^-1 K^-1 (already defined above)
Tref <- 56.79 + 273.15 # choose a reference temperature not off in space (better for fitting)

# Make Keq a function:
Keq <- function(Temp) {
  T <- Temp + 273.15
  Tm_K <- Tm + 273.15
  (4/Ct) * exp(((-DH*1000)/(R*T)) * (1 - T/Tm_K))
}

# Shift the time data to the first point being time=0: (might as well reset the clock)
CC$Time <- CC$Time - CC$Time[1]

#Create a function for the cooling curve temperature vs. time data:
#Interpolate the data using R's very convenient interpolation function:
Temp <- function(time) { # this function won't crap out outside time limits.
  time <- pmin(pmax(time, min(CC$Time)), max(CC$Time))
  approx(CC$Time, CC$Temperature, xout = time)$y
}

#Sanity check:
Temp(0) # This should be the highest (first) temperature
Temp(CC$Time[length(CC$Time)]) # This should be the lowest (last) temperature
Temp(5461) # one second after last data point

# This will look "perfect" because the cooling curve was forced to be linear,
# but with actual observed data, it won't (and doesn't need to) be a straight line.
plot(
  CC$Time,
  Temp(CC$Time),
  main="Temperature vs. Time (Cooling Curve)",
  xlab="Time (s)",
  ylab ="Temperature (°C)",
  col="blue"
)

#Write a function to calculate ln(koff) as a function of temperature.
#This is our Arrhenius relationship for the rate constant.
ln_koff <- function(Temp,ln_kref,Ea_off) {
  T <- Temp + 273.15
  ln_kref - (Ea_off*1000/R)*((1/T) - (1/Tref))
}

#Now calculate our starting point at the starting temperature (70 degC?):
Tstart <- CC$Temperature[1] # get starting temperature (=CC$Temperature[1])
Kt0 <- Keq(Tstart) # calculate starting Keq
a <- 1 # funky coefficients that are dependent on rxn mechanism.
b <- -2 - 2/(Kt0*Ct)
c <- 1
roots <- polyroot(c(c, b, a))
roots

#Now find the real root:
alpha0 <- findReal(roots)
alpha0 # This needs to be a number between 0 and 1 for this fit (it's the first alpha).

#If there is some loose imaginary change left over, get rid of it (visual inspection):
#alpha0 <- Re(roots[1]) # Take the real component of the root that's between 0 and 1.

#Here we define our ordinary diffeq for the deSolve routine to use:
DNA_ODE <- function(time, state, parameters) {
  alpha <- state[1] # this means alpha=alpha(t).
  T <- Temp(time) # This means at whatever time you are currently evaluating, go find the corresponding temperature.
  ln_kref <- parameters["ln_kref"] # ln_kref is a parameter
  Ea_off <- parameters["Ea_off"] # Ea_off is a parameter
  ln_koff_T <- ln_koff(T, ln_kref, Ea_off) # calculate koff_t at this temperature
  dalpha_dt <- exp(ln_koff_T) * (Keq(T)*Ct/2*(1-alpha)^2 - alpha) # This is our model
  list(c(dalpha_dt)) # this is how we hand the derivative back to deSolve.
}

#DNA_ODE(0, c(alpha0), parms) # calculate the ODE at time zero. (System is at eq'm, should be zero)

#This is the structure of simulating a curve:
#create a series of times (0 to 900 seconds)
#sampletimes <- seq(0, 5460, by = 180)
#out <- ode(
#  y = c(alpha = alpha0), # gives the solver our initial condition
#  times = sampletimes, # our 15-minute cooling experiment, sampled every second (to 900s)
#  func = DNA_ODE, # gives the solver the differential equation
#  parms = parms # using parameters "parms" defined above
#)
#head(out)
#tail(out)

# Make cooling curve fitting function:
CC_fitAlpha <- function(p1,p2){
  pvector <- c(ln_kref = unname(p1),Ea_off = unname(p2)) # build a parameter vector
  out <- ode(
    y = c(alpha = alpha0), # gives the solver our initial condition
    times = CC$Time, # The times from our cooling curve
    func = DNA_ODE, # gives the solver the differential equation
    parms = pvector # using parameters "parms" defined above
  )
  out[,2] # return only the fitted y-values
}

# rough starting guess for ln_kref from your old fit, evaluated at Tref:
g1 <- -7.9 # first guess for ln_kref
g2 <- 250 # first guess for Ea_off

# Here is our parameter vector, with ln_kref and Ea_off defined starting guesses:
parms <- c(
  ln_kref = g1,   # trial pre-exponential factor, s^-1
  Ea_off = g2       # trial activation energy, kJ/mol
)

#Look at the initial CC_fit1 to adjust ln_k0_off and Ea_off start values:
CC_start <- CC_fitAlpha(g1, g2)
#Now plot normalized alpha vs. Temp for both curves:
plot(
  CC$Temperature,
  CC_Alpha,
  main="Fraction DNA Duplex vs. Temperature",
  xlab="Temperature (°C)",
  ylab = expression(alpha ~ "(fraction duplex)"),
  col="blue"
)
#plot the first guess:
lines(spline(CC$Temperature, CC_start),col="blue")

#####Go back and change the first guesses until you get a reasonable fit by eye!
# Note: This fit is very sensitive to starting conditions.

#Define the chisq function as (Yobs-Ymodel)^2 (we are going to minimize this)
chisq <- function(p) sum((CC_Alpha - CC_fitAlpha(p[1],p[2]))^2)

#Perform the nlm fit:
CC_fit1 <- nlm(chisq, p = parms, hessian = TRUE) # Minimize chisq starting with initial guesses
#To see the results of the fitting routine:
CC_fit1

#We can have a look at our fit here:
#Now plot normalized alpha vs. Temp for both curves:
plot(
  CC$Temperature,
  CC_Alpha,
  main="Fraction DNA Duplex vs. Temperature",
  xlab="Temperature (°C)",
  ylab = expression(alpha ~ "(fraction duplex)"),
  col="blue"
)
#plot the normalized heating curve to see hysteresis:
points(HC$Temperature,HC_Alpha, col="red")
#add our heating curve fit:
lines(spline(HC$Temperature, HC_fittedAlpha),col="red") #add the model fit here
#calculate the fitted cooling curve:
CC_fittedAlpha <- CC_fitAlpha(CC_fit1$estimate[1], CC_fit1$estimate[2])
#add our cooling curve fit:
lines(spline(CC$Temperature, CC_fittedAlpha),col="blue")

#Calculate and report the correlation matrix: 
#(off-diagonal elements > 0.7 --> too many parameters?)
CC_cov.mat <- 2 * CC_fit1$minimum / (length(CC_Alpha) - 2) * solve(CC_fit1$hessian)
CC_cor.mat <- cov2cor(CC_cov.mat)
CC_cor.mat

#Calcualte and report the approximate standard errors of the parameter estimates:
CC_se <- sqrt(diag(CC_cov.mat))
CC_se

#Calculate the rsq:
CC_sstot = sum((CC_Alpha - mean(CC_Alpha))^2)
CC_sserr = sum((CC_Alpha - CC_fittedAlpha)^2)
CC_rsq = 1-CC_sserr/CC_sstot
CC_rsq

#Now calculate and report all remaining dependent parameters:

#Write a function to calculate ln(kon) as a function of temperature.
#This comes from our equilibrium relationship.
ln_kon <- function(Temp,ln_k0_off,Ea_off) {
  T <- Temp + 273.15
  ln_k_off <- (-Ea_off*1000)/(R*T) + ln_k0_off
  k_off=exp(ln_k_off)
  k_on=Keq(Temp)*k_off
  log(k_on) # return ln(k_on)
}

#Calculate and Report Ea_on:
Ea_off = CC_fit1$estimate[2]
Ea_off

#Get ln_kref:
ln_kref <- CC_fit1$estimate[1]
ln_k0_off <- ln_kref + Ea_off*1000/(R*Tref) # This is the true intercept of ln(koff) vs. 1/T

#Calculate and Report ln(k0_on) and Ea_on:
xvals <- 1.0/(CC$Temperature + 273.15)
yvals <- ln_kon(CC$Temperature,ln_k0_off,Ea_off)
Arrhenfit <- lm(yvals~xvals)
Ea_on <- -coef(Arrhenfit)[2]*R / 1000 # This should be Ea_on in kJ/mol
Ea_on
ln_k0_on <- coef(Arrhenfit)[1] # This should be ln_k0_on
ln_k0_on

#Sanity Check: DH = E_on - E_off
DH
Ea_on-Ea_off

# ============================================================
# Summary table of fitting results
# ============================================================

Fit_Summary <- data.frame(
  Parameter = c(
    "Tm",
    "DH",
    "ln(kref)",
    "ln(koff0)",
    "Ea_off",
    "ln(kon0)",
    "Ea_on"
  ),
  
  Value = c(
    Tm,
    DH,
    CC_fit1$estimate[1],
    ln_k0_off,
    CC_fit1$estimate[2],
    ln_k0_on,
    Ea_on
  ),
  
  Units = c(
    "°C",
    "kJ/mol",
    "dimensionless",
    "dimensionless",
    "kJ/mol",
    "dimensionless",
    "kJ/mol"
  )
)

# Format the numbers nicely
Fit_Summary$Value <- round(Fit_Summary$Value, 3)

# Print the table
Fit_Summary

##################################################################################

# Finding the best Tref:
Tref_scan <- seq(Tm - 15, Tm + 15, by = 2) + 273.15
r_scan <- sapply(Tref_scan, function(Tr) {
  Tref <<- Tr   # ln_koff() reads Tref from the global env
  fit <- nlm(chisq, p = parms, hessian = TRUE)
  cov.mat <- 2 * fit$minimum / (length(CC_Alpha) - 2) * solve(fit$hessian)
  cov2cor(cov.mat)[1,2]
})
plot(Tref_scan - 273.15, r_scan, xlab = "Tref (°C)", ylab = "correlation(ln_kref, Ea_off)")

#Interpolate the data using R's very convenient interpolation function:
y <- approxfun(Tref_scan - 273.15, r_scan)
# Find out what temperature gets you a y of zero. Then make that Tref!
y(43) 
y(57)
y(56.79)