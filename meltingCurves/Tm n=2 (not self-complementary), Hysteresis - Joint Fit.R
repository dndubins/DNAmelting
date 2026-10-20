#Tm n=2 (not self-complementary), Hysteresis - Joint Fit.R
#This script fits TM data assuming D -> SS1 + SS2 (not self-complementary)
#This script was written with assistance from Claude.AI and Perplexity.AI
#The bottom of the script allows you to simulate different heating and cooling
#curve rates to predict the change in hysteresis.
#Date: Jan 17, 2026
#Last Updated: 19-Oct-26
#Platform: R-Studio 2026.01.0 Build 392, on R version 4.5.2
#With Assistance From: ChatGPT Luna 5.6, Claude.AI Sonnet 5
#External Libraries used: deSolve
#(Tools -> Install Packages -> select deSolve)
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
#Uncomment the appropriate lines below for your platform.
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
Lstop <- 48  # stop temperature of lower baseline fit
Ustart <- 67 # start temperature of upper baseline fit
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
#points(CC$Temperature,CC$Abs, col="blue")

###### END OF HEATING CURVE BASELINE FITTING HERE ######

# Cooling Curve Baselines next:
#Plot the data to decide start and end temperatures for baseline fitting.
plot(CC$Temperature,CC$Abs,main="Absorbance vs. Temperature", xlab="Temperature (°C)", ylab = "OD", col="blue")

Lstart <- 30 # start temperature of lower baseline fit
Lstop <- 44  # stop temperature of lower baseline fit
Ustart <- 67 # start temperature of upper baseline fit
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
#points(HC$Temperature,HC$Abs, col="red")

###### END OF COOLING CURVE BASELINE FITTING HERE ######

#NOTE: FINALIZE THE BASELINES before proceeding.
#Visually examine the baselines. Are they nice? If not, change limits above.
#Also considering trimming out garbage data at temperature extremes.
#This plot needs to be pretty!

#Now calculate Alpha vs. Temperature curves
#Alpha = (Upper baseline - Alpha) / (upper baseline - lower baseline)
HC_Alpha <- ((HC_Uint+HC_Uslope*HC$Temperature)-HC$Abs)/((HC_Uint+HC_Uslope*HC$Temperature)-(HC_Lint+HC_Lslope*HC$Temperature))
CC_Alpha <- ((CC_Uint+CC_Uslope*CC$Temperature)-CC$Abs)/((CC_Uint+CC_Uslope*CC$Temperature)-(CC_Lint+CC_Lslope*CC$Temperature))

#Define model parameters
Ct <- 2.0e-6  #2 µM in molar (Should be initial concentration of dsDNA)
R <- 8.314    #gas constant, J mol^-1 K^-1
Tref <- 60 + 273.15 # choose a reference temperature not off in space (better for fitting). Make this close to or equal to the Tm.

library(deSolve) #load the deSolve library, we are going to use it

# Tm and ln_kref fit in natural units; DH and Ea_off are fit in MJ/mol
# (real kJ/mol / 1000) purely for numerical conditioning in nlm — see
# the summary table section below for the conversion back to kJ/mol.

# Generalize Keq to take trial (Tm, DH) 
Keq_gen <- function(Temp, Tm, DH_MJ) {
  T <- Temp + 273.15
  Tm_K <- Tm + 273.15
  (4/Ct) * exp(((-DH_MJ*1000*1000)/(R*T)) * (1 - T/Tm_K))
}

# Bounded version to prevent overflow
#Keq_gen <- function(Temp, Tm, DH) {
#  T <- Temp + 273.15
#  Tm_K <- Tm + 273.15
#  expo <- ((-DH*1000)/(R*T)) * (1 - T/Tm_K)
#  expo <- pmin(pmax(expo, -700), 700)   # exp() overflows/underflows outside roughly ±709
#  (4/Ct) * exp(expo)
#}

#Write a function to calculate ln(koff) as a function of temperature.
#This is our Arrhenius relationship for the rate constant.
ln_koff <- function(Temp,ln_kref,Ea_off_MJ) {
  T <- Temp + 273.15
  ln_kref - (Ea_off_MJ*1000*1000/R)*((1/T) - (1/Tref))
}

# Given trial (Tm, DH), compute the equilibrium starting alpha at a given start temp
# This uses the first root only, which should be the correct one
alpha0_at <- function(Tstart, Tm, DH_MJ) {
  Kt0 <- Keq_gen(Tstart, Tm, DH_MJ)
  b <- -2 - 2/(Kt0*Ct)
  disc <- max(b^2 - 4, 0)   # guard tiny negative disc from floating point at Kt0 -> Inf
  r1 <- (-b - sqrt(disc)) / 2   # this is always the physical root in [0,1]
  r2 <- (-b + sqrt(disc)) / 2   # its reciprocal partner, always >= 1
  r1
}

# Generalized ODE, all 4 kinetic/thermo params passed in explicitly
ODE_gen <- function(time, state, parameters, TempFun) {
  alpha <- state[1]
  T <- TempFun(time)
  Tm <- parameters["Tm"]
  DH <- parameters["DH"]
  ln_kref <- parameters["ln_kref"]; Ea_off <- parameters["Ea_off"]
  ln_koff_T <- ln_koff(T, ln_kref, Ea_off)
  ln_koff_T <- min(ln_koff_T, 700)
  dalpha_dt <- exp(ln_koff_T) * (Keq_gen(T, Tm, DH)*Ct/2*(1-alpha)^2 - alpha)
  if (!is.finite(dalpha_dt)) dalpha_dt <- 0   # last-resort net: any NaN/Inf just stalls, doesn't crash
  list(c(dalpha_dt))
}

simulate_curve <- function(p, TimeVec, TempFun) {
  Tstart <- TempFun(TimeVec[1])
  a0 <- alpha0_at(Tstart, p["Tm"], p["DH"])
  out <- ode(y = c(alpha = a0), times = TimeVec,
             func = function(t,y,parms) ODE_gen(t,y,parms,TempFun), parms = p)
  out[,2]
}

#Create a function for the heating curve temperature vs. time data:
#Interpolate the data using R's very convenient interpolation function:
HC_Temp <- function(time) { # this function won't crap out outside time limits.
  time <- pmin(pmax(time, min(HC$Time)), max(HC$Time))
  approx(HC$Time, HC$Temperature, xout = time)$y
}

CC_Temp <- function(time) { # this function won't crap out outside time limits.
  time <- pmin(pmax(time, min(CC$Time)), max(CC$Time))
  approx(CC$Time, CC$Temperature, xout = time)$y
}

# Joint objective: sum of squared residuals from BOTH curves, same 4 parameters
joint_chisq <- function(p_vec) {
  p <- c(Tm = p_vec[1], DH = p_vec[2], ln_kref = p_vec[3], Ea_off = p_vec[4])
  HC_pred <- simulate_curve(p, HC$Time, HC_Temp)
  CC_pred <- simulate_curve(p, CC$Time, CC_Temp)   # your existing cooling Temp()
  sum((HC_Alpha - HC_pred)^2) + sum((CC_Alpha - CC_pred)^2)
}

# Quick equilibrium-only pre-fit, just to get good starting guesses for Tm, DH (use results from other fitting routines?)
Eq_alpha_simple <- function(Temp, Tm, DH) { #DH in kJ/mol
  Kt <- Keq_gen(Temp, Tm, DH)
  b <- -2 - 2/(Kt*Ct)
  disc <- pmax(b^2 - 4, 0)
  (-b - sqrt(disc)) / 2
}
Eq_chisq <- function(p) sum((HC_Alpha - Eq_alpha_simple(HC$Temperature, p[1], p[2]))^2)
g1 <- 58.0 # First guess for Tm
g2 <- -0.710 # First guess for DH, in kJ/mol
Eq_guess <- nlm(Eq_chisq, p = c(g1, g2))
Eq_guess$estimate   # use these as Tm0, DH0 for p0

# Option 1: Use the quick fit results above
p0 <- c(Tm = Eq_guess$estimate[1], DH = Eq_guess$estimate[2],
        ln_kref = -6.294445 , Ea_off = 0.33)

# Option 2:
# Starting guesses: Brute Force Guessing (for the DIY'er)
p0 <- c(Tm = 57.6, DH = -0.588,
        ln_kref = -5.1 , Ea_off = 0.322)

# ---- Visualize the starting guesses before fitting ----
plot_guess <- function(p_guess) {
  p_g <- c(Tm = unname(p_guess[1]), DH = unname(p_guess[2]),
           ln_kref = unname(p_guess[3]), Ea_off = unname(p_guess[4]))
  
  HC_guess <- simulate_curve(p_g, HC$Time, HC_Temp)
  CC_guess <- simulate_curve(p_g, CC$Time, CC_Temp)
  
  # Pure equilibrium curve for reference (what you'd see with infinitely slow ramps)
  Tseq   <- seq(min(HC$Temperature), max(HC$Temperature), length.out = 300)
  Eq_ref <- Eq_alpha_simple(Tseq, p_g["Tm"], p_g["DH"])
  
  plot(HC$Temperature, HC_Alpha, col = "red",
       xlab = "Temperature (°C)", ylab = expression(alpha),
       main = "Starting guesses (before fitting)", ylim = c(0, 1))
  points(CC$Temperature, CC_Alpha, col = "blue")
  lines(Tseq, Eq_ref, col = "gray50", lty = 2, lwd = 1.5)
  lines(spline(HC$Temperature, HC_guess), col = "red",  lwd = 2)
  lines(spline(CC$Temperature, CC_guess), col = "blue", lwd = 2)
  legend("bottomleft",
         legend = c("Heating data", "Cooling data",
                    "Heating (guess)", "Cooling (guess)", "Equilibrium (guess)"),
         col = c("red", "blue", "red", "blue", "gray50"),
         pch = c(1, 1, NA, NA, NA), lty = c(NA, NA, 1, 1, 2),
         lwd = c(NA, NA, 2, 2, 1.5), bty = "n")
  
  # Numerical score for the guess
  ss_HC <- sum((HC_Alpha - HC_guess)^2)
  ss_CC <- sum((CC_Alpha - CC_guess)^2)
  cat(sprintf("SSR heating: %.5f | SSR cooling: %.5f | total: %.5f\n",
              ss_HC, ss_CC, ss_HC + ss_CC))
  invisible(list(HC = HC_guess, CC = CC_guess))
}
plot_guess(p0)

# After the starting guesses above give a reasonable fit, perform the actual fit here:
joint_fit <- nlm(
  joint_chisq,
  p = p0,
  hessian = TRUE,
  ndigit = 10,
  gradtol = 1e-8,
  iterlim = 200
)

joint_fit$estimate   # your 4 globally-fit parameters
joint_fit # see the results of the joint fit

#Calculate Correlation Matrix: Method 1
joint_cov.mat <- 2 * joint_fit$minimum / (length(HC_Alpha) + length(CC_Alpha) - 4) * solve(joint_fit$hessian)
joint_cor.mat <- cov2cor(joint_cov.mat)
rownames(joint_cor.mat) <- colnames(joint_cor.mat) <- c("Tm","DH","ln_kref","Ea_off")
joint_cor.mat

#calculate Correlation Matrix: Method 2
N <- length(HC_Alpha) + length(CC_Alpha)
P <- length(joint_fit$estimate)
rss <- joint_fit$minimum
sigma2_hat <- rss / (N - P)
H <- (joint_fit$hessian + t(joint_fit$hessian)) / 2
joint_cov.mat <- 2 * sigma2_hat * solve(H)
joint_se <- sqrt(diag(joint_cov.mat))
joint_cor.mat <- cov2cor(joint_cov.mat)
joint_cor.mat

joint_se <- sqrt(diag(joint_cov.mat))
names(joint_se) <- c("Tm","DH","ln_kref","Ea_off")
joint_se

p <- joint_fit$estimate
names(p) <- c("Tm","DH","ln_kref","Ea_off")

HC_joint_pred <- simulate_curve(p, HC$Time, HC_Temp)
CC_joint_pred <- simulate_curve(p, CC$Time, CC_Temp)

# Equilibrium curve from the fitted Tm and DH (infinitely slow ramp)
Tseq   <- seq(min(HC$Temperature), max(HC$Temperature), length.out = 300)
Eq_fit <- Eq_alpha_simple(Tseq, p["Tm"], p["DH"])

plot(HC$Temperature, HC_Alpha, col="red", main="Joint kinetic fit: both curves",
     xlab="Temperature (°C)", ylab=expression(alpha), ylim=c(0, 1))
points(CC$Temperature, CC_Alpha, col="blue")
lines(Tseq, Eq_fit, col="gray50", lty=2, lwd=1.5)
lines(spline(HC$Temperature, HC_joint_pred), col="red", lwd=2)
lines(spline(CC$Temperature, CC_joint_pred), col="blue", lwd=2)
legend("topright",
       legend = c("Heating data", "Cooling data",
                  "Heating (fit)", "Cooling (fit)", "Equilibrium (fit)"),
       col = c("red", "blue", "red", "blue", "gray50"),
       pch = c(1, 1, NA, NA, NA), lty = c(NA, NA, 1, 1, 2),
       lwd = c(NA, NA, 2, 2, 1.5), bty = "n")


#############################################
#   Joint Kinetic Fit — Summary Parameters  #
#############################################

p <- joint_fit$estimate
names(p) <- c("Tm","DH","ln_kref","Ea_off")

Tm      <- p["Tm"]
DH      <- p["DH"]*1000
ln_kref <- p["ln_kref"]
Ea_off  <- p["Ea_off"]*1000

# --- Derived kinetic/thermodynamic parameters ---

# True intercept of ln(koff) at 1/T = 0 (i.e. koff extrapolated to T -> infinity)
ln_k0_off <- unname(ln_kref + Ea_off*1000/(R*Tref))

# Ea_on from DH = Ea_on - Ea_off  (van't Hoff / Arrhenius consistency). Answer is in kJ.
Ea_on <- unname(DH + Ea_off)

# True intercept of ln(kon) at 1/T = 0, derived analytically from
# ln(kon) = ln(Keq) + ln(koff), combining the Keq_gen and koff Arrhenius forms
Tm_K <- Tm + 273.15
ln_k0_on <- unname(log(4/Ct) + DH*1000/(R*Tm_K) + ln_kref + Ea_off*1000/(R*Tref))

# --- Summary table ---

Joint_Fit_Summary <- data.frame(
  Parameter = c(
    "Tm",
    "DH",
    "ln(kref)",
    "ln(koff0)",
    "Ea_off",
    "ln(kon0)",
    "Ea_on",
    "Tref"
  ),
  Value = c(
    Tm,
    DH,
    ln_kref,
    ln_k0_off,
    Ea_off,
    ln_k0_on,
    Ea_on,
    Tref-273.15
  ),
  Units = c(
    "°C",
    "kJ/mol",
    "dimensionless",
    "dimensionless",
    "kJ/mol",
    "dimensionless",
    "kJ/mol",
    "°C"
  )
)

Joint_Fit_Summary$Value <- round(Joint_Fit_Summary$Value, 3)
Joint_Fit_Summary

#None of these values should be negative:
eigen(joint_fit$hessian)$values
# eigenvalues should be comfortably away from zero.
eig <- eigen((H + t(H)) / 2, symmetric = TRUE)$values
eig



##### SIMULATION CORNER (Fitting Routine is Done) ####
T_rate <- 0.5 / 60 # heating (and cooling) rate in °C/min (converted to seconds)

# Create a Linear Ramp Function:
make_ramp <- function(T_start, ramp_rate) {
  force(T_start)
  force(ramp_rate)
  
  function(time) {
    T_start + ramp_rate * time
  }
}

# Now simulate heating and cooling curves:

# Get estimates:
p <- joint_fit$estimate
names(p) <- c("Tm", "DH", "ln_kref", "Ea_off")

# Define temperature ranges:
T_low <- 40
T_high <- 70

heating_rate <- T_rate     # °C/min
cooling_rate <- -T_rate    # °C/min

t_heat <- seq(
  0,
  (T_high - T_low) / heating_rate,
  length.out = 1000
)

t_cool <- seq(
  0,
  (T_low - T_high) / cooling_rate,
  length.out = 1000
)

# Simulate the heating and cooling curves:
HC_sim <- simulate_curve(
  p = p,
  TimeVec = t_heat,
  TempFun = make_ramp(T_low, heating_rate)
)

CC_sim <- simulate_curve(
  p = p,
  TimeVec = t_cool,
  TempFun = make_ramp(T_high, cooling_rate)
)

# Generate the corresponding temperatures:
T_heat <- make_ramp(T_low, heating_rate)(t_heat)
T_cool <- make_ramp(T_high, cooling_rate)(t_cool)

# Plot the curves:
plot(
  T_heat, HC_sim,
  type = "l",
  col = "red",
  lwd = 2,
  xlim = c(T_low, T_high),
  ylim = c(0, 1),
  xlab = "Temperature (°C)",
  ylab = expression(alpha),
  main = "Simulated kinetic hysteresis"
)

lines(
  T_cool, CC_sim,
  col = "blue",
  lwd = 2
)

#Overlay experimental points:
points(HC$Temperature, HC_Alpha, col="red")
points(CC$Temperature, CC_Alpha, col="blue")

legend(
  "topright",
  legend = c(
    paste0("Heating: ", heating_rate*60, " °C/min"),
    paste0("Cooling: ", abs(cooling_rate*60), " °C/min")
  ),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n"
)
