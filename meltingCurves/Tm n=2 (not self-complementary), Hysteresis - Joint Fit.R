#Tm n=2 (not self-complementary), Hysteresis - Joint Fit.R
#This script fits TM data assuming D -> SS1 + SS2 (not self-complementary)
#This script was written with a lot of assistance from Claude.AI Sonnet 5.5 Medium and Perplexity.AI.
#The bottom of the script allows you to simulate different heating and cooling
#curve rates to predict the change in hysteresis.
#Date: Jan 17, 2026
#Last Updated: 08-Nov-26
#Platform: R-Studio 2026.01.0 Build 392, on R version 4.5.2
#With Assistance From: ChatGPT Luna 5.6, Claude.AI Sonnet 5
#External Libraries used: deSolve
#(Tools -> Install Packages -> select deSolve)
#This program expects three columns of data: (the unit of time is in seconds)
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

##################################################
# CHECKPOINT 1: READ IN THE DATA. (USER ACTION)  #
# Note: Only import data you want to fit.        #
# Truncate long, noisy, and drifting baselines.  #
##################################################

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

##################################################
# CHECKPOINT 2: FIT THE BASELINES. (USER ACTION) #
##################################################

# The "lower baseline" (Lint + Lslope*T) is the lower temperature baseline - a linear fit of the data at the lowest temperatures.
# The "upper baseline" (Uint + Uslope*T) is the upper temperature baseline - a linear fit of the data at highest temperatures.

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

Lstart <- 43 # start temperature of lower baseline fit
Lstop <- 47  # stop temperature of lower baseline fit
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

Lstart <- 29 # start temperature of lower baseline fit
Lstop <- 36  # stop temperature of lower baseline fit
Ustart <- 66 # start temperature of upper baseline fit
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

#Examine the heating curve and the cooling curve with baselines together:
plot(HC$Temperature,HC$Abs,main="Absorbance vs. Temperature", xlab="Temperature (°C)", ylab = "OD", col="red")
lines(spline(HC$Temperature, HC_LowerBLfit), col="red")
lines(spline(HC$Temperature, HC_UpperBLfit), col="red")
points(CC$Temperature,CC$Abs, col="blue")
lines(spline(CC$Temperature, CC_LowerBLfit), col="blue")
lines(spline(CC$Temperature, CC_UpperBLfit), col="blue")

#NOTE: FINALIZE THE BASELINES before proceeding. This plot needs to be pretty!
#Do your baselines look like they make sense? If not:
# -Consider going back to CHECKPOINT 1 and truncating the data.
# -Go back to CHECKPOINT 2 and re-fit the baselines.

#################################################################
# CHECKPOINT 3: Calculate and Review Alpha Curves (USER ACTION) #
#################################################################

#Now calculate Alpha vs. Temperature curves
#This equation works for increasing or decreasing absorbance with temperature, because of the definition of "Upper and "Lower".

#Alpha = (Upper baseline - Abs) / (Upper baseline - Lower baseline)
HC_Alpha <- ((HC_Uint+HC_Uslope*HC$Temperature)-HC$Abs)/((HC_Uint+HC_Uslope*HC$Temperature)-(HC_Lint+HC_Lslope*HC$Temperature))
CC_Alpha <- ((CC_Uint+CC_Uslope*CC$Temperature)-CC$Abs)/((CC_Uint+CC_Uslope*CC$Temperature)-(CC_Lint+CC_Lslope*CC$Temperature))

#Examine the normalized curves together:
plot(HC$Temperature,HC_Alpha,main="Alpha vs. Temperature", xlab="Temperature (°C)", ylab = "α (Fraction of SS in Duplex)", col="red")
points(CC$Temperature,CC_Alpha, col="blue")
legend("bottomleft",
       legend = c("Heating data", "Cooling data"),
       col = c("red", "blue"),
       pch = c(1, 1), lty = c(NA, NA),
       lwd = c(NA, NA), bty = "n")

#-Does the normalization look good? 
#-Are the curves roughly parallel? 
#-Are the baselines flat?
#-Visually verify that alpha (fraction of SS in Duplex DNA) goes from 1 to 0 as temperature increases.
#-If for some reason the curves are backwards, change your baseline definitions, or subtract them from 1:
#HC_Alpha <- 1 - (((HC_Uint+HC_Uslope*HC$Temperature)-HC$Abs)/((HC_Uint+HC_Uslope*HC$Temperature)-(HC_Lint+HC_Lslope*HC$Temperature)))  # Invert HC
#CC_Alpha <- 1 - (((CC_Uint+CC_Uslope*CC$Temperature)-CC$Abs)/((CC_Uint+CC_Uslope*CC$Temperature)-(CC_Lint+CC_Lslope*CC$Temperature)))  # Invert CC

###############################################
# CHECKPOINT 4: Model Definition (NO ACTION)  #
###############################################

#Define model parameters
# Ct is the total concentration of strands in the system (it's a constant). 
# At any given time, Ct = [SS1] + [SS2] + 2[D].
# You know Ct, because before you begin, you've added SS1 and SS2 together.
# When DNA is all in duplex form, there should be Ct/2 mol/L of dsDNA.
Ct <- 2.0e-6  #Total concentration of strands, in molar (Ct=(initial [SS1] + [SS2], i.e. 2 × the concentration of each strand).
R <- 8.314    #gas constant, J mol^-1 K^-1
Tref <- NA_real_ # Tref is set automatically in CHECKPOINT 5 (= Tm0, in K)

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

########################################################
# CHECKPOINT 5: Calculate Starting Guesses (NO ACTION) #
########################################################

# The routine will try to estimate reasonable starting guesses.
# ========================-------==============================
# Temperature at which a smoothed alpha curve crosses a given level:
cross_T <- function(Temp, Alpha, level) {
  s  <- smooth.spline(Temp, Alpha)
  Tg <- seq(min(Temp), max(Temp), length.out = 2000)
  Tg[which.min(abs(predict(s, Tg)$y - level))]
}

# Average ramp rate in degC/s:
ramp_rate <- function(D) (max(D$Temperature) - min(D$Temperature)) / (max(D$Time) - min(D$Time))

# 1. Tm: average of the two alpha = 0.5 crossings; gap = hysteresis width
Tm_HC <- cross_T(HC$Temperature, HC_Alpha, 0.5)
Tm_CC <- cross_T(CC$Temperature, CC_Alpha, 0.5)
Tm0   <- (Tm_HC + Tm_CC) / 2
gap   <- max(abs(Tm_HC - Tm_CC), 0.1)   # floor avoids divide-by-zero

# 2. DH (MJ/mol, negative as in the model): width between alpha = 0.75 and 0.25.
#    For a bimolecular melt, ln(27) = 3.296 = |DH|*width/(R*Tm^2).
#    Kinetics broaden curves, so use the narrower (steeper) of the two.
width <- function(D, A) abs(cross_T(D$Temperature, A, 0.25) - cross_T(D$Temperature, A, 0.75))
w     <- min(width(HC, HC_Alpha), width(CC, CC_Alpha))
DH0   <- -log(27) * R * (Tm0 + 273.15)^2 / w / 1e6

# 3. ln_kref: set Tref = Tm so ln_kref = ln(koff) at Tm.
#    Linear response near alpha = 0.5: relaxation rate = 3*koff, lag = rate/relaxation rate,
#    so koff ~ 2*ramp/(3*gap).
Tref     <- Tm0 + 273.15          # ln_koff() reads this global
rate     <- mean(c(ramp_rate(HC), ramp_rate(CC)))
lnkref0  <- log(2 * rate / (3 * gap))

# 4. Ea_off: coarse SSR scan (no fitting). Ea_off is a fraction of |DH|;
#    ln_kref is scanned +/- 1 as well, since the gap formula is only approximate.
grid <- expand.grid(frac = c(0.3, 0.5, 0.7, 0.9), dlnk = c(-1, 0, 1))
ssr  <- apply(grid, 1, function(g)
  tryCatch(joint_chisq(c(Tm0, DH0, lnkref0 + g["dlnk"], g["frac"] * abs(DH0))),
           error = function(e) NA_real_))
best <- grid[which.min(ssr), ]

p0 <- c(Tm = Tm0, DH = DH0,
        ln_kref = lnkref0 + best$dlnk,
        Ea_off  = best$frac * abs(DH0))
p0  # have a look at p0
cbind(grid, ssr)   # if ssr barely changes with frac, Ea_off is poorly determined

# ---> SKIP CHECKPOINT 6 for now (you will go back to it later). Proceed directly to CHECKPOINT 7!

######################################################################
# CHECKPOINT 6: Manual Refinement of Starting Guesses (USER ACTION)  #
######################################################################

# Here are the last guesses I tried from CHECKPOINT 5:
p0
# Tm is the temperature at which alpha=0.5 at equilibrium (should be midpoint)
# DH controls the width (slope) of the alpha curve at each curve at the transition
# ln_kref is ln(koff) at Tref; it controls width of hysteresis (smaller = less hysteresis)
# Ea sets how hysteresis changes with temperature (also controls slope and shape of curve)

# Manually change the guesses inside the function below for Tm, DH, ln_kref, and Ea_off:
p0 <- c(Tm = 57.4, DH = -0.54,
        ln_kref = -5.9, Ea_off = 0.32)
# Now execute CHECKPOINT 7 (visualize the refined guesses).

##########################################################
# CHECKPOINT 7: Visualize Starting Guesses (USER REVIEW) #
##########################################################

# Visualize the starting guesses before fitting:
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

# Now go back to CHECKPOINT 6, and tweak the guesses until you get a closer fit.
# Good starting guesses are very important for the fitting routine to be successful.
# When you are happy with the fit (close enough), proceed to CHECKPOINT 8.

########################################################
# CHECKPOINT 8: RUN THE FITTING ROUTINE (USER REVIEW)  #
########################################################

# Perform the actual fit here:
joint_fit <- nlm(
  joint_chisq,
  p = p0,
  hessian = TRUE,
  ndigit = 10,
  gradtol = 1e-8,
  iterlim = 200
)
#This step may take several seconds to minutes depending on
#the data set and quality of the starting guesses.

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

# Inspect the output of joint_fit. We are hoping for a $code [1] or [2] ([3] is ok but not great) 
# Inspect joint_cor.mat, we are hoping for no NaN's and nothing greater than 1.0.
# If the fit is poor, go back to CHECKPOINT 6, change your starting guesses, then try
# CHECKPOINT 7 and then CHECKPOINT 8 again.
# Sometimes it's better to be close, but not bang on with your guesses.

#####################################################################
# CHECKPOINT 9: ROLL UP AND REPORT THE FINAL RESULTS (USER ACTION)  #
#####################################################################

#Joint Kinetic Fit — Summary Parameters (The "final" answer)
#Run this line to store Tref:
Tref_fit    <- Tref            # Tm0 in K, what the fit used. Do NOT overwrite global Tref.

# The fit used a reference temperature Tref_fit (= Tm0 from CHECKPOINT 5, in K).
# ln_kref is ln(koff) at THAT temperature, so it can't be compared between experiments.
# For comparisons, report koff at a common temperature, Tref_report, the same for all
# experiments. Choose a round number close to the Tms of everything being compared.
# ln(koff0) and ln(kon0) (Arrhenius intercepts) do not depend on either Tref.

# Select ONE of the following (comment out the other):
Tref_report <- 55 + 273.15      # common reporting temperature (K)
#Tref_report <- Tref_fit         # or: report at the temperature used for fitting (Tm0)


p <- joint_fit$estimate
names(p) <- c("Tm","DH","ln_kref","Ea_off")
Tm      <- p["Tm"]
DH      <- p["DH"]*1000
ln_kref <- p["ln_kref"]  # ln(koff) at Tref_fit
Ea_off  <- p["Ea_off"]*1000

# --- Derived kinetic/thermodynamic parameters ---

# Tref-independent intercept, using the Tref the fit was actually done at
ln_k0_off <- unname(ln_kref + Ea_off*1000/(R*Tref_fit))
# ln(koff) at the common reporting temperature
ln_kref_report <- unname(ln_k0_off - Ea_off*1000/(R*Tref_report))

# Ea_on from DH = Ea_on - Ea_off  (van't Hoff / Arrhenius consistency). Answer is in kJ.
Ea_on <- unname(DH + Ea_off)

# True intercept of ln(kon) at 1/T = 0, derived analytically from
# ln(kon) = ln(Keq) + ln(koff), combining the Keq_gen and koff Arrhenius forms
Tm_K <- Tm + 273.15
ln_k0_on <- unname(log(4/Ct) + DH*1000/(R*Tm_K) + ln_k0_off)

# --- Summary table ---

koff_report <- exp(ln_kref_report)   # s^-1 at Tref_report

# --- Standard errors (delta method) ---
# joint_cov.mat is in FIT units: (Tm [°C], DH [MJ/mol], ln_kref, Ea_off [MJ/mol])
# For a derived quantity f, SE = sqrt(g' V g), where g = gradient of f wrt those 4 parameters.
V <- joint_cov.mat
se_of <- function(g) sqrt(drop(t(g) %*% V %*% g))

Tm_K_u <- unname(Tm_K)
DH_J   <- unname(DH) * 1000      # DH is in kJ/mol here, so this is J/mol

se_Tm        <- se_of(c(1, 0, 0, 0))                                   # °C
se_DH        <- se_of(c(0, 1, 0, 0)) * 1000                            # kJ/mol
se_ln_kref_r <- se_of(c(0, 0, 1, 1e6/R * (1/Tref_fit - 1/Tref_report)))
se_koff      <- koff_report * se_ln_kref_r                             # first-order: dk = k * d(ln k)
se_ln_k0_off <- se_of(c(0, 0, 1, 1e6/(R*Tref_fit)))
se_Ea_off    <- se_of(c(0, 0, 0, 1)) * 1000                            # kJ/mol
se_ln_k0_on  <- se_of(c(-DH_J/(R*Tm_K_u^2), 1e6/(R*Tm_K_u), 1, 1e6/(R*Tref_fit)))
se_Ea_on     <- se_of(c(0, 1, 0, 1)) * 1000                            # kJ/mol

# Calculate kon:
# ln(kon) and kon at the common reporting temperature (M^-1 s^-1)
ln_kon_report <- unname(ln_k0_on - Ea_on*1000/(R*Tref_report))
kon_report    <- exp(ln_kon_report)
# Calculate the standard error for kon:
# ln_kon_report = ln(4/Ct) + DH/(R·Tm_K) + ln_kref + Ea_off/(R·Tref_fit) − (DH+Ea_off)/(R·Tref_report)
se_ln_kon_r <- se_of(c(-DH_J/(R*Tm_K_u^2),                    # d/dTm
                       1e6/R * (1/Tm_K_u - 1/Tref_report),    # d/dDH
                       1,                                     # d/dln_kref
                       1e6/R * (1/Tref_fit - 1/Tref_report))) # d/dEa_off
se_kon <- kon_report * se_ln_kon_r                            # first-order



# --- Summary table ---
Joint_Fit_Summary <- data.frame(
  Parameter = c("Tm", "DH",
                "ln(koff) at Tref (report)", "koff",
                "ln(kon) at Tref (report)",  "kon",
                "ln(koff0)", "Ea_off", "ln(kon0)", "Ea_on",
                "Tref (report)", "Tref (fit)"),
  Value = unname(c(Tm, DH,
                   ln_kref_report, koff_report,
                   ln_kon_report,  kon_report,
                   ln_k0_off, Ea_off, ln_k0_on, Ea_on,
                   Tref_report - 273.15, Tref_fit - 273.15)),
  SE    = c(se_Tm, se_DH,
            se_ln_kref_r, se_koff,
            se_ln_kon_r,  se_kon,
            se_ln_k0_off, se_Ea_off, se_ln_k0_on, se_Ea_on,
            NA, NA),                      # Tref values are chosen, not fitted
  Units = c("°C", "kJ/mol",
            "ln(s^-1)", "s^-1",
            "ln(M^-1 s^-1)", "M^-1 s^-1",
            "ln(s^-1)", "kJ/mol", "ln(M^-1 s^-1)", "kJ/mol",
            "°C", "°C")
)

# Format each row individually (decimals per row; NA = scientific notation)
dec <- c(2,   # Tm
         1,   # DH
         3,   # ln(koff) at Tref (report)
         NA,  # koff (scientific)
         3,   # ln(kon) at Tref (report)
         NA,  # kon (scientific)
         3,   # ln(koff0)
         1,   # Ea_off
         3,   # ln(kon0)
         1,   # Ea_on
         2,   # Tref (report)
         2)   # Tref (fit)

fmt <- function(v, d, extra = 0) {
  if (is.na(v)) return("-")
  if (is.na(d)) formatC(v, format = "e", digits = 3)
  else          formatC(v, format = "f", digits = d + extra)
}

Joint_Fit_Summary$Value <- mapply(fmt, Joint_Fit_Summary$Value, dec)
Joint_Fit_Summary$SE    <- mapply(fmt, Joint_Fit_Summary$SE,    dec, MoreArgs = list(extra = 1))

print(Joint_Fit_Summary, row.names = FALSE, right = TRUE)

#None of these values should be negative:
eigen(joint_fit$hessian)$values
# eigenvalues should be comfortably away from zero.
eig <- eigen((H + t(H)) / 2, symmetric = TRUE)$values
eig

#Calculate the global r^2 value:
HC_joint_pred <- simulate_curve(p, HC$Time, HC_Temp)
CC_joint_pred <- simulate_curve(p, CC$Time, CC_Temp)
# Combine observed and predicted values
obs  <- c(HC_Alpha, CC_Alpha)
pred <- c(HC_joint_pred, CC_joint_pred)

#Now calculate the r^2 coefficeint of determination for the global and individual fits:
SS_res <- sum((obs - pred)^2)
SS_tot <- sum((obs - mean(obs))^2)

R2 <- 1 - SS_res/SS_tot
# r^2 value of global fit:
R2

#Calculate the Heating Curve and Cooling Curve r^2 values separately:
R2_HC <- 1 - sum((HC_Alpha - HC_joint_pred)^2) /
  sum((HC_Alpha - mean(HC_Alpha))^2)

R2_CC <- 1 - sum((CC_Alpha - CC_joint_pred)^2) /
  sum((CC_Alpha - mean(CC_Alpha))^2)

# r^2 value of heating curve fit:
R2_HC
# r^2 value of cooling curve fit:
R2_CC

#--------END OF FITTING ROUTINE!----------#


#################################################################
# CHECKPOINT 10: SIMULATING FITS (Not part of fitting routine)  #
#################################################################

# This part of the program takes your fitted parameters above and test drives them.
# You can change the heating/cooling rate and see what happens to hysteresis.
# Your original data will still be plotted on the graph.

T_rate <- 0.5 / 60 # Change the heating (and cooling) rate here in °C/min (converted to seconds)

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

################################################################
#                   ANALYSIS COMPLETE!                         #
################################################################


##########################################################################
# APPENDIX: ARRHENIUS EQUATIONS (Written by Claude.AI Sonnet 5.5 Medium) #
##########################################################################

# Arrhenius form (exponential)
# koff(T) = koff0 * exp(−Ea_off / RT)
# kon(T) = kon0 * exp(−Ea_on / RT)

# Log (linear) form, as used in the script:
# ln koff(T) = ln(koff0) − Ea_off / (R·T)
# ln kon(T) = ln(kon0) − Ea_on / (R·T)

# In the log-linear form, each is a straight line in ln k vs. 1/T. 
# The slope is −Ea/R and the intercept at 1/T = 0 is ln(k0). 
# In the script, ln_k0_off and ln_k0_on are ln(koff0) and ln(kon0), 
# and Ea_off and Ea_on are the activation energies.

# Equivalent form with a reference temperature (what the fit uses)
# ln koff(T) = ln_kref − (Ea_off/R)·(1/T − 1/Tref)

# This is the same line, re-centered at Tref. That is why:
# ln(koff0) = ln_kref + Ea_off / (R·Tref)

# Equilibrium constant and its link to the rates:

# Keq = kon / koff
# ln Keq(T) = ln(4/Ct) − (ΔH/R)·(1/T − 1/Tm)

# This is Keq_gen, with ΔH negative (association enthalpy). 
# Taking the log of kon = Keq·koff gives:

# ln(kon0) = ln(4/Ct) + ΔH/(R·Tm) + ln(koff0)
# Ea_on = ΔH + Ea_off

# Ultimately, we are fitting:
# dα/dt = kon(T)·(Ct/2)(1−α)² − koff(T)·α
# After substitutions, the equation becomes:
# dα/dt = koff(T)·[ Keq(T)·(Ct/2)(1−α)² − α ]
# dα/dt = exp[ ln_kref − (Ea_off/R)(1/T − 1/Tref) ] × ( 2·exp[ −(ΔH/R)(1/T − 1/Tm) ]·(1−α)² − α )
# This is the function ODE_gen() in the program.
