#This script fits a binding curve to data (KCl vs. Molar Elipticity)
#Written by: David Dubins
#This program expects two columns of data:
# Temperature   Abs

#For PC:
Mydata <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#Mydata <- read.table(pipe("pbpaste"), header=TRUE)
attach(Mydata)

nr <- 4  #number of roots (also the order of the polynomial to fit)
points <- dim(Mydata)[1] #number of data points to fit
Ct <- 100e-6  #100 uM in molar

library(polynom) #install the polynom package before loading for the first time

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

# Define the model you are using to fit the data (in this case, melting curve):
kFit <- function(Lint,Lslope,Uint,Uslope,Ktm,DH,Tm){
  #Lint: Lower baseline (intercept)
  #Lslope: Lower baseline (slope)
  #Uint: Upper baseline (intercept)
  #Uslope: Upper baseline (slope)
  #K: Duplex formation equilibrium constant
  #Ct #initial concentration of single-stranded oligonucleotides (Ct) in Mol
  # From solving the model:
  #a = 1
  #b = -4
  #c = 6
  #d = -4 -(1/(4*k*Ct^3))
  #e = 1
  
  for(i in 1:points){
    K <- (Ktm*1e+11)*exp(-(DH*1e+5/8.314)*((1/(Temperature[i]+273.15))-(1/(Tm+273.15))))
    a <- 4*K*Ct^3
    b <- -4*a
    c <- 6*a
    d <- -4*a-1
    e <- a
    
    rawroots <- polynomial (c(e,d,c,b,a)) #polynomial goes UP in order
    foundroots <- solve(rawroots)
    #if(i==4){
    #  temp <- foundroots
    #}
    alpha <- findReal(foundroots)
    absorbance = alpha*((Uint+Uslope*Temperature[i])-(Lint+Lslope*Temperature[i]))+(Lint + Lslope*Temperature[i])
    if(i==1){
      results <- absorbance
    } else {
      results <- c(results,absorbance)
    }
  }
  results
}

# Define the chisq function as sum(Yobs - Ymodel)^2:
chisq <- function(p) sum((Abs-kFit(p[1],p[2],p[3],p[4],p[5],p[6],p[7]))^2)
# Now plot your data:
plot(Temperature,Abs)
# Try your hand at fitting your parameters to get in the right ballpark. 
g1 <- 0.10940037  #first guess for Lint
g2 <- 0.000161479  #first guess for Lslope
g3 <- 0.137396944  #first guess for Uint
g4 <- 7.45937E-05  #first guess for Uslope
g5 <- 10           #first guess for Ktm
g6 <- -8           #first guess for DH
g7 <- 71           #first guess for Tm
Myfitx <- Temperature
Myfity <- kFit(g1,g2,g3,g4,g5,g6,g7)
lines(Myfitx, Myfity,type="p", pch=16)

# Ready for fitting the model with your guesses:
Myfitout <- nlm(chisq, p = c(g1, g2, g3, g4, g5, g6, g7), hessian = TRUE)
# See the results of your fit:
Myfitout

# To obtain the approximate standard errors of the parameter estimates:
sqrt(diag(2*Myfitout$minimum/(length(Abs) - 2) * solve(Myfitout$hessian)))

# Now superimpose our minimized fit on a brand new plot:
plot(Temperature,Abs)
Myfitx <- Temperature
Myfity <- kFit(Myfitout$estimate[1],Myfitout$estimate[2],Myfitout$estimate[3],Myfitout$estimate[4],Myfitout$estimate[5],Myfitout$estimate[6],Myfitout$estimate[7])
lines(spline(Myfitx, Myfity))

# You can also calculate the coefficient of determination (r-squared) for your model:
sstot = sum((Abs-mean(Abs))^2)
sserr = sum((Abs - Myfity)^2)
rsq = 1-sserr/sstot
rsq

#Residual Plot:
plot(Myfitx,Abs-kFit(Myfitout$estimate[1],Myfitout$estimate[2],Myfitout$estimate[3],Myfitout$estimate[4],Myfitout$estimate[5],Myfitout$estimate[6],Myfitout$estimate[7]))

#Final parameter estimates:
#this is Ktm:
Myfitout$estimate[5]*1e+11

#this is DH:
Myfitout$estimate[6]*1e+5

#this is Tm:
Myfitout$estimate[7]