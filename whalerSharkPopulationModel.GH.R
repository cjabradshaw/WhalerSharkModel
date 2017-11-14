############################################
## Whaler shark population projection model
############################################

## Corey J.A. Bradshaw, Flinders University
## corey.bradshaw@flinders.edu.au
## November 2017

## Remove everything
rm(list = ls())

# Set functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc.glm <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]}
  return(AICc.vec)}

k.glm <- function(x) {
  if (length(x$df.residual) == 0) k <- sum(x$dims$ncol) else k <- (length(x$coeff)+1)}

delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
chdev.glm <- function(x) ((( as.numeric(x[12]) - as.numeric(x[10]))/ as.numeric(x[12]))*100) ## % change in deviance, where x is glm object

## source (place in appropriate directory - this is only a suggestion based on RStudio in MacOS)
source("/Applications/RStudio.app/Contents/Resources/R/matrixOperators.r") 

age.vec <- seq(0,39,1)
lage <- length(age.vec)

## construct matrix
stages <- lage
popmat <- matrix(0,nrow=stages,ncol=stages)
colnames(popmat) <- age.vec[1:stages]
rownames(popmat) <- age.vec[1:stages]

## fertility data 
# NB: reproductive cycle of C. brachyurus is probably biennial (potentailly triennial)
# South Australia
f.mat.TL <- c(2880,3020,2950,2800)
f.lit <- c(19,20,26,14)
Lm.fem <- 308; alph.fem <- 74; k.fem <- 0.146
f.mat.age <- (log((f.mat.TL/10)*(Lm.fem-alph.fem))-log(alph.fem)-log(Lm.fem-(f.mat.TL/10)))/k.fem
plot(f.mat.age,f.lit,pch=19)

# youngest mature female 271 cm (16 years old)
# oldest immature female 20 years (250 cm LT)
# proportion mature (Lucifora et al. 2005)
setwd("~/...") # set as appropriate
prop.mat <- read.table("prop.mat.csv",header=T,sep=",")
f.age.maturity <- (log((prop.mat$TL)*(Lm.fem-alph.fem))-log(alph.fem)-log(Lm.fem-(prop.mat$TL)))/k.fem
f.prop.mature <- prop.mat$prop.mat
plot(f.age.maturity,f.prop.mature,pch=19,xlab="age (yrs)", ylab="proportion mature")
mat.dat.out <- data.frame(f.age.maturity, f.prop.mature)

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
mat.dat <- data.frame(f.age.maturity,f.prop.mature)
param.init <- c(1.003150891098860E+00, 1.435062082948057E+01, -3.991451548741554E+01)
fit.logp <- nls(f.prop.mature ~ a / (1+(f.age.maturity/b)^c), 
               data = mat.dat,
               algorithm = "port",
               start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
               trace = TRUE,      
               nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.logp.summ <- summary(fit.logp)
plot(f.age.maturity,f.prop.mature,pch=19,xlab="age (yrs)", ylab="proportion mature")
age.vec.cont <- seq(1,max(age.vec),0.02)
pred.p.mat <- coef(fit.logp)[1] / (1+(age.vec.cont/coef(fit.logp)[2])^coef(fit.logp)[3])
pred.p.matt <- ifelse(pred.p.mat > 1, 1, pred.p.mat)
lines(age.vec.cont,pred.p.matt,lty=2,lwd=3,col="red")
mat.fit.out <- data.frame(age.vec.cont, pred.p.matt)


## modify so that b = new value
b.new <- 16
plot(f.age.maturity,f.prop.mature,pch=19, xlab="age (yrs)", ylab="proportion mature")
age.vec.cont <- seq(1,max(age.vec),0.02)
pred.p.mat <- coef(fit.logp)[1] / (1+(age.vec.cont/b.new)^coef(fit.logp)[3])
pred.p.mat2 <- ifelse(pred.p.mat > 1, 1, pred.p.mat)
lines(age.vec.cont,pred.p.mat,lty=2,lwd=3,col="red")
out.mat <- data.frame(age.vec.cont,pred.p.mat2)
sub.mat <- which(out.mat$age.vec.cont > 15 & out.mat$age.vec.cont < 21)
out.mat[sub.mat,]
mat.fit2.out <- data.frame(age.vec.cont, pred.p.mat2)

# Cliff & Dudley 1992
# TL <- 1.364*PCL - 35.924 # transform precaudal length to TL (Drew et al. 2017)
setwd("~/...") # set as appropriate
size.litter <- read.table("CD92.size.litter.csv",header=T,sep=",")

# convert TL to age to compare to Drew et al. 2017
CD.age <- (log((size.litter$TL)*(Lm.fem-alph.fem))-log(alph.fem)-log(Lm.fem-(size.litter$TL)))/k.fem
plot(CD.age, size.litter$litter.size, pch=19, xlab="age (yrs)", ylab="litter size")

# combine with Drew et al. data
f.age <- c(f.mat.age,CD.age)
f.litter <- c(f.lit, size.litter$litter.size)
plot(f.age, f.litter, pch=19, xlab="age (yrs)", ylab="litter size")
points(f.mat.age,f.lit,col="red",pch=19)
AFR.litter.out <- data.frame(f.age,f.litter)
AUS.litter.out <- data.frame(f.mat.age,f.lit)

# fit exponential model
fit.exp <- lm(log(f.litter) ~ f.age)
fit.exp.summ <- summary(fit.exp)
exp.age.vec <- seq(13,35,0.2)
litt.pred <- exp(coef(fit.exp)[1] + f.age*coef(fit.exp)[2])
litt.up <- as.vector(exp(predict.lm(fit.exp,se.fit=T,interval="prediction")$fit[,3]))
litt.lo <- as.vector(exp(predict.lm(fit.exp,se.fit=T,interval="prediction")$fit[,2]))
litt.se <- as.vector(exp(predict.lm(fit.exp,se.fit=T,interval="prediction")$se.fit))
litt.pred.dat <- data.frame(f.age,litt.pred,litt.lo,litt.up,litt.se)
litt.pred.sort <- litt.pred.dat[order(litt.pred.dat[,1],decreasing=F),]
lines(litt.pred.sort$f.age,litt.pred.sort$litt.pred,col="red",lty=2,lwd=2)
lines(litt.pred.sort$f.age,litt.pred.sort$litt.lo,col="red",lty=3,lwd=2)
lines(litt.pred.sort$f.age,litt.pred.sort$litt.up,col="red",lty=3,lwd=2)

# fit linear model
fit.lin <- lm(f.litter ~ f.age)
fit.lin.summ <- summary(fit.lin)
litt.pred.lin <- (coef(fit.lin)[1] + f.age*coef(fit.lin)[2])
litt.up.lin <- as.vector(predict.lm(fit.lin,se.fit=T,interval="prediction")$fit[,3])
litt.lo.lin <- as.vector(predict.lm(fit.lin,se.fit=T,interval="prediction")$fit[,2])
litt.se.lin <- as.vector(predict.lm(fit.lin,se.fit=T,interval="prediction")$se.fit)
litt.pred.lin.dat <- data.frame(f.age,litt.pred.lin,litt.lo.lin,litt.up.lin,litt.se.lin)
litt.pred.lin.sort <- litt.pred.lin.dat[order(litt.pred.lin.dat[,1],decreasing=F),]
lines(litt.pred.lin.sort$f.age,litt.pred.lin.sort$litt.pred,col="red",lty=2,lwd=2)
lines(litt.pred.lin.sort$f.age,litt.pred.lin.sort$litt.lo,col="red",lty=3,lwd=2)
lines(litt.pred.lin.sort$f.age,litt.pred.lin.sort$litt.up,col="red",lty=3,lwd=2)

LL.vec <- as.vector(c(logLik(fit.exp)[1],logLik(fit.lin)[1]))
k.vec <- c(2,2)
AICc.vec <- -2*LL.vec + ((2*k.vec*length(f.litter))/(length(f.litter)-k.vec-1))
AIC.vec <- c(AIC(fit.exp), AIC(fit.lin))
dAICc.vec <- delta.IC(AICc.vec)
dAIC.vec <- delta.IC(AIC.vec)
wAICc.vec <- weight.IC(dAICc.vec)
wAIC.vec <- weight.IC(dAIC.vec)
AIC.table.out <- data.frame(k.vec,LL.vec,AIC.vec,dAIC.vec,wAIC.vec,AICc.vec,dAICc.vec,wAICc.vec)
colnames(AIC.table.out) <- c("k","LL","AIC","dAIC","wAIC","AICc","dAICc","wAICc")
rownames(AIC.table.out) <- c("exponential","linear")
AIC.table.out

## construct average fertility vector for matrix
pred.p.mat3 <- coef(fit.logp)[1] / (1+(age.vec/b.new)^coef(fit.logp)[3])
pred.p.mat4 <- ifelse(pred.p.mat3 > 1, 1, pred.p.mat3)
pred.p.mat5 <- ifelse(pred.p.mat4 < 0.001, 0, pred.p.mat4)
litt.pred2 <- exp(coef(fit.exp)[1] + age.vec*coef(fit.exp)[2])
f.fert.vec <- 0.5 * (pred.p.mat5*litt.pred2) # biennial breeding, so * 0.5 in all cases


## populate matrix
popmat[1, ] <- 0.5 * f.fert.vec # * 0.5 for female offspring only

# import Sx data
setwd("~/...") # set as appropriate
surv.dat <- read.table("Sx.lin.csv",header=T,sep=",")
surv.vec <- surv.dat$Sx
#surv.vec <- surv.dat$Sx.mod # reduce survival of older age classes following Snat from 18 years on
diag(popmat[2:(stages), ]) <- surv.vec[-stages]
popmat[stages,stages] <- 0
popmat.orig <- popmat ## save original matrix
popmat[1,]
popmat[1:10,1:10]
popmat[31:40,31:40]

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat,stages) # reproductive value
gen.l <- G.val(popmat,stages) # mean generation length
cat.pr <- 0.14/gen.l # probability of catastrophe (Reed et al. 2003)
  
## initial population vector
init.vec <- 1000*stable.stage.dist(popmat)
plot(age.vec,init.vec,xlab="age (yrs)", ylab="N", type="l")


#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 2017
#************************
yr.end <- 2030 # set projection end date
#************************
t <- (yr.end - yr.now)
popmat <- popmat.orig

## set population storage matrices
n.mat <- matrix(0,nrow=stages,ncol=(t+1))
n.mat[,1] <- init.vec

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i] 
}

## year projection vector
yrs <- seq(yr.now,yr.end,1)

# plot
plot(yrs,as.vector(colSums(n.mat)),type="l",xlab="year",ylab="N",xlim=c(yr.now,yr.end))



#####################################################################################
## stochastic projection 
## resampling proportion breeding & litter size for iterated fertility vector values
## set x% SD on Sx values (beta distribution)
## catastrophe sampler
#####################################################################################

## initial population vector
start.pop <- 1000 # start population size
init.vec <- start.pop*stable.stage.dist(popmat.orig)

## set time limit for projection in 1-yr increments
yr.now <- 2017
#************************
yr.end <- 2100 # set projection end date
#************************

setwd("~/...") # set as appropriate
catch <- read.table("structure.csv",header=T,sep=",")
catch.f <- subset(catch, sex="F")
catch.m <- subset(catch, sex="M")

par(mfrow=c(1,2))
hist(catch.f$TL,main="",xlab="length (mm)",br=30,col="black",border="white")
hist(catch.m$TL,main="",xlab="length (mm)",br=30,col="black",border="white")
par(mfrow=c(1,1))

## growth analysis
# 3-parameter logistic
Lm.mal <- 317; Lm.mal.mx <- 357; Lm.mal.mn <- 290
Lm.fem <- 308; Lm.fem.mx <- 321; Lm.fem.mn <- 298
Lm <- 306; Lm.mx <- 318; Lm.mn <- 296
alph.mal <- 78; alph.mal.mx <- 82; alph.mal.mn <- 75
alph.fem <- 74; alph.fem.mx <- 76; alph.fem.mn <- 72
alph <- 75; alph.mx <- 77; alph.mn <- 73
k.mal <- 0.126; k.mal.mx <- 0.14; k.mal.mn <- 0.11
k.fem <- 0.146; k.fem.mx <- 0.15; k.fem.mn <- 0.14
k <- 0.143; k.mx <- 0.15; k.mn <- 0.13

## age equations (from Drew et al. 2016 MFR)
par(mfrow=c(1,2))
age.f <- (log((catch.f$TL/10)*(Lm.fem-alph.fem))-log(alph.fem)-log(Lm.fem-(catch.f$TL/10)))/k.fem
hist(age.f, main="female")
age.m <- (log((catch.m$TL/10)*(Lm.mal-alph.mal))-log(alph.mal)-log(Lm.mal-(catch.m$TL/10)))/k.mal
hist(age.m, main="male")
par(mfrow=c(1,1))

## frequency of individuals per (yearly) age class
age.yr.f <- table(as.integer(age.f))
age.yr.m <- table(as.integer(age.m))

## proportional take for females
p.age.yr.f <- age.yr.f/sum(age.yr.f)
ages.caught <- as.numeric(names(p.age.yr.f))
prop.caught <- as.numeric(p.age.yr.f)

## set take function based on age structure of caught females
plot(ages.caught,prop.caught,pch=19,xlab="age",ylab="proportion caught")

## year projection vector
t <- (yr.end - yr.now)
yrs <- seq(yr.now,yr.end,1)

## density-feedback function
K <- 3*start.pop # set carrying capacity
# reduction in fertility
fert.min.mult <- 0.7
i.K <- start.pop:K
fert.mult <- 1 - cumsum(rep((1 - fert.min.mult)/length(i.K), length(i.K)))
#plot(i.K, fert.mult, pch=19, type="l")
i.K.fit <- lm(fert.mult ~ i.K)

## set population storage matrices
n.mat <- matrix(0,nrow=stages,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.orig

## iterate projection
iter <- 1000
itdiv <- iter/100

# set storage matrices & vectors
n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))
r.mat <- matrix(data = 0, nrow = iter, ncol = t)
age.wm.vec <- TL.wm.vec <- rep(0, iter)

for (e in 1:iter) {
  
  ## set value storage vectors
  r.stoch <- rep(0,t)
  
  ## reset popmat to original values
  popmat <- popmat.orig
  
  ## set up projection loop
  for (i in 1:t) {
  
    ## fertility sampler
    # litter size (rnorm sampler)
    litt.pred.stoch <- round(exp(rnorm(1,fit.exp.summ$coeff[1],fit.exp.summ$coefficients[3]) + age.vec*rnorm(1,fit.exp.summ$coeff[2],fit.exp.summ$coefficients[4])),0) 
    
    # proportion mature
    p.mat.stoch.a <- rnorm(1, fit.logp.summ$coefficients[1], fit.logp.summ$coefficients[4]) / (1+(age.vec/rnorm(1, fit.logp.summ$coefficients[2], fit.logp.summ$coefficients[5]))^rnorm(1, fit.logp.summ$coefficients[3], fit.logp.summ$coefficients[6]))
    p.mat.stoch.b <- ifelse(p.mat.stoch.a > 1, 1, p.mat.stoch.a)
    p.mat.stoch <- ifelse(p.mat.stoch.b < 0.001, 0, p.mat.stoch.b)
    p.mat.stoch
    
    # density feedback in fertility vector
    fert.multiplier <- ifelse(sum(n.mat[,i]) >= start.pop, as.numeric(coef(i.K.fit)[1] + coef(i.K.fit)[2]*sum(n.mat[,i])), 1)
    
    # construct stochastic fertility vector
    f.fert.stoch <- 0.5 * (p.mat.stoch*litt.pred.stoch) * fert.multiplier # biennial breeding, so *0.5 in all cases
    f.fert.stoch
    
    ## survival beta sampler
    # set SD for Sx
    Sx.sd <- 0.05 # can set to any value
    Sx.alpha <- estBetaParams(surv.vec, Sx.sd^2)$alpha
    Sx.beta <- estBetaParams(surv.vec, Sx.sd^2)$beta
    Sx.stoch <- rep(0,stages)
    for (x in 1:stages) {
      Sx.stoch[x] <- rbeta(1,Sx.alpha[x],Sx.beta[x])
    }
    
    ## reconstruct popmat with stochastic elements
    # fertility
    popmat[1, ] <- f.fert.stoch
    
    #survival (+ catastrophic mortality at 50%)
    catastrophe <- rbinom(1, 1, cat.pr)
    if (catastrophe == 1) {
      diag(popmat[2:(stages), ]) <- (0.5*Sx.stoch[-stages])}
    if (catastrophe == 0) {
      diag(popmat[2:(stages), ]) <- Sx.stoch[-stages]}
    popmat[stages,stages] <- 0
    
    ## project over interval
    n.mat[,i+1] <- popmat %*% n.mat[,i] 
    
    ## save r for this iteration' stochastic matrix
    r.running <- log(sum(n.mat[,i+1], na.rm=T) / sum(n.mat[,i], na.rm=T))
    r.stoch[i] <- ifelse(r.running == -Inf, NA, r.running)
    
    #print(i)
  }
  
  r.mat[e,] <- r.stoch
  n.sums.mat[e,] <- as.vector(colSums(n.mat))
  
  # median age & size of final population
  age.wm.vec[e] <- weighted.mean(x=age.vec, w=n.mat[,(t+1)], na.rm=T)
  TL.wm.vec[e] <- 308 / (1 + (3.162162 * exp(-0.146 * age.wm.vec[e]))) # predict TL from weighted mean age

  if (e %% itdiv==0) print(e) 
}

# N confidence limits
par(mfrow=c(1,2))
n.up <- n.lo <- n.mn <- rep(0,(t+1))
for (q in 1:(t+1)) {
  n.mn[q] <- mean(n.sums.mat[,q])
  n.up[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.975))
  n.lo[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.025))
}
plot(yrs, n.mn, type="l", xlab="year",ylab="N",xlim=c(yr.now,yr.end),ylim=c(min(n.lo),max(n.up)))
lines(yrs, n.up, lty=2, col="red")
lines(yrs, n.lo, lty=2, col="red")

# r confidence limits
r.up <- r.lo <- r.mn <- rep(0,(t))
for (q in 1:(t)) {
  r.mn[q] <- mean(r.mat[,q])
  r.up[q] <- as.vector(quantile(r.mat[,q], probs=0.975))
  r.lo[q] <- as.vector(quantile(r.mat[,q], probs=0.025))
}
plot(yrs[-1], r.mn, type="l", xlab="year",ylab="r",xlim=c(yr.now,yr.end),ylim=c(min(r.lo),max(r.up)))
lines(yrs[-1], r.up, lty=2, col="red")
lines(yrs[-1], r.lo, lty=2, col="red")
abline(h=0,lty=3,lwd=2,col="grey")
par(mfrow=c(1,1))

# weighted mean age confidence limits
age.wm.mn <- mean(age.wm.vec, na.rm=T)
age.wm.lo <- quantile(age.wm.vec, probs=0.025, na.rm=T)
age.wm.up <- quantile(age.wm.vec, probs=0.975, na.rm=T)
print(c(age.wm.lo, age.wm.mn, age.wm.up))

# weighted mean TL confidence limits
TL.wm.mn <- mean(TL.wm.vec, na.rm=T)
TL.wm.lo <- quantile(TL.wm.vec, probs=0.025, na.rm=T)
TL.wm.up <- quantile(TL.wm.vec, probs=0.975, na.rm=T)
print(c(TL.wm.lo, TL.wm.mn, TL.wm.up))






## natural mortality
## age-independent equations proposed by Pauly (1980), Chen and Yuan (2006), Jensen’s (1996) k-invariant method
## and Jensen’s (1996) maturity method as modified for use with sharks by Hisano et al. (2011).
## Two age-dependentmethods were also applied following Chen and Watanabe (1989) and Peterson and Wroblewski (1984)

# 3-param von Bertalanffy
VB3.L.inf.fem <- 457
VB3.k.fem <- 0.034
VB3.T.fem <- 17 # average annual sea surface temperature for stock

# Pauly (1980)
# ln(M) = −0.0066 − 0.279 ∗ ln(L∞) + 0.6743 ∗ ln(k) + 0.4634 ∗ ln(T)
mort.nat1 <- exp(-0.0066 - 0.279 * log(VB3.L.inf.fem) + 0.6743 * log(VB3.k.fem) + 0.4634 * log(VB3.T.fem))
surv.nat1 <- 1 - mort.nat1

# Jensen (1996)
# M = 1.6*k
mort.nat3 <- 1.6 * VB3.k.fem
surv.nat3 <- 1 - mort.nat3

L.pc.fem <- (catch.f$TL + 35.924)/1.364
fem.wt <- 6.71*10^-6 * L.pc.fem^3.14
plot(age.f, fem.wt/1000, pch=19, type="l", xlab="age (years)", ylab="mass (kg)")

fem.drywt <- 0.2*fem.wt
fem.Mwt <- 1.92*fem.drywt^-0.25 
plot(age.f, fem.Mwt, pch=19, type="l", xlab="age (years)", ylab="mortality")
fem.St <- exp(-fem.Mwt) # Smart et al. 2017
out.dat <- data.frame(age.f,fem.St)
plot(age.f, fem.St, pch=19, cex=0.4, xlab="age (years)", ylab="survival",ylim=c(0.5,1))
abline(h=surv.nat1, lty=2, lwd=1)
abline(h=surv.nat3, lty=2, lwd=1)

# logarithmic growth function y = a x^b 
age.St.dat <- na.omit(data.frame(age.f, fem.St))
age.St.dat <- age.St.dat[-dim(age.St.dat[1]),]
age.St.dat2 <- subset(age.St.dat, age.f > 0)
colnames(age.St.dat2) <- c("age","St")
fit.logglin <- lm(log(age.St.dat2$St) ~ log(age.St.dat2$age))
summary(fit.logglin)
param.init <- c(as.numeric(coef(fit.logglin)[1]), exp(as.numeric(coef(fit.logglin)[2])))
fit.logg <- nls((St ~ a * (age^b)), 
                data = age.St.dat2,
                algorithm = "port",
                start = c(a = param.init[1], b = param.init[2]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.logg.summ <- summary(fit.logg)
fit.logg.summ
age.vec.cont <- seq(1,max(age.vec),0.02)
pred.fem.St <- (coef(fit.logg)[1] * age.vec.cont^coef(fit.logg)[2])
lines(age.vec.cont,pred.fem.St,lty=2,lwd=3,col="red")


## model formula (after McCool)
## Von Bertalanffy growth function: St = Smax(1-b*exp(-kx))
## Smax = maximum survival
## b = flexibility parameter
## k = rate constant per year

## Von Bertalanffy growth 
fit.surv.growth <- nls(St ~ Smax*(1 - b * (exp(-k * age))), 
               data = age.St.dat2,
               algorithm = "port",
               start = c(Smax = max(age.St.dat2$St, na.rm=T), b = 1, k = 0.2),
               lower = c(0.5,0,0),
               trace = TRUE,      
               nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
pred.fem.St.vb <- coef(fit.surv.growth)[1] * (1 - coef(fit.surv.growth)[2] * (exp(-coef(fit.surv.growth)[3] * age.vec.cont)))
lines(age.vec.cont,pred.fem.St.vb,lty=2,lwd=3,col="blue")

# new max St based on average of all natural mortality maxs
St.max.mean <- mean(surv.nat1, surv.nat3, max(age.St.dat2$St))
age.vec.x <- seq(0, max(age.vec), 1)
pred.fem.St.mean <- St.max.mean * (1 - coef(fit.surv.growth)[2] * (exp(-coef(fit.surv.growth)[3] * age.vec.x)))
lines(age.vec.x,pred.fem.St.mean,lty=2,lwd=3,col="grey")
length(pred.fem.St.mean)
nat.s.out <- data.frame(age.vec.x,pred.fem.St.mean)

## new projection matrix based on mean, age-based natural mortality vector
# import Sx data
popmat.nat.orig <- popmat.orig

diag(popmat.nat.orig[2:(stages), ]) <- pred.fem.St.mean[-stages]
popmat.nat.orig[stages,stages] <- 0
popmat.nat.orig[1,]
popmat.nat.orig[1:10,1:10]
popmat.nat.orig[31:40,31:40]

## matrix properties
max.lambda(popmat.nat.orig) ## 1-yr lambda
max.r(popmat.nat.orig) # rate of population change, 1-yr
stable.stage.dist(popmat.nat.orig) ## stable stage distribution
R.val(popmat.nat.orig,stages) # reproductive value
gen.l <- G.val(popmat.nat.orig,stages) # mean generation length
cat.pr <- 0.14/gen.l # probability of catastrophe (Reed et al. 2003)

## initial population vector
init.vec <- 1000*stable.stage.dist(popmat.nat.orig)
plot(age.vec,init.vec,xlab="age (yrs)", ylab="N", type="l")

#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 2017
#************************
yr.end <- 2030 # set projection end date
#************************
t <- (yr.end - yr.now)
popmat <- popmat.orig

## set population storage matrices
n.mat <- matrix(0,nrow=stages,ncol=(t+1))
n.mat[,1] <- init.vec

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat.nat.orig %*% n.mat[,i] 
}

## year projection vector
yrs <- seq(yr.now,yr.end,1)

# plot
plot(yrs,as.vector(colSums(n.mat)),type="l",xlab="year",ylab="N",xlim=c(yr.now,yr.end))


#####################################################################################
## stochastic projection with natural mortality survival vector
## resampling proportion breeding & litter size for iterated fertility vector values
## set x% SD on Sx values (beta distribution)
## catastrophe sampler
#####################################################################################

## initial population vector
start.pop <- 1000 # start population size
init.vec <- start.pop*stable.stage.dist(popmat.nat.orig)

## set time limit for projection in 1-yr increments
yr.now <- 2017
#************************
yr.end <- 2100 # set projection end date
#************************

## year projection vector
t <- (yr.end - yr.now)
yrs <- seq(yr.now,yr.end,1)

## density-feedback function
K <- 3*start.pop # set carrying capacity
# reduction in fertility
fert.min.mult <- 0.7
i.K <- start.pop:K
fert.mult <- 1 - cumsum(rep((1 - fert.min.mult)/length(i.K), length(i.K)))
i.K.fit <- lm(fert.mult ~ i.K)

## set population storage matrices
n.mat <- matrix(0,nrow=stages,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.nat.orig

## iterate projection
iter <- 1000
itdiv <- iter/100

# set storage matrices & vectors
n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))
r.mat <- matrix(data = 0, nrow = iter, ncol = t)
age.wm.vec <- TL.wm.vec <- rep(0, iter)

for (e in 1:iter) {
  
  ## set value storage vectors
  r.stoch <- rep(0,t)
  
  ## reset popmat to original values
  popmat <- popmat.nat.orig
  
  ## set up projection loop
  for (i in 1:t) {
    
    ## fertility sampler
    # litter size (rnorm sampler)
    litt.pred.stoch <- round(exp(rnorm(1,fit.exp.summ$coeff[1],fit.exp.summ$coefficients[3]) + age.vec*rnorm(1,fit.exp.summ$coeff[2],fit.exp.summ$coefficients[4])),0) 
    
    # proportion mature
    p.mat.stoch.a <- rnorm(1, fit.logp.summ$coefficients[1], fit.logp.summ$coefficients[4]) / (1+(age.vec/rnorm(1, fit.logp.summ$coefficients[2], fit.logp.summ$coefficients[5]))^rnorm(1, fit.logp.summ$coefficients[3], fit.logp.summ$coefficients[6]))
    p.mat.stoch.b <- ifelse(p.mat.stoch.a > 1, 1, p.mat.stoch.a)
    p.mat.stoch <- ifelse(p.mat.stoch.b < 0.001, 0, p.mat.stoch.b)
    p.mat.stoch
    
    # density feedback in fertility vector
    fert.multiplier <- ifelse(sum(n.mat[,i]) >= start.pop, as.numeric(coef(i.K.fit)[1] + coef(i.K.fit)[2]*sum(n.mat[,i])), 1)
    
    # construct stochastic fertility vector
    f.fert.stoch <- 0.5 * (p.mat.stoch*litt.pred.stoch) * fert.multiplier # biennial breeding, so *0.5 in all cases
    
    ## survival beta sampler
    # set SD for Sx
    Sx.sd <- 0.05 # can set to any value
    Sx.alpha <- estBetaParams(pred.fem.St.mean, Sx.sd^2)$alpha
    Sx.beta <- estBetaParams(pred.fem.St.mean, Sx.sd^2)$beta
    Sx.stoch <- rep(0,stages)
    for (x in 1:stages) {
      Sx.stoch[x] <- rbeta(1,Sx.alpha[x],Sx.beta[x])
    }
    
    ## reconstruct popmat with stochastic elements
    # fertility
    popmat[1, ] <- f.fert.stoch
    
    #survival (+ catastrophic mortality at 50%)
    catastrophe <- rbinom(1, 1, cat.pr)
    if (catastrophe == 1) {
      diag(popmat[2:(stages), ]) <- (0.5*Sx.stoch[-stages])}
    if (catastrophe == 0) {
      diag(popmat[2:(stages), ]) <- Sx.stoch[-stages]}
    popmat[stages,stages] <- 0
    
    ## project over interval
    n.mat[,i+1] <- popmat %*% n.mat[,i] 
    
    ## save r for this iteration' stochastic matrix
    r.stoch[i] <- max.r(popmat)
    
    #print(i)
  }
  
  r.mat[e,] <- r.stoch
  n.sums.mat[e,] <- as.vector(colSums(n.mat))
  
  # median age & size of final population
  age.wm.vec[e] <- weighted.mean(x=age.vec, w=n.mat[,(t+1)], na.rm=T)
  TL.wm.vec[e] <- 308 / (1 + (3.162162 * exp(-0.146 * age.wm.vec[e]))) # predict TL from weighted mean age
  
  if (e %% itdiv==0) print(e) 
}

# N confidence limits
par(mfrow=c(1,2))
n.up <- n.lo <- n.mn <- rep(0,(t+1))
for (q in 1:(t+1)) {
  n.mn[q] <- mean(n.sums.mat[,q])
  n.up[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.975))
  n.lo[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.025))
}
plot(yrs, n.mn, type="l", xlab="year",ylab="N",xlim=c(yr.now,yr.end),ylim=c(min(n.lo),max(n.up)))
lines(yrs, n.up, lty=2, col="red")
lines(yrs, n.lo, lty=2, col="red")

# lambda confidence limits
r.up <- r.lo <- r.mn <- rep(0,(t))
for (q in 1:(t)) {
  r.mn[q] <- mean(r.mat[,q])
  r.up[q] <- as.vector(quantile(r.mat[,q], probs=0.975))
  r.lo[q] <- as.vector(quantile(r.mat[,q], probs=0.025))
}
plot(yrs[-1], r.mn, type="l", xlab="year",ylab="r",xlim=c(yr.now,yr.end),ylim=c(min(r.lo),max(r.up)))
lines(yrs[-1], r.up, lty=2, col="red")
lines(yrs[-1], r.lo, lty=2, col="red")
abline(h=0,lty=3,lwd=2,col="grey")
par(mfrow=c(1,1))

# weighted mean age confidence limits
age.wm.mn <- mean(age.wm.vec, na.rm=T)
age.wm.lo <- quantile(age.wm.vec, probs=0.025, na.rm=T)
age.wm.up <- quantile(age.wm.vec, probs=0.975, na.rm=T)
print(c(age.wm.lo, age.wm.mn, age.wm.up))

# weighted mean TL confidence limits
TL.wm.mn <- mean(TL.wm.vec, na.rm=T)
TL.wm.lo <- quantile(TL.wm.vec, probs=0.025, na.rm=T)
TL.wm.up <- quantile(TL.wm.vec, probs=0.975, na.rm=T)
print(c(TL.wm.lo, TL.wm.mn, TL.wm.up))




####################################################
## harvest scenarios in increments of current
####################################################

#####################################################################################
## stochastic projection 
## resampling proportion breeding & litter size for iterated fertility vector values
## set x% SD on Sx values (beta distribution)
## catastrophe sampler
#####################################################################################

## initial population vector
start.pop <- 30000 # start population size
init.vec <- start.pop*stable.stage.dist(popmat.orig)

## use natural-mortality popmat for set-up parameters
gen.l <- G.val(popmat.nat.orig,stages) # mean generation length
cat.pr <- 0.14/gen.l # probability of catastrophe (Reed et al. 2003)

## set time limit for projection in 1-yr increments
yr.now <- 2017
#************************
# run for 3 generations
gen3 <- round((3*gen.l),0)
yr.end <- yr.now + gen3
#yr.end <- 2100 # set projection end date
#************************

## do not count r/n in projections less than 1 generation (burn-in)
gen1 <- round(gen.l,0)

setwd("~/...")
catch <- read.table("structure.csv",header=T,sep=",")
catch.f <- subset(catch, sex="F")
catch.m <- subset(catch, sex="M")

par(mfrow=c(1,2))
hist(catch.f$TL,main="",xlab="length (mm)",br=30,col="black",border="white")
hist(catch.m$TL,main="",xlab="length (mm)",br=30,col="black",border="white")
par(mfrow=c(1,1))

## growth analysis
# 3-parameter logistic
Lm.mal <- 317; Lm.mal.mx <- 357; Lm.mal.mn <- 290
Lm.fem <- 308; Lm.fem.mx <- 321; Lm.fem.mn <- 298
Lm <- 306; Lm.mx <- 318; Lm.mn <- 296
alph.mal <- 78; alph.mal.mx <- 82; alph.mal.mn <- 75
alph.fem <- 74; alph.fem.mx <- 76; alph.fem.mn <- 72
alph <- 75; alph.mx <- 77; alph.mn <- 73
k.mal <- 0.126; k.mal.mx <- 0.14; k.mal.mn <- 0.11
k.fem <- 0.146; k.fem.mx <- 0.15; k.fem.mn <- 0.14
k <- 0.143; k.mx <- 0.15; k.mn <- 0.13

## age equations (from Drew et al. 2016 MFR)
par(mfrow=c(1,2))
age.f <- (log((catch.f$TL/10)*(Lm.fem-alph.fem))-log(alph.fem)-log(Lm.fem-(catch.f$TL/10)))/k.fem
hist(age.f, main="female")
age.m <- (log((catch.m$TL/10)*(Lm.mal-alph.mal))-log(alph.mal)-log(Lm.mal-(catch.m$TL/10)))/k.mal
hist(age.m, main="male")
par(mfrow=c(1,1))

## frequency of individuals per (yearly) age class
age.yr.f <- table(as.integer(age.f))
age.yr.m <- table(as.integer(age.m))

## proportional take for females
p.age.yr.f <- age.yr.f/sum(age.yr.f)
ages.caught <- as.numeric(names(p.age.yr.f))
prop.caught <- as.numeric(p.age.yr.f)

# Gamma function
# y = A*x^k * exp(−x*c)
caught.dat <- data.frame(ages.caught, prop.caught)
param.init <- c(0.05, 5, 0.5)
fit.gam <- nls(prop.caught ~ A * ages.caught^k * exp(-ages.caught * c), 
                data = caught.dat,
                algorithm = "port",
                start = c(A = param.init[1], k = param.init[2], c = param.init[3]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
plot(ages.caught, prop.caught, pch=19, ylim=c(0,0.17))
age.vec.cont <- seq(0,40,0.02)
pred.pcaught <- coef(fit.gam)[1] * age.vec.cont^coef(fit.gam)[2] * exp(-age.vec.cont * coef(fit.gam)[3])
pred.pcaught <- 0.142 * age.vec.cont^1.2 * exp(-age.vec.cont * 0.419)
lines((age.vec.cont),pred.pcaught,lty=2,lwd=3,col="red")
pred.pcaught1.vec <- coef(fit.gam)[1] * (0:39)^1.2 * exp(-(0:39) * coef(fit.gam)[3])
pred.pcaught1.vec[1] <- prop.caught[1]
# standardise to sum to 1
pred.pcaught.vec <- pred.pcaught1.vec/sum(pred.pcaught1.vec)
plot(ages.caught, prop.caught, pch=19, ylim=c(0,0.17))
lines(0:39, pred.pcaught.vec, lty=2, lwd=3, col="red")
harv.prop.dat <- data.frame(ages.caught, prop.caught)
harv.prop.fit <- data.frame(age.vec, pred.pcaught.vec)


## year projection vector
t <- (yr.end - yr.now)
yrs <- seq(yr.now,yr.end,1)

## density-feedback function
K <- start.pop # set carrying capacity (assume in this case we are there)
# reduction in fertility
fert.min.mult <- 0.7
i.K <- (start.pop/3):start.pop
fert.mult <- 1 - cumsum(rep((1 - fert.min.mult)/length(i.K), length(i.K)))
plot(i.K, fert.mult, pch=19, type="l")
i.K.fit <- lm(fert.mult ~ i.K)

## set population storage matrices
n.mat <- matrix(0,nrow=stages,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.nat.orig

## iterate projection
iter <- 1000
itdiv <- iter/100

# set harvest loop
# harvest as multiples of age.yr.f, where 0*prop.caught = current harvest, 0.5*prop.caught = 50% more than now, 1*prop.caught = double now
harv.p.vec <- seq(0,1,0.01)
lharv <- length(harv.p.vec)
r.mn.harv <- r.mn.up <- r.mn.lo <- n.mn.harv <- n.mn.up <- n.mn.lo <- n.min.harv <- n.min.up <- n.min.lo <- age.wm.mn <- age.wm.lo <- age.wm.up <- TL.wm.mn <- TL.wm.lo <- TL.wm.up <- rep(0,lharv)

for (h in 1:lharv) {

  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))
  r.mat <- matrix(data = 0, nrow = iter, ncol = t)
  age.wm.vec <- TL.wm.vec <- rep(0, iter)
  
  for (e in 1:iter) {
    
    ## set value storage vectors
    r.stoch <- rep(0,t)
    
    ## reset popmat to original values
    popmat <- popmat.orig
    
    ## set up projection loop
    for (i in 1:t) {
      
      ## fertility sampler
      # litter size (rnorm sampler)
      litt.pred.stoch <- round(exp(rnorm(1,fit.exp.summ$coeff[1],fit.exp.summ$coefficients[3]) + age.vec*rnorm(1,fit.exp.summ$coeff[2],fit.exp.summ$coefficients[4])),0) 
      
      # proportion mature
      p.mat.stoch.a <- rnorm(1, fit.logp.summ$coefficients[1], fit.logp.summ$coefficients[4]) / (1+(age.vec/rnorm(1, fit.logp.summ$coefficients[2], fit.logp.summ$coefficients[5]))^rnorm(1, fit.logp.summ$coefficients[3], fit.logp.summ$coefficients[6]))
      p.mat.stoch.b <- ifelse(p.mat.stoch.a > 1, 1, p.mat.stoch.a)
      p.mat.stoch <- ifelse(p.mat.stoch.b < 0.001, 0, p.mat.stoch.b)
      p.mat.stoch
      
      # density feedback in fertility vector
      #fert.multiplier <- ifelse(sum(n.mat[,i]) >= start.pop, as.numeric(coef(i.K.fit)[1] + coef(i.K.fit)[2]*sum(n.mat[,i])), 1)
      fert.multiplier <- as.numeric(coef(i.K.fit)[1] + coef(i.K.fit)[2]*sum(n.mat[,i]))
      
      # construct stochastic fertility vector
      f.fert.stoch <- 0.5 * 0.5 * (p.mat.stoch*litt.pred.stoch) * fert.multiplier # biennial breeding, so *0.5 in all cases; daughters only, so *0.5
      f.fert.stoch
      
      ## survival beta sampler
      # set SD for Sx
      Sx.sd <- 0.05 # can set to any value
      Sx.alpha <- estBetaParams(surv.vec, Sx.sd^2)$alpha
      Sx.beta <- estBetaParams(surv.vec, Sx.sd^2)$beta
      Sx.stoch <- rep(0,stages)
      for (x in 1:stages) {
        Sx.stoch[x] <- rbeta(1,Sx.alpha[x],Sx.beta[x])
      }
      
      ## reconstruct popmat with stochastic elements
      # fertility
      popmat[1, ] <- f.fert.stoch
      
      #survival (+ catastrophic mortality at 50%)
      catastrophe <- rbinom(1, 1, cat.pr)
      if (catastrophe == 1) {
        diag(popmat[2:(stages), ]) <- (0.5*Sx.stoch[-stages])}
      if (catastrophe == 0) {
        diag(popmat[2:(stages), ]) <- Sx.stoch[-stages]}
      popmat[stages,stages] <- 0
      
      ## project over interval
      n.mat[,i+1] <- popmat %*% (n.mat[,i])
      
      ## harvest
      n.mat[,i+1] <- n.mat[,i+1] - (harv.p.vec[h] * pred.pcaught.vec * n.mat[,i]) # substitute i for 1 in last element for proportional vs. constant harvest, respectively
      n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
      
      ## save r for this iteration' stochastic matrix
      r.running <- log(sum(n.mat[,i+1], na.rm=T) / sum(n.mat[,i], na.rm=T))
      r.stoch[i] <- ifelse(r.running == -Inf, NA, r.running)
      
      #print(i)
    }
    
    r.mat[e,] <- r.stoch
    n.sums.mat[e,] <- as.vector(colSums(n.mat))
    
    # median age & size of final population
    age.wm.vec[e] <- weighted.mean(x=age.vec, w=n.mat[,(t+1)], na.rm=T)
    TL.wm.vec[e] <- 308 / (1 + (3.162162 * exp(-0.146 * age.wm.vec[e]))) # predict TL from weighted mean age
    
    # plot
    plot(yrs,as.vector(colSums(n.mat)),type="l",xlab="year",ylab="N",xlim=c(yr.now,yr.end))
    plot(yrs[-1],r.mat[e,],type="l",ylab="r", xlab="years")
    abline(h=0, col="red", lty=2)
    
    if (e %% itdiv==0) print(e) 
  }
  
  # 1 generation burn-in
  n.mn <- apply(n.sums.mat[,(gen1+1):(t+1)], MARGIN=2, mean, na.rm=T)
  n.up <- apply(n.sums.mat[,(gen1+1):(t+1)], MARGIN=2, quantile, probs=0.975, na.rm=T)
  n.lo <- apply(n.sums.mat[,(gen1+1):(t+1)], MARGIN=2, quantile, probs=0.025, na.rm=T)
  
  plot(yrs[(gen1+1):(t+1)], n.mn, type="l", xlab="year",ylab="N",xlim=c(yrs[(gen1+1)],yrs[(t+1)]),ylim=c(min(n.lo),max(n.up)))
  lines(yrs[(gen1+1):(t+1)], n.up, lty=2, col="red")
  lines(yrs[(gen1+1):(t+1)], n.lo, lty=2, col="red")
  
  # r confidence limits with burn-in
  r.mn <- apply(r.mat[,(gen1+1):t], MARGIN=2, mean, na.rm=T)
  r.up <- apply(r.mat[,(gen1+1):t], MARGIN=2, quantile, probs=0.975, na.rm=T)
  r.lo <- apply(r.mat[,(gen1+1):t], MARGIN=2, quantile, probs=0.025, na.rm=T)
  
  # plot with burn-in
  plot(yrs[(gen1+1):(t)], r.mn, type="l", xlab="year",ylab="r",xlim=c(yrs[(gen1)],yrs[(t)]),ylim=c(min(r.lo),max(r.up)))
  lines(yrs[(gen1+1):(t)], r.up, lty=2, col="red")
  lines(yrs[(gen1+1):(t)], r.lo, lty=2, col="red")
  abline(h=0,lty=3,lwd=2,col="grey")
  
  # store values per h iteration
  r.mn.harv[h] <- median(r.mn, na.rm=T)
  r.mn.up[h] <- quantile(r.mn, probs=0.975, na.rm=T)
  r.mn.lo[h] <- quantile(r.mn, probs=0.025, na.rm=T)
  n.mn.harv[h] <- median(n.mn, na.rm=T)
  n.mn.up[h] <- quantile(n.mn, probs=0.975, na.rm=T)
  n.mn.lo[h] <- quantile(n.mn, probs=0.025, na.rm=T)
  n.min.harv[h] <- median(n.lo, na.rm=T)
  n.min.up[h] <- quantile(n.lo, probs=0.975, na.rm=T)
  n.min.lo[h] <- quantile(n.lo, probs=0.025, na.rm=T)
  age.wm.mn[h] <- mean(age.wm.vec, na.rm=T)
  age.wm.up[h] <- quantile(age.wm.vec, probs=0.975, na.rm=T)
  age.wm.lo[h] <- quantile(age.wm.vec, probs=0.025, na.rm=T)
  TL.wm.mn[h] <- mean(TL.wm.vec, na.rm=T)
  TL.wm.up[h] <- quantile(TL.wm.vec, probs=0.975, na.rm=T)
  TL.wm.lo[h] <- quantile(TL.wm.vec, probs=0.025, na.rm=T)
  
  print(paste("harvest increment = ", h, sep=""))
}

par(mfrow=c(3,1))
plot(harv.p.vec, r.mn.harv, type="l", lwd=2, ylim=c(-.1,.03), xlab="", ylab="mean long-term r")
abline(h=0, lty=2, col="red")
lines(harv.p.vec, r.mn.lo, lty=2, lwd=1)
lines(harv.p.vec, r.mn.up, lty=2, lwd=1)

plot(harv.p.vec, n.mn.harv, type="l", lwd=2, ylim=c(min(n.mn.lo),max(n.mn.up)), xlab="", ylab="mean long-term N")
lines(harv.p.vec, n.mn.lo, lty=2, lwd=1)
lines(harv.p.vec, n.mn.up, lty=2, lwd=1)

plot(harv.p.vec, n.min.harv, type="l", lwd=2, ylim=c(min(n.min.lo),max(n.min.up)), xlab="harvest multplier on current rate", ylab="mean long-term N min")
lines(harv.p.vec, n.min.lo, lty=2, lwd=1)
lines(harv.p.vec, n.min.up, lty=2, lwd=1)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
plot(harv.p.vec, age.wm.mn, type="l", lwd=2, ylim=c(min(age.wm.lo),max(age.wm.up)), ylab="wmn age (yrs) of final pop", xlab="harvest multiplier")
lines(harv.p.vec, age.wm.lo, lty=2, lwd=1)
lines(harv.p.vec, age.wm.up, lty=2, lwd=1)

plot(harv.p.vec, TL.wm.mn, type="l", lwd=2, ylim=c(min(TL.wm.lo),max(TL.wm.up)), ylab="wmn TL (cm) of final pop", xlab="harvest multiplier")
lines(harv.p.vec, TL.wm.lo, lty=2, lwd=1)
lines(harv.p.vec, TL.wm.up, lty=2, lwd=1)
par(mfrow=c(1,1))


