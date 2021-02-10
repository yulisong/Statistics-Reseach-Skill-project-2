#### Load packages
# Estimate probability densities using Gaussian finite mixture model 
library(mclust)
# Generate two dimensional gaussian distribution
library(mvtnorm)
# Compute two dimensional KDE
library(MASS)
# Compute two dimensional histogram and plot it
library(gplots)
# Estimate probability densities using smoothing spline ANOVA models
library(gss)
# For numeric integration
library(sfsmisc)




#### Define three customized function
# https://www.stat.cmu.edu/~larry/all-of-statistics/=Rprograms/density.examples.r
cv.hist.fun <- function(x){
  ### histogram cross validation function
  n <- length(x)
  a <- min(x)
  b <- max(x)
  k <- 100
  nbins <- seq(1,n,length=k)  ###number of bins
  nbins <- round(nbins)
  h <- (b-a)/nbins            ###width of bins
  risk <- rep(0,k)
  for(i in 1:k){
    ###get counts N_j
    br <- seq(a,b,length=nbins[i]+1)
    N <- hist(x,breaks=br,plot=F)$counts
    risk[i] <- sum(N^2)/(n^2*h[i])  - (2/(h[i]*n*(n-1)))*sum(N*(N-1))
  }
  hbest <- h[risk==min(risk)]
  hbest <- hbest[1]  ###in case of tie take first (smallest) one
  mbest <- (b-a)/hbest   ###optimal number of bins
  list(risk=risk, nbins=nbins, h=h, mbest=mbest)
}

dens.hist <- function(hist) {
  # Density of sample points inside hist() function
  dens = c()
  len = length(hist$density)
  for (i in 1:len) {
    dens = c(dens, rep(hist$density[i], hist$counts[i]))
  }
  return(dens)
}

f <- function (x) {
  0.7*dnorm(x, -10, 1) + 0.3*dnorm(x, 5, 0.5)
}




#### Data Generating Processes and Preliminary Experiments
### Standard Normal distribution
set.seed(101145)
# Number of samples
siz.sample <- 1000
# Generating sample data
sim.data <- rnorm(siz.sample, 0, 1)

# Calculate the optimal number of bins in the histogram
siz.hist <- cv.hist.fun(sim.data)$mbest
# Calculate the optimal bandwidth in the KDE
bw.kde <- bw.SJ(sim.data)
# Calculate the kernel density function with Gassian kernel used
sim.dens <- density(sim.data, bw = bw.kde)
# Calculate the real density function
real.x <- seq(-5, 5, length=1000)
real.y <- dnorm(real.x, mean=0, sd=1)
# Produce a density estimate using a Gaussian finite mixture model 
mod <- densityMclust(sim.data)
summary(mod)
# Check how the number of models influences BIC
plot(mod, what = "BIC")
title("How the number of models influences BIC")
plot(mod, what = "diagnostic", type = "cdf")
plot(mod, what = "diagnostic", type = "qq")
# Estimate probability densities using smoothing spline ANOVA models
fit <- ssden(~ sim.data)
# The sides of the hyper cube are specified by domain
domain <- fit$domain$sim.data
# Create a regular sequence inside the domain
xx <- seq(domain[1], domain[2], length = 101)

# Parametric estimation
para.est.mean <- mean(sim.data)
para.est.sd <- sd(sim.data)

hist(sim.data, probability=T, col="#EAECEE",  
     breaks= seq(min(sim.data), max(sim.data),l=siz.hist), 
     xlab="Data", ylab="Density", main="Normal Dist",
     cex.axis=1.5, cex.lab=1.5)
lines(sim.dens, lwd=2, col="#16A085")
lines(real.x, real.y, type="l", lwd=2, col="#884EA0")
lines(xx, dssden(fit, xx), type = "l",lwd=2, col="#FFBF00")
lines(sort(mod$data), mod$density[order(mod$data)], type="l", lwd=2, col="#6495ED")
legend(1.5, 0.4, legend=c("Hist", "KDE", "Original Dist", "Smoothing Spline", "Mix Gaussian"),
       col=c("#EAECEE", "#16A085", "#884EA0", "#FFBF00", "#6495ED"), lty=1, cex=0.6)



### Gamma(2,3) distribution
set.seed(101145)
siz.sample <- 1000
sim.data <- rgamma(siz.sample, shape = 2, scale = 3)

siz.hist <- cv.hist.fun(sim.data)$mbest
bw.kde <- bw.SJ(sim.data)
sim.dens <- density(sim.data, bw = bw.kde)
real.x <- seq(0, 30, length=1000)
real.y <- dgamma(real.x, shape = 2, scale = 3)
para.est.shape <- mean(sim.data)^2/var(sim.data)
para.est.scale <- var(sim.data)/mean(sim.data)
mod <- densityMclust(sim.data)
summary(mod)
plot(mod, what = "BIC")
plot(mod, what = "diagnostic", type = "cdf")
plot(mod, what = "diagnostic", type = "qq")
fit <- ssden(~ sim.data)
domain <- fit$domain$sim.data
xx <- seq(domain[1], domain[2], length = 101)

hist(sim.data, probability=T, col="#EAECEE",  
     breaks= seq(min(sim.data), max(sim.data),l=siz.hist), 
     xlab="Data", ylab="Density", main="Gamma Dist",
     cex.axis=1.5, cex.lab=1.5)
lines(sim.dens, lwd=2, col="#16A085")
lines(real.x, real.y, type="l", lwd=2, col="#884EA0")
lines(xx, dssden(fit, xx), type = "l",lwd=2, col="#FFBF00")
lines(sort(mod$data), mod$density[order(mod$data)], type="l", lwd=2, col="#6495ED")
legend(20, 0.14, legend=c("Hist", "KDE", "Original Dist", "Smoothing Spline", "Mix Gaussian"),
       col=c("#EAECEE", "#16A085", "#884EA0", "#FFBF00", "#6495ED"), lty=1, cex=0.6)



### A Mix distribution - 0.7*Norm(-10,1)+0.3*Norm(5,0.5)
set.seed(101145)
siz.sample <- 1000
sim.data <- c(rnorm(700, -10, 1), rnorm(300, 5, 0.5))

siz.hist <- cv.hist.fun(sim.data)$mbest
bw.kde <- bw.SJ(sim.data)
sim.dens <- density(sim.data, bw = bw.kde)
# Here we also calculate KDE under default bandwidth
# to check its performance of bw.SJ
sim.dens.default <- density(sim.data)
# Note that here we can still do the parametric estimation
# but it's a little bit complicated, so we don't code here
real.x <- seq(-20, 10, length=1000)
## density matrix
real.y <- f(real.x)
mod <- densityMclust(sim.data)
summary(mod)
plot(mod, what = "BIC")
plot(mod, what = "diagnostic", type = "cdf")
plot(mod, what = "diagnostic", type = "qq")
fit <- ssden(~ sim.data)
domain <- fit$domain$sim.data
xx <- seq(domain[1], domain[2], length = 101)

hist(sim.data, probability=T, col="#EAECEE",  
     breaks= seq(min(sim.data), max(sim.data),l=siz.hist), 
     xlab="Data", ylab="Density", main="Mixture Dist",
     cex.axis=1.5, cex.lab=1.5)
lines(sim.dens, lwd=2, col="#16A085")
lines(real.x, real.y, type="l", lwd=2, col="#884EA0")
lines(xx, dssden(fit, xx), type = "l",lwd=2, col="#FFBF00")
lines(sim.dens.default, lwd=2, col="#D5D8DC")
lines(sort(mod$data), mod$density[order(mod$data)], type="l", lwd=2, col="#6495ED")
legend(-2, 0.3, legend=c("Hist", "KDE", "KDE Default", "Original Dist", "Smoothing Spline", "Mix Gaussian"),
       col=c("#EAECEE", "#16A085", "#D5D8DC", "#884EA0", "#FFBF00", "#6495ED"), lty=1, cex=0.6)



### Some interesting Multivariate plot
### Generating samples inside a annulus
# https://stats.stackexchange.com/questions/406637/how-to-generate-a-data-sample-along-a-annulus-between-fixed-radii
set.seed(101145)
outer_radius <- 1
inner_radius <- 0.7
n <- 1000
rho <- sqrt(runif(n,inner_radius^2,outer_radius^2))
theta <- runif(n, 0, 2*pi)
x <- rho * cos(theta)
y <- rho * sin(theta)
# Draw the point plot in 2 dimension
plot(x, y, pch=19, cex=0.6, col="#00000020")

# Fit KDE to two dimensional data
f1 <- kde2d(x,y)
# Creates a grid of colored rectangles with colors corresponding
# to the number of samples inside the rectangles
image(f1)
# Generate 3D KDE plot of samples
persp(f1,phi = 30, theta = 30, expand=0.5, shade=0.5)
# Fit Histogram to two dimensional data
h2d <- hist2d(x,y, nbins=c(20,20))
# Generate 3D histogram of samples
persp( h2d$x, h2d$y, h2d$counts,
       ticktype="detailed", theta=30, phi=30,
       expand=0.5, shade=0.5, col="#EAECEE", ltheta=-30)
# Produce the density estimate using a Gaussian finite mixture model
mod <- densityMclust(data.frame(x,y))
plot(mod, what = "density", type = "hdr")
plot(mod, what = "density", type = "hdr",
     data = data.frame(x,y), points.cex = 0.5)
plot(mod, what = "density", col="#16A085", type = "persp", shade=0.1, theta=30, 
     phi=30, expand=0.5)



### Generating samples that have two peaks
set.seed(101145)
sim <- rbind(rmvnorm(n=500, mean=c(1,1)), rmvnorm(n=500, mean=c(10,10)))
x <- sim[, 1]
y <- sim[, 2]
# Draw the point plot in 2 dimension
plot(x, y, pch=19, cex=0.6, col="#00000020")

# Fit KDE to two dimensional data
f1 <- kde2d(x,y)
# Creates a grid of colored rectangles with colors corresponding
# to the number of samples inside the rectangles
image(f1)
# Generate 3D KDE plot of samples
persp(f1,phi = 30, theta = 30, expand=0.5, shade=0.5)
# Fit Histogram to two dimensional data
h2d <- hist2d(x,y, nbins=c(20,20))
# Generate 3D histogram of samples
persp( h2d$x, h2d$y, h2d$counts,
       ticktype="detailed", theta=30, phi=30,
       expand=0.5, shade=0.5, col="#EAECEE", ltheta=-30)
# Produce the density estimate using a Gaussian finite mixture model
mod <- densityMclust(data.frame(x,y))
plot(mod, what = "density", type = "hdr")
plot(mod, what = "density", type = "hdr",
     data = data.frame(x,y), points.cex = 0.5)
plot(mod, what = "density", type = "persp", col="#16A085", shade=0.1, theta=30, 
     phi=30, expand=0.5)




#### Monte Carlo Simulation
# https://stats.stackexchange.com/questions/390777/how-to-compute-integrated-squared-error-for-kernel-density-estimation-in-r
# Sample size=1000
sim = 1000
siz.sample = 1000



### Normal(0,1)
## KDE
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = rnorm(siz.sample, 0, 1) 
  #Computing of bandwidth
  bw.kde <- bw.SJ(sim.data)
  #Estimate
  sim.dens <- density(sim.data, bw = bw.kde)
  
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = sim.dens$x, (sim.dens$y - dnorm(sim.dens$x))^2)
}
c1 = sum(ise)/sim
d1 = sd(ise)

## Histogram
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = rnorm(siz.sample, 0, 1) 
  siz.hist <- cv.hist.fun(sim.data)$mbest
  his = hist(sim.data, probability=T, breaks= seq(min(sim.data), max(sim.data),l=siz.hist))
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = sort(sim.data), (dens.hist(his) - dnorm(sort(sim.data)))^2)
}
c2 = sum(ise)/sim
d2 = sd(ise)

## Mixture Gaussian 
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = rnorm(siz.sample, 0, 1) 
  mod <- densityMclust(sim.data)
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = sort(mod$data), (mod$density[order(mod$data)] - dnorm(sort(mod$data)))^2)
}
c3 = sum(ise)/sim
d3 = sd(ise)


## Smooth Spline
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = rnorm(siz.sample, 0, 1) 
  fit <- ssden(~ sim.data)
  domain <- fit$domain$sim.data
  xx <- seq(domain[1], domain[2], length = 101)
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = xx, (dssden(fit, xx) - dnorm(xx))^2)
}
c4 = sum(ise)/sim
d4 = sd(ise)



### Gamma(2,3)
## KDE 
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data <- rgamma(siz.sample, shape = 2, scale = 3)
  #Computing of bandwidth
  bw.kde <- bw.SJ(sim.data)
  #Estimate
  sim.dens <- density(sim.data, bw = bw.kde)
  
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = sim.dens$x, (sim.dens$y - dgamma(sim.dens$x, shape = 2, scale = 3))^2)
}
c5 = sum(ise)/sim
d5 = sd(ise)


## Histogram
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = rgamma(siz.sample, shape = 2, scale = 3)
  siz.hist <- cv.hist.fun(sim.data)$mbest
  his = hist(sim.data, probability=T, breaks= seq(min(sim.data), max(sim.data),l=siz.hist))
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = sort(sim.data), (dens.hist(his) - dgamma(sort(sim.data), shape = 2, scale = 3))^2)
}
c6 = sum(ise)/sim
d6 = sd(ise)


## Mixture Gaussian 
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = rgamma(siz.sample, shape = 2, scale = 3)
  mod <- densityMclust(sim.data)
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = sort(mod$data), (mod$density[order(mod$data)] - dgamma(sort(mod$data), shape = 2, scale = 3))^2)
}
c7 = sum(ise)/sim
d7 = sd(ise)


## Smooth Spline
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = rgamma(siz.sample, shape = 2, scale = 3)
  fit <- ssden(~ sim.data)
  domain <- fit$domain$sim.data
  xx <- seq(domain[1], domain[2], length = 101)
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = xx, (dssden(fit, xx) - dgamma(xx, shape = 2, scale = 3))^2)
}
c8 = sum(ise)/sim
d8 = sd(ise)

### Mixture distribution
## KDE
set.seed(101145)
ise = c()
for (i in 1:250) {
  sim.data = c(rnorm(siz.sample*0.7, -10, 1), rnorm(siz.sample*0.3, 5, 0.5))
  #Computing of bandwidth
  bw.kde <- bw.SJ(sim.data)
  #Estimate
  sim.dens <- density(sim.data, bw = bw.kde)
  
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = sim.dens$x, (sim.dens$y - (0.7*dnorm(sim.dens$x, -10, 1) + 0.3*dnorm(sim.dens$x, 5, 0.5)))^2)
}
c9 = sum(ise)/sim
d9 = sd(ise)


## Histogram
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = c(rnorm(siz.sample*0.7, -10, 1), rnorm(siz.sample*0.3, 5, 0.5))
  siz.hist <- cv.hist.fun(sim.data)$mbest
  his = hist(sim.data, probability=T, breaks= seq(min(sim.data), max(sim.data),l=siz.hist))
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = sort(sim.data), (dens.hist(his) - (0.7*dnorm(sort(sim.data), -10, 1) + 0.3*dnorm(sort(sim.data), 5, 0.5)))^2)
}
c10 = sum(ise)/sim
d10 = sd(ise)


## Mixture Gaussian 
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = c(rnorm(siz.sample*0.7, -10, 1), rnorm(siz.sample*0.3, 5, 0.5))
  mod <- densityMclust(sim.data)
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = sort(mod$data), (mod$density[order(mod$data)] - (0.7*dnorm(sort(mod$data), -10, 1) + 0.3*dnorm(sort(mod$data), 5, 0.5)))^2) # res 0.0006682952
}
c11 = sum(ise)/sim
d11 = sd(ise)


## Smooth Spline
set.seed(101145)
ise = c()
for (i in 1:sim) {
  sim.data = c(rnorm(siz.sample*0.7, -10, 1), rnorm(siz.sample*0.3, 5, 0.5))
  fit <- ssden(~ sim.data)
  domain <- fit$domain$sim.data
  xx <- seq(domain[1], domain[2], length = 101)
  # Simple numerical integration in pkg sfsminsc
  ise[i] = sfsmisc::integrate.xy(x = xx, (dssden(fit, xx) - (0.7*dnorm(xx, -10, 1) + 0.3*dnorm(xx, 5, 0.5)))^2)
}
c12 = sum(ise)/sim
d12 = sd(ise)
