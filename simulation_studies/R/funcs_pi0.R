## ##############################################################################
## START CODE FROM EXTERNAL SOURCE
## source: https://github.com/SiminaB/Fdr-regression/blob/71b123f/functions.R
## license: not listed
## date: 2017/10/17
## modified: Patrick Kimes
## ##############################################################################

##------Functions of covariates-------##
pi0_f1 <- function(x) {
    p2 <- -0.2
    p1 <- 1.2
    a <- 4/(p1-p2)^2
    
    y <- -a*(x-p1)*(x-p2)
    y[x >= 0.7] <- -a*(0.7-p1)*(0.7-p2)
    y[x <= (p1+p2)/2] <- 1  
    y
}
pi0_f2 <- function(x) {
    y <- rep(0, length=length(x))
    y[x >= 0.7] <- -2.5*(x[x >= 0.7]-0.7)^2  
    y
}
pi0f3 <- function(x) {
    y <- rep(0, length=length(x))
    y[x < 0.7] <- -(x[x < 0.7]-0.1)^2
    y[x >= 0.7] <- -(min(x[x >= 0.7])-0.1)^2
    y[x<=0.1] <- 0
    y
}

##smooth function of one covariate for different levels of second covariate
pi0_smooth2 <- function(x1, x2) {
    y1 <- pi0_f1(x1)
    y2 <- pi0_f2(x1)
    y3 <- pi0_f3(x1)
    y <- rep(0, length(x1))
    y[x2 == 1] <- y1[x2 == 1] + y2[x2 == 1] + 0.12*y3[x2 == 1]
    y[x2 == 2] <- y1[x2 == 2] + 0.5*y2[x2 == 2] + 0.06*y3[x2 == 2]
    y[x2 == 3] <- y1[x2 == 3] + 0.3*y2[x2 == 3] 
    y
}

##smooth function of a single variable
pi0_smooth1 <- function(x) {
    y1 <- pi0_f1(x)
    y2 <- pi0_f2(x)
    y3 <- pi0_f3(x)
    y <- rep(0, length(x))
    y <- y1 + y2 + 0.12*y3
    y
}

## ##############################################################################
## END CODE FROM EXTERNAL SOURCE
## ##############################################################################
