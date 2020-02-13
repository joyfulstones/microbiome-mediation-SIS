library(stats)
library(Compositional)
library(nloptr)
library(glmnet)
library(ncvreg)
library(MASS)
library(Matrix)
#library(ccmm)
library(hommel)
library(microHIMA)
#source("microHIMA.R")
# ########### generation data
# ### Sample size and number of mediators
n <- 200
p <- 25


sample.num=n
otu.num=p

Treatment = rbinom(n,1,0.2)

######Two covariates
covariates=cbind(sample(c(1,0),sample.num,replace = TRUE),rnorm(sample.num))
#
# ### parameters
beta0 <- matrix(0,1,p)

beta0 <- as.numeric(beta0)
betaT=rep(0,otu.num)
betaT[c(1,2,3)]=c(1,1.2,1.5) # the first three are non-zero
betaX=matrix(0,otu.num,2)
#
alpha0=0
alphaT=1
alphaZ=alphaC=rep(0,otu.num)
alphaZ[c(1,2,3)]=c(1.3,-0.7, -0.6) # the first three are non-zero for response
alphaX=c(0.5,0.5)
#
#
# ############Microbiome data
library(dirmult)
X=cbind(rep(1,sample.num),covariates,Treatment) #n*(1+q+p)
b=cbind(beta0,betaX,betaT) #num.otu*(1+q+p)
gamma.simu=exp(X %*% t(b)) # n * num.otu
otu.com=t(apply(gamma.simu,1,rdirichlet,n=1))

##################Outcome data
X=cbind(rep(1,sample.num),Treatment,covariates,log(otu.com),log(otu.com)*Treatment)
b=c(alpha0,alphaT,alphaX,alphaZ,alphaC)

outcome=c(b%*%t(X)+rnorm(sample.num,mean = 0, sd =1))
exposure <- t(t(Treatment))

############  example


mhima.fit <- mhima(exposure,covariates,otu.com,outcome)

print(mhima.fit)




