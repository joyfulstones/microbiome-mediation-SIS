mhima <- function(X,M,Y){
library(hdi)
library("glmnet")
library(ncvreg)
library(scalreg)
library(MASS)
library(Matrix)
library(stats)
##
# #X <- ffq.raw$aofib / ffq.raw$calor
# #X <- t(t(sqrt(X)))  # Make it more normal 
# X <- ffq.raw$tfat
# #X <- scale(X)
# Y <- bmi.c  # BMI
# #Y <- scale(X)
# #############  from the Genus level
# # M <- t(abund.list$Genus)
# # n <- dim(M)[1]
# # p <- dim(M)[2]
# # M <- M+1
# # MN <- apply(M,2,mean)
# # for (i in 1:p){
# #   M[,i] <- log(M[,i]/MN[i])
# # }
# # M <- scale(M)
# #######  
# M <- t(otu.tab)
# #M <- t(abund.list$Genus)
# n <- dim(M)[1]
# d <- dim(M)[2]
# 
# MI <- matrix(0,n,d)
# for (i in 1:n){
#   MI[i,] <- t(as.numeric(M[i,] > 0))
# }
# 
# # SM <- length(which(M!= 0))
# # n <- dim(M)[1]
# Index <- matrix(0,1,d)
# for ( i in 1:d){
#    Index[i] <- length(which(M[,i] != 0))
# }
# IDS <- which(Index > n*0.08)
# 
# M <- M[,IDS]

n <- dim(M)[1]
d <- dim(M)[2]
##
for (i in 1:n){
  for (j in 1:d){
    if (M[i,j] == 0){
      M[i,j] <- 0.5
    }
  }
}

##
for ( i in 1:n){
  GM <- sum(M[i,])
  M[i,] <- M[i,]/GM   # composition
}




#################### the following is the SIS+MCP filter step

Index_S <- matrix(0,1,d)
P_b_raw <-  matrix(1,1,d)
P_a_raw <-  matrix(1,1,d)
alpha_est <- matrix(0,1,d)
alpha_SE  <- matrix(0,1,d)
beta_est <- matrix(0,1,d)
beta_SE <-  matrix(0,1,d)
M1 <- t(t(M[,1]))
for (k in 1:d){
  M[,1] <- M[,k]
  M[,k] <- M1 
  MT <- matrix(0,n,d-1)
  for (i in 1:n){
    for (j in 1:(d-1)){
      C_1 <- sqrt((d-j)/(d-j+1))
      C_2 <- prod(M[i,(j+1):d]^(1/(d-j)))
      MT[i,j] <- C_1*log(M[i,j]/C_2)
    }
  }
  MT <- matrix(as.numeric(scale(MT)),nrow(scale(MT)),ncol(scale(MT)))
  X <- matrix(as.numeric(scale(X)),nrow(scale(X)),ncol(scale(X)))
  MX <- cbind(MT,X)
  CN <- round(n/log(n))
  ## the SIS step
  P_S <- matrix(0,1,d)
  for (j in 1:d){
     lm.fit <- lm(Y~MX[,j])
     lm.out <- summary(lm.fit)
     a_est <- lm.out$coefficients[2,1]
     a_se  <- lm.out$coefficients[2,2]
     P_S[j]  <- 2*(1-pnorm(abs(a_est/a_se),0,1))
    
  }
  ID <- which(P_S <= sort(P_S)[CN])
  L = length(ID)
  if ((ID[1] == 1)&&(L > 0)){
    fit <- ncvreg(MX[,ID],Y,family = "gaussian",penalty="MCP",alpha=1)
    lam <- fit$lambda[which.min(BIC(fit))]
    Coefficients <- coef(fit, lambda=lam)
    EST <- Coefficients[2:(L+1)]
    if (EST[1] != 0){
      Index_S[k] = 1
      IS <- which(EST != 0)
      XF <- MX[,ID[IS]]
      lm.fit <- lm(Y~XF)
      lm.out <- summary(lm.fit)
      a_est <- lm.out$coefficients[2,1]
      a_se  <- lm.out$coefficients[2,2]
      P_b_raw[1,k]  <- 2*(1-pnorm(abs(a_est/a_se),0,1))
      beta_est[k] <- a_est
      beta_SE[k]  <- a_se
      ##
      lm.fit <- lm(MT[,1]~X)
      lm.out <- summary(lm.fit)
      a_est <- lm.out$coefficients[2,1]
      a_se  <- lm.out$coefficients[2,2]
      P_a_raw[1,k]  <- 2*(1-pnorm(abs(a_est/a_se),0,1))
      alpha_est[k] <- a_est
      alpha_SE[k] <- a_se
      
    }
    
  }
  
  ##
#print(k)
  
}

LF <- sum(Index_S)

P_joint_corr <- apply(rbind(P_a_raw,P_b_raw),2,max)*LF # the corrected p-values

id <- which(P_joint_corr <= 0.05)
##
# print(n)
# print(d)
# print(LF)
# print(id)
# print(P_joint_corr[id])
##
# print(alpha_est[id])
# print(alpha_SE[id])
# print(beta_est[id])
# print(beta_SE[id])
### the out-put results

out_result <- list(ID = id, alpha_hat = alpha_est[id], alpha_hat_SE= alpha_SE[id], beta_hat=beta_est[id], beta_hat_SE =beta_SE[id],Pval=P_joint_corr[id])
return(out_result)

}

