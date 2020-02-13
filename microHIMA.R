mhima <- function(exposure, covariates, otu.com, outcome){
  library(stats)
  library(Compositional)
  library(nloptr)
  library(glmnet)
  library(ncvreg)
  library(MASS)
  library(Matrix)
  library(stats)
  library(ccmm)
  library(hommel)
#   source(DLASSO_fun)
  #############################   begin the R package %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DLASSO_fun <- function(X,Y){
    library(glmnet)
    n <- dim(X)[1]
    p <- dim(X)[2]
    fit = glmnet(X, Y, alpha = 1)
    cv.fit <- cv.glmnet(X, Y, alpha = 1)
    beta_0 <- coef(fit, s = cv.fit$lambda.min)[2:(p+1)]
    #
    fit <- glmnet(X[,2:p], X[,1], alpha = 1)
    cv.fit <- cv.glmnet(X[,2:p], X[,1], alpha = 1)
    phi_hat <- coef(fit, s = cv.fit$lambda.min)[2:p]
    ##
    R <- X[,1] - X[,2:p]%*%t(t(phi_hat))
    E <- Y - X%*%t(t(beta_0))
    beta_1_hat <- beta_0[1] + sum(R*E)/sum(R*X[,1]) #  The de-biased lasso estimator
    ##
    sigma_e2 <- sum(E^2)/(n-length(which(beta_0!=0)))
    
    sigma_beta1_hat <- sqrt(sigma_e2)*sqrt(sum(R^2))/abs(sum(R*X[,1]))
    
    results <- c(beta_1_hat,sigma_beta1_hat)
    return(results)
    
  }
  
  
  ########################  the DLASSO method
  M_raw <- otu.com
  X <- cbind(exposure,covariates)
  Y <- outcome 
  Y <- Y - mean(Y)
  ###
  M <- M_raw
  n <- dim(M)[1]
  d <- dim(M)[2]
  Index_S <- matrix(0,1,d)
  P_b_raw <-  matrix(0,1,d)
  P_a_raw <-  matrix(0,1,d)
  alpha_EST <- matrix(0,1,d)
  alpha_SE  <- matrix(0,1,d)
  beta_EST <- matrix(0,1,d)
  beta_SE <-  matrix(0,1,d)
  P_raw_DLASSO <-  matrix(0,1,d)
  M1 <- t(t(M_raw[,1]))
  
  for (k in 1:d){
    M <- M_raw
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
    #CN <- round(n/log(n))
    
    ####################
    MXZ <- cbind(MT,X)
    fit.dlasso  <- DLASSO_fun(MXZ, Y)
    beta_est <- fit.dlasso[1]
    beta_se  <- fit.dlasso[2]
    P_b <-  2*(1-pnorm(abs(beta_est/beta_se),0,1))
    beta_EST[k] <-  beta_est
    beta_SE[k]  <- beta_se
    ##
    XZ <- X
    lm.fit <- lm(MT[,1]~XZ)
    lm.out <- summary(lm.fit)
    alpha_est <- lm.out$coefficients[2,1]
    alpha_se <- lm.out$coefficients[2,2]
    P_a <-  2*(1-pnorm(abs(alpha_est/alpha_se),0,1))
    #ab_est_DLASSO[times,k] <- alpha_est*beta_est
    P_raw_DLASSO[k] <- max(P_a,P_b)
    alpha_EST[k] <- alpha_est
    alpha_SE[k] <- alpha_se
    
  } #the end of k
  
  ##
  P_adj_DLASSO <- as.numeric(P_raw_DLASSO)
  ### The FDP method
  set <- which(P_adj_DLASSO < 0.05)
  hom <- hommel(P_adj_DLASSO, simes = FALSE)
  N1 <- discoveries(hom, set,incremental = TRUE, alpha=0.05)
  if (length(set) > 0){
    L <- length(set)
    N2 <- matrix(0,1,L)
    N2[2:L] <- N1[1:(L-1)]
  }
  
  N0 <- N1 - N2
  
  ID_FDR <- set[which(N0 > 0)]
  
  #print(ID_FDR)
  out_result <- list(ID = ID_FDR, alpha_hat = alpha_EST[ID_FDR], alpha_hat_SE= alpha_SE[ID_FDR], beta_hat=beta_EST[ID_FDR], beta_hat_SE =beta_SE[ID_FDR])
  return(out_result)
  
  }