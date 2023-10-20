# rho (SAR) estimation functions ----

sar.lag.mixed.f <- function(rho, env) {
  # Get the beta coefficients from the environment
  beta_0= get("beta_0", envir=env) # Ridge coefficient estimated using beta_0 = (t(X)%*%X+gamma_0*I_p)%*%t(X)Y
  beta_l= get("beta_l", envir=env) # Ridge coefficient estimated using beta_l = (t(X)%*%X+gamma_l*I_p)%*%t(X)WY
  
  
  # Calculate the residuals for  y=beta_0*X
  e.lm.null <- get("y", envir=env) - get("x", envir=env)%*%beta_0
  
  # Calculate the residuals for  Wy=beta_l*X
  e.lm.w <- get("wy", envir=env) - get("x", envir=env)%*%beta_l
  
  # Calculate the sum of squared errors (SSE)
  SSE <-crossprod(e.lm.null-rho*e.lm.w)
  
  # Get the sample size (number of observations)
  n <- get("n", envir=env)
  
  # Calculate the sample variance (s2)
  s2 <- SSE/n
  
  # Calculate the log determinant of (I-rho*W)
  ldet <- spatialreg::do_ldet(rho, env)
  
  # Calculate the log-likelihood value (ret) based on the likelihood function
  ret <- (ldet - ((n/2)*log(2*(3.141593))) - (n/2)*log(s2))
  
  # Return the log-likelihood value
  ret
}

# Function to estimate spatial dependence parameter rho using maximum likelihood estimation
sar_estimation<-function(data,formula, beta_0, beta_l){
  verbose = FALSE
  
  
  method<-"eigen"  
  
  
  
  # Create spatial neighborhood weights
  Neigh <- cell2nb(30, 30, type="rook", torus=FALSE, legacy=FALSE)
  lags <- nblag(Neigh, 2)
  listw  <- nb2listw(lags[[1]],style="W", zero.policy = TRUE)
  mat_w <- listw2mat(listw)
  
  
  interval=NULL
  quiet=FALSE
  control=list()
  con <- list(tol.opt=.Machine$double.eps^0.5, returnHcov=TRUE,
              pWOrder=250, fdHess=NULL, optimHess=FALSE,
              optimHessMethod="optimHess", LAPACK=FALSE,
              compiled_sse=FALSE, Imult=2, cheb_q=5, MC_p=16, MC_m=30,
              super=NULL, spamPivot="MMD", in_coef=0.1, type="MC",
              correct=TRUE, trunc=TRUE, SE_method="LU", nrho=200,
              interpn=2000, small_asy=TRUE, small=1500, SElndet=NULL,
              LU_order=FALSE, pre_eig=NULL)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  
  
  # Extract dependent variable and covariates
  mt <- terms(formula, data = data)
  mf <- lm(formula, data, 
           method="model.frame")
  
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  n <- NROW(x)
  m <- NCOL(x)
  xcolnames <- colnames(x)
  K = 1
  wy <- lag.listw(listw, y, zero.policy=TRUE)
  
  
  
  # Assign variables to the environment
  env <- new.env(parent=globalenv())
  env <- new.env()
  assign("y", y, envir=env)
  assign("beta_0", beta_0, envir=env)
  assign("beta_l", beta_l, envir=env)
  assign("wy", wy, envir=env)
  assign("x", x, envir=env)
  assign("n", n, envir=env)
  assign("m", m, envir=env)
  assign("K", K, envir=env)
  assign("verbose", !quiet, envir=env)
  assign("family", "SAR", envir=env)
  assign("listw", listw, envir=env)
  assign("similar", FALSE, envir=env)
  assign("f_calls", 0L, envir=env)
  assign("verbose", FALSE, envir=env)
  interval <- spatialreg::jacobianSetup(method, env, con, pre_eig=con$pre_eig,
                                        trs=trs, interval=interval)
  assign("interval", interval, envir=env)
  
  # Obtain the SAR lag parameter by maximasing the log likelihood defined in  sar.lag.mixed.f
  opt <- optimize(sar.lag.mixed.f, interval=interval, 
                  maximum=TRUE, tol=con$tol.opt, env=env)
  
  rho <- opt$maximum
  
  
  # Compute sample variance
  e.lm.null <- y - x%*%beta_0
  
  e.lm.w <- wy- x%*%beta_l
  
  SSE <-crossprod(e.lm.null-rho*e.lm.w)
  
  s2 <- SSE/n
  
  ldet <- spatialreg::do_ldet(rho, env)
  
  #Output parameters
  y <- y - rho*wy
  x <- x
  ldet <- ldet
  sar_arguments <- list(rho=rho,
                        y = y, 
                        x = x, 
                        W=mat_w, 
                        ldet=ldet,
                        s2=s2 
  )
  return(sar_arguments)
}




# Function for Ridge Regression in Spatial Lag Model (RRSAR)----
#

# Purpose:
# This function performs ridge regression for spatial lag model 
# using an alternating minimization algorithm.

# Inputs:
# - data: Dataframe containing the data
# - model: Linear regression formula
# - buffer: Buffer size for Spatial Leave One Out (SLOO)
rrsar <- function(data, model, buffer){
  
  # Store number of observations in n
  n<-nrow(data)
  
  # Extract dependent variable and covariates
  mt <- terms(model, data = data)
  mf <- lm(model, data, 
           method="model.frame")
  
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  
  # Create spatial neighborhood weights
  Neigh <- cell2nb(30, 30, type="rook", torus=FALSE, legacy=FALSE)
  lags <- nblag(Neigh, 2)
  listw  <- nb2listw(lags[[1]],style="W", zero.policy = TRUE)
  wy <- lag.listw(listw, y, zero.policy=TRUE)
  
  # Step 1: Ridge Estimation----
  # 
  # Estimate gamma_0, gamma_l, beta_0 and beta_l with  beta_0 = (t(X)%*%X+gamma_0*I_p)%*%t(X)Y and  beta_l = (t(X)%*%X+gamma_l*I_p)%*%t(X)WY
  
  # Obtain interval to search for gamma_0, gamma_l
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda
  
  # Create matrices to store logLik values for SLOO
  logLik_beta_0 <- matrix(NA, n, length(gamma_int))
  logLik_beta_l <- matrix(NA, n, length(gamma_int))
  
  
  #Perform SLOO to select best gamma_0 and gamma_l
  # Loop through data points
  for (k in 1:n){
    
    
    # Assign 'Train', 'Test' and 'DeadZone' labels to partition the dataset 
    # to achieve independence between training and testing sets
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[k]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[k]]))]<-'DeadZone'
      data$SLOO[k]<-'Test'
      
      
      # Store number of covariates in p
      p <- ncol(x)
      
      
      # Separate data into training and testing sets based on 'SLOO' labels
      x_train <- x[data$SLOO=='Train',]
      y_train <- y[data$SLOO=='Train']
      wy_train  <- wy[data$SLOO=='Train']
      
      
      x_test <-x[data$SLOO=='Test',]
      y_test <-y[data$SLOO=='Test']
      wy_test <- wy[data$SLOO=='Test']
      
      # Log-likelihood estimation function for Y= beta_0*X+epsilon
      logLik_conditional_0<-function(gamma){
        
        
        # estimate regression coefficients in the training set
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        
        y_pred = x_train%*%coefs
        
        SSE <- crossprod(y_pred-y_train)
        s2 <- SSE/nrow(x_train)
        
        # compute log-likelihood on the testing set
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        
        res <- logLik
        
        return(res)
      }
      
      
      # Compute log-likelihood for all possible gamma_0 values
      logLik_beta_0[k,] <- sapply(gamma_int, logLik_conditional_0)
      
      
      # Log-likelihood estimation function for WY= beta_l*X+epsilon
      logLik_conditional_l<-function(gamma){
        
        
        
        # estimate coefficients on the training set
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, wy_train)))
        
        
        y_pred = x_train%*%coefs
        
        SSE <- crossprod(y_pred-wy_train)
        s2 <- SSE/nrow(x_train)
        
        # compute log-likelihood on the testing set
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(wy_test-x_test%*%coefs)
        
        
        res <- logLik
        
        return(res)
      }
      
      # Compute log-likelihood for all gamma values
      logLik_beta_l[k,] <- sapply(gamma_int, logLik_conditional_l)
      
  }
  
  # Compute mean log-likelihood and select gamma_0 and gamma_l with highest log-likelihood
  mlogLik_beta_0= colMeans(logLik_beta_0)
  mlogLik_beta_l= colMeans(logLik_beta_l)
  gamma_beta_0<- gamma_int[which.max(mlogLik_beta_0)]
  gamma_beta_l<- gamma_int[which.max(mlogLik_beta_l)]
  
  
  
  # Estimate coefficients for gamma_0 and gamma_l selected through SLOO
  coef_0_ridge <-  drop(solve(crossprod(x) + diag(gamma_beta_0, p), crossprod(x, y)))
  coef_l_ridge <-  drop(solve(crossprod(x) + diag(gamma_beta_l, p), crossprod(x, wy)))
  
  
  # Assign names to coefficients according to covariates' names
  predictorsnames<-all.vars(model[[3]])
  names(coef_0_ridge)<- c(predictorsnames )
  names(coef_l_ridge)<- c(predictorsnames )
  
  # Step 2: Estimate rho ----
  # considering 
  # e.lm.null <- y - x%*%beta_0
  # e.lm.w <- wy - x%*%beta_l
  ridge_sar <- sar_estimation(data, model,beta_0=coef_0_ridge, beta_l=coef_l_ridge)
  
  
  
  # Step 3: Update rho, filtered dependent variable and perform ridge regression----
  y_filtered <- ridge_sar$y
  rho <- ridge_sar$rho
  s2 <- ridge_sar$s2
  
  # Obtain new interval to search for gamma
  fit <- glmnet(x, y_filtered, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda
  
  # Create matrix to store loglikelihood values
  logLik_sar <- matrix(NA, n, length(gamma_int))
  
  #Perform SLOO to select best gamma
  # Loop through each data point 
  for (k in 1:n){
    
    
    # Assign 'Train', 'Test' and 'DeadZone' labels to partition the dataset 
    # to achieve independence between training and testing sets
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[k]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[k]]))]<-'DeadZone'
      data$SLOO[k]<-'Test'
      
      # Separate data into training and testing sets based on 'SLOO' labels
      x_train <- x[data$SLOO=='Train',]
      y_train <- y_filtered[data$SLOO=='Train']
      
      x_test <-x[data$SLOO=='Test',]
      y_test <-y_filtered[data$SLOO=='Test']
      
      
      # Define function to compute log-likelihood
      logLik_conditional<-function(gamma){
        
        
        
        # estimate regression coefficients in the training set
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        
        # compute log-likelihood in the testing set
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        res <- logLik
        
        return(res)
      }
      
      
      
      
      # Compute log-likelihood for different gamma values and store in matrix
      logLik_sar[k,] <- sapply(gamma_int, logLik_conditional)
      
  }
  
  # Compute mean log-likelihood and select gamma with highest log-likelihood
  mlogLik_sar= colMeans(logLik_sar)
  gamma_sar<- gamma_int[which.max(mlogLik_sar)]
  
  # Update coefficients using selected gamma
  coef_sar <-  drop(solve(crossprod(x) + diag(gamma_sar, p), crossprod(x, y_filtered)))
  
  
  
  results <- list(Coefficients=coef_sar,
                  rho=ridge_sar$rho)
  
  
  return(results)
}



# lambda (SEM) estimation functions ----

# Function to calculate the log-likelihood for lambda estimation
sem.error.f <- function(lambda, env) {
  # Get the beta coefficients from the environment
  beta_ridge <- get("beta_ridge", envir=env)
  
  # filter dependent variable
  yl <- get("y", envir=env) - lambda * get("wy", envir=env)
  
  # Get the sample size (number of observations)
  n <-get("n", envir=env)
  
  #filter covariates
  xl <- get("x", envir=env)%*%beta_ridge - lambda * get("WX", envir=env)%*%beta_ridge
  
  # compute the sum of squared errors (SSE)
  SSE <-  crossprod(yl - xl)
  
  # Calculate the sample variance (s2)
  s2 <- SSE/n
  
  # Calculate the log determinant of (I-lambda*W)
  ldet <- spatialreg::do_ldet(lambda, env)
  
  
  # Calculate the log-likelihood value (ret) based on the likelihood function
  ret <- (ldet - ((n/2)*log(2*(3.141593))) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
  
  ret
}

# Function to estimate spatial dependence parameter lambda using maximum likelihood estimation
sem_estimation<-function(data,formula, beta_ridge){
  etype<-"error"
  zero.policy=TRUE
  method<-"eigen" 
  
  # Create spatial neighborhood weights
  Neigh <- cell2nb(30, 30, type="rook", torus=FALSE, legacy=FALSE)
  lags <- nblag(Neigh, 2)
  listw  <- nb2listw(lags[[1]],style="W", zero.policy = TRUE)
  mat_w <- listw2mat(listw)
  
  
  
  interval=NULL
  quiet=FALSE
  control=list()
  con <- list(tol.opt=.Machine$double.eps^0.5, returnHcov=TRUE,
              pWOrder=250, fdHess=NULL, optimHess=FALSE,
              optimHessMethod="optimHess", LAPACK=FALSE,
              compiled_sse=FALSE, Imult=2, cheb_q=5, MC_p=16, MC_m=30,
              super=NULL, spamPivot="MMD", in_coef=0.1, type="MC",
              correct=TRUE, trunc=TRUE, SE_method="LU", nrho=200,
              interpn=2000, small_asy=TRUE, small=1500, SElndet=NULL,
              LU_order=FALSE, pre_eig=NULL)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  
  mt <- terms(formula, data = data)
  mf <- lm(formula, data,  method="model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  n <- nrow(x)
  
  m <- NCOL(x)
  xcolnames <- colnames(x)
  K <- 1
  wy <- lag.listw(listw, y, zero.policy=TRUE)
  wx1 <- as.double(rep(1, n))
  wx <- lag.listw(listw, wx1, zero.policy=zero.policy)
  if (m > 1 || (m == 1 && K == 1)) {
    WX <- matrix(nrow=n,ncol=(m-(K-1)))
    for (k in K:m) {
      wx <- lag.listw(listw, x[,k], zero.policy=zero.policy)
      if (any(is.na(wx)))
        stop("NAs in lagged independent variable")
      WX[,(k-(K-1))] <- wx
    }
  }
  
  colnames(WX) <- xcolnames
  rm(wx)
  
  # Assign variables to the environment
  env <- new.env()
  assign("y", y, envir=env)
  assign("x", x, envir=env)
  assign("beta_ridge", beta_ridge, envir=env)
  assign("wy", wy, envir=env)
  assign("WX", WX, envir=env)
  assign("n", n, envir=env)
  assign("p", m, envir=env)
  assign("verbose", !quiet, envir=env)
  assign("family", "SAR", envir=env)
  assign("compiled_sse", con$compiled_sse, envir=env)
  assign("first_time", TRUE, envir=env)
  assign("LAPACK", con$LAPACK, envir=env)
  assign("listw", listw, envir=env)
  assign("similar", FALSE, envir=env)
  assign("f_calls", 0L, envir=env)
  assign("hf_calls", 0L, envir=env)
  assign("verbose", FALSE, envir = env)
  interval <- spatialreg::jacobianSetup(method, env, con, pre_eig=con$pre_eig,
                                        trs=trs, interval=interval)
  assign("interval", interval, envir=env)
  
  
  # Obtain the SEM spatial dependence parameter by maximasing the log likelihood defined in  sem.error.f
  opt <- optimize(sem.error.f, interval=interval, 
                  maximum=TRUE, tol=con$tol.opt, env=env)
  
  lambda <- opt$maximum
  
  
   " Compute sample variance"
  yl <- y - lambda *wy
  xl <- x%*%beta_ridge - lambda *WX%*%beta_ridge
  SSE <-  crossprod(yl - xl)
  s2 <- SSE/n
  
  
  
  #Output parameters
  y <- y - lambda*wy
  x <- x - lambda*WX
  
  
  sem_arguments <- list(lambda=lambda, y = y, x = x, W = mat_w,  s2=s2)
  return(sem_arguments)
}



# Function for Ridge Regression in Spatial Error Model (RRSEM) ----
# 

# Purpose:
# This function performs ridge regression for spatial error model 
# using an alternating minimization algorithm.

# Inputs:
# - data: Dataframe containing the data
# - model: Linear regression formula
# - buffer: Buffer size for Spatial Leave One Out (SLOO)
rrsem <- function(data, model, buffer){
  
  # Extrat dependent variable and covariates from model
  mf <- lm(model, data,  method="model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  
  Neigh <- cell2nb(30, 30, type="rook", torus=FALSE, legacy=FALSE) # to define buffer for SLOO
  
  n <- nrow(x)
  
  # Step 1: Ridge Estimation----
  # 
  # Estimate gamma and beta_ridge for Y=X*beta_ridge + epsilon
  
  # Obtain interval to search for gamma
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda
  
  # Create matrix to store logLik values for SLOO
  logLik_sem <- matrix(NA, n, length(gamma_int))
  
  
  #Perform SLOO to select best gamma
  # Loop through data points
  for (k in 1:n){
    
    
    # Assign 'Train', 'Test' and 'DeadZone' labels to partition the dataset 
    # to achieve independence between training and testing sets
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[k]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[k]]))]<-'DeadZone'
      data$SLOO[k]<-'Test'
      
      p <- ncol(x)
      
      # Separate data into training and testing sets based on 'SLOO' labels
      x_train <- x[data$SLOO=='Train',]
      y_train <- y[data$SLOO=='Train']

      
      x_test <- x[data$SLOO=='Test',]
      y_test <- y[data$SLOO=='Test']
      
      # Define function to compute log-likelihood
      logLik_conditional<-function(gamma){
        
        
        # estimate regression coefficients in the training set
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        # Compute sample variance
        y_pred = x_train%*%coefs
        SSE <- crossprod(y_pred-y_train)
        s2 <- SSE/nrow(x_train)
        
        # compute log-likelihood on the testing set
        logLik <- -(1/2)*log(2*pi) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        
        res <- logLik
        
        return(res)
      }
      
      
      # Compute log-likelihood for different gamma values and store in matrix
      logLik_sem[k,] <- sapply(gamma_int, logLik_conditional)
      
      
  }

  # Compute mean log-likelihood and select gamma with highest log-likelihood
  mlogLik_sem= colMeans(logLik_sem)
  gamma_sem<-  gamma_int[which.max(mlogLik_sem)]
  

  
  # Initialize beta coefficients using the ridge estimates with the best gamma from the previous step
  beta_ridge <- drop(solve(crossprod(x) + diag(gamma_sem, (p)), crossprod(x, y)))
  predictorsnames<-all.vars(model[[3]])
  names(beta_ridge)<- c(predictorsnames )
  
  
  # Step 2----
  #  Estimate lambda  for y= X*beta_ridge + u, u = lambda*W*u + epsilon
  param_ridge_sem <- sem_estimation(data, model, beta_ridge) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
  
  
  #Step 3----
  # Get filtered  dependent variable and filtered covariates
  x_filtered<-param_ridge_sem$x
  y_filtered<-param_ridge_sem$y
  
  lambda<-param_ridge_sem$lambda
  s2<-param_ridge_sem$s2
  
  
  # Obtain new interval to search for gamma
  fit <- glmnet(x_filtered, y_filtered, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda

  # Initialize storage for log-likelihood
  logLik_sem <- matrix(NA, n, length(gamma_int))
  
  #Perform SLOO to select best gamma
  # Loop through each data point 
  for (i in 1:n){
    
    # Assign 'Train', 'Test' and 'DeadZone' labels to partition the dataset 
    # to achieve independence between training and testing sets
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[i]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[i]]))]<-'DeadZone'
      data$SLOO[i]<-'Test'
      
      
      # Separate data into training and testing sets based on 'SLOO' labels
      x_train <- x_filtered[data$SLOO=='Train',]
      y_train <- y_filtered[data$SLOO=='Train']
    
      x_test <- x_filtered[data$SLOO=='Test',]
      y_test <- y_filtered[data$SLOO=='Test']
      
      # Log-likelihood estimation function for (Y-lambdaWY)= beta_ridge*(X-lambdaWX)+epsilon
      logLik_conditional<-function(gamma){
        
        
        # estimate regression coefficients in the training set
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        
        # compute log-likelihood on the testing set
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        
        
        res <- logLik
        
        return(res)
      }
      
      
      # Compute log-likelihood for different gamma values and store in matrix
      logLik_sem[i,] <- sapply(gamma_int, logLik_conditional)
      
      
  }
  
  # Compute mean log-likelihood and select gamma with highest log-likelihood
  mlogLik_sem= colMeans(logLik_sem)
  gamma_sem<-  gamma_int[which.max(mlogLik_sem)]
  
  # Update beta coefficients using the ridge estimate withe best gamma from the previous step
  beta_ridge <- drop(solve(crossprod(x_filtered) + diag(gamma_sem, (p)), crossprod(x_filtered, y_filtered)))
  predictorsnames<-all.vars(model[[3]])
  names( beta_ridge)<- c(predictorsnames )


  
  
  # Step 4: Estimate lambda using the new ridge coefficients----
  param_ridge_sem <- sem_estimation(data, model, beta_ridge) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
  
  
  # Get the new filtered variables
  x_filtered <-param_ridge_sem$x
  y_filtered <-param_ridge_sem$y
  lambda <- param_ridge_sem$lambda
  W <- param_ridge_sem$W
  
  
  
  
  #Setp 5: Perform ridge for the new filtered variables----
  
  # Obtain new interval to search for gamma
  fit <- glmnet(x_filtered, y_filtered, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda

  # Initialize storage for log-likelihood
  logLik_sem <- matrix(NA, n, length(gamma_int))
  
  #Perform SLOO to select best gamma
  # Loop through each data point 
  for (i in 1:n){
    
    
    # Assign 'Train', 'Test' and 'DeadZone' labels to partition the dataset 
    # to achieve independence between training and testing sets
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[i]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[i]]))]<-'DeadZone'
      data$SLOO[i]<-'Test'
      
      
      # Separate data into training and testing sets based on 'SLOO' labels
      x_train <- x_filtered[data$SLOO=='Train',]
      y_train <- y_filtered[data$SLOO=='Train']
      
      x_test <- x_filtered[data$SLOO=='Test',]
      y_test <- y_filtered[data$SLOO=='Test']
      
      # Log-likelihood estimation function for (Y-lambdaWY)= beta_ridge*(X-lambdaWX)+epsilon
      logLik_conditional<-function(gamma){
        
        
        # estimate regression coefficients in the training set
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        # compute log-likelihood on the testing set
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        
        
        res <- logLik
        
        return(res)
      }
      
      
      # Compute log-likelihood for different gamma values and store in matrix
      logLik_sem[i,] <- sapply(gamma_int, logLik_conditional)
      
      
  }
  
  # Compute mean log-likelihood and select gamma with highest log-likelihood
  mlogLik_sem= colMeans(logLik_sem)
  gamma_sem<-  gamma_int[which.max(mlogLik_sem)]
  
  
  # Update beta coefficients using the ridge estimate withe best gamma from the previous step
  coef_sem <- drop(solve(crossprod(x_filtered) + diag(gamma_sem, (p)), crossprod(x_filtered, y_filtered)))
  
  
  #Output parameters
  results<-list(
    lambda = param_ridge_sem$lambda,
    Coefficients = coef_sem
  )
  
  
  return(results)
} 

# Function to perform Ridge Regression
ridgeRegression <- function(y, x) {
  
  p <- ncol(x)
  
  
  # Calculate bR for gamma in parameter
  # bR is the standardized regression coefficients
  # c is gamma constant
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  parameter <- fit$lambda
  bR <- matrix(0, nrow = length(parameter), ncol = p)
  counter <- 1
  for (c in parameter) {
    bR[counter, ] <-  drop(solve(crossprod(x) + diag(c, (p)), crossprod(x, y)))
    counter <- counter + 1
  }
  
  # Calculate difference between successive elements in bR
  g <- matrix(0, nrow = p, ncol = length(parameter) - 1)
  tb <- t(bR)
  for (i in 1:p) {
    g[i, ] <- diff(tb[i, ])
  }
  
  # Find the index of first-order difference in g for each regression coefficient
  pi <- c()
  for (i in 1:p) {
    for (j in 1:(length(parameter) - 1)) {
      if (abs(g[i, j]) < 0.00105) {
        pi <- c(pi, j)
        break
      }
    }
  }
  gamma <- parameter[max(pi)]
  
  # Calculate bR using the selected gamma constant
  bRk <- matrix(0, nrow = p, ncol = 1)
  bRk <- drop(solve(crossprod(x) + diag(gamma, (p)), crossprod(x, y)))
  
  
  
  
  
  return(list(beta = bRk, gamma = gamma))
}




