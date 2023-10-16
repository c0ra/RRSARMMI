# rho (SAR) estimation functions ----

sar.lag.mixed.f <- function(rho, env) {
  beta_0= get("beta_0", envir=env)
  beta_l= get("beta_l", envir=env)
  
  e.lm.null <- get("y", envir=env) - get("x", envir=env)%*%beta_0
  e.lm.w <- get("wy", envir=env) - get("x", envir=env)%*%beta_l
  
  SSE <-crossprod(e.lm.null-rho*e.lm.w)
  n <- get("n", envir=env)
  s2 <- SSE/n
  ldet <- spatialreg::do_ldet(rho, env)
  ret <- (ldet - ((n/2)*log(2*(3.141593))) - (n/2)*log(s2)-(n/2))
  
  ret
}


sar_estimation<-function(data,formula, beta_0, beta_l, scale= FALSE){
  verbose = FALSE
  
  
  method<-"eigen"  
  
  
  
  
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
  mf <- lm(formula, data, 
           method="model.frame")
  
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  n <- NROW(x)
  m <- NCOL(x)
  xcolnames <- colnames(x)
  K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
  wy <- lag.listw(listw, y, zero.policy=TRUE)
  
  can.sim <- spatialreg::can.be.simmed(listw)
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
  assign("can.sim", can.sim, envir=env)
  assign("listw", listw, envir=env)
  assign("similar", FALSE, envir=env)
  assign("f_calls", 0L, envir=env)
  assign("verbose", FALSE, envir=env)
  interval <- spatialreg::jacobianSetup(method, env, con, pre_eig=con$pre_eig,
                                        trs=trs, interval=interval)
  assign("interval", interval, envir=env)
  
  
  opt <- optimize(sar.lag.mixed.f, interval=interval, 
                  maximum=TRUE, tol=con$tol.opt, env=env)
  
  rho <- opt$maximum
  
  
  # To compute rho p-value using LR
  
  
  e.lm.null <- y - x%*%beta_0
  
  e.lm.w <- wy- x%*%beta_l
  
  SSE <-crossprod(e.lm.null-rho*e.lm.w)
  
  s2 <- SSE/n
  
  # beta <- beta_0-rho*beta_l
  
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
                        s2=s2 #,
                        # beta=beta
  )
  return(sar_arguments)
}




# 2.2 Alternating minimization algorithm  (SAR) ----
#Maximizing likelihood
rrsar <- function(data, model, buffer, scale=FALSE){
  
  
  n<-nrow(data)
  mt <- terms(model, data = data)
  mf <- lm(model, data, 
           method="model.frame")
  
  y <- model.extract(mf, "response")
  
  x <- model.matrix(mt, mf)
  Neigh <- cell2nb(30, 30, type="rook", torus=FALSE, legacy=FALSE)
  lags <- nblag(Neigh, 2)
  listw  <- nb2listw(lags[[1]],style="W", zero.policy = TRUE)
  wy <- lag.listw(listw, y, zero.policy=TRUE)
  
  # Step 1----
  # Estimate gamma and beta_ridge for Y=X*beta_ridge + epsilon
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda
  n<-nrow(x)
  logLik_beta_0 <- matrix(NA, n, length(gamma_int))
  logLik_beta_l <- matrix(NA, n, length(gamma_int))
  # Define training and testing data
  for (k in 1:n){
    
    
    #Define buffer
    
    
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[k]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[k]]))]<-'DeadZone'
      data$SLOO[k]<-'Test'
      
      x_train <- x[data$SLOO=='Train',]
      y_train <- y[data$SLOO=='Train']
      wy_train  <- wy[data$SLOO=='Train']
      
      p <- ncol(x)
      
      
      
      x_test <-x[data$SLOO=='Test',]
      y_test <-y[data$SLOO=='Test']
      wy_test <- wy[data$SLOO=='Test']
      
      
      logLik_conditional_0<-function(gamma){
        
        
        
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        
        y_pred = x_train%*%coefs
        
        SSE <- crossprod(y_pred-y_train)
        s2 <- SSE/nrow(x_train)
        
        
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        
        res <- logLik
        
        return(res)
      }
      
      
      
      logLik_beta_0[k,] <- sapply(gamma_int, logLik_conditional_0)
      
      
      
      logLik_conditional_l<-function(gamma){
        
        
        
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, wy_train)))
        
        
        y_pred = x_train%*%coefs
        
        SSE <- crossprod(y_pred-wy_train)
        s2 <- SSE/nrow(x_train)
        
        
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(wy_test-x_test%*%coefs)
        
        
        res <- logLik
        
        return(res)
      }
      
      
      logLik_beta_l[k,] <- sapply(gamma_int, logLik_conditional_l)
      
  }
  
  #Compute MSE and select best lambda sar
  mlogLik_beta_0= colMeans(logLik_beta_0)
  mlogLik_beta_l= colMeans(logLik_beta_l)
  # plotdat_sar<-tibble(logLik=mlogLik_sar, lambda=gamma_int)
  gamma_beta_0<- gamma_int[which.max(mlogLik_beta_0)]
  gamma_beta_l<- gamma_int[which.max(mlogLik_beta_l)]
  
  
  
  # coef_sar <-  drop(solve(crossprod(x) + diag(gamma, p), crossprod(x, y)))
  
  coef_0_ridge <-  drop(solve(crossprod(x) + diag(gamma_beta_0, p), crossprod(x, y)))
  coef_l_ridge <-  drop(solve(crossprod(x) + diag(gamma_beta_l, p), crossprod(x, wy)))
  
  
  
  predictorsnames<-all.vars(model[[3]])
  names(coef_0_ridge)<- c(predictorsnames )
  names(coef_l_ridge)<- c(predictorsnames )
  
  # Step 2----
  ridge_sar <- sar_estimation(data, model,beta_0=coef_0_ridge, beta_l=coef_l_ridge, scale=scale)
  
  
  
  #Step 3----
  y_filtered <- ridge_sar$y
  rho <- ridge_sar$rho
  s2 <- ridge_sar$s2
  fit <- glmnet(x, y_filtered, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda
  
  
  logLik_sar <- matrix(NA, n, length(gamma_int))
  # Define training and testing data
  for (k in 1:n){
    
    
    #Define buffer
    
    
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[k]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[k]]))]<-'DeadZone'
      data$SLOO[k]<-'Test'
      
      
      x_train <- x[data$SLOO=='Train',]
      y_train <- y_filtered[data$SLOO=='Train']
      
      x_test <-x[data$SLOO=='Test',]
      y_test <-y_filtered[data$SLOO=='Test']
      
      logLik_conditional<-function(gamma){
        
        
        
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        res <- logLik
        
        return(res)
      }
      
      
      
      logLik_sar[k,] <- sapply(gamma_int, logLik_conditional)
      
  }
  #Compute MSE and select best lambda sar
  mlogLik_sar= colMeans(logLik_sar)
  
  gamma_sar<- gamma_int[which.max(mlogLik_sar)]
  coef_sar <-  drop(solve(crossprod(x) + diag(gamma_sar, p), crossprod(x, y_filtered)))
  
  
  
  results <- list(Coefficients=coef_sar,
                  rho=ridge_sar$rho)
  
  
  return(results)
}


# lambda (SEM) estimation functions ----


sem.error.f <- function(lambda, env) {
  beta_ridge <- get("beta_ridge", envir=env)
  yl <- get("y", envir=env) - lambda * get("wy", envir=env)
  n <-get("n", envir=env)
  xl <- get("x", envir=env)%*%beta_ridge - lambda * get("WX", envir=env)%*%beta_ridge
  SSE <-  crossprod(yl - xl)
  n <- get("n", envir=env)
  s2 <- SSE/n
  ldet <- spatialreg::do_ldet(lambda, env)
  ret <- (ldet - ((n/2)*log(2*(3.141593))) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
  if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret, " Jacobian:", ldet, " SSE:", SSE, "\n")
  assign("f_calls", get("f_calls", envir=env)+1L, envir=env)
  ret
}


sem_estimation<-function(data,formula, beta_ridge, scale = FALSE){
  etype<-"error"
  zero.policy=TRUE
  method<-"eigen" 
  
  
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
  K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
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
  if (K == 2) {
    # modified to meet other styles, email from Rein Halbersma
    wx1 <- as.double(rep(1, n))
    wx <- lag.listw(listw, wx1, zero.policy=zero.policy)
    if (m > 1) WX <- cbind(wx, WX)
    else WX <- matrix(wx, nrow=n, ncol=1)
  }
  colnames(WX) <- xcolnames
  rm(wx)
  
  can.sim <- FALSE
  if (listw$style %in% c("W", "S")) 
    can.sim <- spatialreg::can.be.simmed(listw)
  
  
  
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
  assign("can.sim", can.sim, envir=env)
  assign("listw", listw, envir=env)
  assign("similar", FALSE, envir=env)
  assign("f_calls", 0L, envir=env)
  assign("hf_calls", 0L, envir=env)
  assign("verbose", FALSE, envir = env)
  interval <- spatialreg::jacobianSetup(method, env, con, pre_eig=con$pre_eig,
                                        trs=trs, interval=interval)
  assign("interval", interval, envir=env)
  
  
  opt <- optimize(sem.error.f, interval=interval, 
                  maximum=TRUE, tol=con$tol.opt, env=env)
  
  lambda <- opt$maximum
  
  
  
  
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




# 2.2 Alternating minimization algorithm  (SEM) ----
#maximizing likelihood
rrsem <- function(data, model, buffer, scale=FALSE){
  
  
  mf <- lm(model, data,  method="model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  
  Neigh <- cell2nb(30, 30, type="rook", torus=FALSE, legacy=FALSE)
  n <- nrow(x)
  
  # Step 1----
  # Estimate gamma and beta_ridge for Y=X*beta_ridge + epsilon
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda
  logLik_sem <- matrix(NA, n, length(gamma_int))
  
  
  # Define training and testing data
  for (k in 1:n){
    
    
    #Define buffer
    
    
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[k]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[k]]))]<-'DeadZone'
      data$SLOO[k]<-'Test'
      
      
      
      x_train <- x[data$SLOO=='Train',]
      y_train <- y[data$SLOO=='Train']
      p <- ncol(x)
      
      x_test <- x[data$SLOO=='Test',]
      y_test <- y[data$SLOO=='Test']
      
      
      logLik_conditional<-function(gamma){
        
        
        x<-x_train
        
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        
        y_pred = x%*%coefs
        
        SSE <- crossprod(y_pred-y_train)
        s2 <- SSE/nrow(x_train)
        
        
        x<-x_test
        
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x%*%coefs)
        
        res <- logLik
        
        return(res)
      }
      
      
      logLik_sem[k,] <- sapply(gamma_int, logLik_conditional)
      
      
  }
  #Compute MSE and select best lambda sem
  mlogLik_sem= colMeans(logLik_sem)
  
  gamma_sem<-  gamma_int[which.max(mlogLik_sem)]
  
  p<- ncol(x)
  
  beta_ridge <- drop(solve(crossprod(x) + diag(gamma_sem, (p)), crossprod(x, y)))
  predictorsnames<-all.vars(model[[3]])
  names(beta_ridge)<- c(predictorsnames )
  
  
  # Step 2----
  param_ridge_sem <- sem_estimation(data, model, beta_ridge, scale=scale) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
  
  
  #Step 3----
  x_filtered<-param_ridge_sem$x
  
  
  y_filtered<-param_ridge_sem$y
  
  lambda<-param_ridge_sem$lambda
  s2<-param_ridge_sem$s2
  fit <- glmnet(x_filtered, y_filtered, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda
  n<-nrow(x)
  logLik_sem <- matrix(NA, n, length(gamma_int))
  
  # Define training and testing data
  for (i in 1:n){
    
    
    #Define buffer
    
    
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[i]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[i]]))]<-'DeadZone'
      data$SLOO[i]<-'Test'
      
      
      
      x_train <- x_filtered[data$SLOO=='Train',]
      y_train <- y_filtered[data$SLOO=='Train']
      p<-ncol(x_filtered)
      x_test <- x_filtered[data$SLOO=='Test',]
      y_test <- y_filtered[data$SLOO=='Test']
      
      
      logLik_conditional<-function(gamma){
        
        
        
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        
        
        res <- logLik
        
        return(res)
      }
      
      
      
      logLik_sem[i,] <- sapply(gamma_int, logLik_conditional)
      
      
  }
  
  mlogLik_sem= colMeans(logLik_sem)
  gamma_sem<-  gamma_int[which.max(mlogLik_sem)]
  
  coef_sem <- drop(solve(crossprod(x_filtered) + diag(gamma_sem, (p)), crossprod(x_filtered, y_filtered)))
  
  predictorsnames<-all.vars(model[[3]])
  names(coef_sem)<- c(predictorsnames )
  W <- param_ridge_sem$W
  
  beta_ridge = coef_sem
  
  
  # Step 4----
  param_ridge_sem <- sem_estimation(data, model, beta_ridge, scale=scale) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
  
  x_filtered <-param_ridge_sem$x
  y_filtered <-param_ridge_sem$y
  lambda <- param_ridge_sem$lambda
  W <- param_ridge_sem$W
  
  
  
  #Setp 5 ----
  
  fit <- glmnet(x_filtered, y_filtered, alpha = 0, standardize = FALSE, intercept=FALSE)
  gamma_int<-fit$lambda
  n<-nrow(x)
  logLik_sem <- matrix(NA, n, length(gamma_int))
  
  # Define training and testing data
  for (i in 1:n){
    
    
    #Define buffer
    
    
    data$SLOO<-'Train'
    if(buffer>1){
      data.lags <- nblag(Neigh, buffer)
      for (b in 1:buffer){
        data$SLOO[(unlist(data.lags[[b]][[i]]))]<-'DeadZone'
      }}else data$SLOO[(unlist(Neigh[[i]]))]<-'DeadZone'
      data$SLOO[i]<-'Test'
      
      
      
      x_train <- x_filtered[data$SLOO=='Train',]
      y_train <- y_filtered[data$SLOO=='Train']
      p<-ncol(x_filtered)
      x_test <- x_filtered[data$SLOO=='Test',]
      y_test <- y_filtered[data$SLOO=='Test']
      
      
      logLik_conditional<-function(gamma){
        
        
        
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        
        
        res <- logLik
        
        return(res)
      }
      
      
      
      logLik_sem[i,] <- sapply(gamma_int, logLik_conditional)
      
      
  }
  
  mlogLik_sem= colMeans(logLik_sem)
  gamma_sem<-  gamma_int[which.max(mlogLik_sem)]
  
  coef_sem <- drop(solve(crossprod(x_filtered) + diag(gamma_sem, (p)), crossprod(x_filtered, y_filtered)))
  
  
  
  results<-list(
    lambda = param_ridge_sem$lambda,
    Coefficients = coef_sem
  )
  
  
  return(results)
} 


ridgeRegression <- function(y, x) {
  
  p <- ncol(x)
  
  
  # Calculate bR for c in the range from 0 to 1 with the step of 0.01
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




