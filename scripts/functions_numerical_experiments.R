


# rho (SAR) estimation functions ----

sar.lag.mixed.f <- function(rho, env) {
  beta_ridge = get("beta_ridge", envir=env)
  if (is.null(beta_ridge)){
    y.l <- get("y", envir = env)- rho*get("wy", envir = env)
    SSE <- crossprod(y.l)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- spatialreg::do_ldet(rho, env)
    ret <- (ldet - ((n/2)*log(2*(3.141593))) - (n/2)*log(s2)-(n/2)) # the term n/2 comes from (1/(2*s2))*crossprod(y.l)
  }else{
    e.lm.null <- get("y", envir=env) - get("x", envir=env)%*%beta_ridge
    e.lm.w <- get("wy", envir=env) - get("x", envir=env)%*%beta_ridge
   
    SSE <-crossprod(e.lm.null-rho*e.lm.w)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- spatialreg::do_ldet(rho, env)
    ret <- (ldet - ((n/2)*log(2*(3.141593))) - (n/2)*log(s2)-(n/2))
  }
  ret
}


for_ridge_sar<-function(data,formula, beta_ridge, scale= FALSE){
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
  x <- x[,-1]
  n <- NROW(x)
  m <- NCOL(x)
  xcolnames <- colnames(x)
  K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
  wy <- lag.listw(listw, y, zero.policy=TRUE)
  
  can.sim <- spatialreg::can.be.simmed(listw)
  env <- new.env(parent=globalenv())
  env <- new.env()
  assign("y", y, envir=env)
  assign("beta_ridge", beta_ridge, envir=env)
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

    
    e.lm.null <- y - x%*%beta_ridge
    
    e.lm.w <- wy- x%*%beta_ridge
    
    SSE <-crossprod(e.lm.null-rho*e.lm.w)
    
    s2 <- SSE/n
    
  
  
  ldet <- spatialreg::do_ldet(rho, env)
  
  #Output parameters
 
  y <- y - rho*wy
  x <- x
  ldet <- ldet
  sar_arguments <- list(rho=rho, y = y, x = x, W=mat_w, ldet=ldet, s2=s2)
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
  x <- model.matrix(mt, mf)[,-1]
  Neigh <- cell2nb(30, 30, type="rook", torus=FALSE, legacy=FALSE)
  
  
  # Step 1----
  # Estimate gamma and beta_ridge for Y=X*beta_ridge + epsilon
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  lambda_int<-fit$lambda
  n<-nrow(x)
  logLik_sar <- matrix(NA, n, length(lambda_int))
  
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
      
      
      
      x_test <-x[data$SLOO=='Test',]
      y_test <-y[data$SLOO=='Test']
      
      ridge_coef<-function(lambda){
      
        
        coefs<-  drop(solve(crossprod(x_train) + diag(lambda, (p)), crossprod(x_train, y_train)))
        
        
        y_pred = x_train%*%coefs
        
        SSE <- crossprod(y_pred-y_train)
        s2 <- SSE/nrow(x_train)
      
        
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        

        res <- logLik
        
        return(res)
      }
      
      
      
      logLik_sar[k,] <- sapply(lambda_int, ridge_coef)
      
      
  }
  
  #Compute MSE and select best lambda sar
  mlogLik_sar= colMeans(logLik_sar)
  plotdat_sar<-tibble(logLik=mlogLik_sar, lambda=lambda_int)
  lL_sar<- max(plotdat_sar$logLik)
  lambda_sar<- plotdat_sar$lambda[plotdat_sar$logLik==lL_sar]
  
  
  coef_sar <-  drop(solve(crossprod(x) + diag(lambda_sar, p), crossprod(x, y)))
  predictorsnames<-all.vars(model[[3]])
  names(coef_sar)<- c(predictorsnames )

  
  
  j=1
  p <- ncol(x)
  
  # set or initialize some parameters
  vector_rho <- numeric(100) # to store the rho values
  vector_gamma <- numeric(100) # to store the gamma values
  matrix_beta <- matrix(NA,100,p) # to store the beta coefficients
  tol=0.000001 # tolerance to assess convergence of the algorithm
  delta_rho=1 # initialize the element related to rho convergence in the condition
  delta_beta=1 # initialize the element related to beta convergence in the condition
  
  beta_ridge=coef_sar
  
  # Step 2----
  param_ridge_sar <- for_ridge_sar(data, model, beta_ridge, scale=scale) # rho estimation for y=rho*W*y+ X*beta_ridge + epsilon
  
  while ((delta_rho > tol)||(delta_beta > tol)) {
    
    
    
    #Step 3----
    x_filtered <-param_ridge_sar$x
    
    
    y_filtered <-param_ridge_sar$y
    
    rho <- param_ridge_sar$rho
    s2 <- param_ridge_sar$s2
    ldet <- param_ridge_sar$ldet
    fit <- glmnet(x_filtered, y_filtered, alpha = 0, standardize = FALSE, intercept=FALSE)
    lambda_int<-fit$lambda
    
    
    logLik_sar <- matrix(NA, n, length(lambda_int))
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
        
      
        x_train <- x_filtered[data$SLOO=='Train',]
        y_train <- y_filtered[data$SLOO=='Train']
        
        x_test <-x_filtered[data$SLOO=='Test',]
        y_test <-y_filtered[data$SLOO=='Test']
        
        ridge_coef<-function(lambda){
          
          
          
          coefs<-  drop(solve(crossprod(x_train) + diag(lambda, (p)), crossprod(x_train, y_train)))
          
        
          logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
          
          res <- logLik
          
          return(res)
        }
        
        
        
        logLik_sar[k,] <- sapply(lambda_int, ridge_coef)
        
    }
    #Compute MSE and select best lambda sar
    mlogLik_sar= colMeans(logLik_sar)
    
    plotdat_sar<-tibble(logLik=mlogLik_sar, lambda=lambda_int)
    lL_sar<- max(plotdat_sar$logLik)
    
    lambda_sar<- plotdat_sar$lambda[plotdat_sar$logLik==lL_sar]
    coef_sar <-  drop(solve(crossprod(x_filtered) + diag(lambda_sar, p), crossprod(x_filtered, y_filtered)))
    predictorsnames<-all.vars(model[[3]])
    names(coef_sar)<- c(predictorsnames )
    
    beta_ridge=coef_sar
    
    # Step 4----
    param_ridge_sar <- for_ridge_sar(data, model, beta_ridge,scale=scale) # rho estimation for y=rho*W*y+ X*beta_ridge + epsilon
    rho <-  param_ridge_sar$rho
    x_filtered <- param_ridge_sar$x
    W <- param_ridge_sar$W
    
    beta_ridge=coef_sar
    
    #Stop criterion
    if(j==1){
      rho_new=rho
      delta_rho=rho_new
      
      beta_new=coef_sar
      delta_beta=beta_new
      
    }else{
      rho_old=rho_new
      rho_new=rho
      delta_rho=rho_old-rho_new
      
      beta_old=beta_new
      beta_new=coef_sar
      delta_beta=sqrt(sum((beta_old-beta_new)^2))
    }
    
    vector_rho[j]<-rho_new
    matrix_beta[j,]<-beta_new
    vector_gamma[j]<-lambda_sar
    
    j=j+1
  }
  
  
  
  vector_rho <- vector_rho[1:j-1] 
  matrix_beta <- matrix_beta[1:j-1,]
  results<-list(rho = rho_new,
                Coefficients = coef_sar)
  
  return(results)
}


# lambda (SEM) estimation functions ----

sem_error_sse <- function(lambda, env) {
  beta_ridge = get("beta_ridge", envir=env)
  if (is.null(beta_ridge)){
    yl <- get("y", envir=env) - lambda * get("wy", envir=env)
    SSE <- crossprod(yl) } else { 
      yl <- get("y", envir=env) - lambda * get("wy", envir=env)
      n <-get("n", envir=env)
      xl <- get("x", envir=env)%*%beta_ridge - lambda * get("WX", envir=env)%*%beta_ridge
      SSE <-  crossprod(yl - xl)}
  SSE}


sem.error.f <- function(lambda, env) {
  SSE <- sem_error_sse(lambda, env)
  n <- get("n", envir=env)
  s2 <- SSE/n
  ldet <- spatialreg::do_ldet(lambda, env)
  ret <- (ldet - ((n/2)*log(2*(3.141593))) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
  if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret, " Jacobian:", ldet, " SSE:", SSE, "\n")
  assign("f_calls", get("f_calls", envir=env)+1L, envir=env)
  ret
}


for_ridge_sem<-function(data,formula, beta_ridge, scale = FALSE){
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
  x <- x[,-1]
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
  
  
  # Loglikelihood test for lambda significance 

    yl <- y - lambda *wy
    xl <- x%*%beta_ridge - lambda *WX%*%beta_ridge
    SSE <-  crossprod(yl - xl)
    s2 <- SSE/n
    ldet <- spatialreg::do_ldet(lambda, env)
    ret <- (ldet - ((n/2)*log(2*(3.141593))) -  (n/2)*log(s2))
    ret_0 <- ( - ((n/2)*log(2*(3.141593))) - (n/2)*log(crossprod(y-x%*%beta_ridge)/n))

  
  
  #Output parameters
  LR <- 2*(ret-ret_0)
  p_value <- pchisq(LR, df=1, lower.tail = FALSE)
  y <- y - lambda*wy
  x <- x - lambda*WX
  ldet <- ldet
  
  sem_arguments <- list(lambda=lambda, y = y, x = x, W = mat_w, pvalue=p_value, statistic=LR, ldet = ldet, s2=s2)
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
  x <- x[,-1]
  
  Neigh <- cell2nb(30, 30, type="rook", torus=FALSE, legacy=FALSE)
  n <- nrow(x)
  
  # Step 1----
  # Estimate gamma and beta_ridge for Y=X*beta_ridge + epsilon
   fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  lambda_int<-fit$lambda
  logLik_sem <- matrix(NA, n, length(lambda_int))
  
  
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
      
      ldet=0
      
      ridge_coef<-function(lambda){
        
        
        x<-x_train
        
        coefs<-  drop(solve(crossprod(x_train) + diag(lambda, (p)), crossprod(x_train, y_train)))
        
        
        y_pred = x%*%coefs
        
        SSE <- crossprod(y_pred-y_train)
        s2 <- SSE/nrow(x_train)
        
        
        x<-x_test
        
        logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x%*%coefs)
        
        res <- logLik
        
        return(res)
      }
      
      
      logLik_sem[k,] <- sapply(lambda_int, ridge_coef)
      
      
  }
  #Compute MSE and select best lambda sem
  mlogLik_sem= colMeans(logLik_sem)
  
  plotdat_sem<-tibble(logLik=mlogLik_sem, lambda=lambda_int)
  lL_sem<- max(plotdat_sem$logLik)
  lambda_sem<- plotdat_sem$lambda[plotdat_sem$logLik==lL_sem]
  p<- ncol(x)
  coef_sem <- drop(solve(crossprod(x) + diag(lambda_sem, (p)), crossprod(x, y)))
  predictorsnames<-all.vars(model[[3]])
  names(coef_sem)<- c(predictorsnames )
  j=1
  vector_lambda <- numeric(200) # to store the lambda values
  matrix_beta <- matrix(NA,200,p) # to store the beta coefficients
  tol=0.000001 # tolerance to assess convergence of the algorithm
  delta_lambda=1 # initialize the element related to lambda convergence in the condition
  delta_beta=1 # initialize the element related to beta convergence in the condition
  
  beta_ridge = coef_sem
  
  # Step 2----
  param_ridge_sem <- for_ridge_sem(data, model, beta_ridge, scale=scale) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
  
  
  while ((delta_lambda > tol)||(delta_beta > tol)) {
    
    
    
    #Step 3----
    x_filtered<-param_ridge_sem$x
    
    
    y_filtered<-param_ridge_sem$y
    
    lambda<-param_ridge_sem$lambda
    ldet<-param_ridge_sem$ldet
    s2<-param_ridge_sem$s2
    fit <- glmnet(x_filtered, y_filtered, alpha = 0, standardize = FALSE, intercept=FALSE)
    lambda_int<-fit$lambda
    n<-nrow(x)
    logLik_sem <- matrix(NA, n, length(lambda_int))
    
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
        
        
        ridge_coef<-function(lambda){
          
          
          
          coefs<-  drop(solve(crossprod(x_train) + diag(lambda, (p)), crossprod(x_train, y_train)))
          
          
          logLik <- -(1/2)*log(2*(3.141593)) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
          
       
          
          res <- logLik
          
          return(res)
        }
        
        
        
        logLik_sem[i,] <- sapply(lambda_int, ridge_coef)
        
        
    }
    
    mlogLik_sem= colMeans(logLik_sem)
    
    plotdat_sem<-tibble(logLik=mlogLik_sem, lambda=lambda_int)
    lL_sem<- max(plotdat_sem$logLik)
    lambda_sem<- plotdat_sem$lambda[plotdat_sem$logLik==lL_sem]
    coef_sem <- drop(solve(crossprod(x_filtered) + diag(lambda_sem, (p)), crossprod(x_filtered, y_filtered)))
    predictorsnames<-all.vars(model[[3]])
    names(coef_sem)<- c(predictorsnames )
    W <- param_ridge_sem$W
    
    beta_ridge = coef_sem
    
    
    # Step 4----
    param_ridge_sem <- for_ridge_sem(data, model, beta_ridge, scale=scale) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
    
    x_filtered <-param_ridge_sem$x
    lambda <- param_ridge_sem$lambda
    W <- param_ridge_sem$W
    
    y_predicted <- solve((diag(1,n)-lambda*W),(x_filtered%*%coef_sem))
    
    #Stop criterion
    if(j==1){
      lambda_new=lambda
      delta_lambda=lambda_new
      
      beta_new=coef_sem
      delta_beta=beta_new
      
    }else{
      lambda_old=lambda_new
      lambda_new=lambda
      delta_lambda=lambda_old-lambda_new
      
      beta_old=beta_new
      beta_new=coef_sem
      delta_beta=sqrt(sum((beta_old-beta_new)^2))
    }
    
    vector_lambda[j]<-lambda_new
    matrix_beta[j,]<-beta_new
    j=j+1
  }
  
  
  
  
  
  vector_lambda <- vector_lambda[1:j-1] 
  matrix_beta <- matrix_beta[1:j-1,]
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
  parameter <- seq(0.01, 1, 0.01)
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


