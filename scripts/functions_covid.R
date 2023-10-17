adjust_ridge<-function(lambda){
  
  coefs <-  drop(solve(crossprod(x_train) + diag(lambda, p-1), crossprod(x_train, y_train)))
  
  
  y_pred = x_test%*%coefs
  
  return(y_pred)
}


beta =  function(l){
  coefs <- drop(solve(crossprod(x_ridge) + diag(l, (p-1)), crossprod(x_ridge, y_ridge)))
  return(coefs)
}

ridge_coef<-function(lambda){
  
  
  x<-x_train
  
  coefs<-  drop(solve(crossprod(x_train) + diag(lambda, (p-1)), crossprod(x_train, y_train)))
  
  
  y_pred = x%*%coefs
  
  SSE <- crossprod(y_pred-y_train)
  s2 <- SSE/nrow(x_train)
  
  
  x<-x_test
  
  logLik <- -(1/2)*log(2*pi) -(1/2)*log(s2)+ ldet -(1/(2*s2))*crossprod(y_test-x%*%coefs)
  
 
  
  res <- logLik
  
  return(res)
}

# rho (SAR) estimation functions ----

sar.lag.mixed.f <- function(rho, env) {
  beta_ridge = get("beta_ridge", envir=env)
  if (is.null(beta_ridge)){
    y.l <- get("y", envir = env)- rho*get("wy", envir = env)
    SSE <- crossprod(y.l)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- spatialreg::do_ldet(rho, env)
    ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (n/2))
  }else{
    e.lm.null <- get("y", envir=env) - get("x", envir=env)%*%beta_ridge
    e.lm.w <- get("wy", envir=env) - get("x", envir=env)%*%beta_ridge
    e.a <- crossprod(e.lm.null)
    e.b <- crossprod(e.lm.null,e.lm.w)
    e.c <- crossprod(e.lm.w)
    #SSE <- e.a - 2*rho*e.b + rho*rho*e.c
    SSE <-crossprod(e.lm.null-rho*e.lm.w)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- spatialreg::do_ldet(rho, env)
    ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2)- (n/2))
  }
  ret
}


for_ridge_sar<-function(data,formula, beta_ridge, scale= FALSE){
  verbose = FALSE
  
  
  method<-"eigen"   
  
  data_sd<-as(data, "Spatial")
  coords = coordinates(data_sd)

  
  
  Neigh <- poly2nb(data_sd, row.names = NULL, queen=TRUE)
  covid.lags <- nblag(Neigh, 2)
  dlist_1 <- nbdists(covid.lags[[1]], coordinates(coords), longlat=TRUE)
  inv_dlist_1 <- lapply(dlist_1, function(x) 1/(x))
  listw_1  <- nb2listw(covid.lags[[1]],style="B",glist =inv_dlist_1, zero.policy = TRUE)
  mat_w <- listw2mat(listw_1)
  mat_w <- mat_w/eigen(mat_w)$values[1]
  listw <-mat2listw(mat_w)
  
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
  y <- scale(y, scale=scale) #center response
  x <- model.matrix(mt, mf)
  x<-  scale(x, scale=scale) #center variables
  x <- x[,-1]
  n <- NROW(x)
  m <- NCOL(x)
  xcolnames <- colnames(x)
  K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
  wy <- lag.listw(listw, y, zero.policy=TRUE)
  
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
  

  if (is.null(beta_ridge)){
    y.l <- get("y", envir = env)- rho*get("wy", envir = env)
    SSE <- crossprod(y.l)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- spatialreg::do_ldet(rho, env)
    
    ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (n/2)) # loglikelihood 
    
    ret_0 <- ( - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (n/2)) # loglikelihood when rho = 0
    
  }else{
    
    e.lm.null <- y - x%*%beta_ridge
    
    e.lm.w <- wy- x%*%beta_ridge
    
    SSE <-crossprod(e.lm.null-rho*e.lm.w)
    
    s2 <- SSE/n
    
    ldet <- spatialreg::do_ldet(rho, env)
    
    ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2)- (n/2)) # loglikelihood 
    
    ret_0 <- -((n/2)*log(2*pi)) - (n/2)*log(crossprod(e.lm.null)/n)- (n/2) # loglikelihood when rho = 0
  }
  
  ldet <- spatialreg::do_ldet(rho, env)
  
  #Output parameters
  LR <- 2*(ret-ret_0)
  p_value <- pchisq(LR, df=1, lower.tail = FALSE)
  y <- y - rho*wy
  x <- x
  ldet <- ldet
  sar_arguments <- list(rho=rho, y = y, x = x, W=mat_w, pvalue= p_value, statistic= LR, ldet=ldet, s2=s2)
  return(sar_arguments)
}


# 2.2 Alternating minimization algorithm  (SAR) ----
#Maximizing likelihood
# Inizialization with beta
rrsar <- function(data, model, buffer, plot = FALSE, scale=FALSE){
  
  
  n<-nrow(data)
  data.sp<-as(data,'Spatial') #transform to sp object for some functions
  
  mt <- terms(model, data = data)
  mf <- lm(model, data, 
           method="model.frame")
  
  y <- model.extract(mf, "response")
  y <- scale(y, scale)  #center or scale the response; depends on the parameter scale
  x <- model.matrix(mt, mf)[,-1]
  x<-  scale(x, scale)  #center or scale the covariates; depends on the parameter scale
  coords<-coordinates(data.sp)
  Neigh <- poly2nb(data.sp, row.names = NULL, queen=TRUE)
  
  # Step 1----
  # Estimate gamma and beta_ridge for Y=X*beta_ridge + epsilon
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  lambda_int<-n*fit$lambda
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
      
      
      train.data <- data[data$SLOO=='Train',]
      test.data  <- data[data$SLOO=='Test',]
      
      st_geometry(test.data) <- NULL
      predictorvarnames<-all.vars(model[[3]])
      
      x_train <- x[data$SLOO=='Train',]
      y_train <- y[data$SLOO=='Train']
      p <- ncol(test.data %>% dplyr::select(predictorvarnames))
      
      fit <- glmnet(x_train, y_train, alpha = 0, lambda=lambda_int, standardize = FALSE)
      
      
      
      x_test <-x[data$SLOO=='Test',]
      y_test <-y[data$SLOO=='Test']
      
      ridge_coef<-function(lambda){
        
        
        x<-x_train
        
        coefs<-  drop(solve(crossprod(x_train) + diag(lambda, (p)), crossprod(x_train, y_train)))
        
        
        y_pred = x%*%coefs
        
        SSE <- crossprod(y_pred-y_train)
        s2 <- SSE/nrow(x_train)
        
        
        x<-x_test
        
        logLik <- -(1/2)*log(2*pi) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x%*%coefs)
        
        
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
  y_predicted <- predict(fit, s = lambda_sar, newx = x)
  
  
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
  
  x_ridge<-param_ridge_sar$x
  
  
  y_ridge<-param_ridge_sar$y
  
  rho <- param_ridge_sar$rho
  s2 <- param_ridge_sar$s2
  ldet <- param_ridge_sar$ldet
  fit <- glmnet(x_ridge, y_ridge, alpha = 0, standardize = FALSE, intercept=FALSE)
  lambda_int<-n*fit$lambda
  
  while ((delta_rho > tol)||(delta_beta > tol)) {
    
    
    
    #Step 3----
    x_ridge<-param_ridge_sar$x
    
    
    y_ridge<-param_ridge_sar$y
    
    rho <- param_ridge_sar$rho
    s2 <- param_ridge_sar$s2
    ldet <- param_ridge_sar$ldet
    fit <- glmnet(x_ridge, y_ridge, alpha = 0, standardize = FALSE, intercept=FALSE)
    lambda_int<-n*fit$lambda
    
    
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
        
        
        train.data <- data[data$SLOO=='Train',]
        test.data  <- data[data$SLOO=='Test',]
        
        x_train <- x_ridge[data$SLOO=='Train',]
        y_train <- y_ridge[data$SLOO=='Train']
        
        x_test <-x_ridge[data$SLOO=='Test',]
        y_test <-y_ridge[data$SLOO=='Test']
        
        ridge_coef<-function(lambda){
          
          
          x<-x_train
          
          coefs<-  drop(solve(crossprod(x_train) + diag(lambda, (p)), crossprod(x_train, y_train)))
          
          
          y_pred = x%*%coefs
          
          
          
          
          x<-x_test
          
          logLik <- -(1/2)*log(2*pi) -(1/2)*log(s2)+ ldet -(1/(2*s2))*crossprod(y_test-x%*%coefs)
          
          
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
    coef_sar <-  drop(solve(crossprod(x_ridge) + diag(lambda_sar, p), crossprod(x_ridge, y_ridge)))
    predictorsnames<-all.vars(model[[3]])
    names(coef_sar)<- c(predictorsnames )
    
    beta_ridge=coef_sar
    
    # Step 4----
    param_ridge_sar <- for_ridge_sar(data, model, beta_ridge,scale=scale) # rho estimation for y=rho*W*y+ X*beta_ridge + epsilon
    rho <-  param_ridge_sar$rho
    x_ridge <- param_ridge_sar$x
    W <- param_ridge_sar$W
    y_predicted <- solve((diag(1,n)-rho*W),(x_ridge%*%coef_sar))
    
    
    beta_ridge=coef_sar
    
    #Stop criterion
    if(j==1){
      rho_new=param_ridge_sar$rho
      delta_rho=rho_new
      
      beta_new=coef_sar
      delta_beta=beta_new
      
    }else{
      rho_old=rho_new
      rho_new=param_ridge_sar$rho
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
  
  if (plot == TRUE){
    
    beta =  function(l){
      coefs <- drop(solve(crossprod(x_ridge) + diag(l, (p)), crossprod(x_ridge, y_ridge)))
      return(coefs)
    }
    plotdat <-sapply(lambda_int, beta)
    plotdat<-melt(plotdat, id.vars=c("Variable"))
    
    
    colnames(plotdat)<-c("Variable","Model","Coefficients")
    plotdat$lambda<-rep(lambda_int, each=(ncol(x)))
    
    
    
    p_coef_path<-ggplot(plotdat,aes(x=log(lambda), y=Coefficients,color=Variable) )+
      geom_line()+
      theme_classic()+
      ylab(NULL)+
      geom_vline(aes(xintercept = log(lambda_sar)), color="black")+
      geom_text(data = plotdat %>% filter(lambda == last(lambda)), aes(label = Variable, 
                                                                       x = log(lambda), 
                                                                       y = Coefficients, 
                                                                       color = Variable)) + 
      guides(color = "none") +
      xlab( expression(paste("log(", gamma,")")))
    
    
    p_logLik<-ggplot(plotdat_sar, aes(y=logLik, x=log(lambda)))+
      geom_line()+
      theme_classic()+
      geom_vline(aes(xintercept = log(lambda_sar), color="red"))+
      guides(color = "none") +
      xlab( expression(paste("log(", gamma,")")))
    
    p_coef_path
    
    p_logLik
    
  }else{plotdat=plotdat_sar}
  
  vector_rho <- vector_rho[1:j-1] 
  matrix_beta <- matrix_beta[1:j-1,]
  results<-list(Model = model,
                rho = param_ridge_sar$rho,
                Coefficients = coef_sar,
                gamma = lambda_sar,
                predicted_values = y_predicted,
                convergence_rho = vector_rho,
                convergence_gamma = vector_gamma,
                convergence_beta = matrix_beta,
                pvalue=param_ridge_sar$pvalue,
                statistic=param_ridge_sar$statistic,
                iterations=j-1,
                plotdat=plotdat,
                plotdat_sar=plotdat_sar)
  
  return(results)
}


# lambda (SEM) estimation functions ----

sem_error_sse <- function(lambda, env) {
  beta_ridge = get("beta_ridge", envir=env)
  if (is.null(beta_ridge)){
    yl <- get("y", envir=env) - lambda * get("wy", envir=env)
    yl <- get("sw", envir=env) * yl
    SSE <- crossprod(yl) } else { 
      yl <- get("y", envir=env) - lambda * get("wy", envir=env)
      yl <- get("sw", envir=env) * yl
      n <-get("n", envir=env)
      xl <- get("x", envir=env)%*%beta_ridge - lambda * get("WX", envir=env)%*%beta_ridge
      xl <- get("sw", envir=env) * xl
      SSE <-  crossprod(yl - xl)}
  SSE}


sem.error.f <- function(lambda, env) {
  SSE <- sem_error_sse(lambda, env)
  n <- get("n", envir=env)
  s2 <- SSE/n
  ldet <- spatialreg::do_ldet(lambda, env)
  ret <- (ldet + (1/2)*get("sum_lw", envir=env) - ((n/2)*log(2*pi)) - 
            (n/2)*log(s2) - (1/(2*(s2)))*SSE)
  if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret, " Jacobian:", ldet, " SSE:", SSE, "\n")
  assign("f_calls", get("f_calls", envir=env)+1L, envir=env)
  ret
}


for_ridge_sem<-function(data,formula, beta_ridge, scale = FALSE){
  etype<-"error"
  Durbin<-FALSE
  zero.policy=TRUE
  method<-"eigen"   #Also "eigenw" and other options see ML_models.Rd
  
  data_sd<-as(data, "Spatial")
  coords = coordinates(data_sd)
  
  Neigh <- poly2nb(data_sd, row.names = NULL, queen=TRUE)
  covid.lags <- nblag(Neigh, 2)
  dlist_1 <- nbdists(covid.lags[[1]], coordinates(coords), longlat=TRUE)
  inv_dlist_1 <- lapply(dlist_1, function(x) 1/(x))
  listw_1  <- nb2listw(covid.lags[[1]],style="B",glist =inv_dlist_1, zero.policy = TRUE)
  mat_w <- listw2mat(listw_1)
  mat_w <- mat_w/eigen(mat_w)$values[1]
  listw <-mat2listw(mat_w)
  
  
  
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
  y <- scale(y, scale=scale) #center response
  x <- scale(x, scale=scale) #center variables
  x <- x[,-1]
  n <- nrow(x)
  weights <- rep(as.numeric(1), n)
  
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
  
 
  sum_lw <- sum(log(weights))
  sw <- sqrt(weights)
  
  
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
  assign("sum_lw", sum_lw, envir=env)
  assign("sw", sw, envir=env)
  assign("verbose", FALSE, envir = env)
  interval <- spatialreg::jacobianSetup(method, env, con, pre_eig=con$pre_eig,
                                        trs=trs, interval=interval)
  assign("interval", interval, envir=env)
  
  
  opt <- optimize(sem.error.f, interval=interval, 
                  maximum=TRUE, tol=con$tol.opt, env=env)
  
  lambda <- opt$maximum
  
  
  # Loglikelihood test for lambda significance 
  if (is.null(beta_ridge)){
    yl <- y - lambda * wy
    yl <- sw * yl
    SSE <- crossprod(yl) 
    s2 <- SSE/n
    ldet <- spatialreg::do_ldet(lambda, env)
    ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2) - n/2)
    ret_0 <- ( - ((n/2)*log(2*pi)) - (n/2)*log(s2) - n/2)
  } else { 
    yl <- y - lambda *wy
    yl <- sw * yl
    xl <- x%*%beta_ridge - lambda *WX%*%beta_ridge
    xl <- sw * xl
    SSE <-  crossprod(yl - xl)
    s2 <- SSE/n
    ldet <- spatialreg::do_ldet(lambda, env)
    ret <- (ldet - ((n/2)*log(2*pi)) -  (n/2)*log(s2) - n/2)
    ret_0 <- ( - ((n/2)*log(2*pi)) - (n/2)*log(crossprod(y-x%*%beta_ridge)/n) - n/2)
  }
  
  
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
rrsem <- function(data, model, buffer, plot = FALSE, scale=FALSE){
  
  
  
  data.sp<-as(data,'Spatial') #transform to sp object for some functions
  coords<-coordinates(data.sp)
  mf <- lm(model, data,  method="model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.extract(mf, "response")
  y <- scale(y, scale=scale) #center response
  x <- model.matrix(mt, mf)
  x<-  scale(x, scale=scale) #center variables
  x <- x[,-1]
  
  Neigh <- poly2nb(data.sp, row.names = NULL, queen=TRUE)
  n <- nrow(x)
  
  # Step 1----
  # Estimate gamma and beta_ridge for Y=X*beta_ridge + epsilon
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  lambda_int<-n*fit$lambda
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
      
      
      train.data <- data[data$SLOO=='Train',]
      test.data  <- data[data$SLOO=='Test',]
      st_geometry(test.data) <- NULL
      predictorvarnames<-all.vars(model[[3]])
      
      x_train <- x[data$SLOO=='Train',]
      y_train <- y[data$SLOO=='Train']
      p <- ncol(test.data %>% dplyr::select(predictorvarnames))+1
      
      x_test <- x[data$SLOO=='Test',]
      y_test <- y[data$SLOO=='Test']
      
      ldet=0
      
      ridge_coef<-function(lambda){
        
        
        x<-x_train
        
        coefs<-  drop(solve(crossprod(x_train) + diag(lambda, (p-1)), crossprod(x_train, y_train)))
        
        
        y_pred = x%*%coefs
        
        SSE <- crossprod(y_pred-y_train)
        s2 <- SSE/nrow(x_train)
        
        
        x<-x_test
        
        logLik <- -(1/2)*log(2*pi) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x%*%coefs)
        
        #res <- list(coefs=coefs, s2 = s2, lambda = lambda, logLik = logLik)
        
        res <- logLik
        
        return(res)
      }
      
      
      logLik_sem[k,] <- sapply(lambda_int, ridge_coef)
      
      
  }
  #Compute MSE and select best lambda sem
  #cvraw_sem=(y-predmat_sem)^2 #antes era
  mlogLik_sem= colMeans(logLik_sem)
  
  plotdat_sem<-tibble(logLik=mlogLik_sem, lambda=lambda_int)
  lL_sem<- max(plotdat_sem$logLik)
  lambda_sem<- plotdat_sem$lambda[plotdat_sem$logLik==lL_sem]
  # coef_sem <- predict(fit,x,type = "coefficients", s=lambda_sem)
  # coef_sem <- coef_sem@x
  p<- ncol(x)+1
  coef_sem <- drop(solve(crossprod(x) + diag(lambda_sem, (p-1)), crossprod(x, y)))
  predictorsnames<-all.vars(model[[3]])
  names(coef_sem)<- c(predictorsnames )
  
  
  
  
  
  j=1
  vector_lambda <- numeric(100) # to store the lambda values
  matrix_beta <- matrix(NA,100,p-1) # to store the beta coefficients
  tol=0.000001 # tolerance to assess convergence of the algorithm
  delta_lambda=1 # initialize the element related to lambda convergence in the condition
  delta_beta=1 # initialize the element related to beta convergence in the condition
  
  beta_ridge = coef_sem
  
  # Step 2----
  param_ridge_sem <- for_ridge_sem(data, model, beta_ridge, scale=scale) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
  
  x_ridge<-param_ridge_sem$x
  
  
  y_ridge<-param_ridge_sem$y
  
  lambda<-param_ridge_sem$lambda
  ldet<-param_ridge_sem$ldet
  s2<-param_ridge_sem$s2
  fit <- glmnet(x_ridge, y_ridge, alpha = 0, standardize = FALSE, intercept=FALSE)
  lambda_int<-n*fit$lambda
  
  while ((delta_lambda > tol)||(delta_beta > tol)) {
    
    
    
    #Step 3----
    x_ridge<-param_ridge_sem$x
    
    
    y_ridge<-param_ridge_sem$y
    
    lambda<-param_ridge_sem$lambda
    ldet<-param_ridge_sem$ldet
    s2<-param_ridge_sem$s2
    fit <- glmnet(x_ridge, y_ridge, alpha = 0, standardize = FALSE, intercept=FALSE)
    lambda_int<-n*fit$lambda
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
        
        
        train.data <- data[data$SLOO=='Train',]
        test.data  <- data[data$SLOO=='Test',]
        
        x_train <- x_ridge[data$SLOO=='Train',]
        y_train <- y_ridge[data$SLOO=='Train',]
        p<-length(predictorvarnames)+1
        x_test <- x_ridge[data$SLOO=='Test',]
        y_test <- y_ridge[data$SLOO=='Test']
        
        # fit <- glmnet(x_train, y_train, alpha = 0, lambda = lambda_int, standardize = FALSE)
        # 
        # st_geometry(test.data) <- NULL
        # predictorvarnames<-all.vars(model[[3]])
        # p<-length(predictorvarnames)
        # x_test <- t(as.matrix(x_ridge[data$SLOO=='Test',]))
        # y_test <-y_ridge[data$SLOO=='Test',]
        # colnames(x_test) <- colnames(test.data %>% select(predictorvarnames))
        # 
        # predmat_sem[i,] <- predict(fit, newx = x_test,
        #                            type="response")
        # 
        
        ridge_coef<-function(lambda){
          
          
          x<-x_train
          
          coefs<-  drop(solve(crossprod(x_train) + diag(lambda, (p-1)), crossprod(x_train, y_train)))
          
          
          y_pred = x%*%coefs
          
          
          x<-x_test
          
          logLik <- -(1/2)*log(2*pi) -(1/2)*log(s2)+ ldet -(1/(2*s2))*crossprod(y_test-x%*%coefs)
          
          #res <- list(coefs=coefs, s2 = s2, lambda = lambda, logLik = logLik)
          
          res <- logLik
          
          return(res)
        }
        
        
        
        logLik_sem[i,] <- sapply(lambda_int, ridge_coef)
        
        
    }
    
    #cvraw_sem=(y_ridge-predmat_sem)^2 # antes era 
    mlogLik_sem= colMeans(logLik_sem)
    
    plotdat_sem<-tibble(logLik=mlogLik_sem, lambda=lambda_int)
    lL_sem<- max(plotdat_sem$logLik)
    lambda_sem<- plotdat_sem$lambda[plotdat_sem$logLik==lL_sem]
    # coef_sem <- predict(fit,x_ridge,type = "coefficients", s=lambda_sem)
    # coef_sem <- coef_sem@x
    coef_sem <- drop(solve(crossprod(x_ridge) + diag(lambda_sem, (p-1)), crossprod(x_ridge, y_ridge)))
    predictorsnames<-all.vars(model[[3]])
    names(coef_sem)<- c(predictorsnames )
    W <- param_ridge_sem$W
    
    # lL_sem_min<- min(plotdat_sem$logLik)
    # band <- (lL_sem-lL_sem_min)/4
    # new_max_lambda <-plotdat_sem%>%filter(logLik >lL_sem-band)%>%
    #   filter(lambda==max((lambda)))%>%select(lambda)%>%as.numeric()
    # new_min_lambda <-plotdat_sem%>%filter(logLik >lL_sem-band)%>%
    #   filter(lambda==min((lambda)))%>%select(lambda)%>%as.numeric()
    # lambda_int=seq(from=new_min_lambda, to=new_max_lambda, length.out = 100)
    
    # # Sum of Squares Total and Error
    # sst <- sum((y - mean(y))^2)
    # sse <- sum((y_predicted - y)^2)
    # 
    # # R squared
    # rsq <- 1 - sse / sst
    # pseudoR2<-cor(y,y_predicted)^2
    
    beta_ridge = coef_sem
    
    
    # Step 4----
    param_ridge_sem <- for_ridge_sem(data, model, beta_ridge, scale=scale) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
    
    x_ridge <-param_ridge_sem$x
    lambda <- param_ridge_sem$lambda
    W <- param_ridge_sem$W
    
    y_predicted <- solve((diag(1,n)-lambda*W),(x_ridge%*%coef_sem))
    
    #Stop criterion
    if(j==1){
      lambda_new=param_ridge_sem$lambda
      delta_lambda=lambda_new
      
      beta_new=coef_sem
      delta_beta=beta_new
      
    }else{
      lambda_old=lambda_new
      lambda_new=param_ridge_sem$lambda
      delta_lambda=lambda_old-lambda_new
      
      beta_old=beta_new
      beta_new=coef_sem
      delta_beta=sqrt(sum((beta_old-beta_new)^2))
    }
    
    vector_lambda[j]<-lambda_new
    matrix_beta[j,]<-beta_new
    j=j+1
  }
  
  
  if (plot == TRUE){
    
    beta =  function(l){
      coefs <- drop(solve(crossprod(x_ridge) + diag(l, (p-1)), crossprod(x_ridge, y_ridge)))
      return(coefs)
    }
    
    plotdat <-sapply(lambda_int, beta)
    plotdat<-melt(plotdat, id.vars=c("Variable"))
    
    
    
    colnames(plotdat)<-c("Variable","Model","Coefficients")
    plotdat$lambda<-rep(lambda_int, each=(ncol(x)))
    
    
    
    p_coef_path<-ggplot(plotdat,aes(x=log(lambda), y=Coefficients,color=Variable) )+
      geom_line()+
      theme_classic()+
      ylab(NULL)+
      geom_vline(aes(xintercept = log(lambda_sem)), color="black")+
      geom_text(data = plotdat %>% filter(lambda == last(lambda)), aes(label = Variable, 
                                                                       x = log(lambda), 
                                                                       y = Coefficients, 
                                                                       color = Variable)) + 
      guides(color = "none") +
      xlab( expression(paste("log(", gamma,")")))#+
    # ggtitle( bquote("SEM" ~ beta[ridge] ~ " coefficients paths  vs regularization parameter"))
    
    
    p_logLik<-ggplot(plotdat_sem, aes(y=logLik, x=log(lambda)))+
      geom_line()+
      theme_classic()+
      geom_vline(aes(xintercept = log(lambda_sem), color="red"))+
      guides(color = "none") +
      xlab( expression(paste("log(", gamma,")")))#+
    # ggtitle("SEM ridge logLikelihood vs regularization parameter")
    
    p_coef_path
    
    p_logLik
    
  }else{plotdat=plotdat_sem}
  
  
  vector_lambda <- vector_lambda[1:j-1] 
  matrix_beta <- matrix_beta[1:j-1,]
  results<-list(Model = model,
                lambda = param_ridge_sem$lambda,
                Coefficients = coef_sem,
                gamma = lambda_sem,
                x=x_ridge,
                predicted_values = y_predicted,
                convergence_lambda = vector_lambda,
                convergence_beta = matrix_beta,
                iterations=j-1,
                pvalue=param_ridge_sem$pvalue,
                statistic=param_ridge_sem$statistic,
                plotdat=plotdat,
                plotdat_sem=plotdat_sem)
  
  return(results)
} 

# Ridge regression SLOO----
ridge.sloo.lr<- function(data, model, buffer, plot){
  
  
  n<-nrow(data)
  data.sp<-as(data,'Spatial') #transform to sp object for some functions
  
  coords = coordinates(data.sp)
  
  Neigh <- poly2nb(data.sp, row.names = NULL, queen=TRUE)
  
  
  st_geometry(data) <- NULL
  predictorvarnames<-all.vars(model[[3]])
  responsevarname<-all.vars(model[[2]])
  x <- as.matrix(data %>% dplyr::select(predictorvarnames))
  
  x<- scale(x, scale=FALSE) #center variables
  y <- as.matrix(data %>% dplyr::select(responsevarname))
  
  y<- scale(y, scale=FALSE) #center response
  
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE, intercept=FALSE)
  lambda_int<-n*fit$lambda
  predmat_lr<-predict(fit,x)
  
  
  
  
  
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
      
      
      train.data <- data[data$SLOO=='Train',]
      test.data  <- data[data$SLOO=='Test',]
      
      
      
      #Ridge linear regression
      
      x_train <- as.matrix(train.data %>% dplyr::select(predictorvarnames))
      y_train <- as.matrix(train.data %>% dplyr::select(responsevarname))
      
      x_train <- scale(x_train, scale=FALSE) #center variables
      y_train <- scale(y_train, scale=FALSE) #center response
      fit <- glmnet(x_train, y_train, alpha = 0, standardize = FALSE, intercept=FALSE)
      
      
      predmat_lr[k,] <- predict(fit, newx = as.matrix(test.data %>% dplyr::select(predictorvarnames)),
                                type="response")
      
      
  }
  #Compute MSE and select best lambda lr
  fit <- glmnet(x, y, alpha = 0, standardize = FALSE)
  cvraw_lr=sweep(predmat_lr, 1, y)^2
  rmse_lr=apply(cvraw_lr,2,function(x) sqrt(mean(x)))
  plotdat_lr<-tibble(RMSE=rmse_lr, lambda=lambda_int)
  RMSE_lr<- min(plotdat_lr$RMSE)
  lambda_lr<- plotdat_lr$lambda[plotdat_lr$RMSE==min(plotdat_lr$RMSE)]
  coef_lr <- drop(solve(crossprod(x), crossprod(x, y)))
  y_predicted <- x%*%coef_lr
  predictorsnames<-all.vars(model[[3]])
  names(coef_lr)<- c(predictorsnames )
  
  if(plot){
    
    plotdat<- as.data.frame(as.matrix(fit[["beta"]]))
    plotdat$Variable<-predictorsnames
    plotdat<-melt(plotdat, id.vars=c("Variable"))
    
    
    colnames(plotdat)<-c("Variable","Model","Coefficients")
    plotdat$lambda<-rep(fit[["lambda"]], each=ncol(x))
    plotdat$dev.ratio<-rep(fit[["dev.ratio"]], each=ncol(x))
    
    
    
    p11<-ggplot(plotdat,aes(x=log(lambda), y=Coefficients,color=Variable) )+
      geom_line()+
      theme_classic()+
      ylab(NULL)+
      geom_vline(aes(xintercept = log(lambda_lr)), color="black")+
      geom_text(data = plotdat %>% filter(lambda == last(lambda)), aes(label = Variable, 
                                                                       x = log(lambda) - 0.5, 
                                                                       y = Coefficients, 
                                                                       color = Variable)) + 
      guides(color = FALSE) +
      xlab( expression(paste("log(", gamma,")")))
    
    
    p2<-ggplot(plotdat_lr, aes(y=RMSE, x=log(lambda)))+
      geom_line()+
      theme_classic()+
      geom_vline(aes(xintercept = log(lambda_lr), color="red"))+
      guides(color = FALSE) +
      xlab( expression(paste("log(", gamma,")")))
    
    p11
    p2}
  results<-list(Model = model,
                Coefficients = coef_lr,
                gamma = lambda_lr,
                predicted_values=y_predicted)
  
  return(results)
}

#Permutated F-test ----
#Function to perform a permutated F-test for ridge SAR, SEM and linear regression model
# data sf object
# model  formula object with the model to analyse
# B number of permutations to perform for each variable
#Output
#Dataframe with estimated p-values for each variable in each model

perm_f_test_parallel <- function(data, model, B, buffer){
  
  set.seed(1L)
  n <- nrow(data) #Number of data points
  
  v <- length(all.vars(model[[3]])) #Number of variables
  
  p_f_test <- data.frame(SAR = rep(NA,v), SEM = rep(NA,v), OLS = rep(NA,v),
                         row.names = all.vars(model[[3]]))
  
  predictorvarnames<-all.vars(model[[3]])
  
  responsevarname<-all.vars(model[[2]])
  
  mt <- terms(model, data = data)
  mf <- lm(model, data, 
           method="model.frame")
  x_lr <- model.matrix(mt, mf)[,-1]
  y_lr <- model.extract(mf, "response") 
  
  #SAR
  scale=FALSE
  ridge_sar <- rrsar(data,model,buffer, scale=scale, plot=FALSE) #Do ridge
  
  
  y_sar_est <- ridge_sar$predicted_values #Obtain estimated values y_1
  
  #SEM
  scale=FALSE
  ridge_sem <- rrsem(data,model,buffer, scale=scale, plot=FALSE) #Do ridge
  
  
  y_sem_est <- ridge_sem$predicted_values #Obtain estimated values y_1
  
  cores<-detectCores()
  
  cl <- makeCluster(cores)
  
  env <- new.env(parent=globalenv())
  env <- new.env()
  assign("data", data, envir=env)
  assign("model", model, envir=env)
  assign("y_sar_est", y_sar_est, envir=env)
  assign("y_sem_est", y_sem_est, envir=env)
  # assign("y_lr_est", y_lr_est, envir=env)
  assign("buffer", buffer, envir=env)
  assign("y_lr", y_lr, envir=env)
  assign("n", n, envir=env)
  assign("B", B, envir=env)
  assign("x_lr", x_lr, envir=env)
  assign("verbose", FALSE, envir=env)
  
  clusterExport(cl, c("rrsar", "rrsem", "for_ridge_sar", 
                      "for_ridge_sem", "ridge.sloo.lr", "sar.lag.mixed.f","perm_sar",
                      "perm_sem","perm_lr",
                      "sem_error_sse", "sem.error.f", "data", "model",
                      "y_sar_est", "y_sem_est", #"y_lr_est", 
                      "buffer", "y_lr", "n", "B", "x_lr"),
                envir = env)
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl,library("grid"))
  clusterEvalQ(cl,library("lctools"))
  clusterEvalQ(cl,library("FactoMineR"))
  clusterEvalQ(cl,library("factoextra"))
  clusterEvalQ(cl,library("ggrepel"))
  clusterEvalQ(cl,library("reshape"))
  clusterEvalQ(cl,library("geojsonio"))
  clusterEvalQ(cl,library("gridExtra"))
  clusterEvalQ(cl,library("sp"))
  clusterEvalQ(cl,library("Hmisc"))
  clusterEvalQ(cl,library("tidyverse"))
  clusterEvalQ(cl,library("sf"))
  clusterEvalQ(cl,library("spatialreg"))
  clusterEvalQ(cl,library("glmnet"))
  clusterEvalQ(cl,library("StatMeasures"))
  clusterEvalQ(cl,library("spdep"))
  clusterEvalQ(cl, library("mctest"))
  clusterEvalQ(cl,library("GGally"))
  
  p_f_test <- parLapply(cl, predictorvarnames, fun=permute_variable_j)
  
  stopCluster(cl)
  
  res <- data.frame(do.call(rbind, p_f_test))
  
  return(res)
}


permute_variable_j <- function(j){
  
  p_f_test <- data.frame(SAR = NA, SEM = NA, OLS = NA)
  rownames(p_f_test)<-j
  
  #Compute standard F statistic
  data_j<-data%>%dplyr::select(-j) #Eliminate variable j from data set
  model_j<-update(model, paste("~ . -",j)) #Eliminate variable j from formula
  
  #SAR
  scale=FALSE
  ridge_sar_j <- rrsar(data_j,model_j,buffer, scale=scale) #Do ridge
  
  
  y_sar_pred_j <- ridge_sar_j$predicted_values #Obtain estimated values  y_0
  
  F_sar_j   <-  sqrt(sum((y_sar_est-y_sar_pred_j)^2))/sqrt(sum((y_lr-y_sar_est)^2))
  
  
  
  
  #SEM
  scale=FALSE
  ridge_sem_j <-  rrsem(data_j,model_j,buffer, scale=scale) #Do ridge
  
  
  y_sem_pred_j <- ridge_sem_j$predicted_values #Obtain estimated values y_0
  
  
  F_sem_j   <- sqrt(sum((y_sem_est-y_sem_pred_j)^2))/sqrt(sum((y_lr-y_sem_est)^2))
  
  
  #Obtain B permutations of variable j

  data_p <- data
  st_geometry(data_p) <- NULL
  N <- t(replicate(B+B*0.1, sample(as.matrix(data_p%>%dplyr::select(j)),n, replace = TRUE)))
  out <- N[!(duplicated(N)), ][1:B,]
  out <- split(out, row(out))#transform to list
  x_lr_perm<-x_lr
  
  F_j_sar <- sapply(out, perm_sar, varname=j,  y_sar_pred_j, y_lr, F_sar_j)
  
  F_j_sem <- sapply(out, perm_sem, varname=j,  y_sem_pred_j, y_lr, F_sem_j)
  
  # F_j_lr <- sapply(out, perm_lr, varname=j,  y_lr_pred_j, y_lr, F_lr_j)
  
  p_f_test[j,"SAR"] <- format(sum(F_j_sar)/B, digits = 4)
  p_f_test[j,"SEM"] <- format(sum(F_j_sem)/B, digits = 4)
  # p_f_test[j,"OLS"]  <- format(sum(F_j_lr)/B, digits = 4)
  return(p_f_test)
}


perm_sar <-function(x, varname, y_sar_pred_j, y_lr, F_sar_j){
  j=varname
  data_perm <- data
  data_perm[,j]<-x
  #SAR
  scale=FALSE
  ridge_sar_perm <- rrsar(data_perm,model,buffer, scale=scale) #Do ridge
  
  
  y_sar_perm <- ridge_sar_perm$predicted_values #Obtain estimated values
  
  f<-(sqrt(sum((y_lr-y_sar_pred_j)^2))-sqrt(sum((y_sar_perm-y_lr)^2)))/sqrt(sum((y_lr-y_sar_perm)^2))
  F_j_sar<- f >= F_sar_j
  
  return(F_j_sar)
  
}

perm_sem <-function(x, varname, y_sem_pred_j, y_lr, F_sem_j){
  j=varname
  data_perm <- data
  data_perm[,j]<-x
  
  #SEM
  scale=FALSE
  ridge_sem_perm <- rrsem(data_perm,model,buffer, scale=scale) #Do ridge
  y_sem_perm <- ridge_sem_perm$predicted_values #Obtain estimated values
  f<-(sqrt(sum((y_lr-y_sem_pred_j)^2))-sqrt(sum((y_sem_perm-y_lr)^2)))/sqrt(sum((y_lr-y_sem_perm)^2))
  F_j_sem<- f >= F_sem_j
  
  return(F_j_sem)
  
}

perm_lr <- function(x, varname, y_lr_pred_j, y_lr, F_lr_j){
  j=varname
  data_perm <- data
  data_perm[,j]<-x
  
  
  #LR
  ridge_lr_perm <- ridge.sloo.lr(data_perm,model,buffer,plot=FALSE) #Do ridge
  
  
  y_lr_perm <- ridge_lr_perm$predicted_values #Obtain estimated values
  
  f<-sqrt((sum((y_lr-y_lr_pred_j)^2))-sqrt(sum((y_lr_perm-y_lr)^2)))/sqrt(sum((y_lr-y_lr_perm)^2))
  F_j_lr<-f>=F_lr_j
  
  return(F_j_lr)
  
}
