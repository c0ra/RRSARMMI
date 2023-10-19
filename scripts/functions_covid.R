
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
sar_estimation<-function(data, formula, beta_0, beta_l, scale= FALSE){
  
  verbose = FALSE
  
  method<-"eigen"   
  
  data_sd<-as(data, "Spatial")
  coords = coordinates(data_sd)
  
  
  # Create spatial neighborhood weights
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
  control=list() # Initialize a control list for optimization
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
  
  
  
  # Compute loglikelihood rate for test of significance of the spatial dependence parameter
  e.lm.null <- y - x%*%beta_0
  e.lm.w <- wy- x%*%beta_l
  
  SSE <-crossprod(e.lm.null-rho*e.lm.w)
  
  s2 <- SSE/n
  
  ldet <- spatialreg::do_ldet(rho, env)
  
  ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2)) # loglikelihood 
  
  ret_0 <- -((n/2)*log(2*pi)) - (n/2)*log(crossprod(e.lm.null)/n) # loglikelihood when rho = 0
  
  ldet <- spatialreg::do_ldet(rho, env)
  
  #Output parameters
  LR <- 2*(ret-ret_0) # likelihood rate
  p_value <- pchisq(LR, df=1, lower.tail = FALSE) # p-value
  y <- y - rho*wy # filtered dependent variable
  x <- x
  ldet <- ldet
  sar_arguments <- list(rho=rho, y = y, x = x, W=mat_w, pvalue= p_value, statistic= LR, ldet=ldet, s2=s2)
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
# - plot: Wheter to return data sets and plot coefficients paths and loglikelihood curve
# - scale: Logical value indicating whether to scale the variables
rrsar <- function(data, model, buffer, plot = FALSE, scale=FALSE){
  
  # Store number of observations in n
  n<-nrow(data)
  
  
  data.sp<-as(data,'Spatial') #transform to sp object
  
  
  # Extract dependent variable and covariates
  mt <- terms(model, data = data) 
  mf <- lm(model, data, 
           method="model.frame")
  
  y <- model.extract(mf, "response")
  
  x <- model.matrix(mt, mf)
  
  # Create spatial neighborhood weights
  coords<-coordinates(data.sp)
  Neigh <- poly2nb(data.sp, row.names = NULL, queen=TRUE)
  covid.lags <- nblag(Neigh, 2)
  dlist_1 <- nbdists(covid.lags[[1]], coordinates(coords), longlat=TRUE)
  inv_dlist_1 <- lapply(dlist_1, function(x) 1/(x))
  listw_1  <- nb2listw(covid.lags[[1]],style="B",glist =inv_dlist_1, zero.policy = TRUE)
  mat_w <- listw2mat(listw_1)
  mat_w <- mat_w/eigen(mat_w)$values[1]
  listw <-mat2listw(mat_w)
  
  # Obtain lagged dependent variable
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
  plotdat_sar_beta_0<-tibble(logLik=mlogLik_beta_0, gamma=gamma_int)
  plotdat_sar_beta_l<-tibble(logLik=mlogLik_beta_l, gamma=gamma_int)
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
  ridge_sar <- sar_estimation(data, model, beta_0=coef_0_ridge, beta_l=coef_l_ridge, scale=scale)
  
  # Step 3: Update rho, filtered dependent variable and perform ridge regression----
  y_filtered <- ridge_sar$y
  rho <- ridge_sar$rho
  s2 <- ridge_sar$s2
  W <- ridge_sar$W
  
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
  
  plotdat_sar <-tibble(logLik=mlogLik_sar, gamma=gamma_int)
  
  # Store predicted values with final parameters
  y_predicted <- solve((diag(1,n)-rho*W),(x%*%coef_sar))
  
  
  # Plot coefficients paths and log-likelihood curve
  if (plot == TRUE){
    
    beta =  function(l){
      coefs <-  drop(solve(crossprod(x) + diag(l, p), crossprod(x, y_filtered)))
      return(coefs)
    }
    plotdat <-sapply(gamma_int, beta)
    plotdat<-melt(plotdat, id.vars=c("Variable"))
    
    
    colnames(plotdat)<-c("Variable","Model","Coefficients")
    plotdat$gamma<-rep(gamma_int, each=(ncol(x)))
    
    
    
    p_coef_path<-ggplot(plotdat,aes(x=log(gamma), y=Coefficients,color=Variable) )+
      geom_line()+
      theme_classic()+
      ylab(NULL)+
      geom_vline(aes(xintercept = log(gamma_sar)), color="black")+
      geom_text(data = plotdat %>% filter(gamma == last(gamma)), aes(label = Variable, 
                                                                     x = log(gamma), 
                                                                     y = Coefficients, 
                                                                     color = Variable)) + 
      guides(color = "none") +
      xlab( expression(paste("log(", gamma,")")))
    
    
    p_logLik<-ggplot(plotdat_sar, aes(y=logLik, x=log(gamma)))+
      geom_line()+
      theme_classic()+
      geom_vline(aes(xintercept = log(gamma_sar), color="red"))+
      guides(color = "none") +
      xlab( expression(paste("log(", gamma,")")))
    
    p_coef_path
    
    p_logLik
    
  }else{plotdat=plotdat_sar}
  
  # Create output
  results<-list(Model = model,
                rho = ridge_sar$rho,
                Coefficients = coef_sar,
                coef_0_ridge = coef_0_ridge,
                coef_l_ridge = coef_l_ridge,
                gamma = gamma_sar,
                gamma_beta_0 = gamma_beta_0,
                gamma_beta_l = gamma_beta_l,
                predicted_values = y_predicted,
                pvalue=ridge_sar$pvalue,
                statistic=ridge_sar$statistic,
                plotdat=plotdat,
                plotdat_sar=plotdat_sar,
                plodat_sar_beta_0 = plotdat_sar_beta_0,
                plotdat_sar_beta_l = plotdat_sar_beta_l)
  
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
estimation_sem<-function(data,formula, beta_ridge, scale = FALSE){
  
  etype<-"error"
  Durbin<-FALSE
  zero.policy=TRUE
  method<-"eigen" 
  
  data_sd<-as(data, "Spatial")
  coords = coordinates(data_sd)
  
  # Create spatial neighborhood weights
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
  
  # Extract dependent variable and covariates
  mt <- terms(formula, data = data)
  mf <- lm(formula, data,  method="model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  n <- nrow(x)
  weights <- rep(as.numeric(1), n)
  
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
  
  
  ### Likelihood rate test for lambda significance 
  # Compute filtered variables
  yl <- y - lambda *wy
  xl <- x%*%beta_ridge - lambda *WX%*%beta_ridge
  
  
  # Compute squared sum of errors
  SSE <-  crossprod(yl - xl)
  
  # Compute sample variance
  s2 <- SSE/n
  
  
  ldet <- spatialreg::do_ldet(lambda, env)
  
  
  ret <- (ldet - ((n/2)*log(2*pi)) -  (n/2)*log(s2))  # loglikelihood 
  ret_0 <- ( - ((n/2)*log(2*pi)) - (n/2)*log(crossprod(y-x%*%beta_ridge)/n)) # loglikelihood when lambda = 0
  
  
  
  #Output parameters
  LR <- 2*(ret-ret_0) # likelihood rate
  p_value <- pchisq(LR, df=1, lower.tail = FALSE) #p value
  y <- y - lambda*wy # filtered dependent variable
  x <- x - lambda*WX # filtered covariates
  ldet <- ldet
  
  sem_arguments <- list(lambda=lambda, y = y, x = x, W = mat_w, pvalue=p_value, statistic=LR, ldet = ldet, s2=s2)
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
# - plot: Wheter to return data sets and plot coefficients paths and loglikelihood curve
# - scale: Logical value indicating whether to scale the variables
rrsem <- function(data, model, buffer, plot = FALSE, scale=FALSE){
  
  
  # Extrat dependent variable and covariates from model
  mf <- lm(model, data,  method="model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  
  
  data.sp<-as(data,'Spatial') #transform to sp object to obtain neighbourhood
  coords<-coordinates(data.sp)
  Neigh <- poly2nb(data.sp, row.names = NULL, queen=TRUE) # to define buffer for SLOO
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
  param_ridge_sem <-estimation_sem(data, model, beta_ridge, scale=scale) 
  
  # Get filtered  dependent variable and filtered covariates
  x_filtered<-param_ridge_sem$x
  y_filtered<-param_ridge_sem$y
  
  lambda<-param_ridge_sem$lambda
  ldet<-param_ridge_sem$ldet
  s2<-param_ridge_sem$s2
  
  
  
  
  
  #Step 3: Perform ridge for the filtered variables----
  
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
        logLik <- -(1/2)*log(2*pi) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        
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
  beta_ridge  <- drop(solve(crossprod(x_filtered) + diag(gamma_sem, (p)), crossprod(x_filtered, y_filtered)))
  predictorsnames<-all.vars(model[[3]])
  names(beta_ridge )<- c(predictorsnames )
  
  
  
  
  
  
  # Step 4: Estimate lambda using the new ridge coefficients----
  param_ridge_sem <- estimation_sem(data, model, beta_ridge, scale=scale) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
  
  # Get the new filtered variables
  x_filtered<-param_ridge_sem$x
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
  plotdat_sem <-tibble(logLik=mlogLik_sem, gamma=gamma_int)
  
  # Update beta coefficients using the ridge estimate withe best gamma from the previous step
  coef_sem <- drop(solve(crossprod(x_filtered) + diag(gamma_sem, (p)), crossprod(x_filtered, y_filtered)))
  y_predicted <- solve((diag(1,n)-lambda*W),(x_filtered%*%coef_sem))
  
  
  # Plot coefficients paths and log-likelihood curve
  if (plot == TRUE){
    
    beta =  function(l){
      coefs <- drop(solve(crossprod(x_filtered) + diag(l, (p)), crossprod(x_filtered, y_filtered)))
      return(coefs)
    }
    
    plotdat <-sapply(gamma_int, beta)
    plotdat<-melt(plotdat, id.vars=c("Variable"))
    
    
    
    colnames(plotdat)<-c("Variable","Model","Coefficients")
    plotdat$gamma<-rep(gamma_int, each=(ncol(x)))
    
    
    
    p_coef_path<-ggplot(plotdat,aes(x=log(gamma), y=Coefficients,color=Variable) )+
      geom_line()+
      theme_classic()+
      ylab(NULL)+
      geom_vline(aes(xintercept = log(gamma_sem)), color="black")+
      geom_text(data = plotdat %>% filter(gamma == last(gamma)), aes(label = Variable, 
                                                                     x = log(gamma), 
                                                                     y = Coefficients, 
                                                                     color = Variable)) + 
      guides(color = "none") +
      xlab( expression(paste("log(", gamma,")")))#+
    
    
    p_logLik<-ggplot(plotdat_sem, aes(y=logLik, x=log(gamma)))+
      geom_line()+
      theme_classic()+
      geom_vline(aes(xintercept = log(gamma_sem), color="red"))+
      guides(color = "none") +
      xlab( expression(paste("log(", gamma,")")))#+
    
    p_coef_path
    
    p_logLik
    
  }else{plotdat=plotdat_sem}
  
  
  results<-list(Model = model,
                lambda = param_ridge_sem$lambda,
                Coefficients = coef_sem,
                gamma = gamma_sem,
                x=x_filtered,
                predicted_values = y_predicted,
                pvalue=param_ridge_sem$pvalue,
                statistic=param_ridge_sem$statistic,
                plotdat=plotdat,
                plotdat_sem=plotdat_sem)
  
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
  
  
  n <- nrow(data) #Number of data points
  
  v <- length(all.vars(model[[3]])) #Number of variables
  
  p_f_test <- data.frame(SAR = rep(NA,v), SEM = rep(NA,v), 
                         row.names = all.vars(model[[3]]))
  
  predictorvarnames<-all.vars(model[[3]])
  
  responsevarname<-all.vars(model[[2]])
  
  mt <- terms(model, data = data)
  mf <- lm(model, data, 
           method="model.frame")
  x_lr <- model.matrix(mt, mf)
  y_lr <- model.extract(mf, "response") 
  
  #SAR
  scale=FALSE
  ridge_sar <- rrsar(data,model,buffer, scale=scale, plot=FALSE) #Do ridge regresion for SAR
  
  
  y_sar_est <- ridge_sar$predicted_values #Obtain estimated values y_1
  
  #SEM
  scale=FALSE
  ridge_sem <- rrsem(data,model,buffer, scale=scale, plot=FALSE) #Do ridge regresion for SAR
  
  
  y_sem_est <- ridge_sem$predicted_values #Obtain estimated values y_1
  
  # Paralellize
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
  
  clusterExport(cl, c("rrsar", "rrsem", "sar_estimation", 
                      "estimation_sem", 
                      "sar.lag.mixed.f","perm_sar",
                      "perm_sem",
                      "sem.error.f", "data", "model",
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
  
  # Apply function that permutes each covariate B times  and obtains p-value estimate
  p_f_test <- parLapply(cl, predictorvarnames, fun=permute_variable_j)
  
  stopCluster(cl)
  
  res <- data.frame(do.call(rbind, p_f_test))
  
  return(res)
}

# Function to obtain p value of permutation test for variable j
permute_variable_j <- function(j){
  
  # Initialize dataframe to store p value for RRSAR and RRSEM
  p_f_test <- data.frame(SAR = NA, SEM = NA)
  rownames(p_f_test)<-j
  
  #Compute standard F statistic
  data_j<-data%>%select(-j) #Eliminate variable j from data set
  model_j<-update(model, paste("~ . -",j)) #Eliminate variable j from formula
  
  #SAR
  scale=FALSE
  ridge_sar_j <- rrsar(data_j,model_j,buffer, scale=scale) #Do RRSAR without the variable j
  
  
  y_sar_pred_j <- ridge_sar_j$predicted_values #Obtain estimated values  y_0
  
  F_sar_j   <-  (sqrt(crossprod(y_lr-y_sar_pred_j))^2/sqrt(crossprod(y_lr-y_sar_est))^2)-1 # Compute F-statistic
  
  
  
  
  #SEM
  scale=FALSE
  ridge_sem_j <-  rrsem(data_j,model_j,buffer, scale=scale) #Do RRSEM without the variable j
  
  
  y_sem_pred_j <- ridge_sem_j$predicted_values #Obtain estimated values y_0
  
  
  F_sem_j   <- (sqrt(crossprod(y_lr-y_sem_pred_j))^2/sqrt(crossprod(y_lr-y_sem_est))^2)-1 # Compute F-statistic
  
  
  #Obtain B permutations of variable j
  
  data_p <- data
  st_geometry(data_p) <- NULL
  N <- t(replicate(B+B*0.1, sample(as.matrix(data_p%>%select(j)),n, replace = TRUE)))
  out <- N[!(duplicated(N)), ][1:B,]
  out <- split(out, row(out))#transform to list
  x_lr_perm<-x_lr
  
  #Do RRSAR B times, one for each permutation of variable j,
  # compute approximation of F-statistic and compare with F_sar_j 
  F_j_sar <- sapply(out, perm_sar, varname=j,  y_sar_pred_j, y_lr, F_sar_j)
  
  #Do RRSEM B times, one for each permutation of variable j,
  # compute approximation of F-statistic and compare with F_sar_j 
  F_j_sem <- sapply(out, perm_sem, varname=j,  y_sem_pred_j, y_lr, F_sem_j)
  
  # Compute p value
  p_f_test[j,"SAR"] <- format(sum(F_j_sar)/B, digits = 4)
  p_f_test[j,"SEM"] <- format(sum(F_j_sem)/B, digits = 4)
  return(p_f_test)
}

# Function to compare f-statistic aproximation by permutation with its value  F_sar_j:
# x: permutation of variable of interest
# varname: name of the variable of interest
# y_sar_pred_j:(y_0) predicted values of model estimated without the variable of interest
# y_lr: (y_1) predicted model estimated with all the variables
# F_sar_j: F-statistic
perm_sar <-function(x, varname, y_sar_pred_j, y_lr, F_sar_j){
  
  # Substitute variable by its permutation
  j=varname
  data_perm <- data
  data_perm[,j]<-x
  #SAR
  scale=FALSE
  ridge_sar_perm <- rrsar(data_perm,model,buffer, scale=scale) #Do RRSAR with j variable permutated
  
  
  y_sar_perm <- ridge_sar_perm$predicted_values #Obtain estimated values
  
  f<-(sqrt(crossprod(y_lr-y_sar_pred_j))^2/sqrt(crossprod(y_lr-y_sar_perm))^2)-1 # compute approximation of F-statistic
  
  F_j_sar<- f >= F_sar_j # Compare approximation with F-statistic
  
  return(F_j_sar)
  
}

# Function to compare f-statistic aproximation by permutation with its value  F_sem_j:
# x: permutation of variable of interest
# varname: name of the variable of interest
# y_sem_pred_j:(y_0) predicted values of model estimated without the variable of interest
# y_lr: (y_1) predicted model estimated with all the variables
# F_sar_j: F-statistic
perm_sem <-function(x, varname, y_sem_pred_j, y_lr, F_sem_j){
  
  # Substitute variable by its permutation
  j=varname
  data_perm <- data
  data_perm[,j]<-x
  
  #SEM
  scale=FALSE
  ridge_sem_perm <- rrsem(data_perm,model,buffer, scale=scale) #Do RRSEM with j variable permutated
  y_sem_perm <- ridge_sem_perm$predicted_values #Obtain estimated values
  f<-(sqrt(crossprod(y_lr-y_sem_pred_j))^2/sqrt(crossprod(y_lr-y_sem_perm))^2)-1 # compute approximation of F-statistic
  F_j_sem<- f >= F_sem_j # Compare approximation with F-statistic
  
  return(F_j_sem)
  
}
