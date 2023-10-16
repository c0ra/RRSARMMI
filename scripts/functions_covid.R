adjust_ridge<-function(gamma){
  
  coefs <-  drop(solve(crossprod(x_train) + diag(gamma, p-1), crossprod(x_train, y_train)))
  
  
  y_pred = x_test%*%coefs
  
  return(y_pred)
}


beta =  function(l){
  coefs <- drop(solve(crossprod(x_filtered) + diag(l, (p-1)), crossprod(x_filtered, y_filtered)))
  return(coefs)
}

logLik_conditional<-function(gamma){
  
  
  x<-x_train
  
  coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p-1)), crossprod(x_train, y_train)))
  
  
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


sar_estimation<-function(data,formula,  beta_0, beta_l, scale= FALSE){
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
  
  

    
    e.lm.null <- y - x%*%beta_0
    
    e.lm.w <- wy- x%*%beta_l
    
    SSE <-crossprod(e.lm.null-rho*e.lm.w)
    
    s2 <- SSE/n
    
    ldet <- spatialreg::do_ldet(rho, env)
    
    ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2)) # loglikelihood 
    
    ret_0 <- -((n/2)*log(2*pi)) - (n/2)*log(crossprod(e.lm.null)/n) # loglikelihood when rho = 0

  
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
  
  # std_deviations <- st_drop_geometry(data) %>%
  #   select(LnPop, C3, SLivR, Work, Inac, A65pls, Emer, FDoc) %>%
  #   summarise_all(sd)
  # 
  # 
  # data <- data %>%
  #   mutate(
  #     # Standardize the variables
  #     LnPop = scale(LnPop),
  #     C3 = scale(C3),
  #     SLivR = scale(SLivR),
  #     Work = scale(Work),
  #     Inac = scale(Inac),
  #     A65pls = scale(A65pls),
  #     Emer = scale(Emer),
  #     FDoc = scale(FDoc),
  #     
  #     # Center the variable
  #     LnHosp = scale(LnHosp, scale = FALSE)
  #   )
  # 

  
  n<-nrow(data)
  data.sp<-as(data,'Spatial') #transform to sp object for some functions
  
  mt <- terms(model, data = data) 
  mf <- lm(model, data, 
           method="model.frame")
  
  y <- model.extract(mf, "response")
  
  x <- model.matrix(mt, mf)[,-1]
  
  coords<-coordinates(data.sp)
  Neigh <- poly2nb(data.sp, row.names = NULL, queen=TRUE)
  covid.lags <- nblag(Neigh, 2)
  dlist_1 <- nbdists(covid.lags[[1]], coordinates(coords), longlat=TRUE)
  inv_dlist_1 <- lapply(dlist_1, function(x) 1/(x))
  listw_1  <- nb2listw(covid.lags[[1]],style="B",glist =inv_dlist_1, zero.policy = TRUE)
  mat_w <- listw2mat(listw_1)
  mat_w <- mat_w/eigen(mat_w)$values[1]
  listw <-mat2listw(mat_w)
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
  plotdat_sar_beta_0<-tibble(logLik=mlogLik_beta_0, gamma=gamma_int)
  plotdat_sar_beta_l<-tibble(logLik=mlogLik_beta_l, gamma=gamma_int)
  gamma_beta_0<- gamma_int[which.max(mlogLik_beta_0)]
  gamma_beta_l<- gamma_int[which.max(mlogLik_beta_l)]
  

  coef_0_ridge <-  drop(solve(crossprod(x) + diag(gamma_beta_0, p), crossprod(x, y)))
  coef_l_ridge <-  drop(solve(crossprod(x) + diag(gamma_beta_l, p), crossprod(x, wy)))
  
  
  predictorsnames<-all.vars(model[[3]])
  names(coef_0_ridge)<- c(predictorsnames )
  names(coef_l_ridge)<- c(predictorsnames )
  
  # Step 2----
  ridge_sar <- sar_estimation(data, model, beta_0=coef_0_ridge, beta_l=coef_l_ridge, scale=scale)
  
  #Step 3----
  y_filtered <- ridge_sar$y
  rho <- ridge_sar$rho
  s2 <- ridge_sar$s2
  W <- ridge_sar$W
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
  coef_sar <- coef_sar
  plotdat_sar <-tibble(logLik=mlogLik_sar, gamma=gamma_int)
  
  y_predicted <- solve((diag(1,n)-rho*W),(x%*%coef_sar))
  
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


estimation_sem<-function(data,formula, beta_ridge, scale = FALSE){
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
  
  can.sim <- FALSE
  if (listw$style %in% c("W", "S")) 
    can.sim <- spatialreg::can.be.simmed(listw)
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
  assign("can.sim", can.sim, envir=env)
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

    yl <- y - lambda *wy
    yl <- sw * yl
    xl <- x%*%beta_ridge - lambda *WX%*%beta_ridge
    xl <- sw * xl
    SSE <-  crossprod(yl - xl)
    s2 <- SSE/n
    ldet <- spatialreg::do_ldet(lambda, env)
    ret <- (ldet - ((n/2)*log(2*pi)) -  (n/2)*log(s2))
    ret_0 <- ( - ((n/2)*log(2*pi)) - (n/2)*log(crossprod(y-x%*%beta_ridge)/n))

  
  
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
  
  # std_deviations <- data %>%
  #   select(LnPop, C3, SLivR, Work, Inac, A65pls, Emer, FDoc) %>%
  #   summarise_all(sd)
  # 
  # 
  # data <- data %>%
  #   mutate(
  #     # Standardize the variables
  #     LnPop = scale(LnPop),
  #     C3 = scale(C3),
  #     SLivR = scale(SLivR),
  #     Work = scale(Work),
  #     Inac = scale(Inac),
  #     A65pls = scale(A65pls),
  #     Emer = scale(Emer),
  #     FDoc = scale(FDoc),
  #     
  #     # Center the variable
  #     LnHosp = scale(LnHosp, scale = FALSE)
  #   )
  
  
  data.sp<-as(data,'Spatial') #transform to sp object for some functions
  coords<-coordinates(data.sp)
  mf <- lm(model, data,  method="model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  x <- x[,-1]
  
  Neigh <- poly2nb(data.sp, row.names = NULL, queen=TRUE)
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
        
      
        
        coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
        
        
        y_pred = x_train%*%coefs
        
        SSE <- crossprod(y_pred-y_train)
        s2 <- SSE/nrow(x_train)
        
        
        logLik <- -(1/2)*log(2*pi) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x_test%*%coefs)
        
        
        res <- logLik
        
        return(res)
      }
      
      
      logLik_sem[k,] <- sapply(gamma_int, logLik_conditional)
      
      
  }
  #Compute MSE and select best lambda sem
  #cvraw_sem=(y-predmat_sem)^2 #antes era
  mlogLik_sem= colMeans(logLik_sem)
  gamma_sem<-  gamma_int[which.max(mlogLik_sem)]
  


  p<- ncol(x)
  beta_ridge <- drop(solve(crossprod(x) + diag(gamma_sem, (p)), crossprod(x, y)))
  predictorsnames<-all.vars(model[[3]])
  names(  beta_ridge)<- c(predictorsnames )
  
  
  
  # Step 2----
  param_ridge_sem <-estimation_sem(data, model, beta_ridge, scale=scale) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
  
  x_filtered<-param_ridge_sem$x
  
  
  y_filtered<-param_ridge_sem$y
  
  lambda<-param_ridge_sem$lambda
  ldet<-param_ridge_sem$ldet
  s2<-param_ridge_sem$s2
 

    
    
    
    #Step 3----
    fit <- glmnet(x_filtered, y_filtered, alpha = 0, standardize = FALSE, intercept=FALSE)
    gamma_int<-fit$lambda
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
          
          
          x<-x_train
          
          coefs<-  drop(solve(crossprod(x_train) + diag(gamma, (p)), crossprod(x_train, y_train)))
          
          
          y_pred = x%*%coefs
          
          
          x<-x_test
          
          logLik <- -(1/2)*log(2*pi) -(1/2)*log(s2) -(1/(2*s2))*crossprod(y_test-x%*%coefs)
          
          
          res <- logLik
          
          return(res)
        }
        
        
        
        logLik_sem[i,] <- sapply(gamma_int, logLik_conditional)
        
        
    }
    
  
    mlogLik_sem= colMeans(logLik_sem)
    gamma_sem<-  gamma_int[which.max(mlogLik_sem)]
  

    beta_ridge  <- drop(solve(crossprod(x_filtered) + diag(gamma_sem, (p)), crossprod(x_filtered, y_filtered)))
    predictorsnames<-all.vars(model[[3]])
    names(beta_ridge )<- c(predictorsnames )
    W <- param_ridge_sem$W
    

    
    
    
    # Step 4----
    param_ridge_sem <- estimation_sem(data, model, beta_ridge, scale=scale) # lambda estimation for y=lambda*W*y+ X*beta_ridge + epsilon
    
    x_filtered<-param_ridge_sem$x
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
    plotdat_sem <-tibble(logLik=mlogLik_sem, gamma=gamma_int)
    coef_sem <- drop(solve(crossprod(x_filtered) + diag(gamma_sem, (p)), crossprod(x_filtered, y_filtered)))
    y_predicted <- solve((diag(1,n)-lambda*W),(x_filtered%*%coef_sem))
  
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
    # ggtitle( bquote("SEM" ~ beta[ridge] ~ " coefficients paths  vs regularization parameter"))
    
    
    p_logLik<-ggplot(plotdat_sem, aes(y=logLik, x=log(gamma)))+
      geom_line()+
      theme_classic()+
      geom_vline(aes(xintercept = log(gamma_sem), color="red"))+
      guides(color = "none") +
      xlab( expression(paste("log(", gamma,")")))#+
    # ggtitle("SEM ridge logLikelihood vs regularization parameter")
    
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
  
  p_f_test <- parLapply(cl, predictorvarnames, fun=permute_variable_j)
  
  stopCluster(cl)
  
  res <- data.frame(do.call(rbind, p_f_test))
  
  return(res)
}


permute_variable_j <- function(j){
  
  p_f_test <- data.frame(SAR = NA, SEM = NA)
  rownames(p_f_test)<-j
  
  #Compute standard F statistic
  data_j<-data%>%select(-j) #Eliminate variable j from data set
  model_j<-update(model, paste("~ . -",j)) #Eliminate variable j from formula
  
  #SAR
  scale=FALSE
  ridge_sar_j <- rrsar(data_j,model_j,buffer, scale=scale) #Do ridge
  
  
  y_sar_pred_j <- ridge_sar_j$predicted_values #Obtain estimated values  y_0
  
  F_sar_j   <-  (sqrt(crossprod(y_lr-y_sar_pred_j))^2/sqrt(crossprod(y_lr-y_sar_est))^2)-1
  
  
  
  
  #SEM
  scale=FALSE
  ridge_sem_j <-  rrsem(data_j,model_j,buffer, scale=scale) #Do ridge
  
  
  y_sem_pred_j <- ridge_sem_j$predicted_values #Obtain estimated values y_0
  
  
  F_sem_j   <- (sqrt(crossprod(y_lr-y_sem_pred_j))^2/sqrt(crossprod(y_lr-y_sem_est))^2)-1
  
  
  #Obtain B permutations of variable j
  
  data_p <- data
  st_geometry(data_p) <- NULL
  N <- t(replicate(B+B*0.1, sample(as.matrix(data_p%>%select(j)),n, replace = TRUE)))
  out <- N[!(duplicated(N)), ][1:B,]
  out <- split(out, row(out))#transform to list
  x_lr_perm<-x_lr
  
  F_j_sar <- sapply(out, perm_sar, varname=j,  y_sar_pred_j, y_lr, F_sar_j)
  
  F_j_sem <- sapply(out, perm_sem, varname=j,  y_sem_pred_j, y_lr, F_sem_j)

  
  p_f_test[j,"SAR"] <- format(sum(F_j_sar)/B, digits = 4)
  p_f_test[j,"SEM"] <- format(sum(F_j_sem)/B, digits = 4)
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
  
  f<-(sqrt(crossprod(y_lr-y_sar_pred_j))^2/sqrt(crossprod(y_lr-y_sar_perm))^2)-1
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
  f<-(sqrt(crossprod(y_lr-y_sem_pred_j))^2/sqrt(crossprod(y_lr-y_sem_perm))^2)-1
  F_j_sem<- f >= F_sem_j
  
  return(F_j_sem)
  
}
