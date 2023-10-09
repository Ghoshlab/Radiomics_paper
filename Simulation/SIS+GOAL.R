#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%% SIS + GOAL: designed for ultra-high dimension data (p >>n)                                             %%%%%
#%%% R code implementing the algorithm SIS + GOAL                                                           %%%%%
#%%% Note 1: GOAL is the genaralized outcome-adaptve lasso proposed by Baldé et. (2023)                     %%%%%
#%%% Note 2: GOAL genaralized OAL (Shortreed and Ertefaie, 2017)) for high dimension data                   %%%%%
#%%% and   SIS is the sure independence screening proposed  by Tang et al. (20222)                          %%%%%
#%%% Example: Scenario 1; rho=0; n=300; p=1000; threshold for screening= floor(n/log(n)) (Fan and Lv, 2008) %%%%%
#%%% Writen for R version 3.6.2 (Author: Ismaila Baldé)                                                     %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% Install and load the packages: "devtools", "lqa", "Ball", "glmnet", "grpreg" and "MASS" and         %%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# install.packages("devtools") # install "devtool" package first
library(devtools)
# install_github("cran/lqa")  # install "lqa" package from GitHub
library(lqa)
# install.packages("Ball")
library(Ball)
# install.packages("glmnet")
library(glmnet)
# install.packages("grpreg")
library(grpreg)
# install.packages("MASS")
library(MASS) # version 3.3.1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list=ls())    # to remove all objects from the current workspace
ls()             # to check object list
set.seed(2015)   # to create repeatable data sets

# set information for simulating coviariates
mean_x = 0
sig_x = 1

# pairwise correlation between covariates
rho = 0          

# set number of monte carlo (MC) simulation
S=1000

# sample size
n = naug = 300

# total number of predictors
p = 1000

# threshold for screening (taken from Fan and Lv (2008))
d.n=floor(n/log(n)) 
threshold = min(d.n,p)  

# set information for data augmentation
paug=threshold

# note:  pC, pP and pI are number of confounders, pure predictors of outcome and pure predictors of exposure, respectively
pC = pP = pI  = 2

# pS number of spurious covariates
pS = p - (pC+pP+pI)

# list of all p variables
var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""),paste("Xi",1:pI,sep=""),paste("Xs",1:pS,sep=""))

# list of threshold variables
var.list_Ball = c(paste("X",1:threshold,sep=""))

# set strength of relationship between covariates and outcome
beta_v =  c( 0.6, 0.6, 0.6, 0.6, 0, 0, rep(0,p-6) )

# Set strength of relationship between covariates and treatment
alpha_v = c( 1, 1, 0, 0, 1, 1,  rep(0,p-6) )
names(beta_v) = names(alpha_v) = var.list

# set true average treatment effect (taken from Shortreed and Ertefaie (2017))
bA = 0  

# set vector of possible lambda's to try (taken from Shortreed and Ertefaie (2017))
lambda_vec = c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
names(lambda_vec) = as.character(lambda_vec)
# lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
gamma_convergence_factor = 2
# get the gamma value for each value in the lambda vector that corresponds to convergence factor
gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
names(gamma_vals) = names(lambda_vec)

# define some functions for generating data, ATE estimates, and the wAMD,
expit = function(x){
  pr = ( exp(x) / (1+exp(x)) )
  return(pr)
}
ATE_est = function(fY,fw,fA){
  t_ATE = fY*fw
  tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
  return(tt_ATE)
}
create_weights = function(fp,fA,fw){
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  return(fw)
}
wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
  trt = untrt = diff_vec = rep(NA,length(beta))
  names(trt) = names(untrt) = names(diff_vec) = varlist
  for(jj in 1:length(varlist)){
    this.var = paste("w",varlist[jj],sep="")
    DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt]
    trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt])
    untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt])
    diff_vec[jj] = abs( trt[jj] - untrt[jj] )
  }
  wdiff_vec = diff_vec * abs(beta)
  wAMD = c( sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret)
}

#  SIS procedure (taken from Tang et al. (20222))
# Causal.cor
# FUNCTION:  calculate conditional ball covariance 
#           between x and y given z.
#
#INPUT:     x, a n*p matrix, samples of variable X.  
#            y, a n*q matrix, samples of variable y.
#            z, a n*1 binary vector, samples of   z . 
#           distance, if distance = TRUE, the elements
#           of x and y are considered as distance matrices.
#
# OUTPUT:    double; conditional ball covariance 
#              between x and y given z.

Causal.cor <- function(x, y, z, distance = FALSE) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  index0 <- which(z == 0)
  index1 <- which(z == 1)
  alpha = length(index0)/length(z)
  if (distance == TRUE) {
    x0 <- x[index0, index0]
    y0 <- y[index0, index0]
    x1 <- x[index1, index1]
    y1 <- y[index1, index1]
    #see definition 4
    Causal.cor <- alpha*bcov(x0, x0, distance = TRUE)^2 + (1-alpha)*bcov(x1, y1, distance = TRUE)^2
  } else {
    x0 <- x[index0, ]
    y0 <- y[index0, ]
    x1 <- x[index1, ]
    y1 <- y[index1, ]
    #see definition 4
    Causal.cor <- alpha*bcov(x0, y0)^2 + (1-alpha)*bcov(x1, y1)^2
  }
  return(Causal.cor)
}

#########################################################################################
##################### modified lqa.update2 #############################################
#########################################################################################

lqa.upd <-
  function (x, y, family = NULL, penalty = NULL, intercept = TRUE, weights = rep (1, nobs), control = lqa.control (), initial.beta, mustart, eta.new, gamma1 = 1, ...)
  {
    gamma <- gamma1
    
    if (is.null (family))
      stop ("lqa.update: family not specified")
    
    if (is.null (penalty))
      stop ("lqa.update: penalty not specified")
    
    if (!is.null (dim (y)))
      stop ("lqa.update: y must be a vector")
    
    x <- as.matrix (x)
    converged <- FALSE
    n.iter <- control$max.steps
    eps <- control$conv.eps
    c1 <- control$c1
    
    stop.at <- n.iter
    pcol <- ncol (x)
    nobs <- nrow (x)
    converged <- FALSE
    beta.mat <- matrix (0, nrow = n.iter, ncol = pcol)    # to store the coefficient updates
    
    if (missing (initial.beta))
      initial.beta <- rep (0.01, pcol)
    else
      eta.new <- drop (x %*% initial.beta)
    
    if (missing (mustart))
    {
      etastart <- drop (x %*% initial.beta)
      eval (family$initialize)
    }
    
    if (missing (eta.new))
      eta.new <- family$linkfun (mustart)    # predictor
    
    
    
    for (i in 1 : n.iter)
    {
      beta.mat[i,] <- initial.beta  
      mu.new <- family$linkinv (eta.new)      # fitted values
      d.new <- family$mu.eta (eta.new)        # derivative of response function
      v.new <- family$variance (mu.new)       # variance function of the response
      weights <- as.vector(d.new / sqrt (v.new)) # decomposed elements (^0.5) of weight matrix W, see GLM notation
      weights <- c(weights[1:naug],rep(1,paug))
      x.star <- weights * x  
      y.tilde.star <- as.vector(weights * (eta.new  + (y - mu.new) / d.new))
      y.tilde.star <- c( y.tilde.star[1:naug],rep(0,paug))  
      A.lambda <- get.Amat (initial.beta = initial.beta, penalty = penalty, intercept = intercept, c1 = c1, x = x)
      p.imat.new <- crossprod (x.star) + A.lambda       # penalized information matrix
      
      chol.pimat.new <- chol (p.imat.new)               # applying cholesky decomposition for matrix inversion
      inv.pimat.new <- chol2inv (chol.pimat.new)        # inverted penalized information matrix
      beta.new <- gamma * drop (inv.pimat.new %*% t (x.star) %*% y.tilde.star) + (1 - gamma) * beta.mat[i,]  # computes the next iterate of the beta vector
      
      
      if ((sum (abs (beta.new - initial.beta)) / sum (abs (initial.beta)) <= eps))    # check convergence condition
      {
        converged <- TRUE
        stop.at <- i
        if (i < n.iter)
          break
      }
      else
      {
        initial.beta <- beta.new    # update beta vector
        eta.new <- drop (x %*% beta.new)      
      }
    }
    
    
    Hatmat <- x.star %*% inv.pimat.new %*% t (x.star)
    tr.H <- sum (diag (Hatmat))
    dev.m <- sum (family$dev.resids (y, mu.new, weights))
    
    aic.vec <- dev.m + 2 * tr.H
    bic.vec <- dev.m + log (nobs) * tr.H
    
    if (!converged & (stop.at == n.iter))
      cat ("lqa.update with ", penalty$penalty, ": convergence warning! (lambda = ", penalty$lambda, ")\n")
    
    
    fit <- list (coefficients = beta.new, beta.mat = beta.mat[1 : stop.at,], tr.H = tr.H, fitted.values = mu.new, family = family, Amat = A.lambda, converged = converged, stop.at = stop.at, m.stop = stop.at, linear.predictors = eta.new, weights = weights^2, p.imat = p.imat.new, inv.pimat = inv.pimat.new, x.star = x.star, v.new = v.new)
  }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%       modified lqa.default         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lqa.def <- function (x, y, family = binomial (), penalty = NULL, weights = rep (1, nobs),
                     start = NULL, etastart = NULL, mustart = NULL, offset = rep (0, nobs),
                     control = lqa.control (), intercept = TRUE, standardize = TRUE,method="lqa.update2")
{
  call <- match.call ()
  
  ### Check for exponential family and link function:
  ### -----------------------------------------------
  
  if (is.character (family))
    family <- get (family, mode = "function", envir = parent.frame ())
  
  if (is.function (family))
    family <- family ()
  
  if (is.null (family$family))
  {
    print (family)
    stop ("'family' not recognized")
  }
  
  
  ### Check for quadratic penalty:
  ### ----------------------------
  
  if (! (method == "nng.update"))
  {
    
    if (is.null (penalty))
      stop ("penalty not specified \n")
    
    if (is.character (penalty))
      penalty <- get (penalty, mode = "function", envir = parent.frame ())
    
    if (is.function (penalty))
      penalty <- penalty ()
    
    if (is.null (penalty$penalty))
    {
      print (penalty)
      stop ("'penalty' not recognized")
    }
  }
  
  
  ### Check for existence of 'method':
  ### --------------------------------
  
  if (is.null (method))
    stop ("method not specified")
  
  
  ### Check for column of ones in x if intercept == TRUE:
  ### ---------------------------------------------------
  
  if (intercept & (var (x[,1]) > control$var.eps))
    x <- cbind (1, x)
  
  
  ### Standardization:
  ### ----------------
  
  x <- as.matrix (x)
  xnames <- dimnames (x)[[2L]]
  
  ynames <- if (is.matrix (y))
    rownames(y)
  else
    names(y)
  
  nobs <- nrow (x)  
  nvars <- ncol (x)    # number of coefficients in the predictor (including an intercept, if present)
  ones <- rep (1, nobs)
  mean.x <- drop (ones %*% x) / nobs       # computes the vector of means
  
  if (intercept)    # if an intercept is included in the model its corresponding element of mean.x is set to zero
    mean.x[1] <- 0  # (such that x[,1] is not getting centered (and also not standardized later on ...))
  
  x.std <- scale (x, mean.x, FALSE)   # centers the regressor matrix
  
  norm.x <- if (standardize)
  {
    norm.x <- sqrt (drop (ones %*% (x.std^2)))   # computes the euclidean norm of the regressors
    nosignal <- apply (x, 2, var) < control$var.eps
    if (any (nosignal))    # modify norm.x for variables with too small a variance (e.g. the intercept)
      norm.x[nosignal] <- 1
    
    norm.x
  }
  else
    rep (1, nvars)
  
  x.std <- scale (x.std, FALSE, norm.x)  # standardizes the centered regressor matrix
  
  
  ### Call and get the (estimation) method:
  ### -------------------------------------
  
  
  
  
  fit=lqa.upd(x = x.std, y = y, family = family, penalty = penalty, intercept = intercept,
              control = control)
  #fit=lqa.update2(x = x.std, y = y, family = family, penalty = penalty, intercept = intercept,
  #  control = control)
  
  #fit <- do.call (method, list (x = x.std, y = y, family = family, penalty = penalty, intercept = intercept,
  #                             control = control))
  
  
  ### Back-Transformation of estimated coefficients:
  ### ----------------------------------------------
  
  coef <- fit$coefficients
  
  if (intercept)
  {
    coef[1] <- coef[1] - sum (mean.x[-1] * coef[-1] / norm.x[-1])
    coef[-1] <- coef[-1] / norm.x[-1]
  }
  else
    coef <- coef / norm.x
  
  
  ### Computation of some important statistics:
  ### -----------------------------------------
  
  # Remark: The predictors are identical no matter whether we use standardized or unstandardized values, hence all statistics
  #         based on the predictor eta are also equal
  
  eta <- drop (x %*% coef)
  mu <- family$linkinv (eta)
  mu.eta.val <- family$mu.eta (eta)
  wt <- sqrt ((weights * mu.eta.val^2) / family$variance (mu))
  dev <- sum (family$dev.resids (y, mu, weights))
  wtdmu <- sum (weights * y) / sum (weights)
  nulldev <- sum (family$dev.resids (y, wtdmu, weights))
  n.ok <- nobs - sum (weights == 0)
  nulldf <- n.ok - as.integer (intercept)
  residuals <- (y - mu) / mu.eta.val
  
  xnames <- colnames (x)
  ynames <- names (y)
  names (residuals) <- names (mu) <- names (eta) <- names (weights) <- names (wt) <- ynames
  names (coef) <- xnames
  
  Amat <- fit$Amat
  Amat <- t (norm.x * t (norm.x * Amat))   # must be (quadratically) transformed in order to cope with the transformed parameter space
  
  if (is.null (fit$tr.H))
    stop ("quadpen.fit: Element 'tr.H' has not been returned from 'method'")
  
  model.aic <- dev + 2 * fit$tr.H
  model.bic <- dev + log (nobs) * fit$tr.H
  resdf <- n.ok - fit$tr.H
  
  dispersion <- ifelse (!((family$family == "binomial") | (family$family == "poisson")), sum ((y - mu)^2 / family$variance (mu)) / (nobs - fit$tr.H), 1)
  
  fit <- list (coefficients = coef, residuals = residuals, fitted.values = mu, family = family, penalty = penalty,
               linear.predictors = eta, deviance = dev, aic = model.aic, bic = model.bic, null.deviance = nulldev, n.iter = fit$stop.at,
               best.iter = fit$m.stop, weights = wt, prior.weights = weights, df.null = nulldf, df.residual = resdf, converged = fit$converged, mean.x = mean.x,
               norm.x = norm.x, Amat = Amat, method = method, rank = fit$tr.H, x = x, y = y, fit.obj = fit, call = call, dispersion = dispersion)
  
  class (fit) <- c ("lqa", "glm", "lm")
  fit
}

# GOAL of  Baldé et al. (2023) with SIS procedure

GOAL<- function(X,A,Y){
  
  #set threshold for screening.
  ballcor<-rep(NA, p)
  for (j in 1:p){
    # calculate conditional ball covariance for each variable.
    ballcor[j]<-Causal.cor(X[,j],Y,A)
  }
  
  # screening procedure 
  ballorder<-order(ballcor, decreasing=TRUE)
  # select the top threshold number of variables
  ballsisselectedindex<-ballorder[1:threshold]
  ballsisselectedindex = ballsisselectedindex[order(ballsisselectedindex)]
  weight = ballcor[ballsisselectedindex]
  
  # the data matrix after screening
  Data=NULL
  Data = X[,ballsisselectedindex]
  Data = as.data.frame(Data)
  names(Data) = var.list_Ball
  Data$A = A
  Data$Y = Y
  
  # normlize coviarates to have mean 0 and standard deviation 1
  temp.mean = colMeans(Data[,var.list_Ball])
  Temp.mean = matrix(temp.mean,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list_Ball] = Data[,var.list_Ball] - Temp.mean
  temp.sd = apply(Data[var.list_Ball],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[var.list_Ball] = Data[,var.list_Ball] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  
  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list_Ball,collapse="+")))
  lm.Y = lm(y.form,data=Data)
  betaXY = coef(lm.Y)[var.list_Ball] 
  summary(lm.Y)
  betaXY
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%                 create Augmented data      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # set the possible lambda2 value (taken from Zou and Hastie (2005))
  S_lam=c(0,10^c(-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1))
  WM_N=WM_P=M_N=M_P=S_wamd=rep(NA,length(S_lam))
  M_mat=matrix(NA,length(S_lam),threshold)
  for (m in 1:length(S_lam)) {
    # augmented A and X 
    lambda2=S_lam[m]
    q=threshold
    I=diag(1,q,q)
    Iq=sqrt(lambda2)*I
    Anq=c(Data$A,rep(0,q))
    Xnq=matrix(NA,n+q,q)
    X_th=Data[,var.list_Ball]
    for (j in 1:q){
      Xnq[,j]=c(X_th[,j],Iq[,j])
    }
    newData=as.data.frame(Xnq)
    names(newData)=var.list_Ball
    newData$A=Anq
    n.q=n+q
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%% run GOAL based on PIRLS                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # want to save ATE, wAMD and propensity score coefficients for each lambda value
    ATE = wAMD_vec = rep(NA, length(lambda_vec))
    names(ATE) = names(wAMD_vec) = names(lambda_vec)
    coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list_Ball),ncol=length(lambda_vec)))
    names(coeff_XA) = names(lambda_vec)
    rownames(coeff_XA) = var.list_Ball
    
    # weight model with all possible covariates included, this is passed into lasso function
    for( lil in names(lambda_vec) ){
      il = lambda_vec[lil]
      ig = gamma_vals[lil]
      
      # create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
      oal_pen = adaptive.lasso(lambda=n.q^(il),al.weights = abs(betaXY)^(-ig) )
      # run outcome-adaptive lasso model with appropriate penalty
      X=as.matrix(newData[var.list_Ball]);y=as.vector(newData$A);
      lq22=lqa.def(x=X, y=y, family = binomial,penalty = oal_pen)
      # generate propensity score for ATE
      Data[,paste("f.pA",lil,sep="")]=expit(as.matrix(cbind(rep(1,n),Data[var.list_Ball]))%*%as.matrix((1+lambda2)*coef(lq22)))
      # create inverse probability of treatment weights for ATE
      Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
      #save propensity score coef
      coeff_XA[var.list_Ball,lil] = (1+lambda2)*coef(lq22)[var.list_Ball]
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list_Ball,trt.var="A",
                                    wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
      
      # save ATE estimate for this lambda value
      ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)
    } # close loop through lambda values
    
    # print out wAMD for all the lambda values tried
    wAMD_vec
    # find the lambda value that creates the smallest wAMD
    tt = which.min(wAMD_vec)
    # print out ATE corresponding to smallest wAMD value
    ATE[tt]
    # print out the coefficients for the propensity score that corresponds with smalles wAMD value
    GOAL.PIRLS=coeff_XA[,tt]
    GOAL.PIRLS.ate=ATE[tt][[1]]
    
    M_P[m]=GOAL.PIRLS.ate
    WM_P[m]=wAMD_vec[tt][[1]]
    M_mat[m,]=ifelse(abs(coeff_XA[,tt])> 10^(-8),1,0)
    
  } # end for for lambda2
  # final ATE and coef
  
  ptt= which.min(WM_P)
  GOAL=M_P[ptt]
  mGOAL=ifelse(abs(M_mat[ptt,])> 10^(-8),1,0)
  
  return(c(GOAL,mGOAL))
}

# save the S monte carlo results (MCResults)
MCResults=matrix(rep(NA,(threshold+1)*S), ncol = threshold+1,byrow = TRUE)

# Start the clock!
ptm <- proc.time()

for (s in 1:S) {
  
  Data=NULL
  ### simulate data
  Sigma_x = matrix(rho*sig_x^2,nrow=length(var.list),ncol=length(var.list))
  diag(Sigma_x) = sig_x^2
  Mean_x = rep(mean_x,length(var.list))
  Data = as.data.frame(mvrnorm(n = n,mu=Mean_x,Sigma = Sigma_x,empirical = FALSE))
  names(Data) = var.list
  gA_x = rowSums(Data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))
  pA = expit( gA_x )
  Data$A = as.numeric( runif(n=length(pA)) < pA) # simulate A
  gY_xA = rowSums(Data[,var.list]*matrix(beta_v,nrow=n,ncol=length(var.list),byrow=TRUE))  
  Data$Y = gY_xA + rnorm(n=n,sd=sig_x)
  Data$Y = Data$Y + Data$A*bA
  Data.OAL=Data
  X=as.matrix(Data[,var.list])
  A=Data$A
  Y=Data$Y
  
  # save the ATE and the coefficients estimated 
  MCResults[s,]=GOAL(X,A,Y)
  
  # print the current iteration
  print(s)
  
}
# Stop the clock
proc.time() - ptm

# mean estimate from the S replications of ATE with SIS+ GOAL.
mean_ate=mean(MCResults[,1])
# bias (mean_ate-bA)  from the S replications of ATE with SIS+ GOAL.
Bias=mean_ate-bA

# standard error (SE)  from the S replications of ATE with SIS+ GOAL.
SE=sd(MCResults[,1])

# mean squared error (MSE) from the S replications of ATE with SIS+ GOAL.
MSE=(1/S)*sum((MCResults[,1]-rep(bA,S))^2);

# print tabulated results from the S replications of ATE with SIS+ GOAL
SIS_GOAL=c(Bias,SE,MSE)
SIS_GOAL=rbind(SIS_GOAL)
colnames(SIS_GOAL)=c("Bias","SE","MSE") 
SIS_GOAL=as.data.frame(round(SIS_GOAL,2))
rownames(SIS_GOAL)="SIS + GOAL"
View(SIS_GOAL)
print(SIS_GOAL)

# print proportion of times covariate selected from the S replications when estimating the ATE with SIS+ OAL
coef_mat=MCResults[,2:(threshold+1)]
coef_prop=cbind(colMeans(coef_mat))
rownames(coef_prop)=var.list_Ball 
colnames(coef_prop)="Proportion of times covariate selected"
coef_prop=as.data.frame(round(coef_prop,3))
View(coef_prop)
print(coef_prop)

# plot the proportion of times covariate selected
plot(1:threshold,colMeans(MCResults[,2:(threshold+1)]),ylim=c(0,1), type="l",lty=2, col="red", axes=TRUE, ann=FALSE,lwd=3)
# label the x and y axes 
title(xlab="Covariate index")
title(ylab="Proportion of times covariate selected ")
legend("topright",
       "SIS + GOAL",
       col="red", 
       lwd=3, 
       lty=3 
)
