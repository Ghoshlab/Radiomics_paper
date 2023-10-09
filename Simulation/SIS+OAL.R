#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%% SIS + OAL: designed for ultra-high dimension data (p >>n)                                              %%%%%
#%%% R code implementing the algorithm SIS + OAL                                                            %%%%%
#%%% Note: OAL is the outcome-adaptve lasso proposed by Shortreed and Ertefaie (2017)                       %%%%%
#%%% and   SIS is the sure independence screening proposed  by Tang et al. (20222)                          %%%%%
#%%% Example: Scenario 1; rho=0; n=300; p=1000; threshold for screening= floor(n/log(n)) (Fan and Lv, 2008) %%%%%
#%%% Writen for R version 3.6.2 (Author: Ismaila Baldé)                                                     %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% Install and load the packages: "devtools", "lqa", "Ball", "glmnet", "grpreg" and "MASS" and         %%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
n = 300

# total number of predictors
p = 1000

# threshold for screening (taken from Fan and Lv (2008))
d.n=floor(n/log(n)) 
threshold = min(d.n,p) 

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

# OAL of  Shortreed and Ertefaie (2017) with SIS procedure

OAL<- function(X,A,Y){
  
  # set threshold for screening.
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
  
  # want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list_Ball),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list_Ball
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%  Run outcome adaptive lasso for each lambda value                             %%%%%%%%                  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(var.list_Ball,collapse="+")))
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    
    # create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
    oal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig) )
    # run outcome-adaptive lasso model with appropriate penalty
    X=as.matrix(Data[var.list_Ball]);y=as.vector(Data$A);
    logit_oal = lqa.default( x=X, y=y, penalty=oal_pen, family=binomial(logit))
    # generate propensity score for ATE
    Data[,paste("f.pA",lil,sep="")]=expit(as.matrix(cbind(rep(1,n),Data[var.list_Ball]))%*%as.matrix(coef(logit_oal)))
    # save propensity score coefficients
    coeff_XA[var.list_Ball,lil] = coef(logit_oal)[var.list_Ball]
    # create inverse probability of treatment weights
    Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
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
  OAL=ATE[tt][[1]]
  # print out the coefficients for the propensity score that corresponds with smalles wAMD value 
  coeff_XA[,tt]
  mOAL=ifelse(abs(coeff_XA[,tt])> 10^(-8),1,0)
  
  return(c(OAL,mOAL))
}

# save the S monte carlo results (MCResults)
MCResults=matrix(rep(NA,(threshold+1)*S), ncol = threshold+1,byrow = TRUE)

# Start the clock!
ptm <- proc.time()

for (s in 1:S) {
  
  Data=NULL
  # simulate data
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
  MCResults[s,]=OAL(X,A,Y)
  
  # print the current iteration
  print(s)
}

# Stop the clock
proc.time() - ptm

# mean estimate from the S replications of ATE with SIS+ OAL.
mean_ate=mean(MCResults[,1])
# bias (mean_ate-bA)  from the S replications of ATE with SIS+ OAL.
Bias=mean_ate-bA

# standard error (SE)  from the S replications of ATE with SIS+ OAL.
SE=sd(MCResults[,1])

# mean squared error (MSE) from the S replications of ATE with SIS+ OAL.
MSE=(1/S)*sum((MCResults[,1]-rep(bA,S))^2);

# print tabulated results from the S replications of ATE with SIS+ OAL
SIS_OAL=c(Bias,SE,MSE)
SIS_OAL=rbind(SIS_OAL)
colnames(SIS_OAL)=c("Bias","SE","MSE") 
SIS_OAL=as.data.frame(round(SIS_OAL,2))
rownames(SIS_OAL)="SIS + OAL"
View(SIS_OAL)
print(SIS_OAL)

# print proportion of times covariate selected from the S replications when estimating the ATE with SIS+ OAL
coef_mat=MCResults[,2:(threshold+1)]
coef_prop=cbind(colMeans(coef_mat))
rownames(coef_prop)=var.list_Ball 
colnames(coef_prop)="Proportion of times covariate selected"
coef_prop=as.data.frame(round(coef_prop,3))
View(coef_prop)
print(coef_prop)

# plot the proportion of times covariate selected
plot(1:threshold,colMeans(MCResults[,2:(threshold+1)]),ylim=c(0,1), type="l",lty=2, col="blue", axes=TRUE, ann=FALSE,lwd=3)
# label the x and y axes 
title(xlab="Covariate index")
title(ylab="Proportion of times covariate selected ")
legend("topright",
       "SIS + OAL",
       col="blue", 
       lwd=3, 
       lty=2 
)