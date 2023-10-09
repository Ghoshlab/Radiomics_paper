#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% CBS (Tang et al., 2022)                                                                                     %%
#% To faciliate comparison with SIS + OAL and  SIS + GOAL                                                      %%
#% We used the same simulated data,lambda and gamma values, as in Shortreed and Ertefaie (2017)                %%   
#% Example: Scenario 1; rho=0; n=300; p=1000; threshold for screening= floor(n/log(n)) (Fan and Lv, 2008)      %%
#% This R (version 3.6.2) code was prepared by Ismaïla Baldé for comparison purposes with SIS+OAL and SIS+GOAL %%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
lambda_vec = n^c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
names(lambda_vec) = paste0("lambda",1:9)
# lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
gamma_convergence_factor = 2
# get the gamma value for each value in the lambda vector that corresponds to convergence factor
gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
names(gamma_vals) = paste0("gamma",gamma_vals)


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

### Causal.cor
### FUNCTION:  calculate conditional ball covariance 
###            between x and y given z.
###
### INPUT:     x, a n*p matrix, samples of variable X.  
###            y, a n*q matrix, samples of variable y.
###            z, a n*1 binary vector, samples of   z . 
###            distance, if distance = TRUE, the elements
###            of x and y are considered as distance matrices.
###
### OUTPUT:    double; conditional ball covariance 
###               between x and y given z.


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

### DR double robust estimator
### FUNCTION: calculate DR estimator for each sample(which means averaged), for future analysis.
###           
###
### INPUT:     ps   a n*1  vector: estimated propensity score for each samples.
###            or1  a n*1  vector: estimated outcome regression: \hat{b}_1(X_i) for each sample.
###            or0  a n*1  vector: estimated outcome regression: \hat{b}_0(X_i) for each sample.
###            D    a n*1  vector: treatment =1 if treated; = 0 under control.
###            Y    a n*1  vector: observed outcome.
###
### Output:    n*1 DR estimates   
###
###

DR =  function(ps,or1,or0,A,Y){
  (A*Y-(A-ps)*or1)/ps -  (   (1-A)*Y +(A-ps)*or0)/(1-ps)
}


CBS<- function(X,A,Y,alpha=0.05){
  
  p = ncol(X)
  n = nrow(X)
  d.n=floor(n/log(n))
  threshold = min(d.n,p) 
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
  Data = X[,ballsisselectedindex]
  Data = as.data.frame(Data)
  names(Data) = var.list_Ball
  Data$A = A
  Data$Y = Y
  
  
  # centerlize, standardlize
  temp.mean = colMeans(Data[,var.list_Ball])
  Temp.mean = matrix(temp.mean,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list_Ball] = Data[,var.list_Ball] - Temp.mean
  temp.sd = apply(Data[var.list_Ball],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[var.list_Ball] = Data[,var.list_Ball] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%                    CBS                                                  %%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # weight for each variable for refined selection
  betaXY = weight
  betaXY = weight/max(weight)
  
  ## Want to save wAMD and propensity score coefficients for
  ## each lambda and gamma value
  
  wAMD_vec = rep(NA, length(lambda_vec)*length(gamma_vals))
  names(wAMD_vec) = paste0( rep(names(lambda_vec),each = length(gamma_vals)),names(gamma_vals))
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list_Ball),ncol=length(wAMD_vec)))
  names(coeff_XA) = names(wAMD_vec)
  rownames(coeff_XA) = var.list_Ball
  
  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda and gamma value ####################
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  
  
  for( lil in names(lambda_vec) )
    for(mim in names(gamma_vals)){
      il = lambda_vec[lil]
      ig = gamma_vals[mim]
      fitx = as.matrix(Data[,1:threshold])
      alasso <- glmnet(x = fitx, y = Data$A,
                       type.measure = "class",
                       family = "binomial",
                       alpha = 1,
                       penalty.factor = c(abs(betaXY)^(-ig)),
                       lambda = il)
      
      # calculate propensity score 
      prob = predict(alasso,newx = fitx)
      prob = exp(prob)/(1+exp(prob))
      Data[,paste("f.pA",paste0( lil,mim),sep="")] = prob
      # save propensity score coefficients
      coeff_XA[var.list_Ball,paste0( lil,mim)] = coef.glmnet(alasso,s = alasso$lambda.min)[var.list_Ball,]
      # create inverse probability of treatment weights
      Data[,paste("w",paste0( lil,mim),sep="")] = create_weights(fp=Data[,paste("f.pA",paste0( lil,mim),sep="")],fA=Data$A)
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_vec[paste0( lil,mim)] = wAMD_function(DataM=Data,varlist=var.list_Ball,trt.var="A",
                                                 wgt=paste("w",paste0( lil,mim),sep=""),beta=betaXY)$wAMD
      # save ATE estimate for this lambda value
    } # close loop through lambda values
  
  # find the target (gamma,lambda) with smallest wAMD score.
  tt = which.min(wAMD_vec)
  #save the estimated propensity score model
  fitted.ps = Data[,paste("f.pA",names(tt),sep="")]
  #outcome regression
  {
    # use lasso to fit the treated group.
    X1 = X[A==1,]
    Y1 = Y[A==1]
    X0 = X[A==0,]
    Y0 = Y[A==0]
    
    lasso<-cv.glmnet(x = X1, y = Y1,
                     type.measure = "mse",
                     ## K = 10 is the default.
                     nfold = 10,
                     alpha = 1)
    coeflasso1 <- coef.glmnet(lasso,lasso$lambda.min)
    index1 = which(coeflasso1[-1]!=0)
    X1.fit = X1[,index1]
    
    data1 = data.frame(Y1,X1.fit)
    names(data1) = c("Y",paste0("v",1:length(index1)))[1:ncol(data1)]
    
    ormodel1 <- lm(Y~.,data1)
    
    data.fit = data.frame(Y,X[,index1])
    names(data.fit) = c("Y",paste0("v",1:length(index1))) [1:ncol(data1)]
    # save predict for the treated group.
    orfit1 = predict.lm(ormodel1, newdata = data.fit)
    
    # use lasso to fit the controlled group.
    lasso<-cv.glmnet(x = X0, y = Y0,
                     type.measure = "mse",
                     ## K = 10 is the default.
                     nfold = 10,
                     alpha = 1)
    coeflasso0 <- coef.glmnet(lasso,lasso$lambda.min)
    index0 = which(coeflasso0[-1]!=0)
    X0.fit = X0[,index0]
    
    data0 = data.frame(Y0,X0.fit)
    names(data0) = c("Y",paste0("v",1:length(index0)))[1:ncol(data0)]
    
    ormodel0 <- lm(Y~.,data0)
    data.fit = data.frame(Y,X[,index0])
    names(data.fit) = c("Y",paste0("v",1:length(index0))) [1:ncol(data0)]
    #save the estimated data for the controlled group. 
    orfit0 = predict.lm(ormodel0, newdata = data.fit)
  }
  # get the double robust estimate.
  result = DR(fitted.ps,orfit1,orfit0,A,Y) 
  
  # get the point estimate and variance of the resulting estimator. 
  result = c(mean(result),var(result))
  c("point estimate" = result[1],
    "lower bound" = result[1]-qnorm(1-alpha/2)*sqrt(result[2]/n),
    "upper bound" = result[1]+qnorm(1-alpha/2)*sqrt(result[2]/n),
    "variance" =  sqrt(result[2]/n))
  
  
  CBS=result[1]
  mCBS=ifelse(abs(coeff_XA[var.list_Ball,tt])> 10^(-8),1,0)
  
  return(c(CBS,mCBS))  
  
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
  MCResults[s,]=CBS(X,A,Y)
  
  # print the current iteration
  print(s)
  
}

# Stop the clock
proc.time() - ptm

# mean estimate from the S replications of ATE with CBS
mean_ate=mean(MCResults[,1])
# bias (mean_ate-bA)  from the S replications of ATE with CBS
Bias=mean_ate-bA

# standard error (SE)  from the S replications of ATE with CBS
SE=sd(MCResults[,1])

# mean squared error (MSE) from the S replications of ATE with CBS
MSE=(1/S)*sum((MCResults[,1]-rep(bA,S))^2);

# print tabulated results from the S replications of ATE with CBS
CBS=c(Bias,SE,MSE)
CBS=rbind(CBS)
colnames(CBS)=c("Bias","SE","MSE") 
CBS=as.data.frame(round(CBS,2))
rownames(CBS)="CBS"
View(CBS)
print(CBS)

# print proportion of times covariate selected from the S replications when estimating the ATE with SIS+ OAL
coef_mat=MCResults[,2:(threshold+1)]
coef_prop=cbind(colMeans(coef_mat))
rownames(coef_prop)=var.list_Ball 
colnames(coef_prop)="Proportion of times covariate selected"
coef_prop=as.data.frame(round(coef_prop,3))
View(coef_prop)
print(coef_prop)
# plot the proportion of times covariate selected
plot(1:threshold,colMeans(MCResults[,2:(threshold+1)]),ylim=c(0,1), type="l",lty=1, col="black", axes=TRUE, ann=FALSE,lwd=3)
# label the x and y axes 
title(xlab="Covariate index")
title(ylab="Proportion of times covariate selected ")
legend("topright",
       "CBS",
       col="black", 
       lwd=3, 
       lty=1 
)