#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%% CBS (Tang et al., 2022)                                                                      %%%%%
#%%% To faciliate comparison with SIS + OAL and  SIS + GOAL                                       %%%%%
#%%% We used the same simulated data,lambda and gamma values, as in Shortreed and Ertefaie (2017) %%%%%
#%%% Example: Osteosarcoma study; threshold for screening = floor(n/log(n)) (Fan and Lv, 2008)    %%%%%
#%%% This R (version 3.6.2) code  was prepared  by Ismaïla Baldé for comparison purposes          %%%%%
#%%% with SIS+OAL and SIS+GOAL                                                                    %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% Install and load the packages: "devtools", "lqa", "Ball", "glmnet", "grpreg", "MASS" and %%%%%%%%%
#%%%% "cytominer"                                                                              %%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
# install.packages("cytominer")
library(cytominer)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expit = function(x){ 
  pr = ( exp(x) / (1+exp(x)) ) 
  return(pr)
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

### create_weights
### FUNCTION: create weight for n samples, this function is use for 
###           tuning parameter selection
###           
###
### INPUT:     ps, a n*1 vector, estimated propensity score for each sample
###            D , a n*1 vector, treatment =1 if treated; = 0 under control
###
### OUTPUT:    weight, a n*1 vector weight for each sample.


create_weights = function(ps,A){
  weight = ps^(-1)
  weight[A==0] = (1 - ps[A==0])^(-1)
  return(weight)
}


### wAMD_function
### FUNCTION: calculate weighted absolute mean difference for each choice of 
### tuning parameter. This function is use for selecting tuning parameter.
###           
###
### INPUT:     DataM   , a n*(p+2) matrix, contains both covariate, treatment and outcome.
###            varlist , a p*1  vector,names for each covariate.
###            trt.var , name of treatment, in our simulation it is "D".
###            beta    , a p*1 vector, coefficient for each covariate.
###            

### OUTPUT: a list (a) diff_vec , mean difference for each covaiate.
###                (b) wdiff_vec, weighted mean difference for each covaiate.
###                   (diff_vec*|\hat{beta}_j|)  
###                (c) wAMD     , sum of weighted mean difference
###
###
###
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


### CBS
### FUNCTION: Proposed CBS method for estimating average treatment effect.
###           
###
### INPUT:     
###            X:   a n*p  matrix: n*p covariate.
###            D    a n*1  vector: treatment =1 if treated; = 0 under control.
###            Y    a n*1  vector: observed outcome.
###            alpha: double; significant level for confidence interval construction.

### OUTPUT: a list (a) point estimate, the doubly debiased lasso estimator.
###                (b) lower bound, lower bound of the interval.
###                (b) upper bound, upper bound of the interval.
###                (c) Variance, a double; the variance of the dblasso estimator.
###
###
###

# calculate results
f=function(Estim, bVect){
  sd=sd(bVect)
  me=1.96*sd
  CI_n=c(Estim-me,Estim+me)
  bias=xbar-Estim
  LN=CI_n[2]-CI_n[1]
  MSE=(1/length(bVect))*sum((bVect-rep(Estim,length(bVect)))^2)
  return(c(ATE=Estim,Mean=xbar,Bias=bias, SE=sd,MSE=MSE,NCI=CI_n, Length=LN))
}

funCBS <- function(dat,X=X, alpha = 0.05, q=12, n=102){
  #dat=Data.boot
  var.list_Ball = c(paste("X",1:q,sep=""))
  Data=NULL
  Data=dat
  X=as.matrix(X)
  X_Sel=as.matrix(Data[,var.list_Ball])
  
  #set threshold for screening.
  ballcorSel<-rep(NA, q)
  for (j in 1:q){
    # calculate conditional ball covariance for each variable.
    ballcorSel[j]<-Causal.cor(X_Sel[,j],Y,A)
  }
  
  
  weightm = ballcorSel
  betaXY = weightm/max(weightm)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  #var.list_Ball = c(paste("X",1:threshold,sep=""))
  #var.list = c(paste("X",1:p,sep=""))
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # set vector of possible lambda's to try
  lambda_val = n^c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_val) = paste0("lambda",1:9)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*( gamma_convergence_factor - lambda_val + 1 )
  names(gamma_vals) = paste0("gamma",gamma_vals)
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  # centerlize, standardlize
  temp.mean = colMeans(Data[,var.list_Ball])
  Temp.mean = matrix(temp.mean,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list_Ball] = Data[,var.list_Ball] - Temp.mean
  temp.sd = apply(Data[var.list_Ball],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list_Ball),nrow=nrow(Data),byrow=TRUE)
  Data[var.list_Ball] = Data[,var.list_Ball] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  
  
  
  ## Want to save wAMD and propensity score coefficients for
  ## each lambda and gamma value
  
  wAMD_vec = rep(NA, length(lambda_val)*length(gamma_vals))
  names(wAMD_vec) = paste0( rep(names(lambda_val),each = length(gamma_vals)),names(gamma_vals))
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list_Ball),ncol=length(wAMD_vec)))
  names(coeff_XA) = names(wAMD_vec)
  rownames(coeff_XA) = var.list_Ball
  
  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda and gamma value ####################
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  
  
  for( lil in names(lambda_val) )
    for(mim in names(gamma_vals)){
      il = lambda_val[lil]
      ig = gamma_vals[mim]
      fitx = as.matrix(Data[,1:q])
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
      Data[,paste("w",paste0( lil,mim),sep="")] = create_weights(ps=Data[,paste("f.pA",paste0( lil,mim),sep="")],A=Data$A)
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_vec[paste0( lil,mim)] = wAMD_function(DataM=Data,varlist=var.list_Ball,trt.var="A",
                                                 wgt=paste("w",paste0( lil,mim),sep=""),beta=betaXY)$wAMD
      # save ATE estimate for this lambda value
    } # close loop through lambda values
  
  # find the target (gamma,lambda) with smallest wAMD score.
  tt = which.min(wAMD_vec)
  
  #save the estimated propensity score model
  fitted.ps = Data[,paste("f.pA",names(tt),sep="")]
  coeff_XA[var.list_Ball,tt]
  mCBS=ifelse(abs(coeff_XA[var.list_Ball,tt])> 10^(-8),1,0)
  #outcome regression
  {
    # use lasso to fit the treated group.
    X1 = X[A==1,]
    X1=as.matrix(X1)
    Y1 = Y[A==1]
    X0 = X[A==0,]
    X0=as.matrix(X0)
    Y0 = Y[A==0]
    
    lasso<-cv.glmnet(x = X1, y = Y1,
                     type.measure = "mse",
                     family = "binomial",
                     ## K = 10 is the default.
                     nfold = 10,
                     alpha = 1)
    coeflasso1 <- coef.glmnet(lasso,lasso$lambda.min)
    index1 = which(coeflasso1[-1]!=0)
    X1.fit = X1[,index1]
    
    data1 = data.frame(Y1,X1.fit)
    names(data1) = c("Y",paste0("v",1:length(index1)))[1:ncol(data1)]
    
    #ormodel1 <- lm(Y~.,data1)
    ormodel1 <- glm(Y~.,family = "binomial",data1)
    
    data.fit = data.frame(Y,X[,index1])
    names(data.fit) = c("Y",paste0("v",1:length(index1))) [1:ncol(data1)]
    # save predict for the treated group.
    #orfit1 = predict.lm(ormodel1, newdata = data.fit)
    orfit1 = predict.glm(ormodel1, newdata = data.fit, type = "response")
    
    # use lasso to fit the controlled group.
    lasso<-cv.glmnet(x = X0, y = Y0,
                     type.measure = "mse",
                     family = "binomial",
                     ## K = 10 is the default.
                     nfold = 10,
                     alpha = 1)
    coeflasso0 <- coef.glmnet(lasso,lasso$lambda.min)
    index0 = which(coeflasso0[-1]!=0)
    X0.fit = X0[,index0]
    
    data0 = data.frame(Y0,X0.fit)
    names(data0) = c("Y",paste0("v",1:length(index0)))[1:ncol(data0)]
    
    #ormodel0 <- lm(Y~.,data0)
    ormodel0 <- glm(Y~.,family = "binomial",data0)
    
    data.fit = data.frame(Y,X[,index0])
    names(data.fit) = c("Y",paste0("v",1:length(index0))) [1:ncol(data0)]
    #save the estimated data for the controlled group. 
    #orfit0 = predict.lm(ormodel0, newdata = data.fit)
    orfit0 = predict.glm(ormodel0, newdata = data.fit, type = "response")
  }
  # get the double robust estimate.
  result = DR(fitted.ps,orfit1,orfit0,A,Y) 
  
  # get the point estimate and variance of the resulting estimator. 
  result = c(mean(result),var(result))
  c("point estimate" = result[1],
    "lower bound" = result[1]-qnorm(1-alpha/2)*sqrt(result[2]/n),
    "upper bound" = result[1]+qnorm(1-alpha/2)*sqrt(result[2]/n),
    "variance" =  sqrt(result[2]/n))
  
  return(c(result[1],mCBS))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Application to radiomics data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(2015)
setwd("~/Desktop/Ismaila_Recherche_A2022/Radiomic_Paper")
RData <- read.csv("osteosarcoma.csv")
# outcome variable
Y=as.numeric(RData$tumor.volume=='effective group')
summary(as.factor(Y))

# Exposure variable
A=as.numeric(RData$Surgical.staging=="Ⅲ")
summary(as.factor(A))

# Covariates: 1409 Radiomic features
Conf=RData[,10:1418]
char_columns <- sapply(Conf, is.character)                     # Identify character columns
data_chars_as_num <- Conf                                      # Replicate data
data_chars_as_num[ , char_columns] <- as.data.frame(           # Recode characters as numeric
  apply(data_chars_as_num[ , char_columns], 2, as.numeric))
invisible(sapply(data_chars_as_num, class))  
X=data_chars_as_num
X_or=as.matrix(X)
dim(X)

# sample size
n = nrow(X)

# total number of radiomics predictors
p = ncol(X)

# threshold for screening (taken from Fan and Lv (2008))
d.n=floor(n/log(n)) 
threshold = min(d.n,p) 

#set threshold for screening.
ballcor<-rep(NA, p)
for (j in 1:p){
  # calculate conditional ball covariance for each variable.
  ballcor[j]<-Causal.cor(X[,j],Y,A)
}

# list of threshold variables
var.list_Ball = c(paste("X",1:threshold,sep=""))
# list of all p variables
var.list = c(paste("X",1:p,sep=""))

# screening procedure 
ballorder<-order(ballcor, decreasing=TRUE)
# select the top threshold number of variables
ballsisselectedindex<-ballorder[1:threshold]
ballsisselectedindex = ballsisselectedindex[order(ballsisselectedindex)]
weight = ballcor[ballsisselectedindex]

# the data matrix after screening: X_SIS
X_SIS = X[,ballsisselectedindex]
X_SIS=as.data.frame(X_SIS)

# Remove redundant variables
lv_vt=ls(X_SIS)
index_rr=correlation_threshold(lv_vt, X_SIS, cutoff = 0.95, method = "pearson")
Data_CT <- X_SIS[ , ! names(X_SIS) %in% index_rr]                   # Apply %in%-operator

Xq=Data_CT
head(Xq)

# final radiomics variable selected for the osteosarcoma study
var_selected=c("SmallDependenceLowGrayLevelEmphasis.7", "SumAverage.8", "ZonePercentage.9", "ClusterTendency.10", 
               "DependenceEntropy.10", "LargeDependenceLowGrayLevelEmphasis.10", "SmallDependenceLowGrayLevelEmphasis.7",
               "ShortRunLowGrayLevelEmphasis.10", "LongRunLowGrayLevelEmphasis.10", "Variance.12", "Kurtosis.12", "ZoneEntropy.12")
q=ncol(Xq)
Xq=as.matrix(Xq)
Data = as.data.frame(Xq)
var.list_Ball = c(paste("X",1:q,sep=""))
var.list=var.list_Ball
names(Data) = var.list_Ball
Data$A = A
Data$Y = Y

# save data for the bootstrap
Data.boot=Data

head(Data.boot)
dim(Data.boot)
# ATE estimate before the bootstrap
ResIni=funCBS(Data.boot, X=X_or)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%% Boostrap to calculate the SE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=10000
unifnum=rep(NA,n)
Boodata=list()
# save the S monte carlo results (MCResults)
Resvec=matrix(rep(NA,(q+1)*B), ncol = q+1,byrow = TRUE)
# Start the clock!
library(parallel)
ptm <- proc.time()
for (b in 1:B) {
  unifnum = sample(c(1:n),n,replace = T)
  Boodata[[b]]=Data.boot[unifnum,]
  X_CBS=as.matrix(X_or[unifnum,])
  Resvec[b,]=as.numeric(funCBS(Boodata[[b]],X_CBS)) # estimate for iteration b
}

# Stop the clock
proc.time() - ptm
## print tabulated results from the B bootstrap replications of ATE with CBS
AA=Resvec[,1]
# Calculate 10th & 90th percentiles
A4=NA
A4=AA
A4r= quantile(A4, c(0.10, 0.90))  
A4r
A4sub<- A4[A4> A4r[1] &   # Drop rows below/above percentiles
             A4< A4r[2]]

CBSvec=A4sub
xbar=mean(CBSvec)
SIS_CBS=f(ResIni[1],CBSvec)

SIS_CBS=as.data.frame(rbind(round(SIS_CBS,3)))
tabres=c("ATE", "Mean","Bias", "SE", "MSE","95%NCIL", "95%NCIU", "Length")
colnames(SIS_CBS)=tabres
rownames(SIS_CBS)="CBS*"
View(SIS_CBS)
print(SIS_CBS)

# Variable selction for CBS
CBSVS=Resvec[,2:13] 
# print proportion of times covariate selected from the S replications when estimating the ATE with CBS
CBSVS_prop=cbind(colMeans(CBSVS))
rownames(CBSVS_prop)=var_selected
colnames(CBSVS_prop)="Proportion of times covariate selected (%)"
CBSVS_prop=as.data.frame(round(CBSVS_prop,3)*100)
View(CBSVS_prop)
print(CBSVS_prop)
