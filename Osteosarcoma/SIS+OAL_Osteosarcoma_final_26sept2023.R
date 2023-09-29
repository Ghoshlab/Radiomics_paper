#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% SIS + OAL: designed for ultra-high dimension data (p >>n)                                 %%%%%%
#%%%% R code implementing the algorithm SIS + OAL                                               %%%%%%
#%%%% Note: OAL is the outcome-adaptve lasso proposed by Shortreed and Ertefaie (2017)          %%%%%%
#%%%% and   SIS is the sure independence screening proposed  by Tang et al. (20222)             %%%%%%
#%%%% Example: Osteosarcoma study; threshold for screening = floor(n/log(n)) (Fan and Lv, 2008) %%%%%%
#%%%% Writen for R version 3.6.2 (Author: Ismaila Baldé)                                        %%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% Install and load the packages: "devtools", "lqa", "Ball", "glmnet", "grpreg", "MASS" and %%%%%%%
#%%%% "cytominer"                                                                              %%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list=ls())    # to remove all objects from the current workspace
ls()             # to check object list
set.seed(2015)   # to create repeatable data sets

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

# OAL of  Shortreed and Ertefaie (2017) with SIS procedure

funOAL <- function(dat,q=12,n=102){
  
  # list of q variables
  var.list_Ball = c(paste("X",1:q,sep=""))
  var.list=var.list_Ball
  Data=NULL
  Data = as.data.frame(as.matrix(dat[var.list]))
  Data$A=dat$A
  Data$Y=dat$Y
  
  #set vector of possible lambda's to try
  lambda_vec = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_vec) = as.character(lambda_vec)
  
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*(gamma_convergence_factor - lambda_vec + 1)
  names(gamma_vals) = names(lambda_vec)
  
  # Normlize coviarates to have mean 0 and standard deviation 1
  temp.mean = colMeans(Data[,var.list])
  Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list] = Data[,var.list] - Temp.mean
  temp.sd = apply(Data[var.list],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[var.list] = Data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  
  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
  glm.Y = glm(y.form,data=Data,family = "binomial" )
  betaXY = coef(glm.Y)[var.list] 
  summary(glm.Y)
  betaXY
  
  # want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%  Run outcome adaptive lasso for each lambda value                 %%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    
    # create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
    oal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig) )
    # run outcome-adaptive lasso model with appropriate penalty
    X=as.matrix(Data[var.list]);y=as.vector(Data$A);
    logit_oal = lqa.default( x=X, y=y, penalty=oal_pen, family=binomial(logit))
    # generate propensity score for ATE
    Data[,paste("f.pA",lil,sep="")]=expit(as.matrix(cbind(rep(1,n),Data[var.list]))%*%as.matrix(coef(logit_oal)))
    # save propensity score coefficients
    coeff_XA[var.list,lil] = coef(logit_oal)[var.list]
    # create inverse probability of treatment weights
    Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
    # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
    wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
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


# ATE estimate before the bootstrap
ResIni=funOAL(Data)

# save data for the bootstrap
Data.boot=Data
dim(Data.boot)
head(Data.boot)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%   Bootstrap to calculate point and confidence intervals estimates %%%%%%%%%%%%%%%%%%%%%% 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Number of boostrap iterations
B=10000
unifnum=rep(NA,n)
Boodata=list()
# save the B boostrap results (Resvec)
Resvec=matrix(rep(NA,(q+1)*B), ncol = q+1,byrow = TRUE)

# Start the clock!
library(parallel)
ptm <- proc.time()
for (b in 1:B) {
  unifnum = sample(c(1:n),n,replace = T)
  Boodata[[b]]=Data.boot[unifnum,]
  Resvec[b,]=invisible(funOAL(Boodata[[b]])) # estimate for iteration b
}

# Stop the clock
proc.time() - ptm

OALvec=Resvec[,1]
OALVS=Resvec[,2:13]

# print tabulated results from the B bootstrap replications of ATE with SIS + OAL
xbar=mean(OALvec)
SIS_OAL=f(ResIni[1],OALvec)

SIS_OAL=as.data.frame(rbind(round(SIS_OAL,3)))
tabres=c("ATE", "Mean","Bias", "SE", "MSE","95%NCIL", "95%NCIU", "Length")
colnames(SIS_OAL)=tabres
rownames(SIS_OAL)="SIS + OAL"
View(SIS_OAL)
print(SIS_OAL)

# print proportion of times covariate selected from the S replications when estimating the ATE with SIS + OAL

OALVS_prop=cbind(colMeans(OALVS))
rownames(OALVS_prop)=var_selected
colnames(OALVS_prop)="Proportion of times covariate selected (%)"
OALVS_prop=as.data.frame(round(OALVS_prop,3)*100)
View(OALVS_prop)
print(OALVS_prop)