#############################################################
####### code for the proposed method for handling completely missing binary covariate
#############################################################
library(BB)
library(MASS)
source("functions_binary.R")
str1 = "setting 1: two nuisance models are correctly specified."
str2 = "setting 2: only density ratio model is correct, imputation model is mispecified."
str3 = "setting 3: only imputation model is correct, density ratio model is mispecified."
strs = c(str1,str2,str3)
## choose setting 1 or 2 or 3
setting = 3 
# 1. generate data --------
M = 2000 # total sample size of source and target data
## parameter values
gamma = 4*c(-0.3, 0.2, 0.8, -0.8) # gamma.y, gamma.z
eta = c(0.8, 0.3, -0.5, 0.3) # eta.y, eta.z
p = 2; mu_z = rep(0,p); Sigma_z = matrix(NA,p,p)
for(j in 1:p){
  for(l in 1:p){
    Sigma_z[j,l] = 0.3^(abs(j-l))
  }
}
## generate dataset
set.seed(123)
Y = rnorm(M, mean = 1, sd = 1.5)
Z = mvrnorm(M,mu_z,Sigma_z)
if(setting == 1){
  Px = expit(gamma[1] + gamma[2]*Y + Z%*%gamma[3:(p+2)])
  Ps = expit(eta[1] + eta[2]*Y + Z%*%eta[3:(p+2)])
}
if(setting == 2){
  Px = expit( gamma[1] + gamma[2]*Y + Z%*%gamma[3:(p+2)] + 3*Y*as.vector(Z[,1]) )
  Ps = expit(eta[1] + eta[2]*Y + Z%*%eta[3:(p+2)])
}
if(setting == 3){
  Px = expit(gamma[1] + gamma[2]*Y + Z%*%gamma[3:(p+2)])
  Ps = expit( eta[1] + eta[2]*Y + Z%*%eta[3:(p+2)] + 2*Y*as.vector(Z[,1]) )
}

X = rbinom(M,size=1,prob=Px) # binary
S = rbinom(M,size=1,prob=Ps) 
ns = sum(S) # sample size for source data
nt = M - ns # sample size for target data
data = cbind(X,Y,Z,S) 

# observed source data
Xs = data[which(data[,p+3]==1),1]
Ys = data[which(data[,p+3]==1),2]
Zs = data[which(data[,p+3]==1),3:(p+2)]
# observed target data
Yt = data[which(data[,p+3]==0),2]
Zt = data[which(data[,p+3]==0),3:(p+2)]

# 2. estimate nuisance models --------
## estimated correct density ratio model: omega(Ys,Zs)
fit = glm(S ~ Y + Z, family = "binomial")
eta_hat = fit$coefficients
omega_hat = exp(-eta_hat[1] - eta_hat[2]*Ys - Zs%*%eta_hat[3:(p+2)])

## estimated correct imputation model: m_i(Ys,Zs), m_i(Yt,Zt), i = 1,2
fit = glm(Xs ~ Ys + Zs, family = "binomial")
gamma_hat = fit$coefficients
m.s = expit(gamma_hat[1] + gamma_hat[2]*Ys + Zs%*%gamma_hat[3:(p+2)]) # m(Ys,Zs)
m.t = expit(gamma_hat[1] + gamma_hat[2]*Yt + Zt%*%gamma_hat[3:(p+2)]) # m(Yt,Zt)

# 3. estimate parameters using the proposed methods --------
## input
data = data # data frame of (X,Y,Z,S): source data (S=1) and target data (S=0)
omega_hat = omega_hat # estimated density ratio
m.s = m.s # imputed X for source data
m.t = m.t # imputed X for target data
## run our method 
result = fn_dr(data,omega_hat,m.s,m.t)
## output
print(strs[setting])
round(result,3) # point estimators and standard errors for three methods: IW, IM, and DR (proposed doubly robust method)




