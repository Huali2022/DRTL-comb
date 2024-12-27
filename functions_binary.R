## functions for completely missing binary X: dim(Z)>=1
library(BB) # dfsane: Derivative-Free Spectral Approach for solving nonlinear systems of equations
expit<-function(x){
  1/(1+exp(-x))
}
# Comparable method 1: Importance Weighting
ImportanceWeight_Only<-function(omega, beta_init=NULL){
  # weighted estimating equations
  augeq<-function(beta){ 
    temp = t(cbind(Xs,1,Zs)) %*% (omega*(Ys - beta[1] - Xs*beta[2] - as.vector(Zs%*%beta[3:(p+2)]) ))/ns
    temp[,1]
  }
  if(is.null(beta_init)){
    beta_init = rep(0,p+2)
  }
  beta_iw = dfsane(beta_init,augeq,control=list(maxit=2500,trace = FALSE))$par
  names(beta_iw) = c("Intercept", "Xt", paste0("Z",1:p,"t"))
  return(beta_iw)
}
# Comparable method 2: Imputation 
Imputation_Only<-function(m.t, beta_init=NULL){
  augeq<-function(beta){
    temp1 = t(m.t) %*% (Yt - beta[1] - as.vector(Zt%*%beta[3:(p+2)]) )/nt - mean(m.t)*beta[2]
    temp2 = t(cbind(1,Zt)) %*% (Yt - beta[1] - m.t*beta[2] - as.vector(Zt%*%beta[3:(p+2)]) )/nt
    temp = c(temp1, temp2)
    return(temp)
  }
  if(is.null(beta_init)){
    beta_init = rep(0,p+2)
  }
  beta_imp = dfsane(beta_init,augeq,control=list(maxit=2500,trace = FALSE))$par
  names(beta_imp) = c("Intercept", "Xt", paste0("Z",1:p,"t"))
  return(beta_imp)
}
# Proposed Doubly Robust method
Doubly_Robust <- function(omega, m.s, m.t, beta_init=NULL){
  augeq<-function(beta){
    temp1 = t(Xs - m.s) %*% (omega*(Ys - beta[1] - as.vector(Zs%*%beta[3:(p+2)]))) / ns + t(m.s - Xs) %*% omega * beta[2] / ns 
    temp1 = temp1 + t(m.t) %*% (Yt - beta[1] - as.vector(Zt%*%beta[3:(p+2)]) ) / nt - mean(m.t)*beta[2]
    temp2 = t(cbind(1,Zs)) %*% (omega*(m.s - Xs)) * beta[2] / ns
    temp2 = temp2 + t(cbind(1,Zt)) %*% (Yt - beta[1] - m.t*beta[2] - as.vector(Zt%*%beta[3:(p+2)]) ) / nt
    temp = c(temp1, temp2)
    return(temp)
  }
  if(is.null(beta_init)){
    beta_init = rep(0,p+2) 
  }
  beta_dr = dfsane(beta_init,augeq,control=list(maxit=2500,trace = FALSE))$par
  names(beta_dr) = c("Intercept", "Xt", paste0("Z",1:p,"t"))
  return(beta_dr)
}

# functions for estimating parameters using above three proposed methods
fn_dr = function(data,omega_hat,m.s,m.t){
  B = 500  # 500 Bootstrap for variance estimation
  p = dim(data)[2] - 3 # dim of Z
  beta_naive = beta_iw = beta_im = beta_dr = rep(NA,p+2)
  b.se.naive = b.se.iw = b.se.im = b.se.dr = rep(NA,p+2) # intercept, Xt, Zt
  b.beta.naive = b.beta.iw= b.beta.im = b.beta.dr = matrix(NA,B,p+2) # estimator from bootstrap
  
  X = data[,1]; Y = data[,2]; Z = data[,3:(p+2)]; S = data[,p+3]
  ## source data
  Xs = data[which(data[,p+3]==1),1]
  Ys = data[which(data[,p+3]==1),2]
  Zs = data[which(data[,p+3]==1),3:(p+2)]
  ## target data
  Xt = data[which(data[,p+3]==0),1] # artificially missing 
  Yt = data[which(data[,p+3]==0),2]
  Zt = data[which(data[,p+3]==0),3:(p+2)]
  
  X0 = X; Y0 = Y; Z0 = Z; S0 = S
  Xs0 = Xs; Ys0 = Ys; Zs0 = Zs
  Xt0 = Xt; Yt0 = Yt; Zt0 = Zt
  
  ## naive method: using source data
  fit <- lm(Ys ~ Xs + Zs)
  beta_naive = fit$coefficients
  
  ## IW method
  beta_iw <- tryCatch({
    ImportanceWeight_Only(omega_hat, beta_init=beta_naive)
  }, warning = function(w) {
    message("Warning：", conditionMessage(w))
    return(NA)
  })
  ## IM method
  beta_im <- tryCatch({
    Imputation_Only(m.t, beta_init=beta_naive)
  }, warning = function(w) {
    message("Warning：", conditionMessage(w))
    return(NA)
  })
  ## Proposed method
  beta_dr <- tryCatch({
    Doubly_Robust(omega_hat, m.s, m.t, beta_init=beta_naive) 
  }, warning = function(w) {
    message("Warning：", conditionMessage(w))
    return(NA)
  })
  
  ## bootstrap for variance estimation------------
  for(b in 1:B){
    b.s.data = matrix(0,ns,p+2) # source: X,Y,Z,S
    b.s.data[,1] = Xs0
    b.s.data[,2] = Ys0
    b.s.data[,3:(p+2)] = Zs0
    index = sample(1:ns,size = ns,replace = TRUE)
    b.s.data = b.s.data[index,]
    Xs = b.s.data[,1]
    Ys = b.s.data[,2]
    Zs = b.s.data[,3:(p+2)]
    
    b.t.data = matrix(0,nt,p+2) # target: X (artificially missing),Y,Z,S
    b.t.data[,1] = Xt0
    b.t.data[,2] = Yt0
    b.t.data[,3:(p+2)] = Zt0
    index = sample(1:nt,size = nt,replace = TRUE)
    b.t.data = b.t.data[index,]
    Xt = b.t.data[,1] # artificially missing
    Yt = b.t.data[,2]
    Zt = b.t.data[,3:(p+2)]
    
    X = c(Xs,Xt); Y = c(Ys,Yt); Z = rbind(Zs,Zt); S = c(rep(1,ns),rep(0,nt))
    
    ## nuisance models--------------
    # estimated correct density ratio model: omega(Ys,Zs)
    fit = glm(S ~ Y + Z, family = "binomial")
    eta_hat = fit$coefficients
    omega_hat = exp(-eta_hat[1] - eta_hat[2]*Ys - Zs%*%eta_hat[3:(p+2)])
    
    # estimated correct imputation model: m(Ys,Zs), m(Yt,Zt)
    fit = lm(Xs ~ Ys + Zs)
    s_eps_hat = sd(fit$residuals)
    gamma_hat = fit$coefficients
    m.s = expit(gamma_hat[1] + gamma_hat[2]*Ys + Zs%*%gamma_hat[3:(p+2)]) # m(Ys,Zs)
    m.t = expit(gamma_hat[1] + gamma_hat[2]*Yt + Zt%*%gamma_hat[3:(p+2)]) # m(Yt,Zt)
    
    ## different methods--------------
    fit <- lm(Ys ~ Xs + Zs)
    b.beta.naive[b,] = fit$coefficients
    
    b.beta.iw[b,] <- tryCatch({
      ImportanceWeight_Only(omega_hat, beta_init=b.beta.naive[b,])
    }, warning = function(w) {
      #message("Warning：", conditionMessage(w))
      return(NA)
    })
    b.beta.im[b,] <- tryCatch({
      Imputation_Only(m.t, beta_init=b.beta.naive[b,]) 
    }, warning = function(w) {
      #message("Warning：", conditionMessage(w))
      return(NA)
    })
    b.beta.dr[b,] <- tryCatch({
      Doubly_Robust(omega_hat, m.s, m.t, beta_init=b.beta.naive[b,]) 
    }, warning = function(w) {
      #message("Warning：", conditionMessage(w))
      return(NA)
    })
  }
  for(j in 1:(p+2)){
    b.se.naive[j] = sd(na.omit(b.beta.naive[,j]))
    b.se.iw[j]    = sd(na.omit(b.beta.iw[,j]))
    b.se.im[j]    = sd(na.omit(b.beta.im[,j]))
    b.se.dr[j]    = sd(na.omit(b.beta.dr[,j]))
  }
  
  names(beta_iw) = c("Intercept", "Xt", paste0("Z",1:p,"t"))
  result= data.frame(estimate.iw = beta_iw, estimate.im = beta_im, estimate.dr = beta_dr, 
                     se.iw = b.se.iw, se.im = b.se.iw, se.dr = b.se.dr)
  return(result)
}
