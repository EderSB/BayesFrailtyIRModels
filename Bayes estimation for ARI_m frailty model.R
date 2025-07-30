library(tidyverse)
library(coda)
library(msm)

#--------------------------------------------#
#   Estimation for the Frailty ARI_m model   #
#       Hierarchical Bayesian Methods        #
#--------------------------------------------#


# Obs.: The "data" must be in data.frame, matrix or tible format 
# containing at least three necessary information named by: 
# System (numbering of observed systems from 1 to the maximum number of systems), 
# Time (times until recorded failures) and 
# Event (1 for failure and 0 for truncation).

# n is the number of iterations of the chain
# The "parameter"-start are the initial values for the parameters
# m is the chosen memory of the frailty ARI model


#--------------------------------------------#
#             Auxiliary functions            #
#--------------------------------------------#

v.age_components <- function(data, theta, beta, m) {
  K = max(data$System)
  W = NULL
  n_k = NULL
  t.ari.obs=NULL
  for(k in 1:K){
    data.k=data %>% filter(System==k)
    time.k=data.k$Time
    time.obs.k=time.k[-length(time.k)]
    time.k.aux=c(0,time.k)
    t.k.aux=c(rep(0,(m-1)),time.k.aux)
    n_k[k]=length(time.k[-1])
    #
    t.ari.k=NULL
    t.cum.ari.k=NULL
    for(j in 1:length(time.k)){
      t.ari.k[j]=t.k.aux[j+m]^(beta-1) - (1-theta[k])*sum((theta[k]^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
      t.cum.ari.k[j]=t.k.aux[j+m]^(beta) - t.k.aux[j+m-1]^(beta) - 
        beta*(t.k.aux[j+m]-t.k.aux[j+m-1])*(1-theta[k])*sum((theta[k]^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
    }
    t.ari.obs.k=t.ari.k[-length(t.ari.k)]
    
    W[k] = sum(t.cum.ari.k)
    
    t.ari.obs=c(t.ari.obs,t.ari.obs.k)
  }
  
  return(list(W = W, t.ari.obs = t.ari.obs, n_k = n_k))
}


log.pi.theta=function(theta,beta,omega,alpha,g,h,m,time.k){
  time.k.aux=c(0,time.k)
  t.k.aux=c(rep(0,(m-1)),time.k.aux)
  n_k=length(time.k[-1])
  t.ari.k=NULL
  t.cum.ari.k=NULL
  
  for(j in 1:length(time.k)){
    t.ari.k[j]=t.k.aux[j+m]^(beta-1) - (1-theta)*sum((theta^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
    t.cum.ari.k[j]=t.k.aux[j+m]^(beta) - t.k.aux[j+m-1]^(beta) - 
      beta*(t.k.aux[j+m]-t.k.aux[j+m-1])*(1-theta)*sum((theta^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
  }
  t.ari.obs.k=t.ari.k[-length(t.ari.k)]
  #
  resu =sum(log(t.ari.obs.k)) + (g-1)*log(theta) + (h-1)*log(1-theta) - 
    sum((1/alpha+n_k)*log(1+alpha*omega*(t.cum.ari.k)))
  return(resu)
}


log.pi.beta=function(omega,beta,alpha,theta,m,data,a,b){ 
  N=sum(data$Event)
  K=max(data$System)
  
  calc=v.age_components(data=data, theta=theta, beta=beta, m=m)
  
  resu.aux= N*log(beta)+sum(log(calc$t.ari.obs)) -sum((1/alpha+calc$n_k)*log(1+alpha*omega*calc$W))
  resu = ifelse(between(beta,a,b)==TRUE, resu.aux, 0)
  return(resu)
}


log.pi.omega=function(omega,beta,alpha,theta,m,data,c,d){
  N=sum(data$Event)
  K=max(data$System)
  
  calc=v.age_components(data=data, theta=theta, beta=beta, m=m)
  
  resu = N*log(omega) + (c-1)*log(omega) - d*omega - sum((1/alpha+calc$n_k)*log(1+alpha*omega*calc$W))
  return(resu)
}


log.pi.alpha=function(alpha,beta,omega,theta,m,data,e,f){
  N=sum(data$Event)
  K=max(data$System)
  
  calc=v.age_components(data=data, theta=theta, beta=beta, m=m)
  
  resu =  N*log(alpha) - K*lgamma(1/alpha) + sum(lgamma(1/alpha + calc$n_k)) -
    sum((1/alpha+calc$n_k)*log(1+alpha*omega*calc$W)) + (e-1)*log(alpha) - f*alpha
  return(resu)
}


log.pi.g=function(theta,g){
  (g-1)*sum(log(theta))-(g-1)
}


log.pi.h=function(theta,h){
  (h-1)*sum(log(1-theta))-(h-1)
}

#--------------------------------------------#
#               MCMC function                #
#--------------------------------------------#

metr0 <- function(n, data, betastart, omegastart, thetastart, alphastart, gstart, hstart, m,
                      a, b, c, d, e, f) {
  data=data
  m=m
  N=sum(data$Event)
  K=max(data$System)
  #
  alpha = c(alphastart, rep(NA, n - 1))
  beta = c(betastart, rep(NA, n-1))
  omega = c(omegastart, rep(NA, n-1))
  theta = matrix(NA, K, n); theta[ ,1]=thetastart
  g = c(gstart, rep(NA, n-1))
  h = c(hstart, rep(NA, n-1))
  taxa.alpha = 0
  taxa.beta = 0
  taxa.omega = 0
  taxa.theta = c(rep(0,K))
  taxa.g = 0
  taxa.h = 0
  #
  for (i in 2:n) {
    
    # estimate the thetas
    for(k in 1:K){
      data.k=data %>% filter(System==k)
      time.k=data.k$Time
      
      theta.star = rtnorm(1, theta[k,i-1], 0.2, lower=0, upper=1)
      
      log.R.theta = log.pi.theta(theta=theta.star,beta=beta[i-1],omega=omega[i-1],alpha=alpha[i-1],g=g[i-1],h=h[i-1],m=m,time.k=time.k) -
        log.pi.theta(theta=theta[k,i-1],beta=beta[i-1],omega=omega[i-1],alpha=alpha[i-1],g=g[i-1],h=h[i-1],m=m,time.k=time.k)
      
      prob.theta = min(1, exp(log.R.theta))
      
      if(is.finite(prob.theta)){
        if (runif(1) < prob.theta) {
          theta[k,i] = theta.star
          taxa.theta[k] = taxa.theta[k] + 1
        }
        else {
          theta[k,i] = theta[k,i-1]
        }
      }  else { theta[k,i] = theta[k,i-1]}
      
    }
    
    # estimate beta
    s_d=0.5
    beta.star = rtnorm(1, mean=beta[i-1], sd=s_d, lower=1, upper=3) 
    
    log.R.beta = log.pi.beta(beta=beta.star,omega=omega[i-1],alpha=alpha[i-1],theta=theta[,i],m=m,data=data,a=a,b=b) - 
      log.pi.beta(beta=beta[i-1],omega=omega[i-1],alpha=alpha[i-1],theta=theta[,i],m=m,data=data,a=a,b=b)
    
    prob.beta = min(1, exp(log.R.beta))
    
    if(is.finite(prob.beta)){
      if (runif(1) < prob.beta) {
        beta[i] = beta.star
        taxa.beta = taxa.beta + 1
      }
      else {
        beta[i] = beta[i - 1]
      }
    }  else { beta[i] = beta[i-1]}
    
    
    # estimate omega
    omega.star = rtnorm(1, omega[i-1], 0.01, lower=0, upper=0.2)
    
    log.R.omega = log.pi.omega(omega=omega.star,beta=beta[i],alpha=alpha[i-1],theta=theta[,i],m=m,data=data,c=c,d=d) - 
      log.pi.omega(omega=omega[i-1],beta=beta[i],alpha=alpha[i-1],theta=theta[,i],m=m,data=data,c=c,d=d)
    
    prob.omega = min(1, exp(log.R.omega))
    
    if(is.finite(prob.omega)){
      if (runif(1) < prob.omega) {
        omega[i] = omega.star
        taxa.omega = taxa.omega + 1
      }
      else {
        omega[i] = omega[i-1]
      }
    }  else { omega[i] = omega[i-1]}
    
    
    # estimate alpha
    alpha.star = rtnorm(1, alpha[i-1], 1, lower=0, upper=0.5)
    
    log.R.alpha = log.pi.alpha(alpha=alpha.star,beta=beta[i],omega=omega[i],theta=theta[,i],m=m,data=data,e=e,f=f) - 
      log.pi.alpha(alpha=alpha[i-1],beta=beta[i],omega=omega[i],theta=theta[,i],m=m,data=data,e=e,f=f)
    
    prob.alpha = min(1, exp(log.R.alpha))
    
    if(is.finite(prob.alpha)){
      if (runif(1) < prob.alpha) {
        alpha[i] = alpha.star
        taxa.alpha = taxa.alpha + 1
      }
      else {
        alpha[i] = alpha[i-1]
      }
    }  else { alpha[i] = alpha[i-1]}
    
    
    # estimate the hyperparameter "g"
    s_d=1
    g.star = rtnorm(1, mean=g[i-1], sd=s_d, lower=1, upper=3) 
    
    log.R.g = log.pi.g(g=g.star,theta=theta[,i])-log.pi.g(g=g[i-1],theta=theta[,i])
    
    prob.g = min(1, exp(log.R.g))
    
    if(is.finite(prob.g)){
      if (runif(1) < prob.g) {
        g[i] = g.star
        taxa.g = taxa.g + 1
      }
      else {
        g[i] = g[i - 1]
      }
    }  else { g[i] = g[i-1]}
    
    
    # estimate the hyperparameter "h"
    s_d=1
    h.star = rtnorm(1, mean=h[i-1], sd=s_d, lower=1, upper=3) 
    
    
    log.R.h = log.pi.h(h=h.star,theta=theta[,i])-
      log.pi.h(h=h[i-1],theta=theta[,i])
    
    prob.h = min(1, exp(log.R.h))
    
    if(is.finite(prob.h)){
      if (runif(1) < prob.h) {
        h[i] = h.star
        taxa.h = taxa.h + 1
      }
      else {
        h[i] = h[i - 1]
      }
    }  else { h[i] = h[i-1]}
  }
  #
  return(list(beta=beta, omega=omega, theta=theta, alpha=alpha, g=g, h=h, 
              tx.beta=taxa.beta/n, tx.omega=taxa.omega/n, tx.theta=taxa.theta/n, tx.alpha=taxa.alpha/n,
              tx.g=taxa.g/n, tx.h=taxa.h/n))
}


#------------------------------------------------------------------
# An example
# (article application)
#------------------------------------------------------------------

data=read.table("trucks.txt", head=TRUE)

memo=32

mm = metr0(n=50000, data=data, betastart=1.5, omegastart=0.05, thetastart=c(rep(0.5,5)), 
               alphastart=0.2, gstart=1.5, hstart=1.5, m=memo, a=1, b=3, c=0.01, d=1, e=1, f=1)


x.inf=5000; x.sup=50000
jump=10

alpha1=mm$alpha[seq(x.inf,x.sup,jump)]
beta1=mm$beta[seq(x.inf,x.sup,jump)]
omega1=mm$omega[seq(x.inf,x.sup,jump)]
g1=mm$g[seq(x.inf,x.sup,jump)]
h1=mm$h[seq(x.inf,x.sup,jump)]

theta1.1=mm$theta[1,seq(x.inf,x.sup,jump)]
theta1.2=mm$theta[2,seq(x.inf,x.sup,jump)]
theta1.3=mm$theta[3,seq(x.inf,x.sup,jump)]
theta1.4=mm$theta[4,seq(x.inf,x.sup,jump)]
theta1.5=mm$theta[5,seq(x.inf,x.sup,jump)]

chain = data.frame(beta=beta1,omega=omega1,alpha=alpha1,theta1=theta1.1,theta2=theta1.2,
                   theta3=theta1.3,theta4=theta1.4,theta5=theta1.5)

z = as.mcmc(cbind(beta1,omega1,alpha1, g1, h1, theta1.1,theta1.2,theta1.3,theta1.4,theta1.5))
colnames(z) = c("Beta", "Omega", "Alpha", "g1", "h1", "Theta.1","Theta.2","Theta.3","Theta.4","Theta.5")

chain_summary = summary(window(z, start = 1))
print(chain_summary, digits = 3)

beta.est=chain_summary$statistics[1,1]
omega.est=chain_summary$statistics[2,1]
alpha.est=chain_summary$statistics[3,1]
g.est=chain_summary$statistics[4,1]
h.est=chain_summary$statistics[5,1]
theta.est=as.numeric(chain_summary$statistics[6:10,1])

estimates = c(beta.est, omega.est, alpha.est, g.est, h.est, theta.est)
estimates


#------------------------------------------------------------------
# Complementary analyzes
# (article application)
#------------------------------------------------------------------

#-------------------------------------------
# PLP functions
#-------------------------------------------

lambda.est=function(t, beta, omega){
  resu=(beta*omega)*(t^(beta-1))
  return(resu)
}


Lambda.est=function(t, beta, omega){
  resu=omega*(t^beta)
  return(resu)
}



#-------------------------------------------
# ARI_m functions
#-------------------------------------------

lambda.arim.est=function(time, m, beta, omega, theta){
  time.aux=c(0,time)
  lt=c(rep(0,(m-1)),lambda.est(t=time.aux, beta=beta, omega=omega))
  value=NULL
  for(i in 1:length(time)){
    value[i]= lt[i+m]-(1-theta)*sum((theta^(0:(m-1)))*(rev(lt[i:(i+m-1)])))
  }
  return(value)
}


Lambda.arim.est<-function(time,m, beta, omega, theta){
  time.aux=c(0,time)
  time.dif=diff(time.aux)
  time.aux2=time.aux[-length(time.aux)]
  lt=c(rep(0,(m-1)),lambda.est(t=time.aux2, beta=beta, omega=omega))
  sum1=NULL
  for(i in 1:length(time)){
    sum1[i]=time.dif[i]*sum((theta^(0:(m-1)))*(rev(lt[i:(i+m-1)])))
  }
  resu=Lambda.est(t=time, beta=beta, omega=omega)-(1-theta)*cumsum(sum1)
  return(resu)
}

#-------------------------------------------
# Log-likelihood function
#-------------------------------------------

LogLik.shared.ari <- function(par, par.fix=NULL, data, m=1) {
  k=max(data$System)
  
  beta = par[1]
  omega = par[2]
  alpha = par[3]
  theta =par[4:(3+k)]
  
  l.k=NULL
  for(j in 1:k){
    data.k=data %>% filter(System==j)
    time=data.k$Time
    time.aux=c(0,time)
    time.obs=time
    time.ant=time.aux[-length(time.aux)]
    n_i=length(time.ant[-1])
    
    lambda.aux = 
      lambda.arim.est(time=time.ant, omega=omega, beta=beta, theta=theta[j], m=m)
    Lambda.aux = 
      Lambda.arim.est(t=time.obs, omega=omega, beta=beta, theta=theta[j], m=m) - 
      Lambda.arim.est(t=time.ant, omega=omega, beta=beta, theta=theta[j], m=m)
    
    lambda.aux = ifelse(lambda.aux<0, 1, lambda.aux)
    log.lambda.aux = log(lambda.aux)
    log.lambda.aux = ifelse(!is.finite(log.lambda.aux), 0, log.lambda.aux)
    #
    
    l.k[j] = -( sum( log.lambda.aux ) + n_i*log(alpha) + lgamma(1/alpha+n_i) - lgamma(1/alpha)
                -(1/alpha+n_i)*log(1+alpha*sum(Lambda.aux)))
  }
  l=sum(l.k)
  #l=-( sum( log.lambda.aux ) + sum( log.R.aux ) )
  if( !is.finite(l) ) l = 1e50
  return(l %>% as.numeric())
  #
  
}



#-------------------------------------------
# DIC
#-------------------------------------------

N.efect=nrow(chain)

l.chain=NULL
for(i in 1:N.efect){
  l.chain[i] = -LogLik.shared.ari(par=c(chain[i,1], chain[i,2], chain[i,3],chain[i,4],chain[i,5],chain[i,6],
                                        chain[i,7],chain[i,8]), data=data,m=memo)
}

m.l.chain=mean(l.chain)

vero= -LogLik.shared.ari(par=c(beta.est, omega.est, alpha.est, 
                               theta.est[1],theta.est[2],theta.est[3],
                               theta.est[4],theta.est[5]),data=data,m=memo)

DIC = -2*vero + 4*(vero - m.l.chain)
DIC



#-------------------------------------------
# Bayesian diagnostics
#-------------------------------------------

plot(z)
geweke.diag(z)
heidel.diag(z)

#Gelman-Rubin
beta = runif(1,1,3)
omega = runif(1,0,0.1)
alpha = runif(1,0,2)
e=runif(1,1,3)
f=runif(1,1,3)
theta1=runif(1,0,1)
theta2=runif(1,0,1)
theta3=runif(1,0,1)
theta4=runif(1,0,1)
theta5=runif(1,0,1)

mm1=metr0(n=50000, data=data, betastart=beta, omegastart=omega, alphastart=alpha,
          thetastart=c(theta1,theta2,theta3,theta4,theta5), estart=e, fstart=f, m=32)
x.inf=5000
x.sup=50000
jump=10

alpha1=mm1$alpha[seq(x.inf,x.sup,jump)]
beta1=mm1$beta[seq(x.inf,x.sup,jump)]
omega1=mm1$omega[seq(x.inf,x.sup,jump)]
e1=mm1$e[seq(x.inf,x.sup,jump)]
f1=mm1$f[seq(x.inf,x.sup,jump)]

theta1.1=mm1$theta[1,seq(x.inf,x.sup,jump)]
theta1.2=mm1$theta[2,seq(x.inf,x.sup,jump)]
theta1.3=mm1$theta[3,seq(x.inf,x.sup,jump)]
theta1.4=mm1$theta[4,seq(x.inf,x.sup,jump)]
theta1.5=mm1$theta[5,seq(x.inf,x.sup,jump)]

zz = as.mcmc(cbind(beta1,omega1,alpha1, e1, f1, theta1.1,theta1.2,theta1.3,theta1.4,theta1.5))
colnames(zz) = c("beta", "omega", "Alpha", "e1", "f1", "theta.1","theta.2","theta.3","theta.4","theta.5")

combinedchains = mcmc.list(z, zz)
gelman.diag(combinedchains)
