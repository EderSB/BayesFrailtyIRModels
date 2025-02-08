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
# Th "parameter"-start are the initial values for the parameters
# m is the chosen memory of the frailty ARI model
# "eta" in this code refers to the parameter "omega"

metr0 <- function(n, data, betastart, etastart, thetastart, alphastart, estart, fstart, m) {
  data=data
  m=m
  N=sum(data$Event)
  K=max(data$System)
  #
  alpha = c(alphastart, rep(NA, n - 1))
  beta = c(betastart, rep(NA, n-1))
  eta = c(etastart, rep(NA, n-1))
  theta = matrix(NA, K, n); theta[ ,1]=thetastart
  e = c(estart, rep(NA, n-1))
  f = c(fstart, rep(NA, n-1))
  taxa.alpha = 0
  taxa.beta = 0
  taxa.eta = 0
  taxa.theta = c(rep(0,K))
  taxa.e = 0
  taxa.f = 0
  #
  for (i in 2:n) {
    
    # estimate the thetas
    for(k in 1:K){
      data.k=data %>% filter(System==k)
      time.k=data.k$Time
      time.obs.k=time.k[-length(time.k)]
      time.k.aux=c(0,time.k)
      time.prev.k=c(0,time.obs.k)
      t.k.aux=c(rep(0,(m-1)),time.k.aux)
      n_k=length(time.k[-1])
      
      log.pi.theta=function(theta,beta,eta,alpha,e,f){
        t.ari.k=NULL
        t.cum.ari.k=NULL
        for(j in 1:length(time.k)){
          t.ari.k[j]=t.k.aux[j+m]^(beta-1) - (1-theta)*sum((theta^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
          t.cum.ari.k[j]=t.k.aux[j+m]^(beta) - t.k.aux[j+m-1]^(beta) - 
            beta*(t.k.aux[j+m]-t.k.aux[j+m-1])*(1-theta)*sum((theta^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
        }
        t.ari.obs.k=t.ari.k[-length(t.ari.k)]
        #
        resu =sum(log(t.ari.obs.k)) + (e-1)*log(theta) + (f-1)*log(1-theta) - 
          sum((1/alpha+n_k)*log(1+alpha*eta*(t.cum.ari.k)))
        return(resu)
      }
      
      theta.star = rtnorm(1, theta[k,i-1], 0.2, lower=0, upper=1)
      
      log.R.theta = log.pi.theta(theta=theta.star,beta=beta[i-1],eta=eta[i-1],alpha=alpha[i-1],e=e[i-1],f=f[i-1]) -
        log.pi.theta(theta=theta[k,i-1],beta=beta[i-1],eta=eta[i-1],alpha=alpha[i-1],e=e[i-1],f=f[i-1])
    
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
    log.pi.beta=function(eta,beta,alpha){ 
      W=NULL
      n_k=NULL
      t.ari.obs=NULL
      for(k in 1:K){
        data.k=data %>% filter(System==k)
        time.k=data.k$Time
        time.obs.k=time.k[-length(time.k)]
        time.k.aux=c(0,time.k)
        time.prev.k=c(0,time.obs.k)
        t.k.aux=c(rep(0,(m-1)),time.k.aux)
        n_k[k]=length(time.k[-1])
        #
        t.ari.k=NULL
        t.cum.ari.k=NULL
        for(j in 1:length(time.k)){
          t.ari.k[j]=t.k.aux[j+m]^(beta-1) - (1-theta[k,i])*sum((theta[k,i]^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
          t.cum.ari.k[j]=t.k.aux[j+m]^(beta) - t.k.aux[j+m-1]^(beta) - 
            beta*(t.k.aux[j+m]-t.k.aux[j+m-1])*(1-theta[k,i])*sum((theta[k,i]^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
        }
        t.ari.obs.k=t.ari.k[-length(t.ari.k)]
        
        W[k] = sum(t.cum.ari.k)
        
        t.ari.obs=c(t.ari.obs,t.ari.obs.k)
      }
      
      resu.aux= N*log(beta)+sum(log(t.ari.obs)) -sum((1/alpha+n_k)*log(1+alpha*eta*W))
      resu = ifelse(between(beta,a,b)==TRUE, resu.aux, 0)
      return(resu)
    }
    
    s_d=0.5
    beta.star = rtnorm(1, mean=beta[i-1], sd=s_d, lower=1, upper=3) 
    
    log.R.beta = log.pi.beta(beta=beta.star,eta=eta[i-1],alpha=alpha[i-1]) - 
      log.pi.beta(beta=beta[i-1],eta=eta[i-1],alpha=alpha[i-1])
    
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
    
    
    # estimate eta
    log.pi.eta=function(eta,beta,alpha){
      W=NULL
      n_k=NULL
      t.ari.obs=NULL
      for(k in 1:K){
        data.k=data %>% filter(System==k)
        time.k=data.k$Time
        time.obs.k=time.k[-length(time.k)]
        time.k.aux=c(0,time.k)
        time.prev.k=c(0,time.obs.k)
        t.k.aux=c(rep(0,(m-1)),time.k.aux)
        n_k[k]=length(time.k[-1])
        #
        t.ari.k=NULL
        t.cum.ari.k=NULL
        for(j in 1:length(time.k)){
          t.ari.k[j]=t.k.aux[j+m]^(beta-1) - (1-theta[k,i])*sum((theta[k,i]^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
          t.cum.ari.k[j]=t.k.aux[j+m]^(beta) - t.k.aux[j+m-1]^(beta) - 
            beta*(t.k.aux[j+m]-t.k.aux[j+m-1])*(1-theta[k,i])*sum((theta[k,i]^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
        }
        t.ari.obs.k=t.ari.k[-length(t.ari.k)]
        
        W[k] = sum(t.cum.ari.k)
        
        t.ari.obs=c(t.ari.obs,t.ari.obs.k)
      }
      
      resu = N*log(eta) + (c-1)*log(eta) - d*eta - sum((1/alpha+n_k)*log(1+alpha*eta*W))
      return(resu)
    }
    
    eta.star = rtnorm(1, eta[i-1], 0.01, lower=0, upper=0.2)
    
    log.R.eta = log.pi.eta(eta=eta.star,beta=beta[i],alpha=alpha[i-1]) - 
      log.pi.eta(eta=eta[i-1],beta=beta[i],alpha=alpha[i-1])
    
    prob.eta = min(1, exp(log.R.eta))
    
    if(is.finite(prob.eta)){
      if (runif(1) < prob.eta) {
        eta[i] = eta.star
        taxa.eta = taxa.eta + 1
      }
      else {
        eta[i] = eta[i-1]
      }
    }  else { eta[i] = eta[i-1]}
    
    
    # estimate alpha
    log.pi.alpha=function(alpha,beta,eta){
      W=NULL
      n_k=NULL
      t.ari.obs=NULL
      soma.obs=NULL
      for(k in 1:K){
        data.k=data %>% filter(System==k)
        time.k=data.k$Time
        time.obs.k=time.k[-length(time.k)]
        time.prev.k=c(0,time.obs.k)
        time.k.aux=c(0,time.k)
        t.k.aux=c(rep(0,(m-1)),time.k.aux)
        n_k[k]=length(time.k[-1])
        #
        t.ari.k=NULL
        t.cum.ari.k=NULL
        for(j in 1:length(time.k)){
          t.ari.k[j]=t.k.aux[j+m]^(beta-1) - (1-theta[k,i])*sum((theta[k,i]^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
          t.cum.ari.k[j]=t.k.aux[j+m]^(beta) - t.k.aux[j+m-1]^(beta) - 
            beta*(t.k.aux[j+m]-t.k.aux[j+m-1])*(1-theta[k,i])*sum((theta[k,i]^(0:(m-1)))*(rev((t.k.aux[j:(j+m-1)])^(beta-1))))
        }
        t.ari.obs.k=t.ari.k[-length(t.ari.k)]
        
        W[k] = sum(t.cum.ari.k)
        
        t.ari.obs=c(t.ari.obs,t.ari.obs.k)
      }
      
      resu =  N*log(alpha) - K*lgamma(1/alpha) + sum(lgamma(1/alpha + n_k)) -
        sum((1/alpha+n_k)*log(1+alpha*eta*W)) + (r-1)*log(alpha) - s*alpha
      return(resu)
    }
    
    alpha.star = rtnorm(1, alpha[i-1], 1, lower=0, upper=2)
    
    log.R.alpha = log.pi.alpha(alpha=alpha.star,beta=beta[i],eta=eta[i]) - 
      log.pi.alpha(alpha=alpha[i-1],beta=beta[i],eta=eta[i])
    
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
    
    # estimate the hyperparameter "e"
    s_d=1
    e.star = rtnorm(1, mean=e[i-1], sd=s_d, lower=1, upper=5) 
    
    log.pi.e=function(theta,e){
      (e-1)*sum(log(theta))-(e-1)
    }
    
    log.R.e = log.pi.e(e=e.star,theta=theta[,i]) - 
      log.pi.e(e=e[i-1],theta=theta[,i])
    
    prob.e = min(1, exp(log.R.e))
    
    if(is.finite(prob.e)){
      if (runif(1) < prob.e) {
        e[i] = e.star
        taxa.e = taxa.e + 1
      }
      else {
        e[i] = e[i - 1]
      }
    }  else { e[i] = e[i-1]}
    
    
    # estimate the hyperparameter "f"
    s_d=1
    f.star = rtnorm(1, mean=f[i-1], sd=s_d, lower=1, upper=5) 
    
    log.pi.f=function(theta,f){
      (f-1)*sum(log(1-theta))-(f-1)
    }
    
    log.R.f = log.pi.f(f=f.star,theta=theta[,i])-
      log.pi.f(f=f[i-1],theta=theta[,i])
    
    prob.f = min(1, exp(log.R.f))
    
    if(is.finite(prob.f)){
      if (runif(1) < prob.f) {
        f[i] = f.star
        taxa.f = taxa.f + 1
      }
      else {
        f[i] = f[i - 1]
      }
    }  else { f[i] = f[i-1]}
  }
  #
  return(list(beta=beta, eta=eta, theta=theta, alpha=alpha, e=e, f=f, 
              tx.beta=taxa.beta/n, tx.eta=taxa.eta/n, tx.theta=taxa.theta/n, tx.alpha=taxa.alpha/n,
              tx.e=taxa.e/n, tx.f=taxa.f/n))
}



#------------------------------------------------------------------
# An example
# (article application)
#------------------------------------------------------------------

data=read.table("trucks.txt", head=TRUE)

a=1; b=3; c=0.01; d=1; g=1; h=1; r=1; s=1

memo=32

mm = metr0(n=50000, data=data, betastart=1.5, etastart=0.01, thetastart=c(rep(0.1,5)), 
           alphastart=0.5, estart=1, fstart=1, m=memo)

x.inf=5000; x.sup=50000
jump=10

alpha1=mm$alpha[seq(x.inf,x.sup,jump)]
beta1=mm$beta[seq(x.inf,x.sup,jump)]
eta1=mm$eta[seq(x.inf,x.sup,jump)]
e1=mm$e[seq(x.inf,x.sup,jump)]
f1=mm$f[seq(x.inf,x.sup,jump)]

theta1.1=mm$theta[1,seq(x.inf,x.sup,jump)]
theta1.2=mm$theta[2,seq(x.inf,x.sup,jump)]
theta1.3=mm$theta[3,seq(x.inf,x.sup,jump)]
theta1.4=mm$theta[4,seq(x.inf,x.sup,jump)]
theta1.5=mm$theta[5,seq(x.inf,x.sup,jump)]

chain = data.frame(beta=beta1,eta=eta1,alpha=alpha1,theta1=theta1.1,theta2=theta1.2,
                   theta3=theta1.3,theta4=theta1.4,theta5=theta1.5)

z = as.mcmc(cbind(beta1,eta1,alpha1, e1, f1, theta1.1,theta1.2,theta1.3,theta1.4,theta1.5))
colnames(z) = c("Beta", "Eta", "Alpha", "e1", "f1", "Theta.1","Theta.2","Theta.3","Theta.4","Theta.5")

chain_summary = summary(window(z, start = 1))
print(chain_summary, digits = 3)

beta.est=chain_summary$statistics[1,1]
eta.est=chain_summary$statistics[2,1]
alpha.est=chain_summary$statistics[3,1]
e.est=chain_summary$statistics[4,1]
f.est=chain_summary$statistics[5,1]
theta.est=as.numeric(chain_summary$statistics[6:10,1])

estimates = c(beta.est, eta.est, alpha.est, e.est, f.est, theta.est)
estimates



#------------------------------------------------------------------
# Complementary analyzes
# (article application)
#------------------------------------------------------------------

#-------------------------------------------
# PLP functions
#-------------------------------------------

lambda.est=function(t, beta, eta){
  resu=(beta*eta)*(t^(beta-1))
  return(resu)
}


Lambda.est=function(t, beta, eta){
  resu=eta*(t^beta)
  return(resu)
}



#-------------------------------------------
# ARI_m functions
#-------------------------------------------

lambda.arim.est=function(time, m, beta, eta, theta){
  time.aux=c(0,time)
  lt=c(rep(0,(m-1)),lambda.est(t=time.aux, beta=beta, eta=eta))
  value=NULL
  for(i in 1:length(time)){
    value[i]= lt[i+m]-(1-theta)*sum((theta^(0:(m-1)))*(rev(lt[i:(i+m-1)])))
  }
  return(value)
}


Lambda.arim.est<-function(time,m, beta, eta, theta){
  time.aux=c(0,time)
  time.dif=diff(time.aux)
  time.aux2=time.aux[-length(time.aux)]
  lt=c(rep(0,(m-1)),lambda.est(t=time.aux2, beta=beta, eta=eta))
  sum1=NULL
  for(i in 1:length(time)){
    sum1[i]=time.dif[i]*sum((theta^(0:(m-1)))*(rev(lt[i:(i+m-1)])))
  }
  resu=Lambda.est(t=time, beta=beta, eta=eta)-(1-theta)*cumsum(sum1)
  return(resu)
}

#-------------------------------------------
# Log-likelihood function
#-------------------------------------------

LogLik.shared.ari <- function(par, par.fix=NULL, data, m=1) {
  k=max(data$System)
  
  beta = par[1]
  eta = par[2]
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
      lambda.arim.est(time=time.ant, eta=eta, beta=beta, theta=theta[j], m=m)
    Lambda.aux = 
      Lambda.arim.est(t=time.obs, eta=eta, beta=beta, theta=theta[j], m=m) - 
      Lambda.arim.est(t=time.ant, eta=eta, beta=beta, theta=theta[j], m=m)
    
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

vero= -LogLik.shared.ari(par=c(beta.est, eta.est, alpha.est, 
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

#Gelman-Rubi
beta = runif(1,1,3)
eta = runif(1,0,0.1)
alpha = runif(1,0,2)
e=runif(1,1,3)
f=runif(1,1,3)
theta1=runif(1,0,1)
theta2=runif(1,0,1)
theta3=runif(1,0,1)
theta4=runif(1,0,1)
theta5=runif(1,0,1)

mm1=metr0(n=50000, data=data, betastart=beta, etastart=eta, alphastart=alpha,
          thetastart=c(theta1,theta2,theta3,theta4,theta5), estart=e, fstart=f, m=32)
x.inf=5000
x.sup=50000
jump=10

alpha1=mm1$alpha[seq(x.inf,x.sup,jump)]
beta1=mm1$beta[seq(x.inf,x.sup,jump)]
eta1=mm1$eta[seq(x.inf,x.sup,jump)]
e1=mm1$e[seq(x.inf,x.sup,jump)]
f1=mm1$f[seq(x.inf,x.sup,jump)]

theta1.1=mm1$theta[1,seq(x.inf,x.sup,jump)]
theta1.2=mm1$theta[2,seq(x.inf,x.sup,jump)]
theta1.3=mm1$theta[3,seq(x.inf,x.sup,jump)]
theta1.4=mm1$theta[4,seq(x.inf,x.sup,jump)]
theta1.5=mm1$theta[5,seq(x.inf,x.sup,jump)]

zz = as.mcmc(cbind(beta1,eta1,alpha1, e1, f1, theta1.1,theta1.2,theta1.3,theta1.4,theta1.5))
colnames(zz) = c("Beta", "Eta", "Alpha", "e1", "f1", "Theta.1","Theta.2","Theta.3","Theta.4","Theta.5")

combinedchains = mcmc.list(z, zz)
gelman.diag(combinedchains)
