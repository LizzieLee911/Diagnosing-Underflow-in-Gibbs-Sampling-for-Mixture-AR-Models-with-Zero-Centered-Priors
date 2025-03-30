rm(list=ls(all=TRUE))
library(MCMCpack)
library(mvtnorm)
library(rstudioapi)
library(ggplot2)
library(scales)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
script_path <- getSourceEditorContext()$path
script_dir <- dirname(script_path)

file_path <- file.path(script_dir, "GoogleSearchIndex.txt")
dat <- read.delim(file_path)

dat$Week <- as.Date(dat$Week)
y=as.numeric(dat$covid)
print(length(y))

n.all = length(y)
p=1 ## order of AR process

#### Model setup
model_comp_mix=function(tot_num_comp){
  
  ## hyperparameters
  m0=matrix(rep(0,p),ncol=1) ## p is the order of AR process
  C0=10*diag(p)
  C0.inv=0.1*diag(p)
  n0=0.02
  d0=0.02
  K=tot_num_comp ## let the number of mixing component to vary by input
  a=rep(1,K)
  
  Y=matrix(y[(p+1):n.all],ncol=1) ## y_{p+1:T} n.all is the value of T
  
  Fmtx = sapply(1:p, function(k) {
    y[(p + 1 - k):(n.all - k)]
  })
  Fmtx = t(Fmtx)
  
  n=length(Y)
  
  
  ## The code below is used to obtain posterior samples of model parameters
  ## Just copy from the last lecture
  
  sample_omega=function(L.cur){
    n.vec=sapply(1:K, function(k){sum(L.cur==k)})
    rdirichlet(1,a+n.vec)
  }
  
  sample_L_one=function(beta.cur,omega.cur,nu.cur,y.cur,Fmtx.cur){
    prob_k=function(k){
      beta.use=beta.cur[((k-1)*p+1):(k*p)]
      omega.cur[k]*dnorm(y.cur,mean=sum(beta.use*Fmtx.cur),sd=sqrt(nu.cur))
    }
    prob.vec=sapply(1:K, prob_k)
    
    s.pv = sum(prob.vec)
    if (is.na(s.pv) || s.pv <= 0) {
      # 若真的发生概率向量全为0、或出现NA/NaN
      # 可以做一些“fallback”处理，避免报错
      prob.vec = rep(1, K)  # 或者只给某个成分 1、其它 0
      s.pv = K
    }
    
    L.sample=sample(1:K,1,prob=prob.vec/sum(prob.vec))
    return(L.sample)
  }
  
  sample_L=function(y,x,beta.cur,omega.cur,nu.cur){
    L.new=sapply(1:n, function(j){sample_L_one(beta.cur,omega.cur,nu.cur,y.cur=y[j,],Fmtx.cur=x[,j])})
    return(L.new)
  }
  
  sample_nu=function(L.cur,beta.cur){
    n.star=n0+n+p*K
    err.y=function(idx){
      L.use=L.cur[idx]
      beta.use=beta.cur[((L.use-1)*p+1):(L.use*p)]
      err=Y[idx,]-sum(Fmtx[,idx]*beta.use)
      return(err^2)
    }
    err.beta=function(k.cur){
      beta.use=beta.cur[((k.cur-1)*p+1):(k.cur*p)]
      beta.use.minus.m0=matrix(beta.use,ncol=1)-m0
      t(beta.use.minus.m0)%*%C0.inv%*%beta.use.minus.m0
    }
    
    d.star=d0+sum(sapply(1:n,err.y))+sum(sapply(1:K,err.beta))
    1/rgamma(1,shape=n.star/2,rate=d.star/2)
  }
  
  
  sample_beta = function(k, L.cur, nu.cur) {
    idx.select = (L.cur == k)
    n.k = sum(idx.select)
    if (n.k == 0) {
      m.k = m0
      C.k = C0
    } else {
      y.tilde.k = Y[idx.select, , drop = FALSE]  # 保持矩阵结构
      Fmtx.tilde.k = Fmtx[, idx.select, drop = FALSE]  # 保持矩阵结构
      e.k = y.tilde.k - t(Fmtx.tilde.k) %*% m0
      Q.k = t(Fmtx.tilde.k) %*% C0 %*% Fmtx.tilde.k + diag(n.k)
      Q.k.inv = chol2inv(chol(Q.k))
      A.k = C0 %*% Fmtx.tilde.k %*% Q.k.inv
      m.k = m0 + A.k %*% e.k
      C.k = C0 - A.k %*% Q.k %*% t(A.k)
    }
    
    rmvnorm(1, m.k, nu.cur * C.k)
  }
  
  nsim=10000
  
  ## store parameters
  
  beta.mtx=matrix(0,nrow=p*K,ncol=nsim)
  L.mtx=matrix(0,nrow=n,ncol=nsim)
  omega.mtx=matrix(0,nrow=K,ncol=nsim)
  nu.vec=rep(0,nsim)
  
  ## initial value
  
  beta.cur=rep(0,p*K)
  L.cur=rep(1,n)
  omega.cur=rep(1/K,K)
  nu.cur=1
  
  ## Gibbs Sampler
  
  for (i in 1:nsim) {
    set.seed(i)
    
    ## sample omega
    omega.cur=sample_omega(L.cur)
    omega.mtx[,i]=omega.cur
    
    ## sample L
    L.cur=sample_L(Y,Fmtx,beta.cur,omega.cur,nu.cur)
    L.mtx[,i]=L.cur
    
    ## sample nu
    nu.cur=sample_nu(L.cur,beta.cur)
    nu.vec[i]=nu.cur
    
    ## sample beta
    beta.cur=as.vector(sapply(1:K, function(k){sample_beta(k,L.cur,nu.cur)}))
    beta.mtx[,i]=beta.cur
    
    if(i%%1000==0){
      print(i)
    }
    
  }
  
  ## Now compute DIC for mixture model
  ## Somilar as the calculation of DIC in Module 2
  
  cal_log_likelihood_mix_one=function(idx,beta,nu,omega){
    norm.lik=rep(0,K)
    for (k.cur in 1:K) {
      mean.norm=sum(Fmtx[,idx]*beta[((k.cur-1)*p+1):(k.cur*p)])
      norm.lik[k.cur]=dnorm(Y[idx,1],mean.norm,sqrt(nu),log=FALSE)
    }
    log.lik=log(sum(norm.lik*omega))
    return(log.lik)
  }
  
  cal_log_likelihood_mix=function(beta,nu,omega){
    sum(sapply(1:n, function(idx){cal_log_likelihood_mix_one(idx=idx,beta=beta,nu=nu,omega=omega)}))
  }
  
  sample.select.idx=seq(5001,10000,by=1)
  
  beta.mix=rowMeans(beta.mtx[,sample.select.idx])
  nu.mix=mean(nu.vec[sample.select.idx])
  omega.mix=rowMeans(omega.mtx[,sample.select.idx])
  
  log.lik.bayes.mix=cal_log_likelihood_mix(beta.mix,nu.mix,omega.mix)
  
  post.log.lik.mix=sapply(sample.select.idx, function(k){cal_log_likelihood_mix(beta.mtx[,k],nu.vec[k],omega.mtx[,k])})
  E.post.log.lik.mix=mean(post.log.lik.mix)
  
  p_DIC.mix=2*(log.lik.bayes.mix-E.post.log.lik.mix)
  
  DIC.mix=-2*log.lik.bayes.mix+2*p_DIC.mix
  
  return(DIC.mix)
}


mix.model.all=sapply(2:5,model_comp_mix)

print(mix.model.all)
