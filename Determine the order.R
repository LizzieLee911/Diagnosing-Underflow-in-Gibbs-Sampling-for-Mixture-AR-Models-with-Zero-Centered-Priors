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

y.sample=as.numeric(dat$covid)
print(length(y.sample))

plot(y.sample,type='l',xlab='time',ylab='')


par(mfrow=c(1,2))

acf(y.sample,main="",xlab='Lag')
pacf(y.sample,main="",xlab='Lag')

n.all = length(y.sample)       # dataset length
p.star = 15                    # max order
Y = matrix(y.sample[(p.star+1):n.all], ncol=1)
sample.all = matrix(y.sample, ncol=1)
n = length(Y)                  # （100 - 15 = 85）
p = seq(1, p.star, by=1)       # p = 1 ---- p = 15

#Design Matrix construction
design.mtx=function(p_cur){
  Fmtx=matrix(0,ncol=n,nrow=p_cur)
  for (i in 1:p_cur) {
    start.y=p.star+1-i
    end.y=start.y+n-1
    Fmtx[i,]=sample.all[start.y:end.y,1]
  }
  return(Fmtx)
}

criteria.ar=function(p_cur){
  Fmtx=design.mtx(p_cur)
  beta.hat=chol2inv(chol(Fmtx%*%t(Fmtx)))%*%Fmtx%*%Y
  R=t(Y-t(Fmtx)%*%beta.hat)%*%(Y-t(Fmtx)%*%beta.hat)
  sp.square=R/(n-p_cur)
  aic=2*p_cur+n*log(sp.square)
  bic=log(n)*p_cur+n*log(sp.square)
  result=c(aic,bic)
  return(result)
}


criteria=sapply(p,criteria.ar)


plot(p,criteria[1,],type='p',pch='a',col='red',xlab='AR order p',ylab='Criterion',main='',
     ylim=c(min(criteria)-10,max(criteria)+10))
points(p,criteria[2,],pch='b',col='blue')