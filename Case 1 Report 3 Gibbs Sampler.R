library(pscl)
setwd("~/Documents/Fall 2016/Sta 440")
Brain_network <- load('Brain_network.RData')
samples<-5
n<-114
theta_k<-0
psi_k<-1
lambda_j<-0
phi_j<-1
a_j<-1
b_j<-1
sigma_squared<-1
D<-array(0,dim=c(n,2016,samples))
Beta<-matrix(nrow=2, ncol=samples)
Beta[,1]<-rnorm(2,theta_k,psi_k) #first Beta iteration 
mu<-matrix(nrow=3, ncol=samples)
mu[,1]<-rnorm(3,lambda_j,phi_j) #first mu iteration
Sigma<-matrix(nrow=3, ncol=samples)
Sigma[,1]<-rigamma(3,a_j,b_j) #first Sigma iteration
b<-array(0,dim=c(n,3,samples))
b[,,1]<-rnorm(3*n,mu[,1],Sigma[,1]) #first b iteration
z<-matrix(nrow=2016*n, ncol=samples)

for(s in 2:samples){
  #update D (not done)
  for(i in 1:n){
    uv<-matrix(nrow=1, ncol=2016)
    for(u in 1:64){
      for(v in (u+1):64){
        
      }
    }
  }
  #update mu
  for(j in 1:3){
    mu[j,s]<-rnorm(1,(lambda_j*Sigma[j,s-1]+phi_j*sum(b[,j,s-1]))/(Sigma[j,s-1]+n*phi_j),(phi_j*Sigma[j,s-1])/(Sigma[j,s-1]+n*phi_j))
  }
  #update Sigma
  for(j in 1:3){
    Sigma[j,s]<-rigamma(1,a_j+(n/2),b_j+(1/2)*sum((b[,j,s-1]-mu[j,s])^2))
  }
  #update b (not done)
  for(j in 1:3){
    for(i in 1:n){
      b[i,j,s]<-rnorm(1,(mu[j,s]*sigma_squared+Sigma[j,s]*sum)
    }
  }
}