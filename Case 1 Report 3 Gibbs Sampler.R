library(pscl)
library(truncnorm)
setwd("~/Documents/Fall 2016/Sta 440")
Brain_network <- load('Brain_network.RData')
na_sub <- as.vector(unique(which(is.na(covariate), arr.ind = TRUE)[,1]))
covariate <- covariate[-na_sub, ] # remove subjects with NA in covariate
W <- W[ , , -na_sub] # remove connectivity matrices associated with those subjects
n_reg <- 68 # total number of brain regions
n_sub <- nrow(covariate) # total number of subjects

# calculate indicator for same hem
same_hem <- function(u, v){
  if((u<=34) & (v<=34)) ind = 1
  else if ((u>34) & (v>34)) ind = 1
  else ind = 0
  return (ind)
}
ind <- mapply(same_hem, rep(1:n_reg, each = n_reg), rep(1:n_reg, n_reg))

# reshape dataframe
con <- matrix(0, n_reg^2*n_sub, 9)
colnames(con) <- c('i', 'u', 'v', 'y', 'sex', 'age', 'open', 'd', 'ind')
con[, 'i'] <- rep(1:n_sub, each = n_reg^2)
con[, 'u'] <- rep(1:n_reg, each = n_reg)
con[, 'v'] <- rep(1:n_reg, n_reg)
con[, 'y'] <- as.vector(W)
con[con[, 'y'] > 0, 'y'] = 1 # binary connectivity
con[, 'sex'] <- rep(covariate$Sex, each = n_reg^2)
con[, 'age'] <- rep(covariate$Age, each = n_reg^2)
con[, 'open'] <- rep(covariate$Openness, each = n_reg^2)
con[, 'd'] <- rep(d, n_sub)
con[, 'ind'] <- rep(ind, n_sub)
con <- con[-which(con[,'u']<=con[,'v']), ] # remove upper triangles of matrices of W

#initialize variables
samples<-5
n<-112
theta_k<-0
psi_k<-1
lambda_j<-0
phi_j<-1
a_j<-1
b_j<-1
sigma_squared<-1
M<-255136
mu<-matrix(nrow=3, ncol=samples)
mu[,1]<-rnorm(3,lambda_j,phi_j) #first mu iteration
Sigma<-matrix(nrow=3, ncol=samples)
Sigma[,1]<-rigamma(3,a_j,b_j) #first Sigma iteration
Beta<-matrix(nrow=2, ncol=samples)
Beta[,1]<-rnorm(2,theta_k,psi_k) #first Beta iteration 
b<-array(0,dim=c(n,3,samples))
b[,,1]<-rnorm(3*n,mu[,1],Sigma[,1]) #first b iteration
#nu values
nu_array<-matrix(c(seq(1,M)))
nu<-function(m){
  i<-con[m,1]; sex<-con[m,5]; age<-con[m,6]; open<-con[m,7]; ind<-con[m,9]
  nu_iuv<-b[i,1,1]+b[i,2,1]*age+b[i,3,1]*open+Beta[1,1]*sex+Beta[2,1]*ind
  return(nu_iuv)
}
nu_array<-apply(nu_array,1,function(x) nu(x))
#z values
z_array<-matrix(c(seq(1,M)))
z<-function(m){
  return(rnorm(1,nu_array[m],1))
}
z_array<-apply(z_array,1,function(x) z(x))
x_array<-con[,5]

#Gibbs sampler
for(s in 2:samples){
  #update mu
  for(j in 1:3){
    mu[j,s]<-rnorm(1,(lambda_j*Sigma[j,s-1]+phi_j*sum(b[,j,s-1]))/(Sigma[j,s-1]+n*phi_j),(phi_j*Sigma[j,s-1])/(Sigma[j,s-1]+n*phi_j))
  }
  #update Sigma
  for(j in 1:3){
    Sigma[j,s]<-rigamma(1,a_j+(n/2),b_j+(1/2)*sum((b[,j,s-1]-mu[j,s])^2))
  }
  #update Beta
  for(k in 1:2){
    D_m_array<-z_array-(nu_array-x_array)
    Beta[k,s]<-rnorm(1,(theta_k+psi_k*sum(x_array*D_m_array))/(1+psi_k*sum(x_array^2)),psi_k/(1+psi_k*sum(x_array^2)))
  }
  #update b
  #for(j in 1:3){
  #  for(i in 1:n){
  #    b[i,j,s]<-rnorm(1,(mu[j,s]*sigma_squared+Sigma[j,s]*sum)
  #  }
  #}
  #update nu
  nu_array<-matrix(c(seq(1,M)))
  nu<-function(m){
    i<-con[m,1]; sex<-con[m,5]; age<-con[m,6]; open<-con[m,7]; ind<-con[m,9]
    nu_iuv<-b[i,1,s]+b[i,2,s]*age+b[i,3,s]*open+Beta[1,s]*sex+Beta[2,s]*ind
    return(nu_iuv)
  }
  nu_array<-apply(nu_array,1,function(x) nu(x))
  #update z
  z_array<-matrix(c(seq(1,M)))
  z<-function(m){
    if(y==1){
      return(rtruncnorm(1,a=0,mean=nu_array[m],sd=1))
    }
    if(y==0){
      return(rtruncnorm(1,b=0,mean=nu_array[m],sd=1))
    }
  }
  z_array<-apply(z_array,1,function(x) z(x))
}