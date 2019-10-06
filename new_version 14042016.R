# Packages
library(tmvtnorm)

take.sample = function(A,B,indx,indy){
  s = sapply( 1:length(indx), function(i) A[indx[i],]%*%B[indy[i],])
  return(s)
}

# Simulation of the data
n= m1 = 100
p= m2= 100
n = 0.8*m1*m2
Kreal = 2
Ureal = matrix(data=rnorm(Kreal*m1,0,20/sqrt(m1)),nr=m1,nc=Kreal)
Vreal = matrix(data=rnorm(Kreal*m2,0,20/sqrt(m2)),nr=m2,nc=Kreal)
Mreal = x=  Ureal%*%t(Vreal)

ix=seq(m1*m2)
missfrac=0.2
imiss=sample(ix,m1*m2*missfrac,replace=FALSE)
xna=Mreal+matrix(rnorm(m1*m2, sd=1),m1,m2)
xna[imiss]=NA

I=row(xna)[!is.na(xna)]
J=col(xna)[!is.na(xna)]
Y=xna[!is.na(xna)]

# Parameters

MCMC = 200
L = 50
tau = 1/2
lambda = n/4

# Preparation of the computations
K = 10
U = list()
V = list()
for (l in 1:K){
  U[[l]] <- matrix(data=rnorm(l*m1,0,20/sqrt(m1)),nr=m1,nc=l)
  V[[l]] <- matrix(data=rnorm(l*m2,0,20/sqrt(m2)),nr=m2,nc=l)
}
M = matrix(0, nr = m1, nc = m2)
select = 3


for (mcmc in 1:MCMC)  {
  for (rank in 1:K){
    for (i in 1:m1){
      Vi = V[[rank]][J[I==i],]
      Yi = Y[I==i]
      Preci = t(Vi)%*%Vi
      Vari = solve(Preci) 
      Mi = c(Yi%*%Vi%*%Vari)
      U[[rank]][i,] = rtmvnorm(n=1,mean=Mi,sigma=Vari*n/2/lambda,
                               lower=rep(-L,rank),upper=rep(L,rank),algorithm="gibbs",
                               burn.in.samples=10,start.value=U[[rank]][i,])
    }
    for (j in 1:m2){
      Uj = U[[rank]][I[J==j],]
      Yj = Y[J==j]
      Precj = t(Uj)%*%Uj
      Varj = solve(Precj) 
      Mj = c(Yj%*%Uj%*%Varj)
      V[[rank]][j,] = rtmvnorm(n=1,mean=Mj,sigma=Varj*n/2/lambda,
                               lower=rep(-L,rank),upper=rep(L,rank),algorithm="gibbs",
                               burn.in.samples=10,start.value=V[[rank]][j,])
    }
  }
  #update rank
  if (select ==1){
    select.can = sample(c(select, select +1),1,replace = T, prob = c(0.5,0.5))
    if (select.can != select){
      sum.g = sum( (Y - take.sample(U[[select.can]], V[[select.can]], I, J) )^2 
                   - (Y - take.sample(U[[select]], V[[select]], I, J) )^2)  
      ap.select = (2/3)*exp( -lambda*sum.g)*tau^(-select + select.can)
      if (runif(1) <= ap.select) select = select.can
    }
  } else if (select == K){
    select.can = sample(c(select -1, select),1,replace = T, prob = c(0.5,0.5))
  if (select.can != select){
      sum.g = sum( (Y - take.sample(U[[select.can]], V[[select.can]], I, J) )^2 
                   - (Y - take.sample(U[[select]], V[[select]], I, J) )^2)  
      ap.select =  (3/2)*exp( -lambda*sum.g)*tau^(-select + select.can)
      if (runif(1) <= ap.select) select = select.can
    }
  } else { 
    select.can = sample(c(select -1,select, select +1),1,replace = T, prob = c(1/3,1/3,1/3))
    if (select.can != select){
      sum.g = sum( (Y - take.sample(U[[select.can]], V[[select.can]], I, J) )^2 
                   - (Y - take.sample(U[[select]], V[[select]], I, J) )^2)  
      ap.select =  exp( -lambda*sum.g)*tau^( -select + select.can)
      if (runif(1) <= ap.select) select = select.can
    }}
	M = U[[select]]%*%t(V[[select]])*1/mcmc + M*(1-1/mcmc)
  print(select)
}
mean((Mreal-M)^2)



