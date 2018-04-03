#######################################
## correlated particle run for SV
#
# N: number of particles
# y: data (length L)
# U: pool of RVs [2N x L]
# par: parameters (phi,sigma,gamma)
# prior: priors (prior mn, prior sd) for x_1
#

corr_filterSV = function(N, y, U, par, prior){
  phi=par[1];sig=par[2];gam=par[3]
  L = length(y)
  M = vector(); V = vector()
  particles = (sig/sqrt(1-phi^2))*U[1:N,1] # initial sample from stationary distribution
  
  
  particles = sort(particles)
  
  w = dnorm(y[1],0,exp(gam+particles))
  
  logl = log(mean(w))
  
  w = w/sum(w)
  M[1] = sum(w*particles)
  V[1] = sum(w*particles^2)-M[1]^2
  
  
  ind = rep(0, N)
  for(t in 2:L){
    
    inv.dist = cumsum(w/sum(w))
    samp = pnorm(U[(N+1):(2*N),t-1])      # resampling using U variables
    for(i in 1:N){
      ind[i]= min(which(inv.dist>samp[i]))
    }
    
    
    particles = particles[ind]
    
    particles = phi*particles+sig*U[1:N,t] # propagation
    
    
    w = dnorm(y[t], 0, exp(gam+particles))
    logl = logl+log(mean(w))
    
    w = w/sum(w)
    M[t] = sum(w*particles)
    V[t] = sum(w*particles^2)-M[t]^2
  }
  return(list(mean = M, var = V, l = logl))
}





###########################################################################
#####
#SV simulation
#   n is number of time-points
#   par = (phi,sigma,gamma)
#
#   uses stationary distribution of state-process as the prior
#
###

SVsim=function(n,par){
  
  x=rep(0,n)
  x[1]=rnorm(1,0,par[2]/sqrt(1-par[1]^2)) ##stationary distribution
  for(i in 2:n) x[i]=x[i-1]*par[1]+rnorm(1,0,par[2])
  y=rnorm(n,0,exp(par[3]+x))
  return(list(x=x,y=y))
}

