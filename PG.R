## Particle Gibbs (PG)
#   M: number of Gibbs iterations
#   y: data
#   prior: Inverse-Gamma(a,b)
#   N: number of particles
#   sig2.init: initial sigma^2




pg = function(M, y, prior, N, sig2.init){
  alpha = prior[1];beta = prior[2]
  L = length(y)
  sig2 = rep(0,M)
  sig2[1] = sig2.init
  X = matrix(0,M,L)
  
  init.state = cpf_as(y, sig2[1], N, cond = F)
  particles = init.state$x;w=init.state$w[,L]
  J = sample(1:N, prob = w, size = 1)
  X[1,] = particles[J,]
  
  for(k in 2:M){
    a.p = alpha+L/2
    b.p = beta +0.5*( (1-phi^2)*X[k-1,1]+sum((X[k-1,2:L]-phi*X[k-1,1:(L-1)])^2))
    sig2[k]=rinvgamma(1,a.p,b.p)
    cur.state = cpf_as(y, sig2[k], N,x.c = X[k-1,], cond = T)
    particles = cur.state$x;w=cur.state$w[,L]
    J = sample(1:N, prob = w, size = 1)
    X[k,] = particles[J,]
    #print(k)
  }
  return(list('sig2' = sig2, 'X'=X ))
}

## Conditional particle filter
# x.c:  conditional reference path
# cond: conditioning (T/F)

cpf = function(y, sig2, N, x.c, cond){
  L = length(y)
  phi = 0.9
  sig = sqrt(sig2)
  tau = sig/sqrt(1+phi^2)
  gam = 1
  x = matrix(0,N,L)
  a = matrix(0,N,L)
  w = matrix(0,N,L)
  x[,1] = rnorm(N, 0, tau)
  if(cond){
    x[N,1] = x.c[1]
  }
  
  for(t in 1:L){
    if(t!=1){
      particles = phi*x[,t-1]+sig*rnorm(N)
      ind = sample(1:N, prob = w[,t-1], size = N, replace = T)
      x[,t] = particles[ind]
      if(cond){
        x[N,t] = x.c[t]
        m = exp(-1/(2*sig2)*(x.c[t]-particles)^2)
        ind[N] = N #reference trajectory index
      }
      
      a[,t]=ind;
    }
    logweights = dnorm(y[t], 0, exp(gam+x[,t]), log = T)
    const = max(logweights)
    weights = exp(logweights-const)
    w[,t] =weights/sum(weights)
  }
  ind = a[,L]
  for(t in (L-1):1){
    x[,t] = x[ind,t]
    ind = a[ind,t]
  }
  return(list('x'=x, 'w'=w))
}
